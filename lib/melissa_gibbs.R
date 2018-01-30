##---------------------------------------
# Load libraries
##---------------------------------------
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))

#'
#' EM algorithm
#'
.scbpr_EM <- function(x, H, reg_ind, K = 2, pi_k, w, basis, lambda = 1/6,
                      em_max_iter = 100, epsilon_conv = 1e-05, opt_method = "CG",
                      opt_itnmax = 50, is_parallel = TRUE, no_cores = NULL,
                      is_verbose = FALSE){

    #
    # Optimize a promoter regions across cells, which are weighted by the
    # responsibilities of belonging to each cluster.
    #
    optim_regions <- function(x, H, w, K, opt_method = opt_method, opt_itnmax,
                              post_prob, lambda){
        covered_ind <- which(!is.na(H))
        if (is.vector(w)) { w <- matrix(w, ncol = K) }
        for (k in 1:K) {  # For each cluster k
            # TODO: How to handle empty regions???
            w[, k] <- optim(par = w[, k], fn = BPRMeth::sum_weighted_bpr_lik,
                            gr = BPRMeth::sum_weighted_bpr_grad,
                            method = opt_method, control = list(maxit = opt_itnmax),
                            X_list = x[covered_ind], H_list = H[covered_ind],
                            r_nk = post_prob[covered_ind, k], lambda = lambda,
                            is_nll = TRUE)$par
        }
        return(w)
    }

    I <- length(x)      # Number of cells
    N <- length(x[[1]]) # Number of regions
    M <- basis$M + 1    # Number of basis functions
    NLL <- 1e+100       # Initialize and store NLL for each EM iteration
    n <- 0

    # Matrices / Lists for storing results
    w_pdf     <- matrix(0, nrow = I, ncol = K)  # Store weighted PDFs
    post_prob <- matrix(0, nrow = I, ncol = K)  # Hold responsibilities
    w_tmp     <- array(data = 0, dim = c(N, M, K))

    if (is_parallel) {
        # Create cluster object
        cl <- parallel::makeCluster(no_cores)
        doParallel::registerDoParallel(cl)
    }

    # Run EM algorithm until convergence
    for (t in 1:em_max_iter) {
        # TODO: Handle empty clusters!!!
        # TODO: Handle empty clusters!!!
        ## ---------------------------------------------------------------
        # Compute weighted pdfs for each cluster
        for (k in 1:K) {
            # Apply to each cell and only to regions with CpG coverage
            w_pdf[, k] <- log(pi_k[k]) + vapply(X = 1:I, FUN = function(i)
                sum(vapply(X = reg_ind[[i]], FUN = function(y)
                    bpr_log_likelihood(w = w[y, , k], X = x[[i]][[y]],
                                       H = H[[i]][[y]], lambda = lambda,
                                       is_nll = FALSE),
                    FUN.VALUE = numeric(1), USE.NAMES = FALSE)),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)
        }
        # Use the logSumExp trick for numerical stability
        Z <- apply(w_pdf, 1, log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        post_prob <- exp(w_pdf - Z)
        NLL  <- c(NLL, (-1) * sum(Z))   # Evaluate NLL
        # M-Step -----------------------------------------------
        #
        # Compute sum of posterior probabilities for each cluster
        I_k <- colSums(post_prob)
        # Update mixing proportions for each cluster
        pi_k <- I_k / I
        # Update basis function coefficient matrix w
        # If parallel mode is ON
        if (is_parallel) {
            # Parallel optimization for each region n
            res_out <- foreach::"%dopar%"(obj = foreach::foreach(n = 1:N),
              ex  = {out <- optim_regions(x = lapply(x, "[[", n),
                                          H = lapply(H, "[[", n), w = w[n, , ], K = K,
                                          opt_method = opt_method, opt_itnmax = opt_itnmax,
                                          post_prob = post_prob, lambda = lambda) })
            for (k in 1:K) {
                tmp <- sapply(res_out, function(x) x[, k])
                if (is.matrix(tmp)) { w_tmp[, , k] <- t(tmp) }
                else {w_tmp[, 1, k] <- tmp }
            }
            w <- w_tmp
        }else{
            # Sequential optimization for each region n
            res_out <- foreach::"%do%"(obj = foreach::foreach(n = 1:N),
               ex  = {out <- optim_regions(x = lapply(x, "[[", n),
                                           H = lapply(H, "[[", n), w = w[n, , ], K = K,
                                           opt_method = opt_method, opt_itnmax = opt_itnmax,
                                           post_prob = post_prob, lambda = lambda) })
            for (k in 1:K) {
                tmp <- sapply(res_out, function(x) x[, k])
                if (is.matrix(tmp)) { w_tmp[, , k] <- t(tmp) }
                else {w_tmp[, 1, k] <- tmp }
            }
            w <- w_tmp
        }

        if (is_verbose) {
            cat("\r", "It: ", t, "NLL:\t", NLL[t + 1], "\tDiff:\t", NLL[t] - NLL[t + 1])
        }
        if (NLL[t + 1] > NLL[t]) { message("NLL increases!\n"); break; }
        # Check for convergence
        if (NLL[t] - NLL[t + 1] < epsilon_conv) { break }
    }
    if (is_parallel) {
        parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
    }
    # Check if EM converged in the given maximum iterations
    if (t == em_max_iter) { warning("EM did not converge!\n") }
    obj <- structure(list(K = K, N = N, w = w, pi_k = pi_k, lambda = lambda,
                          em_max_iter = em_max_iter, opt_method = opt_method,
                          opt_itnmax = opt_itnmax, NLL = NLL,
                          basis = basis, post_prob = post_prob),
                     class = "scbpr_EM")
    return(obj)
}



# Internal function to make all the appropriate type checks.
.do_scEM_checks <- function(x, H, reg_ind, K, pi_k = NULL, w = NULL, basis,
                            lambda = 1/6, use_kmeans = TRUE, em_init_nstart = 5,
                            em_init_max_iter = 10, epsilon_conv = 1e-04,
                            opt_method = "CG", opt_itnmax = 30, init_opt_itnmax = 20,
                            is_parallel = TRUE, no_cores = NULL, is_verbose = TRUE){
    I <- length(x)
    N <- length(x[[1]])
    M <- basis$M + 1

    if (is.null(w)) {
        ww <- array(data = rnorm(N*M*I, 0, 0.01), dim = c(N, M, I))
        w_init <- rep(0.5, M)
        for (i in 1:I) {
            # Compute regression coefficients using MLE
            ww[reg_ind[[i]], ,i] <- infer_profiles_mle(X = x[[i]][reg_ind[[i]]],
               model = "bernoulli", basis = basis, H = H[[i]][reg_ind[[i]]], w = w_init,
               lambda = lambda, opt_method = opt_method, opt_itnmax = init_opt_itnmax,
               is_parallel = FALSE, no_cores = no_cores)$W
        }

        # Transform to long format to perform k-means
        W_opt <- matrix(0, nrow = I, ncol = N * M)
        for (i in 1:I) { W_opt[i, ] <- as.vector(ww[,,i]) }
        w <- array(data = 0, dim = c(N, M, K))
        NLL_prev <- 1e+120
        optimal_w = optimal_pi_k <- NULL

        # Run 'mini' EM algorithm to find optimal starting points
        if (is_parallel) {
            em_res <- mclapply(X = 1:em_init_nstart, FUN = function(t){
                if (use_kmeans) {
                    # Use Kmeans with random starts
                    cl <- stats::kmeans(W_opt, K, nstart = 1)
                    # Get the mixture components
                    C_n <- cl$cluster
                    # TODO: Check that k-means does not return empty clusters..
                    # Sample randomly one point from each cluster as initial centre
                    for (k in 1:K) { w[, ,k] <- ww[, , sample(which(C_n == k), 1)] }
                    # Mixing proportions
                    if (is.null(pi_k)) { pi_k <- as.vector(table(C_n) / I ) }
                }else{
                    w <- array(data = ww[, ,sample(I, K)], dim = c(N, M, K))
                    if (is.null(pi_k)) { pi_k <- rep(1/K, K) }
                }
                # Run mini EM
                em <- .scbpr_EM(x = x, H = H, reg_ind = reg_ind, K = K, pi_k = pi_k,
                                w = w, basis = basis, lambda = lambda,
                                em_max_iter = em_init_max_iter, epsilon_conv = epsilon_conv,
                                opt_method = opt_method,
                                opt_itnmax = opt_itnmax, is_parallel = FALSE,
                                no_cores = NULL, is_verbose = is_verbose)

                return(em)
            }, mc.cores = no_cores)
        }else{
            em_res <- lapply(X = 1:em_init_nstart, FUN = function(t){
                if (use_kmeans) {
                    # Use Kmeans with random starts
                    cl <- stats::kmeans(W_opt, K, nstart = 1)
                    # Get the mixture components
                    C_n <- cl$cluster
                    # TODO: Check that k-means does not return empty clusters..
                    # Sample randomly one point from each cluster as initial centre
                    for (k in 1:K) { w[, ,k] <- ww[, , sample(which(C_n == k), 1)] }
                    # Mixing proportions
                    if (is.null(pi_k)) { pi_k <- as.vector(table(C_n) / I ) }
                }else{
                    w <- array(data = ww[, ,sample(I, K)], dim = c(N, M, K))
                    if (is.null(pi_k)) { pi_k <- rep(1/K, K) }
                }
                # Run mini EM
                em <- .scbpr_EM(x = x, H = H, reg_ind = reg_ind, K = K, pi_k = pi_k,
                                w = w, basis = basis, lambda = lambda,
                                em_max_iter = em_init_max_iter, epsilon_conv = epsilon_conv,
                                opt_method = opt_method,
                                opt_itnmax = opt_itnmax, is_parallel = FALSE,
                                no_cores = NULL, is_verbose = is_verbose)
                return(em)
            } )
        }
        for (t in 1:em_init_nstart) {
            # Check if NLL is lower and keep the optimal params
            NLL_cur <- utils::tail(em_res[[t]]$NLL, n = 1)
            if (NLL_cur < NLL_prev) {
                # TODO:: Store optimal pi_k from EM
                optimal_w <- em_res[[t]]$w
                optimal_pi_k <- em_res[[t]]$pi_k
                NLL_prev <- NLL_cur
            }
        }
    }
    if (is.null(pi_k)) { optimal_pi_k <- rep(1 / K, K) }
    if (length(w[1,,1]) != (basis$M + 1) ) {
        stop("Coefficient vector should be M+1, M: number of basis functions!")
    }
    return(list(w = optimal_w, basis = basis, pi_k = optimal_pi_k))
}



#' Gibbs sampling algorithm for Melissa model
#'
#' \code{melissa_gibbs} implements the Gibbs sampling algorithm for
#' performing clustering of single cells based on their DNA methylation
#' profiles, where the observation model is the Bernoulli distributed Probit
#' Regression likelihood.
#'
#' @param x A list of length I, where I are the total number of cells. Each
#'   element of the list contains another list of length N, where N is the total
#'   number of genomic regions. Each element of the inner list is an L x 2
#'   matrix of observations, where 1st column contains the locations and the 2nd
#'   column contains the methylation level of the corresponding CpGs.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
#' @param w A N x M x K array, where each column contains the basis function
#'   coefficients for the corresponding cluster.
#' @param basis A 'basis' object. E.g. see \code{\link{create_rbf_object}}
#' @param w_0_mean The prior mean hyperparameter for w
#' @param w_0_cov The prior covariance hyperparameter for w
#' @param dir_a The Dirichlet concentration parameter, prior over pi_k
#' @param lambda The complexity penalty coefficient for penalized regression.
#' @param gibbs_nsim Argument giving the number of simulations of the Gibbs
#'   sampler.
#' @param gibbs_burn_in Argument giving the burn in period of the Gibbs sampler.
#' @param inner_gibbs Logical, indicating if we should perform Gibbs sampling to
#'   sample from the augmented BPR model.
#' @param gibbs_inner_nsim Number of inner Gibbs simulations.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param is_verbose Logical, print results during EM iterations
#'
#' @importFrom stats rmultinom rnorm
#' @importFrom MCMCpack rdirichlet
#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm mvtnorm::rmvnorm
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
melissa_gibbs <- function(x, K = 2, pi_k = rep(1/K, K), w = NULL, basis = NULL,
                          w_0_mean = NULL, w_0_cov = NULL, dir_a = rep(1, K),
                          lambda = 1/2, gibbs_nsim = 1000, gibbs_burn_in = 200,
                          inner_gibbs = FALSE, gibbs_inner_nsim = 50,
                          is_parallel = TRUE, no_cores = NULL, is_verbose = FALSE){

    # Check that x is a list object
    assertthat::assert_that(is.list(x))
    assertthat::assert_that(is.list(x[[1]]))
    I <- length(x)      # Extract number of cells
    N <- length(x[[1]]) # Extract number of promoter regions
    # Number of parallel cores
    no_cores <- BPRMeth:::.parallel_cores(no_cores = no_cores,
                                          is_parallel = is_parallel,
                                          max_cores = N)
    if (is.null(basis)) { basis <- create_rbf_object(M = 3) }
    M <- basis$M + 1    # Number of coefficient parameters

    # Initialize priors over the parameters
    if (is.null(w_0_mean)) { w_0_mean <- rep(0, M) }
    if (is.null(w_0_cov)) { w_0_cov <- diag(4, M) }
    prec_0 <- solve(w_0_cov)          # Invert covariance matrix to get the precision matrix
    w_0_prec_0 <- prec_0 %*% w_0_mean # Compute product of prior mean and prior precision matrix

    # Matrices / Lists for storing results
    w_pdf        <- matrix(0, nrow = I, ncol = K)  # Store weighted PDFs
    post_prob    <- matrix(0, nrow = I, ncol = K)  # Hold responsibilities
    C            <- matrix(0, nrow = I, ncol = K)  # Mixture components
    C_prev       <- C                              # Keep previous components
    C_matrix     <- matrix(0, nrow = I, ncol = K)  # Total mixture components
    NLL          <- vector(mode = "numeric", length = gibbs_nsim)
    NLL[1]       <- 1e+100
    H = y = z = V <- list()
    for (k in 1:K) {
        H[[k]] <- vector("list", N) # List of concatenated design matrices
        y[[k]] <- vector("list", N) # List of observed methylation data
        z[[k]] <- vector("list", N) # List of auxiliary latent variables
        V[[k]] <- vector("list", N) # List of posterior variances
    }
    len_y <- matrix(0, nrow = K, ncol = N) # Total CpG observations per region
    sum_y <- matrix(0, nrow = K, ncol = N) # Total methylated CpGs per region

    # Store mixing proportions draws
    pi_draws <- matrix(NA_real_, nrow = gibbs_nsim, ncol = K)
    # Store BPR coefficient draws
    w_draws <- array(data = 0, dim = c(gibbs_nsim - gibbs_burn_in, N, M , K))

    # List of genes with no coverage for each cell
    region_ind <- lapply(X = 1:I, FUN = function(n) which(!is.na(x[[n]])))
    # Pre-compute the design matrices H for efficiency: each entry is a cell
    if (is_parallel) { des_mat <- parallel::mclapply(X = 1:I, FUN = function(n)
        init_design_matrix(basis = basis, X = x[[n]], coverage_ind = region_ind[[n]]),
        mc.cores = no_cores)
    }else {des_mat <- lapply(X = 1:I, FUN = function(n)
        init_design_matrix(basis = basis, X = x[[n]], coverage_ind = region_ind[[n]]))
    }

    # TODO: Initialize w in a sensible way via mini EM
    if (is.null(w)) {
        # Perform checks for initial parameter values
        no_cores <- BPRMeth:::.parallel_cores(no_cores = no_cores,
                                              is_parallel = is_parallel,
                                              max_cores = 5)

        out <- .do_scEM_checks(x = x, H = des_mat, reg_ind = region_ind, K = K, pi_k = NULL,
                               w = w, basis = basis, lambda = lambda, em_init_nstart = 2,
                               em_init_max_iter = 10, opt_itnmax = 15,
                               init_opt_itnmax = 10, is_parallel = is_parallel,
                               no_cores = no_cores, is_verbose = FALSE)
        w <- out$w; pi_k <- out$pi_k
    }
    pi_draws[1, ] <- pi_k

    if (is_verbose) { message("Starting Gibbs sampling...") }
    # Show progress bar
    pb <- txtProgressBar(min = 1, max = gibbs_nsim, style = 3)
    ##---------------------------------
    # Start Gibbs sampling
    ##---------------------------------
    for (t in 2:gibbs_nsim) {
        empty_C <- vector("integer", K)
        ## ---------------------------------------------------------------
        # Compute weighted pdfs for each cluster
        for (k in 1:K) {
            # Apply to each cell and only to regions with CpG coverage
            w_pdf[, k] <- log(pi_k[k]) + vapply(X = 1:I, FUN = function(i)
                sum(vapply(X = region_ind[[i]], FUN = function(y)
                    bpr_log_likelihood(w = w[y, , k], X = x[[i]][[y]],
                                       H = des_mat[[i]][[y]], lambda = lambda,
                                       is_nll = FALSE),
                    FUN.VALUE = numeric(1), USE.NAMES = FALSE)),
                FUN.VALUE = numeric(1), USE.NAMES = FALSE)
        }
        # Use the logSumExp trick for numerical stability
        Z <- apply(w_pdf, 1, log_sum_exp)
        # Get actual posterior probabilities, i.e. responsibilities
        post_prob <- exp(w_pdf - Z)
        NLL[t] <- -sum(Z)   # Evaluate NLL

        ## -------------------------------------------------------------------
        # Draw mixture components for ith simulation
        # Sample one point from a Multinomial i.e. ~ Discrete
        for (i in 1:I) { C[i, ] <- rmultinom(n = 1, size = 1, post_prob[i, ]) }
        # ## -------------------------------------------------------------------
        # TODO: Should we keep all data
        if (t > gibbs_burn_in) { C_matrix <- C_matrix + C }

        ## -------------------------------------------------------------------
        # Update mixing proportions using updated cluster component counts
        Ci_k <- colSums(C)
        if (is_verbose) {cat("\r", Ci_k) }
        pi_k <- as.vector(MCMCpack::rdirichlet(n = 1, alpha = dir_a + Ci_k))
        pi_draws[t, ] <- pi_k

        # Matrix to keep promoters with no CpG coverage
        empty_region <- matrix(0, nrow = N, ncol = K)
        for (k in 1:K) {
            # Which cells are assigned to cluster k
            C_idx <- which(C[, k] == 1)
            # TODO: Handle cases when we have empty clusters...
            if (length(C_idx) == 0) {
                if (is_verbose) { message("Warning: Empty cluster...") }
                empty_C[k] <- 1
                next
            }
            # Check if current clusters ids are not equal to previous ones
            if (!identical(C[, k], C_prev[, k])) {
                if (is_verbose) { message(t, ": Not identical in cluster ", k) }
                # Iterate over each promoter region
                for (n in 1:N) {
                    # Initialize empty vector for observed methylation data
                    yy <- vector(mode = "integer")
                    # Concatenate the nth promoter from all cells in cluster k
                    tmp <- lapply(des_mat, "[[", n)[C_idx]

                    # TODO: Is this NULL or NA???
                    tmp <- do.call(rbind, tmp[!is.na(tmp)])
                    # TODO: Check when we have empty promoters....
                    if (is.null(tmp)) {
                        H[[k]][[n]] <- NA
                        empty_region[n, k] <- 1
                    }else{
                        H[[k]][[n]] <- tmp
                        # Obtain the corresponding methylation levels
                        for (cell in C_idx) {
                            obs <- x[[cell]][[n]]
                            if (length(obs) > 1) { yy <- c(yy, obs[, 2]) }
                        }
                        # Precompute for faster computations
                        len_y[k, n] <- length(yy)
                        sum_y[k, n] <- sum(yy)
                        y[[k]][[n]] <- yy
                        z[[k]][[n]] <- rep(NA_real_, len_y[k, n])
                        # Compute posterior variance of w_nk
                        V[[k]][[n]] <- solve(prec_0 + crossprod(H[[k]][[n]], H[[k]][[n]]))
                    }
                }
            }
            for (n in 1:N) {
                # In case we have no CpG data in this promoter
                if (is.vector(H[[k]][[n]])) { next }
                # Perform Gibbs sampling on the augmented BPR model
                if (inner_gibbs & t > 4) {
                    w_inner <- matrix(0, nrow = gibbs_inner_nsim, ncol = M)
                    w_inner[1, ] <- w[n, , k]
                    for (tt in 2:gibbs_inner_nsim) {
                        # Update Mean of z
                        mu_z <- H[[k]][[n]] %*% w_inner[tt - 1, ]
                        # Draw latent variable z from z | w, y, X
                        if (sum_y[k, n] == 0) {
                            z[[k]][[n]] <- truncnorm::rtruncnorm(len_y[k, n], mean = mu_z,
                                                                 sd = 1, a = -Inf, b = 0)
                        }else if (sum_y[k, n] == len_y[k, n]) {
                            z[[k]][[n]] <- truncnorm::rtruncnorm(len_y[k, n], mean = mu_z,
                                                                 sd = 1, a = 0, b = Inf)
                        }else{
                            z[[k]][[n]][y[[k]][[n]] == 1] <- truncnorm::rtruncnorm(sum_y[k, n],
                                mean = mu_z[y[[k]][[n]] == 1], sd = 1, a = 0, b = Inf)
                            z[[k]][[n]][y[[k]][[n]] == 0] <- truncnorm::rtruncnorm(len_y[k, n] - sum_y[k, n],
                                mean = mu_z[y[[k]][[n]] == 0], sd = 1, a = -Inf, b = 0)
                        }
                        # Compute posterior mean of w
                        Mu <- V[[k]][[n]] %*% (w_0_prec_0 + crossprod(H[[k]][[n]], z[[k]][[n]]))
                        # Draw variable \w from its full conditional: \w | z, X
                        if (M == 1) { w_inner[tt, ] <- c(rnorm(n = 1, mean = Mu, sd = V[[k]][[n]])) }
                        else {w_inner[tt, ] <- c(mvtnorm::rmvnorm(n = 1, mean = Mu, sigma = V[[k]][[n]])) }
                    }
                    if (M == 1) { w[n, , k] <- mean(w_inner[-(1:(gibbs_inner_nsim/2)), ]) }
                    else {w[n, , k] <- colMeans(w_inner[-(1:(gibbs_inner_nsim/2)), ]) }
                }else{
                    ##-------------=============-=-=-=
                    # TODO:: Should we run this twice to update the z parameter!!!
                    for (l in 1:3) {
                        # Update Mean of z
                        mu_z <- H[[k]][[n]] %*% w[n, , k]
                        # Draw latent variable z from its full conditional: z | w, y, X
                        if (sum_y[k, n] == 0) {
                            z[[k]][[n]] <- truncnorm::rtruncnorm(len_y[k, n], mean = mu_z, sd = 1, a = -Inf, b = 0)
                        }else if (sum_y[k, n] == len_y[k, n]) {
                            z[[k]][[n]] <- truncnorm::rtruncnorm(len_y[k, n], mean = mu_z, sd = 1, a = 0, b = Inf)
                        }else{
                            z[[k]][[n]][y[[k]][[n]] == 1] <- truncnorm::rtruncnorm(sum_y[k, n],
                                    mean = mu_z[y[[k]][[n]] == 1], sd = 1, a = 0, b = Inf)
                            z[[k]][[n]][y[[k]][[n]] == 0] <- truncnorm::rtruncnorm(len_y[k, n] - sum_y[k, n],
                                    mean = mu_z[y[[k]][[n]] == 0], sd = 1, a = -Inf, b = 0)
                        }
                        # Compute posterior mean of w
                        Mu <- V[[k]][[n]] %*% (w_0_prec_0 + crossprod(H[[k]][[n]], z[[k]][[n]]))
                        # Draw variable \w from its full conditional: \w | z, X
                        if (M == 1) { w[n, , k] <- c(rnorm(n = 1, mean = Mu, sd = V[[k]][[n]])) }
                        else{w[n, , k] <- c(mvtnorm::rmvnorm(n = 1, mean = Mu, sigma = V[[k]][[n]])) }
                    }
                }
            }
        }

        # For each empty promoter region, take the methyation profile of the
        # promoter regions that belong to another cluster
        for (n in 1:N) {
            clust_empty_ind <- which(empty_region[n, ] == 1)
            # No empty promoter regions
            if (length(clust_empty_ind) == 0) { next }
            # Case that should never happen with the right preprocessing step
            else if (length(clust_empty_ind) == K) {
                for (k in 1:K) { w[n, , k] <- c(mvtnorm::rmvnorm(1, w_0_mean, w_0_cov)) }
            }else{
                # TODO: Perform a better imputation approach...
                cover_ind <- which(empty_region[n, ] == 0)
                # Randomly choose a cluster to obtain the methylation profiles
                k_imp <- sample(length(cover_ind), 1)
                for (k in seq_along(clust_empty_ind)) {
                    w[n, , clust_empty_ind[k]] <- w[n, , k_imp]
                }
            }
        }

        # Handle empty clusters
        dom_C <- which.max(Ci_k)
        for (k in 1:K) {
            if (empty_C[k] == 1) {
                w[, , k] <- w[, , dom_C]
            }
        }
        C_prev <- C # Make current cluster indices same as previous
        if (t > gibbs_burn_in) {w_draws[t - gibbs_burn_in, , ,] <- w}
        setTxtProgressBar(pb,t)
    }
    close(pb)
    if (is_verbose) { message("Finished Gibbs sampling...") }

    ##-----------------------------------------------
    if (is_verbose) { message("Computing summary statistics...") }
    # Compute summary statistics from Gibbs simulation
    if (K == 1) { pi_post <- mean(pi_draws[gibbs_burn_in:gibbs_nsim, ]) }
    else {pi_post <- colMeans(pi_draws[gibbs_burn_in:gibbs_nsim, ]) }
    C_post <- C_matrix / (gibbs_nsim - gibbs_burn_in)
    w_post <- array(0, dim = c(N, M, K))
    for (k in 1:K) { w_post[, , k] <- colSums(w_draws[, , , k]) /
        (gibbs_nsim - gibbs_burn_in)  }

    # Object to keep input data
    dat <- list(K = K, N = N, I = I, M = M, basis = basis, dir_a = dir_a,
                lambda = lambda, w_0_mean = w_0_mean, w_0_cov = w_0_cov,
                gibbs_nsim = gibbs_nsim, gibbs_burn_in = gibbs_burn_in)
    # Object to hold all the Gibbs draws
    draws <- list(pi = pi_draws, w = w_draws, C = C_matrix, NLL = NLL)
    # Object to hold the summaries for the parameters
    summary <- list(pi = pi_post, w = w_post, C = C_post)
    # Create sc_bayes_bpr_fdmm object
    obj <- structure(list(summary = summary, draws = draws, dat = dat,
                          r_nk = C_post, W = w_post, pi_k = pi_post, basis = basis),
                     class = "melissa_gibbs")
    return(obj)
}
