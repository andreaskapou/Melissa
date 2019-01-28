##---------------------------------------
# Load libraries
##---------------------------------------
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(matrixcalc))

#' @name melissa_vb
#' @rdname melissa_vb
#' @aliases melissa_vb_cluster
#'
#' @title Cluster single cells methylomes using VB
#'
#' @description \code{melissa_vb} clusters single cells based on their
#'   methylation profiles on specific genomic regions, e.g. promoters, using the
#'   Variational Bayes (VB) EM-like algorithm.
#'
#' @param X The input data, which has to be a \link[base]{list} of elements of
#'   length N, where N are the total number of cells. Each element in the list
#'   contains another list of length M, where M is the total number of genomic
#'   regions, e.g. promoters. Each element in the inner list is an \code{I X 2}
#'   matrix, where I are the total number of observations. The first column
#'   contains the input observations x (i.e. CpG locations) and the 2nd columns
#'   contains the corresponding methylation level.
#' @param K Integer denoting the total number of clusters K.
#' @param basis A 'basis' object. E.g. see \code{\link{create_basis}}. If NULL,
#'   will an RBF object will be created.
#' @param delta_0 Parameter vector of the Dirichlet prior on the mixing
#'   proportions pi.
#' @param w Optional, an Mx(D)xK array of the initial parameters, where first
#'   dimension are the genomic regions M, 2nd the number of covariates D (i.e.
#'   basis functions), and 3rd are the clusters K. If NULL, will be assigned
#'   with default values.
#' @param alpha_0 Hyperparameter: shape parameter for Gamma distribution. A
#'   Gamma distribution is used as prior for the precision parameter tau.
#' @param beta_0 Hyperparameter: rate parameter for Gamma distribution. A Gamma
#'   distribution is used as prior for the precision parameter tau.
#' @param vb_max_iter Integer denoting the maximum number of VB iterations.
#' @param epsilon_conv Numeric denoting the convergence threshold for VB.
#' @param is_kmeans Logical, use Kmeans for initialization of model parameters.
#' @param vb_init_nstart Number of VB random starts for finding better
#'   initialization.
#' @param vb_init_max_iter Maximum number of mini-VB iterations.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param is_verbose Logical, print results during VB iterations.
#'
#' @return An object of class \code{melissa_vb} with the following elements:
#'   \itemize{ \item{ \code{W}: An (M+1) X K matrix with the optimized parameter
#'   values for each cluster, M are the number of basis functions. Each column
#'   of the matrix corresponds a different cluster k.} \item{ \code{W_Sigma}: A
#'   list with the covariance matrices of the posterior parmateter W for each
#'   cluster k.} \item{ \code{r_nk}: An (N X K) responsibility matrix of each
#'   observations being explained by a specific cluster. } \item{ \code{delta}:
#'   Optimized Dirichlet paramter for the mixing proportions. } \item{
#'   \code{alpha}: Optimized shape parameter of Gamma distribution. } \item{
#'   \code{beta}: Optimized rate paramter of the Gamma distribution } \item{
#'   \code{basis}: The basis object. } \item{\code{lb}: The lower bound vector.}
#'   \item{\code{labels}: Cluster assignment labels.} \item{ \code{pi_k}:
#'   Expected value of mixing proportions.} }
#'
#' @section Details: The modelling and mathematical details for clustering
#'   profiles using mean-field variational inference are explained here:
#'   \url{http://rpubs.com/cakapourani/} . More specifically: \itemize{\item{For
#'   Binomial/Bernoulli observation model check:
#'   \url{http://rpubs.com/cakapourani/vb-mixture-bpr}} \item{For Gaussian
#'   observation model check: \url{http://rpubs.com/cakapourani/vb-mixture-lr}}}
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
melissa_vb <- function(X, K = 3, basis = NULL, delta_0 = NULL, w = NULL,
                       alpha_0 = .5, beta_0 = 10, vb_max_iter = 100,
                       epsilon_conv = 1e-5, is_kmeans = TRUE,
                       vb_init_nstart = 10, vb_init_max_iter = 20,
                       is_parallel = TRUE, no_cores = 3, is_verbose = TRUE){
    assertthat::assert_that(is.list(X))  # Check that X is a list object...
    assertthat::assert_that(is.list(X[[1]]))  # and its element is a list object
    # Create RBF basis object by default
    if (is.null(basis)) { basis <- create_rbf_object(M = 3) }

    N <- length(X)         # Total number of cells
    M <- length(X[[1]])    # Total number of genomic regions
    D <- basis$M + 1       # Number of basis functions

    if (is.null(delta_0)) { delta_0 <- rep(1e-5, K) }
    # Number of parallel cores
    no_cores <- BPRMeth:::.parallel_cores(no_cores = no_cores,
                                          is_parallel = is_parallel,
                                          max_cores = N)
    # List of genes with no coverage for each cell
    region_ind <- lapply(X = 1:N, FUN = function(n) which(!is.na(X[[n]])))
    # List of cells with no coverage for each genomic region
    cell_ind <- lapply(X = 1:M, FUN = function(m) which(!is.na(lapply(X, "[[", m))))
    # Pre-compute the design matrices H for efficiency: each entry is a cell
    if (is_parallel) { H <- parallel::mclapply(X = 1:N, FUN = function(n)
        init_design_matrix(basis = basis, X = X[[n]], coverage_ind = region_ind[[n]]),
        mc.cores = no_cores)
    }else {H <- lapply(X = 1:N, FUN = function(n)
        init_design_matrix(basis = basis, X = X[[n]], coverage_ind = region_ind[[n]]))
    }
    # Extract responses y_{n}
    y <- lapply(X = 1:N, FUN = function(n) extract_y(X = X[[n]], coverage_ind = region_ind[[n]]))

    # If no initial values
    if (is.null(w)) {
        # Infer MLE profiles for each cell and region
        if (is_parallel) { w_mle <- parallel::mclapply(X = 1:N, FUN = function(n)
            infer_profiles_mle(X = X[[n]][region_ind[[n]]], model = "bernoulli",
                               basis = basis, H = H[[n]][region_ind[[n]]],
                               lambda = .5, opt_itnmax = 15)$W, mc.cores = no_cores)
        }else{w_mle <- lapply(X = 1:N, FUN = function(n)
            infer_profiles_mle(X = X[[n]][region_ind[[n]]], model = "bernoulli",
                               basis = basis, H = H[[n]][region_ind[[n]]],
                               lambda = .5, opt_itnmax = 15)$W)
        }
        # Transform to long format to perform k-means
        W_tmp <- matrix(0, nrow = N, ncol = M * D)
        ww <- array(data = rnorm(M*D*N, 0, 0.01), dim = c(M, D, N))
        for (n in 1:N) { # Iterate over each cell
            # Store optimized w to an array object (genes x basis x cells)
            ww[region_ind[[n]],,n] <- w_mle[[n]]
            W_tmp[n, ] <- c(ww[,,n])
        }
        rm(w_mle)

        # Run mini VB
        lb_prev <- -1e+120
        w <- array(data = 0, dim = c(M, D, K))
        for (t in 1:vb_init_nstart) {
            if (is_kmeans) {
                # Use Kmeans for initial clustering of cells
                cl <- stats::kmeans(W_tmp, K, nstart = 1)
                # Get the mixture components
                C_n <- cl$cluster
                # TODO: Check that k-means does not return empty clusters..
                # Sample randomly one point from each cluster as initial centre
                for (k in 1:K) { w[, ,k] <- ww[, , sample(which(C_n == k), 1)] }
            }else{# Sample randomly data centers
                w <- array(data = ww[, ,sample(N, K)], dim = c(M, D, K))
            }
            # Run mini-VB
            mini_vb <- melissa_vb_inner(H = H, y = y, region_ind = region_ind,
                                        cell_ind = cell_ind, K = K, basis = basis,
                                        w = w, delta_0 = delta_0, alpha_0 = alpha_0,
                                        beta_0 = beta_0, vb_max_iter = vb_init_max_iter,
                                        epsilon_conv = epsilon_conv, is_parallel = is_parallel,
                                        no_cores = no_cores, is_verbose = is_verbose)
            # Check if NLL is lower and keep the optimal params
            lb_cur <- utils::tail(mini_vb$lb, n = 1)
            # TODO: Best way to obtain initial params after mini-restarts
            # if (lb_cur > lb_prev) { optimal_w <- mini_vb$W; lb_prev <- lb_cur; }
            if (lb_cur > lb_prev) { optimal_w <- w; lb_prev <- lb_cur; }
        }
    }
    # Run final VB
    obj <- melissa_vb_inner(H = H, y = y, region_ind = region_ind, cell_ind = cell_ind,
                            K = K, basis = basis, w = optimal_w, delta_0 = delta_0,
                            alpha_0 = alpha_0, beta_0 = beta_0, vb_max_iter = vb_max_iter,
                            epsilon_conv = epsilon_conv, is_parallel = is_parallel,
                            no_cores = no_cores, is_verbose = is_verbose)
    # Add names to the estimated parameters for clarity
    names(obj$delta) <- paste0("cluster", 1:K)
    names(obj$pi_k) <- paste0("cluster", 1:K)
    # colnames(obj$W) <- paste0("cluster", 1:K)
    colnames(obj$r_nk) <- paste0("cluster", 1:K)
    # Get hard cluster assignments for each observation
    obj$labels <- apply(X = obj$r_nk, MARGIN = 1,
                        FUN = function(x) which(x == max(x, na.rm = TRUE)))
    # Subtract \ln(K!) from lower bound to get a better estimate
    obj$lb <- obj$lb - log(factorial(K))
    # Return object
    return(obj)
}

# Compute E[z]
.update_Ez <- function(E_z, mu, y_1, y_0){
    E_z[y_1] <- mu[y_1] + dnorm(-mu[y_1]) / (1 - pnorm(-mu[y_1]))
    E_z[y_0] <- mu[y_0] - dnorm(-mu[y_0]) / pnorm(-mu[y_0])
    return(E_z)
}


##------------------------------------------------------

# Cluster single cells
melissa_vb_inner <- function(H, y, region_ind, cell_ind, K, basis, w, delta_0,
                             alpha_0, beta_0, vb_max_iter, epsilon_conv,
                             is_parallel, no_cores, is_verbose){
    assertthat::assert_that(is.list(H)) # Check that H is a list object
    assertthat::assert_that(is.list(H[[1]])) # Check that H[[1]] is a list object
    N <- length(H)          # Number of cells
    M <- length(H[[1]])     # Number of genomic regions, i.e. genes
    D <- basis$M + 1        # Number of covariates
    assertthat::assert_that(N > 0)
    assertthat::assert_that(M > 0)
    LB <- c(-Inf)           # Store the lower bounds
    r_nk = log_r_nk = log_rho_nk <- matrix(0, nrow = N, ncol = K)
    E_ww <- vector("numeric", length = K)

    # Compute H_{n}'H_{n}
    HH = y_1 = y_0 <- vector(mode = "list", length = N)
    for (n in 1:N) {
        HH[[n]] = y_1[[n]] = y_0[[n]] <- vector(mode = "list", length = M)
        HH[[n]] <- lapply(X = HH[[n]], FUN = function(x) x <- NA)
        y_1[[n]] <- lapply(X = y_1[[n]], FUN = function(x) x <- NA)
        y_0[[n]] <- lapply(X = y_0[[n]], FUN = function(x) x <- NA)

        HH[[n]][region_ind[[n]]] <- lapply(X = H[[n]][region_ind[[n]]],
                                           FUN = function(h) crossprod(h))
        y_1[[n]][region_ind[[n]]] <- lapply(X = y[[n]][region_ind[[n]]],
                                            FUN = function(yy) which(yy == 1))
        y_0[[n]][region_ind[[n]]] <- lapply(X = y[[n]][region_ind[[n]]],
                                            FUN = function(yy) which(yy == 0))
    }
    E_z = E_zz <- y

    # Compute shape parameter of Gamma
    alpha_k <- rep(alpha_0 + M*D/2, K)
    # Mean for each cluster
    m_k <- w
    # Covariance of each cluster
    S_k <- lapply(X = 1:K, function(k) lapply(X = 1:M, FUN = function(m) solve(diag(50, D))))
    # Scale of precision matrix
    if (is.null(beta_0)) { beta_0 <- sqrt(alpha_k[1]) }
    beta_k   <- rep(beta_0, K)
    # Dirichlet parameter
    delta_k  <- delta_0
    # Expectation of log Dirichlet
    e_log_pi <- digamma(delta_k) - digamma(sum(delta_k))
    mk_Sk    <- lapply(X = 1:M, FUN = function(m) lapply(X = 1:K,
                    FUN = function(k) tcrossprod(m_k[m,,k]) + S_k[[k]][[m]]))
    # Update \mu
    mu <- vector(mode = "list", length = N)
    for (n in 1:N) {
        mu[[n]] <- vector(mode = "list", length = M)
        mu[[n]] <- lapply(X = mu[[n]], FUN = function(x) x <- NA)
        if (D == 1) {
            mu[[n]][region_ind[[n]]] <- lapply(X = region_ind[[n]],
                    FUN = function(m) c(H[[n]][[m]] %*% mean(m_k[m,,])))
        }else {
            mu[[n]][region_ind[[n]]] <- lapply(X = region_ind[[n]],
                    FUN = function(m) c(H[[n]][[m]] %*% rowMeans(m_k[m,,])))
        }
    }
    # Update E[z] and E[z^2]
    E_z  <- lapply(X = 1:N, FUN = function(n) {
        l <- E_z[[n]]; l[region_ind[[n]]] <- lapply(X = region_ind[[n]],
        FUN = function(m) .update_Ez(E_z = E_z[[n]][[m]], mu = mu[[n]][[m]],
        y_0 = y_0[[n]][[m]], y_1 = y_1[[n]][[m]])); return(l)})
    E_zz <- lapply(X = 1:N, FUN = function(n) {
        l <- E_zz[[n]]; l[region_ind[[n]]] <- lapply(X = region_ind[[n]],
        FUN = function(m) 1 + mu[[n]][[m]] * E_z[[n]][[m]]); return(l)})

    # Show progress bar
    if (is_verbose) { pb <- utils::txtProgressBar(min = 2, max = vb_max_iter, style = 3) }
    # Run VB algorithm until convergence
    for (i in 2:vb_max_iter) {
        ##-------------------------------
        # Variational E-Step
        ##-------------------------------
        if (is_parallel) {
            for (k in 1:K) {
                log_rho_nk[,k] <- e_log_pi[k] + unlist(parallel::mclapply(X = 1:N, function(n)
                    sum(sapply(X = region_ind[[n]], function(m)
                    -0.5*crossprod(E_zz[[n]][[m]]) + m_k[m,,k] %*% t(H[[n]][[m]]) %*%
                    E_z[[n]][[m]] - 0.5*matrix.trace(HH[[n]][[m]] %*% mk_Sk[[m]][[k]]))),
                    mc.cores = no_cores))
            }
        }else {
            for (k in 1:K) {
                log_rho_nk[,k] <- e_log_pi[k] + sapply(X = 1:N, function(n)
                    sum(sapply(X = region_ind[[n]], function(m)
                    -0.5*crossprod(E_zz[[n]][[m]]) + m_k[m,,k] %*% t(H[[n]][[m]]) %*%
                    E_z[[n]][[m]] - 0.5*matrix.trace(HH[[n]][[m]] %*% mk_Sk[[m]][[k]]))))
            }
        }
        # Calculate responsibilities using logSumExp for numerical stability
        log_r_nk <- log_rho_nk - apply(log_rho_nk, 1, log_sum_exp)
        r_nk     <- exp(log_r_nk)
        ##-------------------------------
        # Variational M-Step
        ##-------------------------------
        # Update Dirichlet parameter
        delta_k <- delta_0 + colSums(r_nk)
        # TODO: Compute expected value of mixing proportions
        pi_k <- (delta_0 + colSums(r_nk)) / (K * delta_0 + N)

        # Iterate over each cluster
        for (k in 1:K) {
            if (is_parallel) {
                tmp <- parallel::mclapply(X = 1:M, function(m) {
                    # Extract temporary objects
                    tmp_HH <- lapply(HH, "[[", m)
                    tmp_H  <- lapply(H, "[[", m)
                    tmp_Ez <- lapply(E_z, "[[", m)
                    # Update covariance for Gaussian
                    S_k <- solve(diag(alpha_k[k]/beta_k[k] + 1e-7, D) +
                        BPRMeth:::.add_func(lapply(X = cell_ind[[m]], FUN = function(n) tmp_HH[[n]]*r_nk[n,k])))
                    # Update mean for Gaussian
                    m_k <- S_k %*% BPRMeth:::.add_func(lapply(X = cell_ind[[m]],
                        FUN = function(n) t(tmp_H[[n]]) %*% tmp_Ez[[n]]*r_nk[n,k]))
                    # Compute E[w^Tw]
                    E_ww <- crossprod(m_k) + matrixcalc::matrix.trace(S_k)
                    return(list(S_k = S_k, m_k = m_k, E_ww = E_ww))
                }, mc.cores = no_cores)
            } else{
                tmp <- lapply(X = 1:M, function(m) {
                    # Extract temporary objects
                    tmp_HH <- lapply(HH, "[[", m)
                    tmp_H  <- lapply(H, "[[", m)
                    tmp_Ez <- lapply(E_z, "[[", m)
                    # Update covariance for Gaussian
                    # TODO: How to make this numerically stable? Does the addition have strong effects?
                    S_k <- solve(diag(alpha_k[k]/beta_k[k] + 1e-7, D) +
                        BPRMeth:::.add_func(lapply(X = cell_ind[[m]], FUN = function(n) tmp_HH[[n]]*r_nk[n,k])))
                    # Update mean for Gaussian
                    m_k <- S_k %*% BPRMeth:::.add_func(lapply(X = cell_ind[[m]],
                        FUN = function(n) t(tmp_H[[n]]) %*% tmp_Ez[[n]]*r_nk[n,k]))
                    # Compute E[w^Tw]
                    E_ww <- crossprod(m_k) + matrixcalc::matrix.trace(S_k)
                    return(list(S_k = S_k, m_k = m_k, E_ww = E_ww))
                })
            }
            # Update covariance
            S_k[[k]] <- lapply(tmp, "[[", "S_k")
            # Update mean
            m_k[,,k] <- t(sapply(tmp, "[[", "m_k"))
            # Update E[w^Tw]
            E_ww[k] <- sum(sapply(tmp, "[[", "E_ww"))
            # Update \beta_k parameter for Gamma
            beta_k[k]  <- beta_0 + 0.5*E_ww[k]
            # Check beta parameter for numerical issues
            # if (is.nan(beta_k[k]) | is.na(beta_k[k]) ) { beta_k[k] <- 1e10}
            # TODO: Does this affect model penalisation??
            if (beta_k[k] > 10*alpha_k[k]) { beta_k[k] <- 10*alpha_k[k] }
        }

        # If parallel mode
        if (is_parallel) {
            # Iterate over cells
            tmp <- parallel::mclapply(X = 1:N, FUN = function(n) {
                # Update \mu
                tmp_mu <- mu[[n]]
                tmp_mu[region_ind[[n]]] <- lapply(X = region_ind[[n]], FUN = function(m)
                    c(H[[n]][[m]] %*% rowSums(matrix(sapply(X = 1:K, FUN = function(k) r_nk[n,k]*m_k[m,,k]), ncol = K))))
                # Update E[z]
                tmp_E_z <- E_z[[n]]
                tmp_E_z[region_ind[[n]]] <- lapply(X = region_ind[[n]], FUN = function(m)
                    .update_Ez(E_z = E_z[[n]][[m]], mu = tmp_mu[[m]], y_0 = y_0[[n]][[m]], y_1 = y_1[[n]][[m]]))
                # Update E[z^2]
                tmp_E_zz <- E_zz[[n]]
                tmp_E_zz[region_ind[[n]]] <- lapply(X = region_ind[[n]], FUN = function(m) 1 + tmp_mu[[m]] * tmp_E_z[[m]])
                # Return list with results
                return(list(tmp_mu = tmp_mu, tmp_E_z = tmp_E_z, tmp_E_zz = tmp_E_zz))

                # tmp_inner <- lapply(X = region_ind[[n]], FUN = function(m) {
                #     # Update \mu
                #     tmp_mu <- c(H[[n]][[m]] %*% rowSums(matrix(sapply(X = 1:K, FUN = function(k) r_nk[n,k]*m_k[m,,k]), ncol = K)))
                #     # Update E[z]
                #     tmp_E_z <- .update_Ez(E_z = E_z[[n]][[m]], mu = tmp_mu, y_0 = y_0[[n]][[m]], y_1 = y_1[[n]][[m]])
                #     # Update E[z^2]
                #     tmp_E_zz <- 1 + tmp_mu * tmp_E_z
                #     # Return list with results
                #     return(list(tmp_mu = tmp_mu, tmp_E_z = tmp_E_z, tmp_E_zz = tmp_E_zz))
                # })
                # # Concatenate results per cell
                # tmp_mu <- mu[[n]]; tmp_mu[region_ind[[n]]] <- lapply(tmp_inner, "[[", "tmp_mu")
                # tmp_E_z <- E_z[[n]]; tmp_E_z[region_ind[[n]]] <- lapply(tmp_inner, "[[", "tmp_E_z")
                # tmp_E_zz <- E_zz[[n]]; tmp_E_zz[region_ind[[n]]] <- lapply(tmp_inner, "[[", "tmp_E_zz")
            }, mc.cores = no_cores)
        } else {
            # Iterate over cells
            tmp <- lapply(X = 1:N, FUN = function(n) {
                # Update \mu
                tmp_mu <- mu[[n]]
                tmp_mu[region_ind[[n]]] <- lapply(X = region_ind[[n]], FUN = function(m)
                    c(H[[n]][[m]] %*% rowSums(matrix(sapply(X = 1:K, FUN = function(k) r_nk[n,k]*m_k[m,,k]), ncol = K))))
                # Update E[z]
                tmp_E_z <- E_z[[n]]
                tmp_E_z[region_ind[[n]]] <- lapply(X = region_ind[[n]], FUN = function(m)
                    .update_Ez(E_z = E_z[[n]][[m]], mu = tmp_mu[[m]], y_0 = y_0[[n]][[m]], y_1 = y_1[[n]][[m]]))
                # Update E[z^2]
                tmp_E_zz <- E_zz[[n]]
                tmp_E_zz[region_ind[[n]]] <- lapply(X = region_ind[[n]], FUN = function(m) 1 + tmp_mu[[m]] * tmp_E_z[[m]])
                # Return list with results
                return(list(tmp_mu = tmp_mu, tmp_E_z = tmp_E_z, tmp_E_zz = tmp_E_zz))
            })
        }
        # Concatenate final results in correct format
        mu <- lapply(tmp, "[[", "tmp_mu")
        E_z <- lapply(tmp, "[[", "tmp_E_z")
        E_zz <- lapply(tmp, "[[", "tmp_E_zz")
        rm(tmp)
        # Update expectations over \ln\pi
        e_log_pi  <- digamma(delta_k) - digamma(sum(delta_k))
        # Compute expectation of E[a]
        E_alpha <- alpha_k / beta_k
        # TODO: Perform model selection using MLE of mixing proportions
        # e_log_pi <- colMeans(r_nk)

        # Compute mm^{T}+S
        mk_Sk <- lapply(X = 1:M, FUN = function(m) lapply(X = 1:K,
                    FUN = function(k) tcrossprod(m_k[m, , k]) + S_k[[k]][[m]]))
        ##-------------------------------
        # Variational lower bound
        ##-------------------------------
        # For efficiency we compute it every 10 iterations and surely on the maximum iteration threshold
        # TODO: Find a better way to not obtain Inf in the pnorm function
        if (i %% 10 == 0 | i == vb_max_iter) {
            if (is_parallel) {
                # TODO: Find a better way to not obtain Inf in the pnorm function
                lb_pz_qz <- sum(unlist(parallel::mclapply(X = 1:N,
                    FUN = function(n) sum(sapply(X = region_ind[[n]],
                    FUN = function(m) 0.5*crossprod(mu[[n]][[m]]) +
                    sum(y[[n]][[m]] * log(1 - (pnorm(-mu[[n]][[m]]) - 1e-10) ) + (1 - y[[n]][[m]]) *
                    log(pnorm(-mu[[n]][[m]]) + 1e-10)) - 0.5*sum(sapply(X = 1:K, FUN = function(k)
                    r_nk[n,k] * matrix.trace(HH[[n]][[m]] %*% mk_Sk[[m]][[k]]) )))),
                    mc.cores = no_cores)))
            }else {
                lb_pz_qz <- sum(sapply(X = 1:N, FUN = function(n)
                    sum(sapply(X = region_ind[[n]],
                    FUN = function(m) 0.5*crossprod(mu[[n]][[m]]) +
                    sum(y[[n]][[m]] * log(1 - (pnorm(-mu[[n]][[m]]) - 1e-10) ) + (1 - y[[n]][[m]]) *
                    log(pnorm(-mu[[n]][[m]]) + 1e-10)) - 0.5*sum(sapply(X = 1:K, FUN = function(k)
                    r_nk[n,k] * matrix.trace(HH[[n]][[m]] %*% mk_Sk[[m]][[k]]) )) ))))
            }
            lb_p_w   <- sum(-0.5*M*D*log(2*pi) + 0.5*M*D*(digamma(alpha_k) -
                log(beta_k)) - 0.5*E_alpha*E_ww)
            lb_p_c   <- sum(r_nk %*% e_log_pi)
            lb_p_pi  <- sum((delta_0 - 1)*e_log_pi) + lgamma(sum(delta_0)) -
                sum(lgamma(delta_0))
            lb_p_tau <- sum(alpha_0*log(beta_0) + (alpha_0 - 1)*(digamma(alpha_k) -
                log(beta_k)) - beta_0*E_alpha - lgamma(alpha_0))
            lb_q_c   <- sum(r_nk*log_r_nk)
            lb_q_pi  <- sum((delta_k - 1)*e_log_pi) + lgamma(sum(delta_k)) -
                sum(lgamma(delta_k))
            lb_q_w   <- sum(sapply(X = 1:M, FUN = function(m) sum(-0.5*log(sapply(X = 1:K,
                FUN = function(k) det(S_k[[k]][[m]]))) - 0.5*D*(1 + log(2*pi)))))
            lb_q_tau <- sum(-lgamma(alpha_k) + (alpha_k - 1)*digamma(alpha_k) +
                                log(beta_k) - alpha_k)
            # Sum all parts to compute lower bound
            LB <- c(LB, lb_pz_qz + lb_p_c + lb_p_pi + lb_p_w + lb_p_tau - lb_q_c -
                       lb_q_pi - lb_q_w - lb_q_tau)
            iter <- length(LB)
            # Check if lower bound decreases
            if (LB[iter] < LB[iter - 1]) { warning("Warning: Lower bound decreases!\n") }
            # Check for convergence
            if (abs(LB[iter] - LB[iter - 1]) < epsilon_conv) { break }
            # Show VB difference
            if (is_verbose) {
                if (i %% 50 == 0) {
                    cat("\r","It:\t",i,"\tLB:\t",LB[iter],"\tDiff:\t",LB[iter] - LB[iter - 1],"\n")
                }
            }
            # tmp <- c(rep(0, K))
            # for (k in 1:K) { tmp[k] <- length(which(r_nk[,k] > 0.8)) }
            # print(colSums(r_nk))
            # print(tmp)
            # print(beta_k)
            # print(LB[iter])
        }
        # Check if VB converged in the given maximum iterations
        if (i == vb_max_iter) {warning("VB did not converge!\n")}
        if (is_verbose) { utils::setTxtProgressBar(pb, i) }
    }
    if (is_verbose) { close(pb) }

    # Store the object
    obj <- list(W = m_k, W_Sigma = S_k, r_nk = r_nk, delta = delta_k,
                alpha = alpha_k, beta = beta_k, basis = basis, pi_k = pi_k,
                lb = LB, delta_0 = delta_0, alpha_0 = alpha_0, beta_0 = beta_0)
    class(obj) <- "melissa_vb"
    return(obj)
}
