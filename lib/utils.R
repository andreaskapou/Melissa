##---------------------------------------
# Load libraries
##---------------------------------------
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(clues))


#'
#' Generate encode synthetic data for the Melissa model
#'
#' @param basis Basis function object
#' @param encode_w Optimized weights of the basis functions for each cluster K
#' @param N Number of cells
#' @param M Number of genomic regions per cell
#' @param K Number of clusters
#' @param C_true Ground truth cluster labels
#' @param max_cov Maximum CpG coverage
#' @param cluster_var Cell variability across clusters
generate_encode_synth_data <- function(basis, encode_w, N = 50, M = 300, K = 2,
                                       pi_k = rep(1/K, K), C_true = NULL,
                                       max_cov = 25, cluster_var = 0.5){
    D <- basis$M + 1                       # Number of covariates
    X <- vector(mode = "list", length = N) # Keep cells
    w_nk <- array(0, dim = c(M, D, K))     # Weights for each region and cluster

    # Generate M synthetic genomic regions
    # TODO: Create more S and inverse S-shape profiles!
    region_profs <- sample(1:NCOL(encode_w), size = M, replace = TRUE,
                           prob = c(.30, .15, .20, .30, .05))
    # Create K shuffles of the data (i.e. different clusters)
    cell_clusters <- matrix(0, nrow = M, ncol = K)
    for (k in 1:K) { cell_clusters[, k] <- sample(region_profs) }

    # Generate overall clusterings C_n from mixing proportions
    if (is.null(C_true)) { C_true <- t(rmultinom(N, 1, pi_k))
    }else if (NROW(C_true) != N || NCOL(C_true) != K) {
        stop("C_true matrix dimensions do not match!")
    }
    # Which cells belong in cluster k
    idx <- list()
    # Extract cells that belong to the same cluster
    for (k in 1:K) { idx[[k]] <- which(C_true[, k] == 1) }
    # Iterate over each genomic region
    for (m in 1:M) {
        # Sample from Bernoulli to either consider this region similar or not
        # across clusters based on the variability across cells
        is_variable <- rbinom(n = 1, size = 1, prob = cluster_var)
        # Generate random weights
        if (!is_variable) { w_sample <- encode_w[, sample(K, 1)] }
        # Iterate over each cluster
        for (k in 1:K) {
            # If we have similar profile across cells for all clusters, just
            # add a small noise in profiles
            if (!is_variable) {
                w_nk[m, , k] <- w_sample + rnorm(D, mean = 0, sd = 0.05)
                # Otherwise we create a new profile for each cluster
            } else{w_nk[m, , k] <- encode_w[, cell_clusters[m, k]] }
            # Simulate data for cells belonging to cluster k
            for (n in idx[[k]]) {
                X[[n]][[m]] <- generate_bpr_data(basis = basis,w = w_nk[m,,k],
                                                 max_cov = max_cov)
            } }
    }
    return(list(X = X, W = w_nk, C_true = C_true))
}


#'
#' Generating single-cell methylation data for a given region
#'
#' @param basis Basis function object
#' @param w Weights of basis functions
#' @param max_cov Maximum CpG coverage
#' @param xmin Minimum x location relative to TSS
#' @param xmax Maximum x location relative to TSS
generate_bpr_data <- function(basis, w, max_cov = 25, xmin = -1000, xmax = 1000){
    require(BPRMeth)
    # L is the number of CpGs found in the ith region
    L <- rbinom(n = 1, size = max_cov, prob = .8)
    x <- matrix(0, nrow = L, ncol = 2)
    # Randomly sample locations for the CpGs
    obs <- sort(sample(xmin:xmax, L))
    # Scale to (-1,1)
    x[, 1] <- BPRMeth:::.minmax_scaling(data = obs, xmin = xmin, xmax = xmax,
                                        fmin = -1,fmax = 1)
    H      <- design_matrix(basis, x[, 1])$H
    p_suc  <- pnorm(H %*% w) + rnorm(NROW(H), mean = 0, sd = 0.05)
    p_suc[which(p_suc > (1 - 1e-10))] <- 1 - 1e-10
    p_suc[which(p_suc < 1e-10)] <- 1e-10
    x[, 2] <- rbinom(NROW(H), 1, p_suc)
    return(x)
}


#' @title Compute stable log-sum-exp
#'
#' @description \code{.log_sum_exp} computes the log sum exp trick for avoiding
#'   numeric underflow and have numeric stability in computations of small
#'   numbers.
#'
#' @param x A vector of observations
#'
#' @return The logs-sum-exp value
#'
#' @references
#' \url{https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/}
#'
log_sum_exp <- function(x) {
    # Computes log(sum(exp(x))
    offset <- max(x)
    return(log(sum(exp(x - offset))) + offset)
}



#' @title Initialise design matrices
#'
#' @description Give a list of observations, initialise design matrices H for
#'   computational efficiency.
#'
#' @param basis Basis object.
#' @param X Observations
#' @param coverage_ind Which observations have coverage
#'
#' @return The design matrix H
#'
init_design_matrix <- function(basis, X, coverage_ind){
    # Design matrix for each genomic region -> i.e. gene
    H <- vector(mode = "list", length = length(X))
    # Set all values to NA
    H <- lapply(H, function(x) x <- NA)
    if (length(coverage_ind) < 1) {
        stop("No coverage across all regions. This cell should be discarded!")
    }
    H[coverage_ind] <- lapply(X = X[coverage_ind], FUN = function(y)
        design_matrix(obj = basis, obs = y[,1])$H)
    return(H)
}


#' @title Extract responses y
#'
#' @description Give a list of observations, extract responses y
#'
#' @param basis Basis object.
#' @param X Observations
#' @param coverage_ind Which observations have coverage
#'
#' @return The design matrix H
#'
extract_y <- function(X, coverage_ind){
    # Responses y for each genomic region -> i.e. gene
    y <- vector(mode = "list", length = length(X))
    # Set all values to NA
    y <- lapply(y, function(x) x <- NA)
    if (length(coverage_ind) < 1) {
        stop("No coverage across all regions. This cell should be discarded!")
    }
    y[coverage_ind] <- lapply(X = X[coverage_ind], FUN = function(x) x[,2])
    return(y)
}


#'
#' Generate synthetic data for the Melissa model
#'
#' @param basis Basis function object
#' @param N Number of cells
#' @param M Number of genomic regions per cell
#' @param K Number of clusters
#' @param C_true Ground truth cluster labels
#' @param max_cov Maximum CpG coverage
#' @param cluster_var Cell variability across clusters
generate_synth_data <- function(basis, N = 50, M = 300, K = 2, pi_k = rep(1/K, K),
                                C_true = NULL, max_cov=25, cluster_var=0.5){
    D <- basis$M + 1   # Number of covariates
    X <- vector(mode = "list", length = N) # Keep cells
    w_nk <- array(0, dim = c(M, D, K))

    # Generate overall clusterings C_n from mixing proportions
    if (is.null(C_true)) { C_true <- t(rmultinom(N, 1, pi_k))
    }else if (NROW(C_true) != N || NCOL(C_true) != K) {
        stop("C_true matrix dimensions do not match!")
    }
    # Which cells belong in cluster k
    idx <- list()
    # Extract cells that belong to the same cluster
    for (k in 1:K) { idx[[k]] <- which(C_true[, k] == 1) }
    # Iterate over each genomic region
    for (m in 1:M) {
        # Sample from Bernoulli to either consider this region similar or not
        # across clusters based on the variability across cells
        is_variable <- rbinom(n = 1, size = 1, prob = cluster_var)
        # Generate random weights
        if (!is_variable) { w_sample <- rnorm(D, mean = 0, sd = 1.2) }
        # Iterate over each cluster
        for (k in 1:K) {
            # If we have similar profile across cells for all clusters, just
            # add a small noise in profiles
            if (!is_variable) {
                w_nk[m, , k] <- w_sample + rnorm(D, mean = 0, sd = 0.05)
                # Otherwise we create a new profile for each cluster
            } else{w_nk[m, , k] <- rnorm(D, mean = 0, sd = 1.2) }
            # Simulate data for cells belonging to cluster k
            for (n in idx[[k]]) {
                X[[n]][[m]] <- generate_bpr_data(basis = basis,w = w_nk[m,,k],
                                                 max_cov = max_cov)
            } }
    }
    return(list(X = X, W = w_nk, C_true = C_true))
}


#'
#' Partition synthetic dataset to training and test set
#'
#' @param X List of all observations.
#' @param data_train_prcg Percentage of data that will be fully used for
#'   training, without any CpGs missing.
#' @param region_train_prcg Fraction of promoters to keep for training set, i.e.
#'   some promoter regions will have no coverage.
#' @param cpg_train_prcg Fraction of CpGs in each promoter region to keep for
#'   training set.
#' @param is_synth Logical, whether we have synthetic data or not.
#'
partition_dataset <- function(X, data_train_prcg = 0.5, region_train_prcg = 0.8,
                              cpg_train_prcg = 0.5, is_synth = FALSE){
    train = test <- X
    N <- length(X)       # Number of cells
    M <- length(X[[1]])  # Number of genomic regions
    for (n in 1:N) {     # Iterate over each cell
        if (!is_synth) { # If we have real data
            # Keep indices of genomic regions that have CpG coverage
            cov_gen_ind <- which(!is.na(X[[n]]))
            # Compute the number of those promoters
            N_cov <- length(cov_gen_ind)
            if (N_cov < 10) { message("Low genomic coverage..."); return(1) }
            pivot <- region_train_prcg * N_cov
            # These are the training data
            train_ind <- cov_gen_ind[sort(sample(N_cov, round(pivot)))]
        }else {
            # Genomic regions with CpG coverage
            pivot <- region_train_prcg * M
            train_ind <- sort(sample(M, round(pivot)))
        }
        train[[n]][-train_ind] <- NA
        test[[n]][train_ind] <- NA
        # Iterate over each covered genomic region
        for (m in train_ind) {
            # Will this sample be used fully for training
            is_train <- rbinom(1, 1, data_train_prcg)
            if (!is_train) {
                # Get fraction of CpGs covered
                pivot <- cpg_train_prcg * NROW(X[[n]][[m]])
                idx <- sort(sample(NROW(X[[n]][[m]]), round(pivot)))
                train[[n]][[m]] <- X[[n]][[m]][idx,,drop = FALSE]
                test[[n]][[m]] <- X[[n]][[m]][-idx,,drop = FALSE]
            }
        }
    }
    return(list(train = train, test = test))
}

#' @title Evaluate imputation model performance
#'
#' @description Evaluate model performance for imputating missing methylation
#'   states.
#'
#' @param obj Output of the BPR-FDMM (joint) or BPR MLE (separate) model
#' @param test Test data to evaluate performance.
#' @param basis Basis object, if NULL we perform joint imputation otherwise
#'   separate
#' @param is_predictive Logical, use predictive distribution for imputation, or
#'   choose the cluster label with the highest responsibility.
#' @param return_test Whether or not to return predictions list
#'
eval_performance <- function(obj, test, basis = NULL, is_predictive = TRUE,
                             return_test = FALSE){
    N         <- length(test)                      # Number of cells
    M         <- length(test[[1]])                 # Number of genomic regions
    test_pred <- test                              # Copy test data
    act_obs   <- vector(mode = "list", length = N) # Keep actual CpG states
    pred_obs  <- vector(mode = "list", length = N) # Keep predicted CpG states
    # Iterate over each cell...
    for (n in 1:N) {
        # Regions that have CpG observations for testing
        idx <- which(!is.na(test[[n]]))
        cell_act_obs  <- vector("list", length = length(idx))
        cell_pred_obs <- vector("list", length = length(idx))
        iter <- 1
        # Iterate over each genomic region
        for (m in idx) {
            # When joint imputation method
            if (is.null(basis)) {
                K <- NCOL(obj$r_nk) # Number of clusters
                # If we use the predictive density for imputation
                if (is_predictive) {
                    tmp_mixt <- vector("numeric", length = length(test[[n]][[m]][, 1]))
                    for (k in 1:K) {
                        # Evalute profile from weighted predictive
                        tmp_mixt <- tmp_mixt + obj$r_nk[n,k] *
                            eval_probit_function(obj = obj$basis,
                                                 obs = test[[n]][[m]][, 1],
                                                 w = obj$W[m,,k])
                    }
                    # Evaluate the methylation profile
                    test_pred[[n]][[m]][,2] <- tmp_mixt
                }else{
                    # Get cluster assignment
                    k <- which.max(obj$r_nk[n,])
                    test_pred[[n]][[m]][,2] <- eval_probit_function(obj = obj$basis,
                            obs = test[[n]][[m]][, 1], w = obj$W[m,,k])
                }
            }else{
                test_pred[[n]][[m]][,2] <- eval_probit_function(obj = basis,
                            obs = test[[n]][[m]][, 1], w = obj[m,,n])
            }
            # TODO
            # Actual CpG states
            cell_act_obs[[iter]]  <- test[[n]][[m]][, 2]
            # Predicted CpG states (i.e. function evaluations)
            cell_pred_obs[[iter]] <- test_pred[[n]][[m]][, 2]
            iter <- iter + 1
        }
        # Combine data to a big vector
        act_obs[[n]] <- do.call("c", cell_act_obs)
        pred_obs[[n]] <- do.call("c", cell_pred_obs)
    }
    if (return_test) {
        return(list(test_pred = test_pred, act_obs = do.call("c", act_obs),
                    pred_obs = do.call("c", pred_obs)))
    }else{
        return(list(act_obs = do.call("c", act_obs),
                    pred_obs = do.call("c", pred_obs)))
    }
}


#' Compute clustering error
#'
#' \code{cluster_error} computes the clustering assignment error, i.e. the
#' average number of incorrect cluster assignments: \deqn{OE =
#' \sum_{n=1}^{N}(I(LT_{n} \neq LP_{n})) / N}
#'
#' @param C_true True cluster assignemnts.
#' @param C_post Posterior mean of predicted cluster assignemnts.
#'
#' @return The clustering assignment error
#'
cluster_error <- function(C_true, C_post){
    # Obtain the total number of objects
    N <- NROW(C_post)
    # Align cluster indices
    C_post <- align_cluster(Z1 = C_true, Z2 = C_post, type = "mat")$Z2
    # Number of correct assignments
    C_match <- sum((C_post == C_true))
    # Compute the error
    error <- 1 - C_match / (NCOL(C_post) * N)
    return(error)
}


#' @title Compute clustering ARI
#'
#' @description \code{cluster_ari} computes the clustering Adjusted Rand Index.
#'
#' @param C_true True cluster assignemnts.
#' @param C_post Posterior responsibilities of predicted cluster assignemnts.
#'
#' @return The clustering ARI.
#'
cluster_ari <- function(C_true, C_post){
    # Obtain labels from 1-hot-encoding data
    C_true_lab <- unlist(apply(C_true, 1,
                               function(x) which(x == max(x, na.rm = TRUE))[1]))
    C_est_lab <- unlist(apply(C_post, 1,
                              function(x) which(x == max(x, na.rm = TRUE))[1]))
    # Compute the overall clustering ARI
    ari <- clues::adjustedRand(C_true_lab, C_est_lab, randMethod = "HA")
    return(ari)
}


#' @title Align cluster labels
#'
#' @description Align cluster labels when we have (soft) clusterings, i.e.
#'   responsibilities from mixture models.
#'
#' @param Z1 True cluster assignments.
#' @param Z2 Estimate cluster assignments
#' @param type Object type of the cluster assignemnts, either 'array',
#'   'mat_array', 'array_mat', 'mat' or 'vec'.
#'
#' @return The aligned labels.
#'
align_cluster <- function(Z1, Z2, params = NULL, type = "mat"){
    if (type == "array") {
        L1 <- apply(Z1, c(1,3), sum) # Source 1 cluster assignments
        Lm <- apply(Z2, c(1,3), sum) # Source m cluster assignments
    }else if (type == "mat_array") {
        L1 <- Z1
        Lm <- apply(Z2, c(1,3), sum) # Source m cluster assignments
    }else if (type == "array_mat") {
        L1 <- apply(Z1, c(1,3), sum) # Source 1 cluster assignments
        Lm <- Z2
    }else if (type == "mat") {
        L1 <- Z1
        Lm <- Z2
    }else if (type == "vec") {
        # For each cluster k in previous Cluster
        for (k in 1:length(unique(Z1))) {
            # Find Max
            Max <- sum(Z1 == k & Z2 == k)/(.01 + sum(Z1 == k) + sum(Z2 == k))
            # For each cluster k in current Cluster
            for (tempk in 1:length(unique(Z2))) {
                # Check if the proportions are higher than Max
                if ( (sum(Z1 == k & Z2 == tempk)/(.01 + sum(Z1 == k) +
                                                  sum(Z2 == tempk))) > Max) {
                    # Get proportion that the two cluster indices are the same
                    Max <- sum(Z1 == k & Z2 == tempk)/(.01 + sum(Z1 == k) +
                                                           sum(Z2 == tempk))
                    dummy <- (Z2 == k)    # Keep indices that do not match
                    Z2[Z2 == tempk] <- k  # Swap the incorrect indices
                    Z2[dummy] <- tempk    # Swap the incorrect indices
                }
            }
        }
        return(list(Z2 = Z2, params = params))
    }else {stop("Type of alignment not implemented yet.") }

    ##-------------------------------------
    # When input is any type except "vec"
    K <- NCOL(Z1) # Number of clusters
    N <- NROW(Z1) # Number of objects
    # TODO: Check if this holds in the general case
    # Convert probs to hard clusterings
    L1_k <- max.col(L1)
    Lm_k <- max.col(Lm)

    L1 <- matrix(0, ncol = K, nrow = N)
    Lm <- matrix(0, ncol = K, nrow = N)
    for (k in 1:K) {
        cl_L1 <- which(L1_k == k)
        cl_Lm <- which(Lm_k == k)
        if (length(cl_L1) > 0) { L1[cl_L1, k] <- 1 }
        if (length(cl_Lm) > 0) { Lm[cl_Lm, k] <- 1 }
    }

    for (k in 1:NCOL(L1)) {            # For each cluster k in L1 source
        for (tempk in 1:NCOL(Lm)) {    # For each cluster k in Lm source
            Max <- sum(L1 == Lm)       # Matches between the cluster indices
            Lm_dummy <- Lm             # Keep the Lm indices in a dummy variable
            Lm_dummy[,k] = Lm[,tempk]  # Swap the incorrect indices
            Lm_dummy[,tempk] = Lm[,k]  # Swap the incorrect indices
            # If the swaps make a better alignment, update indices
            if (sum(L1 == Lm_dummy) > Max) {
                Lm <- Lm_dummy
                if (type == "array" || type == "mat_array") {
                    tmp <- Z2
                    Z2[,,k] <- tmp[,,tempk]
                    Z2[,,tempk] <- tmp[,,k]
                }
                # Align model parameters using the new cluster indices
                if (!is.null(params)) {
                    params <- align_params(obj = params, curr_label = k,
                                           new_label = tempk)
                }
            }
        }
    }
    if (type == "mat") { Z2 <- Lm}
    return(list(Z2 = Z2, params = params)) # Return the aligned parameters
}



#' @title Filter regions across cells
#'
#' @description Filter genomic regions that have low coverage across many cells.
#'   So we keep regions from which we can share information across cells.
#'
#' @param dt Input data object
#' @param opts Options, having at least N = cells, M = regions, filt_region_cov
#'   = threshold on % cells that have coverage for each region.
#'
#' @return The filtered methylation data
#'
filter_regions_across_cells <- function(dt, opts) {
    ##----------------------------------------------------------------------
    message("Filtering low covered regions across cells...")
    ##----------------------------------------------------------------------
    met <- dt$met
    non_cov_reg <- vector(mode = "numeric")
    for (m in 1:opts$M) { # Iterate over each region
        # Number of cells covered in each source
        cov_cells <- length(which(!is.na(lapply(met, "[[", m))))
        # If no coverage
        if (length(cov_cells) == 0) {
            non_cov_reg <- c(non_cov_reg, m)
        }else{
            # If coverage does not pass threshold, filter again
            if (cov_cells/opts$N < opts$filt_region_cov) {
                non_cov_reg <- c(non_cov_reg, m)
            }
        }
    }
    dt$anno_region <- dt$anno_region[-non_cov_reg]  # Filter anno regions
    dt$annos       <- dt$annos[-non_cov_reg]        # Filter annotation data
    for (n in 1:opts$N) {                           # Filter met data
        met[[n]] <- met[[n]][-non_cov_reg]
        met[[n]] <- unname(met[[n]])
    }
    dt$met <- met
    return(dt)
}


#' @title Melissa model imputation analysis
#'
#'
melissa_imputation_analysis <- function(X, opts){
    # Partition to training and test sets
    dt <- partition_dataset(X = X, data_train_prcg = opts$data_train_prcg,
                            region_train_prcg = opts$region_train_prcg,
                            cpg_train_prcg = opts$cpg_train_prcg, is_synth = FALSE)
    rm(X)
    # Run CC EM algorithm
    print("Starting Melissa model - profile")
    scvb_prof <- tryCatch({
        melissa_obj <- melissa_vb(X = dt$train, K = opts$K, basis = opts$basis_prof,
                                  delta_0 = opts$delta_0, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                                  vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                                  is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                                  vb_init_max_iter = opts$vb_init_max_iter, is_parallel = opts$is_parallel,
                                  no_cores = opts$no_cores, is_verbose = TRUE)
    }, error = function(err){
        # error handler picks up where error was generated
        print(paste("ERROR:  ", err)); return(NULL)
    }) # END tryCatch

    # Using methylation profiles
    eval_prof <- eval_performance(obj = scvb_prof, test = dt$test)
    ##----------------------------------------------------------------------
    message("Computing AUC...")
    ##----------------------------------------------------------------------
    pred_prof <- prediction(eval_prof$pred_obs, eval_prof$act_obs)
    auc_prof <- performance(pred_prof, "auc")
    message(unlist(auc_prof@y.values))

    print("Starting Melissa model - rate")
    scvb_mean <- tryCatch({
        melissa_obj <- melissa_vb(X = dt$train, K = opts$K, basis = opts$basis_mean,
                                  delta_0 = opts$delta_0, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                                  vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                                  is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                                  vb_init_max_iter = opts$vb_init_max_iter, is_parallel = opts$is_parallel,
                                  no_cores = opts$no_cores, is_verbose = TRUE)
    }, error = function(err){
        # error handler picks up where error was generated
        print(paste("ERROR:  ", err)); return(NULL)
    }) # END tryCatch

    ##----------------------------------------------------------------------
    message("Computing AUC...")
    ##----------------------------------------------------------------------
    # Using methylation rate
    eval_mean <- eval_performance(obj = scvb_mean, test = dt$test)
    pred_mean <- prediction(eval_mean$pred_obs, eval_mean$act_obs)
    auc_mean <- performance(pred_mean, "auc")
    message(unlist(auc_mean@y.values))

    eval_perf <- list(eval_prof = eval_prof, eval_mean = eval_mean)
    obj <- list(melissa_prof = scvb_prof, melissa_rate = scvb_mean,
                eval_perf = eval_perf, opts = opts)
    return(obj)
}


#' @title Independent model imputation analysis
#'
#'
indep_imputation_analysis <- function(X, opts){
    # Partition to training and test sets
    dt <- partition_dataset(X = X, data_train_prcg = opts$data_train_prcg,
                            region_train_prcg = opts$region_train_prcg,
                            cpg_train_prcg = opts$cpg_train_prcg, is_synth = FALSE)
    rm(X)

    # List of genes with no coverage for each cell
    region_ind <- lapply(X = 1:opts$N, FUN = function(n) which(!is.na(dt$train[[n]])))
    na_ind <- lapply(X = 1:opts$N, FUN = function(n) which(is.na(dt$train[[n]])))
    # List of cells with no coverage for each genomic region
    cell_ind <- lapply(X = 1:opts$M, FUN = function(m) which(!is.na(lapply(dt$train, "[[", m))))
    # Infer MLE methylation profiles for each cell and region
    if (opts$is_parallel) { W_vb_prof <- parallel::mclapply(X = 1:opts$N, FUN = function(n)
        infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                           basis = opts$basis_prof, lambda = 1/10, opt_itnmax = 50)$W,
        mc.cores = opts$no_cores)
    }else{W_vb_prof <- lapply(X = 1:opts$N, FUN = function(n)
        infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                           basis = opts$basis_prof, lambda = 1/10, opt_itnmax = 50)$W)
    }
    # Infer MLE methylation rate for each cell and region
    if (opts$is_parallel) { W_vb_mean <- parallel::mclapply(X = 1:opts$N, FUN = function(n)
        infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                           basis = opts$basis_mean, lambda = 1/10, opt_itnmax = 50)$W,
        mc.cores = opts$no_cores)
    }else{W_vb_mean <- lapply(X = 1:opts$N, FUN = function(n)
        infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                           basis = opts$basis_mean, lambda = 1/10, opt_itnmax = 50)$W)
    }
    # Store data in array objects
    W_prof <- array(data = 0, dim = c(opts$M, opts$basis_prof$M + 1, opts$N))
    W_mean <- array(data = 0, dim = c(opts$M, 1, opts$N))
    for (n in 1:opts$N) { # Iterate over each cell
        # Store optimized W to an array object (genes x basis x cells)
        W_prof[region_ind[[n]],,n] <- W_vb_prof[[n]]
        W_mean[region_ind[[n]],,n] <- W_vb_mean[[n]]
    }
    # Cases when we have empty genomic regions, sample randomly from other cell
    # methylation profiles
    for (cell in 1:opts$N) {
        if (length(na_ind[[cell]]) > 0) {
            for (m in na_ind[[cell]]) {
                ind_cell <- sample(cell_ind[[m]], 1)
                W_prof[m, , cell] <- W_prof[m, , ind_cell]
                W_mean[m, , cell] <- W_mean[m, , ind_cell]
            }
        }
    }
    indep_prof <- W_prof
    indep_mean <- W_mean
    # Evalute model performance
    eval_prof <- eval_performance(obj = W_prof, test = dt$test, basis = opts$basis_prof)
    eval_mean <- eval_performance(obj = W_mean, test = dt$test, basis = opts$basis_mean)
    eval_perf <- list(eval_prof = eval_prof, eval_mean = eval_mean)
    ##----------------------------------------------------------------------
    message("Computing AUC...")
    ##----------------------------------------------------------------------
    pred_prof <- prediction(eval_prof$pred_obs, eval_prof$act_obs)
    auc_prof <- performance(pred_prof, "auc")
    message(unlist(auc_prof@y.values))
    pred_mean <- prediction(eval_mean$pred_obs, eval_mean$act_obs)
    auc_mean <- performance(pred_mean, "auc")
    message(unlist(auc_mean@y.values))

    obj <- list(indep_prof = indep_prof, indep_mean = indep_mean,
                eval_perf = eval_perf, opts = opts)
    return(obj)
}


#' @title RF model imputation analysis
#'
#'
rf_indep_imputation_analysis <- function(X, opts){
    # Partition to training and test sets
    dt <- partition_dataset(X = X, data_train_prcg = opts$data_train_prcg,
                            region_train_prcg = opts$region_train_prcg,
                            cpg_train_prcg = opts$cpg_train_prcg, is_synth = FALSE)

    # List of genes with no coverage for each cell
    train_region_ind <- lapply(X = 1:opts$N, FUN = function(n) which(!is.na(dt$train[[n]])))
    test_region_ind <- lapply(X = 1:opts$N, FUN = function(n) which(!is.na(dt$test[[n]])))
    # List of cells with no coverage for each genomic region
    cell_ind <- lapply(X = 1:opts$M, FUN = function(m) which(!is.na(lapply(dt$train, "[[", m))))
    # Use RF for prediction

    act_obs = pred_obs <- vector("numeric")
    for (n in 1:opts$N) { # Iterate over the cells
        for (m in test_region_ind[[n]]) { # Iterate over genomic regions
            if (m %in% train_region_ind[[n]]) {
                y <- as.factor(dt$train[[n]][[m]][,2])
                if (length(levels(y)) == 1) {
                    if (as.numeric(levels(y)) == 1) {
                        dt$train[[n]][[m]] <- rbind(dt$train[[n]][[m]], c(0.1, 0))
                    }else{
                        dt$train[[n]][[m]] <- rbind(dt$train[[n]][[m]], c(0.1, 1))
                    }
                }
                model <- randomForest::randomForest(x = dt$train[[n]][[m]][,1, drop = FALSE],
                                                    y = as.factor(dt$train[[n]][[m]][,2]),
                                                    ntree = 60, nodesize = 2)
                pred_obs <- c(pred_obs, predict(object = model, newdata = dt$test[[n]][[m]][,1, drop = FALSE],
                                                type = "prob")[,2])
                act_obs <- c(act_obs, dt$test[[n]][[m]][,2])
            } else{
                # Randomly sample a different cell that has coverage and predict from its profile
                ind_cell <- sample(cell_ind[[m]], 1)
                y <- as.factor(dt$train[[ind_cell]][[m]][,2])
                if (length(levels(y)) == 1) {
                    if (as.numeric(levels(y)) == 1) {
                        dt$train[[ind_cell]][[m]] <- rbind(dt$train[[ind_cell]][[m]], c(0.1, 0))
                    }else{
                        dt$train[[ind_cell]][[m]] <- rbind(dt$train[[ind_cell]][[m]], c(0.1, 1))
                    }
                }
                model <- randomForest::randomForest(x = dt$train[[ind_cell]][[m]][,1, drop = FALSE],
                                                    y = as.factor(dt$train[[ind_cell]][[m]][,2]),
                                                    ntree = 60, nodesize = 2)
                pred_obs <- c(pred_obs, predict(object = model, newdata = dt$test[[n]][[m]][,1, drop = FALSE],
                                                type = "prob")[,2])
                act_obs <- c(act_obs, dt$test[[n]][[m]][,2])
            }
        }
    }

    # Store evaluated performance
    eval_perf <- list(act_obs = act_obs, pred_obs = pred_obs)
    ##----------------------------------------------------------------------
    message("Computing AUC...")
    ##----------------------------------------------------------------------
    pred_rf <- prediction(pred_obs, act_obs)
    auc_rf <- performance(pred_rf, "auc")
    message(unlist(auc_rf@y.values))

    obj <- list(eval_perf = eval_perf, opts = opts)
    return(obj)
}


#' @title DeepCpG model imputation analysis performance
deepcpg_imputation_analysis <- function(X, opts) {
    # Partition to training and test sets
    dt_split <- partition_dataset(X = X, data_train_prcg = opts$data_train_prcg,
                                  region_train_prcg = opts$region_train_prcg,
                                  cpg_train_prcg = opts$cpg_train_prcg, is_synth = FALSE)

    # List of genes with no coverage for each cell
    test_region_ind <- lapply(X = 1:opts$N, FUN = function(n) which(!is.na(dt_split$test[[n]])))
    act_obs = pred_obs <- vector("numeric")

    for (n in 1:opts$N) { # Iterate over the cells
        for (m in test_region_ind[[n]]) { # Iterate over genomic regions
            act_obs <- c(act_obs, dt_split$test[[n]][[m]][,3])
            pred_obs <- c(pred_obs, dt_split$test[[n]][[m]][,2])
        }
    }

    ##----------------------------------------------------------------------
    message("Computing AUC...")
    ##----------------------------------------------------------------------
    pred_rf <- prediction(pred_obs, act_obs)
    auc_rf <- performance(pred_rf, "auc")
    message(unlist(auc_rf@y.values))

    # Store evaluated performance
    eval_perf <- list(act_obs = act_obs, pred_obs = pred_obs)
    obj <- list(eval_perf = eval_perf, opts = opts)
    return(obj)
}
