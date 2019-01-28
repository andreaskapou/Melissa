#' @title Generate encode synthetic data for the Melissa model
#'
#' @param basis Basis function object
#' @param encode_w Optimized weights of the basis functions for each cluster K
#' @param N Number of cells
#' @param M Number of genomic regions per cell
#' @param K Number of clusters
#' @param C_true Ground truth cluster labels
#' @param max_cov Maximum CpG coverage
#' @param cluster_var Cell variability across clusters
#'
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

  # Parameter options
  opts <- list(C_true = C_true, W = w_nk, basis_obj = basis, N = N, M = M,
               K = K, pi_k = pi_k, max_cov = max_cov, cluster_var = cluster_var,
               cell_names = paste0("cell_", seq(1:N)))
  # Add cell names to list
  names(X) <- opts$cell_names
  # Store the object
  obj <- structure(list(met = X, anno_region = NULL, opts = opts),
                   class = "melissa_data_obj")
  return(obj)
}


#' Generate synthetic data for the Melissa model
#'
#' @param basis Basis function object
#' @param N Number of cells
#' @param M Number of genomic regions per cell
#' @param K Number of clusters
#' @param C_true Ground truth cluster labels
#' @param max_cov Maximum CpG coverage
#' @param cluster_var Cell variability across clusters
#'
generate_synth_data <- function(basis, N = 200, M = 100, K = 2, pi_k = rep(1/K, K),
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

  # Parameter options
  opts <- list(C_true = C_true, W = w_nk, basis_obj = basis, N = N, M = M,
               K = K, pi_k = pi_k, max_cov = max_cov, cluster_var = cluster_var,
               cell_names = paste0("cell_", seq(1:N)))
  # Add cell names to list
  names(X) <- opts$cell_names
  # Store the object
  obj <- structure(list(met = X, anno_region = NULL, opts = opts),
                   class = "melissa_data_obj")
  return(obj)
}


#' @title Generating single-cell methylation data for a given region
#'
#' @param basis Basis function object
#' @param w Weights of basis functions
#' @param max_cov Maximum CpG coverage
#' @param xmin Minimum x location relative to TSS
#' @param xmax Maximum x location relative to TSS
#'
generate_bpr_data <- function(basis, w, max_cov = 25, xmin = -1000, xmax=1000){
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

library(Melissa)
###-----------------------------------
# Generate ENCODE synthetic data
###-----------------------------------
set.seed(1)
# Load ENCODE cluster profile weights
coef_file <- system.file("extdata", "encode-coef.rds", package = "Melissa")
encode_w <- readRDS(coef_file)
basis_prof <- create_rbf_object(M = 4) # Basis function profiles
# Create synthetic data
melissa_encode_dt <- generate_encode_synth_data(basis = basis_prof,
                        encode_w = encode_w, N = 200, M = 100, K = 4,
                        pi_k = c(.2,.4,.15,.25), cluster_var = 0.5)
usethis::use_data(melissa_encode_dt, overwrite = TRUE)


###-----------------------------------
# Generate synthetic data
###-----------------------------------
set.seed(1)
# Basis function profiles
basis_prof <- create_rbf_object(M = 4)
# Create synthetic data
melissa_synth_dt <- generate_synth_data(basis = basis_prof, N = 100, M = 100,
                                        K = 3, pi_k = c(.2,.4,.4))
usethis::use_data(melissa_synth_dt, overwrite = TRUE)
