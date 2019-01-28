#' @title Compute stable log-sum-exp
#'
#' @description \code{log_sum_exp} computes the log sum exp trick for avoiding
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
#' @description Given a list of observations, initialise design matrices H for
#'   computational efficiency.
#'
#' @param basis Basis object.
#' @param X Observations
#' @param coverage_ind Which observations have coverage
#'
#' @return The design matrix H
#'
#' @importFrom BPRMeth design_matrix
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
#' @description Given a list of observations, extract responses y
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



#' @title Partition synthetic dataset to training and test set
#'
#' @param dt_obj Melissa data object
#' @param data_train_prcg Percentage of genomic regions that will be fully used
#'   for training, i.e. across the whole region we will have no CpGs missing.
#' @param region_train_prcg Fraction of genomic regions to keep for training
#'   set, i.e. some genomic regions will have no coverage at all during
#'   training.
#' @param cpg_train_prcg Fraction of CpGs in each genomic region to keep for
#'   training set.
#' @param is_synth Logical, whether we have synthetic data or not.
#'
#' @return The Melissa object with the following changes. The `met` element will
#'   now contain only the `training` data. An additional element called
#'   `met_test` will store the data that will be used during testing to evaluate
#'   the imputation performance. These data will not be seen from Melissa during
#'   inference.
#'
#' @export
#'
partition_dataset <- function(dt_obj, data_train_prcg = 0.5,
                              region_train_prcg = 0.95, cpg_train_prcg = 0.5,
                              is_synth = FALSE){
  # If `met_test` element already exists, stop
  if (!is.null(dt_obj$met_test)) { stop("Stopping. Test data already exist!")}
  train = test <- dt_obj$met
  N <- length(dt_obj$met)       # Number of cells
  M <- length(dt_obj$met[[1]])  # Number of genomic regions
  for (n in 1:N) {     # Iterate over each cell
    if (!is_synth) { # If we have real data
      # Keep indices of genomic regions that have CpG coverage
      cov_gen_ind <- which(!is.na(dt_obj$met[[n]]))
      # Compute the number of those regions
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
        pivot <- cpg_train_prcg * NROW(dt_obj$met[[n]][[m]])
        idx <- sort(sample(NROW(dt_obj$met[[n]][[m]]), round(pivot)))
        train[[n]][[m]] <- dt_obj$met[[n]][[m]][idx,,drop = FALSE]
        test[[n]][[m]] <- dt_obj$met[[n]][[m]][-idx,,drop = FALSE]
      }
    }
  }
  # Set the training data as the `met` element
  dt_obj$met <- train;
  # Set the test data as the `met_test` element
  dt_obj$met_test <- test
  return(dt_obj)
}

#' @title Impute/predict methylation states
#'
#' @description Make predictions of missing methylation states, i.e. perfrom
#'   imputation using Melissa. This requires keepin a subset of data as a held
#'   out test set during Melissa inference.
#'
#' @param obj Output of Melissa inference object.
#' @param test Test data to evaluate performance.
#' @param basis Basis object, if NULL we perform imputation using Melissa,
#'   otherwise using BPRMeth.
#' @param is_predictive Logical, use predictive distribution for imputation, or
#'   choose the cluster label with the highest responsibility.
#' @param return_test Whether or not to return a list with the predictions.
#'
#' @return A list containing two vectors, the true methylation state and the
#'   predicted/imputed methylation states.
#'
#' @importFrom BPRMeth eval_probit_function
#'
#' @export
impute_met_state <- function(obj, test, basis = NULL, is_predictive = TRUE,
                              return_test = FALSE){
  N         <- length(test)                      # Number of cells
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


#' @title Evaluate imputation performance
#'
#' @description \code{eval_imputation_performance} is a wrapper function for
#'   computing imputation/clustering performance in terms of different metrics,
#'   such as AUC and precision recall curves.
#'
#' @param obj Output of Melissa inference object.
#' @param imputation_obj List containing two vectors of equal length,
#'   corresponding to true methylation states and predicted/imputed methylation
#'   states.
#'
#' @return The `melissa` object, with an additional slot named `imputation`,
#'   containing the AUC, F-measure, True Positive Rate (TPR) and False Positive
#'   Rate (FPR), and Precision Recall (PR) curves.
#'
#' @export
#'
eval_imputation_performance <- function(obj, imputation_obj) {
  # Create object of predictions
  pred_obj <- ROCR::prediction(round(imputation_obj$pred_obs, 2),
                               imputation_obj$act_obs)
  obj$imputation <- list()
  # Store AUC
  obj$imputation$auc <- ROCR::performance(pred_obj, "auc")@y.values[[1]]
  # Store F-measure
  f <- ROCR::performance(pred_obj, "f")
  obj$imputation$f_measure <- f@y.values[[1]][min(which(f@x.values[[1]] <= 0.5))]
  # Store TPR/FPR
  obj$imputation$tpr_fpr <- ROCR::performance(pred_obj, "tpr", "fpr")
  # Store PR
  obj$imputation$pr <- ROCR::performance(pred_obj, "prec", "rec")
  return(obj)
}



#' @title Evaluate clustering performance
#'
#' @description \code{eval_cluster_performance} is a wrapper function for
#'   computing clustering performance in terms of ARI and clustering assignment
#'   error.
#'
#' @param obj Output of Melissa inference object.
#' @param C_true True cluster assignemnts.
#'
#' @return The `melissa` object, with an additional slot named `clustering`,
#'   containing the ARI and clustering assignment error performance.
#'
#' @export
#'
eval_cluster_performance <- function(obj, C_true) {
  # Create new element
  obj$clustering <- list()
  # Compute clustering assignment error
  obj$clustering$error <- cluster_error(C_true, obj$r_nk)
  # Compute ARI for clustering performance
  obj$clustering$ari <- cluster_ari(C_true, obj$r_nk)
  return(obj)
}


#' @title Compute clustering assignment error
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
  C_post <- .align_cluster(Z1 = C_true, Z2 = C_post, type = "mat")$Z2
  # Number of correct assignments
  C_match <- sum((C_post == C_true))
  # Compute the error
  error <- 1 - C_match / (NCOL(C_post) * N)
  return(error)
}


#' @title Compute clustering ARI
#'
#' @description \code{cluster_ari} computes the clustering performance in terms
#'   of the Adjusted Rand Index (ARI) metric.
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


# @title Align cluster labels
#
# @description Align cluster labels when we have (soft) clusterings, i.e.
#   responsibilities from mixture models.
#
# @param Z1 True cluster assignments.
# @param Z2 Estimate cluster assignments
# @param type Object type of the cluster assignemnts, either 'array',
#   'mat_array', 'array_mat', 'mat' or 'vec'.
#
# @return The aligned labels.
#
.align_cluster <- function(Z1, Z2, params = NULL, type = "mat"){
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



#' @title Filter regions by coverage
#'
#' @description Filter genomic regions that have low coverage within that
#'   region. This step does not remove the region, it only sets NA to those
#'   regions.
#'
#' @param obj Melissa data object.
#' @param min_cpgcov Minimum CpG coverage for each genomic region.
#'
#' @return The same object with the regions that do not pass threshold set to
#'   NA.
#'
#' @export
filter_by_cpg_coverage <- function(obj, min_cpgcov = 10) {
  # Consider only regions with enough CpG coverage, the rest are set to NA
  obj$met <- lapply(obj$met, function(x)
    lapply(x, function(y){ if (NROW(y) < min_cpgcov) return(NA) else return(y) }))
  return(obj)
}


#' @title Filter regions by coverage across cells
#'
#' @description Filter genomic regions that have no coverage across many cells.
#'   This way we keep regions from which we can share information across cells.
#'
#' @param obj Melissa data object.
#' @param min_cell_cov_prcg Threshold on the proportion of cells that have
#'   coverage for each region.
#'
#' @return Melissa object without the filtered genomic regions.
#'
#' @export
filter_by_coverage_across_cells <- function(obj, min_cell_cov_prcg = 0.5) {
  N <- length(obj$met)      # Number of cells
  M <- length(obj$met[[1]]) # Number of genomic regions
  non_cov_reg <- vector(mode = "numeric")
  for (m in 1:M) { # Iterate over each region
    # Number of cells covered in each region
    cov_cells <- length(which(!is.na(lapply(obj$met, "[[", m))))
    # If no coverage
    if (length(cov_cells) == 0) {
      non_cov_reg <- c(non_cov_reg, m)
    }else{
      # If coverage does not pass threshold, filter again
      if (cov_cells/N < min_cell_cov_prcg) {
        non_cov_reg <- c(non_cov_reg, m)
      }
    }
  }
  obj$anno_region <- obj$anno_region[-non_cov_reg]  # Filter anno regions
  for (n in 1:N) {                                  # Filter met data
    obj$met[[n]] <- obj$met[[n]][-non_cov_reg]
    obj$met[[n]] <- unname(obj$met[[n]])
  }
  obj$met <- met
  return(obj)
}


#' @title Filter regions by variability
#'
#' @description Filter genomic regions that have low mean methylation variability
#'   across cells; since they are not informative for cell subtype
#'   identification.
#'
#' @param obj Melissa data object.
#' @param min_cov Minimum variability of mean methylation across cells, measured
#'   in terms of standard deviation.
#'
#' @return Melissa object without the filtered genomic regions.
#'
#' @export
filter_by_variability <- function(obj, min_cov = 10) {
  # NUmber of genomic regions
  M <- NROW(obj$anno_region)
  cell_region_sd <- vector("numeric", length = M)
  for (m in 1:M) {
    # Extract all cell observations for specific region m
    tmp <- lapply(obj$met, "[[", m)
    # Keep only regions that have observations
    tmp <- tmp[!is.na(tmp)]
    if (length(tmp) == 0) { cell_region_sd[m] <- 0
    # Compute the standard deviation of region across cells
    } else {cell_region_sd[m] <- sd(sapply(tmp, function(x) mean(x[,2]))) }
  }
  # Keep only highly varying genomic regions
  ind <- which(cell_region_sd > min_cov)

  # Subset according to methylation variability
  obj$anno_region <- obj$anno_region[ind,]
  obj$met <- lapply(obj$met, function(x) x[ind])
  return(obj)
}