# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ROCR))

indep_var_analysis <- function(opts, sim){
  # Initialize lists
  indep_prof = indep_mean = eval_perf <- vector("list", length = length(opts$N))
  i <- 1
  # Iterate
  for (cl_iter in opts$N) {
    # Load synthetic data
    io <- list(data_file = paste0("encode_data_", cl_iter, "_", sim, ".rds"),
               data_dir = "../../local-data/melissa/synthetic/imputation/cells/raw/data-sims/")
    obj <- readRDS(paste0(io$data_dir, io$data_file))
    # Partition to training and test sets
    dt <- partition_dataset(X = obj$synth_data$X, region_train_prcg = opts$region_train_prcg,
                            cpg_train_prcg = opts$cpg_train_prcg, is_synth = TRUE)
    # List of genes with no coverage for each cell
    region_ind <- lapply(X = 1:cl_iter, FUN = function(n) which(!is.na(dt$train[[n]])))
    na_ind <- lapply(X = 1:cl_iter, FUN = function(n) which(is.na(dt$train[[n]])))
    # List of cells with no coverage for each genomic region
    cell_ind <- lapply(X = 1:opts$M, FUN = function(m) which(!is.na(lapply(dt$train, "[[", m))))
    # Infer MLE methylation profiles for each cell and region
    if (opts$is_parallel) { W_vb_prof <- parallel::mclapply(X = 1:cl_iter, FUN = function(n)
      infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                         basis = opts$basis_prof, lambda = 1/10, opt_itnmax = 40)$W,
      mc.cores = opts$no_cores)
    }else{W_vb_prof <- lapply(X = 1:cl_iter, FUN = function(n)
      infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                         basis = opts$basis_prof, lambda = 1/10, opt_itnmax = 40)$W)
    }
    # Infer MLE methylation rate for each cell and region
    if (opts$is_parallel) { W_vb_mean <- parallel::mclapply(X = 1:cl_iter, FUN = function(n)
      infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                         basis = opts$basis_mean, lambda = 1/10, opt_itnmax = 40)$W,
      mc.cores = opts$no_cores)
    }else{W_vb_mean <- lapply(X = 1:cl_iter, FUN = function(n)
      infer_profiles_mle(X = dt$train[[n]][region_ind[[n]]], model = "bernoulli",
                         basis = opts$basis_mean, lambda = 1/10, opt_itnmax = 40)$W)
    }
    # Store data in array objects
    W_prof <- array(data = 0, dim = c(opts$M, opts$basis_prof$M + 1, cl_iter))
    W_mean <- array(data = 0, dim = c(opts$M, 1, cl_iter))
    for (n in 1:cl_iter) { # Iterate over each cell
      # Store optimized W to an array object (genes x basis x cells)
      W_prof[region_ind[[n]],,n] <- W_vb_prof[[n]]
      W_mean[region_ind[[n]],,n] <- W_vb_mean[[n]]
    }
    # Cases when we have empty genomic regions, sample randomly from other cell
    # methylation profiles
    for (cell in 1:cl_iter) {
      if (length(na_ind[[cell]]) > 0) {
        for (m in na_ind[[cell]]) {
          ind_cell <- sample(cell_ind[[m]], 1)
          W_prof[m, , cell] <- W_prof[m, , ind_cell]
          W_mean[m, , cell] <- W_mean[m, , ind_cell]
        }
      }
    }
    indep_prof[[i]] <- W_prof
    indep_mean[[i]] <- W_mean
    # Evalute model performance
    eval_prof <- eval_performance(obj = W_prof, test = dt$test, basis = opts$basis_prof)
    eval_mean <- eval_performance(obj = W_mean, test = dt$test, basis = opts$basis_mean)
    eval_perf[[i]] <- list(eval_prof = eval_prof, eval_mean = eval_mean)
    ##----------------------------------------------------------------------
    message("Computing AUC...")
    ##----------------------------------------------------------------------
    pred_prof <- prediction(eval_prof$pred_obs, eval_prof$act_obs)
    # roc_prof <- performance(pred_prof, "tpr", "fpr")
    auc_prof <- performance(pred_prof, "auc")
    auc_prof <- unlist(auc_prof@y.values)
    message(auc_prof)
    pred_mean <- prediction(eval_mean$pred_obs, eval_mean$act_obs)
    # roc_mean <- performance(pred_mean, "tpr", "fpr")
    auc_mean <- performance(pred_mean, "auc")
    auc_mean <- unlist(auc_mean@y.values)
    message(auc_mean)
    i <- i + 1 # Increase counter
  }
  obj <- list(indep_prof = indep_prof, indep_mean = indep_mean,
              eval_perf = eval_perf, opts = opts)
  return(obj)
}



##------------------------
# Load synthetic data
##------------------------
io <- list(data_file = paste0("raw/data-sims/encode_data_25_1.rds"),
           out_dir = "../../local-data/melissa/synthetic/imputation/cells/")
obj <- readRDS(paste0(io$out_dir, io$data_file))
opts                   <- obj$opts # Get options
opts$data_train_prcg   <- 0.1      # % of data to keep fully for training
opts$region_train_prcg <- 1        # % of regions kept for training
opts$cpg_train_prcg    <- 0.4      # % of CpGs kept for training in each region
opts$is_parallel       <- TRUE     # Use parallelized version
opts$no_cores          <- 3        # Number of cores
rm(obj)

# Parallel analysis
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = opts$total_sims,
                                          is_parallel = TRUE)
print(date())
obj <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
  indep_var_analysis(opts = opts, sim = sim), mc.cores = no_cores_out)
print(date())

##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
saveRDS(obj, file = paste0(io$out_dir, "encode_indep_K", opts$K,
                           "_rbf", opts$basis_prof$M,
                           "_dataTrain", opts$data_train_prcg,
                           "_regionTrain", opts$region_train_prcg,
                           "_cpgTrain", opts$cpg_train_prcg, ".rds") )
