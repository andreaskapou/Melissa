# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ROCR))
set.seed(123)

melissa_cell_analysis <- function(opts, sim){
  # Initialize lists
  melissa_prof = melissa_mean = eval_perf <- vector("list", length = length(opts$N))
  i <- 1
  # Iterate
  for (m in opts$N) {
    # Load synthetic data
    io <- list(data_file = paste0("encode_data_", m, "_", sim, ".rds"),
               data_dir = "../../local-data/melissa/synthetic/imputation/cells/raw/data-sims/")
    obj <- readRDS(paste0(io$data_dir, io$data_file))
    # Partition to training and test sets
    dt <- partition_dataset(X = obj$synth_data$X, region_train_prcg = opts$region_train_prcg,
                            cpg_train_prcg = opts$cpg_train_prcg, is_synth = TRUE)
    # Using methylation profiles
    prof_obj <- melissa_vb(X = dt$train, K = opts$K, basis = opts$basis_prof,
                           delta_0 = opts$delta_0, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                           vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                           is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                           vb_init_max_iter = opts$vb_init_max_iter, is_parallel = opts$is_parallel,
                           no_cores = opts$no_cores, is_verbose = TRUE)
    mean_obj <- melissa_vb(X = dt$train, K = opts$K, basis = opts$basis_mean,
                           delta_0 = opts$delta_0, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                           vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                           is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                           vb_init_max_iter = opts$vb_init_max_iter, is_parallel = opts$is_parallel,
                           no_cores = opts$no_cores, is_verbose = TRUE)
    melissa_prof[[i]] <- prof_obj
    melissa_mean[[i]] <- mean_obj

    ##----------------------------------------------------------------------
    #message("Evaluating model performance...")
    ##----------------------------------------------------------------------
    # Using methylation profiles
    eval_prof <- eval_performance(obj = prof_obj, test = dt$test)
    # Using methylation rate
    eval_mean <- eval_performance(obj = mean_obj, test = dt$test)
    eval_perf[[i]] <- list(eval_prof = eval_prof, eval_mean = eval_mean)

    ##----------------------------------------------------------------------
    message("Computing AUC...")
    ##----------------------------------------------------------------------
    pred_prof <- prediction(eval_prof$pred_obs, eval_prof$act_obs)
    auc_prof <- performance(pred_prof, "auc")
    auc_prof <- unlist(auc_prof@y.values)
    message(auc_prof)

    pred_mean <- prediction(eval_mean$pred_obs, eval_mean$act_obs)
    auc_mean <- performance(pred_mean, "auc")
    auc_mean <- unlist(auc_mean@y.values)
    message(auc_mean)
    i <- i + 1 # Increase counter
  }
  obj <- list(melissa_prof = melissa_prof, melissa_rate = melissa_mean,
              eval_perf = eval_perf, opts = opts)
  return(obj)
}

##------------------------
# Load synthetic data
##------------------------
io <- list(data_file = paste0("raw/data-sims/encode_data_25_1.rds"),
           out_dir = "../../local-data/melissa/synthetic/imputation/cells/")
obj <- readRDS(paste0(io$out_dir, io$data_file))
opts                   <- obj$opts       # Get options
opts$delta_0           <- rep(3, opts$K) # Dirichlet prior
opts$alpha_0           <- 0.5            # Gamma prior
opts$beta_0            <- NULL           # Gamma prior
opts$data_train_prcg   <- 0.1            # % of data to keep fully for training
opts$region_train_prcg <- 1              # % of regions kept for training
opts$cpg_train_prcg    <- 0.4            # % of CpGs kept for training in each region
opts$is_kmeans         <- TRUE           # Use K-means for initialization
opts$vb_max_iter       <- 300            # Maximum VB iterations
opts$epsilon_conv      <- 1e-4           # Convergence threshold for VB
opts$vb_init_nstart    <- 10             # Mini VB restarts
opts$vb_init_max_iter  <- 20             # Mini VB iteratiions
opts$is_parallel       <- TRUE           # Use parallelized version
opts$no_cores          <- 3              # Number of cores
rm(obj)

# Parallel analysis
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = opts$total_sims,
                                          is_parallel = TRUE)
print(date())
obj <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
  melissa_cell_analysis(opts = opts, sim = sim), mc.cores = no_cores_out)
print(date())

##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
saveRDS(obj, file = paste0(io$out_dir, "encode_melissa_K", opts$K,
                           "_rbf", opts$basis_prof$M,
                           "_dataTrain", opts$data_train_prcg,
                           "_regionTrain", opts$region_train_prcg,
                           "_cpgTrain", opts$cpg_train_prcg, ".rds") )
