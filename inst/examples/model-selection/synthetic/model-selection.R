# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ROCR))
set.seed(123)

joint_dissimilarity_analysis <- function(opts, sim){
    # Initialize lists
    melissa = eval_perf <- vector("list", length = length(opts$cluster_var))
    i <- 1
    # Iterate
    for (m in opts$cluster_var) {
        # Load synthetic data
        io <- list(data_file = paste0("encode_data_", m, "_", sim, ".rds"),
                   data_dir = "../../local-data/melissa/synthetic/imputation/dissimilarity/raw/data-sims/")
        obj <- readRDS(paste0(io$data_dir, io$data_file))
        # Partition to training and test sets
        dt <- partition_dataset(X = obj$synth_data$X, region_train_prcg = opts$region_train_prcg,
                                cpg_train_prcg = opts$cpg_train_prcg, is_synth = TRUE)
        # Using methylation profiles
        bpr_prof_res <- melissa_vb(X = dt$train, K = opts$K, basis = opts$basis_prof,
                                   delta_0 = opts$delta_0, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                                   vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                                   is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                                   vb_init_max_iter = opts$vb_init_max_iter, is_parallel = opts$is_parallel,
                                   no_cores = opts$no_cores, is_verbose = FALSE)
        melissa[[i]] <- bpr_prof_res

        ##----------------------------------------------------------------------
        #message("Evaluating model performance...")
        ##----------------------------------------------------------------------
        # Using methylation profiles
        eval_prof <- eval_performance(obj = bpr_prof_res, test = dt$test)
        eval_perf[[i]] <- list(eval_prof = eval_prof)

        ##----------------------------------------------------------------------
        message("Computing AUC...")
        ##----------------------------------------------------------------------
        pred_prof <- prediction(eval_prof$pred_obs, eval_prof$act_obs)
        # roc_prof <- performance(pred_prof, "tpr", "fpr")
        auc_prof <- performance(pred_prof, "auc")
        auc_prof <- unlist(auc_prof@y.values)
        message(auc_prof)
        i <- i + 1 # Increase counter
    }
    obj <- list(melissa = melissa, eval_perf = eval_perf, opts = opts)
    return(obj)
}


##------------------------
# Load synthetic data
##------------------------
io <- list(data_file = paste0("raw/data-sims/encode_data_0.1_1.rds"),
           data_dir = "../../local-data/melissa/synthetic/imputation/dissimilarity/",
           out_dir = "../../local-data/melissa/synthetic/model-selection/")
obj  <- readRDS(paste0(io$data_dir, io$data_file))
opts <- obj$opts                         # Get options
opts$K                 <- 10             # Number of clusters
opts$data_train_prcg   <- 0.1            # % of data to keep fully for training
opts$region_train_prcg <- 1              # % of regions kept for training
opts$cpg_train_prcg    <- 0.4            # % of CpGs kept for training in each region
opts$vb_max_iter       <- 200            # Maximum VB iterations
opts$epsilon_conv      <- 1e-4           # Convergence threshold for VB
opts$vb_init_nstart    <- 2              # Mini VB restarts
opts$vb_init_max_iter  <- 10             # Mini VB max iterations
opts$is_parallel       <- TRUE           # Use parallelized version
opts$no_cores          <- 2              # Number of cores
opts$total_sims        <- 5              # Total simulations
rm(obj)

# Update prior so the model favours large number of clusters
opts$delta_0           <- rep(3, opts$K) # Dirichlet prior
opts$alpha_0           <- 0.5            # Gamma prior
opts$beta_0            <- 10             # Gamma prior
# # Update prior so the model favours smaller number of clusters
# opts$delta_0           <- rep(1e-2, opts$K) # Dirichlet prior
# opts$alpha_0           <- 2                 # Gamma prior
# opts$beta_0            <- .5                # Gamma prior

# Parallel analysis
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = 1,
                                          is_parallel = TRUE)
print(date())
model <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
    joint_dissimilarity_analysis(opts = opts, sim = sim), mc.cores = no_cores_out)
print(date())


##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
saveRDS(model, file = paste0(io$out_dir, "broad_model_selection_K", opts$K, ".rds") )
