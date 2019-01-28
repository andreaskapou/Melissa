# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ROCR))

rf_indep_cpg_analysis <- function(opts, sim){
    # Initialize lists
    eval_perf <- vector("list", length = length(opts$cpg_train_prcg))
    # Load synthetic data
    io <- list(data_file = paste0("high_noise_encode_data_", sim, ".rds"),
               data_dir = "../../local-data/melissa/synthetic/imputation/coverage/raw/data-sims/")
    obj <- readRDS(paste0(io$data_dir, io$data_file))
    i <- 1
    # Iterate
    for (cpg_iter in opts$cpg_train_prcg) {
        # Partition to training and test sets
        dt <- partition_dataset(X = obj$synth_data$X, region_train_prcg = opts$region_train_prcg,
                                cpg_train_prcg = cpg_iter, is_synth = TRUE)

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
                    if (length(y) > 0) {
                        if (length(levels(y)) == 1) {
                            if (as.numeric(levels(y)) == 1) {
                                dt$train[[n]][[m]] <- rbind(dt$train[[n]][[m]], c(0.1, 0))
                            }else{
                                dt$train[[n]][[m]] <- rbind(dt$train[[n]][[m]], c(0.1, 1))
                            }
                        }
                        model <- randomForest::randomForest(x = dt$train[[n]][[m]][,1, drop = FALSE],
                                                            y = as.factor(dt$train[[n]][[m]][,2]),
                                                            ntree = 50, nodesize = 2)
                        pred_obs <- c(pred_obs, predict(object = model, newdata = dt$test[[n]][[m]][,1, drop = FALSE],
                                                        type = "prob")[,2])
                        act_obs <- c(act_obs, dt$test[[n]][[m]][,2])
                    }
                } else{
                    # Randomly sample a different cell that has coverage and predict from its profile
                    ind_cell <- sample(cell_ind[[m]], 1)
                    y <- as.factor(dt$train[[ind_cell]][[m]][,2])
                    if (length(y) > 0) {
                        if (length(levels(y)) == 1) {
                            if (as.numeric(levels(y)) == 1) {
                                dt$train[[ind_cell]][[m]] <- rbind(dt$train[[ind_cell]][[m]], c(0.1, 0))
                            }else{
                                dt$train[[ind_cell]][[m]] <- rbind(dt$train[[ind_cell]][[m]], c(0.1, 1))
                            }
                        }
                        model <- randomForest::randomForest(x = dt$train[[ind_cell]][[m]][,1, drop = FALSE],
                                                            y = as.factor(dt$train[[ind_cell]][[m]][,2]),
                                                            ntree = 50, nodesize = 2)
                        pred_obs <- c(pred_obs, predict(object = model, newdata = dt$test[[n]][[m]][,1, drop = FALSE],
                                                       type = "prob")[,2])
                        act_obs <- c(act_obs, dt$test[[n]][[m]][,2])
                    }
                }
            }
        }

        # Store evaluated performance
        eval_perf[[i]] <- list(act_obs = act_obs, pred_obs = pred_obs)
        ##----------------------------------------------------------------------
        message("Computing AUC...")
        ##----------------------------------------------------------------------
        pred_rf <- prediction(pred_obs, act_obs)
        # roc_prof <- performance(pred_prof, "tpr", "fpr")
        auc_rf <- performance(pred_rf, "auc")
        auc_rf <- unlist(auc_rf@y.values)
        message(auc_rf)
        i <- i + 1 # Increase counter
    }
    obj <- list(eval_perf = eval_perf, opts = opts)
    return(obj)
}


##------------------------
# Load synthetic data
##------------------------
io <- list(data_file = paste0("raw/data-sims/high_noise_encode_data_1.rds"),
           out_dir = "../../local-data/melissa/synthetic/imputation/coverage/")
obj <- readRDS(paste0(io$out_dir, io$data_file))
opts                   <- obj$opts      # Get options
opts$data_train_prcg   <- 0.1           # % of data to keep fully for training
opts$region_train_prcg <- 1             # % of regions kept for training
opts$cpg_train_prcg    <- seq(.1,.9,.1) # % of CpGs kept for training
opts$is_parallel       <- TRUE          # Use parallelized version
opts$no_cores          <- 3             # Number of cores
rm(obj)

# Parallel analysis
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = opts$total_sims,
                                          is_parallel = TRUE)
print(date())
obj <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
    rf_indep_cpg_analysis(opts = opts, sim = sim), mc.cores = no_cores_out)
print(date())

##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
saveRDS(obj, file = paste0(io$out_dir, "high_noise_encode_rf_indep_K", opts$K,
                           "_rbf", opts$basis_prof$M,
                           "_dataTrain", opts$data_train_prcg,
                           "_regionTrain", opts$region_train_prcg,
                           "_clusterVar", opts$cluster_var, ".rds") )
