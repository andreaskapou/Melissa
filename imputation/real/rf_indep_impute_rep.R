# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(ROCR))

rf_indep_imputation_analysis <- function(X, opts){
    # Partition to training and test sets
    dt <- partition_dataset(X = X, region_train_prcg = opts$region_train_prcg,
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
                                                    ntree = 50, nodesize = 2)
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
                                                    ntree = 50, nodesize = 2)
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
    # roc_prof <- performance(pred_prof, "tpr", "fpr")
    auc_rf <- performance(pred_rf, "auc")
    auc_rf <- unlist(auc_rf@y.values)
    message(auc_rf)

    obj <- list(eval_perf = eval_perf, opts = opts)
    return(obj)
}


##------------------------------------
# Load preprocessed data
##------------------------------------
io <- list(data_file = "prom3k", cov = 10, sd = 0.2)
io$data_dir = "../../local-data/melissa/met/"
io$out_dir = "../../local-data/melissa/real/imputation/"
dt <- readRDS(paste0(io$data_dir, "filtered_met/", io$data_file, "_cov",
                     io$cov, "_sd", io$sd, "_gene_var5.rds"))
met <- dt$met
##------------------------------------
# Initialize parameters
##------------------------------------
opts                  <- dt$opts
opts$N                <- length(met) # Number of cells
opts$M                <- length(met[[1]]) # Number of genomic regions
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.5         # % of CpGs kept for training in each region
opts$is_parallel      <- TRUE        # Use parallelized version
opts$no_cores         <- 6           # Number of cores
opts$total_sims       <- 10          # Number of simulations

##----------------------------------------------------------------------
message("Filtering low covered regions across cells...")
##----------------------------------------------------------------------
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
dt$rna         <- dt$rna[-non_cov_reg,]         # Filter RNA data
for (n in 1:opts$N) {                           # Filter met data
    met[[n]] <- met[[n]][-non_cov_reg]
    met[[n]] <- unname(met[[n]])
}
opts$M         <- NROW(dt$annos)                # Total number of regions
rm(dt)

# Run model
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = opts$total_sims,
                                          is_parallel = TRUE,
                                          max_cores = opts$total_sims)
print(date())
model <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
    rf_indep_imputation_analysis(X = met, opts = opts), mc.cores = no_cores_out)
print(date())


##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
obj <- list(model = model, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "final_rf_indep_sim", opts$total_sims,
                           "_", io$data_file,
                           "_cov", io$cov,
                           "_sd", io$sd,
                           "_M", opts$M,
                           "_dataPrcg", opts$data_train_prcg,
                           "_regionPrcg", opts$region_train_prcg,
                           "_cpgPrcg", opts$cpg_train_prcg,
                           "_filter", opts$filt_region_cov,
                           "_gene_var5.rds") )
##----------------------------------------------------------------------
message("Done ...")
##----------------------------------------------------------------------
