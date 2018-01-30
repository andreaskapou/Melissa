# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ROCR))

indep_imputation_analysis <- function(X, opts){
    # Partition to training and test sets
    dt <- partition_dataset(X = X, region_train_prcg = opts$region_train_prcg,
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
    # roc_prof <- performance(pred_prof, "tpr", "fpr")
    auc_prof <- performance(pred_prof, "auc")
    auc_prof <- unlist(auc_prof@y.values)
    message(auc_prof)
    pred_mean <- prediction(eval_mean$pred_obs, eval_mean$act_obs)
    # roc_mean <- performance(pred_mean, "tpr", "fpr")
    auc_mean <- performance(pred_mean, "auc")
    auc_mean <- unlist(auc_mean@y.values)
    message(auc_mean)

    obj <- list(indep_prof = indep_prof, indep_mean = indep_mean,
                eval_perf = eval_perf, opts = opts)
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
opts$basis_prof       <- create_rbf_object(M = 7) # Profile basis functions
opts$basis_mean       <- create_rbf_object(M = 0) # Rate basis function

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
    indep_imputation_analysis(X = met, opts = opts), mc.cores = no_cores_out)
print(date())


##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
obj <- list(model = model, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "final_indep_sim", opts$total_sims,
                           "_", io$data_file,
                           "_cov", io$cov,
                           "_sd", io$sd,
                           "_M", opts$M,
                           "_basis", opts$basis_prof$M,
                           "_dataPrcg", opts$data_train_prcg,
                           "_regionPrcg", opts$region_train_prcg,
                           "_cpgPrcg", opts$cpg_train_prcg,
                           "_filter", opts$filt_region_cov,
                           "_gene_var5.rds") )
##----------------------------------------------------------------------
message("Done ...")
##----------------------------------------------------------------------
