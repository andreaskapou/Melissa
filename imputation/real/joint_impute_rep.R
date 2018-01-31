# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(ROCR))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(require(Matrix))
set.seed(123)

imputation_analysis <- function(X, opts){
    # Partition to training and test sets
    dt <- partition_dataset(X = X, region_train_prcg = opts$region_train_prcg,
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
                                  no_cores = opts$no_cores, is_verbose = FALSE)
    }, error = function(err){
        # error handler picks up where error was generated
        print(paste("ERROR:  ", err)); return(NULL)
    }) # END tryCatch

    print("Starting Melissa model - rate")
    scvb_mean <- tryCatch({
        melissa_obj <- melissa_vb(X = dt$train, K = opts$K, basis = opts$basis_mean,
                                  delta_0 = opts$delta_0, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                                  vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                                  is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                                  vb_init_max_iter = opts$vb_init_max_iter, is_parallel = opts$is_parallel,
                                  no_cores = opts$no_cores, is_verbose = FALSE)
    }, error = function(err){
        # error handler picks up where error was generated
        print(paste("ERROR:  ", err)); return(NULL)
    }) # END tryCatch

    # Using methylation profiles
    eval_prof <- eval_performance(obj = scvb_prof, test = dt$test)
    # Using methylation rate
    eval_mean <- eval_performance(obj = scvb_mean, test = dt$test)
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

    obj <- list(melissa_prof = scvb_prof, melissa_rate = scvb_mean,
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
opts$K                <- 6           # Number of clusters
opts$N                <- length(met) # Number of cells
opts$M                <- length(met[[1]]) # Number of genomic regions
opts$delta_0          <- rep(3, opts$K)   # Dirichlet prior
opts$alpha_0          <- 0.5         # Gamma prior
opts$beta_0           <- 15          # Gamma prior
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.5         # % of CpGs kept for training in each region
opts$is_kmeans        <- TRUE        # Use K-means for initialization
opts$vb_max_iter      <- 400         # Maximum VB iterations
opts$epsilon_conv     <- 1e-4        # Convergence threshold for VB
opts$vb_init_nstart   <- 5           # Mini VB restarts
opts$vb_init_max_iter <- 20          # Mini VB iteratiions
opts$is_parallel      <- TRUE        # Use parallelized version
opts$no_cores         <- 4           # Number of cores
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
no_cores_out <- 10
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = no_cores_out,
                                          is_parallel = TRUE,
                                          max_cores = opts$total_sims)
print(date())
model <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
    imputation_analysis(X = met, opts = opts), mc.cores = no_cores_out)
print(date())

##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
obj <- list(model = model, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "final_melissa_sim", opts$total_sims,
                           "_", io$data_file,
                           "_cov", io$cov,
                           "_sd", io$sd,
                           "_K", opts$K,
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
