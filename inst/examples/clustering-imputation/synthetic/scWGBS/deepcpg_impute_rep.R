# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(ROCR))
set.seed(123)

##------------------------------------
# Load preprocessed data
##------------------------------------
io <- list(dataset = "encode/scWGBS", data_file = "prom10k", cov = 20, sd = 0.2)
io$data_dir <- "../../../local-data/melissa/"
io$sub_dir <- "/"
io$out_dir <- paste0(io$data_dir, io$dataset, "/imputation/deepcpg/", io$sub_dir)
dt <- readRDS(paste0(io$data_dir, "met/filtered_met/", io$dataset, "/deepcpg/", io$sub_dir, io$data_file,
                     "_cov", io$cov, "_sd", io$sd, ".rds"))

##------------------------------------
# Initialize parameters
##------------------------------------
opts                  <- dt$opts
opts$N                <- length(dt$met) # Number of cells
opts$M                <- length(dt$met[[1]]) # Number of genomic regions
opts$cell_names       <- names(dt$met)  # Cell names
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.5         # % of CpGs kept for training in each region
opts$no_cores         <- 1           # Number of cores
opts$total_sims       <- 10          # Number of simulations

# Run model
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = opts$total_sims,
                                          is_parallel = TRUE,
                                          max_cores = opts$total_sims)
print(date())
model <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
    deepcpg_imputation_analysis(X = dt$met, opts = opts), mc.cores = 3)
print(date())

##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
obj <- list(model = model, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "deepcpg_sim", opts$total_sims,
                           "_", io$data_file,
                           "_cov", io$cov,
                           "_sd", io$sd,
                           "_M", opts$M,
                           "_dataPrcg", opts$data_train_prcg,
                           "_regionPrcg", opts$region_train_prcg,
                           "_cpgPrcg", opts$cpg_train_prcg,
                           "_filter", opts$filt_region_cov, ".rds") )
##----------------------------------------------------------------------
message("Done ...")
##----------------------------------------------------------------------
