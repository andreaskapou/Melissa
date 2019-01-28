# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(ROCR))
set.seed(123)

##------------------------------------
# Load preprocessed data
##------------------------------------
io <- list(dataset = "ENCODE", data_file = "prom5k", cov = 8, sd = 0.05)
io$data_dir = "../../../local-data/melissa/"
io$out_dir = paste0(io$data_dir, io$dataset, "/imputation/")
dt <- readRDS(paste0(io$data_dir, "met/filtered_met/", io$dataset, "/", io$data_file,
                     "_cov", io$cov, "_sd", io$sd, ".rds"))

##------------------------------------
# Initialize parameters
##------------------------------------
opts                  <- dt$opts
opts$N                <- length(dt$met) # Number of cells
opts$M                <- length(dt$met[[1]]) # Number of genomic regions
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.2         # % of CpGs kept for training in each region
opts$is_parallel      <- TRUE        # Use parallelized version
opts$no_cores         <- 3           # Number of cores
opts$total_sims       <- 10          # Number of simulations
opts$basis_prof       <- create_rbf_object(M = 9) # Profile basis functions
opts$basis_mean       <- create_rbf_object(M = 0) # Rate basis function

##----------------------------------------------
# Filtering low covered regions across cells
##----------------------------------------------
dt <- filter_regions_across_cells(dt = dt, opts = opts)
anno_region <- dt$anno_region
annos <- dt$annos
met <- dt$met
opts$cell_names <- names(met)
opts$M <- length(met[[1]])  # Number of genomic regions
print(opts$M)
rm(dt)

# Run model
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = opts$total_sims,
                                          is_parallel = TRUE,
                                          max_cores = opts$total_sims)
print(date())
message(io$data_file)
message(opts$basis_prof$M)
model <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
    indep_imputation_analysis(X = met, opts = opts), mc.cores = no_cores_out)
print(date())


##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
obj <- list(model = model, annos = annos, anno_region = anno_region, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "indep_sim", opts$total_sims,
                           "_", io$data_file,
                           "_cov", io$cov,
                           "_sd", io$sd,
                           "_M", opts$M,
                           "_basis", opts$basis_prof$M,
                           "_dataPrcg", opts$data_train_prcg,
                           "_regionPrcg", opts$region_train_prcg,
                           "_cpgPrcg", opts$cpg_train_prcg,
                           "_filter", opts$filt_region_cov, ".rds") )
##----------------------------------------------------------------------
message(io$data_file)
message(opts$basis_prof$M)
message("Done ...")
##----------------------------------------------------------------------
