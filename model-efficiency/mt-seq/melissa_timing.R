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

##------------------------------------
# Load preprocessed data
##------------------------------------
io <- list(dataset = "mt-seq", data_file = "prom5k", cov = 10, sd = 0.2)
io$data_dir = "../../local-data/melissa/"
io$out_dir = paste0(io$data_dir, io$dataset, "/imputation/")
dt <- readRDS(paste0(io$data_dir, "met/filtered_met/", io$dataset, "/", io$data_file, "_cov",
                     io$cov, "_sd", io$sd, "_gene_var5.rds"))

##------------------------------------
# Initialize parameters
##------------------------------------
opts                  <- dt$opts
opts$K                <- 6           # Number of clusters
opts$N                <- length(dt$met) # Number of cells
opts$M                <- length(dt$met[[1]]) # Number of genomic regions
opts$delta_0          <- rep(.1, opts$K)   # Dirichlet prior
opts$alpha_0          <- 0.5         # Gamma prior
opts$beta_0           <- 15          # Gamma prior
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.5         # % of CpGs kept for training in each region
opts$is_kmeans        <- TRUE        # Use K-means for initialization
opts$vb_max_iter      <- 400         # Maximum VB iterations
opts$epsilon_conv     <- 1e-4        # Convergence threshold for VB
opts$vb_init_nstart   <- 5          # Mini VB restarts
opts$vb_init_max_iter <- 10          # Mini VB iteratiions
opts$is_parallel      <- TRUE        # Use parallelized version
opts$no_cores         <- 7           # Number of cores
opts$basis_prof       <- create_rbf_object(M = 7) # Profile basis functions
opts$basis_mean       <- create_rbf_object(M = 0) # Rate basis function

##----------------------------------------------
# Filtering low covered regions across cells
##----------------------------------------------
dt <- filter_regions_across_cells(dt = dt, opts = opts)
met <- dt$met
annos <- dt$annos
anno_region <- dt$anno_region
opts$M <- length(met[[1]])  # Number of genomic regions
print(opts$M)
rm(dt)


print(date())
message(io$data_file)
message(opts$basis_prof$M)
# Partition to training and test sets
dt <- partition_dataset(X = met, data_train_prcg = opts$data_train_prcg,
                        region_train_prcg = opts$region_train_prcg,
                        cpg_train_prcg = opts$cpg_train_prcg, is_synth = FALSE)
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
print(date())

