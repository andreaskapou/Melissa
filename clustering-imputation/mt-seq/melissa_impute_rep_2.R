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
io <- list(dataset = "mt-seq", data_file = "super_enhancers", cov = 10, sd = 0.2)
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
opts$delta_0          <- rep(3, opts$K) + rbeta(opts$K, 1e-1, 1e1)   # Dirichlet prior
opts$alpha_0          <- .5          # Gamma prior
opts$beta_0           <- NULL        # Gamma prior (if NULL beta_0 := sqrt(alpha_0 + M*D/2))
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.5         # % of CpGs kept for training in each region
opts$is_kmeans        <- TRUE        # Use K-means for initialization
opts$vb_max_iter      <- 400         # Maximum VB iterations
opts$epsilon_conv     <- 1e-4        # Convergence threshold for VB
opts$vb_init_nstart   <- 10          # Mini VB restarts
opts$vb_init_max_iter <- 20          # Mini VB iteratiions
opts$is_parallel      <- TRUE        # Use parallelized version
opts$no_cores         <- 3           # Number of cores
opts$total_sims       <- 10          # Number of simulations
opts$basis_prof       <- create_rbf_object(M = 11) # Profile basis functions
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

# Run model
no_cores_out <- 10
no_cores_out <- BPRMeth:::.parallel_cores(no_cores = no_cores_out,
                                          is_parallel = TRUE,
                                          max_cores = opts$total_sims)
print(date())
message(io$data_file)
message(opts$basis_prof$M)
model <- parallel::mclapply(X = 1:opts$total_sims, FUN = function(sim)
    melissa_imputation_analysis(X = met, opts = opts), mc.cores = no_cores_out)
print(date())

##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
obj <- list(model = model, annos = annos, anno_region = anno_region, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "diffuse_melissa_sim", opts$total_sims,
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
message(io$data_file)
message(opts$basis_prof$M)
message("Done ...")
##----------------------------------------------------------------------
