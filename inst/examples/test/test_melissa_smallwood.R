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
#set.seed(123)

##------------------------------------
# Load preprocessed data
##------------------------------------
io <- list(dataset = "smallwood-2014", data_file = "active_enhancers", cov = 10, sd = 0.2)
io$data_dir = "../../local-data/melissa/"
io$out_dir = paste0(io$data_dir, io$dataset, "/imputation/")
dt <- readRDS(paste0(io$data_dir, "met/filtered_met/", io$dataset, "/", io$data_file,
                     "_cov", io$cov, "_sd", io$sd, ".rds"))

##------------------------------------
# Initialize parameters
##------------------------------------
opts                  <- dt$opts
opts$K                <- 5           # Number of clusters
opts$N                <- length(dt$met) # Number of cells
opts$M                <- length(dt$met[[1]]) # Number of genomic regions
opts$delta_0          <- rep(.5, opts$K) + rbeta(opts$K, 1e-1, 1e1)   # Dirichlet prior
opts$alpha_0          <- 0.5         # Gamma prior
opts$beta_0           <- NULL        # Gamma prior
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.5         # % of CpGs kept for training in each region
opts$is_kmeans        <- TRUE        # Use K-means for initialization
opts$vb_max_iter      <- 100         # Maximum VB iterations
opts$epsilon_conv     <- 1e-4        # Convergence threshold for VB
opts$vb_init_nstart   <- 3          # Mini VB restarts
opts$vb_init_max_iter <- 10          # Mini VB iteratiions
opts$is_parallel      <- TRUE        # Use parallelized version
opts$no_cores         <- 2           # Number of cores
opts$total_sims       <- 1          # Number of simulations
opts$basis_prof       <- create_rbf_object(M = 11) # Profile basis functions
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

# Partition to training and test sets
dt <- partition_dataset(X = met, data_train_prcg = opts$data_train_prcg,
                        region_train_prcg = opts$region_train_prcg,
                        cpg_train_prcg = opts$cpg_train_prcg, is_synth = FALSE)
rm(met)

melissa_obj <- melissa_vb(X = dt$train, K = opts$K, basis = opts$basis_prof,
                          delta_0 = opts$delta_0, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                          vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                          is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                          vb_init_max_iter = opts$vb_init_max_iter, is_parallel = TRUE,
                          no_cores = 2, is_verbose = TRUE)

eval_prof <- eval_performance(obj = melissa_obj, test = dt$test)
pred_prof <- prediction(eval_prof$pred_obs, eval_prof$act_obs)
auc_prof <- performance(pred_prof, "auc")
message(unlist(auc_prof@y.values))

##----------------------------------------------------------------------
message(io$data_file)
message(opts$basis_prof$M)
message("Done ...")
##----------------------------------------------------------------------
