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
io <- list(dataset = "smallwood-2014", data_file = "prom10k", cov = 20, sd = 0.2)
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
opts$beta_0           <- NULL        # Gamma prior (if NULL beta_0 := sqrt(alpha_0 + M*D/2))
opts$filt_region_cov  <- 0.5         # Filter low covered genomic regions
opts$data_train_prcg  <- 0.4         # % of data to keep fully for training
opts$region_train_prcg <- 0.95       # % of regions kept for training
opts$cpg_train_prcg   <- 0.5         # % of CpGs kept for training in each region
opts$is_kmeans        <- TRUE        # Use K-means for initialization
opts$vb_max_iter      <- 500         # Maximum VB iterations
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
anno_region <- dt$anno_region
annos <- dt$annos
met <- dt$met
opts$cell_names <- names(met)
opts$M <- length(met[[1]])  # Number of genomic regions
print(opts$M)
rm(dt)


bulk_dt <- matrix(0, ncol = opts$M, nrow = opts$N)
for (m in 1:opts$M) {
    # Extract all cells for specific region
    cells <- lapply(met, "[[", m)
    # Filter to keep only GpC covered regions
    ind <- which(is.na(cells))
    if (length(ind) > 1) { cells_filt <- cells[-ind] }
    # Concatenate to obtain bulk data
    bulk_cell <- do.call(rbind, cells_filt)
    bulk_cov <- length(unique(bulk_cell[,1]))
    for (n in 1:opts$N) {
        if (is.na(cells[[n]])) {
            bulk_dt[n,m] <- 0 / bulk_cov
        }else{
            bulk_dt[n,m] <- NROW(cells[[n]]) / bulk_cov
        }
    }
}

# Create a long vector of elements
bulk_dt_vec <- c(bulk_dt)
# Index of elements with no coverage
zero_ind <- which(bulk_dt_vec == 0)
# Percentage of non-covered regions
non_cov_prcg <- length(zero_ind) / length(c(bulk_dt))
# Coverage percentage including uncovered regions
cov_incl_zero <- mean(bulk_dt)
# Coverage percentage not included uncovered regions
cov_not_incl_zero <- mean(bulk_dt_vec[-zero_ind])


print(paste0("Region ", io$data_file, "\n"))
print(paste0("Percentage of non-covered regions ", round(non_cov_prcg, 2), "\n"))
print(paste0("Coverage percentage including uncovered regions ", round(cov_incl_zero, 2), "\n"))
print(paste0("Coverage percentage not included uncovered regions ", round(cov_not_incl_zero, 2), "\n"))

