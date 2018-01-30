# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(microbenchmark))

##-----------------------------
# Initialize main variables   #
##-----------------------------
set.seed(15)
opts <- list()
opts$K                <- 3      # Number of clusters
opts$M                <- 200    # Number of genomic regions
opts$cluster_var      <- .6     # % of cluster variability across regions
opts$is_kmeans        <- FALSE  # Use K-means for initialization
opts$vb_max_iter      <- 400    # Maximum VB iterations
opts$epsilon_conv     <- 1e-4   # Convergence threshold for VB
opts$vb_init_nstart   <- 2      # Mini VB restarts
opts$vb_init_max_iter <- 10     # Mini VB max iterations
opts$gibbs_nsim       <- 3000   # Gibbs simulations
opts$gibbs_burn_in    <- 200    # Burn in period
opts$is_parallel      <- TRUE   # Use parallelized version
opts$no_cores         <- 2     # Number of cores
opts$basis            <- create_rbf_object(M = 4) # Basis object

# Number of cells
opts$N <- c(50, 100, 200, 500, 1000, 2000)
res <- vector("list", length = length(opts$N))
iter <- 1
for (n in opts$N) {
    print(n)
    # Generate synthetic data from the generative model
    synth_data <- generate_synth_data(basis = opts$basis, N = n, M = opts$M, K = opts$K,
                                      pi_k = c(0.3, 0.4, 0.3), cluster_var = opts$cluster_var)
    res[[iter]] <- microbenchmark("vb" = {
        vb_obj <- melissa_vb(X = synth_data$X, K = opts$K, basis = opts$basis,
                             vb_max_iter = opts$vb_max_iter, epsilon_conv = opts$epsilon_conv,
                             is_kmeans = opts$is_kmeans, vb_init_nstart = opts$vb_init_nstart,
                             vb_init_max_iter = opts$vb_init_max_iter, is_parallel = opts$is_parallel,
                             no_cores = opts$no_cores, is_verbose = FALSE)
    },"gibbs" = {
        gibbs_obj <- melissa_gibbs(x = synth_data$X, K = opts$K, basis = opts$basis,
                                   lambda = 1/8, gibbs_nsim = opts$gibbs_nsim,
                                   gibbs_burn_in = opts$gibbs_burn_in,
                                   inner_gibbs = FALSE, is_parallel = opts$is_parallel,
                                   no_cores = opts$no_cores)
    }, times = 5L)
    print(res[[iter]])
    iter <- iter + 1
}

out_dir = "../../local-data/melissa/synthetic/model-efficiency/"
saveRDS(res, file = paste0(out_dir, "model_efficiency_K", opts$K, "_M", opts$M, "_N", n, ".rds") )

