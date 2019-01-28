# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) { setwd(dirname(parent.frame(2)$ofile)) }
# Source the 'lib' directory
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(require(BPRMeth))
suppressPackageStartupMessages(require(mvtnorm))
set.seed(1234)

##------------------------
# Initialize variables   #
##------------------------
io <- list(out_dir = "../../local-data/melissa/synthetic/imputation/cells/raw/",
           coef_file = "../../local-data/melissa/synthetic/encode-coef.rds")
opts <- list()
opts$N <- c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250)  # Number of cells
opts$M <- 100               # Number of genomic regions
opts$K <- 4                 # Number of clusters
opts$pi_k <- c(.2,.4,.15,.25)  # Mixing proportions
opts$cluster_var  <- 0.5    # % of regions that are variable across clusters
opts$total_sims   <- 10     # Number of simulations
opts$basis_prof   <- create_rbf_object(M = 4) # Profile basis functions
opts$basis_mean   <- create_rbf_object(M = 0) # Rate basis function

##------------------------
# Create synthetic data
##------------------------
# Load ENCODE cluster profiles
encode_w <- readRDS(io$coef_file)

synth_data = params <- vector("list", length = opts$total_sims)
for (n in opts$N) {
  for (sim in 1:opts$total_sims) {
    synth_data[[sim]] <- generate_encode_synth_data(basis = opts$basis_prof,
                                                    encode_w = encode_w, N = n, M = opts$M, K = opts$K,
                                                    pi_k = opts$pi_k, cluster_var = opts$cluster_var)
  }
  ##-----------------------------
  #   Save synthetic data
  ##-----------------------------
  obj <- list(synth_data = synth_data, io = io, opts = opts)
  saveRDS(obj, file = paste0(io$out_dir, "encode_data_", n, ".rds"))

  # Store each simulation independently due to memory issues when doing inference
  for (sim in 1:opts$total_sims) {
    obj <- list(synth_data = synth_data[[sim]], io = io, opts = opts)
    saveRDS(obj, file = paste0(io$out_dir, "data-sims/encode_data_",
                               n, "_", sim, ".rds"))
  }
}

