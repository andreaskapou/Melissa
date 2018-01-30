# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))

# Initialize parameters
set.seed(10)
K <- 2
I <- 40
N <- 20
gibbs_nsim <- 30
gibbs_burn_in <- 10
basis <- create_rbf_object(M = 3)

# Generate synthetic data from the generative model
synth_data <- generate_synth_data(basis = basis, N = I, M = N, K = K,
                                  pi_k = c(0.3, 0.7), cluster_var = 0.5)

# Compute the posterior distribution using the Baysian BPR FDMM model
melissa_res <- melissa_gibbs(x = synth_data$X, K = K, basis = basis,
                             lambda = 1/8, gibbs_nsim = gibbs_nsim,
                             gibbs_burn_in = gibbs_burn_in, inner_gibbs = FALSE,
                             is_parallel = FALSE, no_cores = NULL)
