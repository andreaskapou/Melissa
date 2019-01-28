# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))

##-----------------------------
# Initialize main variables   #
##-----------------------------
set.seed(15)
K <- 2                 # Number of clusters
N <- 20                # Number of cells
M <- 20                # Number of genomic regions
cluster_var  <- 0.8    # % of cluster variability across regions
is_kmeans    <- TRUE   # Use K-means for initialization
is_parallel  <- TRUE   # Use parallelized version
no_cores     <- 2      # Number of cores
vb_max_iter  <- 100    # Maximum VB iterations
epsilon_conv <- 1e-4   # Convergence threshold for VB
vb_init_nstart <- 2    # Mini VB restarts
vb_init_max_iter <- 20 # Mini VB max iterations
basis <- create_rbf_object(M = 4) # Basis object
# Generate synthetic data from the generative model
synth_data <- generate_synth_data(basis = basis, N = N, M = M, K = K,
                                  pi_k = c(0.3, 0.7), cluster_var = cluster_var)

# Compute the posterior distribution using the VB model
melissa_obj <- melissa_vb(X = synth_data$X, K = K, basis = basis,
                        vb_max_iter = vb_max_iter, epsilon_conv = epsilon_conv,
                        is_kmeans = is_kmeans, vb_init_nstart = vb_init_nstart,
                        vb_init_max_iter = vb_init_max_iter, is_parallel = is_parallel,
                        no_cores = no_cores, is_verbose = TRUE)

# Compute clustering assignment errors
ari <- cluster_ari(C_true = synth_data$C_true, C_post = melissa_obj$r_nk)
error <- cluster_error(C_true = synth_data$C_true, C_post = melissa_obj$r_nk)

cat("\n\nResults\n")
cat("Mixing proportions: ", melissa_obj$pi_k, "\n\n")
cat("Best VB:", melissa_obj$lb[length(melissa_obj$lb)], "\n")
cat("Cluster error: ", error, "\n")
cat("Cluster ARI:", ari, "\n")

xs <- seq(-1, 1, len = 100) # create test X values
# Estimate predictive distribution
pred <- matrix(0, nrow = length(xs), ncol = K)
act <- matrix(0, nrow = length(xs), ncol = K)
X_test <- design_matrix(basis, xs)$H
for (m in 1:5) {
    for (k in 1:K) {
        # Predictive mean
        pred[,k] <- pnorm(X_test %*% melissa_obj$W[m,,k])
        act[, k] <- pnorm(X_test %*% synth_data$W[m,,k])
    }
    # Create plot
    p1 <- draw_predictive(xs = xs, pred = pred, title = "Estimated")
    p2 <- draw_predictive(xs = xs, pred = act, title = "True")
    final_fig <- cowplot::plot_grid(p1, p2, labels = c("a", "b"), label_size = 30,
                                    ncol = 2, nrow = 1, rel_widths = c(1, 1))
    print(final_fig)
}
