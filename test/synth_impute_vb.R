# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(require(ROCR))
set.seed(1)
##-----------------------------
# Initialize main variables   #
##-----------------------------
K <- 3                   # Number of clusters
N <- 50                  # Number of cells
M <- 30                  # Number of genomic regions
cluster_var  <- 0.8      # % of cluster variability across regions
region_train_prcg <- 0.7 # % of regions kept for training
cpg_train_prcg <- 0.1    # % of CpGs kept for training in each region
is_kmeans    <- TRUE     # Use K-means for initialization
vb_max_iter  <- 10       # Maximum VB iterations
epsilon_conv <- 1e-4     # Convergence threshold for VB
vb_init_nstart <- 1      # Mini VB restarts
vb_init_max_iter <- 4    # Mini VB max iterations
is_parallel  <- TRUE     # Use parallelized version
no_cores     <- 2        # Number of cores
basis <- create_rbf_object(M = 4) # Basis object
# Generate synthetic data from the generative model
synth_data <- generate_synth_data(basis = basis, N = N, M = M, K = K,
                                  pi_k = c(0.3, 0.5, 0.2), cluster_var = cluster_var)
# Partition to training and test sets
part_data <- partition_dataset(X = synth_data$X, region_train_prcg = region_train_prcg,
                               cpg_train_prcg = cpg_train_prcg, is_synth = TRUE)


# Compute the posterior distribution using the VB model
melissa_obj <- melissa_vb(X = part_data$train, K = K, basis = basis,
                        vb_max_iter = vb_max_iter, epsilon_conv = epsilon_conv,
                        is_kmeans = is_kmeans, vb_init_nstart = vb_init_nstart,
                        vb_init_max_iter = vb_init_max_iter, is_parallel = is_parallel,
                        no_cores = no_cores, is_verbose = TRUE)

# Imputation performance
imp_perf <- eval_performance(obj = melissa_obj, test = part_data$test)
# Compute AUC
pred_prof <- prediction(imp_perf$pred_obs, imp_perf$act_obs)
roc_prof <- performance(pred_prof, "tpr", "fpr")
auc_prof <- performance(pred_prof, "auc")
auc_prof <- unlist(auc_prof@y.values)
message(auc_prof)

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
