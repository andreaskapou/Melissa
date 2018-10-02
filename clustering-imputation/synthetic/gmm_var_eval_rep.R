# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ROCR))
set.seed(123)

# Use the log sum exp trick for having numeric stability
log_sum_exp <- function(x) {
    # Computes log(sum(exp(x))
    offset <- max(x)
    s <- log(sum(exp(x - offset))) + offset
    i <- which(!is.finite(s))
    if (length(i) > 0) { s[i] <- offset }
    return(s)
}


# Fit VBLR model
vb_gmm <- function(X, K = 3, alpha_0 = 1/K, m_0 = c(colMeans(X)),
                   beta_0 = 1, nu_0 = NCOL(X) + 50,
                   W_0 = diag(100, NCOL(X)), max_iter = 500,
                   epsilon_conv = 1e-4, is_verbose = FALSE){
    # Compute logB function
    logB <- function(W, nu){
        D <- NCOL(W)
        return(-0.5*nu*log(det(W)) - (0.5*nu*D*log(2) + 0.25*D*(D - 1) *
                                          log(pi) + sum(lgamma(0.5 * (nu + 1 - 1:NCOL(W)))) ))  #log of B.79
    }
    X <- as.matrix(X)
    D <- NCOL(X)              # Number of features
    N <- NROW(X)              # Number of observations
    W_0_inv <- solve(W_0)     # Compute W^{-1}
    L <- rep(-Inf, max_iter)  # Store the lower bounds
    r_nk = log_r_nk = log_rho_nk <- matrix(0, nrow = N, ncol = K)
    x_bar_k <- matrix(0, nrow = D, ncol = K)
    S_k = W_k <- array(0, c(D, D, K ) )
    log_pi = log_Lambda <- rep(0, K)

    m_k     <- t(kmeans(X, K, nstart = 25)$centers)  # Mean of Gaussian
    beta_k  <- rep(beta_0, K)                # Scale of precision matrix
    nu_k    <- rep(nu_0, K)                  # Degrees of freedom
    alpha   <- rep(alpha_0, K)               # Dirichlet parameter
    log_pi  <- digamma(alpha) - digamma(sum(alpha))
    for (k in 1:K) {
        W_k[,,k] <-  W_0  # Scale matrix for Wishart
        log_Lambda[k] <- sum(digamma((nu_k[k] + 1 - c(1:D))/2)) +
            D*log(2) + log(det(W_k[,,k]))
    }

    # Iterate to find optimal parameters
    for (i in 2:max_iter) {
        ##-------------------------------
        # Variational E-Step
        ##-------------------------------
        for (k in 1:K) {
            diff <- sweep(X, MARGIN = 2, STATS = m_k[, k], FUN = "-")
            log_rho_nk[, k] <- log_pi[k] + 0.5*log_Lambda[k] -
                0.5*(D/beta_k[k]) - 0.5*nu_k[k] * diag(diff %*%W_k[,,k] %*% t(diff)) # log of 10.67
        }
        # Responsibilities using the logSumExp for numerical stability
        Z        <- apply(log_rho_nk, 1, log_sum_exp)
        log_r_nk <- log_rho_nk - Z              # log of 10.49
        r_nk     <- apply(log_r_nk, 2, exp)     # 10.49

        ##-------------------------------
        # Variational M-Step
        ##-------------------------------
        N_k <- colSums(r_nk) + 1e-10  # 10.51
        for (k in 1:K) {
            x_bar_k[, k] <- (r_nk[ ,k] %*% X) / N_k[k]   # 10.52
            x_cen        <- sweep(X,MARGIN = 2,STATS = x_bar_k[, k],FUN = "-")
            S_k[, , k]   <- t(x_cen) %*% (x_cen * r_nk[, k]) / N_k[k]  # 10.53
        }
        # Update Dirichlet parameter
        alpha <- alpha_0 + N_k  # 10.58
        # # Compute expected value of mixing proportions
        pi_k <- (alpha_0 + N_k) / (K * alpha_0 + N)
        # Update parameters for Gaussia-nWishart distribution
        beta_k <- beta_0 + N_k    # 10.60
        nu_k   <- nu_0 + N_k + 1  # 10.63
        for (k in 1:K) {
            # 10.61
            m_k[, k]   <- (1/beta_k[k]) * (beta_0*m_0 + N_k[k]*x_bar_k[, k])
            # 10.62
            W_k[, , k] <- W_0_inv + N_k[k] * S_k[,,k] +
                ((beta_0*N_k[k])/(beta_0 + N_k[k])) *
                tcrossprod((x_bar_k[, k] - m_0))
            W_k[, , k] <- solve(W_k[, , k])
        }
        # Update expectations over \pi and \Lambda
        # 10.66
        log_pi <- digamma(alpha) - digamma(sum(alpha))
        for (k in 1:K) { # 10.65
            log_Lambda[k] <- sum(digamma((nu_k[k] + 1 - 1:D)/2)) +
                D*log(2) + log(det(W_k[,,k]))
        }

        ##-------------------------------
        # Variational lower bound
        ##-------------------------------
        lb_px = lb_pml = lb_pml2 = lb_qml <- 0
        for (k in 1:K) {
            # 10.71
            lb_px <- lb_px + N_k[k] * (log_Lambda[k] - D/beta_k[k] - nu_k[k] *
                                           matrix.trace(S_k[,,k] %*% W_k[,,k]) - nu_k[k]*t(x_bar_k[,k] -
                                                                                               m_k[,k]) %*% W_k[,,k] %*% (x_bar_k[,k] - m_k[,k]) - D*log(2*pi) )
            # 10.74
            lb_pml <- lb_pml + D*log(beta_0/(2*pi)) + log_Lambda[k] -
                (D*beta_0)/beta_k[k] - beta_0*nu_k[k]*t(m_k[,k] - m_0) %*%
                W_k[,,k] %*% (m_k[,k] - m_0)
            # 10.74
            lb_pml2 <- lb_pml2 + nu_k[k] * matrix.trace(W_0_inv %*% W_k[,,k])
            # 10.77
            lb_qml <- lb_qml + 0.5*log_Lambda[k] + 0.5*D*log(beta_k[k]/(2*pi)) -
                0.5*D - logB(W = W_k[,,k], nu = nu_k[k]) -
                0.5*(nu_k[k] - D - 1)*log_Lambda[k] + 0.5*nu_k[k]*D
        }
        lb_px  <- 0.5 * lb_px             # 10.71
        lb_pml <- 0.5*lb_pml + K*logB(W = W_0,nu = nu_0) + 0.5*(nu_0 - D - 1) *
            sum(log_Lambda) - 0.5*lb_pml2 # 10.74
        lb_pz  <- sum(r_nk %*% log_pi)    # 10.72
        lb_qz  <- sum(r_nk * log_r_nk)    # 10.75
        lb_pp  <- sum((alpha_0 - 1)*log_pi) + lgamma(sum(K*alpha_0)) -
            K*sum(lgamma(alpha_0))        # 10.73
        lb_qp  <- sum((alpha - 1)*log_pi) + lgamma(sum(alpha)) -
            sum(lgamma(alpha)) # 10.76
        # Sum all parts to compute lower bound
        L[i] <- lb_px + lb_pz + lb_pp + lb_pml - lb_qz - lb_qp - lb_qml

        ##-------------------------------
        # Evaluate mixture density for plotting
        ##-------------------------------
        # Show VB difference
        if (is_verbose) { cat("It:\t",i,"\tLB:\t",L[i],
                              "\tLB_diff:\t",L[i] - L[i - 1],"\n")}
        # Check if lower bound decreases
        if (L[i] < L[i - 1]) { message("Warning: Lower bound decreases!\n"); }
        # Check for convergence
        if (abs(L[i] - L[i - 1]) < epsilon_conv) { break }
        # Check if VB converged in the given maximum iterations
        if (i == max_iter) {warning("VB did not converge!\n")}
    }
    obj <- structure(list(X = X, K = K, N = N, D = D, pi_k = pi_k,
                          alpha = alpha, r_nk = r_nk,  W = m_k, W_Sigma = W_k,
                          beta = beta_k, nu = nu_k, lb = L[2:i]), class = "vb_gmm")
    return(obj)
}


gmm_var_analysis <- function(opts, sim){
    # Initialize lists
    gmm = eval_perf <- vector("list", length = length(opts$cluster_var))
    # Load synthetic data
    io <- list(data_file = paste0("encode_data_", sim, ".rds"),
               data_dir = "../../local-data/melissa/synthetic/imputation/coverage/raw/data-sims/")
    obj <- readRDS(paste0(io$data_dir, io$data_file))
    i <- 1
    # Iterate
    for (cluster_var in opts$cluster_var) {
        # Load synthetic data
        io <- list(data_file = paste0("encode_data_", cluster_var, "_", sim, ".rds"),
                   data_dir = "../../local-data/melissa/synthetic/imputation/dissimilarity/raw/data-sims/")
        obj <- readRDS(paste0(io$data_dir, io$data_file))
        # Partition to training and test sets
        dt <- partition_dataset(X = obj$synth_data$X, region_train_prcg = opts$region_train_prcg,
                                cpg_train_prcg = opts$cpg_train_prcg, is_synth = TRUE)
        # Convert list to a big MxN matrix
        X <- matrix(0, nrow = opts$N, ncol = opts$M)
        # Compute mean methylation rate for each cell and region
        for (n in 1:opts$N) {
            X[n, ] <- sapply(dt$train[[n]], function(m) mean(m[,2]))
        }
        # Transform to Gaussian data using M values
        X <- log2( (X + 0.01) / (1 - X + 0.01) )

        # Using GMM
        gmm_obj <- vb_gmm(X = X, K = opts$K, alpha_0 = opts$alpha_0, beta_0 = opts$beta_0,
                          max_iter = opts$max_iter, epsilon_conv = opts$epsilon_conv, is_verbose = FALSE)
        gmm[[i]] <- gmm_obj
        # Convert predicted means to (0,1)
        W <- pnorm(gmm_obj$W)

        test_pred <- dt$test                              # Copy test data
        act_obs   <- vector(mode = "list", length = opts$N) # Keep actual CpG states
        pred_obs  <- vector(mode = "list", length = opts$N) # Keep predicted CpG states
        for (n in 1:opts$N) {
            # Cluster assignment
            k <- which.max(gmm_obj$r_nk[n,])
            # Obtain non NA observations
            idx <- which(!is.na(dt$test[[n]]))
            cell_act_obs  <- vector("list", length = length(idx))
            cell_pred_obs <- vector("list", length = length(idx))
            iter <- 1
            # Iterate over each genomic region
            for (m in idx) {
                test_pred[[n]][[m]][,2] <- W[m,k]
                # TODO
                # Actual CpG states
                cell_act_obs[[iter]]  <- dt$test[[n]][[m]][, 2]
                # Predicted CpG states (i.e. function evaluations)
                cell_pred_obs[[iter]] <- test_pred[[n]][[m]][, 2]
                iter <- iter + 1
            }
            # Combine data to a big vector
            act_obs[[n]] <- do.call("c", cell_act_obs)
            pred_obs[[n]] <- do.call("c", cell_pred_obs)
        }
        eval_perf[[i]] <- list(act_obs = do.call("c", act_obs), pred_obs = do.call("c", pred_obs))

        ##----------------------------------------------------------------------
        message("Computing AUC...")
        ##----------------------------------------------------------------------
        pred_gmm <- prediction(do.call("c", pred_obs), do.call("c", act_obs))
        auc_gmm <- performance(pred_gmm, "auc")
        auc_gmm <- unlist(auc_gmm@y.values)
        message(auc_gmm)
        i <- i + 1 # Increase counter
    }
    obj <- list(gmm = gmm, eval_perf = eval_perf, opts = opts)
    return(obj)
}


##------------------------
# Load synthetic data
##------------------------
io <- list(data_file = paste0("raw/data-sims/encode_data_0.1_1.rds"),
           out_dir = "../../local-data/melissa/synthetic/imputation/dissimilarity/")
obj <- readRDS(paste0(io$out_dir, io$data_file))
opts                   <- obj$opts       # Get options
opts$alpha_0           <- rep(3, opts$K) # Dirichlet prior
opts$beta_0            <- 1             # Precision prior
opts$data_train_prcg   <- 0.1            # % of data to keep fully for training
opts$region_train_prcg <- 1              # % of regions kept for training
opts$cpg_train_prcg    <- 0.4            # % of CpGs kept for training in each region
opts$max_iter          <- 500            # Maximum VB iterations
opts$epsilon_conv      <- 1e-4           # Convergence threshold for VB
rm(obj)

print(date())
obj <- lapply(X = 1:opts$total_sims, FUN = function(sim)
    gmm_var_analysis(opts = opts, sim = sim))
print(date())

##----------------------------------------------------------------------
message("Storing results...")
##----------------------------------------------------------------------
saveRDS(obj, file = paste0(io$out_dir, "encode_gmm_K", opts$K,
                           "_dataTrain", opts$data_train_prcg,
                           "_regionTrain", opts$region_train_prcg,
                           "_cpgTrain", opts$cpg_train_prcg, ".rds") )
