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
io <- list(dataset = "smallwood-2014", data_file = "Nanog", cov = 10, sd = 0.2)
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


# bulk_dt <- matrix(0, ncol = opts$M, nrow = opts$N)
# for (m in 1:opts$M) {
#     # Extract all cells for specific region
#     cells <- lapply(met, "[[", m)
#     # Filter to keep only GpC covered regions
#     ind <- which(is.na(cells))
#     if (length(ind) > 1) { cells_filt <- cells[-ind] }
#     # Concatenate to obtain bulk data
#     bulk_cell <- do.call(rbind, cells_filt)
#     bulk_cov <- length(unique(bulk_cell[,1]))
#     for (n in 1:opts$N) {
#         if (is.na(cells[[n]])) {
#             bulk_dt[n,m] <- 0 / bulk_cov
#         }else{
#             bulk_dt[n,m] <- NROW(cells[[n]]) / bulk_cov
#         }
#     }
# }
#
# # Create a long vector of elements
# bulk_dt_vec <- c(bulk_dt)
# # Index of elements with no coverage
# zero_ind <- which(bulk_dt_vec == 0)
# # Percentage of non-covered regions
# non_cov_prcg <- length(zero_ind) / length(c(bulk_dt))
# # Coverage percentage including uncovered regions
# cov_incl_zero <- mean(bulk_dt)
# # Coverage percentage not included uncovered regions
# cov_not_incl_zero <- mean(bulk_dt_vec[-zero_ind])
#
#
# print(paste0("Region ", io$data_file, "\n"))
# print(paste0("Percentage of non-covered regions ", round(non_cov_prcg, 2), "\n"))
# print(paste0("Coverage percentage including uncovered regions ", round(cov_incl_zero, 2), "\n"))
# print(paste0("Coverage percentage not included uncovered regions ", round(cov_not_incl_zero, 2), "\n"))



###------------------------------
####
###

## Compute coverage of CpGs across the region and create a histogram
cpg_cov_dt <- matrix(0, ncol = opts$M, nrow = opts$N)
for (m in 1:opts$M) {
  # Extract all cells for specific region
  cells <- lapply(met, "[[", m)
  cpg_cov_dt[,m] <- unlist(lapply(cells, function(n) NROW(n)))
}

# Create a vector
cpg_cov_dt_vec <- c(cpg_cov_dt)
cpg_cov_dt_vec <- cpg_cov_dt_vec[cpg_cov_dt_vec > 1 & cpg_cov_dt_vec < 200]
hist(c(cpg_cov_dt_vec), breaks = 50, main = "Nanog regions", col = "lavender",
     xlab = paste0("CpG coverage"))



total_iter <- 11
coverage_regions <- matrix(0, ncol = opts$N, nrow = total_iter - 1)
for (iter in 2:total_iter) {
  if (iter == 2) {
    coverage_regions[iter - 1, ] <- sapply(1:NROW(cpg_cov_dt), function(rows)  sum(cpg_cov_dt[rows, ] <= iter*10 &
                                                                                     cpg_cov_dt[rows, ] != 1))
  }else if (iter == total_iter) {
    coverage_regions[iter - 1, ] <- sapply(1:NROW(cpg_cov_dt), function(rows)  sum(cpg_cov_dt[rows, ] > iter*10 &
                                                                                     cpg_cov_dt[rows, ] != 1))
  }else{
    coverage_regions[iter - 1, ] <- sapply(1:NROW(cpg_cov_dt), function(rows)  sum(cpg_cov_dt[rows, ] > (iter - 1)*10 &
                                                                                     cpg_cov_dt[rows, ] < iter*10 &
                                                                                     cpg_cov_dt[rows, ] != 1))
  }
}

dt <- as.data.table(coverage_regions)
colnames(dt) <- paste0("cell", 1:NCOL(coverage_regions))
dt$cov <- c("N<20", "20<N<30", "30<N<40", "40<N<50", "50<N<60", "60<N<70", "70<N<80", "80<N<90", "90<N<100", "N>100")
dt_melt <- melt(dt, id.vars = "cov")

dt_melt <- dt_melt %>%
  .[, cov := factor(cov, levels = c("N<20", "20<N<30", "30<N<40", "40<N<50", "50<N<60", "60<N<70", "70<N<80", "80<N<90", "90<N<100", "N>100"))]


# Define ggplot2 theme for line plots
boxplot_theme <- function(){
  p <- theme(
    plot.title = element_text(size = 15, face = 'bold',
                              margin = ggplot2::margin(0,0,0,0), hjust = 0.5),
    axis.text = element_text(size = rel(1.15), color = 'black'),
    axis.text.x = element_text(size = rel(0.7), color = 'black'),
    axis.title = element_text(size = rel(1.25), color = 'black'),
    axis.title.y = element_text(margin = ggplot2::margin(0,15,0,0)),
    axis.title.x = element_text(margin = ggplot2::margin(5,0,0,0)),
    axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
    axis.ticks.y = element_blank(),
    legend.position = "left",
    legend.key.size = unit(1.9, 'lines'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 19),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "gainsboro"),
    panel.background = element_blank()
  )
  return(p)
}

if (io$data_file == "prom10k") {title <- "Promoter 10kb regions"
}else if (io$data_file == "prom5k") {title <- "Promoter 5kb regions"
}else if (io$data_file == "prom3k") {title <- "Promoter 3kb regions"
}else if (io$data_file == "active_enhancers") {title <- "Active enhancers regions"
}else if (io$data_file == "super_enhancers") {title <- "Super enhancers regions"
}else if (io$data_file == "Nanog") {title <- "Nanog regions"}

p_boxplot <- ggplot(dt_melt, aes(x = cov, y = value)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, fill = "lavender") +
  # scale_y_continuous(limits = c(0.64, 0.79), breaks = pretty_breaks(n = 4)) +
  labs(title = title, x = "CpG coverage", y = "Number of regions") +
  boxplot_theme()
print(p_boxplot)


pdf(file = paste0("out/boxplots/box-", io$dataset, "-", io$data_file, ".pdf"), width = 7, height = 5, useDingbats = FALSE)
print(p_boxplot)
dev.off()

# pdf(file = paste0("out/", io$dataset, "-", io$data_file, ".pdf"), width = 5, height = 4, useDingbats = FALSE)
# hist(c(cpg_cov_dt_vec), breaks = 50, main = title, col = "lavender",
#      xlab = paste0("CpG coverage"))
# dev.off()
