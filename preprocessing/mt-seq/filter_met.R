# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) { cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir) }
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))

data_file  <- "prom3k"
data_dir <- "../local-data/mt-seq/met/processed/unfiltered/"
out_dir <- "../local-data/melissa/met/filtered_met/"

# Load data
obj <- readRDS(paste0(data_dir, data_file, ".rds"))
met <- obj$met
anno_region <- obj$anno_region
annos <- obj$annos
rna <- obj$rna

# Update options
opts <- obj$opts
opts$cov <- 10          # CpG density at each source
opts$met_sd <- 0.2      # Metylation standard deviation across cells
opts$gene_var <- 5      # Gene variance across cells

# Only for promoter regions, so we can reduce dimensionlity of dataset!
# # Filter non varying genes across cells
# rna <- rna[, var := apply(.SD, 1, mean), .SDcols = 2:62] %>%
#   .[var > opts$gene_var] %>% .[, var := NULL]
# # # Subset annotation data
# annos <- annos %>% subset(id %in% rna$id)
# # # Subset promoter regions
# anno_region <- anno_region %>% subset(id %in% rna$id)
# # # Subset methylation data
# met <- lapply(met, function(x) x[rna[, id]])

# Consider only regions with enough CpG coverage
met <- lapply(met, function(x) lapply(x, function(y){
    if (NROW(y) < opts$cov) return(NA) else return(y) }))
# Number of genomic regions
M <- NROW(anno_region)
cell_region_sd <- vector("numeric", length = M)
for (m in 1:M) {
    # Extract all cell observations for specific region m
    tmp <- lapply(met, "[[", m)
    # Keep only regions that have observations
    tmp <- tmp[!is.na(tmp)]
    if (length(tmp) == 0) { cell_region_sd[m] <- 0
    # Compute the standard deviation of region across cells
    } else {cell_region_sd[m] <- sd(sapply(tmp, function(x) mean(x[,2]))) }
}
# Plot histogram of variability across cells
hist(cell_region_sd, breaks = 50)
# Keep only highly varying sources
ind <- which(cell_region_sd > opts$met_sd)

# Subset according to methylation variability
rna <- rna[ind,]; annos <- annos[ind,]; anno_region <- anno_region[ind,]
met <- lapply(met, function(x) x[ind])

# Create object and store results
obj <- list(met = met, anno_region = anno_region, annos = annos,
            rna = rna, opts = opts, io = obj$io)
saveRDS(obj, file = paste0(out_dir, data_file, "_cov", opts$cov, "_sd",
                           opts$met_sd, "_gene_var", opts$gene_var, ".rds"))
