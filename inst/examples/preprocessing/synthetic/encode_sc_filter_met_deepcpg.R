# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) { cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir) }
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
R.utils::sourceDirectory("../../lib", modifiedOnly = FALSE)

dataset   <- "ENCODE/"
data_file <- "prom5k"
sub_dir   <- "/"
data_dir  <- paste0("../../local-data/", dataset, "scBS-seq/deepcpg/processed/unfiltered/")#, sub_dir)
out_dir   <- paste0("../../local-data/melissa/met/filtered_met/", dataset)

# Load deepcpg data
obj <- readRDS(paste0(data_dir, data_file, ".rds"))

# Update options
opts <- obj$opts
opts$filter_chr <- c("1", "3", "5")
opts$cov <- 8            # CpG density at each source
#opts$cov_deepcpg <- 3   # CpG density at each source
opts$met_sd <- 0.05      # Variability

##------------------
# Filter by chromosomes
##------------------
annos <- as.data.table(obj$annos)
idx_filter <- which(annos$seqnames %in% opts$filter_chr)
obj$anno_region <- obj$anno_region[-idx_filter]
obj$annos <- obj$annos[-idx_filter]
obj$met <- lapply(obj$met, function(x) x[-idx_filter])

# Consider only regions with enough CpG coverage
met <- lapply(obj$met, function(x) lapply(x, function(y){
    if (NROW(y) < opts$cov) return(NA) else return(y) }))

##-------------------
# Filter to match original data
##-------------------
# Load original data to map with DeepCpG output
init_obj <- readRDS(paste0(out_dir, data_file, "_cov", opts$cov, "_sd", opts$met_sd, ".rds"))
opts$N   <- length(init_obj$met)      # Number of cells
opts$M   <- length(init_obj$met[[1]]) # Number of genomic regions
opts$filt_region_cov  <- 0.5          # Filter low covered genomic regions
# Filter by region coverage across cells
init_obj <- filter_regions_across_cells(dt = init_obj, opts = opts)

# Map genomic regions by ID
idx <- which(obj$anno_region$id %in% init_obj$anno_region$id)
met <- lapply(met, function(x) x[idx])
anno_region <- obj$anno_region[idx]
annos <- obj$annos[idx]

# Create object and store results
obj <- list(met = met, anno_region = anno_region, annos = annos, opts = opts, io = obj$io)
saveRDS(obj, file = paste0(out_dir, "deepcpg/", sub_dir, data_file, "_cov", opts$cov, "_sd",
                           opts$met_sd, ".rds"))
