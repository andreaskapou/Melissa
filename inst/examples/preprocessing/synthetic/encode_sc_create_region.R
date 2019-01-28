# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) { cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir) }
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))

#
# # Data files
io <- list(anno_name = "promoter")
io$base_dir   <- "../../local-data/ENCODE"
io$out_dir     <- paste0(io$base_dir, "/scBS-seq/processed/unfiltered/")
io$annos_file  <- paste0(io$base_dir, "/annotations/", io$anno_name, ".bed")
io$raw_met_dir <- paste0(io$base_dir, "/scBS-seq/parsed/binarised")
io$met_files  <- sub(".tsv.gz", "", list.files(io$raw_met_dir, pattern = "*.gz", full.names = FALSE))

#
# # Parameter options
opts <- list()
opts$is_centre  <- FALSE   # Whether genomic region is already pre-centred
opts$is_window  <- TRUE    # Use predefined window region
opts$upstream   <- -2500   # Upstream of centre
opts$downstream <- 2500    # Downstream of centre
opts$chrom_size <- NULL    # Chromosome size file
opts$cov        <- 3       # Regions with at least n CpGs
opts$sd_thresh  <- -1      # Threshold for variance of methylation across region

# Read annotation data
annos <- fread(io$annos_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
    setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
    .[,c("anno", "chr") := list(NULL, as.factor(sub("chr", "", chr)))] %>%
    subset(chr %in% c("1", "2", "3", "4", "5", "6")) %>%
    setkey(chr, start, end)
if (io$anno_name == "super_enhancers") {
    annos <- annos %>% .[ (end - start) >= 1000] %>% .[ (end - start) <= 20000]
} else if (io$anno_name == "active_enhancers") {
    annos <- annos %>% .[ (end - start) >= 1000] %>% .[ (end - start) <= 10000]
} else if (io$anno_name == "Nanog") {
    annos <- annos %>% .[ (end - start) >= 1000] %>% .[ (end - start) <= 10000]
}

# If we have promoter region and want pre-centred data
if (opts$is_centre == TRUE) { annos[, strand := "*"] }
annos <- GRanges(annos)
# Create genomic region
anno_region <- create_anno_region(anno = annos, is_centre = opts$is_centre,
                                  is_window = opts$is_window, upstream = opts$upstream,
                                  downstream = opts$downstream)

# Create methylation regions
cores <- 2
met <- mclapply(X = io$met_files, FUN = function(m_file){
    # Read scBS seq data
    met_dt <- read_met(file = sprintf("zcat < %s/%s.tsv.gz", io$raw_met_dir, m_file), type = "sc_seq",
                       strand_info = FALSE)
    # Create promoter methylation regions
    res <- create_region_object(met_dt = met_dt, anno_dt = anno_region,
                                cov = opts$cov, sd_thresh = opts$sd_thresh,
                                ignore_strand = TRUE, filter_empty_region = FALSE)$met
    names(res) <- NULL
    return(res)
}, mc.cores = cores)

names(met) <- io$met_files

#
# # Store the results
message("Storing results...")
obj <- list(met = met, anno_region = anno_region, annos = annos, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "prom5k.rds"))
