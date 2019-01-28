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
io$base_dir   <- "../../local-data/encode/scWGBS_005/"
io$out_dir     <- paste0(io$base_dir, "/parsed/binarised/deepcpg_chr_filtered/")
io$raw_met_dir <- paste0(io$base_dir, "/parsed/binarised/")
io$met_files  <- sub(".tsv.gz", "", list.files(io$raw_met_dir, pattern = "*.gz", full.names = FALSE))
dir.create(io$out_dir)

# Create methylation regions
cores <- 3
met <- mclapply(X = io$met_files, FUN = function(m_file){
    # Read scBS seq data
    met_dt <- read_met(file = sprintf("zcat < %s/%s.tsv.gz", io$raw_met_dir, m_file), type = "sc_seq",
                       strand_info = FALSE) %>% as.data.table %>%
        setnames(c("chr", "start", "end", "width", "strand", "rate")) %>%
        .[,c("end", "width", "strand") := NULL] %>%
        subset(chr %in% c("1", "2", "3", "4", "5", "6"))
    # Write the filtered file
    fwrite(met_dt, file = paste0(io$out_dir, m_file, ".tsv.gz"), col.names = FALSE)
    # Save results
    fwrite(met_dt, file = paste0(io$out_dir, m_file, ".tsv"), sep = "\t", showProgress = FALSE,
           verbose = FALSE, col.names = FALSE)
    #system(sprintf("pigz -p %d -f %s", 1, paste0(io$out_dir, m_file, ".tsv")))
    system(sprintf("gzip -f %s", paste0(io$out_dir, m_file, ".tsv")))
}, mc.cores = cores)
