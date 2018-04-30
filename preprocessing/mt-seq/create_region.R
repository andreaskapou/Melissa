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
io <- list(anno_name = "active_enhancers")
io$base_dir   <- "../local-data/mt-seq/"
io$sub_dir    <- "subsampled/"
#io$sub_dir    <- "/"
io$out_dir    <- paste0(io$base_dir, "/met/processed/")
io$annos_file <- paste0(io$base_dir, "/annotations/", io$anno_name, ".bed")
io$rna_file   <- paste0(io$base_dir, "/rna/parsed/GSE74534_RNA-seq_normalized_counts.txt.gz")
io$met_dir    <- paste0(io$base_dir, "/met/parsed")
io$met_files  <- list.files(io$met_dir, pattern = "*.gz", full.names = TRUE)

#
# # Parameter options
opts <- list()
opts$upstream   <- -2500   # Upstream of TSS
opts$downstream <- 2500    # Downstream of TSS
opts$chrom_size <- NULL    # Chromosome size file
opts$cov        <- 3       # Promoters with at least n CpGs
opts$sd_thresh  <- -1      # Threshold for variance of methylation across region

#
# # Read scRNA-Seq data
rna <- fread(input = sprintf("zcat < %s", io$rna_file), sep = "\t",
             header = TRUE, showProgress = FALSE, stringsAsFactors = FALSE) %>%
    setnames(c("ens_id"), c("id")) %>% setkey(id)

#
# # Read annotation file
annos <- fread(input = io$annos_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
               showProgress = FALSE) %>%
    setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
    .[, c("anno", "chr") := list(NULL, as.factor(paste0("chr", chr)))] %>%
    .[ (end - start) >= 1000] %>%
    setkey(chr, start, end) #%>% GRanges()

# Create annotation region
anno_region <- copy(annos)
anno_region <- anno_region %>% .[, centre := floor((start + end)/2) ]

# Only for Nanog regions we create a larger genomic region due to sparse CpG coverage
if (io$anno_name == "Nanog") {
    anno_region <- anno_region[, c("start", "end") := list(centre + opts$upstream,
                                                           centre + opts$downstream)]
}
# Create GRanges objects
anno_region <- GRanges(anno_region)
annos <- GRanges(annos)

# # Create methylation regions
met <- list()
for (m_file in io$met_files) {
    # Extract cell ID
    cell <- gsub(".*_([^;]+)\\.CpG.*", "\\1", m_file)
    print(cell)
    # Read scBS seq data
    met_dt <- fread(input = sprintf("zcat < %s", m_file), sep = "\t", header = TRUE,
                    stringsAsFactors = TRUE, showProgress = FALSE) %>%
        .[, c("end", "meth_reads", "unmeth_reads", "strand") :=
              list(start, ifelse((meth_reads + 0.1)/(unmeth_reads + 0.1) > 1, 1, 0), NULL, NULL)] %>%
        .[,c("chr", "start", "end", "meth_reads")] %>%
        setnames(c("chr", "start", "end", "met")) %>% setkey(chr, start, end) %>% GRanges()
    # Create promoter methylation regions
    met[[cell]] <- create_region_object(met_dt = met_dt, anno_dt = anno_region,
                                        cov = opts$cov, sd_thresh = opts$sd_thresh,
                                        ignore_strand = TRUE, filter_empty_region = FALSE)$met
}
rm(met_dt)

#
# # Store the results
message("Storing results...")
obj <- list(met = met, anno_region = anno_region, annos = annos, rna = rna,
            io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "unfiltered/", io$anno_name, ".rds"))
