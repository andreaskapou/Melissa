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
io <- list(dataset = "mt-seq")
io$base_dir   <- paste0("../local-data/", io$dataset)
io$out_dir    <- paste0(io$base_dir, "/met/processed/")
io$annos_file <- paste0(io$base_dir, "/annotations/Mus_musculus.GRCm38.75.protein.coding.bed")
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
# # Read annotation file and create a GRanges object
annos <- fread(input = io$annos_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
               showProgress = FALSE, select = c(1,2,3,4,6,10)) %>%
    setnames(c("chr", "start", "end", "id", "strand", "gene_name")) %>%
    .[, gene_name := gsub(".* gene_name \"([^;]+)\";.*", "\\1", gene_name)] %>%
    .[, chr := as.factor(paste0("chr", chr))] %>% subset(id %in% rna$id) %>%
    setkey(chr, start, end) %>% GRanges()

#
# # Create promoter regions
anno_region <- create_anno_region(anno = annos, chrom_size = opts$chrom_size,
                                  upstream = opts$upstream, downstream = opts$downstream)
#
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
saveRDS(obj, file = paste0(io$out_dir, "unfiltered/prom5k.rds"))
