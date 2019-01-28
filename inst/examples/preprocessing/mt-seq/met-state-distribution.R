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
io$base_dir   <- paste0("../../local-data/", io$dataset)
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
for (m_file in io$met_files[1:10]) {
    # Extract cell ID
    cell <- gsub(".*_([^;]+)\\.CpG.*", "\\1", m_file)
    print(cell)
    # Read scBS seq data
    met_dt <- fread(input = sprintf("zcat < %s", m_file), sep = "\t", header = TRUE,
                    stringsAsFactors = TRUE, showProgress = FALSE) %>%
        .[, c("rate") := list(meth_reads / (meth_reads + unmeth_reads))] %>%
        setkey(chr, start)
    # Create promoter methylation regions
    met[[cell]] <- met_dt
}
rm(met_dt)

met_dt <- rbindlist(met)

met_dt$rate <- met_dt$rate * 100
# breaks = 30,

hist(met_dt$rate, main = "Digital output of DNA methylation",
     xlab = "% of methylation", ylab = "Number of CpGs", col = "cornflowerblue")

#' Define ggplot2 theme
#'
gg_theme <- function(){
    p <- theme(
        plot.title = element_text(size = 16,face = 'bold',
                                  margin = margin(0,0,3,0), hjust = 0.5),
        axis.text = element_text(size = rel(1.05), color = 'black'),
        axis.title = element_text(size = rel(1.45), color = 'black'),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1.4, 'lines'),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 12)#,
        #panel.border = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.background = element_blank()
    )
    return(p)
}

#  "#1B9E77" "#D95F02" "#7570B3"
library(ggplot2)
suppressPackageStartupMessages(library(scales))
dna_binary_plot <- ggplot(met_dt, aes(x = rate)) +
    geom_histogram(color = "#7570B3", fill = "#7570B3", bins = 17, alpha=0.75) +
    # scale_y_continuous(trans = 'log2') +
    theme_classic() +
    gg_theme() +
    scale_y_continuous(breaks = pretty_breaks(n = 4)) +
    labs(title = "Digital output of single cell DNA methylation", x = "% of methylation", y = "Number of CpGs")


pdf(file = paste0("met-binary-state.pdf"), width = 8, height = 4, useDingbats = FALSE)
dna_binary_plot
dev.off()

binary_total <- length(which(met_dt$rate == 0.0)) + length(which(met_dt$rate == 100.0))
total <- binary_total / NROW(met_dt)


hemi_total <- length(which(met_dt$rate == 50.0))
total <- hemi_total / NROW(met_dt)
