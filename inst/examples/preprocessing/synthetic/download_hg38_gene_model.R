library(EnsDb.Hsapiens.v86)
library(data.table)
library(purrr)

## https://bioconductor.org/packages/devel/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

## Making a "short cut"
edb <- EnsDb.Hsapiens.v86
## print some informations for this package
edb

edb_protein_coding_genes <- addFilter(edb, TxBiotypeFilter("protein_coding"))


gene_model <- genes(edb_protein_coding_genes)


gene_model <- as.data.table(gene_model)


gene_model <- gene_model %>% .[,c("seqnames", "start", "end", "strand", "gene_id", "gene_name")] %>%
  subset(seqnames %in% c(seq(1:22), "X", "Y"))

fwrite(gene_model, "promoter_hg38.bed", sep = "\t", row.names = FALSE, col.names = FALSE)
