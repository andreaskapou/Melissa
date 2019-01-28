# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) { cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir) }
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rhdf5))

create_region_object_local <- function(met_dt, anno_dt, cov = 5, sd_thresh = 1e-1,
                                       ignore_strand = TRUE,
                                       filter_empty_region = TRUE,
                                       fmin = -1, fmax = 1){
    message("Creating methylation regions ...")
    assertthat::assert_that(methods::is(met_dt, "GRanges"))
    assertthat::assert_that(methods::is(anno_dt, "GRanges"))
    # Find overlaps between met and anno
    overlaps <- GenomicRanges::findOverlaps(query = anno_dt, subject = met_dt,
                                            ignore.strand = ignore_strand)
    if (length(overlaps) < 2) { stop("Low overlap between anno and met data!")}

    # Convert data in vector format for efficiency
    query_hits <- S4Vectors::queryHits(overlaps)
    subj_hits  <- S4Vectors::subjectHits(overlaps)
    gen_ind    <- unique(query_hits) # Indices of genomic locations
    centre     <- anno_dt$centre     # Central locations
    id         <- anno_dt$id         # (Ensembl) IDs
    strand     <- as.character(GenomicRanges::strand(anno_dt))
    cpg_loc    <- GenomicRanges::ranges(met_dt)@start  # CpG locations
    met        <- met_dt$met         # Methylated read
    # Number of columns for each matrix
    if (is.null(met_dt$total)) { D <- 2 } else {total <- met_dt$total; D <- 3 }

    # Extract upstream and downstream lengths
    N <- NROW(anno_dt)
    up_anno <- vector(mode = "integer", N)    # Start location in chromosome
    down_anno <- vector(mode = "integer", N)  # End location in chromosome
    anno_start <- GenomicRanges::ranges(anno_dt)@start # Start info
    anno_end <- anno_start + GenomicRanges::ranges(anno_dt)@width - 1 # End info
    for (i in 1:N) {
        # Depending on the strand we change regions up or downstream of centre
        if (identical(strand[i], "-")) {
            up_anno[i] <- centre[i] - anno_end[i]
            down_anno[i] <- abs(anno_start[i] - centre[i])
        }else {
            up_anno[i] <- anno_start[i] - centre[i]
            down_anno[i] <- anno_end[i] - centre[i]
        }
    }
    rm(N, anno_start, anno_end)

    # Initialize variables -----------------------------------------
    met_region <- list()
    # Initilize list to NA
    for (j in 1:length(id)) { met_region[[id[j]]] <- NA }
    LABEL     <- FALSE                     # Flag variable
    reg_count <- 0                         # Promoter counter
    cpg_ind   <- vector(mode = "integer")  # Vector of CpG indices
    cpg_ind   <- c(cpg_ind, subj_hits[1])  # Add the first subject hit

    total_iter <- NROW(query_hits)
    for (i in 2:total_iter) {
        # If query hits is the same as the previous one
        if (query_hits[i] == query_hits[i - 1]) {
            cpg_ind <- c(cpg_ind, subj_hits[i])  # Add subject hit
            # In case we have the last region
            if (i == total_iter) { reg_count <- reg_count + 1; LABEL <- TRUE }
        }else {reg_count <- reg_count + 1; LABEL <- TRUE } # Increase counter

        if (LABEL) {
            # Central location for region 'reg_count'
            idx <- id[gen_ind[reg_count]]
            # Only keep regions that have more than 'n' CpGs
            if (length(cpg_ind) > cov) {
                # If sd of the methylation level is above threshold
                #if (D == 2) { obs_var <- stats::sd(met[cpg_ind])
                #}else {obs_var <- stats::sd(met[cpg_ind]/total[cpg_ind]) }
                obs_var <- 5
                if (obs_var > sd_thresh) {
                    # Locations of CpGs in the genome
                    region <- cpg_loc[cpg_ind]
                    # Middle location for region 'reg_count'
                    middle <- centre[gen_ind[reg_count]]
                    # Extract strand information, i.e. direction
                    strand_direction <- strand[gen_ind[reg_count]]
                    # Extract upstream information
                    upstream <- up_anno[gen_ind[reg_count]]
                    # Extract downstream information
                    downstream <- down_anno[gen_ind[reg_count]]
                    # Shift CpG locations relative to TSS
                    center_data  <- BPRMeth:::.do_centre_loc(region = region,
                                                             centre = middle, strand_direction = strand_direction)
                    # In the "-" strand the order of the locations should change
                    Order <- base::order(center_data)
                    met_region[[idx]] <- matrix(0, length(cpg_ind), D)
                    # Store normalized locations of methylated CpGs
                    met_region[[idx]][, 1] <- round(BPRMeth:::.minmax_scaling(
                        data = center_data[Order], xmin = upstream,
                        xmax = downstream, fmin = fmin, fmax = fmax), 4)
                    # Store methylated reads in the corresponding locations
                    if (D == 2) { met_region[[idx]][, 2] <- met[cpg_ind][Order]
                    }else{
                        # Store total reads in the corresponding locations
                        met_region[[idx]][, 2] <- total[cpg_ind][Order]
                        # Store methylated reads in the corresponding locations
                        met_region[[idx]][, 3] <- met[cpg_ind][Order]
                    }
                }
            }
            LABEL   <- FALSE
            cpg_ind <- vector(mode = "integer")
            cpg_ind <- c(cpg_ind, subj_hits[i])
        }
    }
    # Keep only covered genomic regions
    if (filter_empty_region) {
        cov_ind <- which(!is.na(met_region))
        met_region <- met_region[cov_ind]
        anno_dt <- anno_dt[cov_ind, ]
    }
    return(structure(list(met = met_region, anno = anno_dt),
                     class = "region_object"))
}


#
# # Data files
io <- list(anno_name = "promoter_hg38")
io$base_dir   <- "../../local-data/encode/"
io$sub_dir    <- "/"
io$out_dir    <- paste0(io$base_dir, "/scWGBS/deepcpg/processed/unfiltered/", io$sub_dir)
io$annos_file <- paste0(io$base_dir, "/annotations/", io$anno_name, ".bed")
io$met_file   <- paste0(io$base_dir, "/scWGBS/deepcpg/raw/", io$sub_dir, "data.h5")

#
# # Parameter options
opts <- list()
opts$is_centre  <- FALSE   # Whether genomic region is already pre-centred
opts$is_window  <- TRUE    # Use predefined window region
opts$upstream   <- -5000   # Upstream of centre
opts$downstream <- 5000    # Downstream of centre
opts$chrom_size <- NULL    # Chromosome size file
opts$cov        <- 3       # Regions with at least n CpGs
opts$sd_thresh  <- -1      # Threshold for variance of methylation across region

# Read annotation data
annos <- fread(io$annos_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
    setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
    .[,c("anno", "chr") := list(NULL, as.factor(sub("chr", "", chr)))] %>%
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

# Extract h5 contents
listing <- h5ls(io$met_file, all = TRUE)
# Find all data nodes
data_nodes <- grep("GM|H1", listing$name)
# Extract full paths
data_paths = paste(listing$group[data_nodes], listing$name[data_nodes], sep = "/")
# Extract chromo
chromo <- h5read(io$met_file, "/chromo")
# Extract position
pos <- h5read(io$met_file, "/pos")

num_cells <- length(data_paths) / 2
met_list <- list()
for (idx in 1:num_cells) {
    cell_name <- sub(".*/", "", data_paths[idx])
    print(cell_name)
    # Extract actual and predicted methylation
    actual <- h5read(io$met_file, data_paths[idx])
    pred <- round(h5read(io$met_file, data_paths[idx + num_cells]), 2)
    # Remove CpGs that we did not know ground truth label
    idx_remove <- which(actual == -1)
    actual <- actual[-idx_remove]
    pred <- pred[-idx_remove]
    chromo_filt <- chromo[-idx_remove]
    pos_filt <- pos[-idx_remove]

    met_dt <- data.table(chr = chromo_filt, start = pos_filt, end = pos_filt, met = actual, total = pred) %>%
        setkey(chr, start, end) %>% GRanges()

    # Create promoter methylation regions
    res <- create_region_object_local(met_dt = met_dt, anno_dt = anno_region,
                                      cov = opts$cov, sd_thresh = opts$sd_thresh,
                                      ignore_strand = TRUE, filter_empty_region = FALSE)$met
    #names(res) <- NULL
    met_list[[cell_name]] <- res
}

#
# # Store the results
message("Storing results...")
obj <- list(met = met_list, anno_region = anno_region, annos = annos, io = io, opts = opts)
saveRDS(obj, file = paste0(io$out_dir, "prom10k.rds"))
