#' @title Binarise CpG sites
#'
#' @description Script for binarising CpG sites and formatting the coverage file
#'   so it can be directly used from the BPRMeth package. The format of each
#'   file is the following: <chr> <start> <met_level>, where met_level can be
#'   either 0 or 1. To read compressed files, e.g ending in .gz or .bz2, the
#'   R.utils package needs to be installed.
#' @param indir Directory containing the coverage files, output from Bismark.
#' @param outdir Directory to store the output files for each cell with exactly
#'   the same name. If NULL, then a directory called `binarised` inside `indir`
#'   will be create by default.
#' @param format Integer, denoting the format of coverage file. When set to `1`,
#'   the coverage file format is assumed to be: "<chr> <start> <end> <met_prcg>
#'   <met_reads> <unmet_reads>". When set to `2`, then the format is assumed to
#'   be: "<chr> <start> <met_prcg> <met_reads> <unmet_reads>".
#' @param no_cores Number of cores to use for parallel processing. If NULL, no
#'   parallel processing is used.
#'
#' @return No value is returned, the binarised data are stored in the outdir.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' \dontrun{
#' # Met directory
#' met_dir <- "name_of_met_dir"
#'
#' binarise_files(met_dir)
#' }
#'
#' @seealso \code{\link{create_melissa_data_obj}}, \code{\link{melissa}},
#'   \code{\link{filter_regions}}
#'
#' @export
#'
binarise_files <- function(indir, outdir = NULL, format = 1, no_cores = NULL) {
  # Whether or not to run on parallel mode
  is_parallel <- TRUE
  if (is.null(no_cores)) {
    is_parallel <- FALSE
    no_cores <- 1
  }
  # The out directory will be inside `indir/binarised`
  if (is.null(outdir)) {
    outdir <- paste0(indir, "/binarised/")
  }

  # Load cell filenames
  filenames <- setdiff(list.files(path = indir),
                       list.dirs(path = indir, recursive = FALSE,
                                 full.names = FALSE))

  # Create out directory if it doesn't exist
  ifelse(!dir.exists(outdir), dir.create(outdir), FALSE)

  i <- 0 # FOR CMD check to pass
  # Parallelise processing
  if (is_parallel) {
    doParallel::registerDoParallel(no_cores = no_cores)
    invisible(foreach::foreach(i = seq_along(filenames)) %dopar% {
      # Process each file
      .process_bismark_file(filename = filenames[i], outdir = outdir,
                            format = format)
    })
    doParallel::stopImplicitCluster()
  }else {
    for (i in seq_along(filenames)) {
      # Process each file
      .process_bismark_file(filename = filenames[i], indir = indir,
                            outdir = outdir, format = format)
    }
  }
}


# Private function for reading and processing a coverage bismark file
.process_bismark_file <- function(filename, indir, outdir, format) {
  # So we can pass Build without NOTEs
  rate = met_reads = unnmet_reads = chr <- NULL
  cell_id <- sub(".gz|.zip|.bz2|.rar|.7z","", filename)
  outfile <- sprintf("%s/%s", outdir, cell_id)
  if ( file.exists(paste0(outfile)) ||
       file.exists(paste0(outfile, ".gz")) ||
       file.exists(paste0(outfile, ".zip")) ||
       file.exists(paste0(outfile, ".bz2")) ||
       file.exists(paste0(outfile, ".rar")) ||
       file.exists(paste0(outfile, ".7z"))) {
    cat(sprintf("Sample %s already processed, skipping...\n", cell_id))
  } else {
    cat(sprintf("Processing %s...\n", cell_id))
    # Load data
    data <- data.table::fread(file = sprintf("%s/%s", indir, filename),
                              verbose = FALSE, showProgress = FALSE)

    # Input format 1
    if (format == 1) {
      if (NCOL(data) != 6) {
        stop("Wrong file format, it should contain 6 columns!")
      }
      colnames(data) <- c("chr","pos", "pos_end", "met_prcg",
                          "met_reads","unnmet_reads")
      data[,rate := round((met_reads/(met_reads + unnmet_reads)))] %>%
        .[,c("pos_end", "met_prcg", "met_reads","unnmet_reads") := NULL] %>%
        .[, chr := as.factor(sub("chr", "", chr))] %>%
        data.table::setkey(chr, pos)
    } else if (format == 2) { # Input format 2
      if (NCOL(data) != 5) {
        stop("Wrong file format, it should contain 5 columns!")
      }
      colnames(data) <- c("chr","pos", "met_prcg", "met_reads","unnmet_reads")
      data[,rate := round((met_reads/(met_reads + unnmet_reads)))] %>%
        .[,c("met_prcg","met_reads","unnmet_reads") := NULL] %>%
        .[, chr := as.factor(sub("chr", "", chr))] %>%
        data.table::setkey(chr, pos)
    } else {
      stop("File format currently not supported!")
    }

    # Sanity check
    tmp <- sum((max(data$rate) > 1) | (min(data$rate) < 0))
    if (tmp > 0) {
      cat(sprintf("%s: There are %d CpG sites that have
                  methylation rate > 1 or < 0\n", cell_id, tmp))
    }
    # Calculate binary methylation status
    cat(sprintf("%s: There are %0.03f%% of sites with non-binary methylation
                rate\n", cell_id, mean(!data$rate %in% c(0,1))))
    # Save results
    data.table::fwrite(data, file = outfile, showProgress = FALSE,
                       verbose = FALSE, col.names = FALSE, sep = "\t")
    # Removing gz step, since it is not portable. The user needs to gzip files
    # after running the binarise function.
    # system(sprintf("gzip -f %s", outfile))
  }
}


#' @title Create methylation regions for all cells
#'
#' @description Wrapper function for creating methylation regions for all cells,
#'   which is the input object for Melissa prior to filtering.
#'
#' @param met_dir Directory of (binarised) methylation files, each file
#'   corresponds to a single cell.
#' @param anno_file The annotation file with `tab` delimited format:
#'   "chromosome", "start", "end", "strand", "id", "name" (optional). Read the
#'   `BPRMeth` documentation for more details.
#' @param chrom_size_file Optional file name to read genome chromosome sizes.
#' @param chr_discarded Optional vector with chromosomes to be discarded.
#' @param is_centre Logical, whether 'start' and 'end' locations are
#'   pre-centred. If TRUE, the mean of the locations will be chosen as centre.
#'   If FALSE, the 'start' will be chosen as the center; e.g. for genes the
#'   'start' denotes the TSS and we use this as centre to obtain K-bp upstream
#'   and downstream of TSS.
#' @param is_window Whether to consider a predefined window region around
#'   centre. If TRUE, then 'upstream' and 'downstream' parameters are used,
#'   otherwise we consider the whole region from start to end location.
#' @param upstream Integer defining the length of bp upstream of 'centre' for
#'   creating the genomic region. If is_window = FALSE, this parameter is
#'   ignored.
#' @param downstream Integer defining the length of bp downstream of 'centre'
#'   for creating the genomic region. If is_window = FALSE, this parameter is
#'   ignored.
#' @param cov Integer defining the minimum coverage of CpGs that each region
#'   must contain.
#' @param sd_thresh Optional numeric defining the minimum standard deviation of
#'   the methylation change in a region. This is used to filter regions with no
#'   methylation variability.
#' @param no_cores Number of cores to be used for parallel processing of data.
#'
#' @return A \code{melissa_data_obj} object, with the following elements:
#'   \itemize{ \item{ \code{met}: A list of elements of length N, where N are
#'   the total number of cells. Each element in the list contains another list
#'   of length M, where M is the total number of genomic regions, e.g.
#'   promoters. Each element in the inner list is an \code{I X 2} matrix, where
#'   I are the total number of observations. The first column contains the input
#'   observations x (i.e. CpG locations) and the 2nd column contains the
#'   corresponding methylation level.} \item {\code{anno_region}: The annotation
#'   object.} \item {\code{opts}: A list with the parameters that were used for
#'   creating the object. } }
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @examples
#' \dontrun{
#' # Met directory
#' met_dir <- "name_of_met_dir"
#' # Annotation file name
#' anno_file <- "name_of_anno_file"
#'
#' obj <- create_melissa_data_obj(met_dir, anno_file)
#'
#' # Extract annotation regions
#' met <- obj$met
#'
#' # Extract annotation regions
#' anno <- obj$anno_region
#' }
#'
#' @seealso \code{\link{binarise_files}}, \code{\link{melissa}},
#'   \code{\link{filter_regions}}
#'
#' @export
#'
create_melissa_data_obj <- function(met_dir, anno_file, chrom_size_file = NULL,
    chr_discarded = NULL, is_centre = FALSE, is_window = TRUE, upstream = -5000,
    downstream = 5000, cov = 5, sd_thresh = -1, no_cores = NULL) {

  # Parameter options
  opts <- list()
  # Load cell filenames
  opts$met_files <- setdiff(list.files(path = met_dir, full.names = FALSE),
                            list.dirs(path = met_dir, recursive = FALSE,
                                      full.names = FALSE))
  opts$cell_id <- unname(sapply(opts$met_files, function(filename)
    sub(".gz|.zip|.bz2|.rar|.7z", "", filename) ) )
  # opts$met_files <- list.files(met_dir, pattern = "*.gz", full.names = FALSE)
  # opts$cell_id <- opts$met_files
  # opts$cell_id <- sapply(strsplit(opts$met_files,".", fixed = TRUE),`[`,1)
  opts$is_centre  <- is_centre   # Whether genomic region is already pre-centred
  opts$is_window  <- is_window   # Use predefined window region
  opts$upstream   <- upstream    # Upstream of centre
  opts$downstream <- downstream  # Downstream of centre
  opts$chrom_size <- chrom_size_file  # Chromosome size file
  opts$chr_discarded <- chr_discarded # Chromosomes to discard
  opts$cov        <- cov         # Regions with at least n CpGs
  opts$sd_thresh  <- sd_thresh   # Variance of methylation within region

  # Read annotation file and create annotation regions
  anno_region <- BPRMeth::read_anno(file = anno_file,
        chrom_size_file = opts$chrom_size, chr_discarded = opts$chr_discarded,
        is_centre = opts$is_centre, is_window = opts$is_window,
        upstream = opts$upstream, downstream = opts$downstream,
        is_anno_region = TRUE, delimiter = "\t")

  # Create methylation regions
  if (is.null(no_cores)) {
    met <- lapply(X = opts$met_files, FUN = function(n){
      # Read scBS seq data
      met_dt <- BPRMeth::read_met(file = sprintf("%s/%s", met_dir, n),
                                  type = "sc_seq", strand_info = FALSE)
      # Create promoter methylation regions
      res <- BPRMeth::create_region_object(met_dt = met_dt, anno_dt = anno_region,
                  cov = opts$cov, sd_thresh = opts$sd_thresh,
                  ignore_strand = TRUE, filter_empty_region = FALSE)$met
      names(res) <- NULL
      return(res)
    })
  } else{
    met <- parallel::mclapply(X = opts$met_files, FUN = function(n){
      # Read scBS seq data
      met_dt <- BPRMeth::read_met(file = sprintf("%s/%s", met_dir, n),
                                  type = "sc_seq", strand_info = FALSE)
      # Create promoter methylation regions
      res <- BPRMeth::create_region_object(met_dt = met_dt, anno_dt = anno_region,
                  cov = opts$cov, sd_thresh = opts$sd_thresh,
                  ignore_strand = TRUE, filter_empty_region = FALSE)$met
      names(res) <- NULL
      return(res)
    }, mc.cores = no_cores)
  }

  # Add cell IDs to list
  names(met) <- opts$cell_id
  # Store the object
  obj <- structure(list(met = met, anno_region = anno_region, opts = opts),
                   class = "melissa_data_obj")
  return(obj)
}

