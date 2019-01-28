#' @title \code{Melissa}: Bayesian clustering and imputation of single cell
#'   methylomes
#' @description Bayesian clustering and imputation of single cell methylomes
#' @docType package
#' @name Melissa
#'
#' @return Melissa main package documentation.
#'
#' @author C.A.Kapourani \email{kapouranis.andreas@@gmail.com}
#'
#' @rawNamespace importFrom(magrittr,"%>%")
#' @rawNamespace importFrom(data.table,":=")
#' @import GenomicRanges ggplot2
#' @importFrom stats pnorm dbinom dnorm
#' @importFrom matrixcalc matrix.trace
#' @importFrom cowplot plot_grid
#'
.datatable.aware <- TRUE
NULL
#> NULL


.onLoad <- function(libname = find.package("Melissa"), pkgname = "Melissa"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(
      # sample file names from taxstats
      c(# we use the magrittr pipe
        "."
      )
    )
  invisible()
}
