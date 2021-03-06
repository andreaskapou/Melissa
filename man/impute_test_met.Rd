% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{impute_test_met}
\alias{impute_test_met}
\title{Impute/predict test methylation states}
\usage{
impute_test_met(
  obj,
  test,
  basis = NULL,
  is_predictive = TRUE,
  return_test = FALSE
)
}
\arguments{
\item{obj}{Output of Melissa inference object.}

\item{test}{Test data to evaluate performance.}

\item{basis}{Basis object, if NULL we perform imputation using Melissa,
otherwise using BPRMeth.}

\item{is_predictive}{Logical, use predictive distribution for imputation, or
choose the cluster label with the highest responsibility.}

\item{return_test}{Whether or not to return a list with the predictions.}
}
\value{
A list containing two vectors, the true methylation state and the
  predicted/imputed methylation states.
}
\description{
Make predictions of missing methylation states, i.e. perfrom
  imputation using Melissa. This requires keepin a subset of data as a held
  out test set during Melissa inference. If you want to impute a whole
  directory containing cells (files) with missing methylation levels, see
  \code{\link{impute_met_files}}.
}
\examples{
# Extract synthetic data
dt <- melissa_synth_dt

# Partition to train and test set
dt <- partition_dataset(dt)

# Create basis object from BPRMeth package
basis_obj <- BPRMeth::create_rbf_object(M = 3)

# Run Melissa
melissa_obj <- melissa(X = dt$met, K = 2, basis = basis_obj, vb_max_iter=10,
   vb_init_nstart = 1, is_parallel = FALSE, is_verbose = FALSE)

imputation_obj <- impute_test_met(obj = melissa_obj,
                                   test = dt$met_test)

}
\seealso{
\code{\link{create_melissa_data_obj}}, \code{\link{melissa}},
  \code{\link{filter_regions}}, \code{\link{eval_imputation_performance}},
  \code{\link{eval_cluster_performance}}
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
