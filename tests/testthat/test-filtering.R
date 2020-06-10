context("Filtering correctly genomic regions")

test_that("Filter by CpG coverage works", {
  obj <- melissa_encode_dt
  filt_obj <- filter_by_cpg_coverage(obj, min_cpgcov = 20)
  expect_equal(sum(sapply(filt_obj$met[[1]], function(n) is.na(n[[1]]))), 31)
  expect_equal(filt_obj$met[[1]][[94]], NA)
  expect_error(filter_by_cpg_coverage(obj, min_cpgcov = -1))

  filt_obj <- filter_by_cpg_coverage(obj, min_cpgcov = 1000)
  expect_equal(sum(sapply(filt_obj$met[[1]], function(n) is.na(n[[1]]))),
               NROW(obj$met[[1]]))
})

test_that("Filter by coverage across cells works", {
  obj <- melissa_encode_dt
  filt_obj <- filter_by_coverage_across_cells(obj, min_cell_cov_prcg = 0.2)
  expect_equal(sum(sapply(filt_obj$met[[1]], function(n) is.na(n[[1]]))), 0)
  expect_error(filter_by_coverage_across_cells(obj, min_cell_cov_prcg = -1))
  expect_error(filter_by_coverage_across_cells(obj, min_cell_cov_prcg = 2))
})

test_that("Filter by variability works", {
  obj <- melissa_encode_dt
  filt_obj <- filter_by_variability(obj, min_var = 0.1)
  expect_equal(length(filt_obj$met[[1]]), 99)

  filt_obj <- filter_by_variability(obj, min_var = 0.15)
  expect_equal(length(filt_obj$met[[1]]), 3)

  expect_error(filter_by_variability(obj, min_var = -1))

  filt_obj <- filter_by_variability(obj, min_var = 1000)
  expect_equal(length(filt_obj$met[[1]]), 0)
})

test_that("log sum exp function works", {
  x <- c(0.001, 0.5, 2, 1.4, 1.5)
  expect_gt(log_sum_exp(x), 2.92)
  expect_lt(log_sum_exp(x), 2.93)
})
