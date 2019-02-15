test_that("log sum exp function works", {
    x <- c(0.001, 0.5, 2, 1.4, 1.5)
    expect_gt(log_sum_exp(x), 2.92)
    expect_lt(log_sum_exp(x), 2.93)
})
