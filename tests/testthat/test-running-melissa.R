context("Checking correct functionality of Melissa")

test_that("Melissa runs and returns correct results", {
  set.seed(1)
  obj <- melissa_synth_dt
  # Partition to train and test set
  expect_error(partition_dataset(obj, data_train_prcg = 1.4))
  expect_error(partition_dataset(obj, data_train_prcg = -1))
  expect_error(partition_dataset(obj, data_train_prcg = -1))
  expect_error(partition_dataset(obj, region_train_prcg = 1.4))
  expect_error(partition_dataset(obj, region_train_prcg = -1))

  obj <- partition_dataset(obj, data_train_prcg = 0.5, region_train_prcg = 0.95,
                           cpg_train_prcg = 0.5, is_synth = FALSE)
  expect_equal(NROW(obj$met[[1]][[1]]), 20)
  expect_equal(NROW(obj$met[[1]][[12]]), 21)

  # Create basis object from BPRMeth package
  basis_obj <- BPRMeth::create_rbf_object(M = 3)
  # Run Melissa VB model
  melissa_obj <- melissa(X = obj$met, K = 2, basis = basis_obj,
                         delta_0 = NULL, alpha_0 = 0.5, vb_max_iter = 5,
                         vb_init_nstart = 1, is_parallel = FALSE,
                         is_verbose = FALSE)

  # Check lots of parameters from Melissa output that match expectations
  expect_gt(melissa_obj$W[1,1,1], 1.13)
  expect_lt(melissa_obj$W[1,1,1], 1.136)

  expect_lt(melissa_obj$r_nk[1,2], 1e-20)

  expect_lt(melissa_obj$delta[1], 14)
  expect_lt(melissa_obj$alpha[1], 61)
  expect_lt(melissa_obj$beta[1], 37.99)
  expect_lt(melissa_obj$lb[2], -3426.022)
  expect_gt(melissa_obj$lb[2], -3426.03)

  # Check imputation performance is correct
  imputation_obj <- impute_met_state(obj = melissa_obj, test = obj$met_test)
  expect_lt(imputation_obj$pred_obs[1], 0.515)
  expect_gt(imputation_obj$pred_obs[1], 0.514)
  expect_lt(melissa_obj$lb[2], -3426.022)

  # Check evaluation of imputation performance is correct
  melissa_obj <- eval_imputation_performance(obj = melissa_obj,
                                             imputation_obj = imputation_obj)

  expect_lt(melissa_obj$imputation$auc, 0.9172)
  expect_gt(melissa_obj$imputation$auc, 0.9171)

  expect_lt(melissa_obj$imputation$f_measure, 0.86812)
  expect_gt(melissa_obj$imputation$f_measure, 0.86811)

  expect_lt(melissa_obj$imputation$pr@y.values[[1]][2], 0.85715)
  expect_gt(melissa_obj$imputation$pr@y.values[[1]][2], 0.85714)


  # Check the clustering performance is correct
  melissa_obj <- eval_cluster_performance(melissa_obj, obj$opts$C_true)

  expect_lt(melissa_obj$clustering$ari, 1.00001)
  expect_gt(melissa_obj$clustering$ari, 0.99999)

  expect_lt(melissa_obj$clustering$error, 0.001)
  expect_gt(melissa_obj$clustering$error, -0.00001)


  ##
  # Test for Melissa Gibbs output
  melissa_gibbs <- melissa_gibbs(X = obj$met, K = 2, basis = basis_obj,
                               gibbs_nsim = 10, gibbs_burn_in = 5,
                               is_parallel = FALSE, is_verbose = FALSE)

  expect_lt(melissa_gibbs$summary$pi[1], 0.668493)
  expect_gt(melissa_gibbs$summary$pi[1], 0.668492)

  expect_lt(melissa_gibbs$r_nk[2,1], 1.001)
  expect_gt(melissa_gibbs$r_nk[2,1], 0.999)

  expect_lt(melissa_gibbs$W[1,1,1], 1.3791)
  expect_gt(melissa_gibbs$W[1,1,1], 1.3790)
})
