context("P1")

fit <- mb_reconstruction(p1Seasonal, pc3seasons, 1750, log.trans = 1:3)
set.seed(24)
cvFolds <- make_Z(1922:2003, nRuns = 50, frac = 0.25, contiguous = TRUE)
cv <- cv_mb(p1Seasonal, pc3seasons, cvFolds, 1750, log.trans = 1:3, return.type = 'metric means')

test_that("P.1 reconstruction produces data.table output", {
  expect_is(fit, 'data.table')
})

test_that("P.1 reconstruction is numerically correct", {
  expect_equal(fit$Q[1], 596.9214, tolerance = 1e-4)
})

test_that("P.1 cross-validation produces data.table output", {
  expect_is(cv, 'data.table')
})

test_that("P.1 cross-validation is numerically correct", {
  expect_equal(cv$R2[1], 0.4588, tolerance = 1e-4)
})
