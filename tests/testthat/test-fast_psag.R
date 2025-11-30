test_that("fast PSAG matches to classical PSAG", {
  set.seed(1)
  X = matrix(runif(40),8,5)
  X = sweep(X, 1, rowSums(X), '/')
  res.fast = fast_psag(X, n_components = 4)
  res = psag(X)
  expect_equal(res.fast$X, res$X)
  expect_equal(res.fast$Vhat, res$Vhat)
  expect_equal(res.fast$Xhat, res$Xhat)
  expect_equal(res.fast$Xhat_reduced, res$Xhat_reduced)
  expect_equal(res.fast$scores, res$scores)
  expect_equal(res.fast$RSS, res$RSS)
  expect_equal(res.fast$backwards_mean, res$backwards_mean)
  expect_equal(res.fast$loadings, res$loadings)
  expect_equal(res.fast$construction_info, res$construction_info)
})


test_that("fast PSAG matches to classical PSAG", {
  set.seed(1)
  X = cbind(matrix(runif(200,0,1),20,10),
            matrix(runif(200,0,0.01),20,10))
  X = sweep(X, 1, rowSums(X), '/')
  res.fast = fast_psag(X, n_components = 19)
  res = psag(X)

  expect_equal(res.fast$X, res$X)
  expect_equal(res.fast$Vhat, res$Vhat)
  expect_equal(res.fast$Xhat, res$Xhat)
  expect_equal(res.fast$Xhat_reduced, res$Xhat_reduced)
  expect_equal(res.fast$scores, res$scores)
  expect_equal(res.fast$RSS, res$RSS)
  expect_equal(res.fast$backwards_mean, res$backwards_mean)
  expect_equal(res.fast$loadings, res$loadings)
  expect_equal(res.fast$construction_info, res$construction_info)
})

test_that("fast PSAG matches to classical PSAG", {
  set.seed(1)
  X = cbind(matrix(runif(400,0,1),20,20),
            matrix(runif(400,0,0.01),20,20))
  X = sweep(X, 1, rowSums(X), '/')
  res.fast = fast_psag(X, n_components = NULL, percent_appx = 0.999)
  res = psag(X)

  expect_equal(res.fast$scores[,1:10], res$scores[,1:10], tolerance = 1e-3)
  expect_equal(res.fast$RSS, res$RSS, tolerance = 1e-3)
  expect_equal(res.fast$backwards_mean, res$backwards_mean, tolerance = 1e-3)
  expect_equal(res.fast$loadings[,1:10], res$loadings[,1:10], tolerance = 1e-3)
})
