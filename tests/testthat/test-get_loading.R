test_that("test loading vectors for PSA-S", {
  X = matrix(c(0, 0, 1,
               3/8, 1/8, 1/2,
               3/4, 1/4, 0),
             byrow = T, nrow = 3, ncol = 3)
  out = psas(X, testweights = seq(0, 1, length.out = 100))
  loadings = get_loading(out)
  loadings2 = c(1, -1, 0)
  loadings1 = c(3/4, 1/4, -1)
  expect_true(base::min(sum((loadings[,2] - loadings2)^2), sum((loadings[,2] + loadings2)^2)) < 1e-4)
  expect_true(base::min(sum((loadings[,1] - loadings1)^2), sum((loadings[,1] + loadings1)^2)) < 1e-4)
})


test_that("test loading vectors for PSA-O", {
  x1 = c(0, 0, 1)
  x2 = c(cos(pi/6), sin(pi/6), 0)
  x3 = x2/sqrt(2) + c(0, 0, 1/sqrt(2))

  X = rbind(x1, x2, x3)
  X = sweep(X, 1, rowSums(X), '/')
  out = psao(X, testweights = seq(0, 1, length.out = 100))
  loadings = get_loading(out)
  loadings2 = c(1, -1, 0)
  loadings1 = x2 - c(0,0,1)
  expect_true(base::min(sum((loadings[,2] - loadings2)^2), sum((loadings[,2] + loadings2)^2)) < 1e-4)
  expect_true(base::min(sum((loadings[,1] - loadings1)^2), sum((loadings[,1] + loadings1)^2)) < 1e-4)
})
