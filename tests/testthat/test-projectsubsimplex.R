test_that("projectsubsimplex() behaves with simulated data", {
  # Plan: start with some data, represent it in one higher dimension,
  # translate point back and forth according to the merge direction,
  # apply projection to get back to original data with scores given by size of translation exactly
  # starting data
  set.seed(34)
  X_5 <- matrix(runif(5*10), nrow = 10)
  X_5 <- divL1(X_5)

  #raise to higher D
  w <- 0.7
  X <- cbind(w * X_5[, 1], (1-w) * X_5[, 1], X_5[, 2:5]) #assume it is first two components to merge

  #translate points eps - make sure translation doesn't take points out of the simplex
  eps <- 1E-3
  tranvec <- (diag(6)[2, ] - diag(6)[1, ]) / sqrt(2)  #and make tranvec Euclidean size of 1
  tranX <- X * NA
  tranX[2 * (1:5), ] <- t(t(X[2 * (1:5), ]) + eps*tranvec)
  tranX[2 * (1:5) - 1, ] <- t(t(X[2 * (1:5) - 1, ]) - eps*tranvec)
  stopifnot(all(rowSums(tranX) == 1))
  dists <- sqrt(rowSums((X - tranX)^2))

  proj <- projectsubsimplex(X = tranX, v1 = 1, v2 = 2, w = w)
  expect_equal(proj$mergedirection, c(-1, 1, 0, 0, 0, 0))
  expect_equal(abs(proj$scores), rep(eps, nrow(X)), ignore_attr = TRUE)
  expect_equal(sign(proj$scores), c(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1), ignore_attr = TRUE)
  expect_equal(proj$V[1, ], (w * diag(6)[1, ] + (1-w) * diag(6)[2, ]))
  expect_equal(proj$X, X_5, ignore_attr = TRUE)
})

test_that("projectsubsimplex works", {
  X = matrix(c(0.5, 0, 0.5,
               0.2, 0.6, 0.2,
               0.5, 0.3, 0.2),
             byrow = T, nrow = 3, ncol = 3)
  V = diag(3)
  colnames(X) <- colnames(V) <- c('Var1', 'Var2', 'Var3')
  res = projectsubsimplex(X, V, 2, 3, 0.6)

  X.expect = matrix(c(0.5, 0.5,
                      0.8, 0.2,
                      0.5, 0.5),
                    byrow = T, nrow = 3, ncol = 2)
  colnames(X.expect) = c('(Var2,Var3)', 'Var1')
  expect_equal(res$X, X.expect)

  V.expect = matrix(c(0, 0.6, 0.4,
                      1, 0, 0),
                    byrow = T, nrow = 2, ncol = 3)
  colnames(V.expect) = c('Var1', 'Var2', 'Var3')
  expect_equal(res$V, V.expect)

  expect_equal(res$scores, c(sqrt(2)*0.3, sqrt(2)*(-0.12), 0))

  merge.expect = c(0, -1, 1)
  names(merge.expect) = c('Var1', 'Var2', 'Var3')
  expect_equal(res$mergedirection, merge.expect)
})

test_that("projectsubsimplex works when X is reduced to one dimensional", {
  X = matrix(c(0, 1,
               0.75, 0.25,
               0.6, 0.4),
             byrow = T, nrow = 3, ncol = 2)
  V = diag(2)
  colnames(X) <- colnames(V) <- c('Var1', 'Var2')
  res = projectsubsimplex(X, V, 1, 2, 0.6)

  X.expect = matrix(c(1,
                      1,
                      1),
                    byrow = T, nrow = 3, ncol = 1)
  colnames(X.expect) = '(Var1,Var2)'
  expect_equal(res$X, X.expect)

  V.expect = matrix(c(0.6, 0.4),
                    byrow = T, nrow = 1, ncol = 2)
  colnames(V.expect) = c('Var1', 'Var2')
  expect_equal(res$V, V.expect)

  expect_equal(res$scores, c(sqrt(2)*0.6, sqrt(2)*(-0.15), 0))

  merge.expect = c(-1, 1)
  names(merge.expect) = c('Var1', 'Var2')
  expect_equal(res$mergedirection, merge.expect)


})

