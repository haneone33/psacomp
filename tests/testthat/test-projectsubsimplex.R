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
