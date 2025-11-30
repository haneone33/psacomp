test_that("PSAG dimension 2 to 1", {
  # setup
  X = t(matrix(c(1/3, 1/3, 1/3,
                 0,   0,   1,
                 0,   1,   0),
               byrow = T, nrow = 3, ncol = 3,
               dimnames = list(NULL, paste0('e',1:3))))
  V = diag(1, 3)
  dimnames(V) = list(paste0('e',1:3), paste0('V',1:3))

  res = psag_tolowerdim(X, V)

  idx_exp = 1
  Vhat_exp = matrix(c(1/11,  1/11,
                      0,     10/11,
                      10/11, 0),
                    byrow = T, nrow = 3, ncol = 2,
                    dimnames = list(paste0('e',1:3), paste0('V',1:2)))
  Xhat_exp = t(matrix(c(1/11, 5/11,  5/11,
                        1/11, 0,     10/11,
                        1/11, 10/11, 0),
                      byrow = T, nrow = 3, ncol = 3,
                      dimnames = list(NULL, paste0('e',1:3))))
  Xhat_reduced_exp = t(matrix(c(0.5, 0.5,
                        1, 0,
                        0, 1),
                      byrow = T, nrow = 3, ncol = 2,
                      dimnames = list(NULL, paste0('V',1:2))))
  residuals_exp = t(matrix(c(8/33,  -4/33, -4/33,
                             -1/11, 0,     1/11,
                             -1/11, 1/11,  0),
                           byrow = T, nrow = 3, ncol = 3,
                           dimnames = list(NULL, paste0('e',1:3))))
  scores_exp = c(1,-1,-1) * sqrt(colSums(residuals_exp^2))
  RSS_exp = sum(scores_exp^2)

  expect_equal(res$idx, idx_exp)
  expect_equal(res$Vhat, Vhat_exp)
  expect_equal(res$Xhat, Xhat_exp)
  expect_equal(res$Xhat_reduced, Xhat_reduced_exp)
  expect_equal(res$residuals, residuals_exp)
  expect_equal(res$scores, scores_exp)
  expect_equal(res$RSS, RSS_exp)
})

test_that("PSAG dimension 1 to 0", {
  X = matrix(c(1/4, 1/4, 1/4,
               3/4, 0,   3/8,
               0,   3/4, 3/8),
             byrow = T, nrow = 3, ncol = 3,
             dimnames = list(NULL, paste0('e',1:3)))
  V = matrix(c(1/4, 1/4,
               3/4, 0,
               0,   3/4),
             byrow = T, nrow = 3, ncol = 2)
  dimnames(V) = list(paste0('e',1:3), paste0('V',1:2))

  res = psag_tozerodim(X, V, type = '2')

  idx_exp = 2
  Vhat_exp = matrix(c(1/4,3/8,3/8), nrow = 3, ncol = 1, dimnames = list(NULL, 'V1'))
  Xhat_exp = matrix(c(1/4,3/8,3/8,
                      1/4,3/8,3/8,
                      1/4,3/8,3/8),byrow = F,3,3)
  Xhat_reduced_exp = matrix(c(1,1,1), nrow = 1, ncol = 3, dimnames = list('V1', NULL))
  residuals_exp = X - Xhat_exp
  scores_exp = c(-1,1,0)*sqrt(colSums(residuals_exp^2))
  RSS_exp = sum(scores_exp^2)

  expect_equal(res$idx, idx_exp)
  expect_equal(res$Vhat, Vhat_exp)
  expect_equal(res$Xhat, Xhat_exp)
  expect_equal(res$Xhat_reduced, Xhat_reduced_exp)
  expect_equal(res$residuals, residuals_exp)
  expect_equal(res$scores, scores_exp)
  expect_equal(res$RSS, RSS_exp)
})



