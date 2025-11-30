test_that("push_vertex_inwards unit simplex", {
  X = diag(1,4)
  dimnames(X) = list(paste0('e',1:4), NULL)
  V = diag(1,4)
  dimnames(V) = list(paste0('e',1:4), paste0('V',1:4))
  idx = 2
  margin = 1e-4
  res = push_vertex_inwards(X, V, idx, margin)
  X_exp = matrix(c(1, margin/3, 0, 0,
                   0, 1-margin, 0, 0,
                   0, margin/3, 1, 0,
                   0, margin/3, 0, 1),
                 byrow = T, nrow = 4, ncol = 4,
                 dimnames = list(paste0('e',1:4), NULL))
  expect_equal(res, X_exp)
})


test_that("push_vertex_inwards subsimplex", {
  X = matrix(c(0.3, 0,   0,   0,
               0.1, 0.3, 0,   0.3,
               0.1, 0,   0.3, 0,
               0.5, 0.7, 0.7, 0.7),
             byrow = T, nrow = 4, ncol = 4,
             dimnames = list(paste0('e',1:4), NULL))
  V = matrix(c(0.9, 0,   0,
               0,   0.3, 0,
               0,   0,   0.3,
               0.1, 0.7, 0.7),
             byrow = T, nrow = 4, ncol = 3,
             dimnames = list(paste0('e',1:4), paste0('V',1:3)))
  idx = 2
  margin = 1e-4
  res = push_vertex_inwards(X, V, idx, margin)
  X_exp = V %*% matrix(c(1/3, margin/2, 0, margin/2,
                         1/3, 1-margin, 0, 1-margin,
                         1/3, margin/2, 1, margin/2),
                       byrow = T, nrow = 3, ncol = 4)
  expect_equal(res, X_exp)
})
