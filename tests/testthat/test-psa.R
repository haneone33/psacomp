test_that('check classes for PSA-S output', {
  nr = 3
  nc = 4
  X = matrix(1:(nr*nc), byrow = T, nrow = nr, ncol = nc)
  X = sweep(X, 1, rowSums(X), '/')
  colnames(X) = paste0('Var',1:4)
  res = suppressMessages(psa('s', X))

  expect_type(res$vertices, 'list')
  expect_type(res$pts, 'list')
  expect_type(res$pts.approx, 'list')
  expect_type(res$modes, 'list')
  expect_true(is.matrix(res$scores))
  expect_true(is.matrix(res$loadings))
  expect_true(is.numeric(res$rss))
  expect_true(is.data.frame(res$const.info))

})

test_that("check dimensions for PSA-S output", {
  nr = 3
  nc = 4
  X = matrix(1:(nr*nc), byrow = T, nrow = nr, ncol = nc)
  X = sweep(X, 1, rowSums(X), '/')
  colnames(X) = paste0('Var',1:4)
  res = suppressMessages(psa('s', X))

  expect_equal(names(res$vertices), paste0('r=',1:nc))
  expect_equal(names(res$pts), paste0('r=',1:nc))
  expect_equal(names(res$pts.approx), paste0('r=',1:nc))
  expect_equal(names(res$modes), paste0('r=',1:nc))
  expect_equal(colnames(res$scores), paste0('Comp.',1:nc))
  expect_equal(colnames(res$loadings), paste0('Comp.',1:nc))
  expect_equal(names(res$rss), paste0('Comp.',1:nc))
  expect_equal(rownames(res$const.info), paste0('r=',nc:1))

})

test_that('check classes for PSA-O output', {
  nr = 3
  nc = 4
  X = matrix(1:(nr*nc), byrow = T, nrow = nr, ncol = nc)
  X = sweep(X, 1, rowSums(X), '/')
  colnames(X) = paste0('Var',1:4)
  res = suppressMessages(psa('o', X))

  expect_type(res$vertices, 'list')
  expect_type(res$pts, 'list')
  expect_type(res$pts.approx, 'list')
  expect_type(res$modes, 'list')
  expect_true(is.matrix(res$scores))
  expect_true(is.matrix(res$loadings))
  expect_true(is.numeric(res$rss))
  expect_true(is.data.frame(res$const.info))

})

test_that("check dimensions for PSA-O output", {
  nr = 3
  nc = 4
  X = matrix(1:(nr*nc), byrow = T, nrow = nr, ncol = nc)
  X = sweep(X, 1, rowSums(X), '/')
  colnames(X) = paste0('Var',1:4)
  res = suppressMessages(psa('o', X))

  expect_equal(names(res$vertices), paste0('r=',1:nc))
  expect_equal(names(res$pts), paste0('r=',1:nc))
  expect_equal(names(res$pts.approx), paste0('r=',1:nc))
  expect_equal(names(res$modes), paste0('r=',1:nc))
  expect_equal(colnames(res$scores), paste0('Comp.',1:nc))
  expect_equal(colnames(res$loadings), paste0('Comp.',1:nc))
  expect_equal(names(res$rss), paste0('Comp.',1:nc))
  expect_equal(rownames(res$const.info), paste0('r=',nc:1))

})
