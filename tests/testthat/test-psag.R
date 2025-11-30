test_that("2 dimensional case", {
  # setup
  X = matrix(c(1/3, 1/3, 1/3,
               0,   1,   0,
               0,   0,   1),
             byrow = T, nrow = 3, ncol = 3)
  res = psag(X)

  Vhat1_exp = matrix(c(1/11,  1/11,
                       0,     10/11,
                       10/11, 0),
                     byrow = T, nrow = 3, ncol = 2,
                     dimnames = list(paste0('e',1:3), paste0('V',1:2)))
  Xhat1_exp = matrix(c(1/11, 1/11, 1/11,
                       5/11, 10/11, 0,
                       5/11, 0,    10/11),
                     byrow = T, nrow = 3, ncol = 3,
                     dimnames = list(paste0('e',1:3), NULL))
  Xhat_reduced1_exp = matrix(c(0.5, 0, 1,
                               0.5, 1, 0),
                             byrow = T, nrow = 2, ncol = 3,
                             dimnames = list(paste0('V',1:2), NULL))
  residual1_exp = t(X) - Xhat1_exp
  scores1_exp = c(1, -1, -1)*sqrt(colSums(residual1_exp**2))
  RSS1_exp = sum(residual1_exp**2)

  Vhat0_exp = matrix(c(1/11, 5/11, 5/11), nrow = 3, ncol = 1,
                     dimnames = list(paste0('e',1:3), 'V1'))
  Xhat0_exp = matrix(c(1/11, 1/11, 1/11,
                       5/11, 5/11, 5/11,
                       5/11, 5/11, 5/11),
                     byrow = T, nrow = 3, ncol = 3,
                     dimnames = list(paste0('e',1:3), NULL))
  Xhat_reduced0_exp = matrix(c(1, 1, 1), nrow = 1, ncol = 3,
                             dimnames = list('V1', NULL))
  residual0_exp = Xhat1_exp - Xhat0_exp
  scores0_exp = c(0, 1, -1)*sqrt(colSums(residual0_exp**2))
  RSS0_exp = sum(residual0_exp**2)

  expect_equal(res$Vhat$`r=1`, Vhat1_exp)
  expect_equal(res$Xhat$`r=1`, t(Xhat1_exp))
  expect_equal(res$Xhat_reduced$`r=1`, t(Xhat_reduced1_exp))

  expect_equal(res$Vhat$`r=0`, Vhat0_exp)
  expect_equal(res$Xhat$`r=0`, t(Xhat0_exp))
  expect_equal(res$Xhat_reduced$`r=0`, t(Xhat_reduced0_exp))

  backwards_mean_exp = Vhat0_exp
  scores_exp = matrix(c(scores0_exp, scores1_exp),
                      byrow = F, ncol = 2,
                      dimnames = list(NULL, paste0('PC',1:2)))
  RSS_exp = c(RSS0_exp, RSS1_exp)
  names(RSS_exp) = paste0('PC',1:2)
  loadings_exp = matrix(c(0,     10/11,
                          5/11,  -5/11,
                          -5/11, -5/11),
                        byrow = T, ncol = 2,
                        dimnames = list(paste0('e',1:3), paste0('PC',1:2)))
  expect_equal(res$backwards_mean, backwards_mean_exp)
  expect_equal(res$scores, scores_exp)
  expect_equal(res$RSS, RSS_exp)
  expect_equal(res$loadings, loadings_exp)

})

test_that("3 dimensional case", {
  # setup
  X = t(matrix(c(0.3, 0,   0,
                 0.1, 0.2, 0.1,
                 0.1, 0.1, 0.2,
                 0.5, 0.7, 0.7),
               byrow = T, nrow = 4, ncol = 3))
  res = psag(X)

  Vhat2_exp = matrix(c(0.9, 0,   0,
                       0,   0.3, 0,
                       0,   0,   0.3,
                       0.1, 0.7, 0.7),
                     byrow = T, nrow = 4, ncol = 3,
                     dimnames = list(paste0('e',1:4), paste0('V',1:3)))
  Xhat2_exp = matrix(c(0.3, 0,   0,
                       0.1, 0.2, 0.1,
                       0.1, 0.1, 0.2,
                       0.5, 0.7, 0.7),
                     byrow = T, nrow = 4, ncol = 3,
                     dimnames = list(paste0('e',1:4), NULL))
  Xhat_reduced2_exp = matrix(c(1/3, 0,   0,
                               1/3, 2/3, 1/3,
                               1/3, 1/3, 2/3),
                             byrow = T, nrow = 3, ncol = 3,
                             dimnames = list(paste0('V',1:3), NULL))
  residual2_exp = t(X) - Xhat2_exp
  scores2_exp = c(0, 0, 0)*sqrt(colSums(residual2_exp**2))
  RSS2_exp = sum(residual2_exp**2)

  Vhat1_exp = matrix(c(0.9, 0,
                       0,   0.15,
                       0,   0.15,
                       0.1, 0.7),
                     byrow = T, nrow = 4, ncol = 2,
                     dimnames = list(paste0('e',1:4), paste0('V',1:2)))
  Xhat1_exp = matrix(c(0.3, 0,    0,
                       0.1, 0.15, 0.15,
                       0.1, 0.15, 0.15,
                       0.5, 0.7,  0.7),
                     byrow = T, nrow = 4, ncol = 3,
                     dimnames = list(paste0('e',1:4), NULL))
  Xhat_reduced1_exp = matrix(c(1/3, 0, 0,
                               2/3, 1, 1),
                             byrow = T, nrow = 2, ncol = 3,
                             dimnames = list(paste0('V',1:2), NULL))
  residual1_exp = t(X) - Xhat1_exp
  scores1_exp = c(0, 1, -1)*sqrt(colSums(residual1_exp**2))
  RSS1_exp = sum(residual1_exp**2)

  Vhat0_exp = matrix(c(0.1, 0.4/3, 0.4/3, 1.9/3), nrow = 4, ncol = 1,
                     dimnames = list(paste0('e',1:4), 'V1'))
  Xhat0_exp = matrix(c(0.1,   0.1,   0.1,
                       0.4/3, 0.4/3, 0.4/3,
                       0.4/3, 0.4/3, 0.4/3,
                       1.9/3, 1.9/3, 1.9/3),
                     byrow = T, nrow = 4, ncol = 3,
                     dimnames = list(paste0('e',1:4), NULL))
  Xhat_reduced0_exp = matrix(c(1, 1, 1), nrow = 1, ncol = 3,
                             dimnames = list('V1', NULL))
  residual0_exp = Xhat1_exp - Xhat0_exp
  scores0_exp = c(-1, 1, 1)*sqrt(colSums(residual0_exp**2))
  RSS0_exp = sum(residual0_exp**2)

  expect_equal(res$Vhat$`r=2`, Vhat2_exp)
  expect_equal(res$Xhat$`r=2`, t(Xhat2_exp))
  expect_equal(res$Xhat_reduced$`r=2`, t(Xhat_reduced2_exp))

  expect_equal(res$Vhat$`r=1`, Vhat1_exp)
  expect_equal(res$Xhat$`r=1`, t(Xhat1_exp))
  expect_equal(res$Xhat_reduced$`r=1`, t(Xhat_reduced1_exp))

  expect_equal(res$Vhat$`r=0`, Vhat0_exp)
  expect_equal(res$Xhat$`r=0`, t(Xhat0_exp))
  expect_equal(res$Xhat_reduced$`r=0`, t(Xhat_reduced0_exp))

  backwards_mean_exp = Vhat0_exp
  RSS_exp = c(RSS0_exp, RSS1_exp, RSS2_exp)
  names(RSS_exp) = paste0('PC',1:3)

  if(res$construction_info$`r=2`$idx == 2){
    scores_exp = matrix(c(scores0_exp, scores1_exp, scores2_exp),
                        byrow = F, ncol = 3,
                        dimnames = list(NULL, paste0('PC',1:3)))
    loadings_exp = cbind(Vhat1_exp[,2]-backwards_mean_exp,
                         Vhat2_exp[,2]-backwards_mean_exp,
                         c(0,0,0,1)-backwards_mean_exp)
  }else if(res$construction_info$`r=2`$idx == 3){
    scores_exp = matrix(c(scores0_exp, -scores1_exp, scores2_exp),
                        byrow = F, ncol = 3,
                        dimnames = list(NULL, paste0('PC',1:3)))
    loadings_exp = cbind(Vhat1_exp[,2]-backwards_mean_exp,
                         Vhat2_exp[,3]-backwards_mean_exp,
                         c(0,0,0,1)-backwards_mean_exp)
  }
  dimnames(loadings_exp) = list(paste0('e',1:4), paste0('PC',1:3))
  expect_equal(res$backwards_mean, backwards_mean_exp)
  expect_equal(res$scores, scores_exp)
  expect_equal(res$RSS, RSS_exp)
  expect_equal(res$loadings, loadings_exp)
})


test_that("test signs of scores", {
  X = matrix(c(0.2, 0.7, 0.1,
               0.2, 0.6, 0.2,
               0.1, 0.8, 0.1,
               0.1, 0.7, 0.2,
               0.2, 0.2, 0.6,
               0.2, 0.3, 0.5,
               0.1, 0.2, 0.7,
               0.1, 0.3, 0.6),
             byrow = T, ncol = 3)
  res = psag(X)
  # second mode is always e1 as positive
  # first mode is e2 vs e3, need to check which is positive
  # sign of scores when e2 is positive
  sign_exp = matrix(c(1, 1,
                      1, 1,
                      1, -1,
                      1, -1,
                      -1, 1,
                      -1, 1,
                      -1, -1,
                      -1, -1), byrow = T, ncol = 2,
                    dimnames = list(NULL, paste0('PC',1:2)))
  if(res$loadings['e3','PC1']>0){
    sign_exp[,2] = -sign_exp[,2]
  }
  expect_equal(sign(res$scores), sign_exp)
})
