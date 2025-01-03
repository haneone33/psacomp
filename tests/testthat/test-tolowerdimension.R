toyevaluator <- function(X, v1, v2, w = w){
  out <- (w - 0.2)^2  #best w is 0.2
  if (!((v1 == 2) && (v2 == 3))){ #best pair is 2, 3
    out <- out + 1
  }
  return(out)
}

# A function that takes result of `gridsearch()` and applies a projection. Returning `X`, `V`, `scores`, `mergedirection`.
toyprojector <- function(X, V, v1, v2, w){
  newX <- cbind(X[, v1] + X[, v2], X[, -c(v1, v2), drop = FALSE])
  newV <- rbind(w * V[v1, ] + (1-w) * V[v2, ],
                V[-c(v1, v2), , drop = FALSE])
  return(list(
    X = newX,
    V = newV,
    scores = (v1-v2) * (1:nrow(X)),
    mergedirection = V[v2, ] - V[v1, ]
  ))
}

test_that("tolowerdimension behaves on toy evaluator and projector", {

  set.seed(3)
  X <- matrix(runif(5 * 20), nrow = 20)
  X <- X/rowSums(X)
  lowerD <- tolowerdimension(X, evaluator = toyevaluator, projector = toyprojector, testweights = seq(0, 1, by = 0.1))
  expect_equal(lowerD$v1, 3)
  expect_equal(lowerD$v2, 2)
  expect_equal(lowerD$w, 0.8)
  expect_equal(lowerD$V[1, ], c(0,0.2,0.8,0,0))
  expect_equal(lowerD$scores, 1:nrow(X))
  expect_equal(lowerD$mergedirection,c(0,1,-1,0,0))
  expect_equal(lowerD$X, cbind(X[, 2] + X[, 3], X[, -c(2, 3)]))

  lowerDrepeat <- tolowerdimension(X, evaluator = NULL, projector = toyprojector,
                                   merge = data.frame(v1 = lowerD$v1,
                                                      v2 = lowerD$v2,
                                                      w = lowerD$w))
  expect_equal(lowerD, lowerDrepeat)
})
