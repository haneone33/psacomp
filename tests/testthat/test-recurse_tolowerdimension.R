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

test_that("recurse_tolowerdimension behaves on toy evaluator and projector", {
  set.seed(3)
  X <- matrix(runif(5 * 20), nrow = 20)
  X <- X/rowSums(X)
  lowerD2 <- suppressMessages(recurse_tolowerdimension(X, evaluator = toyevaluator, projector = toyprojector, testweights = seq(0, 1, by = 0.1), maxsteps = 2))
  expect_equal(lowerD2$Vbirth$`r=4`, c(2, -1, -4, -5), ignore_attr = TRUE) #2 stands for reduced
  expect_equal(lowerD2$Vbirth$`r=3`, c(3, 2, -5), ignore_attr = TRUE)
  expect_equal(lowerD2$vertices$`r=3`[1, ], c(0.2,0,0,0.8,0), ignore_attr = TRUE)
  expect_equal(lowerD2$vertices$`r=3`[2, ], c(0,0.2,0.8,0,0), ignore_attr = TRUE)
  expect_equal(lowerD2$scores$`r=3`, 1:nrow(X), ignore_attr = TRUE)
  expect_equal(lowerD2$modes$`r=3`[1, ],c(1, 0, 0, -1, 0), ignore_attr = TRUE)
  expect_equal(lowerD2$pts$`r=3`, cbind(X[, 1] + X[, 4], X[, 2] + X[, 3],  X[, 5]), ignore_attr = TRUE)

  rep <- suppressMessages(recurse_tolowerdimension(X, evaluator = NULL, projector = toyprojector,
                                                   merges = lowerD2$merges))
  expect_equal(lowerD2, rep)
})

# test_that("recurse_tolowerdimension gives good names", {
#   set.seed(3)
#   X <- matrix(runif(5 * 20), nrow = 20)
#   X <- X/rowSums(X)
#   colnames(X) <- c("A", "B", "C", "D", "E")
#   res <- suppressMessages(recurse_tolowerdimension(X, evaluator = toyevaluator, projector = toyprojector, testweights = seq(0, 1, by = 0.1), maxsteps = 2))
#   expect_equal(names(res$Vbirth[[1]]), colnames(X))
# })
