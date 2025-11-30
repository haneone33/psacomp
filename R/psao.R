#' @title Principle Sub Orthant Analysis for the Simplex
#' @inheritParams recurse_tolowerdimension
#' @noRd
psao <- function(X, V = diag(ncol(X)), testweights = seq(0, 1, length.out = 100),
                 merges = NULL, maxsteps = min(ncol(X), nrow(merges)) - 1){
  results <- recurse_tolowerdimension(X, V = V,
                                       evaluator = evaluatesuborthant,
                                       projector = projectsuborthant,
                                       testweights = testweights,
                                       merges = merges,
                                       maxsteps = maxsteps
  )

  results$X = X
  results$Xhat = list()
  for(i in names(results$pts)){
    results$Xhat[[i]] = divL1(divL2(results$pts[[i]]) %*% results$vertices[[i]])
  }

  return(results)
}
