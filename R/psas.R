#' @title Principle Sub Orthant Analysis for the Simplex
#' @inheritParams recurse_tolowerdimension
#' @noRd
psas <- function(X, V = diag(ncol(X)), testweights = seq(0, 1, length.out = 100),
                 merges = NULL, maxsteps = min(ncol(X), nrow(merges)) - 1){

  if(!is.null(colnames(X))) colnames(V) <- colnames(X)
  else if(!is.null(colnames(V))) colnames(X) <- colnames(X)
  else colnames(X) <- colnames(V) <- paste0('V',1:ncol(X))

  results <- recurse_tolowerdimension(X, V = V,
                                      evaluator = evaluatesubsimplex,
                                      projector = projectsubsimplex,
                                      testweights = testweights,
                                      merges = merges,
                                      maxsteps = maxsteps
  )

  results$X = X
  results$Xhat = list()
  for(i in names(results$pts)){
    results$Xhat[[i]] = results$pts[[i]] %*% results$vertices[[i]]
  }

  return(results)
}
