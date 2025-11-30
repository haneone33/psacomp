#' @title Map data to the simplex
#' @description Map data to the simplex while preserving proportions.
#' @param X a data matrix
#' @return a matrix of the same dimension as `X`
#' @export
to_simplex <- function(X){
  if(min(X) < 0){
    warning('`X` has negative values. Replaced by zero and normalized')
    X[X < 0] = 0
    X = sweep(X, 1, rowSums(X), "/")
  }
  if(any(abs(rowSums(X)-1)>1e-12)){
    warning('rows of `X` does not sum to one. Rows are normalized')
    X = sweep(X, 1, rowSums(X), "/")
  }
  return(X)
}
