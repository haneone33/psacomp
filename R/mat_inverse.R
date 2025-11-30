mat_inverse <- function(V){
  x = svd(V)
  V.inv = x$v%*%diag(1/x$d)%*%t(x$u)
  rownames(V.inv) = colnames(V)
  colnames(V.inv) = rownames(V)
  return(V.inv)
}
