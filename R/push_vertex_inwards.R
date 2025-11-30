push_vertex_inwards <- function(X, V, idx, margin = 1e-10){
  X_reduced = mat_inverse(V) %*% X
  d = nrow(X_reduced)-1
  xnew = rep(margin/d, d+1)
  xnew[idx] = 1 - margin
  X_reduced[,X_reduced[idx,] > 1-margin] = xnew
  return(V %*% X_reduced)
}
