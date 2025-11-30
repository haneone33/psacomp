#' @title Projection map for PSA-G from a unit simplex
#' @param a normal vector defining the lower dimensional simplex. A vector of
#' length \eqn{r}.
#' @param X (d+1) times n data matrix. The last row must correspond to the vertex being lost.
#' @return data matrix of the same dimension.
#' @noRd
psag_projection_unit <- function(a, X){

  d = nrow(X) - 1

  atx = t(c(a, -1)) %*% X
  Xhat = X
  Xhat[d+1,] = Xhat[d+1,] + atx
  Xhat = sweep(Xhat, 2, atx+1, "/")
  return(Xhat)
}

#' @title Get PSA-G approximation with respect to a normal vector
#' @param a normal vector defining the lower dimensional simplex of length \eqn{r}.
#' @param X (d+1) times n data matrix.
#' @param V (d+1) times (r+1) basis matrix. The last column must correspond to the vertex to be lost.
#' @return (d+1) times n data matrix.
#' @noRd
get_Xhat <- function(a, X, V){
  # compute coordinates of X in the reduced coordinate system
  Y = mat_inverse(V) %*% X
  # compute coordinates of Yhat in the full coordinate system
  Xhat = V %*% psag_projection_unit(a, Y)
  return(Xhat)
}

#' @title Get the vertex set of the subsimplex with respect to a normal vector
#' @param a normal vector defining the lower dimensional simplex of length \eqn{r}.
#' @param V (d+1) times (r+1) basis matrix.
#' @return (d+1) times r basis matrix
#' @noRd
get_Vhat <- function(a, V){
  Vnew = diag(1, ncol(V))[,1:(ncol(V)-1)]
  #colnames(Vnew) = colnames(V)[1:(ncol(V)-1)]
  colnames(Vnew) = paste0('V',1:(ncol(V)-1))
  Vhat = V %*% psag_projection_unit(a, Vnew)
  return(Vhat)
}

get_Res_y <- function(a, X, V){
  Y = mat_inverse(V) %*% X
  Yhat = psag_projection_unit(a, Y)
  Res = Y - Yhat
  return(Res)
}

get_Res <- function(a, X, V){
  Xhat = get_Xhat(a, X, V)
  Res = X - Xhat
  return(Res)
}

get_score <- function(a, X, V){
  Res = get_Res(a, X, V)
  Res_y = get_Res_y(a, X, V)
  scores = sign(Res_y[nrow(Res_y),]) * sqrt(colSums(Res^2))
  return(scores)
}

## RSS is computed directly from a, X, V for faster optimization
## equivalent to RSS = sum(get_Res(a,X,V)^2)
get_RSS <- function(a, X, V){
  Y = mat_inverse(V) %*% X
  Res = X - V %*% psag_projection_unit(a, Y)
  RSS = sum(Res^2)
  return(RSS)
}


RSS_gradient <- function(a, X, V){
  Y = mat_inverse(V) %*% X
  atY = t(c(a, -1)) %*% Y
  coeff = sweep(X, 1, V[,ncol(V)], '-')
  coeff = colSums(coeff^2)
  gr.mat = sweep(Y, 2, coeff*atY/((atY+1)**3), '*')
  gr = rowSums(gr.mat)
  return(gr[-length(gr)])
}
