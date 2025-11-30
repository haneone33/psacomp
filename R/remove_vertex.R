#' @title Find lower dimensional simplex from given index
#'
#' @param X (d+1) times n data matrix in full coordinate system.
#' @param V (d+1) times (r+1) basis matrix for reduced coordinate system.
#' @param idx the vertex to remove. If character, it is the column name of the vertex.
#'        If integer, it is the column index of the vertex.

#' @return
#' A list of the following components.
#' \item{a}{normal vector for the lower dimensional simplex}
#' \item{Xhat}{(d+1) times n lower dimesional representation in full coordinate system}
#' \item{Vhat}{estimated basis matrix}
#' \item{Xhat_reduced}{lower dimensional representation in reduced coordinate system}
#' \item{residuals}{vector of residuals}
#' \item{scores}{matrix of scores}
#' \item{RSS}{residual sum of squares}
#'
#' @noRd
remove_vertex <- function(X, V, idx){
  r = ncol(V) - 1
  if(!is.numeric(idx)){
    idx = which(colnames(V) == idx)
  }
  removed_vertex = V[,idx]

  ## Step 1: find normal vector a
  X = push_vertex_inwards(X, V, idx)
  X_reduced = mat_inverse(V) %*% X

  # pivoting: move the unique vertex to the last column
  new_idx = 1:(r+1)
  new_idx[c(idx, r+1)] = c(r+1, idx)
  V = V[,new_idx]
  X_reduced = X_reduced[new_idx,]
  X = V %*% X_reduced

  # optimizing
  opt_res = stats::optim(par = rep(1, r), fn = get_RSS, gr = RSS_gradient,
                         X = X, V = V,
                         method = 'L-BFGS-B', lower = rep(0, r),
                         control = list(factr = 1))

  ## Step 2: format output
  a = opt_res$par

  Xhat = get_Xhat(a, X, V)
  Vhat = get_Vhat(a, V)
  Xhat_reduced = mat_inverse(Vhat) %*% Xhat
  residuals = get_Res(a, X, V)
  scores = get_score(a, X, V)
  RSS = sum(scores^2)

  a = c(a, -1)
  a = a[new_idx]

  return(list(Xhat = Xhat, Vhat = Vhat, Xhat_reduced = Xhat_reduced,
              idx = idx, a = a, removed_vertex = removed_vertex,
              residuals = residuals, scores = scores, RSS = RSS))
}



