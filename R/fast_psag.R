#' @title fast PSA-G
#'
#' @param X data matrix. May have column names.
#' @param n_components dimension of compositions to which the data is reduced by PCA
#' @param percent_appx percent of variance to be explained by PCA
#'
#' @return
#' A list of the following components of PSA-G.
#' \item{X}{the input matrix.}
#' \item{Xhat}{a list of lower dimensional representations with respect to the original basis.}
#' \item{Vhat}{a list of matrix representing vertices of the lower dimensional subsimplex.}
#' \item{Xhat_reduced}{a list of lower dimensional representations with respect to the reduced basis `Vhat`}
#' \item{scores}{a matrix of scores.}
#' \item{residuals}{a list of residuals.}
#' \item{scores}{a matrix of scores.}
#' \item{RSS}{a vector of residual sums of squares.}
#' \item{backwards_mean}{the backwards mean. Equal to `Vhat$'r=0'`.}
#' \item{modes}{a matrix of loading vectors.}
#' \item{construction_info}{a list of `idx` and `a`. `idx` contains the index of
#'  the vertex that were chosen at each step. `a` contains the corresponding normal vector.}
#'
#' @export
fast_psag <- function(X, n_components = NULL, percent_appx = 0.99){

  ## initialize
  X = to_simplex(X)
  if(is.null(colnames(X))) colnames(X) = paste0('e',1:ncol(X))
  n = nrow(X)
  d = ncol(X)-1

  V = diag(1, d+1)
  rownames(V) = colnames(X)
  colnames(V) = paste0('V',1:(d+1))

  ## setup result
  res = as.list(rep(NA, d+1))
  names(res) = paste0('r=',0:d)

  Xhat = t(X)
  Vhat = V
  Xhat_reduced = Xhat
  rownames(Xhat_reduced) = colnames(V)

  res[[paste0('r=',d)]] = list(Xhat = Xhat, Vhat = Vhat, Xhat_reduced = Xhat_reduced,
                               idx = NULL, a = NULL, removed_vertex = rep(NA, d+1),
                               residuals = matrix(NA, d+1, n),
                               scores = rep(NA, n),
                               RSS = NA)

  ## Step 1: approximated PSA-G to reduce dimension to n_components
  message('Setup Initialize')
  rss_per_index = sapply(colnames(V), function(i) remove_vertex(Xhat, V, i), simplify = F)
  RSS = sapply(rss_per_index, function(l) l$RSS)
  RSS_idx = sort(RSS, decreasing = T, index.return = T)

  if(is.null(n_components)){
    n_components = which(cumsum(RSS)/sum(RSS) > percent_appx)[1]
  }
  message(paste0('Selected n_components : ', n_components))

  message('Finish initialization')

  if(n_components <= d-1){
    for(r in d:(n_components+1)){
      i = RSS_idx$ix[r+1]
      idx = which(Vhat[i,]>0)
      opt_res = remove_vertex(Xhat, Vhat, idx)
      res[[paste0('r=',r-1)]] = opt_res
      Xhat = opt_res$Xhat
      Vhat = opt_res$Vhat
      message(paste0("Finished approximating dimension ", r-1, "."))
    }
  }

  ## Step 2: exact PSA-G
  for(r in n_components:1){
    opt_res = psag_tolowerdim(Xhat, Vhat)
    res[[paste0('r=',r-1)]] = opt_res
    Xhat = opt_res$Xhat
    Vhat = opt_res$Vhat
    message(paste0("Finished creating dimension ", r-1, "."))
  }

  out = psag_format_output(res)
  out$n_components = n_components

  return(out)
}

