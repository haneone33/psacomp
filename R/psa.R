#' @title Principal Subsimplex Analysis
#' @description Estimate PSA-S or PSA-O of given data matrix.
#'
#' @param type 's' for PSA-S or 'o' for PSA-O.
#' @param X a data matrix.
#' @param testweights a vector of weights for grid search for alpha.
#' @return A list of the following components of PSA.
#' \item{Vhat}{a list of matrix representing vertices of the lower dimensional subsimplex.}
#' \item{Xhat}{a list of lower dimensional representations with respect to the original basis.}
#' \item{Xhat_reduced}{a list of lower dimensional representations with respect to the reduced basis `Vhat`}
#' \item{scores}{a matrix of scores.}
#' \item{X}{the input matrix.}
#' \item{residuals}{a list of residuals.}
#' \item{scores}{a matrix of scores.}
#' \item{RSS}{a vector of residual sums of squares.}
#' \item{backwards_mean}{the backwards mean. Equal to `Vhat$'r=0'`.}
#' \item{loadings}{a matrix of loading vectors.}
#' \item{construction_info}{a data frame of merged vertices and merging weight at each merge.}
#'
#' @export

psa <- function(type, X, testweights = seq(0, 1, length.out = 100)){

  ## setup
  if(! type %in% c('s','o')){
    stop('available types are `s` and `o`')
  }
  if(min(X) < 0){
    warning('`X` has negative values. Replaced by zero and normalized')
    X[X < 0] = 0
    X = sweep(X, 1, rowSums(X), "/")
  }
  if(any(abs(rowSums(X)-1)>1e-12)){
    warning('rows of `X` does not sum to one. Rows are normalized')
    X = sweep(X, 1, rowSums(X), "/")
  }

  n = nrow(X)
  d = ncol(X)-1
  if(is.null(colnames(X))) colnames(X) = paste0('e',1:ncol(X))

  ## core computation
  if(type == 's'){
    res = psas(X, testweights = testweights)
    res$Xhat = list()
    for(i in names(res$pts)){
      res$Xhat[[i]] = res$pts[[i]] %*% res$vertices[[i]]
    }
  }else if(type == 'o'){
    res = psao(X, testweights = testweights)
    res$Xhat = list()
    for(i in names(res$pts)){
      res$Xhat[[i]] = divL1(divL2(res$pts[[i]]) %*% res$vertices[[i]])
    }
  }else{
    stop('type should be either `s` or `o`')
  }

  ## output format
  out = list(X = X,                   # given data
             Vhat = list(),           # vertices of low rank subsimplices
             Xhat = list(),           # low rank approximations w.r.t. vertices
             Xhat_reduced = list(),   # low rank approximations w.r.t. original vertices
             scores = list(),         # scalars
             RSS = list(),            # residual sum of squares
             loadings = NULL,         # loading vectors (matrix)
             backwards_mean = NULL,   # backwards mean
             construction_info = list())     # merge indices and weight

  # Vhat
  out$Vhat = rev(res$vertices)
  names(out$Vhat) = paste0('r=',0:d)
  for(r in 1:(d+1)){
    out$Vhat[[r]] = t(out$Vhat[[r]])
    colnames(out$Vhat[[r]]) = paste0('V',1:ncol(out$Vhat[[r]]))
  }
  # Xhat
  out$Xhat = rev(res$Xhat)
  names(out$Xhat) = paste0('r=',0:d)
  # Xhat_reduced
  out$Xhat_reduced = rev(res$pts)
  names(out$Xhat_reduced) = paste0('r=',0:d)
  for(r in 1:(d+1)){
    colnames(out$Xhat_reduced[[r]]) = paste0('V',1:ncol(out$Xhat_reduced[[r]]))
  }

  out$scores = rev(res$scores)
  out$scores = do.call('cbind', out$scores[1:d])
  colnames(out$scores) = paste0('PC',1:d)
  out$RSS = apply(out$scores, 2, function(x) sum(x^2))
  out$loadings = get_loading(res)

  out$backwards_mean = out$Vhat$`r=0`
  out$construction_info = res$merges

  return(out)
}
