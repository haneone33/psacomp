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

  ## output format
  res = list(Vhat = list(),           # vertices of low rank subsimplices
             Xhat = list(),           # low rank approximations w.r.t. vertices
             Xhat_reduced = list(),   # low rank approximations w.r.t. original vertices
             X = X,                   # given data
             residuals = list(),      # residuals (vectors; in Delta_d scale)
             scores = list(),         # scalars
             RSS = list(),            # residual sum of squares
             backwards_mean = NULL,   # backwards mean
             loadings = NULL,            # loading vectors (matrix)
             construction_info = list())     # merge indices and weight

  if(type == 's'){
    out = psas(X, testweights = testweights)
    out$Xhat = list()
    for(i in names(out$pts)){
      out$Xhat[[i]] = out$pts[[i]] %*% out$vertices[[i]]
    }
  }else if(type == 'o'){
    out = psao(X, testweights = testweights)
    out$Xhat = list()
    for(i in names(out$pts)){
      out$Xhat[[i]] = divL1(divL2(out$pts[[i]]) %*% out$vertices[[i]])
    }
  }else{
    stop('type should be either `s` or `o`')
  }

  res$Vhat = rev(out$vertices)
  res$Xhat = rev(out$Xhat)
  res$Xhat_reduced = rev(out$pts)
  names(res$Vhat) = paste0('r=',0:d)
  names(res$Xhat) = paste0('r=',0:d)
  names(res$Xhat_reduced) = paste0('r=',0:d)

  res$scores = rev(out$scores)
  res$scores = do.call('cbind', res$scores[1:d])
  colnames(res$scores) = paste0('r=',1:d)
  res$RSS = apply(res$scores, 2, function(vec) sum(vec^2))

  res$backwards_mean = res$Vhat$`r=0`
  res$residuals = lapply(1:d, function(r) res$Xhat[[paste0('r=',r)]] - res$Xhat[[paste0('r=',r-1)]])
  names(res$residuals) = paste0('r=',1:d)
  res$construction_info = out$merges

  res$loadings = get_loading(out)

  return(res)
}
