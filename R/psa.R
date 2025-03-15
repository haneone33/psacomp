#' @title Principal Nested Simplices Analysis
#' @description Estimate PSA-S or PSA-O of given data matrix.
#'
#' @param type `s` for PSA-S or `o` for PSA-O.
#' @param X A data matrix.
#' @param testweights A vector of weights to try.
#' @return A list
#' + vertices a list of matrix representing vertices of the lower dimensional subsimplex.
#'The $r$th element of the list corresponds to the rank $r-1$ subsimplex.
#' + pts a list of lower dimensional representations with respect to the reduced basis `vertices`
#' + pts.approx a list of lower dimensional representations with respect to the original basis
#' + scores a matrix of scores.
#' + rss a vector of residual sums of squares.
#' + modes a list of modes of variation. The $r$th element of the list is the difference of
#' rank $r$ approximations to rank $r-1$ approximations.
#' + loadings a matrix of loading vectors.
#' + const.info a data frame of merged vertices and merging weight at each merge.
#' + dendrogram.input .
#'
#' @export

psa <- function(type, X, testweights = seq(0, 1, length.out = 100)){
  n = nrow(X)
  d = ncol(X)-1
  if(is.null(colnames(X))) colnames(X) = paste0('e',1:ncol(X))
  var.names = colnames(X)

  ## output format
  res = list(vertices = list(),       # vertices of low rank subsimplices
             pts = list(),            # low rank approximations w.r.t. vertices
             pts.approx = list(),     # low rank approximations w.r.t. original vertices
             scores = list(),         # scalars
             rss = list(),            # residual sum of squares
             modes = list(),          # residuals (vectors; in Delta_d scale)
             loadings = NULL,         # loading vectors (matrix)
             const.info = list())     # merge indices and weight

  if(type == 's'){
    out = psas(X, testweights = testweights)
    out$pts.approx = list()
    for(i in names(out$pts)){
      out$pts.approx[[i]] = out$pts[[i]] %*% out$vertices[[i]]
    }
  }else if(type == 'o'){
    out = psao(X, testweights = testweights)
    out$pts.approx = list()
    for(i in names(out$pts)){
      out$pts.approx[[i]] = divL1(divL2(out$pts[[i]]) %*% out$vertices[[i]])
    }
  }else{
    stop('type should be either `s` or `o`')
  }

  res$vertices = rev(out$vertices)
  res$pts = rev(out$pts)
  res$pts.approx = rev(out$pts.approx)

  res$scores = rev(out$scores)
  res$scores = do.call('cbind', res$scores)
  colnames(res$scores) = paste0('Comp.',1:ncol(res$scores))
  res$rss = apply(res$scores, 2, function(vec) sum(vec^2))
  res$rss[d+1]=0

  res$modes = rev(out$modes)
  res$const.info = out$merges

  res$loadings = get_loading(out)
  colnames(res$loadings) = paste0('Comp.',1:ncol(res$loadings))

  return(res)
}
