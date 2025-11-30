#' @title Principal Subsimplex Analysis
#' @description Estimate PSA-S or PSA-O of given data matrix.
#'
#' @param X a data matrix.
#' @param type Type of PSA. 's' for PSA-S, 'o' for PSA-O, or 'g' for PSA-G.
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

psa <- function(X, type = c('s','o','g'), testweights = seq(0, 1, length.out = 100)){

  ## setup
  type = match.arg(type)

  X = to_simplex(X)
  if(is.null(colnames(X))) colnames(X) = paste0('e',1:ncol(X))
  n = nrow(X)
  d = ncol(X)-1

  ## core computation
  if(type == 's'){
    out = psas_format_output(psas(X, testweights = testweights))
  }else if(type == 'o'){
    out = psao_format_output(psao(X, testweights = testweights))
  }else if(type == 'g'){
    out = psag(X)
  }

  return(out)
}
