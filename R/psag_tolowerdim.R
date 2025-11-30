#' @title Find lower dimensional simplex
#'
#' @param X data matrix in full coordinate system.
#' @param V basis matrix for reduced coordinate system.
#' @param type a character indicating which algorithm to use. '1' for rescaling
#' after each step. '2' for no rescaling.
#' @return an outcome from `psag_format_optim`
#' @noRd
psag_tolowerdim <- function(X, V){
  r = ncol(V) - 1
  ## if r == 1 return zero dimensional representation and terminate, continue otherwise
  if(r == 1) return(psag_tozerodim(X, V))

  res = lapply(seq(r+1), function(i) remove_vertex(X, V, i))
  RSS = sapply(res, function(l) l$RSS)
  idx = which.min(RSS)

  return(remove_vertex(X, V, idx))
}
