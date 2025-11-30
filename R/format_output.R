psag_format_output <- function(res){
  d = length(res) - 1

  out = list(X = t(res[[paste0('r=',d)]]$Xhat),
             Vhat = NA, Xhat = NA, Xhat_reduced = NA,
             scores = NA, RSS = NA, loadings = NA,
             backwards_mean = NA,
             construction_info = list(idx = NA, a = NA))

  out$Vhat = lapply(res, function(ls) ls$Vhat)
  out$Xhat = lapply(res, function(ls) t(ls$Xhat))
  out$Xhat_reduced = lapply(res, function(ls) t(ls$Xhat_reduced))

  out$construction_info = lapply(res[1:d], function(ls) list(idx = ls$idx, a = ls$a))
  names(out$construction_info) = paste0('r=',1:d)

  out$scores = sapply(res[1:d], function(ls) ls$scores)
  colnames(out$scores) = paste0('PC',1:d)
  out$RSS = sapply(res[1:d], function(ls) ls$RSS)
  names(out$RSS) = paste0('PC',1:d)

  out$backwards_mean = out$Vhat$`r=0`
  out$loadings = sapply(res[1:d], function(ls) ls$removed_vertex - out$backwards_mean)
  colnames(out$loadings) = paste0('PC',1:d)
  rownames(out$loadings) = colnames(out$X)

  return(out)
}

psas_format_output <- function(res){
  d = length(res$scores) - 1
  X = res$X

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
    rownames(out$Vhat[[r]]) = colnames(X)
  }

  # Xhat
  out$Xhat = rev(res$Xhat)
  names(out$Xhat) = paste0('r=',0:d)
  for(r in 1:(d+1)){
    out$Xhat[[r]] = out$Xhat[[r]]
    colnames(out$Xhat[[r]]) = colnames(X)
  }

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
  rownames(out$loadings) = colnames(X)

  backwards_mean = c(out$Vhat$`r=0`)
  names(backwards_mean) = colnames(X)
  out$backwards_mean = backwards_mean
  out$construction_info = res$merges

  return(out)
}

psao_format_output <- function(res){
  return(psas_format_output(res))
}
