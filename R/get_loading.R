get_loading <- function(merge.info){
  max_r = nrow(merge.info$merges)
  loadings.ls = vector(mode = 'list', length = max_r)

  for(r in 1:(max_r-1)){
    v1 = merge.info$vertices[[paste0('r=',r+1)]][merge.info$merges[paste0('r=',r),'v1'],]
    v2 = merge.info$vertices[[paste0('r=',r+1)]][merge.info$merges[paste0('r=',r),'v2'],]
    loadings.ls[[paste0('r=',r)]] = v2 - v1
  }

  loadings <- do.call('cbind', loadings.ls)
  return(loadings)
}
