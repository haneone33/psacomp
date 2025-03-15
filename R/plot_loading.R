#' @title Loading Plot for PSA
#' @description Creating a bar plot summarizing a loading vectors of PSA
#' @import ggplot2
#'
#' @param psa.res output of `psa()`
#' @param k number of loading vectors to display
#' @inheritParams plot_bar
#'
#' @return a ggplot object
#'
#' @export
plot_loading <- function(V, k = 4, max.k = 12){

  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop(
      "Package \"cowplot\" must be installed to use this function.",
      call. = FALSE
    )
  }

  ls = list()
  for(i in 1:k){
    ls[[i]] = plot_bar(V[,i], max.k) + ggtitle(paste0('Comp.',i))
  }
  g = cowplot::plot_grid(plotlist = ls, nrow = 1)
  return(g)
}

