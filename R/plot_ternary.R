#' @title Vertices Plot for PSA
#' @description Creating bar plots representing vertices of PSA
#'
#' @param V a matrix of vertices (each row representing a vertex)
#' @inheritParams plot_bar
#'
#' @return a ggplot object
plot_vertex <- function(V, max.k = 12){

  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop(
      "Package \"cowplot\" must be installed to use this function.",
      call. = FALSE
    )
  }

  k = nrow(V)
  ls = list()
  for(i in 1:k){
    ls[[i]] = plot_bar(V[i,], max.k) +
      ggtitle(paste0('V',i)) +
      scale_fill_manual(values = 'darkgray') +
      scale_y_continuous(breaks = pretty(2))
  }
  cowplot::plot_grid(plotlist = ls, nrow = 1)
}

#' @title Ternary Plot for PSA
#' @description Creating a ternary plot with its vertices for 2-dimensional
#' representation of obtained from PSA
#' @import ggplot2
#'
#' @param psa.res output of `psa()`
#' @param groups a vector for colors of data points if applicable
#' @param group.name title for legend if `groups` exists
#' @inheritParams plot_vertex
#'
#' @return a ggplot object
#'
#' @export
plot_ternary <- function(psa.res, groups = NULL, group.name = 'Groups', max.k = 12){

  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop(
      "Package \"ggtern\" must be installed to use this function.",
      call. = FALSE
    )
  }

  df = as.data.frame(psa.res$pts$`r=3`)
  colnames(df) = c('V1','V2','V3')

  if(is.null(groups)){
    g1 = ggtern::ggtern(df, aes(x = .data$V1, y = .data$V2, z = .data$V3))
  }else{
    df$Groups = Groups
    g1 = ggtern::ggtern(df, aes(x = .data$V1, y = .data$V2, z = .data$V3,
                                col = .data$Groups))
  }

  g1 = g1 + geom_point(shape = 1) +
    theme_bw() +
    theme(legend.margin = margin(0,-60,0,0),
          legend.position = 'left',
          legend.justification = c(0, 0.9),
          plot.margin = margin(0,-30,0,30))

  if(!is.null(groups)) g1 = g1 + guides(fill = guide_legend(title = group.name))

  g2 = plot_vertex(psa.res$vertices$`r=3`, max.k) +
    theme(plot.margin = margin(30,0,30,0))

  g = cowplot::plot_grid(ggtern::ggplot_gtable(ggtern::ggplot_build(g1)),
                         g2,
                         nrow = 1, rel_widths = c(1.5,1))

  return(g)

}
