#' @title Single Bar Plot for PSA Loadings or Vertices
#' @description Creating a bar plot summarizing a loading vector or a vertex vector
#' @import ggplot2
#'
#' @param v a loading vector
#' @param max.k maximum number of elements to display
#'
#' @return a ggplot object
plot_bar <- function(v, max.k = 12){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  v = v[abs(v)>1e-8]
  if(length(v) > max.k){
    threshold = sort(abs(v), decreasing = T)[max.k]
    v = v[abs(v) >= threshold]
  }

  v = sort(v, decreasing = F)
  df = data.frame(variable = factor(names(v), levels = names(v)),
                  value = v,
                  vertex = factor(v>0, levels = c(T,F)))
  rownames(df) = NULL

  g = ggplot(data = df, aes(x = .data$variable, y = .data$value, fill = .data$vertex)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme_bw() +
    labs(x = '', y = '') +
    theme(axis.ticks.y = element_blank()) +
    theme(legend.position = 'none') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.margin = unit(c(2.5,2.5,2.5,2.5), "points"))

  return(g)
}
