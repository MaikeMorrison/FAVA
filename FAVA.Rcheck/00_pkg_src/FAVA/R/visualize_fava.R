# plot_relaband ----------------------------------------
#' Visualize a relative abundance matrix as a stacked bar plot.
#'
#' This function enables graphical visualization of a matrix of compostional data.
#'In the output plot, each vertical bar represents a single vector;
#' the height of each color in the bar corresponds to the abundance of each category
#' in that vector. Because this function produces a
#' ggplot object, its output can be modified using
#' standard ggplot2 syntax.
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param group Optional; a string specifying the name of the column that describes which group each row (sample) belongs to. Use if \code{matrices} is a single matrix containing multiple groups of samples you wish to compare.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param arrange Optional; controls horizontal ordering of samples and vertical ordering of categories.
#'   If \code{arrange = TRUE} or \code{arrange = "both"}, samples are ordered by the categories of greatest
#'   abundance and categories are ordered in decreasing abundance. If \code{arrange = "vertical"}, sample
#'   order is unchanged but categories are ordered in decreasing abundance. If \code{arrange = "horizontal"},
#'   samples are ordered by the most abundant categories, but category order is unchanged. If \code{arrange} is missing
#'   or \code{arrange = FALSE}, neither order is changed.
#' @returns A ggplot object containing a bar plot visualization of the relative abundance matrix.
#' @examples
#'
#' # Make an example matrix of compositional data
#' # Each row is an individual. Rows sum to 1.
#' population_A = matrix(c(
#'     .5, .3, .2,
#'     .4, .2, .4,
#'     .5, .4, .1,
#'     .6, .1, .3,
#'     .2, 0, .8
#'   ),
#'   nrow = 5,
#'   byrow = TRUE
#'   )
#'
#'   plot_relabund(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = FALSE
#'               )
#'   plot_relabund(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = "horizontal"
#'               )
#'   plot_relabund(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = "vertical"
#'               )
#'    plot_relabund(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = TRUE  # could also be "both"
#'               )
#'
#'
#' # You can modify the plot as you would any ggplot2 object
#' plot_relabund(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = TRUE
#'               ) +
#'   # Below are example, optional modifications to the default plot
#'   ggplot2::ggtitle("Population A") +
#'   ggplot2::scale_fill_brewer("Blues") +
#'   ggplot2::scale_color_brewer("Blues") +
#'   ggplot2::xlab("Individuals")
#'   # Note that both scale_fill and scale_color are needed to change the color of the bars.
#'
#'
#'   # Plot a dataset which has 2 populations
#'
#'   population_B = matrix(c(
#'     .9, 0, .1,
#'     .6, .4, 0,
#'     .7, 0, .3,
#'     .3, .4, .3,
#'     .5, .3, .2
#'   ),
#'   nrow = 5,
#'   byrow = TRUE
#'   )
#'
#'
#'   populations_AB = cbind(data.frame(c("A", "A", "A", "A", "A",
#'                                      "B", "B", "B", "B", "B")),
#'                          rbind(population_A, population_B))
#'   colnames(populations_AB) = c("population", "category_1", "category_2", "category_3")
#'
#'
#'  plot_relabund(relab_matrix = populations_AB, group = "population")
#'  plot_relabund(relab_matrix = populations_AB, group = "population", arrange = "vertical")
#'  plot_relabund(relab_matrix = populations_AB, group = "population", arrange = "horizontal")
#'  plot_relabund(relab_matrix = populations_AB, group = "population", arrange = "both")
#'
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
plot_relabund <- function(relab_matrix, group = NULL, time = NULL, w = NULL, K = NULL, arrange = FALSE) {

  relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

  K = ncol(relab_checker_out$relab_matrix)

  # Repeat rows to account for time or weight (w) if provided,
  # otherwise return relab_matrix unaltered
  relab_edited = relab_sample_weighter(relab = relab_matrix, K = K, time = time, w = w, group = group)

  # Re-arrange rows or columns as specified by arrange
  # otherwise return relab_edited unaltered
  relab_edited = arrange_categories(relab_matrix = relab_edited,
                                    arrange = arrange,
                                    K = K, group = group, time = time)


  # Generate the data to plot
  relab_plot = dplyr::mutate(relab_edited, ID = 1:nrow(relab_edited), .before = 1)


  start =  2 + (!is.null(group)) + (!is.null(time))

  relab_plot_long = tidyr::pivot_longer(relab_plot, cols = start:ncol(relab_plot))


  relab_plot_long$ID = factor(relab_plot_long$ID , levels = unique(relab_plot$ID), ordered = TRUE)
  relab_plot_long$name = factor(relab_plot_long$name, levels = colnames(relab_plot)[start:ncol(relab_plot)] %>% rev, ordered = TRUE)

  if(!is.null(group)){

    if(is.null(relab_checker_out$group)){
      stop("The group provided is not a column name in relab_matrix. Please provide a valid group.")
    }

    ggplot2::ggplot(data = relab_plot_long, ggplot2::aes(fill = .data$name,
                                            color = .data$name,
                                            y = .data$value,
                                            x = .data$ID)) +
      ggplot2::geom_bar(position = "stack", stat = "identity",
                        width = 1) + ggplot2::theme_void() + ggplot2::ylab("") +
      ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(size = 12)) +
      ggplot2::facet_wrap(~ group, scales = "free_x")
  }else{

    ggplot2::ggplot(data = relab_plot_long, ggplot2::aes(fill = .data$name,
                                                         color = .data$name,
                                                         y = .data$value,
                                                         x = .data$ID)) +
      ggplot2::geom_bar(position = "stack", stat = "identity",
                        width = 1) + ggplot2::theme_void() + ggplot2::ylab("") +
      ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(size = 12))

  }
}




