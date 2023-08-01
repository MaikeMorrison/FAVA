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
#' @param relab_matrix A matrix with \code{I=nrow(relab_matrix)} rows, each containing \code{K=ncol(relab_matrix)} non-negative entries that sum to 1.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix and the number of entries
#' that sum to 1 (\code{K}) must be specified.
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
#' @return A ggplot object containing a bar plot visualization of the relative abundance matrix.
#' @examples
#' plot_relaband(
#'   # Make an example matrix of compositional data
#'   # Each row is an individual. Rows sum to 1.
#'   relab_matrix = matrix(c(
#'     .4, .2, .4,
#'     .5, .3, .2,
#'     .5, .4, .1,
#'     .6, .1, .3,
#'     .6, .3, .1
#'   ),
#'   nrow = 5,
#'   byrow = TRUE
#'   ),
#'  K = 3, # How many categories per vector?
#'  arrange = TRUE
#' ) +
#'   # Below are example, optional modifications to the default plot
#'   ggplot2::ggtitle("Population A") +
#'   ggplot2::scale_fill_brewer("Blues") +
#'   ggplot2::scale_color_brewer("Blues") +
#'   ggplot2::xlab("Individuals")
#'   # Note that both scale_fill and scale_color are needed to change the color of the bars.
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
plot_relaband <- function(relab_matrix, group = NULL, time = NULL, w = NULL, K = NULL, arrange = FALSE) {

  if(!is.null(group)){

    relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

    if(is.null(relab_checker_out$group)){
      stop("The group provided is not a column name in relab_matrix. Please provide a valid group.")
    }

    # If group was provided but K was not, we assume group is the only
    # descriptive column in the data set
    if(is.null(K)) K = ncol(relab_matrix) - 1

    # Repeat rows to account for time or weight (w) if provided,
    # otherwise return relab_matrix unaltered
    relab_edited = relab_sample_weighter(relab = relab_matrix, K = K, time = time, w = w, group = group)

    # Re-arrange rows or columns as specified by arrange
    # otherwise return relab_matrix unaltered
    relab_edited = arrange_categories(relab_matrix = relab_edited,
                       arrange = arrange,
                       K = K, group = group, time = time)

    # MAIKE WAS HERE

    # Generate the data to plot
    df <- data.frame(cbind(data.frame(Individuals = 1:nrow(relab_matrix)), relab_matrix)) %>%
      tidyr::pivot_longer(cols = 3:(ncol(relab_matrix) + 1))
    df$name <- factor(df$name, levels = unique(df$name) %>% rev())
    df$Individuals = factor(df$Individuals)


    ggplot2::ggplot(data = df, ggplot2::aes(fill = .data$name,
                                            color = .data$name,
                                            y = .data$value,
                                            x = .data$Individuals)) +
      ggplot2::geom_bar(position = "stack", stat = "identity",
                        width = 1) + ggplot2::theme_void() + ggplot2::ylab("") +
      ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(size = 12)) +
      ggplot2::facet_wrap(~ group, scales = "free_x")


  }else{

    # Clean the matrices for plotting:
    relab_matrix <- relab_matrix_checker(relab_matrix = relab_matrix, K = K)

    # Generate the data to plot
    df <- data.frame(cbind(data.frame(Individuals = 1:nrow(relab_matrix)), relab_matrix)) %>%
      tidyr::pivot_longer(cols = 2:(ncol(relab_matrix) + 1))
    df$name <- factor(df$name, levels = unique(df$name) %>% rev())

    # Re-order vertical bars if arrange == TRUE, "both", or "horizontal",
    # Re-order categories if arrange == TRUE, "both" or "vertical"
    if (!missing(arrange) & arrange !=FALSE){
      if(arrange %in% c(TRUE, "both")){
        clustermeans <- colMeans(relab_matrix) %>% sort() %>% rev
        ordernames <- names(clustermeans)
        relab_matrix <- data.frame(relab_matrix) %>%
          dplyr::arrange(dplyr::across({{ ordernames }})) %>%
          dplyr::select(names(which(clustermeans != 0)))

        if(sum(clustermeans != 0) != K){
          warning(paste0("Only plotting the ", sum(clustermeans != 0),
                         " categories with non-zero abundances. If you are manually changing the fill or color of the plot, please provide ",
                         sum(clustermeans != 0), " colors, instead of ", K, "."))
        }

      }else if(arrange == "vertical"){
        clustermeans <- colMeans(relab_matrix) %>% sort() %>% rev
        ordernames <- names(clustermeans)
        relab_matrix <- data.frame(relab_matrix) %>%
          # dplyr::arrange(dplyr::across({{ ordernames }})) %>%
          dplyr::select(names(which(clustermeans != 0)))

        if(sum(clustermeans != 0) != K){
          warning(paste0("Only plotting the ", sum(clustermeans != 0),
                         " categories with non-zero abundances. If you are manually changing the fill or color of the plot, please provide ",
                         sum(clustermeans != 0), " colors, instead of ", K, "."))
        }

      }else if(arrange == "horizontal"){
        clustermeans <- colMeans(relab_matrix) %>% sort() %>% rev
        ordernames <- names(clustermeans)
        relab_matrix <- data.frame(relab_matrix) %>%
          dplyr::arrange(dplyr::across({{ ordernames }}))# %>%
        # dplyr::select(names(which(clustermeans != 0)))
      }else{
        stop("The options for arrange are TRUE, vertical, horizontal, or both. TRUE and both are interchangeable.")
      }

    }

    # Generate the data to plot
    df <- data.frame(cbind(data.frame(Individuals = 1:nrow(relab_matrix)), relab_matrix)) %>%
      tidyr::pivot_longer(cols = 2:(ncol(relab_matrix) + 1))
    df$name <- factor(df$name, levels = unique(df$name) %>% rev())

    # Generate the plot
    ggplot2::ggplot(
      data = df,
      ggplot2::aes(fill = .data$name, color = .data$name, y = .data$value, x = .data$Individuals)
    ) +
      ggplot2::geom_bar(position = "stack", stat = "identity", width = 1) +
      ggplot2::theme_void() +
      ggplot2::ylab("") +
      ggplot2::theme(legend.position = "none")
  }

}




# helper function: arrange

arrange_categories <- function(relab_matrix, arrange, K, group, time){

  relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

  relab_matrix = relab_checker_out$relab_matrix

  if (arrange !=FALSE){
    clustermeans <- colMeans(relab_matrix) %>% sort() %>% rev
    ordernames <- c(names(clustermeans))

    if(sum(clustermeans != 0) != K){
      warning(paste0("Only plotting the ", sum(clustermeans != 0),
                     " categories with non-zero abundances. If you are manually changing the fill or color of the plot, please provide ",
                     sum(clustermeans != 0), " colors, instead of ", K, "."))
    }

    if(arrange %in% c(TRUE, "both")){

      relab_matrix <- data.frame(relab_matrix) %>%
        dplyr::arrange(dplyr::across({{ ordernames }})) %>%
        dplyr::select(c(names(which(clustermeans != 0))))

    }else if(arrange == "vertical"){
      relab_matrix <- data.frame(relab_matrix) %>%
        # dplyr::arrange(dplyr::across({{ ordernames }})) %>%
        dplyr::select(c(names(which(clustermeans != 0))))


    }else if(arrange == "horizontal"){
      relab_matrix <- data.frame(relab_matrix) %>%
        dplyr::arrange(dplyr::across({{ ordernames }}))# %>%
      # dplyr::select(c("group", names(which(clustermeans != 0))))
    }else{
      stop("The options for arrange are FALSE, TRUE, vertical, horizontal, or both. TRUE and both are interchangeable.")
    }
  }

  if(!is.null(group)){
    relab_matrix = relab_matrix %>% data.frame %>% dplyr::mutate(group = relab_checker_out$group, .before = 1)
  }

  if(!is.null(time)){
    relab_matrix = relab_matrix %>% data.frame %>% dplyr::mutate(time = relab_checker_out$time, .before = 1)
  }

  return(relab_matrix %>% data.frame)
}
