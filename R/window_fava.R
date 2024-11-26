# window_list -----------------------------------------------------------------
#' Generate sliding windows of specified length given the maximum number of samples
#'
#' This function generates a list of of sliding windows conditional on two parameters: the length of each window (number of samples) and the total number of samples present in the data.
#'
#' @param window_size An integer number specifying the number of samples per window.
#' @param length An integer number specifying the total number of samples.
#' @param window_step Optional; an integer number specifying the distance between the first entry of adjacent windows. Default is \code{window_step=1}.
#' @returns A list of samples of sample indices. Each list entry represents one window.
#' @examples
#' window_list(window_size = 6, length = 40)
#' window_list(window_size = 6, length = 40, window_step = 2)
#' @export
window_list <- function(window_size, length, window_step = 1){
  if(window_size <= 0){stop("Window_size must be greater than zero")}
  if(length <= 0){stop("length (total number of samples) must be greater than zero")}
  if(window_size > length){stop("window_size must be less than the total number of samples")}
  if(length < window_size + window_step){stop("There are not enough samples to yield multiple windows with the current window_size and window_step parameters.")}

  winlist = list()
  i=1
  for(win in seq(from = 1, to = length-window_size+1, by = window_step)){
    winlist[[i]] <- win:(win + window_size - 1)
    i = i+1
  }
  winlist
}


# window_fava -----------------------------------------------------------------
#' Compute FAVA in sliding windows.
#'
#' This function computes FAVA in sliding window slices of a dataset.
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param window_size An integer number specifying the number of samples per window.
#' @param window_step Optional; an integer specifying the distance between the first entry of adjacent windows. Default is \code{window_step=1}.
#' @param index Optional; a string specifying the name of the column in \code{relab_matrix} containing an index for each sample. For example, if \code{relab_matrix} contains time series data, \code{index} would be the column containing the time of each sample. If \code{index} is not specified but \code{time} is, \code{time} is by default used as the index.
#' @param group Optional; a string specifying the name of the column that describes which group each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param normalized Optional; should normalized FAVA be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#' @param alpha Optional; number between 0 and 1 specifying the opacity of the horizontal
#' lines plotted. Default is \code{alpha = 0.5}.
#' @returns A list of values of FAVA for each window.
#' @importFrom dplyr %>%
#' @examples
#' A = matrix(c(.3,.7,0,.1,0,.9,.2,.5,.3,.1,.8,.1,.3,.4,.3,.6,.4,0,0,.5,.5),
#'            ncol = 3, byrow = TRUE)
#' window_out = window_fava(relab_matrix = A, window_size = 4, normalized = TRUE)
#' @export
window_fava <- function(relab_matrix, window_size, window_step = 1,
                        group = NULL,
                        index = NULL,
                        time = NULL, w = NULL,
                        S = NULL,
                        K = NULL, normalized = FALSE,
                        alpha = 0.5){

  # If time is provided but index is not, make index same as time
  if((!is.null(time)) & (is.null(index))) index = time

  # Process relab to allow for multiple grouping variables:
  if(length(group) > 1){
    relab_grouping_vars = dplyr::select(relab_matrix, dplyr::all_of(group))

    relab_grouping_vars$grouping_var_multiple = apply(relab_grouping_vars, 1, paste, collapse = "_")

    relab_matrix = dplyr::mutate(relab_matrix,
                                 grouping_var_multiple = relab_grouping_vars$grouping_var_multiple,
                                 .before = 1)

    group = "grouping_var_multiple"

    if(is.null(K)){K = ncol(relab_matrix) - (length(group)) - (!is.null(time)) - 1}
  }

  # Should the plot say "Weighted" or "Normalized" on the y-axis?
  y_axis = ifelse(!(is.null(time) & is.null(w) & is.null(S)),
                  "Weighted FAVA",
                  ifelse(normalized,
                         "Normalized FAVA",
                         "FAVA"))

  # Arrange data by time if not already
  if(!is.null(time)){
    if(!is.null(group)){
      relab_matrix = dplyr::arrange(.data = relab_matrix, .data[[group]], .data[[time]])
    } else{
      relab_matrix = dplyr::arrange(.data = relab_matrix, .data[[time]])
    }}


  if(is.null(group)){
    window_indices = window_list(window_size = window_size, length = nrow(relab_matrix),
                                 window_step = window_step)

    window_data = window_fava_sub(relab_matrix = relab_matrix,
                                  window_indices = window_indices, window_size = window_size,
                                  time = time, w = w, K = K, S = S, normalized = normalized)

    # Replace generic index with specified index or time if applicable
    if((!is.null(index))){

      window_data = tidyr::pivot_longer(window_data, cols = paste0("w", 1:window_size),
                                        names_to = "window_number", values_to = "generic_index") %>%
        dplyr::left_join(data.frame(original_index = relab_matrix[,index] %>% unlist,
                                    generic_index = 1:nrow(relab_matrix))) %>%
        dplyr::select(-generic_index) %>%
        tidyr::pivot_wider(names_from = window_number, values_from = original_index)
    }

  } else{

    group_fava_list = list()
    i = 1

    relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

    for(element in unique(unlist(relab_matrix[,group]))){

      relab_matrix_group = relab_matrix[relab_checker_out$group == element,]

      window_indices = window_list(window_size = window_size, length = nrow(relab_matrix_group), window_step = window_step)

      window_data_group = dplyr::mutate(window_fava_sub(relab_matrix = relab_matrix_group, group = group,
                                                        window_indices = window_indices, window_size = window_size,
                                                        time = time, w = w, K = K, S = S, normalized = normalized),
                                        group = element, .before = "FAVA")

      # Replace generic index with specified index or time if applicable
      if((!is.null(index))){

        window_data_group = tidyr::pivot_longer(window_data_group, cols = paste0("w", 1:window_size),
                                          names_to = "window_number", values_to = "generic_index") %>%
          dplyr::left_join(data.frame(original_index = relab_matrix_group[,index] %>% unlist,
                                      generic_index = 1:nrow(relab_matrix_group))) %>%
          dplyr::select(-generic_index) %>%
          tidyr::pivot_wider(names_from = window_number, values_from = original_index)
      }

      group_fava_list[[i]] = window_data_group

      i = i + 1
    }

    window_data = do.call(rbind, group_fava_list)
  }

return(list(window_data = window_data,
            window_plot = window_plot(window_fava = window_data, alpha = alpha) +
              ggplot2::ylab(y_axis)))
}

# To appease R command check
utils::globalVariables(c("generic_index", "original_index", "window_number"))

# helper function: compute FAVA for given windows on a full dataset
window_fava_sub = function(relab_matrix, window_indices, window_size,
                           time = NULL, w = NULL,
                           K = NULL, S = NULL, group = NULL,
                           normalized = FALSE){
  fava_list = c()
  for(window in window_indices){
    if(!is.null(w)){
      w_subset = w[window]/sum(w[window])
    }else{
      w_subset = NULL
    }

    fava_out = fava(relab_matrix = relab_matrix[window,],
                    time = time, K = K, normalized = normalized,
                    w = w_subset, S = S, group = group)

    if(is.data.frame(fava_out)){
      fava_out = unlist(fava_out[,2])
    }

    fava_list = c(fava_list,
                  fava_out)
  }

  return(cbind(data.frame(FAVA = fava_list,
                          window_index = 1:length(fava_list)),
               do.call(rbind, window_indices) %>%
                 `colnames<-`(paste0("w", 1:window_size))))
}


# window_plot -----------------------------------------------------------------
#' Generate a plot of FAVA in sliding windows.
#'
#' This function generates a plot of normalized or unnormalized, weighted or
#' unweighted FAVA computed in sliding windows across samples for one or many
#' groups of samples.
#'
#' @param window_fava The output of \code{window_fava}.
#' @param alpha Optional; number between 0 and 1 specifying the opacity of the horizontal
#' lines plotted. Default is \code{alpha = 0.5}.
#' @returns A ggplot2 object.
#' @examples
#' A = matrix(c(.3,.7,0,.1,0,.9,.2,.5,.3,.1,.8,.1,.3,.4,.3,.6,.4,0,0,.5,.5),
#'            ncol = 3, byrow = TRUE)
#' window_out = window_fava(relab_matrix = A, window_size = 4, normalized = TRUE)
#' window_out$window_data
#' window_out$window_plot
#' @export
window_plot <- function(window_fava, alpha = 0.5){
  # To appease R cmd check
  group <- FAVA <- window_index <- index <- NULL

  # If data has multiple groups ----------------------------------------
  if("group" %in% colnames(window_fava)){
    window_long = window_fava %>%
      tidyr::pivot_longer(-c(group, FAVA, window_index),
                          values_to = "index", names_to = "window_count")

    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = index,
                                      y = FAVA,
                                      color = group,
                                      group = paste0(window_index, group)),
                         window_long, linewidth = 1, alpha = alpha) +
      ggplot2::theme_bw() +
      ggplot2::ylab("FAVA") +
      ggplot2::xlab("Index") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  }else{
    # If data has only one group ----------------------------------------
    window_long = window_fava %>%
      tidyr::pivot_longer(-c(FAVA, window_index),
                          values_to = "index", names_to = "window_count")

    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = index,
                                      y = FAVA,
                                      group = paste0(window_index)),
                         window_long, linewidth = 1, alpha = alpha) +
      ggplot2::theme_bw() +
      ggplot2::ylab("FAVA") +
      ggplot2::xlab("Index") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  }

}
