# THIS FILE CONTAINS:
# windows

# window_list -----------------------------------------------------------------
#' Generate sliding windows of specified length given the maximum number of samples
#'
#' This function generates a list of of sliding windows conditional on two parameters: the length of each window (number of samples) and the total number of samples present in the data.
#'
#' @param window_size An integer number specifying the number of samples per window.
#' @param length An integer number specifying the total number of samples.
#' @param window_step Optional; an integer specifying the distance between the first entry of adjacent windows. Default is \code{window_step=1}.
#' @returns A list of samples of sample indices. Each list entry represents one window.
#' @examples
#' window_list(6,40)
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
#' @param Q A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{Q} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param window_size An integer number specifying the number of samples per window.
#' @param K Optional; number of categories which sum to 1 for each sample. Default is \code{K=ncol(Q)}.
#' @param window_step Optional; an integer specifying the distance between the first entry of adjacent windows. Default is \code{window_step=1}.
#' @param group Optional; a string specifying the name of the column that describes which group each row (sample) belongs to. Use if \code{Q} is a
#' single matrix containing multiple groups of samples you wish to compare, such as samples from different subjects.
#' @param normalized Default is \code{normalized = FALSE}. If \code{normalized = TRUE}, then fava is normalized by the theoretical upper bound conditional
#' on the most abundant category and the number of categories. This setting is recommended when there are fewer than 5 samples per window. Note that it is
#' incompatible with weightings at this time.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} i
#' n the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(Q), nrow(Q))}. If \code{group} is specified,
#' then \code{w} can contain the weights for multiple groups, concatenated in the same order as the rows of \code{Q}. In this case, \code{w} must sum to the
#' number of distinct elements of \code{group}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for
#' \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling
#' 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(Q))}.
#' @param index Optional; the name of the column in \code{Q} containing an index for each sample. For example, if \code{Q} contains time series data,
#' \code{index} would be the column containing the time of each sample.
#' @returns A list of values of fava for each window.
#' @importFrom dplyr %>%
#' @examples
#' A = matrix(c(.3,.7,0,.1,0,.9,.2,.5,.3,.1,.8,.1,.3,.4,.3,.6,.4,0,0,.5,.5),
#'            ncol = 3, byrow = TRUE)
#' window_fava(Q = A, window_size = 4, normalized = TRUE)
#' @export
window_fava <- function(Q, K = ncol(Q), window_size, group,
                       window_step = 1, normalized = FALSE,
                       w = rep(1/nrow(Q), nrow(Q)), S = diag(K),
                       index){
  Q_full = Q
  # UNGROUPED DATA ------------------------------------------------------------
  if(missing(group)){
    Q = Q[,(ncol(Q)-K+1):ncol(Q)]

    window_indices = window_list(window_size = window_size, length = nrow(Q), window_step = window_step)

    fava_list = c()
    for(window in window_indices){
      if(normalized){
        fava_list = c(fava_list, fava_norm(Q[window,]))
      }else{
        fava_list = c(fava_list,
                     fava(Q = Q[window,], w = w[window]/sum(w[window]), S = S))
      }
    }

    output = cbind(data.frame(FAVA = fava_list, window_index = 1:length(fava_list)),
          do.call(rbind, window_indices) %>% `colnames<-`(paste0("w", 1:window_size)))

    if(missing(index)){ return(output) }else{

        tidyr::pivot_longer(output, -c(FAVA, window_index),
                            names_to = "window_number", values_to = "generic_index") %>%
        left_join(data.frame(original_index = Q_full[,index] %>% unlist,
                 generic_index = 1:nrow(Q_full))) %>%
        select(-generic_index) %>%
        pivot_wider(names_from = window_number, values_from = original_index)
    }
  }else{

    # GROUPED DATA -------------------------------------------------------------

    group_fava_list = list()
    i = 1

    for(element in dplyr::select(Q, {{ group }}) %>% unique %>% unlist){
      Q_group = dplyr::filter(Q, get(group) == element)[,(ncol(Q)-K+1):ncol(Q)]
      window_indices = window_list(window_size = window_size, length = nrow(Q_group),
                                   window_step = window_step)

      fava_list = c()
      for(window in window_indices){
        if(normalized){
          fava_list = c(fava_list, fava_norm(Q_group[window,]))
        }else{
          fava_list = c(fava_list,
                       fava(Q = Q_group[window,],
                           w = w[window]/sum(w[window]), S = S))
        }
      }

      group_fava_list[[i]] = cbind(data.frame(group = element,
                                             FAVA = fava_list,
                                             window_index = 1:length(fava_list)),
                                  do.call(rbind, window_indices) %>%
                                    `colnames<-`(paste0("w", 1:window_size)))
      i = i + 1

    }
    output = do.call(rbind, group_fava_list)

    if(missing(index)){ return(output) }else{

      tidyr::pivot_longer(output,
                          -c(FAVA, window_index, group),
                          names_to = "window_number",
                          values_to = "generic_index") %>%
        left_join(data.frame(group = unlist(Q[,group]),
                             original_index = unlist(Q[,index])) %>%
                    group_by(group) %>% mutate(generic_index = 1:n())) %>%
        select(-generic_index) %>%
        pivot_wider(names_from = window_number, values_from = original_index)
    }
  }
}

# window_plot -----------------------------------------------------------------
#' Generate a plot of FAVA in sliding windows.
#'
#' This function generates a plot of normalized or unnormalized, weighted or
#' unweighted FAVA computed in sliding windows for one or many groups of samples.
#'
#' @param window_fava The output of \code{window_fava}.
#' @param alpha Optional; number between 0 and 1 specifying the opacity of the horizontal
#' lines plotted. Default is \code{alpha = 0.5}.
#' @returns A ggplot2 object.
#' @examples
#' A = matrix(c(.3,.7,0,.1,0,.9,.2,.5,.3,.1,.8,.1,.3,.4,.3,.6,.4,0,0,.5,.5),
#'            ncol = 3, byrow = TRUE)
#' window_out = window_fava(Q = A, window_size = 4, normalized = TRUE)
#' window_plot(window_fava = window_out, alpha = 0.8)
#' @export
window_plot <- function(window_fava, alpha = 0.5){

  # If data has multiple groups ----------------------------------------
  if("group" %in% colnames(window_fava)){
    window_long = window_fava %>%
      tidyr::pivot_longer(-c(.data$group, .data$FAVA, .data$window_index),
                          values_to = "index", names_to = "window_count")

    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = .data$index,
                                      y = .data$FAVA,
                                      color = .data$group,
                                      group = paste0(.data$window_index, .data$group)),
                         window_long, size = 1, alpha = alpha) +
      ggplot2::theme_bw() +
      ggplot2::ylab("FAVA") +
      ggplot2::xlab("Index") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  }else{
    # If data has only one group ----------------------------------------
    window_long = window_fava %>%
      tidyr::pivot_longer(-c(.data$FAVA, .data$window_index),
                          values_to = "index", names_to = "window_count")

    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = .data$index,
                                      y = .data$FAVA,
                                      group = paste0(.data$window_index)),
                         window_long, size = 1, alpha = alpha) +
      ggplot2::theme_bw() +
      ggplot2::ylab("FAVA") +
      ggplot2::xlab("Index") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  }

}



