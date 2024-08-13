# bootstrap_fava
#' Statistically compare FAVA values between pairs of relative abundance matrices.
#'
#' \code{bootstrap_fava} uses bootstrapping to statistically compare FAVA values between pairs of relative abundance matrices. \code{bootstrap_fava} takes the same options as \code{fava}, so, as with \code{fava}, you can separately analyze multiple populations or groups of samples (specify \code{group}), and account for similarity among categories (specify \code{S}) or uneven weighting of rows (specify \code{w} or \code{time}). \code{bootstrap_fava} follows the bootstrapping procedure defined by Efron and Tibshirani (1993). Details on the bootstrapping procedure are available in the Methods section of the accompanying paper.
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a category, and each entry represents the abundance of that category in the sample. If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param n_replicates The number of bootstrap replicate matrices to generate. Default is `n_replicates = 1000`.
#' @param group A string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param normalized Optional; should normalized FAVA be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#' @param seed Optional; an integer to be used as a random seed for the simulations.
# #' @param save_replicates Optional; should all of the bootstrap replicate matrices be included in the output? Default is \code{save_replicates = FALSE};  \code{save_replicates = FALSE} save memory when analyzing large datasets.
#' @param alternative Optional; do you want to do a one- or two.sided test? Default is \code{alternative = "two.sided"}. If you wish to do a one-sided test, specify either \code{alternative = "lesser"} or \code{alternative = "greater"}.
#' @returns A named list containing the following entries:
#' \itemize{
#' \item \code{p_values}: The probability of observing the observed difference in variability between each pair of groups if there were no difference between groups. Computed as the fraction of bootstrap differences greater than or equal to the observed difference. Depends on what \code{alternative} is specified ("greater", "lesser", or "two.sided").
#' \item \code{bootstrap_distribution_plot}: The distribution of bootstrap replicate differences in each variability value. The observed differences are shown in red. The further the red points are from 0, the more significant the statistical difference between groups.
#' \item \code{observed_stats}: The observed diversity statistics for the groups.
#' \item \code{bootstrap_stats}: The bootstrap replicate diversity statistics for the groups.
#' \item \code{bootstrap_replicates}: The bootstrap replicate matrices, reported only if  \code{save_replicates = TRUE}.}
#' @examples
#' # Statistically compare values of FAVA between
#' # subjects in the xue_microbiome_sample data:
#'
#'  boot_out = bootstrap_fava(relab_matrix = xue_microbiome_sample,
#'                n_replicates = 100,
#'                seed = 1,
#'                group = "subject",
#'                K = 524,
#'                S = xue_species_similarity)
#'
#' # Table of P-values comparing values of FAVA between group 1 and group 2:
#'  boot_out$P_values
#'
#'  # Plots of the bootstrap distributions of differences in FAVA between each pair of matrices,
#'  # and how the true observed differences (red dots) compare to the distribution.
#'  boot_out$bootstrap_distribution_plot
#' @export
bootstrap_fava <- function(relab_matrix,
                           n_replicates = 1000,
                           group,
                           K = NULL,
                           S = NULL,
                           w = NULL,
                           time = NULL,
                           normalized = FALSE,
                           seed = NULL,
                           # save_replicates = FALSE,
                           alternative = "two.sided"){

  # To appease R cmd check
  P_value <- P_value_numeric <- Comparison <- Difference <- combn <- . <- NULL

  if((!is.null(time)) && (!is.null(w))){
    stop("Please specify either time or w, but not both.")
  }

  # Set random seed (optional)
  if(!is.null(seed)){
    set.seed(seed)
  }

  # If multiple grouping variables are provided, make a new grouping column
  multiple_groups = FALSE
  if(length(group)>1){
    multiple_groups = TRUE
    relab_matrix = dplyr::mutate(relab_matrix,
                                 group =  apply( relab_matrix[ , group ] , 1 , paste , collapse = "_" ),
                                 .before = 1)
    group_multiple = group
    group = "group"

    group_table = dplyr::distinct(dplyr::select(relab_matrix, dplyr::all_of(c("group", group_multiple))))

    relab_matrix = relab_matrix %>%
      dplyr::select(-dplyr::all_of(group_multiple))
  }

  # Repeat rows if time or w are provided
  # relab_weighted = relab_sample_weighter(relab = relab_matrix, K = K, time = time, w = w, group = group)

  relab_check_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

  relab_clean = relab_check_out$relab_matrix
  groups = relab_check_out$group
  times = relab_check_out$time
  K = ncol(relab_clean)
  if(is.null(w) & !is.null(time)){
    w = time_weights(times = times, group = groups)
  }

  # How many groups are there in the data? Do we need to do multiple pairwise comparisons?
  if(length(unique(groups)) < 2){
    stop(paste0("bootstrap_fava must be provided with multiple groups to compare. The grouping column '",
                group,
                "' contains only the following group: '",
                groups, "'\n" ))
  }
  if(length(unique(groups)) == 2){
    bootstrap_list = pairwise_comparison(group_pair = unique(groups), relab_clean = relab_clean,
                              n_replicates = n_replicates, groups = groups, K = K, S = S, w = w, normalized = normalized, alternative = alternative)


    p_values = data.frame(Comparison = bootstrap_list$comparison, P_value = bootstrap_list$P_value) %>%
      dplyr::mutate(P_value_numeric = as.numeric(P_value),
                    P_value = ifelse(P_value_numeric==0, paste0("<", 1/n_replicates), P_value_numeric))

    # 2 - observed_difference
    observed_difference = data.frame(Comparison = bootstrap_list$comparison, Difference = bootstrap_list$observed_difference)



    # 3 - bootstrap_difference
    bootstrap_difference = data.frame(Comparison = bootstrap_list$comparison, Difference = bootstrap_list$bootstrap_difference)

    # 4 - bootstrap_plot
    bootstrap_plot = ggplot2::ggplot(data = bootstrap_difference, mapping = ggplot2::aes(x = Comparison, y = Difference)) +
      ggplot2::geom_violin(fill = "grey", color = NA, alpha = 0.6) +
      ggplot2::geom_boxplot(alpha = 0, width = 0.3) +
      ggplot2::geom_jitter(alpha = 0.3, width = 0.05, size = 2) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Difference in FAVA") +
      ggplot2::geom_point(data = observed_difference, color = "red", size = 5)


    return(list(P_values = p_values,
                bootstrap_distribution_plot = bootstrap_plot,
                observed_difference = observed_difference,
                bootstrap_difference = bootstrap_difference))
    #
    #
    # out$P_values = out$P_values %>% data.frame() %>% t() %>%
    #   cbind(group_1 = groups[[1]], group_2 = groups[[2]],.)
    # rownames( out$P_values ) = NULL
    #
    # if(multiple_groups){
    #   out$observed_stats = dplyr::left_join(group_table, out$observed_stats)
    #   out$bootstrap_stats = dplyr::left_join(group_table, out$bootstrap_stats)
    # }
    # return(out)
  }else{
    # Make a list of all unique pairs of groups
    group_pairs = t(combn(unique(groups), 2))

    # Do the bootstrap comparison procedure for each group:
    bootstrap_list = list()
    for(pair in 1:nrow(group_pairs)){
      group_pair = group_pairs[pair,]
      # relab_pair = relab_matrix[relab_matrix[[group]] %in% group_pair, ]

      bootstrap_list[[pair]] = pairwise_comparison(group_pair = group_pair, relab_clean = relab_clean, n_replicates = n_replicates,
                                                   groups = groups, K = K, S = S, w = w,
                                                   normalized = normalized, alternative = alternative)


    }

    # Combine all elements from each category

    # 1 - P-values
    p_values = lapply(bootstrap_list, function(pair) c(pair$comparison, pair$P_value)) %>%
      do.call(rbind, .) %>%
      data.frame() %>%
      `colnames<-`(c("Comparison", "P_value")) %>%
      dplyr::mutate(P_value_numeric = as.numeric(P_value),
                    P_value = ifelse(P_value_numeric==0, paste0("<", 1/n_replicates), P_value_numeric))

      # cbind(data.frame(t(combn(groups, 2))) %>% `colnames<-`(c("group_1", "group_2")),
      #                lapply(bootstrap_list, function(list) list$P_values) %>%
      #                  do.call(rbind, .))

    # # 2 - bootstrap_distribution_plot
    # bootstrap_distribution_plots = lapply(bootstrap_list, function(list)
    #   list$bootstrap_distribution_plot) %>% `names<-`(apply(group_pairs, 1, paste0, collapse = "--"))

    # 2 - observed_difference
    observed_difference = lapply(bootstrap_list, function(list) c(list$comparison, list$observed_difference)) %>%
      do.call(rbind, .) %>% data.frame %>% `colnames<-`(c("Comparison", "Difference")) %>%
      dplyr::mutate(Difference = as.numeric(Difference))
    observed_difference$Comparison = stringr::str_replace_all(observed_difference$Comparison, ' - ', " -\n")
    # if(multiple_groups){ observed_difference = left_join(group_table, observed_difference) }


    # 3 - bootstrap_difference
    bootstrap_difference = lapply(bootstrap_list, function(list) list$bootstrap_difference) %>%
      do.call(cbind, .) %>% `colnames<-`(observed_difference$Comparison) %>%
      data.frame() %>%
      tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Comparison", values_to = "Difference")
    bootstrap_difference$Comparison = stringr::str_replace_all(bootstrap_difference$Comparison, '\\.\\.\\.', " -\n")

    # if(multiple_groups){ bootstrap_difference = left_join(group_table, bootstrap_difference) }

    # 4 - bootstrap_plot
    bootstrap_plot = ggplot2::ggplot(data = bootstrap_difference, mapping = ggplot2::aes(x = Comparison, y = Difference)) +
      ggplot2::geom_violin(fill = "grey", color = NA, alpha = 0.6) +
      ggplot2::geom_boxplot(alpha = 0, width = 0.3) +
      ggplot2::geom_jitter(alpha = 0.3, width = 0.05, size = 2) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Difference in FAVA") +
      ggplot2::geom_point(data = observed_difference, color = "red", size = 5)

    # # 5 - bootstrap_replicates
    #
    # bootstrap_replicates = lapply(bootstrap_list, function(list) list$bootstrap_replicates) %>%
    #   do.call(rbind, .)

    return(list(P_values = p_values,
                bootstrap_distribution_plot = bootstrap_plot,
                observed_difference = observed_difference,
                bootstrap_difference = bootstrap_difference))

  }
}



pairwise_comparison <- function(group_pair,
                                relab_clean,
                                n_replicates,
                                groups,

                                K,
                                S,
                                w,
                                # time,
                                normalized,

                                # save_replicates,
                                alternative){
  # To appease R cmd check
  . <- NULL

  # Confirm there are only two groups provided
  if(length(group_pair) != 2){
    stop(paste0("There must be exactly 2 groups. There are ", length(groups), " groups in the provided relab_pair matrix."))
  }



  # Split the data into the two groups
  A = relab_clean[groups == group_pair[[1]],]
  B = relab_clean[groups == group_pair[[2]],]
  pooled = rbind(A, B)


  m = nrow(A)
  n = nrow(B)
  N = m+n
  if(!is.null(w)){
    wA = w[groups == group_pair[[1]]]
    wB = w[groups == group_pair[[2]]]
    wPooled = c(wA, wB)
  }else{wPooled = NULL; wA = NULL; wB = NULL}

  # Generate bootstrap replicates of the pooled groups
  rep_list = list()
  for(rep in 1:n_replicates){
    a_samp = sample(1:N, m, replace = TRUE)
    b_samp = sample(1:N, n, replace = TRUE)
    # A_rep = pooled[sample(1:N, m, replace = TRUE),]
    # B_rep = pooled[sample(1:N, n, replace = TRUE),]

    rep_list[[rep]] = list(A = a_samp, B = b_samp)
  }

  # Compute statistics for each replicate
  bootstrap_stats = lapply(rep_list,
                           function(rep){
                             # Compute for A
                             c(fava(relab_matrix = pooled[rep$A,], K = K, S = S, normalized = normalized,
                                    w = `if`(is.null(NULL), NULL, wPooled[rep$A]/sum(wPooled[rep$A]))),
                               # Compute for B
                               fava(relab_matrix = pooled[rep$B,], K = K, S = S, normalized = normalized,
                                    w = `if`(is.null(NULL), NULL, wPooled[rep$B]/sum(wPooled[rep$B]))))
                           }

  ) %>%
    do.call(rbind, .) %>%
    `colnames<-`(group_pair) %>%
    data.frame()

  bootstrap_stats$Difference = bootstrap_stats[,group_pair[[1]]] - bootstrap_stats[,group_pair[[2]]]

  # Compute the difference between the two original populations
  observed_difference = fava(relab_matrix = A, K = K, S = S, normalized = normalized, w = wA) -
    fava(relab_matrix = B, K = K, S = S, normalized = normalized, w = wB)

  diff_label = paste0( #"Difference (",
    group_pair[[1]], " - ", group_pair[[2]]#, ")"
  )

  # COMPUTE P-VALUES
  if(alternative == "greater"){
    p_val = mean( bootstrap_stats$Difference >= observed_difference )
  } else if(alternative == "less"){
    p_val = mean( bootstrap_stats$Difference <= observed_difference )
  } else if(alternative == "two.sided"){
    # two.sided p-value is from pg 212, eqn 15.26 of Efron and Tibshirani book
    p_val =  mean(abs(bootstrap_stats$Difference) >= abs(observed_difference))
  } else {
    stop("Valid options for alternative are 'greater', 'less', and 'two.sided'.")
  }

  # Replace P-values of 0 with less than 1 over the number of replicates
  # p_values[p_values ==  0] = paste0("<", round(1/n_replicates, 10))



  # PLOT
  # boot_diff_long <- bootstrap_difference %>%
  #   tidyr::pivot_longer(cols = 1:3, names_to = "Statistic", values_to = "Difference")
  #
  # obs_diff_long <- data.frame(Statistic = names(observed_difference),
  #                             Difference = observed_difference)

  # plot =
  #   ggplot2::ggplot(data = bootstrap_difference, ggplot2::aes(x = diff_label, y = FAVA)) +
  #   ggplot2::geom_violin(fill = "grey", color = NA, width = 1, alpha = 0.75) +
  #   ggplot2::geom_boxplot(width = 0.5, alpha = 0) +
  #   ggplot2::geom_hline(yintercept = 0) +
  #   ggplot2::geom_jitter(width = 0.05, alpha = 0.5) +
  #   ggplot2::geom_point(ggplot2::aes(y = observed_difference, x = diff_label), color = "red", size = 5, alpha = 0.7) +
  #   ggplot2::xlab(diff_label) +
  #   ggplot2::theme_bw() +
  #   ggplot2::ylab("Difference in FAVA") +
  #   ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10))
  # plot
  #   ggplot2::geom_histogram(fill = "grey", color = NA, width = 1, alpha = 0.75) +
  #   ggplot2::geom_vline(xintercept = 0, linewidth = 2) +
  #   ggplot2::geom_vline(xintercept  = observed_difference, color = "red", linewidth = 2) +
  #   ggplot2::xlab(diff_label) +
  #   ggplot2::theme_bw() #+
    # ggplot2::scale_x_discrete(labels = c("FAVA" = "Across-sample\nheterogeneity",
    #                                      "gini_simpson_mean" = "Mean within-\nsample diversity",
    #                                      "gini_simpson_pooled" =  "Diversity of\npooled samples"))
  # plot
  # bootstrap_replicates = rep_list

  # if(!save_replicates){
  #   bootstrap_replicates = NULL
  # }

  # bootstrap_stats = cbind(group = c(rep(groups[[1]], n_replicates),
  #                                   rep(groups[[2]], n_replicates)),
  #                         rbind(bootstrap_stats_A, bootstrap_stats_B))

  # observed_stats = cbind(group = c(groups[[1]], groups[[2]]),
  #                        rbind(observed_stats_A, observed_stats_B) %>%
  #                          apply(c(1,2), as.numeric) %>% `row.names<-`(NULL)) %>%
  #   data.frame()

  return(list(comparison = diff_label,
              P_value = p_val,
              # bootstrap_distribution_plot = plot,
              observed_difference = observed_difference,
              bootstrap_difference = bootstrap_stats$Difference
              # bootstrap_replicates = bootstrap_replicates
              ))
}











