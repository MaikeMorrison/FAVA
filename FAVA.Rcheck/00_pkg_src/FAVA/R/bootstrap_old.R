#' # bootstrap_fava ----------------------------------------
#' #' Generate and analyze bootstrap replicates of one or more relative abundance matrices
#' #'
#' #' Generates bootstrap replicate relative abundance matrices, computes fava (possibly weighted or normalized) for each bootstrap replicate, produces several plots of the bootstrap distributions of Fst (possibly weighted) and/or Fst/FstMax for each provided Q matrix, and runs two statistical tests comparing these bootstrap distributions. The tests comparing bootstrap distributions of Fst (possibly weighted) and/or Fst/FstMax facilitate statistical comparison of the variability in each of multiple Q matrices.
#' #'
#' #' @param matrices A dataframe, matrix, or array representing a relative abundance matrix or a (possibly named) list of arbitrarily many relative abundance matrices. For each relative abundance matrix, matrix rows represent samples and the last \code{K} columns contain relative abundances of distinct categories (when restricted to the last \code{K} columns, the rows must sum to 1). If the matrices are not named (e.g., \code{matrices = list(matrix1, matrix2)} instead of \code{matrices = list(A = matrix1, B = matrix2)}), the matrices will be numbered in the order they are provided in the list. If \code{matrices} is a single matrix, dataframe, or array and \code{group} is specified, the matrix will be split into multiple relative abundance matrices, one for each distinct value of the column \code{group}, which will each be analyzed separately.
#' #' @param n_replicates The number of bootstrap replicate matrices to generate for each provided relative abundance matrix.
#' #' @param K Optional; the number of categories in each provided relative abundance matrix, or a vector of such K values if the number of categories differs between matrices. If a single K is provided, each sample in every matrix must have \code{K} categories. If a vector of multiple K values is provided, \code{matrices} must be a list and the \eqn{i^{th}} entry of \code{K} must correspond to the \eqn{i^{th}} Q matrix in \code{matrices}. The default value of \code{K} is the number of columns in the matrix, the number of columns in the first matrix if a list is provided, or the number of columns minus 1 if \code{group} is specified but \code{K} is not.
#' #' @param seed Optional; a number to set as the random seed. Use if reproducibility of random results is desired.
#' #' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{matrices} is a single matrix containing multiple groups of samples you wish to compare.
#' #' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' #' @param w Optional; a vector of length \code{nrow(matrices)} (if \code{matrices} is a dataframe, matrix, or array and \code{group} is not provided) or a list of vectors, each corresponding to one group (if \code{matrices} is a dataframe, matrix, or array and \code{group} is provided) or one matrix (if \code{matrices} is a list of matrices). Each vector must contain non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/I, nrow(matrices))} where \code{I} is the number of rows in each matrix/group.
#' #' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(matrices))}.
#' #' @param normalized Optional; should normalized FAVA be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#' #' @param save_replicates Optional; should all of the bootstrap replicate matrices be included in the output? Default is \code{save_replicates = FALSE}; use \code{save_replicates = FALSE} to save memory when analyzing large datasets.
#' #'
#' #'
#' #' @return A named list containing the following entries:
#' #' \itemize{
#' #' \item \code{bootstrap_replicates}: A named list of lists. Each element is named for a relative abundance matrix provided in \code{matrices} and contains a list of \code{n_replicates} bootstrap replicates of the provided matrix. E.g., if \code{n_replicates = 100} and the first relative abundance matrix in \code{matrices} is named \code{A}, then the first element of \code{bootstrap_replicates}, \code{bootstrap_replicates$bootstrap_matrices_A}, is itself a list of 100 matrices, each representing a bootstrap replicate of matrix A.
#' #' \item \code{statistics}: A dataframe containing the desired statistic (FAVA, normalized FAVA, or weighted FAVA), computed for each bootstrap replicate matrix in \code{bootstrap_replicates}. The first column, titled \code{group}, is a factor indicating which provided relative abundance matrix the row corresponds to (the matrix name if \code{matrices} is a named list, or a number otherwise). The row names are of the form \code{stats_matrix.replicate} where \code{matrix} is the name of one of the provided relative abundance matrices (or the entry number if the list elements were not named) and replicate is the number of bootstrap replicate (rep takes values from 1 to \code{n_replicates}).
#' #' \item \code{plot_boxplot}: A ggplot2 box plot depicting the bootstrap distribution of FAVA for each matrix in \code{matrices}.
#' #' \item \code{plot_violin}: A ggplot2 violin plot depicting the bootstrap distribution of FAVA for each matrix in \code{matrices}.
#' #' \item \code{plot_ecdf}: A ggplot2 empirical cumulative distribution function plot depicting the bootstrap distribution of FAVA for each matrix in \code{matrices}.
#' #' \item \code{test_kruskal_wallis}: Results of a Kruskal-Wallis test performed on the bootstrap distributions of FAVA. This test is a non-parametric statistical test of whether all provided bootstrap distributions are identically distributed.
#' #' \item \code{test_pairwise_wilcox}: Results of a Wilcoxon rank-sum test performed on the bootstrap distributions of FAVA. This test is a non-parameteric statistical test of whether \emph{each pairwise combination} of provided bootstrap distributions is identically distributed. The result is a matrix of p-values whose entries correspond to each pair of relative abundance matrices.
#' #' }
#' #' @examples
#' #'
#' #' # Generate 100 bootstrap replicates for each of the
#' #' # 3 subjects included in the example data, xue_microbiome_sample.
#' #' # In this example, we:
#' #' # use the group argument so that each subject is analyzed,
#' #' # provide a species similarity matrix to incorporate phylogenetic information,
#' #' # and use the time argument so that samples are weighted by the amount of
#' #' # time between the samples before and after it.
#' #' boot_out = bootstrap_fava(matrices = xue_microbiome_sample,
#' #'                n_replicates = 100,
#' #'                seed = 1,
#' #'                group = "subject",
#' #'                time = "timepoint",
#' #'                S = xue_species_similarity,
#' #'                save_replicates = TRUE)
#' #'
#' #' # To look at a plot of the distribution of
#' #' # FAVA for each Q matrix:
#' #' boot_out$plot_violin
#' #'
#' #' # To determine if each of the 4 distibutions of
#' #' # FAVA is significantly different from
#' #' # each of the other distributions:
#' #' boot_out$test_pairwise_wilcox
#' #'
#' #'
#' #'
#' #' @importFrom dplyr %>%
#' #' @importFrom dplyr mutate
#' #' @importFrom dplyr filter
#' #' @importFrom rlang .data
#' #' @export
#' bootstrap_fava <- function(matrices,
#'                            n_replicates,
#'                            seed = NULL,
#'                            group = NULL,
#'                            time = NULL,
#'                            w = NULL,
#'                            K = NULL,
#'                            S = NULL,
#'                            normalized = FALSE,
#'                            save_replicates = TRUE) {
#'   . <- NULL # to please R command check
#'
#'   # set seed
#'   if (!is.null(seed)) {
#'     set.seed(seed)
#'   }
#'
#'   if(normalized == TRUE && ((!is.null(w)) | (!is.null(S)))){
#'     stop("FAVA can be either normalized or weighted, but not both. Please specify `normalized = TRUE` if you wish to compute normalized FAVA OR provide the weighting parameters w and/or S.")
#'   }
#'
#'   # If group provided, create a new matrices list
#'   # which converts the long single matrix to a list of matrices. -------------------------------
#'
#'   # If multiple grouping variables are provided, create a new grouping variable
#'   # that is the pairwise combination of all
#'   multiple_groups = FALSE
#'   relab_grouping_vars = NULL
#'   if(length(group) > 1){
#'     relab_grouping_vars = dplyr::select(matrices, dplyr::all_of(group))
#'
#'     relab_grouping_vars$grouping_var_multiple = apply(relab_grouping_vars, 1, paste, collapse = "_")
#'
#'     matrices = dplyr::mutate(matrices,
#'                                  grouping_var_multiple = relab_grouping_vars$grouping_var_multiple,
#'                                  .before = 1)
#'
#'     if(any(table(matrices$grouping_var_multiple)<2)){
#'       ignore = names(which(table(matrices$grouping_var_multiple)<2))
#'       warning("Only analyzing combinations of grouping variables with at least two samples. Ignoring the following combinations of grouping variables: ", paste(ignore, collapse = "  "))
#'       matrices = dplyr::filter(matrices,
#'                                    grouping_var_multiple %in%
#'                                      names(which(table(matrices$grouping_var_multiple) >= 2)))
#'     }
#'
#'     group = "grouping_var_multiple"
#'
#'     multiple_groups = TRUE
#'   }
#'
#'
#'
#'   if(!is.null(group)){
#'     # the true K, if K was not provided, will be nrow(Q) - 1 since one of the columns is group
#'     if(is.null(K)){ K = ncol(matrices) - (length(group)) - (!is.null(time))}
#'
#'     relab_matrix_clean = relab_checker(matrices, K = K, group = group, time = time)
#'
#'     groups <- unique(relab_matrix_clean$group)
#'
#'     matrix_list <- vector("list", length(groups))
#'     i = 1
#'     for(matrix in groups){
#'       matrix_list[[i]] <- dplyr::filter(dplyr::mutate(matrices, "group" = get(group)),
#'                                         group == matrix) %>%
#'         dplyr::select(-"group")
#'       i = i +1
#'     }
#'     names(matrix_list) <- groups
#'
#'     matrices <- matrix_list
#'   }
#'
#'
#'
#'   # Do computations if matrices = a single matrix ---------------------------------------
#'
#'   if ((is.data.frame(matrices) | is.array(matrices) | length(matrices) == 1) & is.null(group)) {
#'     if(is.null(K)){
#'       K = ncol(matrices) - (!is.null(time)) - (!is.null(group))
#'     }
#'
#'     n_matrix <- 1
#'
#'     names <- "Q"
#'
#'     # Clean Q matrix - isolate categories
#'     # matrices <- relab_checker(relab = matrices, K = K, group = group, time = time)$relab_matrix
#'
#'     bootstrap_matrices_Q <- list()
#'     matrix <- matrices
#'
#'     # Repeat rows if time or weight is provided
#'     if((!is.null(w)) | (!is.null(time))){
#'       matrix = relab_sample_weighter(relab = matrix, K = K, time = time, w = w)
#'     }
#'
#'     # Generate bootstrap data sets
#'     for (replicate in 1:n_replicates) {
#'       bootstrap_matrices_Q[[replicate]] <- matrix[purrr::rdunif(
#'         n = nrow(matrix),
#'         a = 1, b = nrow(matrix)
#'       ), ]
#'     }
#'
#'     # Compute statistics for these reps
#'     fava_Q <- lapply(
#'       X = bootstrap_matrices_Q,
#'       FUN = function(matrix) fava(relab_matrix = matrix, K = K, S = S, normalized = normalized)
#'     ) %>%
#'       unlist()
#'
#'     gini_simpson_mean_Q <- lapply(
#'       X = bootstrap_matrices_Q,
#'       FUN = function(matrix) gini_simpson_mean(relab_matrix = matrix, K = K, S = S)
#'     ) %>%
#'       unlist()
#'
#'     gini_simpson_pooled_Q <- lapply(
#'       X = bootstrap_matrices_Q,
#'       FUN = function(matrix) gini_simpson_pooled(relab_matrix = matrix, K = K, S = S)
#'     ) %>%
#'       unlist()
#'
#'
#'     # Do computations if matrices = a list ---------------------------------------------
#'   } else if (is.list(matrices)) {
#'     if(is.null(K)){
#'       K = ncol(matrices[[1]])
#'     }
#'     n_matrix <- length(matrices)
#'     K.list <- K
#'
#'     # List of names of matrices
#'     names <- if (sum(!is.na(names(matrices)))) {
#'       names(matrices)
#'     } else {
#'       1:n_matrix
#'     }
#'
#'     ## For each matrix: ##
#'     for (m in 1:n_matrix) {
#'       if(length(K.list)>1){
#'         K <- K.list[[m]]
#'       }
#'
#'       bs_list <- list()
#'       matrix <- matrices[[m]]
#'
#'       # Repeat rows if time or weight is provided
#'       if(((!is.null(w)) | (!is.null(time))) & nrow(matrix) > 2){
#'         matrix = relab_sample_weighter(relab = matrix, K = K, time = time, w = w)
#'       }
#'
#'
#'       # Generate bootstrap data sets
#'       for (replicate in 1:n_replicates) {
#'         bs_list[[replicate]] <- matrix[purrr::rdunif(
#'           n = nrow(matrix),
#'           a = 1, b = nrow(matrix)
#'         ), ]
#'       }
#'
#'       # Compute statistics for these reps
#'       fava <- lapply(
#'         X = bs_list,
#'         FUN = function(matrix) fava(relab_matrix = matrix, K = K, S = S, normalized = normalized)
#'       ) %>%
#'         unlist()
#'
#'       gini_simpson_mean <- lapply(
#'         X = bs_list,
#'         FUN = function(matrix) gini_simpson_mean(relab_matrix = matrix, K = K, S = S)
#'       ) %>%
#'         unlist()
#'
#'       gini_simpson_pooled <- lapply(
#'         X = bs_list,
#'         FUN = function(matrix) gini_simpson_pooled(relab_matrix = matrix, K = K, S = S)
#'       ) %>%
#'         unlist()
#'
#'
#'       # Name this dataset, based on the name of the matrices in the list or the entry number
#'       assign(paste0("fava_", names[m]), fava, pos = -1)
#'       assign(paste0("gini_simpson_mean_", names[m]), gini_simpson_mean, pos = -1)
#'       assign(paste0("gini_simpson_pooled_", names[m]), gini_simpson_pooled, pos = -1)
#'       assign(paste0("bootstrap_matrices_", names[m]), bs_list, pos = -1)
#'     }
#'   } else {
#'     stop("Error: The entry `matrices` must be a data frame, matrix, or array, or a list of these objects.")
#'   }
#'
#'   # Make a dataset with all matrices' statistics:
#'   all_stats <- data.frame(
#'     group =
#'       names %>%
#'       lapply(function(name) rep(name, n_replicates)) %>%
#'       unlist(),
#'     FAVA = mget(paste0("fava_", names)) %>%
#'       do.call(what = c, args = .),
#'     gini_simpson_mean = mget(paste0("gini_simpson_mean_", names)) %>%
#'       do.call(what = c, args = .),
#'     gini_simpson_pooled = mget(paste0("gini_simpson_pooled_", names)) %>%
#'       do.call(what = c, args = .)
#'   )
#'
#'   stat_name = ifelse(normalized == TRUE, "Normalized FAVA",
#'                      ifelse(any(!sapply(list(time, w, S), is.null)),
#'                             "Weighted FAVA", "FAVA"))
#'
#'   if(multiple_groups){
#'     merge_group_names = dplyr::distinct(relab_grouping_vars)
#'     merge_group_names$group = merge_group_names$grouping_var_multiple
#'     merge_group_names <- dplyr::select(merge_group_names,
#'                                        -grouping_var_multiple)
#'
#'     all_stats = dplyr::right_join(merge_group_names,
#'                                   all_stats,
#'                                   .before = 1)
#'
#'   }
#'
#'   all_stats$group <- factor(all_stats$group, levels = unique(all_stats$group))
#'
#'   plot_ecdf <- ggplot2::ggplot(data = all_stats) +
#'     ggplot2::stat_ecdf(ggplot2::aes(x = .data$FAVA, color = .data$group)) +
#'     ggplot2::xlab(stat_name) +
#'     ggplot2::ylab("Cumulative Probability") +
#'     ggplot2::xlim(0, 1) +
#'     ggplot2::theme_bw() +
#'     ggplot2::scale_color_viridis_d()
#'
#'   plot_boxplot <- ggplot2::ggplot(
#'     data = all_stats,
#'     ggplot2::aes(x = .data$group, y = .data$FAVA)
#'   ) +
#'     ggplot2::geom_boxplot() +
#'     ggplot2::ylab(stat_name) +
#'     ggplot2::xlab("") +
#'     ggplot2::theme_bw()
#'
#'   plot_violin <- ggplot2::ggplot(
#'     data = all_stats,
#'     ggplot2::aes(
#'       x = .data$group,
#'       y = round(.data$FAVA, 5)
#'     )
#'   ) +
#'     ggplot2::geom_violin(scale = "width") +
#'     ggplot2::geom_boxplot(width = 0.3) +
#'     ggplot2::ylab(stat_name) +
#'     ggplot2::xlab("") +
#'     ggplot2::theme_bw() +
#'     ggplot2::geom_jitter(alpha = 0.4, width = 0.05)
#'
#'   if (is.data.frame(matrices) | is.array(matrices)) {
#'     test_kruskal_wallis <- "This statistical test can only be performed if a list of matrices is provided."
#'
#'     test_pairwise_wilcox <- "This statistical test can only be performed if a list of matrices is provided."
#'   } else {
#'     test_kruskal_wallis <- stats::kruskal.test(all_stats$FAVA ~ all_stats$group)
#'
#'     test_pairwise_wilcox <- stats::pairwise.wilcox.test(
#'       x = all_stats$FAVA,
#'       g = all_stats$group,
#'       paired = FALSE
#'     )
#'   }
#'
#'   bootstrap_replicates <- mget(paste0(
#'     "bootstrap_matrices_",
#'     names
#'   ))
#'
#'   if(!save_replicates){
#'     bootstrap_replicates = NULL
#'   }
#'
#'   # Set name to match group name, if applicable
#'   if((!multiple_groups) & (!is.null(group))){
#'     colnames(all_stats)[1] = group
#'   }
#'
#'   return(list(
#'     bootstrap_replicates = bootstrap_replicates,
#'     statistics = all_stats,
#'     plot_boxplot = plot_boxplot,
#'     plot_violin = plot_violin,
#'     plot_ecdf = plot_ecdf,
#'     test_kruskal_wallis = test_kruskal_wallis,
#'     test_pairwise_wilcox = test_pairwise_wilcox
#'   ))
#' }
#'
