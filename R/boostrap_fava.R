# bootstrap_fava ----------------------------------------
#' Generate and analyze bootstrap replicates of one or more relative abundance matrices
#'
#' Generates bootstrap replicate relative abundance matrices, computes fava (possibly weighted or normalized) for each bootstrap replicate, produces several plots of the bootstrap distributions of Fst (possibly weighted) and/or Fst/FstMax for each provided Q matrix, and runs two statistical tests comparing these bootstrap distributions. The tests comparing bootstrap distributions of Fst (possibly weighted) and/or Fst/FstMax facilitate statistical comparison of the variability in each of multiple Q matrices.
#'
#' @param matrices A dataframe, matrix, or array representing a relative abundance matrix or a (possibly named) list of arbitrarily many relative abundance matrices. For each relative abundance matrix, matrix rows represent samples and the last \code{K} columns contain relative abundances of distinct categories (when restricted to the last \code{K} columns, the rows must sum to 1). If the matrices are not named (e.g., \code{matrices = list(matrix1, matrix2)} instead of \code{matrices = list(A = matrix1, B = matrix2)}), the matrices will be numbered in the order they are provided in the list. If \code{matrices} is a single matrix, dataframe, or array and \code{group} is specified, the matrix will be split into multiple relative abundance matrices, one for each distinct value of the column \code{group}, which will each be analyzed separately.
#' @param n_replicates The number of bootstrap replicate matrices to generate for each provided relative abundance matrix.
#' @param K Optional; the number of categories in each provided relative abundance matrix, or a vector of such K values if the number of categories differs between matrices. If a single K is provided, each sample in every matrix must have \code{K} categories. If a vector of multiple K values is provided, \code{matrices} must be a list and the \eqn{i^{th}} entry of \code{K} must correspond to the \eqn{i^{th}} Q matrix in \code{matrices}. The default value of \code{K} is the number of columns in the matrix, the number of columns in the first matrix if a list is provided, or the number of columns minus 1 if \code{group} is specified but \code{K} is not.
#' @param seed Optional; a number to set as the random seed. Use if reproducibility of random results is desired.
#' @param group Optional; a string specifying the name of the column that describes which group each row (sample) belongs to. Use if \code{matrices} is a single matrix containing multiple groups of samples you wish to compare.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param w Optional; a vector of length \code{nrow(matrices)} (if \code{matrices} is a dataframe, matrix, or array and \code{group} is not provided) or a list of vectors, each corresponding to one group (if \code{matrices} is a dataframe, matrix, or array and \code{group} is provided) or one matrix (if \code{matrices} is a list of matrices). Each vector must contain non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/I, nrow(relab_matrix))} where \code{I} is the number of rows in each matrix/group.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param normalized Optional; should normalized FAVA be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#' @param save_replicates Optional; should all of the bootstrap replicate matrices be included in the output? Default is \code{save_replicates = FALSE}; use \code{save_replicates = FALSE} to save memory when analyzing large datasets.
#'
#'
#' @return A named list containing the following entries:
#' \itemize{
#' \item \code{bootstrap_replicates}: A named list of lists. Each element is named for a relative abundance matrix provided in \code{matrices} and contains a list of \code{n_replicates} bootstrap replicates of the provided matrix. E.g., if \code{n_replicates = 100} and the first relative abundance matrix in \code{matrices} is named \code{A}, then the first element of \code{bootstrap_replicates}, \code{bootstrap_replicates$bootstrap_matrices_A}, is itself a list of 100 matrices, each representing a bootstrap replicate of matrix A.
#' \item \code{statistics}: A dataframe containing the desired statistic (FAVA, normalized FAVA, or weighted FAVA), computed for each bootstrap replicate matrix in \code{bootstrap_replicates}. The first column, titled \code{Matrix}, is a factor indicating which provided relative abundance matrix the row corresponds to (the matrix name if \code{matrices} is a named list, or a number otherwise). The row names are of the form \code{stats_matrix.replicate} where \code{matrix} is the name of one of the provided relative abundance matrices (or the entry number if the list elements were not named) and replicate is the number of bootstrap replicate (rep takes values from 1 to \code{n_replicates}).
#' \item \code{plot_boxplot}: A ggplot2 box plot depicting the bootstrap distribution of FAVA for each matrix in \code{matrices}.
#' \item \code{plot_violin}: A ggplot2 violin plot depicting the bootstrap distribution of FAVA for each matrix in \code{matrices}.
#' \item \code{plot_ecdf}: A ggplot2 empirical cumulative distribution function plot depicting the bootstrap distribution of FAVA for each matrix in \code{matrices}.
#' \item \code{test_kruskal_wallis}: Results of a Kruskal-Wallis test performed on the bootstrap distributions of FAVA. This test is a non-parametric statistical test of whether all provided bootstrap distributions are identically distributed.
#' \item \code{test_pairwise_wilcox}: Results of a Wilcoxon rank-sum test performed on the bootstrap distributions of FAVA. This test is a non-parameteric statistical test of whether \emph{each pairwise combination} of provided bootstrap distributions is identically distributed. The result is a matrix of p-values whose entries correspond to each pair of relative abundance matrices.
#' }
#' @examples
#
# # Generate 100 bootstrap replicates for each of the
# bootstrap_2 <- Q_bootstrap(matrices = Q_4,
#                            n_replicates = 100,
#                            K = 2,
#                            seed = 1,
#                            group = "Pop")
#
# # To look at a plot of the distribution of
# # Fst/FstMax for each Q matrix:
# bootstrap_2$plot_violin
#
# # To determine if each of the 4 distibutions of
# # Fst/FstMax is significantly different from
# # each of the other distributions:
# bootstrap_2$test_pairwise_wilcox
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
Q_bootstrap <- function(matrices, n_replicates, K, seed, group,
                        time, w, S, normalized = FALSE, save_replicates = TRUE) {
  . <- NULL # to please R command check

  # set seed
  if (!missing(seed)) {
    set.seed(seed)
  }

  if(!missing(normalized) & (!missing(w) | !missing(S))){
    stop("FAVA can be either normalized or weighted, but not both. Please specify `normalized = TRUE` if you wish to compute normalized FAVA OR provide the weighting parameters w and/or S.")
  }

  # If group provided, create a new matrices list
  # which converts the long single matrix to a list of matrices.

  if(!missing(group)){
    # the true K, if K was not provided, will be nrow(Q) - 1 since one of the columns is group
    if(missing(K)){ K = ncol(matrices) - 1 - !missing(timepoint)}

    groups <- unique(matrices %>% dplyr::select(dplyr::all_of(group)) %>% unlist)

    matrix_list <- vector("list", length(groups))
    i = 1
    for(matrix in groups){
      matrix_list[[i]] <- dplyr::filter(dplyr::mutate(matrices, "group" = get(group)),
                                        group == matrix) %>%
        select(-"group")
      i = i +1
    }
    names(matrix_list) <- groups

    matrices <- matrix_list
  }


  # Do computations if matrices = a single matrix ---------------------------------------

  if ((is.data.frame(matrices) | is.array(matrices) | length(matrices) == 1) & missing(group)) {
    if(missing(K)){
      K = ncol(matrices) - !missing(time)
    }

    n_matrix <- 1

    names <- "Q"

    # Clean Q matrix - isolate categories
    matrices <- Q_checker(Q = matrices, K = K)

    bootstrap_matrices_Q <- list()
    matrix <- matrices
    # Generate bootstrap data sets
    for (replicate in 1:n_replicates) {
      bootstrap_matrices_Q[[replicate]] <- matrix[purrr::rdunif(
        n = nrow(matrix),
        a = 1, b = nrow(matrix)
      ), ]
    }

    # Compute statistics for these reps
    stats_Q <- lapply(
      X = bootstrap_matrices_Q,
      FUN = function(matrix) Q_stat(Q = matrix, K = ncol(matrix))
    ) %>%
      unlist() %>%
      matrix(ncol = 3, byrow = TRUE) %>%
      data.frame() %>%
      `colnames<-`(c("Fst", "FstMax", "ratio"))


    # Do computations if matrices = a list ---------------------------------------------
  } else if (is.list(matrices)) {
    if(missing(K)){
      K = ncol(matrices[[1]])
    }
    n_matrix <- length(matrices)
    K.list <- K

    # List of names of matrices
    names <- if (sum(!is.na(names(matrices)))) {
      names(matrices)
    } else {
      1:n_matrix
    }

    ## For each matrix: ##
    for (m in 1:n_matrix) {
      if(length(K.list)>1){
        K <- K.list[[m]]
      }

      bs_list <- list()
      matrix <- matrices[[m]]

      # Check format of matrix
      matrix <- Q_checker(Q = matrix, K = K, rep = m)

      # Generate bootstrap data sets
      for (replicate in 1:n_replicates) {
        bs_list[[replicate]] <- matrix[purrr::rdunif(
          n = nrow(matrix),
          a = 1, b = nrow(matrix)
        ), ]
      }

      # Compute statistics for these reps
      stats <- lapply(
        X = bs_list,
        FUN = function(matrix) Q_stat(Q = matrix, K = ncol(matrix))
      ) %>%
        unlist() %>%
        matrix(ncol = 3, byrow = TRUE) %>%
        data.frame() %>%
        `colnames<-`(c("Fst", "FstMax", "ratio"))

      # Sometimes, for values of Fst very close to 0 (i.e. order 10^-6, 10^-7), the
      # value of Fst ends up negative due to precision errors.
      # Find matrices for which this is the case, and replace them and their statistics

      while (sum(stats$ratio < 0)) {
        # Which bootstrap matrices have negative values for Fst/FstMax?
        negatives <- which(stats$ratio < 0)

        # Replace those bootstrap replicates with new, random bootstrap replicates
        bs_list[negatives] <- lapply(
          X = 1:length(negatives),
          FUN = function(x) {
            matrix[purrr::rdunif(
              n = nrow(matrix),
              a = 1, b = nrow(matrix)
            ), ]
          }
        )

        # Replace the corresponding entries of the statistics matrix
        stats[negatives, ] <- lapply(
          X = bs_list[negatives],
          FUN = function(matrix) {
            Q_stat(Q = matrix, K = ncol(matrix))
          }
        ) %>%
          unlist() %>%
          matrix(ncol = 3, byrow = TRUE) %>%
          data.frame()
      } # repeat this until there are no more errors

      # Name this dataset, based on the name of the matrices in the list or the entry number
      assign(paste0("stats_", names[m]), stats, pos = -1)
      assign(paste0("bootstrap_matrices_", names[m]), bs_list, pos = -1)
    }
  } else {
    stop("Error: The entry `matrices` must be a data frame, matrix, or array, or a list of these objects.")
  }

  # Make a dataset with all matrices' statistics:
  all_stats <- cbind(
    Matrix =
      names %>%
      lapply(function(name) rep(name, n_replicates)) %>%
      unlist(),
    mget(paste0("stats_", names)) %>%
      do.call(what = rbind, args = .) %>%
      rbind()
  )

  all_stats$Matrix <- factor(all_stats$Matrix, levels = unique(all_stats$Matrix))

  plot_ecdf <- ggplot2::ggplot(data = all_stats) +
    ggplot2::stat_ecdf(ggplot2::aes(x = .data$ratio, color = .data$Matrix)) +
    ggplot2::xlab(expression(F[ST]/F[ST]^{max})) +
    ggplot2::ylab("Cumulative Probability") +
    ggplot2::xlim(0, 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_viridis_d()

  plot_boxplot <- ggplot2::ggplot(
    data = all_stats,
    ggplot2::aes(x = .data$Matrix, y = .data$ratio)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::ylab(expression(F[ST]/F[ST]^{max})) +
    ggplot2::xlab("") +
    ggplot2::theme_bw()

  plot_violin <- ggplot2::ggplot(
    data = all_stats,
    ggplot2::aes(
      x = .data$Matrix,
      y = round(.data$ratio, 5)
    )
  ) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::geom_boxplot(width = 0.3) +
    ggplot2::ylab(expression(F[ST]/F[ST]^{max})) +
    ggplot2::xlab("") +
    ggplot2::theme_bw()

  if (is.data.frame(matrices) | is.array(matrices)) {
    test_kruskal_wallis <- "This statistical test can only be performed if a list of matrices is provided."

    test_pairwise_wilcox <- "This statistical test can only be performed if a list of matrices is provided."
  } else {
    test_kruskal_wallis <- stats::kruskal.test(all_stats$ratio ~ all_stats$Matrix)

    test_pairwise_wilcox <- stats::pairwise.wilcox.test(
      x = all_stats$ratio,
      g = all_stats$Matrix,
      paired = FALSE
    )
  }

  bootstrap_replicates <- mget(paste0(
    "bootstrap_matrices_",
    names
  ))

  return(list(
    bootstrap_replicates = bootstrap_replicates,
    statistics = all_stats,
    plot_boxplot = plot_boxplot,
    plot_violin = plot_violin,
    plot_ecdf = plot_ecdf,
    test_kruskal_wallis = test_kruskal_wallis,
    test_pairwise_wilcox = test_pairwise_wilcox
  ))
}
