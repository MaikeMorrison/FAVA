
# Function to process relab matrix and return correct grouping columns, w, K, S, etc.

process_relab <- function(relab_matrix,
                          K = NULL,
                          S = NULL,
                          w = NULL,
                          time = NULL,
                          group = NULL){
  # To appease R cmd check
  grouping_var_multiple <- NULL


  # Define K if not provided
  if(is.null(K)){
    K = ncol(relab_matrix) - (length(group)) - (!is.null(time))
  }

  # Define S if not provided
  if(is.null(S)){
    S = diag(K)
  }

  # Modify relab_matrix if there are multiple groups
  multiple_groups = FALSE
  relab_grouping_vars = NULL
  if(length(group) > 1){
    relab_grouping_vars = dplyr::select(relab_matrix, dplyr::all_of(group))

    relab_grouping_vars$grouping_var_multiple = apply(relab_grouping_vars, 1,
                                                      paste, collapse = "_")

    relab_matrix = dplyr::mutate(relab_matrix,
                                 grouping_var_multiple =
                                   relab_grouping_vars$grouping_var_multiple,
                                 .before = 1)

    if(any(table(relab_matrix$grouping_var_multiple)<2)){
      ignore = names(which(table(relab_matrix$grouping_var_multiple)<2))
      warning("Only analyzing combinations of grouping variables with at least two samples. Ignoring the following combinations of grouping variables: ",
              paste(ignore, collapse = "  "))
      relab_matrix = dplyr::filter(relab_matrix,
                                   grouping_var_multiple %in%
                                     names(which(table(relab_matrix$grouping_var_multiple) >= 2)))
    }

    group = "grouping_var_multiple"

    multiple_groups = TRUE
  }

  # Arrange data by group and time if not already
  # Need to re-order w to match, if provided
  if(!is.null(time)){
    if(!is.null(group)){
      relab_matrix$index = 1:nrow(relab_matrix)

      relab_matrix = dplyr::arrange(.data = relab_matrix, .data[[group]], .data[[time]])

      w = w[relab_matrix$index]
      relab_matrix$index = NULL
    } else{
      relab_matrix$index = 1:nrow(relab_matrix)

      relab_matrix = dplyr::arrange(.data = relab_matrix, .data[[time]])

      w = w[relab_matrix$index]
      relab_matrix$index = NULL
    }}else if(!is.null(group)){
      relab_matrix$index = 1:nrow(relab_matrix)

      relab_matrix = dplyr::arrange(.data = relab_matrix, .data[[group]])

      w = w[relab_matrix$index]
      relab_matrix$index = NULL
    }

  if(length(group) == 1){
    if(any(table(relab_matrix[[group]])<2)){
      ignore = names(which(table(relab_matrix[[group]])<2))
      warning("Only analyzing groups with at least two samples. Ignoring the following groups: ", paste(ignore, collapse = "  "))
      relab_matrix = relab_matrix[relab_matrix[[group]] %in% names(which(table(relab_matrix[[group]]) >= 2)),]
    }
  }


  relab_matrix_clean = relab_checker(relab_matrix,
                                     K = K,
                                     group = group,
                                     time = time)

  if(any(is.na(relab_matrix_clean$group))){
    stop("There is at least one NA in a grouping column.")
  }


  if(!is.null(time)){
    w = time_weights(times = relab_matrix_clean$time, group = relab_matrix_clean$group)
  }

  w_default = FALSE
  if((!is.null(group)) & (is.null(w))){
    w_default = TRUE
    w = rep(1/table(relab_matrix_clean$group), table(relab_matrix_clean$group))
  }

  if(is.null(w)){
    w = rep(1/nrow(relab_matrix), nrow(relab_matrix))
  }

  S = as.matrix(S)
  S = S_checker(S = S, K = K, relab_matrix = relab_matrix)

  if(!is.null(w) & length(w) != nrow(relab_matrix)){
    stop("Length of w must equal number of rows of relab_matrix.")
  }

  return(list(K = K,
              S = S,
              w = w,
              w_default = w_default,
              time = time,
              group = group,
              multiple_groups = multiple_groups,
              relab_matrix_clean = relab_matrix_clean,
              relab_grouping_vars = relab_grouping_vars))

}



# gini_simpson -----------------------------------------------------------------
#' Compute the Gini-Simpson index of a compositional vector
#'
#' This function computes the Gini-Simpson index, a statistical measure of variability known in population genetics as heterozygosity, of avector of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the variability of the vector. Values of 0 are achieved when the vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param q A vector with \code{K=length(q)} non-negative entries that sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=length(q)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(q))}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute unweighted Gini-Simpson index:
#' gini_simpson(q = c(0.4, 0.3, 0.3))
#'
#' # Compute Gini-Simpson index assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(3)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' gini_simpson(q = c(0.4, 0.3, 0.3), S = similarity_matrix)
#' @export
gini_simpson <- function(q, K = length(q), S = diag(K)){
  S = as.matrix(S)
  S = S_checker(S = S, K = K)
  q = unlist(q)[(length(q)-K+1):length(q)]
  q = as.numeric(q)

  if(round(sum(q), 4) != 1){stop("Vector does not sum to 1")}
  if(length(q) == 1){q = c(1,0)}

  1 - sum(q * c(S %*% q))
}

# gini_simpson_mean -----------------------------------------------------------------
#' Compute the mean Gini-Simpson index of the rows in a matrix of compositional vectors
#'
#' This function computes the mean Gini-Simpson index, a statistical measure of variability known in population genetics as heterozygosity, of a set of vectors of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the mean variability of the vectors. Values of 0 are achieved when each vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vectors are equal to (1/K, 1/K, ..., 1/K).
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # To compute the mean Gini-Simpson index of
#' # the following compositional vectors...
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#'
#' # we could compute the mean manually:
#' mean(sapply(list(q1, q2, q3, q4), gini_simpson))
#'
#' # Or we could use gini_simpson_mean:
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' gini_simpson_mean(relative_abundances)
#'
#'# Incoporating weights:
#'
#' # Compute mean Gini-Simpson index ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' gini_simpson_mean(relative_abundances, w = row_weights)
#'
#' # Compute mean Gini-Simpson index assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' gini_simpson_mean(relative_abundances, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' gini_simpson_mean(relative_abundances, w = row_weights, S = similarity_matrix)
#' @export
gini_simpson_mean <- function(relab_matrix,
                              K = NULL,
                              S = NULL,
                              w = NULL,
                              time = NULL,
                              group = NULL
                              # K = ncol(relab_matrix),
                              # w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                              # S = diag(ncol(relab_matrix))
){

  process_out = process_relab(relab_matrix = relab_matrix, K = K, S = S, w = w, time = time, group = group)

  K = process_out$K
  S = process_out$S
  w = process_out$w
  time = process_out$time
  group = process_out$group

  w_default = process_out$w_default
  multiple_groups = process_out$multiple_groups
  relab_matrix_clean = process_out$relab_matrix_clean
  relab_grouping_vars = process_out$relab_grouping_vars

  if(is.null(group)){
    gini_simpson_mean_fast(relab_matrix_clean$relab_matrix, K, w, S)
  } else{
    gs_list = c()
    for(subgroup in unique(relab_matrix_clean$group)){

      include = relab_matrix_clean$group == subgroup

      relab_sub = relab_matrix_clean$relab_matrix[include,]

      if(w_default){
        w_sub = w[names(w) == subgroup]
      } else{
        w_sub = w[include]
      }

      gs_list = c(gs_list,
                  gini_simpson_mean_fast(relab_sub, K, w_sub, S))
    }
    gs_df = data.frame(unique(relab_matrix_clean$group), gs_list)
    colnames(gs_df) = c(group, "gini_simpson_mean")
    # names(gs_list) = unique(relab_matrix_clean$group)

    if(multiple_groups){
      gs_df = dplyr::left_join(dplyr::distinct(relab_grouping_vars), gs_df)
    }
    return(gs_df)
  }

  # OLD VERSION BELOW
  # # K, S, w, time, and group are optional arguments
  #
  # # S = as.matrix(S)
  # I = nrow(relab_matrix)
  # # S = S_checker(S = S, K = K)
  #
  #
  # if(!missing(w) && length(w) != nrow(relab_matrix)){
  #   stop("Length of w must equal number of rows of relab_matrix.")
  # }
  #
  # relab_matrix = relab_matrix[,(ncol(relab_matrix)-K+1):ncol(relab_matrix)]
  #
  # # Average Gini-Simpson index of each of the I subpopulations
  # sum(w *sapply(1:I, function(i){ gini_simpson(q = unlist(relab_matrix[i,]), S = S) }))

}

# gini_simpson_pooled -----------------------------------------------------------------
#' Compute the pooled Gini-Simpson index of the rows in a matrix of compositional vectors
#'
#' This function computes the Gini-Simpson index of a "pooled" vector equal to \code{colMeans(relab_matrix)}. Values of 0 are achieved when this pooled vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when this pooled vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # To compute the pooled Gini-Simpson index of
#' # the following compositional vectors...
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#'
#' # we could compute the mean manually:
#' qPooled = (q1 + q2 + q3 + q4)/4
#' gini_simpson(qPooled)
#'
#' # Or we could use gini_simpson_pooled:
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' gini_simpson_pooled(relative_abundances)
#'
#'# Incoporating weights:
#'
#' # Compute pooled Gini-Simpson index ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' gini_simpson_pooled(relative_abundances, w = row_weights)
#'
#' # Compute pooled Gini-Simpson index assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' gini_simpson_pooled(relative_abundances, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' gini_simpson_pooled(relative_abundances, w = row_weights, S = similarity_matrix)
#' @export
gini_simpson_pooled <- function(relab_matrix,
                                K = NULL,
                                S = NULL,
                                w = NULL,
                                time = NULL,
                                group = NULL
                                # K = ncol(relab_matrix),
                                # w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                                # S = diag(ncol(relab_matrix))
){

  process_out = process_relab(relab_matrix = relab_matrix, K = K, S = S, w = w, time = time, group = group)

  K = process_out$K
  S = process_out$S
  w = process_out$w
  time = process_out$time
  group = process_out$group

  w_default = process_out$w_default
  multiple_groups = process_out$multiple_groups
  relab_matrix_clean = process_out$relab_matrix_clean
  relab_grouping_vars = process_out$relab_grouping_vars

  if(is.null(group)){
    gini_simpson_pooled_fast(relab_matrix_clean$relab_matrix, K, w, S)
  } else{
    gs_list = c()
    for(subgroup in unique(relab_matrix_clean$group)){

      include = relab_matrix_clean$group == subgroup

      relab_sub = relab_matrix_clean$relab_matrix[include,]

      if(w_default){
        w_sub = w[names(w) == subgroup]
      } else{
        w_sub = w[include]
      }

      gs_list = c(gs_list,
                  gini_simpson_pooled_fast(relab_sub, K, w_sub, S))
    }
    gs_df = data.frame(unique(relab_matrix_clean$group), gs_list)
    colnames(gs_df) = c(group, "gini_simpson_pooled")
    # names(gs_list) = unique(relab_matrix_clean$group)

    if(multiple_groups){
      gs_df = dplyr::left_join(dplyr::distinct(relab_grouping_vars), gs_df)
    }
    return(gs_df)
  }

  # OLD VERSION BELOW
  # # w and S are optional arguments
  #
  # I = nrow(relab_matrix)
  #
  # # S = as.matrix(S)
  # # S = S_checker(S = S, K = K)
  #
  # if(missing(S)){
  #   S = diag(K)
  # }
  #
  # if(missing(w)){
  #   w = rep(1, I)/I
  # }
  #
  #
  # if(!missing(w) && length(w) != nrow(relab_matrix)){
  #   stop("Length of w must equal number of rows of relab_matrix.")
  # }
  #
  # relab_matrix = relab_matrix[,(ncol(relab_matrix)-K+1):ncol(relab_matrix)]
  #
  # gini_simpson(q = colSums(sweep(x = relab_matrix, MARGIN = 1, w, `*`)), S = S)

}



# fava_norm -----------------------------------------------------------------
#' Compute the normalized Fst of a matrix of compositional vectors
#'
#' This function computes the normalized Fst given the number of rows and the mean abundance of the most abundant category.
#' We employ the normalization employed in the [FSTruct package](https://github.com/MaikeMorrison/FSTruct) by
#' Morrison, Alcala, and Rosenberg (2020) \doi{doi:10.1111/1755-0998.13647}.
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute the weighted fava of
#' # the following compositional vectors:
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' fava_norm(relative_abundances)
#' @export
fava_norm <- function(relab_matrix, K = ncol(relab_matrix)){
  relab_matrix = relab_matrix[,(ncol(relab_matrix)-K+1):ncol(relab_matrix)]
  I <- nrow(relab_matrix) # Here, I is the number of individuals (number of subpopulations)
  p <- colSums(relab_matrix) # yields vector of summed allele frequencies across populations
  sig1 <- max(p)
  J <- ceiling(1 / sig1)
  sig1.frac <- sig1 - floor(sig1)

  if (sig1 == I) {
    favaMax <- 0
    fava <- 0
    ratio <- 0
  } else {
    if (sig1 <= 1) {
      favaMax <- ((I - 1) * (1 - sig1 * (J - 1) * (2 - J * sig1))) /
        (I - (1 - sig1 * (J - 1) * (2 - J * sig1)))
    } else {
      favaMax <- (I * (I - 1) - sig1^2 + floor(sig1) - 2 * (I - 1) * sig1.frac + (2 * I - 1) * sig1.frac^2) / (I * (I - 1) - sig1^2 - floor(sig1) + 2 * sig1 - sig1.frac^2)
    }

    fava <-
      (sum(relab_matrix^2) / I - sum(colSums(relab_matrix / I)^2)) /
      (1 - sum(colSums(relab_matrix / I)^2))

    fava / favaMax
  }
}



# fava -----------------------------------------------------------------
#' Compute the Fst of a matrix of compositional vectors
#'
#' This function computes the population-genetic statistic Fst on any matrix with rows that sum to 1. Values of 0 are achieved when each row is a permutation of (1,0,..., 0) and at least two categories have non-zero abundance across all rows. The value equals 1 when each row is identical.
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @param normalized Optional; should normalized FAVA be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute the Fst of
#' # the following compositional vectors:
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' fava(relative_abundances)
#'
#'# Incoporating weights:
#'
#' # Compute fava ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' fava(relative_abundances, w = row_weights)
#'
#' # Compute fava assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' fava(relative_abundances, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' fava(relative_abundances, w = row_weights, S = similarity_matrix)
#' @export
fava <- function(relab_matrix,
                 K = NULL,
                 S = NULL,
                 w = NULL,
                 time = NULL,
                 group = NULL,
                 normalized = FALSE){

  if(normalized == TRUE && any(!sapply(list(time, w, S), is.null))){
    stop("FAVA can be either normalized or weighted, but not both. Please specify `normalized = TRUE` if you wish to compute normalized FAVA OR provide the weighting parameters w or time and/or S.")
  }

  if((!is.null(time)) && (!is.null(w))){
    stop("Please specify either time or w, but not both.")
  }

  process_out = process_relab(relab_matrix = relab_matrix, K = K, S = S, w = w, time = time, group = group)

  K = process_out$K
  S = process_out$S
  w = process_out$w
  time = process_out$time
  group = process_out$group

  w_default = process_out$w_default
  multiple_groups = process_out$multiple_groups
  relab_matrix_clean = process_out$relab_matrix_clean
  relab_grouping_vars = process_out$relab_grouping_vars



  if(is.null(group)){
    if(normalized){
      fava_norm(relab_matrix = relab_matrix_clean$relab_matrix)
    }else{
      (gini_simpson_pooled_fast(relab_matrix_clean$relab_matrix, K, w, S) -
         gini_simpson_mean_fast(relab_matrix_clean$relab_matrix, K, w, S))/
        gini_simpson_pooled_fast(relab_matrix_clean$relab_matrix, K, w, S)
    }} else{
      fava_list = c()
      for(subgroup in unique(relab_matrix_clean$group)){

        include = relab_matrix_clean$group == subgroup

        relab_sub = relab_matrix_clean$relab_matrix[include,]

        if(w_default){
          w_sub = w[names(w) == subgroup]
        } else{
          w_sub = w[include]
        }

        fava_list = c(fava_list,
                      ifelse(normalized,
                             fava_norm(relab_matrix = relab_sub),
                             (gini_simpson_pooled_fast(relab_sub, K, w_sub, S) -
                                gini_simpson_mean_fast(relab_sub, K, w_sub, S))/
                               gini_simpson_pooled_fast(relab_sub, K, w_sub, S)))
      }
      fava_df = data.frame(unique(relab_matrix_clean$group), fava_list)
      colnames(fava_df) = c(group, "FAVA")
      # names(fava_list) = unique(relab_matrix_clean$group)

      if(multiple_groups){
        fava_df = dplyr::left_join(dplyr::distinct(relab_grouping_vars), fava_df)
      }
      return(fava_df)
    }
}




# fast versions of functions to use in fava function:
gini_simpson_fast <- function(q, S = diag(length(q)), K = length(q)){

  1 - sum(q * c(S %*% q))
}

gini_simpson_mean_fast <- function(relab_matrix,
                                   K = ncol(relab_matrix),
                                   w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                                   S = diag(ncol(relab_matrix))){
  I = nrow(relab_matrix)

  # Average Gini-Simpson index of each of the I subpopulations
  sum(w *sapply(1:I, function(i){ gini_simpson_fast(q = unlist(relab_matrix[i,]), S = S) }))

}

gini_simpson_pooled_fast <- function(relab_matrix,
                                     K = ncol(relab_matrix),
                                     w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                                     S = diag(ncol(relab_matrix))){

  I = nrow(relab_matrix)

  gini_simpson_fast(q = colSums(sweep(x = relab_matrix, MARGIN = 1, w, `*`)), S = S)

}
