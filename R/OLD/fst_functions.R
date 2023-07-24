# THIS FILE CONTAINS:
# het - compute the heterozygosity of a compositional vector
# hetMean
# hetPooled
# fava
# fava_norm
# time_weights

# S_checker ------------------------------------------------------------------
# An internal function to check if similarity matrices are up to spec and fix any issues automatically.
# S - a similarity matrix
# K - the number of taxa
S_checker <- function(S, K) {

  if(!isSymmetric(S)){
    warning("S is not symmetric.")
  }
  if(any(diag(S)!=1)){
    stop("All diagonal elements of S must be equal to 1.")
  }
  if(nrow(S) != K | ncol(S) != K){
    stop("S must have K rows and K columns.")
  }
  if(any(S>1)){
    stop("All elements of S must be less than or equal to 1.")
  }
  if(any(S<0)){
    stop("All elements of S must be positive.")
  }
}

# het -----------------------------------------------------------------
#' Compute the heterozygosity of a compositional vector
#'
#' This function computes the heterozygosity, a statistical measure of variability also known as the Gini-Simpson index, of vector of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the variability of the vector. Values of 0 are achieved when the vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param q A vector with \code{K=length(q)} non-negative entries that sum to 1.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(q))}.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=length(q)}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute unweighted heterozygosity:
#' het(q = c(0.4, 0.3, 0.3))
#'
#' # Compute heterozygosity assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(3)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' het(q = c(0.4, 0.3, 0.3), S = similarity_matrix)
#' @export
het <- function(q, S = diag(length(q)), K = length(q)){
  S = as.matrix(S)
  S_checker(S = S, K = K)
  q = unlist(q)[(length(q)-K+1):length(q)]

  if(round(sum(q), 4) != 1){stop("Vector does not sum to 1")}
  if(length(q) == 1){q = c(1,0)}

  1 - sum(q * c(S %*% q))
}

# hetMean -----------------------------------------------------------------
#' Compute the mean heterozygosity of the rows in a matrix of compositional vectors
#'
#' This function computes the mean heterozygosity, a statistical measure of variability also known as the Gini-Simpson index, of a set of vectors of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the mean variability of the vectors. Values of 0 are achieved when each vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vectors are equal to (1/K, 1/K, ..., 1/K).
#'
#' @param Q A matrix with \code{I=nrow(Q)} rows, each containing \code{K=ncol(Q)} non-negative entries that sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(Q)}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean variability across rows. The default value is \code{w = rep(1/nrow(Q), nrow(Q))}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(Q))}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # To compute the mean heterozygosity of
#' # the following compositional vectors...
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#'
#' # we could compute the mean manually:
#' mean(sapply(list(q1, q2, q3, q4), het))
#'
#' # Or we could use hetMean:
#' Q_matrix = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' hetMean(Q_matrix)
#'
#'# Incoporating weights:
#'
#' # Compute mean heterozygosity ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' hetMean(Q_matrix, w = row_weights)
#'
#' # Compute mean heterozygosity assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' hetMean(Q_matrix, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' hetMean(Q_matrix, w = row_weights, S = similarity_matrix)
#' @export
hetMean <- function(Q,
                    K = ncol(Q),
                    w = rep(1/nrow(Q), nrow(Q)),
                    S = diag(ncol(Q))){
  # K, w, and S are optional arguments

  # S = as.matrix(S)
  I = nrow(Q)
  # S_checker(S = S, K = K)


  if(!missing(w) && length(w) != nrow(Q)){
    stop("Length of w must equal number of rows of Q.")
  }

  Q = Q[,(ncol(Q)-K+1):ncol(Q)]

  # Average heterozygosity of each of the I subpopulations
  sum(w *sapply(1:I, function(i){ het(q = unlist(Q[i,]), S = S) }))

}

# hetPooled -----------------------------------------------------------------
#' Compute the pooled heterozygosity of the rows in a matrix of compositional vectors
#'
#' This function computes the heterozygosity of a "pooled" vector equal to \code{colMeans(Q)}. Values of 0 are achieved when this pooled vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when this pooled vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param Q A matrix with \code{I=nrow(Q)} rows, each containing \code{K=ncol(Q)} non-negative entries that sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(Q)}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(Q), nrow(Q))}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(Q))}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # To compute the pooled heterozygosity of
#' # the following compositional vectors...
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#'
#' # we could compute the mean manually:
#' qPooled = (q1 + q2 + q3 + q4)/4
#' het(qPooled)
#'
#' # Or we could use hetPooled:
#' Q_matrix = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' hetPooled(Q_matrix)
#'
#'# Incoporating weights:
#'
#' # Compute pooled heterozygosity ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' hetPooled(Q_matrix, w = row_weights)
#'
#' # Compute pooled heterozygosity assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' hetPooled(Q_matrix, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' hetPooled(Q_matrix, w = row_weights, S = similarity_matrix)
#' @export
hetPooled <- function(Q,
                      K = ncol(Q),
                      w = rep(1/nrow(Q), nrow(Q)),
                      S = diag(ncol(Q))){
  # w and S are optional arguments

  I = nrow(Q)

  # S = as.matrix(S)
  # S_checker(S = S, K = K)

  if(missing(S)){
    S = diag(K)
  }

  if(missing(w)){
    w = rep(1, I)/I
  }


  if(!missing(w) && length(w) != nrow(Q)){
    stop("Length of w must equal number of rows of Q.")
  }

  Q = Q[,(ncol(Q)-K+1):ncol(Q)]

  het(q = colSums(sweep(x = Q, MARGIN = 1, w, `*`)), S = S)

}

# fava -----------------------------------------------------------------
#' Compute the Fst of a matrix of compositional vectors
#'
#' This function computes the population genetic statistic Fst on any matrix with rows that sum to 1. Values of 0 are achieved when each row is a permutation of (1,0,..., 0) and at least two categories have non-zero abundance across all rows. The value equals 1 when each row is identical.
#'
#' @param Q A matrix with \code{I=nrow(Q)} rows, each containing \code{K=ncol(Q)} non-negative entries that sum to 1.
#' If \code{Q} contains any metadata, it must be on the left-hand side of the matrix and the number of entries
#' that sum to 1 (\code{K}) must be specified.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(Q)}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(Q), nrow(Q))}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(Q))}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute the Fst of
#' # the following compositional vectors:
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#
#' Q_matrix = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' fava(Q_matrix)
#'
#'# Incoporating weights:
#'
#' # Compute fava ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' fava(Q_matrix, w = row_weights)
#'
#' # Compute fava assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' fava(Q_matrix, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' fava(Q_matrix, w = row_weights, S = similarity_matrix)
#' @export
fava <- function(Q, w = rep(1/nrow(Q), nrow(Q)), S = diag(ncol(Q)), K = ncol(Q)){
  # S = as.matrix(S)
  Q = Q[,(ncol(Q)-K+1):ncol(Q)]
  (hetPooled(Q, K, w, S) - hetMean(Q, K, w, S))/hetPooled(Q, K, w, S)
}

# fava_norm -----------------------------------------------------------------
#' Compute the normalized Fst of a matrix of compositional vectors
#'
#' This function computes the normalized Fst given the number of rows and the mean abundance of the most abundant category.
#' We employ the normalization employed in the [FSTruct package](https://github.com/MaikeMorrison/FSTruct) by
#' [Morrison, Alcala, and Rosenberg (2020)](https://doi.org/10.1111/1755-0998.13647).
#'
#' @param Q A matrix with \code{I=nrow(Q)} rows, each containing \code{K=ncol(Q)} non-negative entries that sum to 1.
#' If \code{Q} contains any metadata, it must be on the left-hand side of the matrix and the number of entries
#' that sum to 1 (\code{K}) must be specified.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(Q)}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute the weighted fava of
#' # the following compositional vectors:
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#
#' Q_matrix = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' fava_norm(Q_matrix)
#' @export
fava_norm <- function(Q, K = ncol(Q)){
  Q = Q[,(ncol(Q)-K+1):ncol(Q)]
  I <- nrow(Q) # Here, I is the number of individuals (number of subpopulations)
  p <- colSums(Q) # yields vector of summed allele frequencies across populations
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
      (sum(Q^2) / I - sum(colSums(Q / I)^2)) /
      (1 - sum(colSums(Q / I)^2))

    fava / favaMax
  }
}


# time_weights -----------------------------------------------------------------
#' Compute a normalized weighting vector based on a vector of sampling times.
#'
#' This function takes a vector of sampling times, \eqn{t = (t_1, t_2, \ldots, t_I)}{latex}
#' and computes a normalized vector which can be used to weight each sample based on
#' the time between the subsequent and the preceding samples. The weighting vector \eqn{w}
#' is defined such that each entry, \eqn{w_i = d_i / 2T}, where \eqn{T=t_I - t_1} and
#' \eqn{d_i = t_{i+1} - t_{i-1}} for \eqn{i} not equal to 1 or I. \eqn{d_1 = t_2-t_1} and \eqn{d_I = t_I-t_{I-1}}.
#'
#' @param times A numeric vector of sampling times. Each entry must be non-negative and
#' greater than the previous entry.
#' @param group Optional; a character vector specifying the group identity of each
#' sampling time. Use if there are samples from multiple replicates or subjects
#' in one data set.
#' @examples

#' time_vector = c(1, 8, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
#'                 32, 33, 34, 35, 36, 37, 38, 39, 44, 50, 57, 64)
#'
#' time_weights(times = time_vector)
#' @export
time_weights <- function(times, group){

  if(missing(group)){
    I = length(times)
    T = times[I] - times[1]

    if(I<2){
      stop("times must have length greater than 1.")
    }
    if(any(sapply(2:I, function(i) times[i] - times[i-1]) <= 0)){
      stop("times must be increasing. Each entry must be greater than the previous entry.")
    }

    di <- c(times[2] - times[1])
    for(i in 2:(I-1)){
      di[i] = times[i+1] - times[i-1]
    }
    di[I] = times[I] - times[I-1]
    return(di/(2*T))
  }else{

    wi = c()
    for(name in unique(group)){
      time_name = times[which(group == name)]

      I = length(time_name)
      T = time_name[I] - time_name[1]

      if(I<2){
        stop("Within each group, times must have length greater than 1.")
      }
      if(any(sapply(2:I, function(i) time_name[i] - time_name[i-1]) <= 0)){
        stop("Within each group, times must be increasing. Each entry must be greater than the previous entry.")
      }

      wi <- c(wi, (time_name[2] - time_name[1])/(2*T))
      for(i in 2:(I-1)){
        wi = c(wi, (time_name[i+1] - time_name[i-1])/(2*T))
      }
      wi = c(wi, (time_name[I] - time_name[I-1])/(2*T))
    }
    return(wi)
  }

}






# # 0: unweighted functions ---------------------------------------------------------------
#
# # compute the heterozygosity of a vector that sums to 1
# het <- function(q){
#   if(round(sum(q), 4) != 1){stop("Vector does not sum to 1")}
#   if(length(q) == 1){q = c(1,0)}
#   1 - sum(q^2)
# }
#
# # Subpopulation heterozygosity
# Hs <- function(Q){
#   I = nrow(Q)
#
#   # Average heterozygosity of each of the I subpopulations
#   sum(1/I *sapply(1:I, function(i){ het(q = Q[i,]) }))
# }
#
# # Total heterozygosity
# Ht <- function(Q){
#   I = nrow(Q)
#
#   qmean <- colSums(Q)/I
#
#   het(qmean)
# }
#
#
# fava <- function(Q) (Ht(Q) - Hs(Q))/Ht(Q)
#
#
#
# # 1: weighted functions  -----------------------------------------------------------
