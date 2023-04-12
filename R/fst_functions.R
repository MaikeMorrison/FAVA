# THIS FILE CONTAINS:
# het - compute the heterozygosity of a compositional vector


# het -----------------------------------------------------------------
#' Compute the heterozygosity of a compositional vector
#'
#' This function computes the heterozygosity, a statistical measure of variability also known as the Gini-Simpson index, of vector of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the variability of the vector. Values of 0 are achieved when the vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param q A vector with \code{K=length(q)} non-negative entries that sum to 1.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(q))}.
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
het <- function(q, S = diag(length(q))){
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
hetMean <- function(Q, w = rep(1/nrow(Q), nrow(Q)), S = diag(ncol(Q))){
  # w and S are optional arguments

  I = nrow(Q)
  K = ncol(Q)


  if(any(diag(S)!=1)){
    stop("The diagonal elements of S must equal 1.")
  }
  if(any(S>1)){
    stop("All elements of S must be less than or equal to 1.")
  }
  if(any(S<0)){
    stop("All elements of S must be positive.")
  }
  if(!missing(w) && length(w) != nrow(Q)){
    stop("Length of w must equal number of rows of Q.")
  }

  # Average heterozygosity of each of the I subpopulations
  sum(w *sapply(1:I, function(i){ het(q = unlist(Q[i,]), S = S) }))

}

# hetPooled -----------------------------------------------------------------
#' Compute the pooled heterozygosity of the rows in a matrix of compositional vectors
#'
#' This function computes the heterozygosity of a "pooled" vector equal to \code{colMeans(Q)}. Values of 0 are achieved when this pooled vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when this pooled vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param Q A matrix with \code{I=nrow(Q)} rows, each containing \code{K=ncol(Q)} non-negative entries that sum to 1.
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
hetPooled <- function(Q, w, S){
  # w and S are optional arguments

  I = nrow(Q)
  K = ncol(Q)

  if(missing(S)){
    S = diag(K)
  }

  if(missing(w)){
    w = rep(1, I)/I
  }

  if(any(diag(S)!=1)){
    stop("The diagonal elements of S must equal 1.")
  }
  if(any(S>1)){
    stop("All elements of S must be less than or equal to 1.")
  }
  if(any(S<0)){
    stop("All elements of S must be positive.")
  }
  if(!missing(w) && length(w) != nrow(Q)){
    stop("Length of w must equal number of rows of Q.")
  }

  het(q = colSums(sweep(x = Q, MARGIN = 1, w, `*`)), S = S)

}

# fst -----------------------------------------------------------------
#' Compute the Fst of a matrix of compositional vectors
#'
#' This function computes the population genetic statistic Fst on any matrix with rows that sum to 1. Values of 0 are achieved when each row is a permutation of (1,0,..., 0) and at least two categories have non-zero abundance across all rows. The value equals 1 when each row is identical.
#'
#' @param Q A matrix with \code{I=nrow(Q)} rows, each containing \code{K=ncol(Q)} non-negative entries that sum to 1.
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
#' fst(Q_matrix)
#'
#'# Incoporating weights:
#'
#' # Compute Fst ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' fst(Q_matrix, w = row_weights)
#'
#' # Compute Fst assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' fst(Q_matrix, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' fst(Q_matrix, w = row_weights, S = similarity_matrix)
#' @export
fst <- function(Q, w = rep(1/nrow(Q), nrow(Q)), S = diag(ncol(Q))){
  (hetPooled(Q, w, S) - hetMean(Q, w, S))/hetPooled(Q, w, S)
}

# fst_norm -----------------------------------------------------------------
#' Compute the normalized Fst of a matrix of compositional vectors
#'
#' This function computes the normalized Fst given the number of rows and the mean abundance of the most abundant category. We employ the normalization employed in the [FSTruct package](https://github.com/MaikeMorrison/FSTruct) by [Morrison, Alcala, and Rosenberg (2020)](https://doi.org/10.1111/1755-0998.13647).
#'
#' @param Q A matrix with \code{I=nrow(Q)} rows, each containing \code{K=ncol(Q)} non-negative entries that sum to 1.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute the weighted Fst of
#' # the following compositional vectors:
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#
#' Q_matrix = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' fst_norm(Q_matrix)
#' @export
fst_norm <- function(Q){
  I <- nrow(Q) # Here, I is the number of individuals (number of subpopulations)
  p <- colSums(Q) # yields vector of summed allele frequencies across populations
  sig1 <- max(p)
  J <- ceiling(1 / sig1)
  sig1.frac <- sig1 - floor(sig1)

  if (sig1 == I) {
    FstMax <- 0
    Fst <- 0
    ratio <- 0
  } else {
    if (sig1 <= 1) {
      FstMax <- ((I - 1) * (1 - sig1 * (J - 1) * (2 - J * sig1))) /
        (I - (1 - sig1 * (J - 1) * (2 - J * sig1)))
    } else {
      FstMax <- (I * (I - 1) - sig1^2 + floor(sig1) - 2 * (I - 1) * sig1.frac + (2 * I - 1) * sig1.frac^2) / (I * (I - 1) - sig1^2 - floor(sig1) + 2 * sig1 - sig1.frac^2)
    }

    Fst <-
      (sum(Q^2) / I - sum(colSums(Q / I)^2)) /
      (1 - sum(colSums(Q / I)^2))

    Fst / FstMax
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
# fst <- function(Q) (Ht(Q) - Hs(Q))/Ht(Q)
#
#
#
# # 1: weighted functions  -----------------------------------------------------------
