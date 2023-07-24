
# relab_checker ------------------------------------------------------------------
# An internal function to check if relab matrices are up to spec and fix any issues automatically.
# relab - a relative abundance matrix provided as input to another function
# K - Optional; the number of ancestral clusters
# rep - an optional parameter used if relab matrix is one of a list, in order to provide more useful warning messages
# group - Optional; the group that each row of the relative matrix abundance corresponds to
# time - Optional; the timepoint that each row of the relative abundance matrix corresponds to
relab_checker <- function(relab, K, rep = NULL, group = NULL, time = NULL) {
  original_relab = relab

  if(missing(K)){
      K = ncol(relab) - (!is.null(group)) - (!is.null(time))
  }

  # Check if relab matrix is within a list, and extract if needed
  if (is.list(relab) && !is.data.frame(relab) && !is.array(relab)) {
    relab <- relab[[1]]
  }
  # Check if relab matrix is in STRUCTURE/ADMIXTURE output form, or if it only contains columns of ancestry coefficients
  # If it is in STRUCTURE/ADMIXTURE output form, extract the ancestry coefficients--the last K columns.
  if (ncol(relab) > K) {
    relab <- relab[, (ncol(relab) - K + 1):ncol(relab)]
  }

  # convert relab matrix entries to numbers
  relab <- data.matrix(relab)

  # Name relab matrix columns q1, q2, ..., qK
  # colnames(relab) <- paste0("q",1:K)

  # Check if relab matrix has any missing values, and give warning if necessary
  if(any(is.na(relab))){
    # Identify location of missing entries
    na.pos <- sapply(which(is.na(relab)),
                     function(index) c(index %% nrow(relab), ceiling(index/nrow(relab)))) %>%
      t()
    # Format missing entries as a string
    na.pos.format <- list()
    for(row in 1:nrow(na.pos)){
      na.pos.format[row] <- paste0("(", na.pos[row,1], ", ", na.pos[row,2], ")")
    }
    na.pos.format.string <- as.character(na.pos.format) %>% paste(collapse = ", ")

    stop(paste0("There is at least one NA value in your relab matrix. The missing entries are found in the following positions: ",
                na.pos.format.string))
  }

  # check if matrix rows sum to 1, and give useful warnings if rounding is necessary
  sums <- rowSums(relab) %>% round(5)
  if (any(sums != 1)) {
    if (is.null(rep)) {
      warning("At least one relab matrix has rows which do not sum to exactly 1. Rounding the sum of each row to 1 by dividing all entries by the sum of the row.")
    } else {
      warning(paste0(
        "At least one of the rows of relab matrix number ", rep,
        " (restricted to the last K columns) does not sum to 1. Rounding the sum of each row to 1 by dividing all entries by the sum of the row."
      ))
    }
    # Normalize each row of the matrix by dividing by the rowsums
    relab <- relab / sums
  }
  if(!is.null(group) | !is.null(time)){
    if(is.null(group)){
      group_return = NULL
    } else{
      group_return = unlist(original_relab[,group], use.names = FALSE)
    }

    if(is.null(time)){
      time_return = NULL
    } else{
      time_return = unlist(original_relab[,time])
    }


    return(list("relab_matrix" = relab,
                "group" = group_return,
                "time" = time_return))
  }
  return(relab)
}

# relab_checker(relab = xue_microbiome_sample, group = "subject", time = "timepoint")
