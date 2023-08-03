
# relab_checker ------------------------------------------------------------------
# An internal function to check if relab matrices are up to spec and fix any issues automatically.
# relab - a relative abundance matrix provided as input to another function
# K - Optional; the number of ancestral clusters
# rep - an optional parameter used if relab matrix is one of a list, in order to provide more useful warning messages
# group - Optional; the group that each row of the relative matrix abundance corresponds to
# time - Optional; the timepoint that each row of the relative abundance matrix corresponds to
relab_checker <- function(relab, K = NULL, rep = NULL, group = NULL, time = NULL) {
  original_relab = relab

  if(is.null(K)){
      K = ncol(relab) - (!is.null(group)) - (!is.null(time))
  }

  # Check if relab matrix is within a list, and extract if needed
  if (is.list(relab) && !is.data.frame(relab) && !is.array(relab)) {
    relab <- relab[[1]]
  }
  # If the matrix contains meta information, extract the last K columns.
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

# relab_checker(relab = xue_microbiome_sample, group = "subject", time = "timepoint")








# Repeat rows based on time series data, to be used
# when generating boostrap replicates of matrices
relab_sample_weighter = function(relab, K = NULL, time = NULL, w = NULL, group = NULL){

  relab_matrix_clean = relab_checker(relab = relab, K = K, time = time, group = group)

  if(is.null(group)){

    if(is.null(w) & (!is.null(time))){
      t = relab_matrix_clean$time
      w = time_weights(t)
      T = max(t) - min(t)

      return(relab[rep(1:length(w), round(w*T*2)),])
    }else  if(is.null(time) & !is.null(w)){
      return(relab[rep(1:length(w), round(w*800)),])
    } else{
      # warning("Please provide either time or w to relab_sample_weighter function.")
      return(relab)
    }

  }else{
    groups = unique(relab_matrix_clean$group)
    df_list = vector(mode='list', length = length(groups))
    i = 1
    for(g in groups){
      group_sub = relab_matrix_clean$group == g

      if(is.null(w) & (!is.null(time))){

        t_sub = (relab_matrix_clean$time)[group_sub]
        w_sub = time_weights(t_sub)
        T = max(t_sub) - min(t_sub)

        df_list[[i]] = relab[group_sub,][rep(1:length(w_sub), round(w_sub*T*2)),]

      }else  if(is.null(time) & !is.null(w)){
        df_list[[i]] = relab[group_sub,][rep(1:length(w), round(w*800)),]
      }else{
        # warning("Please provide either time or w to relab_sample_weighter function.")
        return(relab)
      }
      i = i + 1
    }

    return(do.call(rbind, df_list))

  }

}


# helper function: arrange

arrange_categories <- function(relab_matrix, arrange, K = NULL, group = NULL, time = NULL){

  relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

  K = ncol(relab_checker_out$relab_matrix)

  relab_matrix = relab_checker_out$relab_matrix

  if(is.null(colnames(relab_matrix))){
    colnames(relab_matrix) = paste0("cat_", 1:K)
  }


  if(any(colnames(relab_matrix) %in% as.character(1:1000))){
    warning("relab_matrix has at least one numeric column name, which can cause errors when re-arranging the columns. The columns are being automatically renamed in order to avoid errors. If you wish to avoid this renaming, please rename the columns of relab_matrix.")
    colnames(relab_matrix) = paste0("cat_", 1:K)
  }

  if (arrange !=FALSE){
    clustermeans <- colMeans(relab_matrix) %>% sort() %>% rev
    ordernames <- c(names(clustermeans))

    if(sum(clustermeans != 0) != K){
      warning(paste0("Only plotting the ", sum(clustermeans != 0),
                     " categories with non-zero abundances. If you are manually changing the fill or color of the plot, you can provide ",
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

