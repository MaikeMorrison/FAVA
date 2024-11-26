# relab_phyloseq ---------------------------------------------------------------
#' Generate a relative abundance matrix with sample metadata and OTU abundances from a  phyloseq object.
#'
#' The R package phyloseq streamlines the storage and analysis of microbiome sequence data. This function takes a phyloseq object and extracts the OTU table and the sample metadata and combines them into one relative abundance matrix with rows corresponding to samples, metadata on the left-hand side, and OTU relative abundances on the right-hand side.
#'
#' @param phyloseq_object A phyloseq object containing both an OTU table (`otu_table`) and sample metadata (`sample_data`).
#' @returns A data frame with rows representing samples and columns representing sample data categories or OTU relative abundances.
#' OTU abundances are automatically normalized so that they sum to 1 for each sample, though a warning will be provided if a
#' renormalization was necessary.
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data(GlobalPatterns, package = "phyloseq")
#'
#' # Make a small phyloseq object for demonstration
#' phyloseq_subset = phyloseq::subset_taxa(phyloseq::subset_samples(GlobalPatterns,
#'                                                                  X.SampleID %in%
#'                                                                  c("CL3", "CC1")),
#'                                         Order == "Cenarchaeales")
#'   otu_table = relab_phyloseq(phyloseq_subset)
#'   otu_table[, 1:10]
#' }
#' @export
relab_phyloseq <- function(phyloseq_object){
  if(is.null(phyloseq_object@sam_data)){
    warning("phyloseq_object does not have sample_data")
    if(phyloseq::taxa_are_rows(phyloseq_object)){
      return(t(phyloseq::otu_table(phyloseq_object)))
    }else{
      return(phyloseq::otu_table(phyloseq_object))
    }
  }
  if(is.null(phyloseq_object@otu_table)){
    stop("phyloseq_object does not include an otu_table")
  }

  if(phyloseq::taxa_are_rows(phyloseq_object)){
    output = cbind(phyloseq::sample_data(phyloseq_object),
          t(phyloseq::otu_table(phyloseq_object)))
  }else{
    output = cbind(phyloseq::sample_data(phyloseq_object),
          (phyloseq::otu_table(phyloseq_object)))
  }

  meta_cols = ncol(phyloseq::sample_data(phyloseq_object))
  data_cols = (meta_cols+1):ncol(output)

  relab_sums = rowSums(output[,data_cols])

  if(any(relab_sums == 0)){
    warning(paste0("The following samples summed to 0 and were excluded: ", paste0(names(which(relab_sums==0)), collapse = ", ")))
    output = output[-which(relab_sums==0),]
    relab_sums = rowSums(output[,data_cols])
  }

  if(!all(round(relab_sums, 7)==1)){
    output[,data_cols] = output[,data_cols]/relab_sums
    warning("Some of the sample abundances do not sum to exactly 1. Rounding the sum of each sample to 1 by dividing all entries by the sum of the sample.")
  }

  return(output)

}



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






# FUNCTIONS FOR WEIGHTINGS

# S_checker ------------------------------------------------------------------
# An internal function to check if similarity matrices are up to spec and fix any issues automatically.
# S - a similarity matrix
# K - the number of taxa
S_checker <- function(S, K, relab_matrix = NULL) {

  if(!is.null(relab_matrix) & !is.null(colnames(relab_matrix)) & !is.null(colnames(S))){

    taxa_names = colnames(relab_matrix)[(ncol(relab_matrix)-K + 1):ncol(relab_matrix)]

    if(length(taxa_names) != ncol(S) & !all(taxa_names %in% colnames(S))){
      stop("Not all of the taxa in relab_matrix are included in S. Your similarity matrix must describe every taxon included in relab_matrix.")
    }

    if(length(taxa_names) < ncol(S)){
      S = S[taxa_names, taxa_names]
      warning("S describes more taxa than are included in relab_matrix.\nI am subsetting the rows and columns of S so that they match the names of the K categories in relab_matrix.")
    }

    else if(any(taxa_names!=colnames(S)) || any(taxa_names!=rownames(S))){
      S = S[taxa_names, taxa_names]
      warning("The row and/or column names of the similarity matrix S do not match the names of the K categories in relab_matrix.\nI am re-ordering the rows and columns of S so that they match the names of the K categories in relab_matrix.")
    }

  }


  if(!isSymmetric(S)){
    warning("S is not symmetric.")
  }
  if(any(round(diag(S),8)!=1)){
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

  return(S)
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
#' @param times A numeric vector of sampling times. Each entry must be
#' greater than the previous entry.
#' @param group Optional; a character vector specifying the group identity of each
#' sampling time. Use if there are samples from multiple replicates or subjects
#' in one dataset.
#' @returns A numeric vector. Each entry provides a weight for each entry in the
#' provided `times` vector. If `group` is not specified, the vector sums to 1. If
#' `group` is specified, the vector sums to the number of distinct groups.
#' @examples

#' time_vector = c(1, 8, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
#'                 32, 33, 34, 35, 36, 37, 38, 39, 44, 50, 57, 64)
#'
#' time_weights(times = time_vector)
#' @export
time_weights <- function(times, group = NULL){

  if(!is.null(group)){
    if(length(group) != length(times)){
      stop(paste0("The times and group vectors must have the same length. times has length ",
           length(times), ". group has length ", length(group), "."))
    }
  }

  # Create an ID for each sample to ensure the returned sample weights match the
  # order that the samples were given in.
  ID = paste0(group, times)


  if(is.null(group)){
    I = length(times)
    T = times[I] - times[1]

    if(I<2){
      stop("times must have length greater than 1.")
    }
    if(any(sapply(2:I, function(i) times[i] - times[i-1]) <= 0)){
      stop("times must be increasing. Each entry must be greater than the previous entry.")
    }

    di <- c(times[2] - times[1])
    if(I>2){
      for(i in 2:(I-1)){
        di[i] = times[i+1] - times[i-1]
      }
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

      wi_addition <- c(name1 = (time_name[2] - time_name[1])/(2*T))
      if(I>2){
        for(i in 2:(I-1)){
          wi_addition = c(wi_addition, (time_name[i+1] - time_name[i-1])/(2*T))
        }
      }
      wi_addition = c(wi_addition, (time_name[I] - time_name[I-1])/(2*T))

      names(wi_addition) = paste0(name, time_name)

      wi = c(wi, wi_addition)
    }
    return(wi[ID])
  }

}








# Repeat rows based on time series data, to be used
# when generating boostrap replicates of matrices

# lcm_pair = function(x, y){
#   if(x>y)
#   {
#     greater=x
#   }
#   else
#   {
#     greater=y
#   }
#   while(TRUE)
#   {
#     if((greater%%x==0)&&(greater%%y==0))
#     {
#       lcm=greater
#       break
#     }
#     greater=greater+1
#   }
#   return(lcm)
# }
#
# lcm = function(integers){Reduce(lcm_pair, integers)}

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

  if(arrange == FALSE){
    relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)
    relab_matrix_clean = relab_checker_out$relab_matrix
    relab_matrix_clean = relab_matrix_clean[,colSums(relab_matrix_clean) > 0]
  }else{

  relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

  K = ncol(relab_checker_out$relab_matrix)

  relab_matrix_clean = relab_checker_out$relab_matrix
  relab_matrix_clean = relab_matrix_clean[,colSums(relab_matrix_clean) > 0]

  if(is.null(colnames(relab_matrix_clean))){
    colnames(relab_matrix_clean) = paste0("cat_", 1:K)
  }


  if(any(colnames(relab_matrix_clean) %in% as.character(1:1000))){
    warning("relab_matrix has at least one numeric column name, which can cause errors when re-arranging the columns. The columns are being automatically renamed in order to avoid errors. If you wish to avoid this renaming, please rename the columns of relab_matrix.")
    colnames(relab_matrix_clean) = paste0("cat_", 1:K)
  }


    clustermeans <- colMeans(relab_matrix_clean) %>% sort() %>% rev
    ordernames <- c(names(clustermeans))

    if(sum(clustermeans != 0) != K){
      warning(paste0("Only plotting the ", sum(clustermeans != 0),
                     " categories with non-zero abundances. If you are manually changing the fill or color of the plot, you can provide ",
                     sum(clustermeans != 0), " colors, instead of ", K, "."))
    }

    if(arrange %in% c(TRUE, "both")){

      relab_matrix_clean <- data.frame(relab_matrix_clean) %>%
        dplyr::arrange(dplyr::across({{ ordernames }})) %>%
        dplyr::select(c(names(which(clustermeans != 0))))

    }else if(arrange == "vertical"){
      relab_matrix_clean <- data.frame(relab_matrix_clean) %>%
        # dplyr::arrange(dplyr::across({{ ordernames }})) %>%
        dplyr::select(c(names(which(clustermeans != 0))))


    }else if(arrange == "horizontal"){
      relab_matrix_clean <- data.frame(relab_matrix_clean) %>%
        dplyr::arrange(dplyr::across({{ ordernames }}))# %>%
      # dplyr::select(c("group", names(which(clustermeans != 0))))
    }else{
      stop("The options for arrange are FALSE, TRUE, vertical, horizontal, or both. TRUE and both are interchangeable.")
    }

  }

  if(!is.null(group)){
    relab_matrix_clean = relab_matrix_clean %>% data.frame %>% dplyr::mutate(group = relab_checker_out$group, .before = 1)
  }

  if(!is.null(time)){
    relab_matrix_clean = relab_matrix_clean %>% data.frame %>% dplyr::mutate(time = relab_checker_out$time, .before = 1)
  }

  return(relab_matrix_clean %>% data.frame)



}

