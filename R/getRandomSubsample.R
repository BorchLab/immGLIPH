#' Draw a stratified random subsample from the reference repertoire
#'
#' Draws a random subset of reference motif regions with the same size as the
#' sample set. When \code{cdr3_len_stratify} and/or \code{vgene_stratify} are
#' enabled, the function preserves the CDR3 length and/or V-gene distribution
#' of the sample in the subsample. This is used internally by the repeated
#' random sampling (RRS) local-similarity method in \code{\link{runGLIPH}}.
#'
#' @param cdr3_len_stratify Whether to preserve the CDR3 length distribution.
#'   **Default:** `FALSE`
#' @param vgene_stratify Whether to preserve the V-gene distribution.
#'   **Default:** `FALSE`
#' @param refseqs_motif_region Character vector of reference motif regions.
#' @param motif_region Character vector of sample motif regions.
#' @param motif_lengths_list Named list mapping CDR3 lengths to their
#'   frequency in \code{motif_region}. Required when
#'   \code{cdr3_len_stratify = TRUE}.
#' @param ref_motif_lengths_id_list Named list mapping CDR3 lengths to
#'   indices in \code{refseqs_motif_region}. Required when
#'   \code{cdr3_len_stratify = TRUE}.
#' @param motif_region_vgenes_list Named list mapping V-genes to their
#'   frequency in \code{motif_region}. Required when
#'   \code{vgene_stratify = TRUE}.
#' @param ref_motif_vgenes_id_list Named list mapping V-genes to indices in
#'   \code{refseqs_motif_region}. Required when \code{vgene_stratify = TRUE}.
#' @param lengths_vgenes_list Nested list mapping CDR3 length x V-gene
#'   combinations to their frequency in the sample. Required when both
#'   stratification flags are \code{TRUE}.
#' @param ref_lengths_vgenes_list Nested list mapping CDR3 length x V-gene
#'   combinations to indices in \code{refseqs_motif_region}. Required when
#'   both stratification flags are \code{TRUE}.
#'
#' @return A character vector of length \code{length(motif_region)} drawn from
#'   \code{refseqs_motif_region}.
#'
#' @export

getRandomSubsample <- function(cdr3_len_stratify = FALSE,
                               vgene_stratify = FALSE,
                               refseqs_motif_region,
                               motif_region,
                               motif_lengths_list,
                               ref_motif_lengths_id_list,
                               motif_region_vgenes_list,
                               ref_motif_vgenes_id_list,
                               ref_lengths_vgenes_list,
                               lengths_vgenes_list){
  random_subsample <- c()

  ### Return an unbiased reference subsample
  if(cdr3_len_stratify == FALSE && vgene_stratify == FALSE){
    random_subsample <- refseqs_motif_region[sample(x = seq_along(refseqs_motif_region),size = length(motif_region),replace = FALSE)]
  }

  ### Return a reference subsample with biased CDR3b length distribution according to the sample
  if(cdr3_len_stratify == TRUE && vgene_stratify == FALSE){

    ## if there are more seqs with specific cdr3 length in sample than in reference database, surplus number of seqs will be taken randomly from remaining seqs with different cdr3 length
    random_length <- 0
    subsample_parts <- vector("list", length(motif_lengths_list))
    idx <- 0L
    for(cdr3_length in names(motif_lengths_list)){
      idx <- idx + 1L
      if(length(ref_motif_lengths_id_list[[cdr3_length]]) < motif_lengths_list[[cdr3_length]]){
        random_length <- random_length + motif_lengths_list[[cdr3_length]] - length(ref_motif_lengths_id_list[[cdr3_length]])
        subsample_parts[[idx]] <- ref_motif_lengths_id_list[[cdr3_length]]
      }else {
        subsample_parts[[idx]] <- sample(x = ref_motif_lengths_id_list[[cdr3_length]],size = motif_lengths_list[[cdr3_length]],replace = FALSE)
      }
    }
    random_subsample <- unlist(subsample_parts)
    if(random_length > 0){
      random_subsample <- c(random_subsample, sample(x = (seq_along(refseqs_motif_region))[-random_subsample],size = random_length,replace = FALSE))
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }

  ### Return a reference subsample with biased V gene distribution according to the sample
  if(vgene_stratify == TRUE && cdr3_len_stratify == FALSE){

    ## if there are more seqs with specific v gene in sample than in reference database, surplus number of seqs will be taken randomly from remaining seqs with different vgenes
    random_length <- 0
    subsample_parts <- vector("list", length(motif_region_vgenes_list))
    idx <- 0L
    for(act_vgene in names(motif_region_vgenes_list)){
      idx <- idx + 1L
      if(length(ref_motif_vgenes_id_list[[act_vgene]]) < motif_region_vgenes_list[[act_vgene]]){
        random_length <- random_length + motif_region_vgenes_list[[act_vgene]] - length(ref_motif_vgenes_id_list[[act_vgene]])
        subsample_parts[[idx]] <- ref_motif_vgenes_id_list[[act_vgene]]
      }else {
        subsample_parts[[idx]] <- sample(x = ref_motif_vgenes_id_list[[act_vgene]],size = motif_region_vgenes_list[[act_vgene]],replace = FALSE)
      }
    }
    random_subsample <- unlist(subsample_parts)
    if(random_length > 0){
      random_subsample <- c(random_subsample, sample(x = (seq_along(refseqs_motif_region))[-random_subsample],size = random_length,replace = FALSE))
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }

  ### Return a reference subsample with biased V gene and CDR3b length distribution according to the sample
  if(vgene_stratify == TRUE && cdr3_len_stratify == TRUE){

    ## if there are more seqs with specific v gene or CDR3b length in sample than in reference database, surplus number of seqs will be taken randomly from remaining seqs
    random_length <- 0
    subsample_parts <- vector("list", length(motif_lengths_list) * length(motif_region_vgenes_list))
    idx <- 0L
    for(cdr3_length in names(motif_lengths_list)){
      for(act_vgene in names(motif_region_vgenes_list)){
        idx <- idx + 1L
        if(length(ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]]) < lengths_vgenes_list[[cdr3_length]][[act_vgene]]){
          random_length <- random_length + lengths_vgenes_list[[cdr3_length]][[act_vgene]] - length(ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]])
          subsample_parts[[idx]] <- ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]]
        }else {
          subsample_parts[[idx]] <- sample(x = ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]],size = lengths_vgenes_list[[cdr3_length]][[act_vgene]],replace = FALSE)
        }
      }
    }
    random_subsample <- unlist(subsample_parts)

    if(random_length > 0){
      random_subsample <- c(random_subsample, sample(x = (seq_along(refseqs_motif_region))[-random_subsample],size = random_length,replace = FALSE))
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }

  ## Closing time!
  return(random_subsample)
}
