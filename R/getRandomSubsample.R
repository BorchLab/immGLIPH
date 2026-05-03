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
#' @examples
#' ref_seqs <- c("ASSG", "ASSD", "ASSE", "ASSF", "ASSK", "ASSL")
#' sample_seqs <- c("ASSG", "ASSF", "ASSL")
#' sub <- getRandomSubsample(
#'   refseqs_motif_region = ref_seqs,
#'   motif_region = sample_seqs
#' )
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
                               lengths_vgenes_list) {
  random_subsample <- c()

  ### Return an unbiased reference subsample
  if (!cdr3_len_stratify && !vgene_stratify) {
    random_subsample <- refseqs_motif_region[
      sample(
        x = seq_along(refseqs_motif_region),
        size = length(motif_region),
        replace = FALSE
      )
    ]
  }

  ### Return a reference subsample with biased CDR3b length distribution
  if (cdr3_len_stratify && !vgene_stratify) {
    random_length <- 0
    subsample_parts <- vector("list", length(motif_lengths_list))
    idx <- 0L
    for (cdr3_length in names(motif_lengths_list)) {
      idx <- idx + 1L
      ref_ids <- ref_motif_lengths_id_list[[cdr3_length]]
      needed <- motif_lengths_list[[cdr3_length]]
      if (length(ref_ids) < needed) {
        random_length <- random_length + needed - length(ref_ids)
        subsample_parts[[idx]] <- ref_ids
      } else {
        subsample_parts[[idx]] <- sample(
          x = ref_ids, size = needed, replace = FALSE
        )
      }
    }
    random_subsample <- unlist(subsample_parts)
    if (random_length > 0) {
      random_subsample <- c(
        random_subsample,
        sample(
          x = seq_along(refseqs_motif_region)[-random_subsample],
          size = random_length,
          replace = FALSE
        )
      )
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }

  ### Return a reference subsample with biased V gene distribution
  if (vgene_stratify && !cdr3_len_stratify) {
    random_length <- 0
    subsample_parts <- vector("list", length(motif_region_vgenes_list))
    idx <- 0L
    for (act_vgene in names(motif_region_vgenes_list)) {
      idx <- idx + 1L
      ref_ids <- ref_motif_vgenes_id_list[[act_vgene]]
      needed <- motif_region_vgenes_list[[act_vgene]]
      if (length(ref_ids) < needed) {
        random_length <- random_length + needed - length(ref_ids)
        subsample_parts[[idx]] <- ref_ids
      } else {
        subsample_parts[[idx]] <- sample(
          x = ref_ids, size = needed, replace = FALSE
        )
      }
    }
    random_subsample <- unlist(subsample_parts)
    if (random_length > 0) {
      random_subsample <- c(
        random_subsample,
        sample(
          x = seq_along(refseqs_motif_region)[-random_subsample],
          size = random_length,
          replace = FALSE
        )
      )
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }

  ### Return a reference subsample with biased V gene and CDR3b length distribution
  if (vgene_stratify && cdr3_len_stratify) {
    random_length <- 0
    subsample_parts <- vector(
      "list",
      length(motif_lengths_list) * length(motif_region_vgenes_list)
    )
    idx <- 0L
    for (cdr3_length in names(motif_lengths_list)) {
      for (act_vgene in names(motif_region_vgenes_list)) {
        idx <- idx + 1L
        ref_ids <- ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]]
        needed <- lengths_vgenes_list[[cdr3_length]][[act_vgene]]
        if (length(ref_ids) < needed) {
          random_length <- random_length + needed - length(ref_ids)
          subsample_parts[[idx]] <- ref_ids
        } else {
          subsample_parts[[idx]] <- sample(
            x = ref_ids, size = needed, replace = FALSE
          )
        }
      }
    }
    random_subsample <- unlist(subsample_parts)

    if (random_length > 0) {
      random_subsample <- c(
        random_subsample,
        sample(
          x = seq_along(refseqs_motif_region)[-random_subsample],
          size = random_length,
          replace = FALSE
        )
      )
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }

  return(random_subsample)
}
