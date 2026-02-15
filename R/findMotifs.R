#' Find continuous and discontinuous sequence motifs
#'
#' Searches a character vector of amino acid sequences for k-mer motifs and
#' returns their frequencies. Both continuous and discontinuous (gapped) motifs
#' are supported. When \pkg{immApex} is installed, the C++-accelerated
#' \code{\link[immApex]{calculateMotif}} backend is used automatically for
#' improved performance; otherwise the function falls back to a pure-R
#' implementation based on \code{\link[stringdist]{qgrams}}.
#'
#' @param seqs A character vector of amino acid sequences in which motifs
#'   will be identified and counted.
#' @param q A numeric vector of motif lengths to search for.
#'   **Default:** \code{2:4}.
#' @param kmer_mindepth The minimum number of times a k-mer must be observed
#'   in \code{seqs} for it to be included in the output.
#'   **Default:** \code{NULL} (no filtering).
#' @param discontinuous Whether to include discontinuous (gapped) motifs in
#'   the search. **Default:** \code{FALSE}.
#'
#' @return A \code{data.frame} with two columns: \code{motif} (the k-mer
#'   string) and \code{V1} (the observed frequency).
#'
#' @examples
#' utils::data("gliph_input_data")
#' sample_seqs <- as.character(gliph_input_data$CDR3b)
#' res <- findMotifs(seqs = sample_seqs)
#'
#' @export
findMotifs <- function(seqs, q = 2:4, kmer_mindepth = NULL,
                       discontinuous = FALSE) {
  seqs <- as.character(seqs)

  ## ---- immApex C++ fast path ------------------------------------------------
  if (requireNamespace("immApex", quietly = TRUE) &&
      exists("calculateMotif", asNamespace("immApex"))) {
    result <- immApex::calculateMotif(
      input.sequences      = seqs,
      motif.lengths        = q,
      min.depth            = if (is.null(kmer_mindepth)) 1L else kmer_mindepth,
      discontinuous        = discontinuous,
      discontinuous.symbol = ".",
      nthreads             = 1L
    )
    ## Rename "frequency" -> "V1" for backward compatibility with callers
    colnames(result)[colnames(result) == "frequency"] <- "V1"
    return(result)
  }

  ## ---- Fallback: stringdist-based implementation ----------------------------
  .find_motifs_stringdist(seqs       = seqs,
                          q          = q,
                          kmer_mindepth = kmer_mindepth,
                          discontinuous = discontinuous)
}

#' @keywords internal
.find_motifs_stringdist <- function(seqs, q = 2:4, kmer_mindepth = NULL,
                                    discontinuous = FALSE) {
  kmer_parts <- list()
  seqs <- as.character(seqs)
  all_q <- q
  if (discontinuous == TRUE) all_q <- sort(unique(c(all_q, q + 1)))

  ### Iterate for all kmer lengths
  for (i in all_q) {
    ## Get all continuous motifs with length i
    cont_kmer <- stringdist::qgrams(seqs, q = i)
    cont_kmer <- as.data.frame(t(cont_kmer))
    cont_kmer$motif <- rownames(cont_kmer)

    ## Exclude motifs with less counts than kmer_mindepth
    if (!is.null(kmer_mindepth))
      cont_kmer <- cont_kmer[which(cont_kmer$V1 >= kmer_mindepth), ]
    cont_kmer <- cont_kmer[, 2:1]

    ## Save motifs and corresponding counts
    if (i %in% q) kmer_parts[[length(kmer_parts) + 1L]] <- cont_kmer

    ## Get all discontinuous motifs with (i-1) fixed positions with i in q
    if (discontinuous == TRUE && i %in% (q + 1)) {
      ## Find all discontinuous motifs with (i-1) fixed positions and
      ## variable position j
      for (j in 2:(i - 1)) {
        ## Replace position j by a dot in all continuous motifs with length i
        disc_kmer <- cont_kmer
        substr(disc_kmer$motif, j, j) <- "."
        disc_kmer <- as.data.frame(table(rep(disc_kmer$motif, disc_kmer$V1)))
        colnames(disc_kmer) <- c("motif", "V1")

        ## Exclude motifs with less counts than kmer_mindepth
        if (!is.null(kmer_mindepth))
          disc_kmer <- disc_kmer[which(disc_kmer$V1 >= kmer_mindepth), ]

        ## Save motifs and corresponding counts
        kmer_parts[[length(kmer_parts) + 1L]] <- disc_kmer
      }
    }
  }

  ## Closing time!
  all_kmers <- do.call(rbind, kmer_parts)
  return(all_kmers)
}
