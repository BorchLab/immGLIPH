#' Local similarity detection using Fisher's exact test
#'
#' Implements the GLIPH2-style Fisher's exact test approach for detecting
#' locally enriched CDR3 motifs in a sample set compared to a reference
#' database.
#'
#' @param motif_region Character vector. Motif regions extracted from sample
#'   sequences (e.g. CDR3b with structural boundaries trimmed).
#' @param refseqs_motif_region Character vector. Motif regions extracted from
#'   reference sequences.
#' @param seqs Character vector. Unique sample CDR3b sequences.
#' @param refseqs Data frame. Reference database containing at least a
#'   \code{CDR3b} column.
#' @param sequences Data frame. Sample data containing at least a \code{CDR3b}
#'   column.
#' @param motif_length Numeric vector. Lengths of motifs to search for (e.g.
#'   \code{2:4}).
#' @param kmer_mindepth Numeric. Minimum number of times a motif must be
#'   observed in the sample set to be considered.
#' @param lcminp Numeric. Maximum p-value threshold for a motif to be
#'   considered significant.
#' @param lcminove Numeric or numeric vector. Minimum fold change threshold(s)
#'   for enrichment filtering. If a vector, each element corresponds to the
#'   matching entry in \code{motif_length}.
#' @param discontinuous_motifs Logical. Whether to include discontinuous motifs
#'   in the search.
#' @param motif_distance_cutoff Numeric. Maximum positional distance for motifs
#'   to be grouped together. Not used directly in this function but kept for
#'   interface consistency.
#' @param no_cores Numeric. Number of cores to use for parallel motif finding.
#' @param verbose Logical. If \code{TRUE}, print status messages via
#'   \code{message()}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{selected_motifs}{Data frame of motifs passing significance and
#'     enrichment filters.}
#'   \item{all_motifs}{Data frame of all motifs found in the sample set with
#'     associated statistics.}
#' }
#'
#' @import foreach
#' @keywords internal
.local_fisher <- function(motif_region,
                          refseqs_motif_region,
                          seqs,
                          refseqs,
                          sequences,
                          motif_length,
                          kmer_mindepth,
                          lcminp,
                          lcminove,
                          discontinuous_motifs,
                          motif_distance_cutoff,
                          no_cores,
                          verbose) {

  ## -----------------------------------------------------------

  ## Step 1 & 2: Find motifs in reference and sample sets
  ## -----------------------------------------------------------
  if (requireNamespace("immApex", quietly = TRUE) &&
      exists("calculateMotif", asNamespace("immApex"))) {
    ## ---- immApex C++ fast path: single-call with multithreading ----
    calcMotif <- get("calculateMotif", asNamespace("immApex"))
    if (verbose) message("Finding motifs in reference sequences (immApex)...")
    ref_motifs_df <- calcMotif(
      input.sequences = refseqs_motif_region,
      motif.lengths   = motif_length,
      min.depth       = 1L,
      discontinuous   = discontinuous_motifs,
      nthreads        = no_cores
    )
    colnames(ref_motifs_df)[colnames(ref_motifs_df) == "frequency"] <- "count"

    if (verbose) message("Finding motifs in sample sequences (immApex)...")
    motifs_df <- calcMotif(
      input.sequences = motif_region,
      motif.lengths   = motif_length,
      min.depth       = 1L,
      discontinuous   = discontinuous_motifs,
      nthreads        = no_cores
    )
    colnames(motifs_df)[colnames(motifs_df) == "frequency"] <- "count"

  } else {
    ## ---- Fallback: foreach chunking with stringdist backend ----
    if (verbose) message("Finding motifs in reference sequences...")

    # Divide the sequences equally among all cores
    overhang <- length(refseqs_motif_region) %% no_cores
    id_list  <- list()
    last_id  <- 0
    next_id  <- 0
    steps    <- (length(refseqs_motif_region) - overhang) / no_cores

    for (i in seq_len(no_cores)) {
      next_id <- last_id + steps
      if (overhang > 0) {
        next_id <- next_id + 1
        overhang <- overhang - 1
      }
      id_list[[i]] <- (last_id + 1):next_id
      last_id <- next_id
    }

    # Receive all motifs in the reference sequences
    ref_motifs_list <- foreach::foreach(i = seq_len(no_cores)) %dopar% {
      return(findMotifs(seqs = refseqs_motif_region[id_list[[i]]],
                        q = motif_length,
                        discontinuous = discontinuous_motifs))
    }

    # Convert the list into a more manageable data frame
    ref_motifs_df <- NULL
    for (i in seq_len(no_cores)) {
      if (i == 1) {
        ref_motifs_df <- ref_motifs_list[[i]]
      } else {
        ref_motifs_df <- merge(x = ref_motifs_df,
                               y = ref_motifs_list[[i]],
                               by = "motif",
                               all = TRUE)
      }
    }
    ref_motifs_df[is.na(ref_motifs_df)] <- 0
    ref_motifs_df <- data.frame(
      motif = ref_motifs_df$motif,
      count = rowSums(data.frame(ref_motifs_df[, -1]))
    )

    if (verbose) message("Finding motifs in sample sequences...")

    # Divide the sequences equally among all cores
    overhang <- length(motif_region) %% no_cores
    id_list  <- list()
    last_id  <- 0
    next_id  <- 0
    steps    <- (length(motif_region) - overhang) / no_cores

    for (i in seq_len(no_cores)) {
      next_id <- last_id + steps
      if (overhang > 0) {
        next_id <- next_id + 1
        overhang <- overhang - 1
      }
      id_list[[i]] <- (last_id + 1):next_id
      last_id <- next_id
    }

    # Receive all motifs in the sample sequences
    motifs_list <- foreach::foreach(i = seq_len(no_cores)) %dopar% {
      return(findMotifs(seqs = motif_region[id_list[[i]]],
                        q = motif_length,
                        discontinuous = discontinuous_motifs))
    }

    # Convert the list into a more manageable data frame
    motifs_df <- NULL
    for (i in seq_len(no_cores)) {
      if (i == 1) {
        motifs_df <- motifs_list[[i]]
      } else {
        motifs_df <- merge(x = motifs_df,
                           y = motifs_list[[i]],
                           by = "motif",
                           all = TRUE)
      }
    }
    motifs_df[is.na(motifs_df)] <- 0
    motifs_df <- data.frame(
      motif = motifs_df$motif,
      count = rowSums(data.frame(motifs_df[, -1]))
    )
  }

  ## -----------------------------------------------------------
  ## Step 3: Compare motifs between sample and reference
  ## -----------------------------------------------------------
  if (verbose) message("Comparing motifs between sample and reference sets...")

  # Summarise motifs found in the sample sequences and their counts in both sets
  motifs_df <- merge(x = motifs_df, y = ref_motifs_df, by = "motif", all.x = TRUE)
  motifs_df[is.na(motifs_df)] <- 0
  colnames(motifs_df) <- c("motif", "counts", "num_in_ref")
  motifs_df$avgRef <- motifs_df$num_in_ref
  motifs_df$topRef <- rep(0, nrow(motifs_df))
  motifs_df$OvE   <- rep(0, nrow(motifs_df))

  # Assess significance by Fisher's Exact Test (one-sided, greater) using
  # the hypergeometric distribution.
  #
  # Contingency table:
  # -----------------------------------------------------------------------
  #   # sample seqs with motif    | # reference seqs with motif
  # -----------------------------------------------------------------------
  #   # sample seqs without motif | # reference seqs without motif
  # -----------------------------------------------------------------------
  n_unique_seqs    <- length(unique(seqs))
  n_unique_refseqs <- length(unique(refseqs))

  motifs_df$p.value <- stats::phyper(
    q          = motifs_df$counts - 1,
    m          = motifs_df$num_in_ref + motifs_df$counts,
    n          = n_unique_refseqs + n_unique_seqs - motifs_df$num_in_ref - motifs_df$counts,
    k          = n_unique_seqs,
    lower.tail = FALSE
  )
  motifs_df$p.value <- as.numeric(formatC(motifs_df$p.value,
                                          digits = 1, format = "e"))

  # Fold change enrichment of motif in sample vs reference
  motifs_df$OvE <- round(
    motifs_df$counts / (motifs_df$num_in_ref + 0.01) /
      n_unique_seqs * n_unique_refseqs,
    digits = 1
  )

  # For comparability with the resampling method, normalise reference counts
  # to the size of the sample set
  motifs_df$avgRef <- round(
    motifs_df$num_in_ref / length(refseqs_motif_region) * length(motif_region),
    digits = 3
  )

  ## -----------------------------------------------------------
  ## Step 4: Determine minimum fold change per motif length
  ## -----------------------------------------------------------
  temp_minove <- c()
  if (length(lcminove) == 1) {
    temp_minove <- rep(lcminove, nrow(motifs_df))
  } else {
    temp_minove <- rep(0, nrow(motifs_df))

    # The number of informative positions in discontinuous motifs is reduced
    # by one compared to continuous motifs of the same length
    nam_len <- nchar(motifs_df$motif)
    disc_idx <- grep(".", motifs_df$motif, fixed = TRUE)
    nam_len[disc_idx] <- nam_len[disc_idx] - 1

    for (i in seq_along(motif_length)) {
      temp_minove[nam_len == motif_length[i]] <- lcminove[i]
    }
  }

  ## -----------------------------------------------------------
  ## Step 5: Filter significantly enriched motifs
  ## -----------------------------------------------------------
  # Criteria:
  #   - p-value <= lcminp (Fisher's exact test)
  #   - fold change >= lcminove (length-dependent)
  #   - sample occurrence >= kmer_mindepth
  sel_mot_df <- motifs_df[motifs_df$OvE >= temp_minove &
                            motifs_df$p.value <= lcminp &
                            motifs_df$counts >= kmer_mindepth, ]

  if (verbose) {
    message(nrow(sel_mot_df), " significantly enriched motif(s) identified.")
  }

  ## -----------------------------------------------------------
  ## Return results
  ## -----------------------------------------------------------
  list(
    selected_motifs = sel_mot_df,
    all_motifs      = motifs_df
  )
}
