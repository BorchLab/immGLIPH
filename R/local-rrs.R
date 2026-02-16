#' Local similarity detection via repeated random sampling (GLIPH1-style)
#'
#' Identifies enriched motifs in TCR CDR3b sequences by comparing motif
#' frequencies in the sample set against repeated random subsamples drawn from
#' a naive reference database.
#'
#' @param motif_region Character vector. Motif regions extracted from sample
#'   sequences.
#' @param refseqs_motif_region Character vector. Motif regions extracted from
#'   reference sequences.
#' @param seqs Character vector. Unique sample CDR3b sequences.
#' @param sequences Data frame. Must contain columns \code{CDR3b} and
#'   \code{TRBV}.
#' @param motif_length Numeric vector. Lengths of motifs to search for.
#' @param sim_depth Numeric. Number of repeated random sampling iterations.
#' @param kmer_mindepth Numeric. Minimum number of times a motif must appear
#'   in the sample set to be considered.
#' @param lcminp Numeric. Maximum p-value threshold for a motif to be
#'   selected.
#' @param lcminove Numeric vector. Minimum fold-change (observed / expected)
#'   threshold(s). If a single value, the same threshold is applied to all
#'   motif lengths.
#' @param discontinuous_motifs Logical. Whether to include discontinuous
#'   motifs.
#' @param cdr3_len_stratify Logical. Whether to stratify random subsamples by
#'   CDR3 length distribution.
#' @param vgene_stratify Logical. Whether to stratify random subsamples by
#'   V-gene usage distribution.
#' @param no_cores Integer. Number of cores for parallel execution.
#' @param verbose Logical. If \code{TRUE}, emit progress messages.
#' @param motif_lengths_list List. Pre-computed CDR3 length counts from the
#'   sample.
#' @param ref_motif_lengths_id_list List. Pre-computed reference indices by
#'   CDR3 length.
#' @param motif_region_vgenes_list List. Pre-computed V-gene counts from the
#'   sample.
#' @param ref_motif_vgenes_id_list List. Pre-computed reference indices by
#'   V-gene.
#' @param lengths_vgenes_list List. Pre-computed joint CDR3-length / V-gene
#'   counts from the sample.
#' @param ref_lengths_vgenes_list List. Pre-computed reference indices by
#'   joint CDR3-length / V-gene.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{sample_log}{Data frame. Rows are \code{"Discovery"} followed by
#'       one row per simulation; columns are motifs. Values are motif counts.}
#'     \item{selected_motifs}{Data frame of motifs passing enrichment filters
#'       with columns \code{motif}, \code{counts}, \code{num_in_ref},
#'       \code{avgRef}, \code{topRef}, \code{OvE}, \code{p.value}.}
#'     \item{all_motifs}{Data frame with the same columns as
#'       \code{selected_motifs} but for every motif found.}
#'   }
#'
#' @import foreach
#' @keywords internal
.local_rrs <- function(motif_region,
                       refseqs_motif_region,
                       seqs,
                       sequences,
                       motif_length,
                       sim_depth,
                       kmer_mindepth,
                       lcminp,
                       lcminove,
                       discontinuous_motifs,
                       cdr3_len_stratify,
                       vgene_stratify,
                       no_cores,
                       verbose,
                       motif_lengths_list,
                       ref_motif_lengths_id_list,
                       motif_region_vgenes_list,
                       ref_motif_vgenes_id_list,
                       lengths_vgenes_list,
                       ref_lengths_vgenes_list) {

  ## ---- Step 1: Find motifs in the sample set ----------------------------
  if (verbose) message("Computing motif frequency in sample set.")
  discovery <- findMotifs(
    seqs         = motif_region,
    q            = motif_length,
    discontinuous = discontinuous_motifs
  )

  ## ---- Step 2: Repeated random sampling from the reference DB -----------
  if (verbose) message("Running ", sim_depth, " random sampling iterations.")

  res <- foreach::foreach(i = seq_len(sim_depth), .inorder = FALSE,
                          .packages = "immGLIPH") %dopar% {
    motif_sample <- getRandomSubsample(
      cdr3_len_stratify        = cdr3_len_stratify,
      vgene_stratify           = vgene_stratify,
      refseqs_motif_region     = refseqs_motif_region,
      motif_region             = motif_region,
      motif_lengths_list       = motif_lengths_list,
      ref_motif_lengths_id_list = ref_motif_lengths_id_list,
      motif_region_vgenes_list = motif_region_vgenes_list,
      ref_motif_vgenes_id_list = ref_motif_vgenes_id_list,
      ref_lengths_vgenes_list  = ref_lengths_vgenes_list,
      lengths_vgenes_list      = lengths_vgenes_list
    )

    sim <- findMotifs(
      seqs         = motif_sample,
      q            = motif_length,
      discontinuous = discontinuous_motifs
    )

    sim <- merge(discovery, sim, by = "motif", all.x = TRUE)
    sim$V1.y[is.na(sim$V1.y)] <- 0
    sim$V1.y
  }

  ## ---- Step 3: Build sample log -----------------------------------------
  if (verbose) message("Building sample log.")

  res    <- unlist(res)
  copies <- matrix(as.numeric(res), ncol = sim_depth, byrow = FALSE)

  sample_log <- cbind(discovery$motif, discovery$V1, copies)
  sample_log <- t(sample_log)

  mot        <- sample_log[1, ]
  sample_log <- sample_log[-1, ]
  sample_log <- matrix(as.numeric(sample_log), nrow = sim_depth + 1)
  colnames(sample_log) <- mot
  sample_log <- as.data.frame(sample_log)
  rownames(sample_log) <- c("Discovery", paste("sim", 0:(sim_depth - 1),
                                               sep = "-"))

  ## ---- Step 4: Evaluate motif enrichment --------------------------------
  if (verbose) message("Evaluating motif enrichment.")

  sample_reads <- sample_log
  actual_reads <- as.data.frame(sample_reads[1, ])
  colnames(actual_reads) <- colnames(sample_reads)
  sample_reads <- as.data.frame(sample_reads[-1, ])
  nam <- colnames(actual_reads)

  ## Build results data frame

  motifs_df <- data.frame(
    matrix(rep(0, ncol(sample_reads) * 7), ncol = 7),
    stringsAsFactors = FALSE
  )
  colnames(motifs_df) <- c("motif", "counts", "num_in_ref",
                           "avgRef", "topRef", "OvE", "p.value")
  motifs_df$motif  <- nam
  motifs_df$counts <- unlist(actual_reads)
  motifs_df$avgRef <- colMeans(sample_reads)
  motifs_df$topRef <- vapply(sample_reads, max, FUN.VALUE = c(1))

  motifs_df$OvE[motifs_df$avgRef > 0] <-
    motifs_df$counts[motifs_df$avgRef > 0] /
    motifs_df$avgRef[motifs_df$avgRef > 0]

  motifs_df$OvE[motifs_df$avgRef == 0] <-
    1 / (sim_depth * length(motif_region))

  motifs_df$p.value <- vapply(
    seq_len(ncol(sample_reads)),
    function(x) sum(sample_reads[, x] >= actual_reads[, x]) / sim_depth,
    FUN.VALUE = c(1)
  )

  motifs_df$avgRef  <- round(motifs_df$avgRef, digits = 2)
  motifs_df$OvE     <- round(motifs_df$OvE,    digits = 3)
  motifs_df$p.value <- round(motifs_df$p.value, digits = 6)
  motifs_df$p.value[motifs_df$p.value == 0] <- 1 / sim_depth

  ## ---- Step 5: Determine per-motif fold-change threshold ----------------
  if (length(lcminove) == 1) {
    temp_minove <- rep(lcminove, ncol(actual_reads))
  } else {
    temp_minove <- rep(0, ncol(actual_reads))

    # For discontinuous motifs the dot does not count as an informative
    # position, so effective length is one less.
    nam_len <- nchar(nam)
    nam_len[grep(".", nam, fixed = TRUE)] <-
      nam_len[grep(".", nam, fixed = TRUE)] - 1

    for (i in seq_along(motif_length)) {
      temp_minove[nam_len == motif_length[i]] <- lcminove[i]
    }
  }

  ## ---- Step 6: Select significantly enriched motifs ---------------------
  # Criteria:
  #   - p-value   <  lcminp

  #   - OvE       >= per-motif lcminove threshold
  #   - counts    >= kmer_mindepth
  sel_mot_df <- motifs_df[
    motifs_df$OvE     >= temp_minove &
    motifs_df$p.value <  lcminp &
    motifs_df$counts  >= kmer_mindepth,
  ]

  if (verbose) {
    if (nrow(sel_mot_df) == 0) {
      message("No significantly enriched motifs found in sample set.")
    } else {
      message(nrow(sel_mot_df),
              " significantly enriched motif(s) found in sample set.")
    }
  }

  ## ---- Return -----------------------------------------------------------
  list(
    sample_log      = sample_log,
    selected_motifs = sel_mot_df,
    all_motifs      = motifs_df
  )
}
