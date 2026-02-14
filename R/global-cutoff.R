#' Global similarity search using Hamming distance cutoff (GLIPH1.0 method)
#'
#' @param seqs character vector. Unique CDR3b sequences filtered by C/F
#'   start/end (if applicable) and minimum length.
#' @param motif_region character vector. The motif region of each sequence in
#'   \code{seqs} (i.e. with boundaries stripped if \code{structboundaries} is
#'   active).
#' @param sequences data.frame. The full input data frame containing at least
#'   columns \code{CDR3b} and optionally \code{TRBV}.
#' @param gccutoff numeric. Maximum Hamming distance for two sequences to be
#'   considered globally similar.
#' @param global_vgene logical. If \code{TRUE}, global connections are
#'   restricted to sequence pairs that share a V gene.
#' @param no_cores numeric. Number of cores registered with the parallel
#'   backend (used only for documentation; the function relies on a
#'   pre-registered \code{foreach} backend).
#' @param verbose logical. If \code{TRUE}, progress messages are printed.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{edges}{A \code{data.frame} with columns \code{V1}, \code{V2}, and
#'     \code{type} (always \code{"global"}). Each row represents a pair of
#'     CDR3b sequences that are globally similar.}
#'   \item{not_in_global_ids}{An integer vector of indices (into \code{seqs})
#'     for sequences that have no global neighbour.}
#' }
#'
#' @keywords internal
#' @import foreach
.global_cutoff <- function(seqs,
                           motif_region,
                           sequences,
                           gccutoff,
                           global_vgene,
                           no_cores,
                           verbose) {

  if (verbose) message("Searching for global similarities (cutoff method).")

  ## ---- immApex C++ fast path ------------------------------------------------
  if (requireNamespace("immApex", quietly = TRUE)) {
    return(.global_cutoff_immapex(seqs, motif_region, sequences,
                                  gccutoff, global_vgene, verbose))
  }

  ## ---- Fallback: stringdist + foreach implementation ------------------------
  .global_cutoff_stringdist(seqs, motif_region, sequences,
                            gccutoff, global_vgene, no_cores, verbose)
}

#' immApex-accelerated global cutoff via buildNetwork()
#' @keywords internal
.global_cutoff_immapex <- function(seqs, motif_region, sequences,
                                   gccutoff, global_vgene, verbose) {

  ## Build input data frame when V-gene filtering is requested
  if (global_vgene) {
    ## Match V genes for each unique CDR3b
    vgene_lookup <- sequences[!duplicated(sequences$CDR3b),
                              c("CDR3b", "TRBV")]
    input_df <- data.frame(
      motif = motif_region,
      CDR3b = seqs,
      stringsAsFactors = FALSE
    )
    input_df <- merge(input_df, vgene_lookup, by = "CDR3b", all.x = TRUE)

    edge_result <- immApex::buildNetwork(
      input.data = input_df,
      seq_col    = "motif",
      v_col      = "TRBV",
      threshold  = gccutoff,
      filter.v   = TRUE,
      ids        = input_df$CDR3b,
      output     = "edges",
      metric     = "hamming"
    )
  } else {
    edge_result <- immApex::buildNetwork(
      input.sequences = motif_region,
      threshold       = gccutoff,
      ids             = seqs,
      output          = "edges",
      metric          = "hamming"
    )
  }

  ## Convert to expected format
  if (!is.null(edge_result) && nrow(edge_result) > 0) {
    edges <- data.frame(
      V1   = edge_result$from,
      V2   = edge_result$to,
      type = "global",
      stringsAsFactors = FALSE
    )
  } else {
    edges <- data.frame(
      V1   = character(0),
      V2   = character(0),
      type = character(0),
      stringsAsFactors = FALSE
    )
  }

  ## Identify isolated sequences (no global neighbour)
  connected <- unique(c(edges$V1, edges$V2))
  not_in_global_ids <- which(!seqs %in% connected)

  if (verbose) message(nrow(edges), " global edges found (cutoff method).")

  list(
    edges             = edges,
    not_in_global_ids = as.integer(not_in_global_ids)
  )
}

#' stringdist + foreach fallback for global cutoff
#' @keywords internal
.global_cutoff_stringdist <- function(seqs, motif_region, sequences,
                                      gccutoff, global_vgene, no_cores,
                                      verbose) {

  ## Prepare V-gene lookup vectors when V-gene filtering is requested
  temp_seqs   <- c()
  temp_vgenes <- c()
  if (global_vgene) {
    temp_seqs   <- seqs
    temp_vgenes <- sequences$TRBV[which(sequences$CDR3b %in% temp_seqs)]
  }

  ## ------------------------------------------------------------------
  ## Parallel loop: for every sequence compute Hamming distance to all
  ## others and identify neighbours within gccutoff
  ## ------------------------------------------------------------------
  res <- foreach::foreach(i = seq_along(seqs)) %dopar% {
    not_in_global_ids <- c()
    global_con        <- c()

    # Hamming distance between motif_region[i] and all other motif regions
    dis <- stringdist::stringdist(
      a      = motif_region[i],
      b      = motif_region,
      method = "hamming",
      nthread = 1
    )
    global_ids <- which(dis <= gccutoff)

    # If requested, restrict to shared V genes
    if (global_vgene) {
      act_vgenes         <- temp_vgenes[temp_seqs == seqs[i]]
      global_ids_seqs    <- seqs[global_ids]
      all_global_ids_seqs <- temp_seqs[temp_seqs %in% global_ids_seqs]
      global_ids <- global_ids[
        global_ids_seqs %in%
          all_global_ids_seqs[temp_vgenes[temp_seqs %in% global_ids_seqs] %in% act_vgenes]
      ]
    }

    # Track sequences that have no neighbour
    if (length(global_ids) == 1) not_in_global_ids <- i

    # Keep only the upper triangle (j > i) to avoid duplicate edges
    global_ids <- global_ids[which(global_ids > i)]

    if (length(global_ids) > 0) {
      global_con <- cbind(
        rep(seqs[i], length(global_ids)),
        seqs[global_ids],
        rep("global", length(global_ids))
      )
    }

    list(not_in_global_ids = not_in_global_ids, global_con = global_con)
  }

  ## ------------------------------------------------------------------
  ## Collect results from parallel workers
  ## ------------------------------------------------------------------
  not_in_global_ids <- c()
  global_con        <- c()
  for (i in seq_along(res)) {
    not_in_global_ids <- c(not_in_global_ids, res[[i]]$not_in_global_ids)
    global_con        <- rbind(global_con, res[[i]]$global_con)
  }

  ## Build edge data frame
  if (!is.null(global_con)) {
    edges <- data.frame(
      V1   = global_con[, 1],
      V2   = global_con[, 2],
      type = global_con[, 3],
      stringsAsFactors = FALSE
    )
  } else {
    edges <- data.frame(
      V1   = character(0),
      V2   = character(0),
      type = character(0),
      stringsAsFactors = FALSE
    )
  }

  if (verbose) message(nrow(edges), " global edges found (cutoff method).")

  list(
    edges             = edges,
    not_in_global_ids = as.integer(not_in_global_ids)
  )
}
