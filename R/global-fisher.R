#' Global similarity search using structural matching and Fisher's test (GLIPH2)
#'
#' Identifies globally similar CDR3b sequences by generating "struct" tags
#' (motif region with one variable position replaced by \code{"\%"}), then
#' testing for enrichment of each struct in the sample set vs. a naive
#' reference database using the hypergeometric distribution (one-sided
#' Fisher's exact test). Optionally filters for BLOSUM62-compatible amino
#' acid substitutions at the variable position.
#'
#' @param seqs Character vector. Unique CDR3b sequences from the sample.
#' @param motif_region Character vector. Motif regions of the sample
#'   sequences (boundaries already stripped if applicable).
#' @param sequences Data frame. Full sample data frame with at least columns
#'   \code{CDR3b} and (if applicable) \code{TRBV}.
#' @param refseqs Data frame. Reference database with at least a \code{CDR3b}
#'   column.
#' @param refseqs_motif_region Character vector. Motif regions of the reference
#'   sequences.
#' @param structboundaries Logical. Whether structural boundaries are applied.
#' @param boundary_size Integer. Number of AAs trimmed from each end.
#' @param global_vgene Logical. If \code{TRUE}, restrict edges to pairs sharing
#'   a V-gene.
#' @param all_aa_interchangeable Logical. If \code{FALSE}, only pairs whose
#'   variable-position amino acids have BLOSUM62 >= 0 are kept.
#' @param no_cores Integer. Number of registered parallel cores.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with elements:
#' \describe{
#'   \item{cluster_list}{Data frame of struct clusters with columns
#'     \code{cluster_tag}, \code{cluster_size}, \code{unique_CDR3b},
#'     \code{num_in_ref}, \code{fisher.score}, \code{aa_at_position},
#'     \code{TRBV}, \code{CDR3b}.}
#'   \item{global_similarities}{Logical indicating whether any global
#'     similarities were found.}
#' }
#'
#' @import foreach
#' @keywords internal
.global_fisher <- function(seqs,
                           motif_region,
                           sequences,
                           refseqs,
                           refseqs_motif_region,
                           structboundaries,
                           boundary_size,
                           global_vgene,
                           all_aa_interchangeable,
                           no_cores,
                           verbose) {

  ## Load BLOSUM62 compatible amino acid pairs
  BlosumVec <- .get_blosum_vec()

  ## Empty result template
  empty_df <- data.frame(
    cluster_tag     = character(0),
    cluster_size    = integer(0),
    unique_CDR3b    = integer(0),
    num_in_ref      = integer(0),
    fisher.score    = numeric(0),
    aa_at_position  = character(0),
    TRBV            = character(0),
    CDR3b           = character(0),
    stringsAsFactors = FALSE
  )

  ## -----------------------------------------------------------------
  ## Step 1: Build struct tags for sample sequences
  ## -----------------------------------------------------------------
  if (verbose) message("Searching for global similarities (Fisher method).")

  sample_seqs <- data.frame(
    seq    = seqs,
    struct = motif_region,
    stringsAsFactors = FALSE
  )
  sample_seqs$nchar      <- nchar(sample_seqs$struct)
  sample_seqs$pos        <- 0L
  sample_seqs$aa_at_pos  <- ""
  sample_seqs$tag        <- sample_seqs$struct

  ## Expand: for each position, replace that position with "%"
  exp_sample_seqs <- foreach::foreach(
    i = seq_len(max(sample_seqs$nchar)),
    .combine = rbind
  ) %dopar% {
    temp_part <- sample_seqs[sample_seqs$nchar >= i, ]
    temp_part$pos <- i
    temp_part$aa_at_pos <- substr(temp_part$struct, i, i)
    substr(temp_part$tag, i, i) <- "%"

    ## Keep only structs that occur at least twice
    dup_check <- duplicated(temp_part$tag) |
                 duplicated(temp_part$tag, fromLast = TRUE)
    temp_part[dup_check, ]
  }

  ## Early exit if nothing duplicated
  if (is.null(exp_sample_seqs) || nrow(exp_sample_seqs) == 0) {
    if (verbose) message("No global similarities found in sample set.")
    return(list(cluster_list = empty_df, global_similarities = FALSE))
  }

  unq_sample_struct <- unique(exp_sample_seqs$tag)

  if (length(unq_sample_struct) == 0) {
    if (verbose) message("No global similarities found in sample set.")
    return(list(cluster_list = empty_df, global_similarities = FALSE))
  }

  ## -----------------------------------------------------------------
  ## Step 2: Build struct tags for reference sequences
  ## -----------------------------------------------------------------
  reference_seqs <- data.frame(
    seq = unique(refseqs$CDR3b),
    stringsAsFactors = FALSE
  )
  if (structboundaries) {
    reference_seqs$struct <- substr(
      reference_seqs$seq,
      boundary_size + 1L,
      nchar(reference_seqs$seq) - boundary_size
    )
  } else {
    reference_seqs$struct <- reference_seqs$seq
  }
  reference_seqs$nchar      <- nchar(reference_seqs$struct)
  reference_seqs$pos        <- 0L
  reference_seqs$aa_at_pos  <- ""
  reference_seqs$tag        <- reference_seqs$struct

  exp_reference_seqs <- foreach::foreach(
    i = seq_len(min(max(reference_seqs$nchar), max(sample_seqs$nchar))),
    .combine = rbind
  ) %dopar% {
    temp_part <- reference_seqs[reference_seqs$nchar >= i, ]
    temp_part$pos <- i
    temp_part$aa_at_pos <- substr(temp_part$struct, i, i)
    substr(temp_part$tag, i, i) <- "%"
    temp_part[temp_part$tag %in% unq_sample_struct, ]
  }

  ## -----------------------------------------------------------------
  ## Step 3: Compute struct frequencies and enrichment
  ## -----------------------------------------------------------------
  sample_stats <- data.frame(table(exp_sample_seqs$tag),
                             stringsAsFactors = FALSE)
  colnames(sample_stats) <- c("tag", "Freq")

  if (is.null(exp_reference_seqs) || nrow(exp_reference_seqs) == 0) {
    ref_stats <- data.frame(tag = character(0), Freq = integer(0),
                            stringsAsFactors = FALSE)
  } else {
    ref_stats <- data.frame(table(exp_reference_seqs$tag),
                            stringsAsFactors = FALSE)
  }
  colnames(ref_stats) <- c("tag", "Freq")

  sample_stats <- merge(sample_stats, ref_stats, by = "tag", all.x = TRUE)
  sample_stats$tag <- as.character(sample_stats$tag)
  sample_stats[is.na(sample_stats)] <- 0
  colnames(sample_stats) <- c("tag", "num_in_sample", "num_in_ref")

  ## -----------------------------------------------------------------
  ## Step 4: Add V-gene info and build edges
  ## -----------------------------------------------------------------
  seqs_w_vgenes <- sequences
  if (!global_vgene) {
    seqs_w_vgenes$TRBV <- rep(" ", nrow(seqs_w_vgenes))
  }
  exp_sample_seqs <- merge(
    exp_sample_seqs,
    seqs_w_vgenes[, c("CDR3b", "TRBV")],
    by.x = "seq", by.y = "CDR3b",
    all.x = TRUE
  )

  edges <- foreach::foreach(
    i = seq_len(nrow(sample_stats)),
    .combine = rbind
  ) %dopar% {
    act_members <- unique(
      exp_sample_seqs[exp_sample_seqs$tag == sample_stats$tag[i], ]
    )
    act_members$summary <- paste(
      act_members$seq, act_members$tag,
      act_members$aa_at_pos, act_members$pos,
      act_members$TRBV,
      sep = "#"
    )

    if (nrow(act_members) >= 2) {
      combn_ids <- t(utils::combn(seq_len(nrow(act_members)), m = 2))
    } else {
      combn_ids <- t(utils::combn(rep(1, 2), m = 2))
    }

    ret <- data.frame(
      V1 = act_members$summary[combn_ids[, 1]],
      V2 = act_members$summary[combn_ids[, 2]],
      stringsAsFactors = FALSE
    )

    keep <- rep(TRUE, nrow(ret))
    if (global_vgene) {
      keep <- keep &
        (act_members$TRBV[combn_ids[, 1]] == act_members$TRBV[combn_ids[, 2]])
    }
    if (!all_aa_interchangeable) {
      keep <- keep &
        (paste0(act_members$aa_at_pos[combn_ids[, 1]],
                act_members$aa_at_pos[combn_ids[, 2]]) %in% BlosumVec)
    }
    ret[keep, ]
  }

  ## -----------------------------------------------------------------
  ## Step 5: Build cluster network
  ## -----------------------------------------------------------------
  if (is.null(edges) || nrow(edges) == 0) {
    if (verbose) message("No global similarities found in sample set.")
    return(list(cluster_list = empty_df, global_similarities = FALSE))
  }

  gr <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
  cm <- igraph::components(gr)

  if (cm$no == 0) {
    if (verbose) message("No global similarities found in sample set.")
    return(list(cluster_list = empty_df, global_similarities = FALSE))
  }

  ## Parse the component membership back into columns
  cm_splitted <- data.frame(
    matrix(unlist(strsplit(names(cm$membership), split = "#")),
           ncol = 5, byrow = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(cm_splitted) <- c("CDR3b", "tag", "aa_at_position",
                              "position", "TRBV")
  cm_splitted$position <- as.numeric(cm_splitted$position)

  cluster_list_raw <- foreach::foreach(i = seq_len(cm$no)) %dopar% {
    csize     <- cm$csize[i]
    member_df <- unique(cm_splitted[which(cm$membership == i), ])
    num_cdr3s <- length(unique(member_df$CDR3b))

    c(member_df$tag[1], csize, num_cdr3s,
      paste(unique(member_df$CDR3b), collapse = " "),
      paste(sort(unique(member_df$aa_at_position)), collapse = ""),
      member_df$TRBV[1])
  }

  cluster_list <- data.frame(
    matrix(unlist(cluster_list_raw), ncol = 6, byrow = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(cluster_list) <- c("cluster_tag", "cluster_size", "unique_CDR3b",
                               "CDR3b", "aa_at_position", "TRBV")
  cluster_list$cluster_size <- as.numeric(cluster_list$cluster_size)
  cluster_list$unique_CDR3b <- as.numeric(cluster_list$unique_CDR3b)

  ## Add reference frequency
  cluster_list <- merge(
    cluster_list,
    sample_stats[, c("tag", "num_in_ref")],
    by.x = "cluster_tag", by.y = "tag",
    all.x = TRUE
  )
  cluster_list$num_in_ref[is.na(cluster_list$num_in_ref)] <- 0

  ## Fisher's exact test (hypergeometric)
  n_uniq_ref    <- length(unique(refseqs$CDR3b))
  n_uniq_sample <- length(unique(seqs))

  cluster_list$fisher.score <- stats::phyper(
    q          = cluster_list$unique_CDR3b - 1,
    m          = cluster_list$num_in_ref + cluster_list$unique_CDR3b,
    n          = n_uniq_ref + n_uniq_sample -
                   cluster_list$num_in_ref - cluster_list$unique_CDR3b,
    k          = n_uniq_sample,
    lower.tail = FALSE
  )

  cluster_list <- cluster_list[, c("cluster_tag", "cluster_size",
                                    "unique_CDR3b", "num_in_ref",
                                    "fisher.score", "aa_at_position",
                                    "TRBV", "CDR3b")]

  if (verbose) {
    message(nrow(cluster_list),
            " clusters based on global similarity found in sample set.")
  }

  list(cluster_list = cluster_list, global_similarities = TRUE)
}
