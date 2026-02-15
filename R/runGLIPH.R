#' Run the GLIPH or GLIPH2 TCR clustering algorithm
#'
#' Unified entry point for the GLIPH/GLIPH2 algorithm for grouping T cell
#' receptors by antigen specificity. The function identifies locally and
#' globally similar CDR3b sequences, clusters them into convergence groups,
#' and scores each group for biological relevance.
#'
#' @param cdr3_sequences Input data containing CDR3b amino acid sequences.
#'   Accepts a character vector, a \code{data.frame} with columns described
#'   below, a \code{Seurat} object, a \code{SingleCellExperiment} object, or
#'   a list returned by \code{scRepertoire::combineTCR()/combineBCR()}.
#'
#'   When a \code{data.frame} is supplied, the following column names are
#'   recognized (alternative names in parentheses are mapped automatically):
#'   \describe{
#'     \item{CDR3b}{(cdr3, cdr3_aa, CDR3.beta, junction_aa) Required.
#'       CDR3 beta-chain amino acid sequences.}
#'     \item{TRBV}{(v_gene, v.gene, Vgene, v_call) Optional. V-gene usage.}
#'     \item{patient}{(sample, donor, sample_id) Optional. Donor index.}
#'     \item{HLA}{(hla, HLA_alleles) Optional. HLA alleles, comma-separated.}
#'     \item{counts}{(frequency, clone_count, cloneCount) Optional. Clone
#'       frequency.}
#'   }
#' @param method Character. Algorithm preset to use.
#'   \itemize{
#'     \item{\code{"gliph2"}}{ Fisher-based local and global similarity,
#'       GLIPH2-style isolated clustering and scoring.}
#'     \item{\code{"gliph1"}}{ Repeated random sampling for local similarity,
#'       Hamming distance cutoff for global similarity, GLIPH1-style
#'       connected-component clustering.}
#'     \item{\code{"custom"}}{ All parameters can be set independently.}
#'   }
#'   **Default:** \code{"gliph2"}
#' @param chains Character. Chain type for extraction from \code{Seurat} or
#'   \code{SingleCellExperiment} objects via \code{immApex::getIR()}.
#'   **Default:** \code{"TRB"}
#' @param result_folder Character. Path to output folder. If \code{""},
#'   results are not saved to disk. **Default:** \code{""}
#' @param refdb_beta Character or \code{data.frame}. Reference database for
#'   motif enrichment testing. Built-in databases include
#'   \code{"human_v1.0_CD4"}, \code{"human_v1.0_CD8"},
#'   \code{"human_v1.0_CD48"}, \code{"human_v2.0_CD4"},
#'   \code{"human_v2.0_CD8"}, \code{"human_v2.0_CD48"},
#'   \code{"mouse_v1.0_CD4"}, \code{"mouse_v1.0_CD8"},
#'   \code{"mouse_v1.0_CD48"}, and the legacy alias
#'   \code{"gliph_reference"} (= \code{"human_v1.0_CD48"}).
#'   Alternatively, supply a \code{data.frame} with CDR3b in the first
#'   column and optional V-gene in the second.
#'   See \code{\link{reference_list}} for details.
#'   **Default:** \code{"human_v2.0_CD48"}
#' @param v_usage_freq \code{data.frame} or \code{NULL}. V-gene frequencies
#'   for scoring. If \code{NULL}, built-in defaults are used.
#'   **Default:** \code{NULL}
#' @param cdr3_length_freq \code{data.frame} or \code{NULL}. CDR3 length
#'   frequencies for scoring. If \code{NULL}, built-in defaults are used.
#'   **Default:** \code{NULL}
#' @param ref_cluster_size Character. Reference cluster size strategy.
#'   \itemize{
#'     \item{\code{"original"}}{ Use the original sample size.}
#'     \item{\code{"simulated"}}{ Use simulated cluster sizes.}
#'   }
#'   **Default:** \code{"original"}
#' @param sim_depth Integer. Simulation depth for repeated random sampling
#'   (local method \code{"rrs"}) or cluster scoring.
#'   **Default:** \code{1000}
#' @param lcminp Numeric. Local convergence maximum p-value threshold.
#'   **Default:** \code{0.01}
#' @param lcminove Numeric vector. Local convergence minimum fold-change
#'   per motif length (lengths 2, 3, and 4 respectively).
#'   **Default:** \code{c(1000, 100, 10)}
#' @param kmer_mindepth Integer. Minimum number of kmer observations
#'   required to consider a motif.
#'   **Default:** \code{3}
#' @param accept_CF Logical. If \code{TRUE}, accept only sequences starting
#'   with C and ending with F.
#'   **Default:** \code{TRUE}
#' @param min_seq_length Integer. Minimum CDR3b sequence length to retain.
#'   **Default:** \code{8}
#' @param gccutoff Numeric or \code{NULL}. Global convergence Hamming distance
#'   cutoff (used when \code{global_method = "cutoff"}). If \code{NULL},
#'   the cutoff is auto-selected based on sample size.
#'   **Default:** \code{NULL}
#' @param structboundaries Logical. If \code{TRUE}, trim structural
#'   boundaries from CDR3b sequences before motif search.
#'   **Default:** \code{TRUE}
#' @param boundary_size Integer. Number of positions to trim from each end
#'   when \code{structboundaries = TRUE}.
#'   **Default:** \code{3}
#' @param motif_length Numeric vector. Motif lengths to search.
#'   **Default:** \code{c(2, 3, 4)}
#' @param local_similarities Logical. If \code{TRUE}, search for locally
#'   similar CDR3b sequences.
#'   **Default:** \code{TRUE}
#' @param global_similarities Logical. If \code{TRUE}, search for globally
#'   similar CDR3b sequences.
#'   **Default:** \code{TRUE}
#' @param local_method Character or \code{NULL}. Method for local similarity
#'   detection. If \code{NULL}, set by the \code{method} preset.
#'   \itemize{
#'     \item{\code{"fisher"}}{ Fisher exact test for motif enrichment.}
#'     \item{\code{"rrs"}}{ Repeated random sampling.}
#'   }
#'   **Default:** \code{NULL}
#' @param global_method Character or \code{NULL}. Method for global similarity
#'   detection. If \code{NULL}, set by the \code{method} preset.
#'   \itemize{
#'     \item{\code{"fisher"}}{ Fisher exact test for struct enrichment.}
#'     \item{\code{"cutoff"}}{ Hamming distance cutoff.}
#'   }
#'   **Default:** \code{NULL}
#' @param clustering_method Character or \code{NULL}. Clustering strategy.
#'   If \code{NULL}, set by the \code{method} preset.
#'   \itemize{
#'     \item{\code{"GLIPH1.0"}}{ Connected-component clustering.}
#'     \item{\code{"GLIPH2.0"}}{ Isolated clustering with merging.}
#'   }
#'   **Default:** \code{NULL}
#' @param scoring_method Character or \code{NULL}. Scoring strategy.
#'   If \code{NULL}, set by the \code{method} preset.
#'   \itemize{
#'     \item{\code{"GLIPH1.0"}}{ GLIPH1-style scoring.}
#'     \item{\code{"GLIPH2.0"}}{ GLIPH2-style scoring.}
#'   }
#'   **Default:** \code{NULL}
#' @param cluster_min_size Integer. Minimum number of unique CDR3b sequences
#'   required to retain a convergence group.
#'   **Default:** \code{2}
#' @param hla_cutoff Numeric. Significance cutoff for HLA enrichment testing.
#'   **Default:** \code{0.1}
#' @param n_cores Integer or \code{NULL}. Number of cores for parallel
#'   processing. If \code{NULL}, the number of available cores is
#'   auto-detected.
#'   **Default:** \code{1}
#' @param motif_distance_cutoff Integer. Maximum positional distance between
#'   shared motifs for two CDR3b sequences to be linked (GLIPH2).
#'   **Default:** \code{3}
#' @param discontinuous_motifs Logical. If \code{TRUE}, allow discontinuous
#'   motif patterns during local similarity search.
#'   **Default:** \code{FALSE}
#' @param all_aa_interchangeable Logical. If \code{FALSE}, BLOSUM62 filtering
#'   is applied to global similarities, restricting substitutions to
#'   biochemically similar amino acids.
#'   **Default:** \code{FALSE}
#' @param boost_local_significance Logical. If \code{TRUE}, boost local
#'   p-values using germline N-nucleotide insertion information.
#'   **Default:** \code{FALSE}
#' @param global_vgene Logical. If \code{TRUE}, restrict global similarity
#'   edges to pairs sharing the same V-gene.
#'   **Default:** \code{FALSE}
#' @param cdr3_len_stratify Logical. If \code{TRUE}, stratify random
#'   subsamples by CDR3 length (used with \code{local_method = "rrs"}).
#'   **Default:** \code{FALSE}
#' @param vgene_stratify Logical. If \code{TRUE}, stratify random subsamples
#'   by V-gene usage (used with \code{local_method = "rrs"}).
#'   **Default:** \code{FALSE}
#' @param public_tcrs Logical or character. Controls cross-donor edge
#'   filtering. For \code{method = "gliph1"} or \code{"gliph2"}: if
#'   \code{FALSE}, restrict edges to same donor. For \code{method = "custom"}:
#'   \itemize{
#'     \item{\code{"all"}}{ Allow cross-donor edges for all similarity types.}
#'     \item{\code{"local"}}{ Allow cross-donor edges for local only.}
#'     \item{\code{"global"}}{ Allow cross-donor edges for global only.}
#'     \item{\code{"none"}}{ Restrict all edges to same donor.}
#'   }
#'   **Default:** \code{TRUE}
#' @param vgene_match Character. V-gene matching requirement for custom
#'   clustering.
#'   \itemize{
#'     \item{\code{"none"}}{ No V-gene matching required.}
#'     \item{\code{"local"}}{ Require V-gene match for local edges.}
#'     \item{\code{"global"}}{ Require V-gene match for global edges.}
#'     \item{\code{"all"}}{ Require V-gene match for all edges.}
#'   }
#'   **Default:** \code{"none"}
#' @param scoring_sim_depth Integer. Simulation depth used specifically for
#'   convergence group scoring.
#'   **Default:** \code{1000}
#' @param verbose Logical. If \code{TRUE}, print progress messages to the
#'   console.
#'   **Default:** \code{TRUE}
#'
#' @return A \code{list} with the following elements:
#' \describe{
#'   \item{\code{sample_log}}{\code{data.frame}. Motif counts per simulation
#'     iteration (only present when \code{local_method = "rrs"}).}
#'   \item{\code{motif_enrichment}}{\code{list} with two elements:
#'     \describe{
#'       \item{\code{selected_motifs}}{\code{data.frame} of significantly
#'         enriched motifs passing all thresholds.}
#'       \item{\code{all_motifs}}{\code{data.frame} of all evaluated motifs
#'         with enrichment statistics.}
#'     }}
#'   \item{\code{global_enrichment}}{\code{list}. Global struct enrichment
#'     results (GLIPH2 only; \code{NULL} otherwise).}
#'   \item{\code{connections}}{\code{data.frame}. Edge list representing the
#'     clone network.}
#'   \item{\code{cluster_properties}}{\code{data.frame}. Convergence group
#'     properties and scores.}
#'   \item{\code{cluster_list}}{Named \code{list} of \code{data.frame}
#'     objects with per-cluster member details.}
#'   \item{\code{parameters}}{\code{list}. All input parameters used for the
#'     run.}
#' }
#'
#' @references
#' Glanville, J. et al. (2017). Identifying specificity groups in the
#' T cell receptor repertoire. \emph{Nature}, 547, 94--98.
#' \doi{10.1038/nature22976}
#'
#' Huang, H. et al. (2020). Analyzing the Mycobacterium tuberculosis immune
#' response by T-cell receptor clustering with GLIPH2 and genome-wide antigen
#' screening. \emph{Nature Biotechnology}, 38, 1194--1202.
#' \doi{10.1038/s41587-020-0505-4}
#'
#' @examples
#' utils::data("gliph_input_data")
#' res <- runGLIPH(
#'   cdr3_sequences = gliph_input_data[seq_len(200), ],
#'   method = "gliph2",
#'   sim_depth = 50,
#'   n_cores = 1
#' )
#'
#' @import foreach
#' @export
runGLIPH <- function(cdr3_sequences,
                     method = c("gliph2", "gliph1", "custom"),
                     chains = "TRB",
                     result_folder = "",
                     refdb_beta = "human_v2.0_CD48",
                     v_usage_freq = NULL,
                     cdr3_length_freq = NULL,
                     ref_cluster_size = "original",
                     sim_depth = 1000,
                     lcminp = 0.01,
                     lcminove = c(1000, 100, 10),
                     kmer_mindepth = 3,
                     accept_CF = TRUE,
                     min_seq_length = 8,
                     gccutoff = NULL,
                     structboundaries = TRUE,
                     boundary_size = 3,
                     motif_length = c(2, 3, 4),
                     local_similarities = TRUE,
                     global_similarities = TRUE,
                     local_method = NULL,
                     global_method = NULL,
                     clustering_method = NULL,
                     scoring_method = NULL,
                     cluster_min_size = 2,
                     hla_cutoff = 0.1,
                     n_cores = 1,
                     motif_distance_cutoff = 3,
                     discontinuous_motifs = FALSE,
                     all_aa_interchangeable = FALSE,
                     boost_local_significance = FALSE,
                     global_vgene = FALSE,
                     cdr3_len_stratify = FALSE,
                     vgene_stratify = FALSE,
                     public_tcrs = TRUE,
                     vgene_match = "none",
                     scoring_sim_depth = 1000,
                     verbose = TRUE) {

  t1 <- Sys.time()

  ## ---- Match method preset ----

  method <- match.arg(method)

  if (is.null(local_method)) {
    local_method <- switch(method,
      gliph1 = "rrs",
      gliph2 = "fisher",
      custom = "rrs"
    )
  }
  if (is.null(global_method)) {
    global_method <- switch(method,
      gliph1 = "cutoff",
      gliph2 = "fisher",
      custom = "cutoff"
    )
  }
  if (is.null(clustering_method)) {
    clustering_method <- switch(method,
      gliph1 = "GLIPH1.0",
      gliph2 = "GLIPH2.0",
      custom = "GLIPH1.0"
    )
  }
  if (is.null(scoring_method)) {
    scoring_method <- switch(method,
      gliph1 = "GLIPH1.0",
      gliph2 = "GLIPH2.0",
      custom = "GLIPH1.0"
    )
  }

  ## ---- Extract input ----
  cdr3_sequences <- .extract_input(cdr3_sequences, chains = chains)

  ## ---- Validate parameters ----
  params <- .validate_params(
    refdb_beta              = refdb_beta,
    v_usage_freq            = v_usage_freq,
    cdr3_length_freq        = cdr3_length_freq,
    ref_cluster_size        = ref_cluster_size,
    sim_depth               = sim_depth,
    lcminp                  = lcminp,
    lcminove                = lcminove,
    kmer_mindepth           = kmer_mindepth,
    accept_CF               = accept_CF,
    min_seq_length          = min_seq_length,
    gccutoff                = gccutoff,
    structboundaries        = structboundaries,
    boundary_size           = boundary_size,
    motif_length            = motif_length,
    local_similarities      = local_similarities,
    global_similarities     = global_similarities,
    cluster_min_size        = cluster_min_size,
    hla_cutoff              = hla_cutoff,
    n_cores                 = n_cores,
    motif_distance_cutoff   = motif_distance_cutoff,
    discontinuous_motifs    = discontinuous_motifs,
    all_aa_interchangeable  = all_aa_interchangeable,
    boost_local_significance = boost_local_significance
  )

  ## Unpack validated params
  min_seq_length   <- params$min_seq_length
  boundary_size    <- params$boundary_size
  sim_depth        <- params$sim_depth
  kmer_mindepth    <- params$kmer_mindepth
  motif_length     <- params$motif_length
  cluster_min_size <- params$cluster_min_size
  n_cores          <- params$n_cores

  ## ---- Prepare result folder ----
  result_folder <- .prepare_result_folder(result_folder)
  save_results <- result_folder != ""

  ## ---- Parse sequences ----
  parsed <- .parse_sequences(
    cdr3_sequences  = cdr3_sequences,
    accept_CF       = accept_CF,
    min_seq_length  = min_seq_length,
    global_vgene    = global_vgene,
    vgene_stratify  = vgene_stratify,
    verbose         = verbose
  )
  sequences    <- parsed$sequences
  vgene.info   <- parsed$vgene.info
  patient.info <- parsed$patient.info
  hla.info     <- parsed$hla.info
  count.info   <- parsed$count.info

  ## ---- Load reference ----
  ref_result <- .load_reference(
    refdb_beta     = refdb_beta,
    accept_CF      = accept_CF,
    min_seq_length = min_seq_length,
    verbose        = verbose
  )
  refseqs <- ref_result$refseqs

  ## ---- Prepare motif regions ----
  seqs <- unique(sequences$CDR3b)
  motif_region <- .prepare_motif_region(seqs, structboundaries, boundary_size)
  refseqs_motif_region <- .prepare_motif_region(
    unique(refseqs), structboundaries, boundary_size
  )
  ## Also need "all" motif region for combined approach
  all_motif_region <- .prepare_motif_region(
    sequences$CDR3b, structboundaries, boundary_size
  )

  if (length(refseqs_motif_region) < length(motif_region)) {
    stop("Reference database must have more sequences than input data.",
         call. = FALSE)
  }

  ## ---- Set up parallel ----
  no_cores <- .setup_parallel(n_cores)

  ## ================================================================
  ## Part 1: Local similarities
  ## ================================================================
  sample_log      <- NULL
  selected_motifs <- NULL
  all_motifs      <- NULL
  local_clone_network <- NULL

  if (local_similarities) {
    if (verbose) message("Part 1: Searching for local similarities.")

    if (local_method == "rrs") {
      ## Prepare stratification lists
      motif_lengths_list       <- list()
      ref_motif_lengths_id_list <- list()
      motif_region_vgenes_list <- list()
      ref_motif_vgenes_id_list <- list()
      lengths_vgenes_list      <- list()
      ref_lengths_vgenes_list  <- list()

      if (vgene_stratify) {
        for (act_vgene in sort(as.character(unique(sequences$TRBV)))) {
          motif_region_vgenes_list[[act_vgene]] <-
            sum(sequences$TRBV == act_vgene)
          ref_motif_vgenes_id_list[[act_vgene]] <-
            which(unique(refseqs$CDR3b) %in%
                    refseqs$CDR3b[refseqs$TRBV == act_vgene])
        }
      }
      if (cdr3_len_stratify) {
        motif_region_lengths     <- nchar(motif_region)
        ref_motif_region_lengths <- nchar(refseqs_motif_region)
        for (cdr3_length in sort(unique(motif_region_lengths))) {
          motif_lengths_list[[as.character(cdr3_length)]] <-
            sum(motif_region_lengths == cdr3_length)
          ref_motif_lengths_id_list[[as.character(cdr3_length)]] <-
            which(ref_motif_region_lengths == cdr3_length)
        }
      }
      if (vgene_stratify && cdr3_len_stratify) {
        for (cdr3_length in names(ref_motif_lengths_id_list)) {
          for (act_vgene in names(ref_motif_vgenes_id_list)) {
            lengths_vgenes_list[[cdr3_length]][[act_vgene]] <-
              sum(sequences$TRBV == act_vgene &
                    nchar(motif_region) == as.numeric(cdr3_length))
            ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]] <-
              ref_motif_lengths_id_list[[cdr3_length]][
                ref_motif_lengths_id_list[[cdr3_length]] %in%
                  ref_motif_vgenes_id_list[[act_vgene]]
              ]
          }
        }
      }

      rrs_result <- .local_rrs(
        motif_region              = motif_region,
        refseqs_motif_region      = refseqs_motif_region,
        seqs                      = seqs,
        sequences                 = sequences,
        motif_length              = motif_length,
        sim_depth                 = sim_depth,
        kmer_mindepth             = kmer_mindepth,
        lcminp                    = lcminp,
        lcminove                  = lcminove,
        discontinuous_motifs      = discontinuous_motifs,
        cdr3_len_stratify         = cdr3_len_stratify,
        vgene_stratify            = vgene_stratify,
        no_cores                  = no_cores,
        verbose                   = verbose,
        motif_lengths_list        = motif_lengths_list,
        ref_motif_lengths_id_list = ref_motif_lengths_id_list,
        motif_region_vgenes_list  = motif_region_vgenes_list,
        ref_motif_vgenes_id_list  = ref_motif_vgenes_id_list,
        lengths_vgenes_list       = lengths_vgenes_list,
        ref_lengths_vgenes_list   = ref_lengths_vgenes_list
      )

      sample_log      <- rrs_result$sample_log
      selected_motifs <- rrs_result$selected_motifs
      all_motifs      <- rrs_result$all_motifs

    } else if (local_method == "fisher") {
      fisher_result <- .local_fisher(
        motif_region           = motif_region,
        refseqs_motif_region   = refseqs_motif_region,
        seqs                   = seqs,
        refseqs                = refseqs,
        sequences              = sequences,
        motif_length           = motif_length,
        kmer_mindepth          = kmer_mindepth,
        lcminp                 = lcminp,
        lcminove               = lcminove,
        discontinuous_motifs   = discontinuous_motifs,
        motif_distance_cutoff  = motif_distance_cutoff,
        no_cores               = no_cores,
        verbose                = verbose
      )

      selected_motifs <- fisher_result$selected_motifs
      all_motifs      <- fisher_result$all_motifs
    }

    ## Build local clone network from selected motifs
    if (!is.null(selected_motifs) && nrow(selected_motifs) > 0) {
      local_clone_network <- foreach::foreach(
        i = seq_len(nrow(selected_motifs)),
        .combine = rbind
      ) %dopar% {
        all_ids <- grep(
          pattern = selected_motifs$motif[i],
          x       = all_motif_region,
          value   = FALSE
        )
        if (length(all_ids) >= 2) {
          combn_ids <- t(utils::combn(all_ids, m = 2))
        } else {
          combn_ids <- t(utils::combn(rep(1, 2), m = 2))
        }
        temp_df <- data.frame(
          V1   = combn_ids[, 1],
          V2   = combn_ids[, 2],
          type = rep("local", nrow(combn_ids)),
          tag  = rep(selected_motifs$motif[i], nrow(combn_ids)),
          stringsAsFactors = FALSE
        )
        temp_df
      }

      ## Apply motif distance cutoff
      if (!is.null(local_clone_network) && nrow(local_clone_network) > 0) {
        diffs <- abs(
          stringr::str_locate(
            sequences$CDR3b[local_clone_network$V1],
            local_clone_network$tag
          )[, 1] -
          stringr::str_locate(
            sequences$CDR3b[local_clone_network$V2],
            local_clone_network$tag
          )[, 1]
        )
        local_clone_network <- local_clone_network[
          diffs < motif_distance_cutoff, ]
      }

      if (verbose) {
        message(nrow(selected_motifs),
                " significantly enriched motif(s) found.")
      }
    } else {
      local_similarities <- FALSE
      if (verbose) message("No significantly enriched motifs found.")
    }
  }

  ## ================================================================
  ## Part 2: Global similarities
  ## ================================================================
  global_res      <- NULL
  global_clone_network <- NULL
  not_in_global_ids <- seq_along(seqs)

  if (global_similarities) {
    if (verbose) message("Part 2: Searching for global similarities.")

    if (global_method == "cutoff") {
      ## Auto-set gccutoff if NULL
      if (is.null(gccutoff)) {
        gccutoff <- if (length(all_motif_region) < 125) 2 else 1
      }

      cutoff_result <- .global_cutoff(
        seqs         = seqs,
        motif_region = motif_region,
        sequences    = sequences,
        gccutoff     = gccutoff,
        global_vgene = global_vgene,
        no_cores     = no_cores,
        verbose      = verbose
      )

      global_clone_network <- cutoff_result$edges
      not_in_global_ids    <- cutoff_result$not_in_global_ids

      if (nrow(global_clone_network) == 0) {
        global_similarities <- FALSE
      }

    } else if (global_method == "fisher") {
      fisher_global_result <- .global_fisher(
        seqs                    = seqs,
        motif_region            = motif_region,
        sequences               = sequences,
        refseqs                 = refseqs,
        refseqs_motif_region    = refseqs_motif_region,
        structboundaries        = structboundaries,
        boundary_size           = boundary_size,
        global_vgene            = global_vgene,
        all_aa_interchangeable  = all_aa_interchangeable,
        no_cores                = no_cores,
        verbose                 = verbose
      )

      global_res <- fisher_global_result$cluster_list
      if (!fisher_global_result$global_similarities) {
        global_similarities <- FALSE
      }
    }
  }

  ## ================================================================
  ## Part 3: Clustering
  ## ================================================================
  if (verbose) message("Part 3: Clustering sequences.")

  cluster_properties  <- NULL
  cluster_list        <- list()
  clone_network       <- NULL
  save_cluster_list_df <- NULL

  if (clustering_method == "GLIPH1.0") {
    ## Build combined clone network
    combined_network <- NULL

    ## Convert local IDs to CDR3b sequences
    if (!is.null(local_clone_network) && nrow(local_clone_network) > 0) {
      local_net <- data.frame(
        V1   = sequences$CDR3b[local_clone_network$V1],
        V2   = sequences$CDR3b[local_clone_network$V2],
        type = local_clone_network$type,
        stringsAsFactors = FALSE
      )
      combined_network <- local_net
    }

    if (!is.null(global_clone_network) && nrow(global_clone_network) > 0) {
      global_net <- global_clone_network[, c("V1", "V2", "type")]
      if (!is.null(combined_network)) {
        combined_network <- rbind(combined_network, global_net)
      } else {
        combined_network <- global_net
      }
    }

    gliph1_result <- .cluster_gliph1(
      clone_network     = combined_network,
      sequences         = sequences,
      not_in_global_ids = not_in_global_ids,
      seqs              = seqs,
      vgene.info        = vgene.info,
      patient.info      = patient.info,
      global_vgene      = global_vgene,
      public_tcrs       = public_tcrs,
      cluster_min_size  = cluster_min_size,
      verbose           = verbose
    )

    cluster_properties   <- gliph1_result$cluster_properties
    cluster_list         <- gliph1_result$cluster_list
    clone_network        <- gliph1_result$clone_network
    save_cluster_list_df <- gliph1_result$save_cluster_list_df

  } else if (clustering_method == "GLIPH2.0") {
    gliph2_result <- .cluster_gliph2(
      local_res                = if (!is.null(selected_motifs) &&
                                     nrow(selected_motifs) > 0) {
                                   ## Add start/stop/members columns needed
                                   ## by .cluster_gliph2
                                   local_res_df <- selected_motifs
                                   if (!("start" %in% colnames(local_res_df))) {
                                     ## Compute start/stop from motif positions
                                     motif_names <- local_res_df$motif
                                     starts <- vapply(motif_names, function(m) {
                                       pos <- stringr::str_locate(
                                         all_motif_region, m
                                       )
                                       min(pos[, 1], na.rm = TRUE)
                                     }, numeric(1))
                                     stops <- vapply(motif_names, function(m) {
                                       pos <- stringr::str_locate(
                                         all_motif_region, m
                                       )
                                       max(pos[, 1], na.rm = TRUE)
                                     }, numeric(1))
                                     local_res_df$start <- starts
                                     local_res_df$stop <- stops
                                   }
                                   if (!("num_in_sample" %in%
                                         colnames(local_res_df))) {
                                     local_res_df$num_in_sample <-
                                       local_res_df$counts
                                   }
                                   if (!("num_fold" %in%
                                         colnames(local_res_df))) {
                                     local_res_df$num_fold <- local_res_df$OvE
                                   }
                                   if (!("fisher.score" %in%
                                         colnames(local_res_df))) {
                                     local_res_df$fisher.score <-
                                       local_res_df$p.value
                                   }
                                   if (!("members" %in%
                                         colnames(local_res_df))) {
                                     local_res_df$members <- vapply(
                                       local_res_df$motif,
                                       function(m) {
                                         ids <- grep(m, all_motif_region)
                                         paste(
                                           unique(sequences$CDR3b[ids]),
                                           collapse = " "
                                         )
                                       },
                                       character(1)
                                     )
                                   }
                                   local_res_df
                                 } else {
                                   NULL
                                 },
      global_res               = global_res,
      sequences                = sequences,
      local_similarities       = local_similarities,
      global_similarities      = global_similarities,
      global_vgene             = global_vgene,
      all_aa_interchangeable   = all_aa_interchangeable,
      structboundaries         = structboundaries,
      boundary_size            = boundary_size,
      motif_distance_cutoff    = motif_distance_cutoff,
      cluster_min_size         = cluster_min_size,
      boost_local_significance = boost_local_significance,
      verbose                  = verbose
    )

    cluster_properties   <- gliph2_result$merged_clusters
    cluster_list         <- gliph2_result$cluster_list
    clone_network        <- gliph2_result$clone_network
    save_cluster_list_df <- gliph2_result$save_cluster_list_df
  }

  ## ================================================================
  ## Part 4: Scoring
  ## ================================================================
  if (verbose) message("Part 4: Scoring convergence groups.")

  if (!is.null(cluster_properties) && length(cluster_list) > 0) {
    sc_gliph_version <- if (scoring_method == "GLIPH1.0") 1L else 2L

    scoring_res <- clusterScoring(
      cluster_list    = cluster_list,
      cdr3_sequences  = cdr3_sequences,
      refdb_beta      = refdb_beta,
      v_usage_freq    = v_usage_freq,
      cdr3_length_freq = cdr3_length_freq,
      ref_cluster_size = ref_cluster_size,
      sim_depth       = scoring_sim_depth,
      gliph_version   = sc_gliph_version,
      hla_cutoff      = hla_cutoff,
      n_cores         = n_cores
    )

    cluster_properties <- cbind(cluster_properties, scoring_res)
    ## Reorder so members is last
    members_col <- "members"
    if (members_col %in% colnames(cluster_properties)) {
      other_cols <- setdiff(colnames(cluster_properties), members_col)
      cluster_properties <- cluster_properties[, c(other_cols, members_col)]
    }
  }

  ## ---- Stop parallel ----
  .stop_parallel()

  ## ================================================================
  ## Save results
  ## ================================================================
  if (save_results) {
    if (!is.null(sample_log)) {
      utils::write.table(
        sample_log,
        file = paste0(result_folder, "kmer_resample_", sim_depth, "_log.txt"),
        quote = FALSE, sep = "\t", row.names = TRUE
      )
    }
    if (!is.null(selected_motifs)) {
      utils::write.table(
        selected_motifs,
        file = paste0(result_folder, "local_similarities.txt"),
        quote = FALSE, sep = "\t", row.names = FALSE
      )
    }
    if (!is.null(all_motifs)) {
      utils::write.table(
        all_motifs,
        file = paste0(result_folder, "all_motifs.txt"),
        quote = FALSE, sep = "\t", row.names = FALSE
      )
    }
    if (!is.null(clone_network)) {
      utils::write.table(
        clone_network,
        file = paste0(result_folder, "clone_network.txt"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
      )
    }
    if (!is.null(save_cluster_list_df)) {
      utils::write.table(
        save_cluster_list_df,
        file = paste0(result_folder, "cluster_member_details.txt"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
      )
    }
    if (!is.null(cluster_properties)) {
      utils::write.table(
        cluster_properties,
        file = paste0(result_folder, "convergence_groups.txt"),
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
      )
    }

    ## Save parameters
    all_params <- list(
      method                   = method,
      local_method             = local_method,
      global_method            = global_method,
      clustering_method        = clustering_method,
      scoring_method           = scoring_method,
      result_folder            = result_folder,
      ref_cluster_size         = ref_cluster_size,
      sim_depth                = sim_depth,
      lcminp                   = lcminp,
      lcminove                 = lcminove,
      kmer_mindepth            = kmer_mindepth,
      accept_CF                = accept_CF,
      min_seq_length           = min_seq_length,
      gccutoff                 = gccutoff,
      structboundaries         = structboundaries,
      boundary_size            = boundary_size,
      motif_length             = motif_length,
      motif_distance_cutoff    = motif_distance_cutoff,
      discontinuous_motifs     = discontinuous_motifs,
      local_similarities       = local_similarities,
      global_similarities      = global_similarities,
      global_vgene             = global_vgene,
      all_aa_interchangeable   = all_aa_interchangeable,
      boost_local_significance = boost_local_significance,
      cluster_min_size         = cluster_min_size,
      hla_cutoff               = hla_cutoff,
      n_cores                  = n_cores
    )
    .save_parameters(all_params, result_folder)
  }

  ## ================================================================
  ## Coerce numeric columns
  ## ================================================================
  if (is.data.frame(selected_motifs)) {
    selected_motifs <- .coerce_numeric_cols(selected_motifs)
  }
  if (is.data.frame(all_motifs)) {
    all_motifs <- .coerce_numeric_cols(all_motifs)
  }
  if (is.data.frame(cluster_properties)) {
    cluster_properties <- .coerce_numeric_cols(cluster_properties)
  }
  if (is.list(cluster_list)) {
    cluster_list <- lapply(cluster_list, .coerce_numeric_cols)
  }

  ## ================================================================
  ## Assemble output
  ## ================================================================
  t2 <- Sys.time()
  dt <- t2 - t1
  if (verbose) message("Total time: ", round(dt, 2), " ", units(dt))

  output <- list(
    sample_log         = sample_log,
    motif_enrichment   = list(
      selected_motifs = selected_motifs,
      all_motifs      = all_motifs
    ),
    global_enrichment  = global_res,
    connections        = clone_network,
    cluster_properties = cluster_properties,
    cluster_list       = cluster_list,
    parameters         = list(
      method                   = method,
      local_method             = local_method,
      global_method            = global_method,
      clustering_method        = clustering_method,
      scoring_method           = scoring_method,
      result_folder            = result_folder,
      ref_cluster_size         = ref_cluster_size,
      sim_depth                = sim_depth,
      lcminp                   = lcminp,
      lcminove                 = lcminove,
      kmer_mindepth            = kmer_mindepth,
      accept_CF                = accept_CF,
      min_seq_length           = min_seq_length,
      gccutoff                 = gccutoff,
      structboundaries         = structboundaries,
      boundary_size            = boundary_size,
      motif_length             = motif_length,
      motif_distance_cutoff    = motif_distance_cutoff,
      discontinuous_motifs     = discontinuous_motifs,
      local_similarities       = local_similarities,
      global_similarities      = global_similarities,
      global_vgene             = global_vgene,
      all_aa_interchangeable   = all_aa_interchangeable,
      boost_local_significance = boost_local_significance,
      cluster_min_size         = cluster_min_size,
      hla_cutoff               = hla_cutoff,
      n_cores                  = n_cores
    )
  )

  if (verbose && save_results) {
    message("Results saved to: ", result_folder)
  }

  output
}


#' @rdname runGLIPH
#' @section Deprecated Functions:
#' \code{turbo_gliph}, \code{gliph2}, and \code{gliph_combined} are
#' deprecated aliases. Use \code{runGLIPH()} instead.
#' @export
turbo_gliph <- function(cdr3_sequences, ...) {
  .Deprecated("runGLIPH", package = "immGLIPH",
              msg = paste("turbo_gliph() is deprecated.",
                          "Use runGLIPH(..., method = 'gliph1') instead."))
  runGLIPH(cdr3_sequences, method = "gliph1", ...)
}

#' @rdname runGLIPH
#' @export
gliph2 <- function(cdr3_sequences, ...) {
  .Deprecated("runGLIPH", package = "immGLIPH",
              msg = paste("gliph2() is deprecated.",
                          "Use runGLIPH(..., method = 'gliph2') instead."))
  runGLIPH(cdr3_sequences, method = "gliph2", ...)
}

#' @rdname runGLIPH
#' @export
gliph_combined <- function(cdr3_sequences, ...) {
  .Deprecated("runGLIPH", package = "immGLIPH",
              msg = paste("gliph_combined() is deprecated.",
                          "Use runGLIPH(..., method = 'custom') instead."))
  runGLIPH(cdr3_sequences, method = "custom", ...)
}
