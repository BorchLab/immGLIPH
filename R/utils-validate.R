#' Validate parameters for runGLIPH
#'
#' Consolidates all parameter validation into one function.
#'
#' @param refdb_beta Reference database
#' @param v_usage_freq V-gene usage frequency data frame
#' @param cdr3_length_freq CDR3 length frequency data frame
#' @param ref_cluster_size Cluster size reference type
#' @param sim_depth Simulation depth
#' @param lcminp Local convergence min p-value
#' @param lcminove Local convergence min OvE
#' @param kmer_mindepth Minimum kmer observations
#' @param accept_CF Accept C/F start/end only
#' @param min_seq_length Minimum sequence length
#' @param gccutoff Global convergence cutoff
#' @param structboundaries Use structural boundaries
#' @param boundary_size Boundary size in AAs
#' @param motif_length Motif lengths to search
#' @param local_similarities Search local similarities
#' @param global_similarities Search global similarities
#' @param cluster_min_size Minimum cluster size
#' @param hla_cutoff HLA significance cutoff
#' @param n_cores Number of cores
#' @param motif_distance_cutoff Motif distance cutoff (GLIPH2)
#' @param discontinuous_motifs Allow discontinuous motifs
#' @param all_aa_interchangeable BLOSUM62 filtering
#' @param boost_local_significance Germline boost
#' @return List of validated (and possibly adjusted) parameter values
#' @keywords internal
.validate_params <- function(refdb_beta = "gliph_reference",
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
                             cluster_min_size = 2,
                             hla_cutoff = 0.1,
                             n_cores = 1,
                             motif_distance_cutoff = 3,
                             discontinuous_motifs = FALSE,
                             all_aa_interchangeable = FALSE,
                             boost_local_significance = FALSE) {

    # refdb_beta
    if (!is.data.frame(refdb_beta)) {
        if (length(refdb_beta) != 1 || !is.character(refdb_beta)) {
            stop("refdb_beta must be a data frame or 'gliph_reference'.",
                 call. = FALSE)
        }
        if (!(refdb_beta %in% c("gliph_reference"))) {
            stop("refdb_beta must be a data frame or 'gliph_reference'.",
                 call. = FALSE)
        }
    }

    # v_usage_freq
    if (!is.null(v_usage_freq)) {
        if (!is.data.frame(v_usage_freq) || ncol(v_usage_freq) < 2 ||
            nrow(v_usage_freq) < 1) {
            stop("v_usage_freq must be a data frame with V-genes in column 1 ",
                 "and frequencies in column 2.", call. = FALSE)
        }
        v_usage_freq[, 2] <- as.numeric(v_usage_freq[, 2])
        if (any(is.na(v_usage_freq[, 2]))) {
            stop("v_usage_freq column 2 must contain numeric frequencies.",
                 call. = FALSE)
        }
    }

    # cdr3_length_freq
    if (!is.null(cdr3_length_freq)) {
        if (!is.data.frame(cdr3_length_freq) || ncol(cdr3_length_freq) < 2 ||
            nrow(cdr3_length_freq) < 1) {
            stop("cdr3_length_freq must be a data frame with CDR3 lengths ",
                 "in column 1 and frequencies in column 2.", call. = FALSE)
        }
        cdr3_length_freq[, 2] <- as.numeric(cdr3_length_freq[, 2])
        if (any(is.na(cdr3_length_freq[, 2]))) {
            stop("cdr3_length_freq column 2 must contain numeric frequencies.",
                 call. = FALSE)
        }
    }

    # ref_cluster_size
    if (!(ref_cluster_size %in% c("original", "simulated"))) {
        stop("ref_cluster_size must be 'original' or 'simulated'.",
             call. = FALSE)
    }

    # sim_depth
    if (!is.numeric(sim_depth) || length(sim_depth) != 1 || sim_depth < 1) {
        stop("sim_depth must be a single number >= 1.", call. = FALSE)
    }
    sim_depth <- round(sim_depth)

    # lcminp
    if (!is.numeric(lcminp) || length(lcminp) != 1 || lcminp <= 0) {
        stop("lcminp must be a single number > 0.", call. = FALSE)
    }

    # lcminove
    if (!is.numeric(lcminove)) {
        stop("lcminove must be numeric.", call. = FALSE)
    }
    if (length(lcminove) > 1 && length(lcminove) != length(motif_length)) {
        stop("lcminove must be length 1 or same length as motif_length.",
             call. = FALSE)
    }
    if (any(lcminove < 1)) {
        stop("lcminove values must be >= 1.", call. = FALSE)
    }

    # kmer_mindepth
    if (!is.numeric(kmer_mindepth) || length(kmer_mindepth) != 1 ||
        kmer_mindepth < 1) {
        stop("kmer_mindepth must be a single number >= 1.", call. = FALSE)
    }
    kmer_mindepth <- round(kmer_mindepth)

    # accept_CF
    if (!is.logical(accept_CF)) {
        stop("accept_sequences_with_C_F_start_end must be logical.",
             call. = FALSE)
    }

    # min_seq_length
    if (!is.numeric(min_seq_length) || length(min_seq_length) != 1 ||
        min_seq_length < 0) {
        stop("min_seq_length must be a single number >= 0.", call. = FALSE)
    }
    min_seq_length <- round(min_seq_length)

    # structboundaries / boundary_size
    if (!is.logical(structboundaries)) {
        stop("structboundaries must be logical.", call. = FALSE)
    }
    if (!is.numeric(boundary_size) || length(boundary_size) != 1 ||
        boundary_size < 0) {
        stop("boundary_size must be a single number >= 0.", call. = FALSE)
    }
    boundary_size <- round(boundary_size)
    if (structboundaries) {
        min_seq_length <- max(min_seq_length, 2L * boundary_size + 1L)
    }

    # motif_length
    if (!is.numeric(motif_length) || any(motif_length < 1)) {
        stop("motif_length must be numeric with values >= 1.", call. = FALSE)
    }
    motif_length <- round(motif_length)

    # similarities
    if (!is.logical(local_similarities)) {
        stop("local_similarities must be logical.", call. = FALSE)
    }
    if (!is.logical(global_similarities)) {
        stop("global_similarities must be logical.", call. = FALSE)
    }
    if (!local_similarities && !global_similarities) {
        stop("At least one of local_similarities or global_similarities ",
             "must be TRUE.", call. = FALSE)
    }

    # gccutoff
    if (!is.null(gccutoff)) {
        if (!is.numeric(gccutoff) || length(gccutoff) != 1 || gccutoff < 0) {
            stop("gccutoff must be NULL or a single number >= 0.",
                 call. = FALSE)
        }
    }

    # cluster_min_size
    if (!is.numeric(cluster_min_size) || length(cluster_min_size) != 1 ||
        cluster_min_size < 1) {
        stop("cluster_min_size must be a single number >= 1.", call. = FALSE)
    }
    cluster_min_size <- round(cluster_min_size)

    # hla_cutoff
    if (!is.numeric(hla_cutoff) || length(hla_cutoff) != 1 ||
        hla_cutoff < 0 || hla_cutoff > 1) {
        stop("hla_cutoff must be between 0 and 1.", call. = FALSE)
    }

    # n_cores
    if (!is.null(n_cores)) {
        if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) {
            stop("n_cores must be NULL or a single number >= 1.",
                 call. = FALSE)
        }
        n_cores <- round(n_cores)
    }

    # motif_distance_cutoff
    if (!is.numeric(motif_distance_cutoff) ||
        length(motif_distance_cutoff) != 1) {
        stop("motif_distance_cutoff must be a single number.", call. = FALSE)
    }

    # discontinuous_motifs
    if (!is.logical(discontinuous_motifs)) {
        stop("discontinuous_motifs must be logical.", call. = FALSE)
    }

    # all_aa_interchangeable
    if (!is.logical(all_aa_interchangeable)) {
        stop("all_aa_interchangeable must be logical.", call. = FALSE)
    }

    # boost_local_significance
    if (!is.logical(boost_local_significance)) {
        stop("boost_local_significance must be logical.", call. = FALSE)
    }

    list(
        refdb_beta                = refdb_beta,
        v_usage_freq              = v_usage_freq,
        cdr3_length_freq          = cdr3_length_freq,
        ref_cluster_size          = ref_cluster_size,
        sim_depth                 = sim_depth,
        lcminp                    = lcminp,
        lcminove                  = lcminove,
        kmer_mindepth             = kmer_mindepth,
        accept_CF                 = accept_CF,
        min_seq_length            = min_seq_length,
        gccutoff                  = gccutoff,
        structboundaries          = structboundaries,
        boundary_size             = boundary_size,
        motif_length              = motif_length,
        local_similarities        = local_similarities,
        global_similarities       = global_similarities,
        cluster_min_size          = cluster_min_size,
        hla_cutoff                = hla_cutoff,
        n_cores                   = n_cores,
        motif_distance_cutoff     = motif_distance_cutoff,
        discontinuous_motifs      = discontinuous_motifs,
        all_aa_interchangeable    = all_aa_interchangeable,
        boost_local_significance  = boost_local_significance
    )
}
