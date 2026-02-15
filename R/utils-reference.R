#' Valid built-in reference database names
#'
#' @return Character vector of valid reference database names.
#' @keywords internal
.valid_reference_names <- function() {
    c("human_v1.0_CD4", "human_v1.0_CD8", "human_v1.0_CD48",
      "human_v2.0_CD4", "human_v2.0_CD8", "human_v2.0_CD48",
      "mouse_v1.0_CD4", "mouse_v1.0_CD8", "mouse_v1.0_CD48",
      "gliph_reference")
}

#' Get or download the immGLIPH reference list
#'
#' Downloads the reference repertoire data from GitHub releases on first use
#' and caches locally via \pkg{BiocFileCache}. Subsequent calls load from the
#' cache without re-downloading.
#'
#' The cached file contains a named \code{list} with entries for each built-in
#' reference database (see \code{\link{.valid_reference_names}}).
#'
#' @param force_download Logical. If \code{TRUE}, re-download even if cached.
#'   **Default:** \code{FALSE}
#' @param verbose Logical. Print messages. **Default:** \code{TRUE}
#'
#' @return A named \code{list} of reference databases. Each element is a list
#'   with \code{refseqs}, \code{vgene_frequencies}, and
#'   \code{cdr3_length_frequencies}.
#'
#' @examples
#' \dontrun{
#' ref <- getGLIPHreference()
#' names(ref)
#' head(ref[["human_v2.0_CD48"]]$refseqs)
#' }
#'
#' @export
getGLIPHreference <- function(force_download = FALSE, verbose = TRUE) {
    if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
        stop("Package 'BiocFileCache' is required to download reference data.\n",
             "Install with: BiocManager::install('BiocFileCache')",
             call. = FALSE)
    }

    cache_dir <- tools::R_user_dir("immGLIPH", which = "cache")
    bfc <- BiocFileCache::BiocFileCache(cache_dir, ask = FALSE)

    rname <- "immGLIPH_reference_list_v1"
    file_url <- paste0(
        "https://github.com/BorchLab/immGLIPH/releases/download/",
        "reference_data/reference_list.RData"
    )

    rid <- BiocFileCache::bfcquery(bfc, rname, "rname", exact = TRUE)$rid

    if (force_download && length(rid)) {
        BiocFileCache::bfcremove(bfc, rid)
        rid <- character(0)
    }

    if (!length(rid)) {
        if (verbose) message("Downloading immGLIPH reference data ...")
        rid <- names(BiocFileCache::bfcadd(bfc, rname, file_url))
    } else {
        if (verbose) message("Loading immGLIPH reference data from cache.")
    }

    file_path <- BiocFileCache::bfcrpath(bfc, rids = rid)

    env <- new.env(parent = emptyenv())
    load(file_path, envir = env)
    env[["reference_list"]]
}

#' Load reference list (internal, with caching)
#'
#' Loads the reference list, using a package-level cache to avoid repeated
#' disk reads or downloads within a single session.
#'
#' @return The reference list.
#' @keywords internal
.get_reference_list <- function() {
    # Session-level cache: store in package namespace environment
    pkg_env <- parent.env(environment())
    if (!is.null(pkg_env$.reference_list_cache)) {
        return(pkg_env$.reference_list_cache)
    }
    ref <- getGLIPHreference(verbose = TRUE)
    pkg_env$.reference_list_cache <- ref
    ref
}

# Session-level cache placeholder (populated on first use)
.reference_list_cache <- NULL

#' Load and prepare reference database
#'
#' Handles both named reference databases and user-provided data frames.
#'
#' @param refdb_beta Character name or data frame
#' @param accept_CF Filter for C/F start/end
#' @param min_seq_length Minimum sequence length
#' @param global_vgene Whether V-gene info is needed
#' @param vgene_stratify Whether V-gene stratification is needed
#' @param structboundaries Whether to use structural boundaries
#' @param boundary_size Boundary size
#' @param verbose Print messages
#' @return A list with refseqs (character vector of CDR3b), ref_vgenes,
#'   refseqs_motif_region, and the full reference data frame
#' @keywords internal
.load_reference <- function(refdb_beta,
                            accept_CF = TRUE,
                            min_seq_length = 8,
                            global_vgene = FALSE,
                            vgene_stratify = FALSE,
                            structboundaries = TRUE,
                            boundary_size = 3,
                            verbose = TRUE) {

    ref_vgenes <- NULL
    refseqs_df <- NULL

    if (is.data.frame(refdb_beta)) {
        refseqs <- refdb_beta
        refseqs[] <- lapply(refseqs, as.character)

        if (nrow(refseqs) <= 0) {
            stop("Reference repertoire is empty.", call. = FALSE)
        }

        if ("CDR3b" %in% colnames(refseqs)) {
            if (verbose) message("Using 'CDR3b' column from reference database.")
        } else {
            if (verbose) message("Using first column as CDR3b sequences.")
            colnames(refseqs)[1] <- "CDR3b"
        }

        if (ncol(refseqs) > 1 && global_vgene) {
            if ("TRBV" %in% colnames(refseqs)) {
                if (verbose) message("Using 'TRBV' column from reference database.")
            } else {
                if (verbose) message("Using second column as V-gene information.")
                colnames(refseqs)[2] <- "TRBV"
            }
        } else if (global_vgene) {
            stop("V-gene restriction requires V-gene info in reference database.",
                 call. = FALSE)
        }

        if (!("TRBV" %in% colnames(refseqs))) {
            refseqs$TRBV <- ""
        }
        if (ncol(refseqs) == 1) {
            refseqs <- cbind(refseqs, data.frame(TRBV = rep("", nrow(refseqs)),
                                                  stringsAsFactors = FALSE))
        }

        refseqs <- refseqs[, c("CDR3b", "TRBV")]

        if (accept_CF) {
            refseqs <- refseqs[grep("^C.*F$", refseqs$CDR3b, perl = TRUE), ]
        }
        if (nrow(refseqs) == 0) {
            stop("No reference sequences matching C/F criteria.", call. = FALSE)
        }

        refseqs <- unique(refseqs)
        refseqs <- refseqs[nchar(refseqs$CDR3b) >= min_seq_length, ]
        if (nrow(refseqs) == 0) {
            stop("No reference sequences >= min_seq_length.", call. = FALSE)
        }

        refseqs <- refseqs[grep("^[ACDEFGHIKLMNPQRSTVWY]*$", refseqs$CDR3b), ]
        if (nrow(refseqs) == 0) {
            stop("No valid amino acid sequences in reference.", call. = FALSE)
        }

        refseqs_df <- refseqs

        if (vgene_stratify) {
            ref_vgenes <- as.character(refseqs$TRBV)
        }

        refseqs_cdr3 <- as.character(refseqs$CDR3b)

    } else {
        # Named reference database - load via BiocFileCache
        reference_list <- .get_reference_list()
        refseqs <- as.data.frame(reference_list[[refdb_beta]]$refseqs)
        refseqs[] <- lapply(refseqs, as.character)

        if (accept_CF) {
            refseqs <- refseqs[grep("^C.*F$", refseqs$CDR3b, perl = TRUE), ]
        }
        refseqs <- unique(refseqs)
        refseqs <- refseqs[nchar(refseqs$CDR3b) >= min_seq_length, ]
        refseqs <- refseqs[grep("^[ACDEFGHIKLMNPQRSTVWY]*$", refseqs$CDR3b), ]

        refseqs_df <- refseqs

        if (vgene_stratify) {
            ref_vgenes <- as.character(refseqs$TRBV)
            if (verbose) message("Using V-gene information from reference database.")
        }

        refseqs_cdr3 <- as.character(refseqs$CDR3b)
    }

    # Prepare motif region
    refseqs_motif_region <- .prepare_motif_region(refseqs_cdr3,
                                                   structboundaries,
                                                   boundary_size)

    list(
        refseqs       = refseqs_cdr3,
        ref_vgenes    = ref_vgenes,
        refseqs_motif = refseqs_motif_region,
        refseqs_df    = refseqs_df
    )
}

#' Get BLOSUM62-compatible amino acid pairs
#'
#' Returns a character vector of two-letter amino acid pairs whose BLOSUM62
#' substitution score is >= 0, derived from the full BLOSUM62 matrix in
#' \code{immApex::immapex_blosum.pam.matrices}.
#'
#' @return Character vector of AA pair strings (e.g. \code{"AA"}, \code{"AS"}).
#' @keywords internal
.get_blosum_vec <- function() {
    env <- new.env(parent = emptyenv())
    utils::data("immapex_blosum.pam.matrices", package = "immApex",
                envir = env)
    mat <- env[["immapex_blosum.pam.matrices"]][["BLOSUM62"]]
    idx <- which(mat >= 0, arr.ind = TRUE)
    paste0(rownames(mat)[idx[, 1]], colnames(mat)[idx[, 2]])
}

#' Prepare motif regions from sequences
#'
#' Extracts the motif region by removing boundary amino acids.
#'
#' @param seqs Character vector of sequences
#' @param structboundaries Whether to trim boundaries
#' @param boundary_size Number of AAs to trim from each end
#' @return Character vector of motif regions
#' @keywords internal
.prepare_motif_region <- function(seqs, structboundaries, boundary_size) {
    if (structboundaries) {
        substr(seqs, boundary_size + 1L, nchar(seqs) - boundary_size)
    } else {
        seqs
    }
}
