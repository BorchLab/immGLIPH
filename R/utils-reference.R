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
        # Named reference database
        reference_list <- NULL
        utils::data("reference_list", envir = environment(), package = "immGLIPH")
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
#' substitution score is >= 0. When \pkg{immApex} is installed, the vector is
#' derived from the full BLOSUM62 matrix in
#' \code{immapex_blosum.pam.matrices}; otherwise the bundled
#' \code{BlosumVec} dataset is loaded as a fallback.
#'
#' @return Character vector of AA pair strings (e.g. \code{"AA"}, \code{"AS"}).
#' @keywords internal
.get_blosum_vec <- function() {
    if (requireNamespace("immApex", quietly = TRUE)) {
        mat <- immApex::immapex_blosum.pam.matrices[["BLOSUM62"]]
        idx <- which(mat >= 0, arr.ind = TRUE)
        return(paste0(rownames(mat)[idx[, 1]], colnames(mat)[idx[, 2]]))
    }
    # Fallback: load from bundled package data
    BlosumVec <- NULL
    utils::data("BlosumVec", envir = environment(), package = "immGLIPH")
    BlosumVec
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
