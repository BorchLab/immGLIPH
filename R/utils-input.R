#' Extract TCR input data from various sources
#'
#' Handles Seurat objects, SingleCellExperiment objects,
#' combineTCR/combineBCR list output, data frames, and character vectors.
#' Uses immApex::getIR() when single-cell objects are detected.
#'
#' @param input Input data (data frame, vector, Seurat, SCE, or list)
#' @param chains Chain type for getIR extraction. Default "TRB".
#' @return A data frame with standardized column names
#' @keywords internal
.extract_input <- function(input, chains = "TRB") {

    # Seurat or SingleCellExperiment object
    if (inherits(input, "Seurat") || inherits(input, "SingleCellExperiment")) {
        if (!requireNamespace("immApex", quietly = TRUE)) {
            stop("The immApex package is required to extract TCR data from ",
                 class(input)[1], " objects.\n",
                 "Install with: devtools::install_github('BorchLab/immApex')",
                 call. = FALSE)
        }
        ir_data <- immApex::getIR(
            input.data = input,
            chains = chains,
            sequence.type = "aa"
        )
        df <- .standardize_colnames(ir_data)
        return(df)
    }

    # combineTCR/combineBCR output (list of data frames, not a data frame itself)
    if (is.list(input) && !is.data.frame(input)) {
        if (!requireNamespace("immApex", quietly = TRUE)) {
            stop("The immApex package is required for list input from ",
                 "combineTCR/combineBCR.\n",
                 "Install with: devtools::install_github('BorchLab/immApex')",
                 call. = FALSE)
        }
        ir_data <- immApex::getIR(
            input.data = input,
            chains = chains,
            sequence.type = "aa"
        )
        df <- .standardize_colnames(ir_data)
        return(df)
    }

    # Character vector of CDR3b sequences
    if (is.atomic(input) && is.character(input)) {
        return(data.frame(CDR3b = input, stringsAsFactors = FALSE))
    }

    # Data frame - standardize column names
    if (is.data.frame(input)) {
        return(.standardize_colnames(input))
    }

    stop("cdr3_sequences must be a character vector, data frame, ",
         "Seurat object, or SingleCellExperiment object.",
         call. = FALSE)
}

#' Standardize column names from various input formats
#'
#' Maps alternative column names to canonical GLIPH names.
#' Supports scRepertoire, immApex getIR, and native GLIPH formats.
#'
#' @param df A data frame
#' @return Data frame with standardized column names
#' @keywords internal
.standardize_colnames <- function(df) {
    col_map <- list(
        CDR3b   = c("CDR3b", "cdr3", "cdr3_aa", "CDR3.beta", "cdr3b",
                     "CDR3_beta", "junction_aa"),
        TRBV    = c("TRBV", "v_gene", "v.gene", "Vgene", "v_call"),
        patient = c("patient", "sample", "donor", "sample_id", "sample.id"),
        HLA     = c("HLA", "hla", "HLA_alleles"),
        counts  = c("counts", "frequency", "clone_count", "Frequency",
                     "cloneCount")
    )

    current_names <- colnames(df)
    new_names <- current_names

    for (canonical in names(col_map)) {
        if (canonical %in% current_names) next
        for (alt in col_map[[canonical]]) {
            idx <- which(current_names == alt)
            if (length(idx) > 0) {
                new_names[idx[1]] <- canonical
                break
            }
        }
    }

    colnames(df) <- new_names

    # If still no CDR3b column and only one column, assume it's CDR3b
    if (!("CDR3b" %in% colnames(df)) && ncol(df) == 1) {
        colnames(df)[1] <- "CDR3b"
    }

    df
}

#' Parse and filter sequences data frame
#'
#' Consolidates the input preparation logic from turbo_gliph, gliph2,
#' and gliph_combined into a single function.
#'
#' @param cdr3_sequences Data frame or vector of sequences
#' @param accept_CF Accept only sequences starting with C and ending with F
#' @param min_seq_length Minimum sequence length
#' @param global_vgene Whether global vgene matching is required
#' @param vgene_stratify Whether vgene stratification is required
#' @param verbose Print notifications
#' @return A list with sequences data frame and info flags
#' @keywords internal
.parse_sequences <- function(cdr3_sequences,
                             accept_CF = TRUE,
                             min_seq_length = 8,
                             global_vgene = FALSE,
                             vgene_stratify = FALSE,
                             verbose = TRUE) {

    vgene.info <- FALSE
    patient.info <- FALSE
    hla.info <- FALSE
    count.info <- FALSE

    if (is.atomic(cdr3_sequences)) {
        sequences <- data.frame(CDR3b = cdr3_sequences, stringsAsFactors = FALSE)
        if (global_vgene) {
            stop("V-gene information required for global_vgene. ",
                 "Provide a data frame with a 'TRBV' column.", call. = FALSE)
        }
        if (vgene_stratify) {
            stop("V-gene information required for vgene_stratify. ",
                 "Provide a data frame with a 'TRBV' column.", call. = FALSE)
        }
    } else if (is.data.frame(cdr3_sequences)) {
        if (ncol(cdr3_sequences) == 1) {
            sequences <- data.frame(CDR3b = cdr3_sequences[, 1],
                                    stringsAsFactors = FALSE)
            if (verbose) message("Using first column as CDR3b sequences.")
        } else {
            if (!("CDR3b" %in% colnames(cdr3_sequences))) {
                stop("cdr3_sequences must contain a column named 'CDR3b'.",
                     call. = FALSE)
            }
            sequences <- data.frame(CDR3b = cdr3_sequences$CDR3b,
                                    stringsAsFactors = FALSE)
        }

        if ("TRBV" %in% colnames(cdr3_sequences)) {
            vgene.info <- TRUE
            sequences$TRBV <- cdr3_sequences$TRBV
        } else {
            if (global_vgene) {
                stop("V-gene information required for global_vgene.",
                     call. = FALSE)
            }
            if (vgene_stratify) {
                stop("V-gene information required for vgene_stratify.",
                     call. = FALSE)
            }
        }

        if ("patient" %in% colnames(cdr3_sequences)) {
            patient.info <- TRUE
            sequences$patient <- cdr3_sequences$patient
        }
        if ("HLA" %in% colnames(cdr3_sequences)) {
            hla.info <- TRUE
            sequences$HLA <- cdr3_sequences$HLA
        }
        if ("counts" %in% colnames(cdr3_sequences)) {
            count.info <- TRUE
            cdr3_sequences$counts <- as.numeric(cdr3_sequences$counts)
            cdr3_sequences$counts[is.na(cdr3_sequences$counts)] <- 1
            sequences$counts <- cdr3_sequences$counts
        }

        # Carry forward any additional columns
        extra_cols <- setdiff(colnames(cdr3_sequences), colnames(sequences))
        if (length(extra_cols) > 0) {
            sequences <- cbind(sequences, cdr3_sequences[, extra_cols,
                                                          drop = FALSE])
        }

        sequences[] <- lapply(sequences, as.character)
    } else {
        stop("cdr3_sequences must be a character vector or data frame.",
             call. = FALSE)
    }

    # Filter to valid amino acid sequences
    valid_aa <- grep("^[ACDEFGHIKLMNPQRSTVWY]*$", sequences$CDR3b)
    if (length(valid_aa) == 0) {
        stop("No valid CDR3b amino acid sequences found.", call. = FALSE)
    }
    sequences <- sequences[valid_aa, , drop = FALSE]

    # Add sequence IDs
    sequences <- cbind(data.frame(seq_ID = seq_len(nrow(sequences))),
                       sequences)

    list(
        sequences    = sequences,
        vgene.info   = vgene.info,
        patient.info = patient.info,
        hla.info     = hla.info,
        count.info   = count.info
    )
}
