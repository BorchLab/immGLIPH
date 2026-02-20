#' Coerce character columns to numeric where possible
#'
#' @param df A data frame
#' @return Data frame with numeric columns where appropriate
#' @keywords internal
.coerce_numeric_cols <- function(df) {
    if (!is.data.frame(df)) return(df)
    for (i in seq_len(ncol(df))) {
        vals <- suppressWarnings(as.numeric(df[, i]))
        if (!any(is.na(vals))) {
            df[, i] <- vals
        }
    }
    df
}

#' Check for existing output files
#'
#' @param result_folder Path to output folder
#' @param filenames Character vector of filenames to check
#' @return Logical indicating whether saving should proceed
#' @keywords internal
.check_existing_files <- function(result_folder, filenames) {
    for (fname in filenames) {
        fpath <- paste0(result_folder, fname)
        if (file.exists(fpath)) {
            warning("File already exists: ", fpath,
                    ". Saving skipped.", call. = FALSE)
            return(FALSE)
        }
    }
    TRUE
}

#' Prepare result folder path
#'
#' @param result_folder Path string
#' @return Normalized path or "" if no saving
#' @keywords internal
.prepare_result_folder <- function(result_folder) {
    if (!is.character(result_folder)) {
        stop("result_folder must be a character string.", call. = FALSE)
    }
    if (length(result_folder) > 1) {
        stop("result_folder must be a single path.", call. = FALSE)
    }
    if (result_folder == "") return("")

    if (substr(result_folder, nchar(result_folder), nchar(result_folder)) != "/") {
        result_folder <- paste0(result_folder, "/")
    }
    if (!dir.exists(result_folder)) dir.create(result_folder, recursive = TRUE)
    result_folder
}

#' Save parameter list to file
#'
#' @param parameters Named list of parameters
#' @param result_folder Path to output folder
#' @return NULL (invisibly). Called for side effect of writing file.
#' @keywords internal
.save_parameters <- function(parameters, result_folder) {
    paras <- data.frame(
        name  = names(parameters),
        value = vapply(parameters, function(x) paste0(x, collapse = ","),
                       character(1)),
        stringsAsFactors = FALSE
    )
    utils::write.table(paras, file = paste0(result_folder, "parameter.txt"),
                        quote = FALSE, sep = "\t", row.names = FALSE,
                        col.names = FALSE)
}
