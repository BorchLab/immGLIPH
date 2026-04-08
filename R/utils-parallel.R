#' Create a BiocParallel backend parameter object
#'
#' Returns a \code{\link[BiocParallel]{BiocParallelParam}} suitable for the
#' current platform and requested number of cores. On Unix-like systems
#' (macOS, Linux) with more than one core, a
#' \code{\link[BiocParallel]{MulticoreParam}} is returned. On Windows or
#' when only one core is requested, a
#' \code{\link[BiocParallel]{SerialParam}} is used instead.
#'
#' @param n_cores Number of cores. \code{NULL} auto-detects.
#' @return A \code{BiocParallelParam} object.
#' @keywords internal
.setup_parallel <- function(n_cores) {
    if (is.null(n_cores)) {
        n_cores <- max(1L, parallel::detectCores() - 2L)
    }
    n_cores <- max(1L, min(n_cores, parallel::detectCores() - 1L))

    if (n_cores <= 1L || .Platform$OS.type == "windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(workers = n_cores)
    }
}
