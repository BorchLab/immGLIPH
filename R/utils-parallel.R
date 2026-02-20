#' Setup parallel backend
#'
#' On Unix-like systems (macOS, Linux) \code{doParallel} uses forked
#' processes that inherit the loaded package namespace. On Windows, PSOCK
#' clusters are used instead; these require loading the package on each
#' worker via \code{library()}, which fails during \code{R CMD check}
#' because the package is not yet installed in the standard library path.
#' To avoid this, we fall back to sequential execution (\code{registerDoSEQ})
#' when only one core is requested or when running on Windows.
#'
#' @param n_cores Number of cores. NULL auto-detects.
#' @return The actual number of cores being used
#' @keywords internal
.setup_parallel <- function(n_cores) {
    if (is.null(n_cores)) {
        n_cores <- max(1L, parallel::detectCores() - 2L)
    }
    n_cores <- max(1L, min(n_cores, parallel::detectCores() - 1L))

    if (n_cores <= 1L || .Platform$OS.type == "windows") {
        foreach::registerDoSEQ()
        n_cores <- 1L
    } else {
        doParallel::registerDoParallel(n_cores)
    }
    n_cores
}

#' Stop parallel backend
#' @return NULL (invisibly). Called for side effect.
#' @keywords internal
.stop_parallel <- function() {
    doParallel::stopImplicitCluster()
}
