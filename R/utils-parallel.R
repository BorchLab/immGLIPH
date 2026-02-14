#' Setup parallel backend
#'
#' @param n_cores Number of cores. NULL auto-detects.
#' @return The actual number of cores being used
#' @keywords internal
.setup_parallel <- function(n_cores) {
    if (is.null(n_cores)) {
        n_cores <- max(1L, parallel::detectCores() - 2L)
    }
    n_cores <- max(1L, min(n_cores, parallel::detectCores() - 1L))
    doParallel::registerDoParallel(n_cores)
    n_cores
}

#' Stop parallel backend
#' @keywords internal
.stop_parallel <- function() {
    doParallel::stopImplicitCluster()
}
