.onAttach <- function(libname, pkgname) {
    msg <- paste0(
        "immGLIPH v", utils::packageVersion("immGLIPH"), "\n",
        "immGLIPH is an R implementation of the GLIPH algorithm.\n",
        "Please cite the original work:\n",
        "  Glanville et al. Nature 547, 94-98 (2017)\n",
        "  Huang et al. Nat Biotechnol 38, 1194-1202 (2020)\n"
    )
    packageStartupMessage(msg)
}
