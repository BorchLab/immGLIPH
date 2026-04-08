# immGLIPH 0.99.3

* Replaced `foreach`/`doParallel` with `BiocParallel` for parallelization
  across all functions per Bioconductor recommendations.
* Replaced `devtools::install_github()` references with
  `BiocManager::install()` in vignette and error messages.
* Updated vignette to use `SingleCellExperiment` instead of `Seurat` for
  the single-cell workflow example.
* Fixed `combineTCR()` example in vignette (removed obsolete `cells`
  argument).
* Made more vignette code chunks evaluable (`clusterScoring()`,
  `deNovoTCRs()`, `plotNetwork()` examples).
* Noted that `scRepertoire` and `immApex` are Bioconductor packages.
* Replaced iterative for-loop list growing with vectorized alternatives
  (`lapply()`, `vapply()`, `Reduce()`).
* Standardized code spacing around operators and after commas per
  Bioconductor coding style.

# immGLIPH 0.99.2

* Fixing roxygen documentation issue creating warnings. 

# immGLIPH 0.99.1

* Moved reference repertoire data hosting from GitHub releases to Zenodo
  (https://zenodo.org/records/18925758) per Bioconductor requirements.
* Updated `getGLIPHreference()` to download from Zenodo.

# immGLIPH 0.99.0

* Initial Bioconductor submission
* R implementation of GLIPH and GLIPH2 algorithms for TCR clustering
* Integration with scRepertoire ecosystem via immApex
* Interactive network visualization with plotNetwork()
* De novo TCR sequence generation with deNovoTCRs()
