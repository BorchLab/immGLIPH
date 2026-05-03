# immGLIPH 0.99.4

* Added a Validation section to the README with concordance metrics
  against the published GLIPH and GLIPH2 cluster vectors from
  Glanville et al. (2017) and Huang et al. (2020). With paper-matched
  parameters, immGLIPH reproduces the published cluster vectors at
  ARI 0.985 (Glanville) and 0.863 (Huang) on the intersection of
  shared CDR3s. Full benchmark code lives at
  [BorchLab/immGLIPH-benchmark](https://github.com/BorchLab/immGLIPH-benchmark).
* Fixed `clusterScoring()` failure in the clonal-expansion-enrichment
  test when a cluster contains more members than the reference pool
  has rows. The null draw now uses `replace = TRUE` (bootstrap),
  matching the V-gene null and the statistically appropriate choice
  for resampling. Surfaced on the Huang 2020 benchmark.

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
