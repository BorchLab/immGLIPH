# immGLIPH

An R implementation of the **GLIPH** and **GLIPH2** algorithms for clustering
T cell receptors (TCRs) predicted to bind the same HLA-restricted peptide
antigen.

immGLIPH identifies specificity groups based on local (motif-based) and global
(sequence-based) CDR3 similarities, then clusters them into convergence groups
and scores each group for biological significance.

**Please cite the original publications:**

- **GLIPH:** Glanville, J. et al. *Identifying specificity groups in the T cell
  receptor repertoire.* Nature 547, 94-98 (2017).
  [doi:10.1038/nature22976](https://doi.org/10.1038/nature22976)

- **GLIPH2:** Huang, H. et al. *Analyzing the Mycobacterium tuberculosis immune
  response by T-cell receptor clustering with GLIPH2 and genome-wide antigen
  screening.* Nature Biotechnology 38, 1194-1202 (2020).
  [doi:10.1038/s41587-020-0505-4](https://doi.org/10.1038/s41587-020-0505-4)

## Installation

```r
devtools::install_github("BorchLab/immGLIPH")
```

## Reference Data

immGLIPH uses naive TCR repertoire reference databases for motif enrichment
testing and cluster scoring. The reference data (~19 MB) is **not** bundled
with the package to keep the install size small. Instead, it is downloaded
automatically on first use and cached locally via
[BiocFileCache](https://bioconductor.org/packages/BiocFileCache/).

### Setup

```r
# Install BiocFileCache (one-time)
BiocManager::install("BiocFileCache")

# Pre-download the reference data (optional -- happens automatically on first runGLIPH() call)
library(immGLIPH)
ref <- getGLIPHreference()
```

### Available References

| Name | Species | Subset | Source |
|:-----|:--------|:-------|:-------|
| `"human_v1.0_CD4"` | Human | CD4 | Glanville et al. (2017) |
| `"human_v1.0_CD8"` | Human | CD8 | Glanville et al. (2017) |
| `"human_v1.0_CD48"` | Human | CD4+CD8 | Glanville et al. (2017) |
| `"human_v2.0_CD4"` | Human | CD4 | Huang et al. (2020) |
| `"human_v2.0_CD8"` | Human | CD8 | Huang et al. (2020) |
| `"human_v2.0_CD48"` | Human | CD4+CD8 | Huang et al. (2020) |
| `"mouse_v1.0_CD4"` | Mouse | CD4 | Huang et al. (2020) |
| `"mouse_v1.0_CD8"` | Mouse | CD8 | Huang et al. (2020) |
| `"mouse_v1.0_CD48"` | Mouse | CD4+CD8 | Huang et al. (2020) |
| `"gliph_reference"` | Human | CD4+CD8 | Legacy alias for `human_v1.0_CD48` |

Select a reference with the `refdb_beta` parameter (default: `"human_v2.0_CD48"`):

```r
# Human (default)
res <- runGLIPH(my_data, refdb_beta = "human_v2.0_CD48")

# Mouse
res <- runGLIPH(mouse_data, refdb_beta = "mouse_v1.0_CD48")

# Custom data frame
res <- runGLIPH(my_data, refdb_beta = my_custom_ref)
```

Each reference includes pre-computed V-gene usage and CDR3 length frequency
distributions that are automatically used during cluster scoring.

### Rebuilding Reference Data

The build script downloads the raw reference files from the GLIPH web server,
processes them, and saves the resulting `reference_list.RData`:

```bash
Rscript data-raw/build_reference_list.R
```

## Quick Start

```r
library(immGLIPH)
data("gliph_input_data")

res <- runGLIPH(
  cdr3_sequences = gliph_input_data,
  method         = "gliph2",
  sim_depth      = 500,
  n_cores        = 2
)

# Convergence group scores
head(res$cluster_properties)

# Enriched motifs
head(res$motif_enrichment$selected_motifs)
```

## Integration with scRepertoire

immGLIPH integrates with the
[scRepertoire](https://github.com/BorchLab/scRepertoire) ecosystem through
[immApex](https://github.com/BorchLab/immApex). `runGLIPH()` can directly
accept Seurat objects, SingleCellExperiment objects, or `combineTCR()` output.

## Documentation

See the package vignette for the full tutorial:

```r
vignette("immGLIPH")
```

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/immGLIPH/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/immGLIPH/issues).

#### [Pull Requests](https://github.com/BorchLab/immGLIPH/pulls) are welcome for bug fixes, new features, or enhancements.
