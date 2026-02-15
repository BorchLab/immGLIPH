#' Example TCR input data
#'
#' A \code{data.frame} of 365 TRB CDR3 sequences with V-gene and patient
#' annotations, derived from the \pkg{scRepertoire} example dataset (Yost
#' et al. 2021). The data were extracted from the \code{\link{gliph_sce}}
#' SingleCellExperiment object using \code{immApex::getIR()}.
#'
#' \describe{
#'   \item{\code{CDR3b}}{Amino acid sequence of the TRB CDR3 region.}
#'   \item{\code{TRBV}}{TRBV gene name (e.g. \code{"TRBV9"}).}
#'   \item{\code{patient}}{Patient/sample identifier
#'     (e.g. \code{"P17B"}, \code{"P19L"}).}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name gliph_input_data
#' @usage data(gliph_input_data)
#' @seealso \code{\link{gliph_sce}} for the parent SingleCellExperiment object,
#'   \code{\link{runGLIPH}} for the main analysis function.
#' @source Yost, K. E. et al. Clonal replacement of tumor-specific T cells
#'   following PD-1 blockade. \emph{Nature Medicine} 25, 1251--1259 (2019).
#'
#'   Built from \pkg{scRepertoire} example data; see
#'   \code{data-raw/build_example_data.R}.
"gliph_input_data"

#' Example SingleCellExperiment with TCR clonal information
#'
#' A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object containing
#' 2,000 genes across 500 cells, with T-cell receptor clonotype information
#' stored in the \code{colData}. Built from the \pkg{scRepertoire} example
#' dataset using \code{combineTCR()} and \code{combineExpression()}.
#'
#' The \code{colData} includes scRepertoire columns such as \code{CTaa}
#' (amino acid clonotype), \code{CTgene} (gene-level clonotype), \code{CTnt}
#' (nucleotide clonotype), \code{CTstrict} (strict clonotype), and clone
#' frequency/proportion columns. These can be parsed by
#' \code{immApex::getIR()} to extract chain-specific TCR data.
#'
#' This object demonstrates how to pass a SingleCellExperiment directly to
#' \code{\link{runGLIPH}}.
#'
#' @docType data
#' @keywords datasets
#' @name gliph_sce
#' @usage data(gliph_sce)
#' @seealso \code{\link{gliph_input_data}} for a plain data.frame extracted
#'   from this object, \code{\link{runGLIPH}} for the main analysis function.
#' @source Yost, K. E. et al. Clonal replacement of tumor-specific T cells
#'   following PD-1 blockade. \emph{Nature Medicine} 25, 1251--1259 (2019).
#'
#'   Built from \pkg{scRepertoire} example data; see
#'   \code{data-raw/build_example_data.R}.
"gliph_sce"

#' Cluster size probabilities in naive reference repertoire
#'
#' A \code{list} with two elements providing expected cluster-size
#' probabilities under the null model (no true convergence):
#' \describe{
#'   \item{\code{original}}{Probabilities from the original GLIPH algorithm,
#'     applied uniformly across all sample sizes.}
#'   \item{\code{simulated}}{Probabilities estimated from 500-step simulations
#'     at sample sizes of 125, 250, 500, 1000, 2000, 4000, 6000, 8000, and
#'     10000 random reference sequences. During scoring the row closest to the
#'     actual sample size is used.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ref_cluster_sizes
#' @usage data(ref_cluster_sizes)
#' @source Glanville, J. et al. Identifying specificity groups in the T cell
#'   receptor repertoire. \emph{Nature} 547, 94--98 (2017).
"ref_cluster_sizes"

#' Germline TCR-beta CDR3 fragments
#'
#' A \code{list} of three \code{data.frame}s containing germline-encoded
#' fragments of V (\code{gTRV}), D (\code{gTRD}), and J (\code{gTRJ}) gene
#' segments that may appear in the CDR3 region. These fragments are used by the
#' GLIPH2 algorithm to identify germline-encoded sequence segments.
#'
#' @docType data
#' @keywords datasets
#' @name gTRB
#' @usage data(gTRB)
#' @source Lefranc, M.-P. IMGT, the international ImMunoGeneTics database.
#'   \emph{Nucl. Acids Res.} 29, 207--209 (2001).
"gTRB"


#' GLIPH reference repertoire list (external data)
#'
#' A named \code{list} of naive TCR repertoire reference databases used for
#' motif enrichment testing and cluster scoring. The data is \strong{not}
#' bundled with the package; it is downloaded on first use from GitHub releases
#' and cached locally via \pkg{BiocFileCache} (see
#' \code{\link{getGLIPHreference}}).
#'
#' Each element is itself a list with three components:
#' \describe{
#'   \item{\code{refseqs}}{A \code{data.frame} with columns \code{CDR3b}
#'     (amino acid sequence) and \code{TRBV} (V-gene name).}
#'   \item{\code{vgene_frequencies}}{A \code{data.frame} with columns
#'     \code{vgene} and \code{freq} giving the relative frequency of each
#'     V gene in the reference repertoire.}
#'   \item{\code{cdr3_length_frequencies}}{A \code{data.frame} with columns
#'     \code{len} and \code{freq} giving the relative frequency of each
#'     CDR3 length in the reference repertoire.}
#' }
#'
#' The following named entries are available:
#' \itemize{
#'   \item \code{"human_v1.0_CD4"}, \code{"human_v1.0_CD8"},
#'         \code{"human_v1.0_CD48"} -- Glanville et al. (2017)
#'   \item \code{"human_v2.0_CD4"}, \code{"human_v2.0_CD8"},
#'         \code{"human_v2.0_CD48"} -- Huang et al. (2020)
#'   \item \code{"mouse_v1.0_CD4"}, \code{"mouse_v1.0_CD8"},
#'         \code{"mouse_v1.0_CD48"} -- Glanville et al. (2017)
#'   \item \code{"gliph_reference"} -- Legacy alias for
#'         \code{"human_v1.0_CD48"}
#' }
#'
#' @seealso \code{\link{getGLIPHreference}} to download or load the data,
#'   \code{\link{runGLIPH}} and \code{\link{clusterScoring}} which use the
#'   reference internally via the \code{refdb_beta} parameter.
#'
#' @source
#' Glanville, J. et al. Identifying specificity groups in the T cell receptor
#'   repertoire. \emph{Nature} 547, 94--98 (2017).
#'
#' Huang, H. et al. Analyzing the \emph{Mycobacterium tuberculosis} immune
#'   response by T-cell receptor clustering with GLIPH2 and genome-wide antigen
#'   screening. \emph{Nature Biotechnology} 38, 1194--1202 (2020).
#'
#' Raw data downloaded from \url{http://50.255.35.37:8080/tools}.
#'
#' @name reference_list
#' @keywords datasets
NULL
