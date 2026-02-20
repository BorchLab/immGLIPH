# Tests for .global_fisher()

# Register sequential backend for %dopar%
foreach::registerDoSEQ()

# ---- Helper: small synthetic data ------------------------------------------

.make_global_fisher_data <- function() {
  # Sequences with similar structures (differ by 1 AA)
  sample_seqs <- c(
    "CASSLAPGATNEKLFF",
    "CASSLTPGATNEKLFF",  # differs at pos 5 (A -> T)
    "CASSLMPGATNEKLFF",  # differs at pos 5 (A -> M)
    "CASSLDRGEVFF",
    "CASSLDRGEVFF",       # duplicate
    "CASSYLAGGRNTLYF"
  )
  sample_seqs_uniq <- unique(sample_seqs)

  sequences <- data.frame(
    seq_ID  = seq_along(sample_seqs),
    CDR3b   = sample_seqs,
    TRBV    = c("TRBV5-1", "TRBV5-1", "TRBV6-2",
                "TRBV5-1", "TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  motif_region <- substr(sample_seqs_uniq, 4,
                          nchar(sample_seqs_uniq) - 3)

  set.seed(42)
  ref_seqs <- vapply(seq_len(100), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1, 3:9, 11:14, 16:20, 23, 25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_seqs <- unique(ref_seqs)
  ref_motif_region <- substr(ref_seqs, 4, nchar(ref_seqs) - 3)

  list(
    sample_seqs      = sample_seqs_uniq,
    sequences        = sequences,
    motif_region     = motif_region,
    ref_seqs         = ref_seqs,
    ref_motif_region = ref_motif_region
  )
}

# ---- Output structure -------------------------------------------------------

test_that(".global_fisher returns list with cluster_list and global_similarities", {
  skip_if_not_installed("immApex")

  d <- .make_global_fisher_data()
  result <- immGLIPH:::.global_fisher(
    seqs                   = d$sample_seqs,
    motif_region           = d$motif_region,
    sequences              = d$sequences,
    refseqs                = d$ref_seqs,
    refseqs_motif_region   = d$ref_motif_region,
    structboundaries       = TRUE,
    boundary_size          = 3,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    no_cores               = 1,
    verbose                = FALSE
  )

  expect_type(result, "list")
  expect_true("cluster_list" %in% names(result))
  expect_true("global_similarities" %in% names(result))
  expect_type(result$global_similarities, "logical")
})

test_that(".global_fisher cluster_list has expected columns when similarities found", {
  skip_if_not_installed("immApex")

  d <- .make_global_fisher_data()
  result <- immGLIPH:::.global_fisher(
    seqs                   = d$sample_seqs,
    motif_region           = d$motif_region,
    sequences              = d$sequences,
    refseqs                = d$ref_seqs,
    refseqs_motif_region   = d$ref_motif_region,
    structboundaries       = TRUE,
    boundary_size          = 3,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    no_cores               = 1,
    verbose                = FALSE
  )

  expect_s3_class(result$cluster_list, "data.frame")
  if (nrow(result$cluster_list) > 0) {
    expected_cols <- c("cluster_tag", "cluster_size", "unique_CDR3b",
                       "num_in_ref", "fisher.score", "aa_at_position",
                       "TRBV", "CDR3b")
    expect_true(all(expected_cols %in% colnames(result$cluster_list)))
    expect_true(is.numeric(result$cluster_list$fisher.score))
  }
})

# ---- No similarities --------------------------------------------------------

test_that(".global_fisher returns FALSE for no structural matches", {
  skip_if_not_installed("immApex")

  # All very different sequences - no two share a struct
  dissimilar_seqs <- c("CASSLAPGATNEKLFF", "CXXXXXXXXXXXXXXF",
                        "CYYYYYYYYYYYYYYYF")
  motif_region <- substr(dissimilar_seqs, 4, nchar(dissimilar_seqs) - 3)

  sequences <- data.frame(
    seq_ID = seq_along(dissimilar_seqs),
    CDR3b  = dissimilar_seqs,
    TRBV   = c("TRBV5-1", "TRBV6-2", "TRBV7-2"),
    stringsAsFactors = FALSE
  )

  ref_seqs <- c("CASSQQQQQQQQQQF", "CASSWWWWWWWWWWF")

  result <- immGLIPH:::.global_fisher(
    seqs                   = dissimilar_seqs,
    motif_region           = motif_region,
    sequences              = sequences,
    refseqs                = ref_seqs,
    refseqs_motif_region   = substr(ref_seqs, 4, nchar(ref_seqs) - 3),
    structboundaries       = TRUE,
    boundary_size          = 3,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    no_cores               = 1,
    verbose                = FALSE
  )

  expect_false(result$global_similarities)
  expect_equal(nrow(result$cluster_list), 0)
})

# ---- V-gene restriction -----------------------------------------------------

test_that(".global_fisher with global_vgene restricts to same V-gene", {
  skip_if_not_installed("immApex")

  # Two similar seqs but different V-genes
  similar_seqs <- c("CASSLAPGATNEKLFF", "CASSLTPGATNEKLFF")
  motif_region <- substr(similar_seqs, 4, nchar(similar_seqs) - 3)

  sequences <- data.frame(
    seq_ID = 1:2,
    CDR3b  = similar_seqs,
    TRBV   = c("TRBV5-1", "TRBV6-2"),  # different V-genes
    stringsAsFactors = FALSE
  )

  ref_seqs <- c("CASSQQQQQQQQQQF", "CASSWWWWWWWWWWF")

  result_vgene <- immGLIPH:::.global_fisher(
    seqs                   = similar_seqs,
    motif_region           = motif_region,
    sequences              = sequences,
    refseqs                = ref_seqs,
    refseqs_motif_region   = substr(ref_seqs, 4, nchar(ref_seqs) - 3),
    structboundaries       = TRUE,
    boundary_size          = 3,
    global_vgene           = TRUE,
    all_aa_interchangeable = TRUE,
    no_cores               = 1,
    verbose                = FALSE
  )

  result_no_vgene <- immGLIPH:::.global_fisher(
    seqs                   = similar_seqs,
    motif_region           = motif_region,
    sequences              = sequences,
    refseqs                = ref_seqs,
    refseqs_motif_region   = substr(ref_seqs, 4, nchar(ref_seqs) - 3),
    structboundaries       = TRUE,
    boundary_size          = 3,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    no_cores               = 1,
    verbose                = FALSE
  )

  # V-gene restriction should yield fewer or equal clusters
  n_vgene    <- if (result_vgene$global_similarities)
    nrow(result_vgene$cluster_list) else 0
  n_no_vgene <- if (result_no_vgene$global_similarities)
    nrow(result_no_vgene$cluster_list) else 0
  expect_true(n_vgene <= n_no_vgene)
})
