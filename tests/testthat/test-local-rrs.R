# Tests for .local_rrs()

# Register sequential backend for %dopar%
foreach::registerDoSEQ()

# ---- Helper: small synthetic data ------------------------------------------

.make_rrs_data <- function() {
  sample_seqs <- c(
    "CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
    "CASSLTGGEETQYF", "CASSLGGRETQYF", "CASSLGQAYEQYF",
    "CASSFSTCSANYGYTF", "CASSPTGGYEQYF"
  )
  sequences <- data.frame(
    seq_ID = seq_along(sample_seqs),
    CDR3b  = sample_seqs,
    TRBV   = rep(c("TRBV5-1", "TRBV6-2"), length.out = length(sample_seqs)),
    stringsAsFactors = FALSE
  )
  motif_region <- substr(sample_seqs, 4, nchar(sample_seqs) - 3)

  set.seed(42)
  ref_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1, 3:9, 11:14, 16:20, 23, 25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_seqs <- unique(ref_seqs)
  ref_motif_region <- substr(ref_seqs, 4, nchar(ref_seqs) - 3)

  list(
    sample_seqs      = sample_seqs,
    sequences        = sequences,
    motif_region     = motif_region,
    ref_seqs         = ref_seqs,
    ref_motif_region = ref_motif_region
  )
}

# ---- Output structure -------------------------------------------------------

test_that(".local_rrs returns list with sample_log, selected_motifs, all_motifs", {
  d <- .make_rrs_data()
  result <- immGLIPH:::.local_rrs(
    motif_region             = d$motif_region,
    refseqs_motif_region     = d$ref_motif_region,
    seqs                     = d$sample_seqs,
    sequences                = d$sequences,
    motif_length             = 2,
    sim_depth                = 5,
    kmer_mindepth            = 1,
    lcminp                   = 1.0,
    lcminove                 = 0,
    discontinuous_motifs     = FALSE,
    cdr3_len_stratify        = FALSE,
    vgene_stratify           = FALSE,
    no_cores                 = 1,
    verbose                  = FALSE,
    motif_lengths_list       = list(),
    ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = list(),
    ref_motif_vgenes_id_list = list(),
    lengths_vgenes_list      = list(),
    ref_lengths_vgenes_list  = list()
  )

  expect_type(result, "list")
  expect_true("sample_log" %in% names(result))
  expect_true("selected_motifs" %in% names(result))
  expect_true("all_motifs" %in% names(result))
})

# ---- sample_log structure ---------------------------------------------------

test_that(".local_rrs sample_log has sim_depth + 1 rows", {
  d <- .make_rrs_data()
  sim_depth <- 5
  result <- immGLIPH:::.local_rrs(
    motif_region             = d$motif_region,
    refseqs_motif_region     = d$ref_motif_region,
    seqs                     = d$sample_seqs,
    sequences                = d$sequences,
    motif_length             = 2,
    sim_depth                = sim_depth,
    kmer_mindepth            = 1,
    lcminp                   = 1.0,
    lcminove                 = 0,
    discontinuous_motifs     = FALSE,
    cdr3_len_stratify        = FALSE,
    vgene_stratify           = FALSE,
    no_cores                 = 1,
    verbose                  = FALSE,
    motif_lengths_list       = list(),
    ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = list(),
    ref_motif_vgenes_id_list = list(),
    lengths_vgenes_list      = list(),
    ref_lengths_vgenes_list  = list()
  )

  expect_equal(nrow(result$sample_log), sim_depth + 1)
  expect_equal(rownames(result$sample_log)[1], "Discovery")
})

# ---- all_motifs structure ---------------------------------------------------

test_that(".local_rrs all_motifs has expected columns", {
  d <- .make_rrs_data()
  result <- immGLIPH:::.local_rrs(
    motif_region             = d$motif_region,
    refseqs_motif_region     = d$ref_motif_region,
    seqs                     = d$sample_seqs,
    sequences                = d$sequences,
    motif_length             = 2,
    sim_depth                = 5,
    kmer_mindepth            = 1,
    lcminp                   = 1.0,
    lcminove                 = 0,
    discontinuous_motifs     = FALSE,
    cdr3_len_stratify        = FALSE,
    vgene_stratify           = FALSE,
    no_cores                 = 1,
    verbose                  = FALSE,
    motif_lengths_list       = list(),
    ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = list(),
    ref_motif_vgenes_id_list = list(),
    lengths_vgenes_list      = list(),
    ref_lengths_vgenes_list  = list()
  )

  expected_cols <- c("motif", "counts", "num_in_ref", "avgRef",
                     "topRef", "OvE", "p.value")
  expect_true(all(expected_cols %in% colnames(result$all_motifs)))
})

# ---- p-values ---------------------------------------------------------------

test_that(".local_rrs p-values are in (0, 1]", {
  d <- .make_rrs_data()
  result <- immGLIPH:::.local_rrs(
    motif_region             = d$motif_region,
    refseqs_motif_region     = d$ref_motif_region,
    seqs                     = d$sample_seqs,
    sequences                = d$sequences,
    motif_length             = 2,
    sim_depth                = 10,
    kmer_mindepth            = 1,
    lcminp                   = 1.0,
    lcminove                 = 0,
    discontinuous_motifs     = FALSE,
    cdr3_len_stratify        = FALSE,
    vgene_stratify           = FALSE,
    no_cores                 = 1,
    verbose                  = FALSE,
    motif_lengths_list       = list(),
    ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = list(),
    ref_motif_vgenes_id_list = list(),
    lengths_vgenes_list      = list(),
    ref_lengths_vgenes_list  = list()
  )

  pvals <- result$all_motifs$p.value
  expect_true(is.numeric(pvals))
  expect_true(all(pvals > 0 & pvals <= 1))
})

# ---- Filtering --------------------------------------------------------------

test_that(".local_rrs kmer_mindepth filters correctly", {
  d <- .make_rrs_data()

  result <- immGLIPH:::.local_rrs(
    motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
    seqs = d$sample_seqs, sequences = d$sequences,
    motif_length = 2, sim_depth = 5, kmer_mindepth = 5,
    lcminp = 1.0, lcminove = 0, discontinuous_motifs = FALSE,
    cdr3_len_stratify = FALSE, vgene_stratify = FALSE,
    no_cores = 1, verbose = FALSE,
    motif_lengths_list = list(), ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = list(), ref_motif_vgenes_id_list = list(),
    lengths_vgenes_list = list(), ref_lengths_vgenes_list = list()
  )

  if (nrow(result$selected_motifs) > 0) {
    expect_true(all(result$selected_motifs$counts >= 5))
  }
})

test_that(".local_rrs strict lcminp yields fewer selected motifs", {
  d <- .make_rrs_data()

  run_rrs <- function(lcminp) {
    immGLIPH:::.local_rrs(
      motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
      seqs = d$sample_seqs, sequences = d$sequences,
      motif_length = 2, sim_depth = 5, kmer_mindepth = 1,
      lcminp = lcminp, lcminove = 0, discontinuous_motifs = FALSE,
      cdr3_len_stratify = FALSE, vgene_stratify = FALSE,
      no_cores = 1, verbose = FALSE,
      motif_lengths_list = list(), ref_motif_lengths_id_list = list(),
      motif_region_vgenes_list = list(), ref_motif_vgenes_id_list = list(),
      lengths_vgenes_list = list(), ref_lengths_vgenes_list = list()
    )
  }

  res_strict  <- run_rrs(0.001)
  res_lenient <- run_rrs(1.0)

  expect_true(nrow(res_strict$selected_motifs) <=
                nrow(res_lenient$selected_motifs))
})
