# Tests for .local_fisher()

# Register sequential backend for %dopar%
foreach::registerDoSEQ()

# ---- Helper: small synthetic data ------------------------------------------

.make_fisher_data <- function() {
  # Sample sequences that share motifs (e.g., "SLA" appears in several)
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

  # Reference sequences (larger, different distribution)
  set.seed(42)
  ref_seqs <- vapply(seq_len(100), function(i) {
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

test_that(".local_fisher returns list with selected_motifs and all_motifs", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 2,
    lcminp                = 1.0,
    lcminove              = c(1, 1),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_type(result, "list")
  expect_true("selected_motifs" %in% names(result))
  expect_true("all_motifs" %in% names(result))
})

test_that(".local_fisher all_motifs has expected columns", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expected_cols <- c("motif", "counts", "num_in_ref", "avgRef",
                     "topRef", "OvE", "p.value")
  expect_true(all(expected_cols %in% colnames(result$all_motifs)))
  expect_s3_class(result$all_motifs, "data.frame")
})

# ---- p-value and fold change ------------------------------------------------

test_that(".local_fisher p-values are numeric and in [0, 1]", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  pvals <- result$all_motifs$p.value
  expect_true(is.numeric(pvals))
  expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that(".local_fisher OvE is numeric", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_true(is.numeric(result$all_motifs$OvE))
})

# ---- Filtering parameters ---------------------------------------------------

test_that(".local_fisher kmer_mindepth filters out rare motifs", {
  d <- .make_fisher_data()

  result_low <- immGLIPH:::.local_fisher(
    motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
    seqs = d$sample_seqs, refseqs = d$ref_seqs, sequences = d$sequences,
    motif_length = 2, kmer_mindepth = 1, lcminp = 1.0, lcminove = 0,
    discontinuous_motifs = FALSE, motif_distance_cutoff = 1,
    no_cores = 1, verbose = FALSE
  )
  result_high <- immGLIPH:::.local_fisher(
    motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
    seqs = d$sample_seqs, refseqs = d$ref_seqs, sequences = d$sequences,
    motif_length = 2, kmer_mindepth = 5, lcminp = 1.0, lcminove = 0,
    discontinuous_motifs = FALSE, motif_distance_cutoff = 1,
    no_cores = 1, verbose = FALSE
  )

  # Higher mindepth should yield fewer or equal selected motifs
  expect_true(nrow(result_high$selected_motifs) <=
                nrow(result_low$selected_motifs))
  # All selected should satisfy mindepth
  if (nrow(result_high$selected_motifs) > 0) {
    expect_true(all(result_high$selected_motifs$counts >= 5))
  }
})

test_that(".local_fisher lcminp filters out high p-value motifs", {
  d <- .make_fisher_data()

  result_strict <- immGLIPH:::.local_fisher(
    motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
    seqs = d$sample_seqs, refseqs = d$ref_seqs, sequences = d$sequences,
    motif_length = 2, kmer_mindepth = 1, lcminp = 1e-10, lcminove = 0,
    discontinuous_motifs = FALSE, motif_distance_cutoff = 1,
    no_cores = 1, verbose = FALSE
  )
  result_lenient <- immGLIPH:::.local_fisher(
    motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
    seqs = d$sample_seqs, refseqs = d$ref_seqs, sequences = d$sequences,
    motif_length = 2, kmer_mindepth = 1, lcminp = 1.0, lcminove = 0,
    discontinuous_motifs = FALSE, motif_distance_cutoff = 1,
    no_cores = 1, verbose = FALSE
  )

  expect_true(nrow(result_strict$selected_motifs) <=
                nrow(result_lenient$selected_motifs))
})

test_that(".local_fisher lcminove filters by fold change per motif length", {
  d <- .make_fisher_data()

  # Very high OvE threshold should eliminate most motifs
  result <- immGLIPH:::.local_fisher(
    motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
    seqs = d$sample_seqs, refseqs = d$ref_seqs, sequences = d$sequences,
    motif_length = c(2, 3), kmer_mindepth = 1, lcminp = 1.0,
    lcminove = c(1e6, 1e6),
    discontinuous_motifs = FALSE, motif_distance_cutoff = 1,
    no_cores = 1, verbose = FALSE
  )

  expect_equal(nrow(result$selected_motifs), 0)
})

# ---- Discontinuous motifs ---------------------------------------------------

test_that(".local_fisher includes discontinuous motifs when enabled", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region = d$motif_region, refseqs_motif_region = d$ref_motif_region,
    seqs = d$sample_seqs, refseqs = d$ref_seqs, sequences = d$sequences,
    motif_length = 2, kmer_mindepth = 1, lcminp = 1.0, lcminove = 0,
    discontinuous_motifs = TRUE, motif_distance_cutoff = 1,
    no_cores = 1, verbose = FALSE
  )

  disc <- result$all_motifs[grep("\\.", result$all_motifs$motif), ]
  expect_true(nrow(disc) > 0)
})

# ---- immApex path -----------------------------------------------------------

test_that(".local_fisher works via immApex when available", {
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(1, 1),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_type(result, "list")
  expect_s3_class(result$all_motifs, "data.frame")
  expect_true(nrow(result$all_motifs) > 0)
  expect_true(all(c("motif", "counts", "num_in_ref", "OvE", "p.value") %in%
                    colnames(result$all_motifs)))
})

test_that(".local_fisher immApex path with discontinuous motifs", {
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = TRUE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  disc <- result$all_motifs[grep("\\.", result$all_motifs$motif), ]
  expect_true(nrow(disc) > 0)
})

# ---- Verbose messaging -------------------------------------------------------

test_that(".local_fisher prints messages when verbose is TRUE", {
  d <- .make_fisher_data()
  expect_message(
    immGLIPH:::.local_fisher(
      motif_region          = d$motif_region,
      refseqs_motif_region  = d$ref_motif_region,
      seqs                  = d$sample_seqs,
      refseqs               = d$ref_seqs,
      sequences             = d$sequences,
      motif_length          = c(2, 3),
      kmer_mindepth         = 1,
      lcminp                = 1.0,
      lcminove              = c(1, 1),
      discontinuous_motifs  = FALSE,
      motif_distance_cutoff = 1,
      no_cores              = 1,
      verbose               = TRUE
    ),
    "motif"
  )
})

# ---- Single lcminove value ---------------------------------------------------

test_that(".local_fisher works with single lcminove value for multiple motif_lengths", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3, 4),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = 1,
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_type(result, "list")
  expect_s3_class(result$all_motifs, "data.frame")
})

# ---- avgRef normalization ----------------------------------------------------

test_that(".local_fisher avgRef is normalized to sample set size", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_true(is.numeric(result$all_motifs$avgRef))
  expect_true(all(result$all_motifs$avgRef >= 0))
})

# ---- topRef column -----------------------------------------------------------

test_that(".local_fisher topRef column is numeric", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_true(is.numeric(result$all_motifs$topRef))
})

# ---- Empty selected_motifs when all below kmer_mindepth ----------------------

test_that(".local_fisher returns empty selected_motifs when all below kmer_mindepth", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 9999,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_equal(nrow(result$selected_motifs), 0)
  # all_motifs should still contain motifs found
  expect_true(nrow(result$all_motifs) > 0)
})

# ---- Multiple motif lengths with matching lcminove vector --------------------

test_that(".local_fisher correctly applies per-length lcminove thresholds", {
  d <- .make_fisher_data()

  # Very high threshold for length 2, very low for length 3
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(1e6, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  if (nrow(result$selected_motifs) > 0) {
    # No 2-mer motifs should pass the high threshold
    two_mers <- result$selected_motifs[nchar(result$selected_motifs$motif) == 2, ]
    expect_equal(nrow(two_mers), 0)
  }
})

# ---- motif length 4 ---------------------------------------------------------

test_that(".local_fisher works with motif_length = 4", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = 4,
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = 0,
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_type(result, "list")
  if (nrow(result$all_motifs) > 0) {
    # All continuous motifs should be 4 characters
    expect_true(all(nchar(result$all_motifs$motif) == 4))
  }
})

# ---- counts column is always >= 1 in all_motifs -----------------------------

test_that(".local_fisher all_motifs counts are positive", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_true(all(result$all_motifs$counts >= 1))
})

# ---- num_in_ref is non-negative in all_motifs --------------------------------

test_that(".local_fisher num_in_ref is non-negative", {
  d <- .make_fisher_data()
  result <- immGLIPH:::.local_fisher(
    motif_region          = d$motif_region,
    refseqs_motif_region  = d$ref_motif_region,
    seqs                  = d$sample_seqs,
    refseqs               = d$ref_seqs,
    sequences             = d$sequences,
    motif_length          = c(2, 3),
    kmer_mindepth         = 1,
    lcminp                = 1.0,
    lcminove              = c(0, 0),
    discontinuous_motifs  = FALSE,
    motif_distance_cutoff = 1,
    no_cores              = 1,
    verbose               = FALSE
  )

  expect_true(all(result$all_motifs$num_in_ref >= 0))
})
