# Tests for .validate_params() in utils-validate.R

# ---- Valid defaults ---------------------------------------------------------

test_that(".validate_params returns list with all expected elements", {
  result <- immGLIPH:::.validate_params()
  expect_type(result, "list")
  expected_names <- c("refdb_beta", "v_usage_freq", "cdr3_length_freq",
                      "ref_cluster_size", "sim_depth", "lcminp", "lcminove",
                      "kmer_mindepth", "accept_CF", "min_seq_length",
                      "gccutoff", "structboundaries", "boundary_size",
                      "motif_length", "local_similarities",
                      "global_similarities", "cluster_min_size",
                      "hla_cutoff", "n_cores", "motif_distance_cutoff",
                      "discontinuous_motifs", "all_aa_interchangeable",
                      "boost_local_significance")
  expect_true(all(expected_names %in% names(result)))
})

# ---- refdb_beta validation --------------------------------------------------

test_that(".validate_params rejects invalid refdb_beta name", {
  expect_error(
    immGLIPH:::.validate_params(refdb_beta = "invalid_name"),
    "refdb_beta"
  )
})

test_that(".validate_params rejects non-character, non-data.frame refdb_beta", {
  expect_error(
    immGLIPH:::.validate_params(refdb_beta = 42),
    "refdb_beta"
  )
})

test_that(".validate_params accepts valid refdb_beta names", {
  result <- immGLIPH:::.validate_params(refdb_beta = "human_v2.0_CD48")
  expect_equal(result$refdb_beta, "human_v2.0_CD48")
})

test_that(".validate_params accepts data frame refdb_beta", {
  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF"),
    TRBV = c("TRBV5-1"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.validate_params(refdb_beta = ref_df)
  expect_s3_class(result$refdb_beta, "data.frame")
})

# ---- v_usage_freq validation ------------------------------------------------

test_that(".validate_params rejects invalid v_usage_freq", {
  expect_error(
    immGLIPH:::.validate_params(v_usage_freq = "not_a_df"),
    "v_usage_freq"
  )
})

test_that(".validate_params rejects v_usage_freq with non-numeric frequencies", {
  bad_freq <- data.frame(vgene = "TRBV5-1", freq = "abc",
                          stringsAsFactors = FALSE)
  expect_error(
    immGLIPH:::.validate_params(v_usage_freq = bad_freq),
    "numeric"
  )
})

# ---- cdr3_length_freq validation -------------------------------------------

test_that(".validate_params rejects invalid cdr3_length_freq", {
  expect_error(
    immGLIPH:::.validate_params(cdr3_length_freq = "not_a_df"),
    "cdr3_length_freq"
  )
})

test_that(".validate_params rejects cdr3_length_freq with non-numeric frequencies", {
  bad_freq <- data.frame(length = "10", freq = "abc",
                          stringsAsFactors = FALSE)
  expect_error(
    immGLIPH:::.validate_params(cdr3_length_freq = bad_freq),
    "numeric"
  )
})

# ---- ref_cluster_size validation -------------------------------------------

test_that(".validate_params rejects invalid ref_cluster_size", {
  expect_error(
    immGLIPH:::.validate_params(ref_cluster_size = "invalid"),
    "ref_cluster_size"
  )
})

# ---- sim_depth validation --------------------------------------------------

test_that(".validate_params rejects non-numeric sim_depth", {
  expect_error(
    immGLIPH:::.validate_params(sim_depth = "abc"),
    "sim_depth"
  )
})

test_that(".validate_params rejects sim_depth < 1", {
  expect_error(
    immGLIPH:::.validate_params(sim_depth = 0),
    "sim_depth"
  )
})

test_that(".validate_params rounds sim_depth", {
  result <- immGLIPH:::.validate_params(sim_depth = 5.7)
  expect_equal(result$sim_depth, 6)
})

# ---- lcminp validation ----------------------------------------------------

test_that(".validate_params rejects non-numeric lcminp", {
  expect_error(
    immGLIPH:::.validate_params(lcminp = "abc"),
    "lcminp"
  )
})

test_that(".validate_params rejects lcminp <= 0", {
  expect_error(
    immGLIPH:::.validate_params(lcminp = 0),
    "lcminp"
  )
})

# ---- lcminove validation ---------------------------------------------------

test_that(".validate_params rejects non-numeric lcminove", {
  expect_error(
    immGLIPH:::.validate_params(lcminove = "abc"),
    "lcminove"
  )
})

test_that(".validate_params rejects lcminove length mismatch", {
  expect_error(
    immGLIPH:::.validate_params(lcminove = c(10, 20), motif_length = c(2, 3, 4)),
    "lcminove"
  )
})

test_that(".validate_params rejects lcminove values < 1", {
  expect_error(
    immGLIPH:::.validate_params(lcminove = c(0.5, 1, 1)),
    "lcminove"
  )
})

# ---- kmer_mindepth validation -----------------------------------------------

test_that(".validate_params rejects non-numeric kmer_mindepth", {
  expect_error(
    immGLIPH:::.validate_params(kmer_mindepth = "abc"),
    "kmer_mindepth"
  )
})

test_that(".validate_params rejects kmer_mindepth < 1", {
  expect_error(
    immGLIPH:::.validate_params(kmer_mindepth = 0),
    "kmer_mindepth"
  )
})

# ---- accept_CF validation ---------------------------------------------------

test_that(".validate_params rejects non-logical accept_CF", {
  expect_error(
    immGLIPH:::.validate_params(accept_CF = "yes"),
    "logical"
  )
})

# ---- min_seq_length validation -----------------------------------------------

test_that(".validate_params rejects non-numeric min_seq_length", {
  expect_error(
    immGLIPH:::.validate_params(min_seq_length = "abc"),
    "min_seq_length"
  )
})

test_that(".validate_params adjusts min_seq_length with structboundaries", {
  result <- immGLIPH:::.validate_params(
    structboundaries = TRUE,
    boundary_size = 3,
    min_seq_length = 1
  )
  # min_seq_length should be at least 2 * boundary_size + 1 = 7
  expect_true(result$min_seq_length >= 7)
})

# ---- structboundaries / boundary_size validation ----------------------------

test_that(".validate_params rejects non-logical structboundaries", {
  expect_error(
    immGLIPH:::.validate_params(structboundaries = "yes"),
    "logical"
  )
})

test_that(".validate_params rejects negative boundary_size", {
  expect_error(
    immGLIPH:::.validate_params(boundary_size = -1),
    "boundary_size"
  )
})

# ---- motif_length validation ------------------------------------------------

test_that(".validate_params rejects motif_length with values < 1", {
  expect_error(
    immGLIPH:::.validate_params(motif_length = c(0, 2)),
    "motif_length"
  )
})

# ---- similarities validation ------------------------------------------------

test_that(".validate_params rejects non-logical local_similarities", {
  expect_error(
    immGLIPH:::.validate_params(local_similarities = "yes"),
    "logical"
  )
})

test_that(".validate_params rejects non-logical global_similarities", {
  expect_error(
    immGLIPH:::.validate_params(global_similarities = "yes"),
    "logical"
  )
})

test_that(".validate_params rejects both similarities as FALSE", {
  expect_error(
    immGLIPH:::.validate_params(local_similarities = FALSE,
                                 global_similarities = FALSE),
    "At least one"
  )
})

# ---- gccutoff validation ----------------------------------------------------

test_that(".validate_params rejects invalid gccutoff", {
  expect_error(
    immGLIPH:::.validate_params(gccutoff = "abc"),
    "gccutoff"
  )
  expect_error(
    immGLIPH:::.validate_params(gccutoff = -1),
    "gccutoff"
  )
})

test_that(".validate_params accepts NULL gccutoff", {
  result <- immGLIPH:::.validate_params(gccutoff = NULL)
  expect_null(result$gccutoff)
})

# ---- cluster_min_size validation --------------------------------------------

test_that(".validate_params rejects cluster_min_size < 1", {
  expect_error(
    immGLIPH:::.validate_params(cluster_min_size = 0),
    "cluster_min_size"
  )
})

# ---- hla_cutoff validation --------------------------------------------------

test_that(".validate_params rejects hla_cutoff outside [0, 1]", {
  expect_error(
    immGLIPH:::.validate_params(hla_cutoff = 1.5),
    "hla_cutoff"
  )
  expect_error(
    immGLIPH:::.validate_params(hla_cutoff = -0.1),
    "hla_cutoff"
  )
})

# ---- n_cores validation -----------------------------------------------------

test_that(".validate_params rejects n_cores < 1", {
  expect_error(
    immGLIPH:::.validate_params(n_cores = 0),
    "n_cores"
  )
})

test_that(".validate_params accepts NULL n_cores", {
  result <- immGLIPH:::.validate_params(n_cores = NULL)
  expect_null(result$n_cores)
})

# ---- motif_distance_cutoff validation ----------------------------------------

test_that(".validate_params rejects non-numeric motif_distance_cutoff", {
  expect_error(
    immGLIPH:::.validate_params(motif_distance_cutoff = "abc"),
    "motif_distance_cutoff"
  )
})

# ---- discontinuous_motifs validation ----------------------------------------

test_that(".validate_params rejects non-logical discontinuous_motifs", {
  expect_error(
    immGLIPH:::.validate_params(discontinuous_motifs = "yes"),
    "logical"
  )
})

# ---- all_aa_interchangeable validation --------------------------------------

test_that(".validate_params rejects non-logical all_aa_interchangeable", {
  expect_error(
    immGLIPH:::.validate_params(all_aa_interchangeable = "yes"),
    "logical"
  )
})

# ---- boost_local_significance validation ------------------------------------

test_that(".validate_params rejects non-logical boost_local_significance", {
  expect_error(
    immGLIPH:::.validate_params(boost_local_significance = "yes"),
    "logical"
  )
})
