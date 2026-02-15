test_that(".validate_params returns all parameters with defaults", {
  result <- immGLIPH:::.validate_params()
  expect_type(result, "list")
  expect_equal(result$refdb_beta, "human_v2.0_CD48")
  expect_equal(result$sim_depth, 1000)
  expect_equal(result$lcminp, 0.01)
  expect_equal(result$motif_length, c(2, 3, 4))
  expect_equal(result$n_cores, 1)
})

test_that(".validate_params rejects bad refdb_beta", {
  expect_error(immGLIPH:::.validate_params(refdb_beta = 42),
               "refdb_beta must be")
  expect_error(immGLIPH:::.validate_params(refdb_beta = "invalid_name"),
               "refdb_beta must be")
})

test_that(".validate_params accepts all built-in reference names", {
  for (nm in immGLIPH:::.valid_reference_names()) {
    result <- immGLIPH:::.validate_params(refdb_beta = nm)
    expect_equal(result$refdb_beta, nm)
  }
})

test_that(".validate_params accepts data frame refdb_beta", {
  df <- data.frame(CDR3b = c("CASSLAPGATNEKLFF"), stringsAsFactors = FALSE)
  result <- immGLIPH:::.validate_params(refdb_beta = df)
  expect_true(is.data.frame(result$refdb_beta))
})

test_that(".validate_params rejects bad sim_depth", {
  expect_error(immGLIPH:::.validate_params(sim_depth = -1), "sim_depth")
  expect_error(immGLIPH:::.validate_params(sim_depth = "abc"), "sim_depth")
})

test_that(".validate_params rounds numeric parameters", {
  result <- immGLIPH:::.validate_params(sim_depth = 100.7, kmer_mindepth = 2.3)
  expect_equal(result$sim_depth, 101)
  expect_equal(result$kmer_mindepth, 2)
})

test_that(".validate_params rejects bad lcminp", {
  expect_error(immGLIPH:::.validate_params(lcminp = 0), "lcminp")
  expect_error(immGLIPH:::.validate_params(lcminp = -1), "lcminp")
})

test_that(".validate_params rejects mismatched lcminove and motif_length", {
  expect_error(
    immGLIPH:::.validate_params(lcminove = c(100, 10), motif_length = c(2, 3, 4)),
    "lcminove"
  )
})

test_that(".validate_params rejects both similarities FALSE", {
  expect_error(
    immGLIPH:::.validate_params(local_similarities = FALSE,
                                global_similarities = FALSE),
    "At least one"
  )
})

test_that(".validate_params adjusts min_seq_length for structboundaries", {
  result <- immGLIPH:::.validate_params(structboundaries = TRUE,
                                         boundary_size = 3,
                                         min_seq_length = 1)
  expect_equal(result$min_seq_length, 7)  # 2*3 + 1
})

test_that(".validate_params rejects bad hla_cutoff", {
  expect_error(immGLIPH:::.validate_params(hla_cutoff = 2), "hla_cutoff")
  expect_error(immGLIPH:::.validate_params(hla_cutoff = -0.1), "hla_cutoff")
})

test_that(".validate_params validates logical parameters", {
  expect_error(immGLIPH:::.validate_params(accept_CF = "yes"),
               "must be logical")
  expect_error(immGLIPH:::.validate_params(structboundaries = 1),
               "must be logical")
  expect_error(immGLIPH:::.validate_params(discontinuous_motifs = "no"),
               "must be logical")
})
