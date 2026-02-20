# Tests for runGLIPH()

# ---- Parameter validation (.validate_params) ---------------------------------

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

test_that(".validate_params rejects bad n_cores", {
  expect_error(immGLIPH:::.validate_params(n_cores = 0), "n_cores")
  expect_error(immGLIPH:::.validate_params(n_cores = -1), "n_cores")
  expect_error(immGLIPH:::.validate_params(n_cores = "abc"), "n_cores")
})

test_that(".validate_params accepts NULL n_cores", {
  result <- immGLIPH:::.validate_params(n_cores = NULL)
  expect_null(result$n_cores)
})

test_that(".validate_params rejects bad motif_length", {
  expect_error(immGLIPH:::.validate_params(motif_length = "abc"),
               "motif_length")
  expect_error(immGLIPH:::.validate_params(motif_length = 0),
               "motif_length")
  expect_error(immGLIPH:::.validate_params(motif_length = -1),
               "motif_length")
})

test_that(".validate_params rounds motif_length", {
  # Must also adjust lcminove to match motif_length length
  result <- immGLIPH:::.validate_params(motif_length = c(2.7, 3.2),
                                         lcminove = c(100, 10))
  expect_equal(result$motif_length, c(3, 3))
})

test_that(".validate_params rejects bad cluster_min_size", {
  expect_error(immGLIPH:::.validate_params(cluster_min_size = "abc"),
               "cluster_min_size")
  expect_error(immGLIPH:::.validate_params(cluster_min_size = 0),
               "cluster_min_size")
})

test_that(".validate_params rejects bad boundary_size", {
  expect_error(immGLIPH:::.validate_params(boundary_size = "abc"),
               "boundary_size")
  expect_error(immGLIPH:::.validate_params(boundary_size = -1),
               "boundary_size")
  # 0 is valid (>= 0)
  result <- immGLIPH:::.validate_params(boundary_size = 0)
  expect_equal(result$boundary_size, 0)
})

test_that(".validate_params validates ref_cluster_size", {
  expect_error(immGLIPH:::.validate_params(ref_cluster_size = "invalid"),
               "ref_cluster_size")
  result <- immGLIPH:::.validate_params(ref_cluster_size = "original")
  expect_equal(result$ref_cluster_size, "original")
  result2 <- immGLIPH:::.validate_params(ref_cluster_size = "simulated")
  expect_equal(result2$ref_cluster_size, "simulated")
})

test_that(".validate_params validates gccutoff", {
  expect_error(immGLIPH:::.validate_params(gccutoff = -1), "gccutoff")
  expect_error(immGLIPH:::.validate_params(gccutoff = "abc"), "gccutoff")
  result <- immGLIPH:::.validate_params(gccutoff = NULL)
  expect_null(result$gccutoff)
})

# ---- Input validation --------------------------------------------------------

test_that("runGLIPH rejects invalid method", {
  expect_error(
    runGLIPH(cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF"),
             method = "invalid_method"),
    "should be one of"
  )
})

test_that("runGLIPH rejects empty input", {
  expect_error(
    runGLIPH(cdr3_sequences = character(0)),
    ""
  )
})

# ---- End-to-end integration tests -------------------------------------------

test_that("runGLIPH method presets set correct sub-methods", {
  # Test that method presets configure internal methods correctly
  # We can verify this via the returned parameters list

  utils::data("gliph_input_data", package = "immGLIPH")
  small_data <- gliph_input_data[seq_len(50), ]

  # Create a small custom reference to avoid needing BiocFileCache
  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF", "CASSLGQAYEQYF",
              "CASSFSTCSANYGYTF", "CASSPTGGYEQYF", "CASRLAGGRNTLYF",
              "CASSLSGGRNTLYF", "CASSLTGGATNEKLFF", "CASSLDRYEVFF",
              "CASSYLAGARNTLYF", "CASSLTGAEETQYF", "CASSLGQRETQYF",
              "CASSLAQAYEQYF", "CASSFSACSANYGYTF", "CASSPAGGYEQYF",
              "CASRLATGRNTLYF", "CASSLASGGRNTLYF",
              paste0("CASS", paste0(sample(LETTERS[c(1,3:9,11:14,16:20,23,25)],
                                           8, replace = TRUE),
                                    collapse = ""), "F")),
    TRBV = rep(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), length.out = 21),
    stringsAsFactors = FALSE
  )
  # Pad to have enough reference sequences
  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1,3:9,11:14,16:20,23,25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  extra_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  ref_df <- rbind(ref_df, extra_df)
  ref_df <- ref_df[!duplicated(ref_df$CDR3b), ]

  # Test gliph1 method
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available (dev version required)")

  res_g1 <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  expect_type(res_g1, "list")
  expect_true("parameters" %in% names(res_g1))
  expect_equal(res_g1$parameters$method, "gliph1")
  expect_equal(res_g1$parameters$local_method, "rrs")
  expect_equal(res_g1$parameters$global_method, "cutoff")
  expect_equal(res_g1$parameters$clustering_method, "GLIPH1.0")
  expect_true("motif_enrichment" %in% names(res_g1))
  expect_true("connections" %in% names(res_g1))
  expect_true("cluster_properties" %in% names(res_g1))
  expect_true("cluster_list" %in% names(res_g1))
})

test_that("runGLIPH gliph2 method runs end-to-end with custom reference", {
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available (dev version required)")

  utils::data("gliph_input_data", package = "immGLIPH")
  small_data <- gliph_input_data[seq_len(50), ]

  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1,3:9,11:14,16:20,23,25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  ref_df <- ref_df[!duplicated(ref_df$CDR3b), ]

  res_g2 <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph2",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  expect_type(res_g2, "list")
  expect_equal(res_g2$parameters$method, "gliph2")
  expect_equal(res_g2$parameters$local_method, "fisher")
  expect_equal(res_g2$parameters$global_method, "fisher")
  expect_equal(res_g2$parameters$clustering_method, "GLIPH2.0")
  expect_true("motif_enrichment" %in% names(res_g2))
  expect_true("cluster_list" %in% names(res_g2))
})

test_that("runGLIPH output structure is complete", {
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available (dev version required)")

  utils::data("gliph_input_data", package = "immGLIPH")
  small_data <- gliph_input_data[seq_len(30), ]

  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1,3:9,11:14,16:20,23,25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  # Check all expected output elements
  expected_names <- c("sample_log", "motif_enrichment", "global_enrichment",
                      "connections", "cluster_properties", "cluster_list",
                      "parameters")
  for (nm in expected_names) {
    expect_true(nm %in% names(res), info = paste("Missing element:", nm))
  }

  # motif_enrichment should have sub-elements
  expect_true("selected_motifs" %in% names(res$motif_enrichment))
  expect_true("all_motifs" %in% names(res$motif_enrichment))

  # parameters should record what was used
  expect_true("method" %in% names(res$parameters))
  expect_true("sim_depth" %in% names(res$parameters))
  expect_true("n_cores" %in% names(res$parameters))
})

test_that("runGLIPH handles character vector input", {
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available (dev version required)")

  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
            "CASSLTGGEETQYF", "CASSLGGRETQYF", "CASSLGQAYEQYF",
            "CASSFSTCSANYGYTF", "CASSPTGGYEQYF")

  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1,3:9,11:14,16:20,23,25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )

  res <- runGLIPH(
    cdr3_sequences = seqs,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  expect_type(res, "list")
  expect_true("motif_enrichment" %in% names(res))
})

test_that("runGLIPH can save and load results via result_folder", {
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available (dev version required)")

  utils::data("gliph_input_data", package = "immGLIPH")
  small_data <- gliph_input_data[seq_len(50), ]

  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1,3:9,11:14,16:20,23,25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )

  tmp_dir <- file.path(tempdir(), "gliph_test_output")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    result_folder = tmp_dir,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  # Verify output files were created
  expect_true(dir.exists(tmp_dir))
  expect_true(file.exists(file.path(tmp_dir, "parameter.txt")))

  # Test loadGLIPH can read them back
  loaded <- loadGLIPH(result_folder = tmp_dir)
  expect_type(loaded, "list")
  expect_true("parameters" %in% names(loaded))
  expect_true("cluster_list" %in% names(loaded))
})
