# Tests for deNovoTCRs()

test_that("deNovoTCRs rejects non-character convergence_group_tag", {
  expect_error(deNovoTCRs(convergence_group_tag = 123), "character")
})

test_that("deNovoTCRs rejects multiple tags", {
  expect_error(deNovoTCRs(convergence_group_tag = c("tag1", "tag2")),
               "single character")
})

test_that("deNovoTCRs rejects non-character result_folder", {
  expect_error(deNovoTCRs(convergence_group_tag = "CRG-1",
                           result_folder = 123),
               "character")
})

test_that("deNovoTCRs requires clustering_output when result_folder is empty", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1", result_folder = ""),
    "clustering_output"
  )
})

test_that("deNovoTCRs rejects non-numeric sims", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               sims = "abc"),
    "numeric"
  )
})

test_that("deNovoTCRs rejects non-numeric num_tops", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               num_tops = "abc"),
    "numeric"
  )
})

test_that("deNovoTCRs rejects non-logical make_figure", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               make_figure = "yes"),
    "logical"
  )
})

test_that("deNovoTCRs rejects non-logical normalization", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               normalization = "yes"),
    "logical"
  )
})

test_that("deNovoTCRs rejects tag not in cluster_list", {
  mock_output <- list(
    cluster_list = list(
      "CRG-CASSLAPGATNEKLFF" = data.frame(
        CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
        TRBV = c("TRBV5-1", "TRBV6-2"),
        stringsAsFactors = FALSE
      )
    )
  )

  expect_error(
    deNovoTCRs(convergence_group_tag = "NONEXISTENT_TAG",
               clustering_output = mock_output,
               sims = 10, num_tops = 5, n_cores = 1),
    "Could not find"
  )
})

test_that("deNovoTCRs rejects sims less than 1", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               sims = 0),
    "at least 1"
  )
})

test_that("deNovoTCRs rejects num_tops less than 1", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               num_tops = 0),
    "at least 1"
  )
})

test_that("deNovoTCRs rejects min_length less than 1", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               min_length = 0),
    "at least 1"
  )
})

# ---- Functional tests -------------------------------------------------------

test_that("deNovoTCRs generates sequences from a pre-computed cluster", {
  skip_on_cran()

  # Build a mock clustering output with a cluster of similar sequences
  cluster_members <- data.frame(
    seq_ID  = 1:5,
    CDR3b   = c("CASSLAPGATNEKLFF", "CASSLAPRATNEKLFF",
                "CASSLAPGETQEKLFF", "CASSLAPQATNEKLFF",
                "CASSLAPGAGNEKLFF"),
    TRBV    = rep("TRBV5-1", 5),
    patient = rep("P1", 5),
    stringsAsFactors = FALSE
  )

  mock_output <- list(
    cluster_list = list("CRG-CASSLAPGATNEKLFF" = cluster_members),
    cluster_properties = data.frame(
      tag = "CRG-CASSLAPGATNEKLFF",
      cluster_size = 5,
      stringsAsFactors = FALSE
    ),
    parameters = list(method = "gliph1")
  )

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  result <- deNovoTCRs(
    convergence_group_tag = "CRG-CASSLAPGATNEKLFF",
    clustering_output = mock_output,
    refdb_beta = ref_df,
    sims = 100,
    num_tops = 10,
    min_length = 10,
    make_figure = FALSE,
    n_cores = 1
  )

  expect_type(result, "list")

  # Check output structure
  expect_true("de_novo_sequences" %in% names(result))
  expect_true("sample_sequences_scores" %in% names(result))
  expect_true("cdr3_length_probability" %in% names(result))
  expect_true("PWM_Scoring" %in% names(result))
  expect_true("PWM_Prediction" %in% names(result))

  # de_novo_sequences should be a data.frame
  expect_s3_class(result$de_novo_sequences, "data.frame")

  # PWM_Scoring should have amino acid-related columns
  expect_s3_class(result$PWM_Scoring, "data.frame")
  expect_true(ncol(result$PWM_Scoring) > 0)

  # PWM_Prediction should be a list of data.frames
  expect_type(result$PWM_Prediction, "list")

  # cdr3_length_probability should have probabilities summing to ~1
  expect_s3_class(result$cdr3_length_probability, "data.frame")
  probs <- result$cdr3_length_probability[, ncol(result$cdr3_length_probability)]
  expect_true(abs(sum(as.numeric(probs)) - 1) < 0.01)
})
