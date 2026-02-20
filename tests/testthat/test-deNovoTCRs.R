# Tests for deNovoTCRs()

# ---- Input validation -------------------------------------------------------

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

test_that("deNovoTCRs rejects non-logical accept_sequences_with_C_F_start_end", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               accept_sequences_with_C_F_start_end = "yes"),
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

test_that("deNovoTCRs rejects non-numeric n_cores", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               clustering_output = list(cluster_list = list()),
               n_cores = "abc"),
    "numeric"
  )
})

test_that("deNovoTCRs rejects multiple result_folder paths", {
  expect_error(
    deNovoTCRs(convergence_group_tag = "CRG-1",
               result_folder = c("a", "b")),
    "single path"
  )
})

# ---- Functional tests -------------------------------------------------------

test_that("deNovoTCRs generates sequences from a pre-computed cluster", {
  skip_on_cran()

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
  expect_true("de_novo_sequences" %in% names(result))
  expect_true("sample_sequences_scores" %in% names(result))
  expect_true("cdr3_length_probability" %in% names(result))
  expect_true("PWM_Scoring" %in% names(result))
  expect_true("PWM_Prediction" %in% names(result))
  expect_s3_class(result$de_novo_sequences, "data.frame")
  expect_s3_class(result$PWM_Scoring, "data.frame")
  expect_true(ncol(result$PWM_Scoring) > 0)
  expect_type(result$PWM_Prediction, "list")
  expect_s3_class(result$cdr3_length_probability, "data.frame")
  probs <- result$cdr3_length_probability[, ncol(result$cdr3_length_probability)]
  expect_true(abs(sum(as.numeric(probs)) - 1) < 0.01)
})

test_that("deNovoTCRs with normalization produces norm_score column", {
  skip_on_cran()

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
    cluster_list = list("CRG-CASSLAPGATNEKLFF" = cluster_members)
  )

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  result <- deNovoTCRs(
    convergence_group_tag = "CRG-CASSLAPGATNEKLFF",
    clustering_output = mock_output,
    refdb_beta = ref_df,
    normalization = TRUE,
    sims = 50,
    num_tops = 5,
    min_length = 10,
    make_figure = FALSE,
    n_cores = 1
  )

  expect_true("norm_score" %in% colnames(result$de_novo_sequences))
  expect_true("norm_scores" %in% colnames(result$sample_sequences_scores))
})

test_that("deNovoTCRs saves output to result_folder", {
  skip_on_cran()

  cluster_members <- data.frame(
    seq_ID  = 1:5,
    CDR3b   = c("CASSLAPGATNEKLFF", "CASSLAPRATNEKLFF",
                "CASSLAPGETQEKLFF", "CASSLAPQATNEKLFF",
                "CASSLAPGAGNEKLFF"),
    TRBV    = rep("TRBV5-1", 5),
    stringsAsFactors = FALSE
  )

  mock_output <- list(
    cluster_list = list("CRG-CASSLAPGATNEKLFF" = cluster_members)
  )

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  tmp_dir <- file.path(tempdir(), paste0("denovo_save_test_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  dir.create(tmp_dir, showWarnings = FALSE)

  result <- deNovoTCRs(
    convergence_group_tag = "CRG-CASSLAPGATNEKLFF",
    clustering_output = mock_output,
    result_folder = tmp_dir,
    refdb_beta = ref_df,
    sims = 50,
    num_tops = 5,
    min_length = 10,
    make_figure = FALSE,
    n_cores = 1
  )

  expected_file <- file.path(tmp_dir, "CRG-CASSLAPGATNEKLFF_de_novo.txt")
  expect_true(file.exists(expected_file))
})

test_that("deNovoTCRs warns when sequences fall below min_length", {
  skip_on_cran()

  cluster_members <- data.frame(
    seq_ID  = 1:3,
    CDR3b   = c("CASSLAPGATNEKLFF", "CASSLAP", "CASSLAPGAGNEKLFF"),
    TRBV    = rep("TRBV5-1", 3),
    stringsAsFactors = FALSE
  )

  mock_output <- list(
    cluster_list = list("CRG-test" = cluster_members)
  )

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    TRBV = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    deNovoTCRs(
      convergence_group_tag = "CRG-test",
      clustering_output = mock_output,
      refdb_beta = ref_df,
      sims = 50,
      num_tops = 5,
      min_length = 10,
      make_figure = FALSE,
      n_cores = 1
    ),
    "excluded"
  )
})

test_that("deNovoTCRs with accept_sequences_with_C_F_start_end = FALSE", {
  skip_on_cran()

  cluster_members <- data.frame(
    seq_ID  = 1:5,
    CDR3b   = c("CASSLAPGATNEKLFF", "CASSLAPRATNEKLFF",
                "CASSLAPGETQEKLFF", "CASSLAPQATNEKLFF",
                "CASSLAPGAGNEKLFF"),
    TRBV    = rep("TRBV5-1", 5),
    stringsAsFactors = FALSE
  )

  mock_output <- list(
    cluster_list = list("CRG-test" = cluster_members)
  )

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    TRBV = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  result <- deNovoTCRs(
    convergence_group_tag = "CRG-test",
    clustering_output = mock_output,
    refdb_beta = ref_df,
    accept_sequences_with_C_F_start_end = FALSE,
    sims = 50,
    num_tops = 5,
    min_length = 10,
    make_figure = FALSE,
    n_cores = 1
  )

  expect_type(result, "list")
  expect_s3_class(result$de_novo_sequences, "data.frame")
})
