# Tests for loadGLIPH()

# ---- Input validation -------------------------------------------------------

test_that("loadGLIPH rejects non-character result_folder", {
  expect_error(loadGLIPH(result_folder = 123), "character")
})

test_that("loadGLIPH rejects multiple paths", {
  expect_error(loadGLIPH(result_folder = c("path1", "path2")), "single path")
})

test_that("loadGLIPH rejects empty string", {
  expect_error(loadGLIPH(result_folder = ""), "path must be specified")
})

test_that("loadGLIPH rejects non-existent path", {
  expect_error(loadGLIPH(result_folder = "/tmp/nonexistent_gliph_path_xyz"),
               "does not exist")
})

# ---- Functional tests: load saved GLIPH1 results ----------------------------

test_that("loadGLIPH loads gliph1 output correctly", {
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  utils::data("gliph_input_data", package = "immGLIPH")
  small_data <- gliph_input_data[seq_len(50), ]

  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1, 3:9, 11:14, 16:20, 23, 25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  ref_df <- ref_df[!duplicated(ref_df$CDR3b), ]

  tmp_dir <- file.path(tempdir(), paste0("gliph1_load_test_", Sys.getpid()))
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

  loaded <- loadGLIPH(result_folder = tmp_dir)
  expect_type(loaded, "list")
  expect_true("parameters" %in% names(loaded))
  expect_true("cluster_list" %in% names(loaded))
  expect_true("cluster_properties" %in% names(loaded))
  expect_true("connections" %in% names(loaded))
  expect_true("motif_enrichment" %in% names(loaded))

  # Parameters should be loadable
  expect_true("method" %in% names(loaded$parameters))
  expect_equal(loaded$parameters$method, "gliph1")

  # cluster_list should be a named list
  if (!is.null(loaded$cluster_list)) {
    expect_type(loaded$cluster_list, "list")
    expect_true(length(names(loaded$cluster_list)) > 0)
  }
})

# ---- Functional tests: load saved GLIPH2 results ----------------------------

test_that("loadGLIPH loads gliph2 output correctly", {
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  utils::data("gliph_input_data", package = "immGLIPH")
  small_data <- gliph_input_data[seq_len(50), ]

  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1, 3:9, 11:14, 16:20, 23, 25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  ref_df <- ref_df[!duplicated(ref_df$CDR3b), ]

  tmp_dir <- file.path(tempdir(), paste0("gliph2_load_test_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph2",
    refdb_beta = ref_df,
    result_folder = tmp_dir,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  loaded <- loadGLIPH(result_folder = tmp_dir)
  expect_type(loaded, "list")
  expect_true("parameters" %in% names(loaded))
  expect_true("cluster_list" %in% names(loaded))
  expect_true("cluster_properties" %in% names(loaded))
  expect_true("motif_enrichment" %in% names(loaded))

  expect_equal(loaded$parameters$method, "gliph2")
})

# ---- Path normalization ------------------------------------------------------

test_that("loadGLIPH adds trailing slash if missing", {
  skip_on_cran()
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  utils::data("gliph_input_data", package = "immGLIPH")
  small_data <- gliph_input_data[seq_len(50), ]

  set.seed(42)
  extra_seqs <- vapply(seq_len(200), function(i) {
    paste0("C", paste0(sample(LETTERS[c(1, 3:9, 11:14, 16:20, 23, 25)],
                               sample(8:14, 1), replace = TRUE),
                        collapse = ""), "F")
  }, character(1))
  ref_df <- data.frame(
    CDR3b = extra_seqs,
    TRBV = sample(c("TRBV5-1", "TRBV6-2", "TRBV7-2"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  ref_df <- ref_df[!duplicated(ref_df$CDR3b), ]

  tmp_dir <- file.path(tempdir(), paste0("gliph_slash_test_", Sys.getpid()))
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

  # Load without trailing slash - should still work
  path_no_slash <- sub("/$", "", tmp_dir)
  loaded <- loadGLIPH(result_folder = path_no_slash)
  expect_type(loaded, "list")
})

test_that("loadGLIPH errors when parameter.txt is missing", {
  tmp_dir <- file.path(tempdir(), paste0("gliph_empty_test_", Sys.getpid()))
  dir.create(tmp_dir, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  expect_error(loadGLIPH(result_folder = tmp_dir), "missing")
})
