# Tests for plotNetwork()

# ---- Input validation -------------------------------------------------------

test_that("plotNetwork rejects non-character result_folder", {
  expect_error(plotNetwork(result_folder = 123), "character")
})

test_that("plotNetwork rejects multiple result_folder paths", {
  expect_error(plotNetwork(result_folder = c("a", "b")), "single path")
})

test_that("plotNetwork rejects non-character color_info", {
  expect_error(plotNetwork(clustering_output = list(), color_info = 123),
               "character")
})

test_that("plotNetwork rejects non-function color_palette", {
  expect_error(plotNetwork(clustering_output = list(), color_palette = "red"),
               "function")
})

test_that("plotNetwork rejects invalid local_edge_color", {
  expect_error(plotNetwork(clustering_output = list(), local_edge_color = "not_a_color"),
               "color")
})

test_that("plotNetwork rejects invalid global_edge_color", {
  expect_error(plotNetwork(clustering_output = list(), global_edge_color = "not_a_color"),
               "color")
})

test_that("plotNetwork rejects non-logical absolute_size", {
  expect_error(plotNetwork(clustering_output = list(), absolute_size = "yes"),
               "logical")
})

test_that("plotNetwork rejects non-numeric cluster_min_size", {
  expect_error(plotNetwork(clustering_output = list(), cluster_min_size = "abc"),
               "numeric")
})

test_that("plotNetwork rejects cluster_min_size < 1", {
  expect_error(plotNetwork(clustering_output = list(), cluster_min_size = 0),
               "at least 1")
})

test_that("plotNetwork returns NULL for clustering output with no clusters", {
  mock_output <- list(
    cluster_list = NULL,
    cluster_properties = NULL,
    parameters = list(gliph_version = 1)
  )
  expect_message(result <- plotNetwork(clustering_output = mock_output),
                 "No clusters found")
  expect_null(result)
})

test_that("plotNetwork returns NULL for empty cluster_list", {
  mock_output <- list(
    cluster_list = list(),
    cluster_properties = data.frame(),
    parameters = list(clustering_method = "GLIPH2.0")
  )
  expect_message(result <- plotNetwork(clustering_output = mock_output),
                 "No clusters found")
  expect_null(result)
})

test_that("plotNetwork rejects non-character size_info", {
  expect_error(plotNetwork(clustering_output = list(), size_info = 123),
               "character")
})

test_that("plotNetwork rejects non-character show_additional_columns", {
  expect_error(
    plotNetwork(clustering_output = list(), show_additional_columns = 42),
    "character"
  )
})

test_that("plotNetwork rejects multiple size_info values", {
  expect_error(
    plotNetwork(clustering_output = list(), size_info = c("a", "b")),
    "single column"
  )
})

test_that("plotNetwork rejects multiple color_info values", {
  expect_error(
    plotNetwork(clustering_output = list(), color_info = c("a", "b")),
    "single column"
  )
})

test_that("plotNetwork rejects non-numeric n_cores", {
  expect_error(
    plotNetwork(clustering_output = list(), n_cores = "abc"),
    "numeric"
  )
})

test_that("plotNetwork rejects n_cores < 1", {
  expect_error(
    plotNetwork(clustering_output = list(), n_cores = 0),
    "at least 1"
  )
})

# ---- Functional tests with GLIPH1 output ------------------------------------

test_that("plotNetwork produces visNetwork from gliph_version=1 output", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    plot_obj <- plotNetwork(
      clustering_output = res,
      cluster_min_size = 2,
      n_cores = 1
    )
    expect_s3_class(plot_obj, "visNetwork")
  }
})

test_that("plotNetwork works with color_info = 'none'", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    plot_obj <- plotNetwork(
      clustering_output = res,
      color_info = "none",
      cluster_min_size = 2,
      n_cores = 1
    )
    expect_s3_class(plot_obj, "visNetwork")
  }
})

# ---- Functional tests with GLIPH2 output ------------------------------------

test_that("plotNetwork produces visNetwork from gliph_version=2 output", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph2",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    plot_obj <- plotNetwork(
      clustering_output = res,
      cluster_min_size = 2,
      n_cores = 1
    )
    expect_s3_class(plot_obj, "visNetwork")
  }
})

# ---- plotNetwork with size_info and various coloring -------------------------

test_that("plotNetwork handles size_info parameter", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    # Test with size_info = "cluster_size"
    plot_obj <- plotNetwork(
      clustering_output = res,
      size_info = "cluster_size",
      absolute_size = FALSE,
      cluster_min_size = 2,
      n_cores = 1
    )
    expect_s3_class(plot_obj, "visNetwork")

    # Test with absolute_size = TRUE
    plot_obj2 <- plotNetwork(
      clustering_output = res,
      size_info = "cluster_size",
      absolute_size = TRUE,
      cluster_min_size = 2,
      n_cores = 1
    )
    expect_s3_class(plot_obj2, "visNetwork")
  }
})

test_that("plotNetwork errors for non-existent size_info column", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    expect_error(
      plotNetwork(
        clustering_output = res,
        size_info = "nonexistent_column",
        cluster_min_size = 2,
        n_cores = 1
      ),
      "not found"
    )
  }
})

test_that("plotNetwork errors for non-existent color_info column", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    expect_error(
      plotNetwork(
        clustering_output = res,
        color_info = "nonexistent_column",
        cluster_min_size = 2,
        n_cores = 1
      ),
      "not found"
    )
  }
})

test_that("plotNetwork with categorical color_info", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    # Use TRBV as categorical color_info
    if ("TRBV" %in% colnames(res$cluster_list[[1]])) {
      plot_obj <- plotNetwork(
        clustering_output = res,
        color_info = "TRBV",
        cluster_min_size = 2,
        n_cores = 1
      )
      expect_s3_class(plot_obj, "visNetwork")
    }
  }
})

test_that("plotNetwork with show_additional_columns", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties) &&
      any(as.numeric(res$cluster_properties$cluster_size) >= 2)) {
    plot_obj <- plotNetwork(
      clustering_output = res,
      show_additional_columns = c("TRBV"),
      cluster_min_size = 2,
      n_cores = 1
    )
    expect_s3_class(plot_obj, "visNetwork")
  }
})

test_that("plotNetwork errors when no clusters meet cluster_min_size", {
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

  res <- runGLIPH(
    cdr3_sequences = small_data,
    method = "gliph1",
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1,
    verbose = FALSE
  )

  if (!is.null(res$cluster_properties)) {
    expect_error(
      plotNetwork(
        clustering_output = res,
        cluster_min_size = 10000,
        n_cores = 1
      ),
      "does not contain any clusters"
    )
  }
})
