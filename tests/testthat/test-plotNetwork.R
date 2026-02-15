# Tests for plotNetwork()

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

test_that("plotNetwork requires clustering output with clusters", {
  mock_output <- list(
    cluster_list = NULL,
    cluster_properties = NULL,
    parameters = list(gliph_version = 1)
  )
  expect_error(plotNetwork(clustering_output = mock_output), "does not contain")
})
