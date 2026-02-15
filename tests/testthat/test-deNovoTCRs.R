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
