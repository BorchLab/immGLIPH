# Tests for getGLIPHreference()

test_that("getGLIPHreference requires BiocFileCache", {
  skip_if(requireNamespace("BiocFileCache", quietly = TRUE),
          "BiocFileCache is installed; skipping missing-dependency test")
  expect_error(getGLIPHreference(verbose = FALSE), "BiocFileCache")
})

test_that("getGLIPHreference returns a list when BiocFileCache is available", {
  skip_if_not_installed("BiocFileCache")
  skip_on_cran()

  result <- getGLIPHreference(verbose = FALSE)
  expect_type(result, "list")
  expect_true(length(result) > 0)

  # All valid reference names should be present
  valid_names <- immGLIPH:::.valid_reference_names()
  for (nm in valid_names) {
    expect_true(nm %in% names(result),
                info = paste("Missing reference:", nm))
  }
})
