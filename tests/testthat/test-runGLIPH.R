test_that("runGLIPH rejects invalid method", {
  expect_error(runGLIPH(c("CASSLAPGATNEKLFF"), method = "invalid"),
               "'arg' should be one of")
})

test_that("runGLIPH rejects empty input", {
  expect_error(runGLIPH(character(0)), "No valid")
})

test_that("deprecated aliases emit warnings", {
  expect_warning(
    tryCatch(
      turbo_gliph(character(0)),
      error = function(e) NULL
    ),
    "deprecated"
  )
  expect_warning(
    tryCatch(
      gliph2(character(0)),
      error = function(e) NULL
    ),
    "deprecated"
  )
})
