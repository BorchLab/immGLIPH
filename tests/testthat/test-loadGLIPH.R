# Tests for loadGLIPH()

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
