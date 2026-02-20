# Tests for utility functions in utils-output.R

# ---- .coerce_numeric_cols ---------------------------------------------------

test_that(".coerce_numeric_cols converts character columns to numeric", {
  df <- data.frame(
    a = c("1", "2", "3"),
    b = c("hello", "world", "foo"),
    c = c("1.5", "2.5", "3.5"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.coerce_numeric_cols(df)
  expect_true(is.numeric(result$a))
  expect_true(is.character(result$b))
  expect_true(is.numeric(result$c))
})

test_that(".coerce_numeric_cols returns non-data.frame input unchanged", {
  expect_equal(immGLIPH:::.coerce_numeric_cols("hello"), "hello")
  expect_equal(immGLIPH:::.coerce_numeric_cols(42), 42)
  expect_null(immGLIPH:::.coerce_numeric_cols(NULL))
})

test_that(".coerce_numeric_cols handles mixed numeric/character columns", {
  df <- data.frame(
    a = c("1", "abc", "3"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.coerce_numeric_cols(df)
  # Should remain character because "abc" can't be coerced
  expect_true(is.character(result$a))
})

# ---- .check_existing_files --------------------------------------------------

test_that(".check_existing_files returns TRUE when no files exist", {
  tmp_dir <- file.path(tempdir(), paste0("check_test_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  dir.create(tmp_dir, showWarnings = FALSE)

  result <- immGLIPH:::.check_existing_files(
    paste0(tmp_dir, "/"),
    c("file1.txt", "file2.txt")
  )
  expect_true(result)
})

test_that(".check_existing_files returns FALSE and warns when files exist", {
  tmp_dir <- file.path(tempdir(), paste0("check_test2_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  dir.create(tmp_dir, showWarnings = FALSE)

  writeLines("test", file.path(tmp_dir, "existing.txt"))

  expect_warning(
    result <- immGLIPH:::.check_existing_files(
      paste0(tmp_dir, "/"),
      c("existing.txt", "other.txt")
    ),
    "already exists"
  )
  expect_false(result)
})

# ---- .prepare_result_folder -------------------------------------------------

test_that(".prepare_result_folder appends trailing slash", {
  tmp_dir <- file.path(tempdir(), paste0("prep_test_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  result <- immGLIPH:::.prepare_result_folder(tmp_dir)
  expect_equal(substr(result, nchar(result), nchar(result)), "/")
  expect_true(dir.exists(result))
})

test_that(".prepare_result_folder returns empty string for empty input", {
  result <- immGLIPH:::.prepare_result_folder("")
  expect_equal(result, "")
})

test_that(".prepare_result_folder errors on non-character input", {
  expect_error(
    immGLIPH:::.prepare_result_folder(123),
    "character"
  )
})

test_that(".prepare_result_folder errors on multiple paths", {
  expect_error(
    immGLIPH:::.prepare_result_folder(c("a", "b")),
    "single path"
  )
})

test_that(".prepare_result_folder preserves existing trailing slash", {
  tmp_dir <- file.path(tempdir(), paste0("prep_test2_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  result <- immGLIPH:::.prepare_result_folder(paste0(tmp_dir, "/"))
  expect_equal(substr(result, nchar(result), nchar(result)), "/")
})

test_that(".prepare_result_folder creates nested directories", {
  tmp_dir <- file.path(tempdir(), paste0("prep_nest_", Sys.getpid()),
                       "sub1", "sub2")
  on.exit(unlink(file.path(tempdir(), paste0("prep_nest_", Sys.getpid())),
                 recursive = TRUE), add = TRUE)

  result <- immGLIPH:::.prepare_result_folder(tmp_dir)
  expect_true(dir.exists(result))
})

# ---- .save_parameters -------------------------------------------------------

test_that(".save_parameters writes parameter file", {
  tmp_dir <- file.path(tempdir(), paste0("save_params_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  dir.create(tmp_dir, showWarnings = FALSE)

  params <- list(
    method = "gliph2",
    sim_depth = 1000,
    motif_length = c(2, 3, 4)
  )

  immGLIPH:::.save_parameters(params, paste0(tmp_dir, "/"))

  param_file <- file.path(tmp_dir, "parameter.txt")
  expect_true(file.exists(param_file))

  content <- readLines(param_file)
  expect_true(any(grepl("method", content)))
  expect_true(any(grepl("gliph2", content)))
  expect_true(any(grepl("motif_length", content)))
  expect_true(any(grepl("2,3,4", content)))
})
