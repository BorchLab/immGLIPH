# Tests for .setup_parallel() and .stop_parallel()

# ---- .setup_parallel with n_cores = 1 ----------------------------------------

test_that(".setup_parallel with 1 core returns 1", {
  result <- immGLIPH:::.setup_parallel(1)
  expect_equal(result, 1L)
})

# ---- .setup_parallel with NULL auto-detects -----------------------------------

test_that(".setup_parallel with NULL returns at least 1", {
  result <- immGLIPH:::.setup_parallel(NULL)
  expect_true(result >= 1L)
  expect_true(result <= parallel::detectCores())
  # Clean up
  immGLIPH:::.stop_parallel()
})

# ---- .setup_parallel clamps to valid range ------------------------------------

test_that(".setup_parallel clamps excessive core count", {
  max_cores <- parallel::detectCores()
  result <- immGLIPH:::.setup_parallel(max_cores + 100)
  expect_true(result <= max_cores)
  expect_true(result >= 1L)
  # Clean up
  immGLIPH:::.stop_parallel()
})

test_that(".setup_parallel clamps negative core count to 1", {
  result <- immGLIPH:::.setup_parallel(-5)
  expect_equal(result, 1L)
})

test_that(".setup_parallel clamps zero to 1", {
  result <- immGLIPH:::.setup_parallel(0)
  expect_equal(result, 1L)
})

# ---- .stop_parallel runs without error ----------------------------------------

test_that(".stop_parallel executes without error", {
  immGLIPH:::.setup_parallel(1)
  expect_no_error(immGLIPH:::.stop_parallel())
})

# ---- Sequential fallback on Windows or 1 core ---------------------------------

test_that(".setup_parallel registers sequential backend for single core", {
  result <- immGLIPH:::.setup_parallel(1)
  expect_equal(result, 1L)
  # After registerDoSEQ, foreach should still work
  res <- foreach::foreach(i = 1:3, .combine = c) %dopar% { i * 2 }
  expect_equal(res, c(2, 4, 6))
})
