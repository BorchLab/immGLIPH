# Tests for .setup_parallel()

# ---- .setup_parallel with n_cores = 1 returns SerialParam --------------------

test_that(".setup_parallel with 1 core returns a SerialParam", {
  result <- immGLIPH:::.setup_parallel(1)
  expect_s4_class(result, "SerialParam")
})

# ---- .setup_parallel with NULL returns a valid param -------------------------

test_that(".setup_parallel with NULL returns a valid BiocParallelParam", {
  result <- immGLIPH:::.setup_parallel(NULL)
  expect_true(is(result, "BiocParallelParam"))
})

# ---- .setup_parallel with multiple cores on non-Windows ----------------------

test_that(".setup_parallel returns MulticoreParam on non-Windows with multiple cores", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 3,
          "Need at least 3 cores to test MulticoreParam")
  result <- immGLIPH:::.setup_parallel(2)
  expect_s4_class(result, "MulticoreParam")
})

# ---- .setup_parallel clamps to valid range -----------------------------------

test_that(".setup_parallel clamps excessive core count", {
  result <- immGLIPH:::.setup_parallel(9999)
  expect_true(is(result, "BiocParallelParam"))
})

test_that(".setup_parallel clamps negative core count to SerialParam", {
  result <- immGLIPH:::.setup_parallel(-5)
  expect_s4_class(result, "SerialParam")
})

test_that(".setup_parallel clamps zero to SerialParam", {
  result <- immGLIPH:::.setup_parallel(0)
  expect_s4_class(result, "SerialParam")
})
