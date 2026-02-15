# Tests for .extract_input() and .parse_sequences()

# ---- .extract_input ----------------------------------------------------------

test_that(".extract_input handles character vector", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  result <- immGLIPH:::.extract_input(seqs)
  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result), "CDR3b")
  expect_equal(nrow(result), 3)
})

test_that(".extract_input handles data frame with CDR3b column", {
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    TRBV  = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.extract_input(df)
  expect_s3_class(result, "data.frame")
  expect_true("CDR3b" %in% colnames(result))
  expect_true("TRBV" %in% colnames(result))
})

test_that(".extract_input rejects invalid input types", {
  expect_error(immGLIPH:::.extract_input(42), "must be a character vector")
  expect_error(immGLIPH:::.extract_input(NULL), "must be a character vector")
})

# ---- .standardize_colnames ---------------------------------------------------

test_that(".standardize_colnames maps alternative names", {
  df <- data.frame(
    cdr3     = c("CASSLAPGATNEKLFF"),
    v_gene   = c("TRBV5-1"),
    sample   = c("P1"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.standardize_colnames(df)
  expect_true("CDR3b" %in% colnames(result))
  expect_true("TRBV" %in% colnames(result))
  expect_true("patient" %in% colnames(result))
})

test_that(".standardize_colnames does not overwrite existing canonical names", {
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF"),
    TRBV  = c("TRBV5-1"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.standardize_colnames(df)
  expect_equal(colnames(result), c("CDR3b", "TRBV"))
})

# ---- .parse_sequences --------------------------------------------------------

test_that(".parse_sequences filters amino acid sequences", {
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "INVALID123", "CASSLDRGEVFF"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_equal(nrow(result$sequences), 2)
  expect_true("seq_ID" %in% colnames(result$sequences))
})

test_that(".parse_sequences detects optional columns", {
  df <- data.frame(
    CDR3b   = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    TRBV    = c("TRBV5-1", "TRBV6-2"),
    patient = c("P1", "P2"),
    HLA     = c("A*02:01", "B*07:02"),
    counts  = c(5, 10),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_true(result$vgene.info)
  expect_true(result$patient.info)
  expect_true(result$hla.info)
  expect_true(result$count.info)
})

test_that(".parse_sequences errors on missing CDR3b column", {
  df <- data.frame(x = c("CASSLAPGATNEKLFF"), stringsAsFactors = FALSE)
  colnames(df) <- "not_cdr3"
  expect_error(immGLIPH:::.parse_sequences(df), "CDR3b")
})
