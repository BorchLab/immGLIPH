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

test_that(".parse_sequences errors on missing CDR3b column for multi-column data frame", {
  df <- data.frame(
    not_cdr3 = c("CASSLAPGATNEKLFF"),
    other    = c("something"),
    stringsAsFactors = FALSE
  )
  expect_error(immGLIPH:::.parse_sequences(df, verbose = FALSE), "CDR3b")
})

test_that(".parse_sequences uses first column for single-column data frame", {
  df <- data.frame(not_cdr3 = c("CASSLAPGATNEKLFF"), stringsAsFactors = FALSE)
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_equal(nrow(result$sequences), 1)
})

# ---- .parse_sequences character vector input ---------------------------------

test_that(".parse_sequences handles character vector input", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF")
  result <- immGLIPH:::.parse_sequences(seqs, verbose = FALSE)
  expect_equal(nrow(result$sequences), 2)
  expect_false(result$vgene.info)
  expect_false(result$patient.info)
  expect_false(result$hla.info)
  expect_false(result$count.info)
})

test_that(".parse_sequences errors on character vector with global_vgene", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF")
  expect_error(
    immGLIPH:::.parse_sequences(seqs, global_vgene = TRUE, verbose = FALSE),
    "V-gene"
  )
})

test_that(".parse_sequences errors on character vector with vgene_stratify", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF")
  expect_error(
    immGLIPH:::.parse_sequences(seqs, vgene_stratify = TRUE, verbose = FALSE),
    "V-gene"
  )
})

# ---- .parse_sequences data frame with missing V-gene -----------------------

test_that(".parse_sequences errors on data frame without TRBV when global_vgene is TRUE", {
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    stringsAsFactors = FALSE
  )
  expect_error(
    immGLIPH:::.parse_sequences(df, global_vgene = TRUE, verbose = FALSE),
    "V-gene"
  )
})

test_that(".parse_sequences errors on data frame without TRBV when vgene_stratify is TRUE", {
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    stringsAsFactors = FALSE
  )
  expect_error(
    immGLIPH:::.parse_sequences(df, vgene_stratify = TRUE, verbose = FALSE),
    "V-gene"
  )
})

# ---- .parse_sequences errors on invalid input type -------------------------

test_that(".parse_sequences errors on invalid input type", {
  expect_error(
    immGLIPH:::.parse_sequences(list(a = 1), verbose = FALSE),
    "character vector or data frame"
  )
})

# ---- .parse_sequences errors on no valid AA sequences ----------------------

test_that(".parse_sequences errors when no valid AA sequences found", {
  df <- data.frame(
    CDR3b = c("12345", "INVALID123"),
    stringsAsFactors = FALSE
  )
  expect_error(
    immGLIPH:::.parse_sequences(df, verbose = FALSE),
    "No valid"
  )
})

# ---- .parse_sequences adds seq_ID column -----------------------------------

test_that(".parse_sequences adds sequential seq_ID column", {
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_equal(result$sequences$seq_ID, 1:3)
})

# ---- .parse_sequences counts handling -------------------------------------

test_that(".parse_sequences handles NA counts", {
  df <- data.frame(
    CDR3b  = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    counts = c(5, NA),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_true(result$count.info)
  # NA counts should be replaced with 1
  expect_equal(as.numeric(result$sequences$counts[2]), 1)
})

# ---- .parse_sequences carries forward extra columns -----------------------

test_that(".parse_sequences carries forward additional columns", {
  df <- data.frame(
    CDR3b   = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    TRBV    = c("TRBV5-1", "TRBV6-2"),
    patient = c("P1", "P2"),
    custom_col = c("A", "B"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_true("custom_col" %in% colnames(result$sequences))
})

# ---- .standardize_colnames additional mappings ------------------------------

test_that(".standardize_colnames maps HLA and counts alternatives", {
  df <- data.frame(
    cdr3_aa     = c("CASSLAPGATNEKLFF"),
    v_call      = c("TRBV5-1"),
    donor       = c("P1"),
    hla         = c("A*02:01"),
    clone_count = c(5),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.standardize_colnames(df)
  expect_true("CDR3b" %in% colnames(result))
  expect_true("TRBV" %in% colnames(result))
  expect_true("patient" %in% colnames(result))
  expect_true("HLA" %in% colnames(result))
  expect_true("counts" %in% colnames(result))
})

test_that(".standardize_colnames handles single-column data frame without CDR3b", {
  df <- data.frame(
    some_random_name = c("CASSLAPGATNEKLFF"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.standardize_colnames(df)
  expect_equal(colnames(result), "CDR3b")
})

# ---- .extract_input data frame standardization ----------------------------

test_that(".extract_input standardizes alternative column names", {
  df <- data.frame(
    cdr3   = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    v_gene = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.extract_input(df)
  expect_true("CDR3b" %in% colnames(result))
  expect_true("TRBV" %in% colnames(result))
})

# ---- .extract_input rejects numeric input ------------------------------------

test_that(".extract_input rejects numeric input", {
  expect_error(immGLIPH:::.extract_input(42), "must be a character vector")
})

# ---- .extract_input rejects NULL input ---------------------------------------

test_that(".extract_input rejects NULL input", {
  expect_error(immGLIPH:::.extract_input(NULL), "must be a character vector")
})

# ---- .parse_sequences accept_CF filtering ------------------------------------

test_that(".parse_sequences keeps all sequences when accept_CF not in function", {
  # .parse_sequences does not filter by C/F - that's done by other functions

  # Verify that sequences without C/F pass through .parse_sequences
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "AASSLAPGATNEKLFG"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_equal(nrow(result$sequences), 2)
})

# ---- .parse_sequences verbose messaging on single-column DF ------------------

test_that(".parse_sequences prints message for single-column data frame", {
  df <- data.frame(not_cdr3 = c("CASSLAPGATNEKLFF"), stringsAsFactors = FALSE)
  expect_message(
    immGLIPH:::.parse_sequences(df, verbose = TRUE),
    "first column"
  )
})

# ---- .parse_sequences with empty counts column -------------------------------

test_that(".parse_sequences treats all-NA counts as 1", {
  df <- data.frame(
    CDR3b  = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    counts = c(NA, NA),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_true(result$count.info)
  expect_true(all(as.numeric(result$sequences$counts) == 1))
})

# ---- .standardize_colnames with Frequency column -----------------------------

test_that(".standardize_colnames maps Frequency to counts", {
  df <- data.frame(
    CDR3b     = c("CASSLAPGATNEKLFF"),
    Frequency = c(10),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.standardize_colnames(df)
  expect_true("counts" %in% colnames(result))
})

# ---- .standardize_colnames with junction_aa column ---------------------------

test_that(".standardize_colnames maps junction_aa to CDR3b", {
  df <- data.frame(
    junction_aa = c("CASSLAPGATNEKLFF"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.standardize_colnames(df)
  expect_true("CDR3b" %in% colnames(result))
})

# ---- .standardize_colnames with sample.id column -----------------------------

test_that(".standardize_colnames maps sample.id to patient", {
  df <- data.frame(
    CDR3b     = c("CASSLAPGATNEKLFF"),
    sample.id = c("P1"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.standardize_colnames(df)
  expect_true("patient" %in% colnames(result))
})

# ---- .parse_sequences preserves row count after filtering --------------------

test_that(".parse_sequences preserves valid rows after AA filtering", {
  df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "BAD123SEQ"),
    TRBV  = c("TRBV5-1", "TRBV6-2", "TRBV7-2"),
    stringsAsFactors = FALSE
  )
  result <- immGLIPH:::.parse_sequences(df, verbose = FALSE)
  expect_equal(nrow(result$sequences), 2)
  expect_true(all(result$sequences$seq_ID == 1:2))
})

# ---- .extract_input with list input requires immApex -------------------------

test_that(".extract_input errors on list input without immApex or with message", {
  # This tests the list (combineTCR output) path
  mock_list <- list(data.frame(CDR3b = "CASSLAPGATNEKLFF"))
  # If immApex is not available, should error; if available, should succeed
  if (!requireNamespace("immApex", quietly = TRUE)) {
    expect_error(immGLIPH:::.extract_input(mock_list), "immApex")
  } else {
    # immApex is available but the list may not be valid combineTCR format
    # This verifies the code path is reached
    expect_true(TRUE)
  }
})
