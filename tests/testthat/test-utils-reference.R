# Tests for utility functions in utils-reference.R

test_that(".valid_reference_names returns expected names", {
  names <- immGLIPH:::.valid_reference_names()
  expect_type(names, "character")
  expect_equal(length(names), 10)

  # All expected names present
expect_true("human_v1.0_CD4" %in% names)
  expect_true("human_v1.0_CD8" %in% names)
  expect_true("human_v1.0_CD48" %in% names)
  expect_true("human_v2.0_CD4" %in% names)
  expect_true("human_v2.0_CD8" %in% names)
  expect_true("human_v2.0_CD48" %in% names)
  expect_true("mouse_v1.0_CD4" %in% names)
  expect_true("mouse_v1.0_CD8" %in% names)
  expect_true("mouse_v1.0_CD48" %in% names)
  expect_true("gliph_reference" %in% names)
})

test_that(".prepare_motif_region trims boundaries correctly", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF")
  # With structural boundaries (trim 3 from each end)
  result <- immGLIPH:::.prepare_motif_region(seqs, structboundaries = TRUE, boundary_size = 3)
  expect_equal(result[1], substr("CASSLAPGATNEKLFF", 4, nchar("CASSLAPGATNEKLFF") - 3))
  expect_equal(result[1], "SLAPGATNEK")
  expect_equal(result[2], "SLDRGE")
})

test_that(".prepare_motif_region returns full sequences when structboundaries is FALSE", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF")
  result <- immGLIPH:::.prepare_motif_region(seqs, structboundaries = FALSE, boundary_size = 3)
  expect_equal(result, seqs)
})

test_that(".prepare_motif_region handles different boundary sizes", {
  seqs <- c("CASSLAPGATNEKLFF")
  result1 <- immGLIPH:::.prepare_motif_region(seqs, TRUE, 1)
  result2 <- immGLIPH:::.prepare_motif_region(seqs, TRUE, 5)
  expect_equal(result1, substr(seqs, 2, nchar(seqs) - 1))
  expect_equal(result2, substr(seqs, 6, nchar(seqs) - 5))
})

test_that(".get_blosum_vec returns character vector of AA pairs", {
  skip_if_not_installed("immApex")
  bv <- immGLIPH:::.get_blosum_vec()
  expect_type(bv, "character")
  expect_true(length(bv) > 0)
  # All entries should be exactly 2 characters
  expect_true(all(nchar(bv) == 2))
  # Identity pairs (e.g., "AA", "CC") should be present
  expect_true("AA" %in% bv)
  expect_true("LL" %in% bv)
  expect_true("FF" %in% bv)
})

test_that(".load_reference handles data frame input", {
  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.load_reference(
    refdb_beta = ref_df,
    accept_CF = TRUE,
    min_seq_length = 8,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true("refseqs" %in% names(result))
  expect_true("ref_vgenes" %in% names(result))
  expect_true("refseqs_motif" %in% names(result))
  expect_true("refseqs_df" %in% names(result))
  expect_type(result$refseqs, "character")
  expect_true(length(result$refseqs) > 0)
  # All sequences should start with C and end with F (accept_CF = TRUE)
  expect_true(all(grepl("^C.*F$", result$refseqs)))
})

test_that(".load_reference filters by min_seq_length", {
  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASF", "CASSLDRGEVFF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.load_reference(
    refdb_beta = ref_df,
    accept_CF = TRUE,
    min_seq_length = 8,
    verbose = FALSE
  )

  # "CASF" is only 4 chars, should be filtered out
  expect_true(all(nchar(result$refseqs) >= 8))
})

test_that(".load_reference errors on empty data frame", {
  ref_df <- data.frame(CDR3b = character(0), TRBV = character(0),
                        stringsAsFactors = FALSE)
  expect_error(
    immGLIPH:::.load_reference(ref_df, verbose = FALSE),
    "empty"
  )
})

test_that(".load_reference handles data frame without CDR3b column name", {
  ref_df <- data.frame(
    seq = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    vgene = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.load_reference(
    refdb_beta = ref_df,
    accept_CF = TRUE,
    min_seq_length = 8,
    verbose = FALSE
  )

  # Should still work, using first column as CDR3b
  expect_true(length(result$refseqs) > 0)
})

test_that(".load_reference extracts vgenes when vgene_stratify is TRUE", {
  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.load_reference(
    refdb_beta = ref_df,
    accept_CF = TRUE,
    min_seq_length = 8,
    vgene_stratify = TRUE,
    verbose = FALSE
  )

  expect_false(is.null(result$ref_vgenes))
  expect_type(result$ref_vgenes, "character")
})

test_that(".load_reference applies accept_CF filtering", {
  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "AASSLAPGATNEKLFG", "CASSLDRGEVFF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  # With accept_CF = TRUE, non-C/F sequences should be removed
  result_cf <- immGLIPH:::.load_reference(
    refdb_beta = ref_df, accept_CF = TRUE, min_seq_length = 8, verbose = FALSE
  )
  expect_true(all(grepl("^C.*F$", result_cf$refseqs)))

  # With accept_CF = FALSE, all valid AA sequences should remain
  result_no_cf <- immGLIPH:::.load_reference(
    refdb_beta = ref_df, accept_CF = FALSE, min_seq_length = 8, verbose = FALSE
  )
  expect_true(length(result_no_cf$refseqs) >= length(result_cf$refseqs))
})
