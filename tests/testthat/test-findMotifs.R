# Tests for findMotifs()

# ---- Core findMotifs tests ---------------------------------------------------

test_that("findMotifs returns data frame with motif and V1 columns", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  result <- findMotifs(seqs)
  expect_s3_class(result, "data.frame")
  expect_true("motif" %in% colnames(result))
  expect_true("V1" %in% colnames(result))
})

test_that("findMotifs finds known 2-mer motifs", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF")
  result <- findMotifs(seqs, q = 2)
  # "CA" appears at position 1 in both sequences
  ca_row <- result[result$motif == "CA", ]
  expect_equal(nrow(ca_row), 1)
  expect_equal(ca_row$V1, 2)
})

test_that("findMotifs respects kmer_mindepth", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  result <- findMotifs(seqs, q = 2, kmer_mindepth = 3)
  # All returned motifs should have count >= 3
  expect_true(all(result$V1 >= 3))
})

test_that("findMotifs handles multiple motif lengths", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  result <- findMotifs(seqs, q = c(2, 3))
  # Should have 2-char and 3-char motifs
  motif_lens <- nchar(result$motif)
  expect_true(any(motif_lens == 2))
  expect_true(any(motif_lens == 3))
})

test_that("findMotifs with discontinuous motifs includes dots", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  result <- findMotifs(seqs, q = 2, discontinuous = TRUE)
  disc <- result[grep("\\.", result$motif), ]
  expect_true(nrow(disc) > 0)
  # All discontinuous motifs should contain at least one dot
  expect_true(all(grepl("\\.", disc$motif)))
})

test_that("findMotifs returns empty data frame for empty input", {
  result <- findMotifs(character(0), q = 2)
  expect_true(is.null(result) || nrow(result) == 0)
})

# ---- stringdist fallback tests -----------------------------------------------

test_that(".find_motifs_stringdist returns data frame with motif and V1", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  result <- immGLIPH:::.find_motifs_stringdist(seqs, q = 2)
  expect_s3_class(result, "data.frame")
  expect_true("motif" %in% colnames(result))
  expect_true("V1" %in% colnames(result))
})

test_that(".find_motifs_stringdist handles discontinuous motifs", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  result <- immGLIPH:::.find_motifs_stringdist(seqs, q = 2, discontinuous = TRUE)
  disc <- result[grep("\\.", result$motif), ]
  expect_true(nrow(disc) > 0)
})

# ---- immApex backend equivalence tests ---------------------------------------

test_that("immApex backend produces same motifs as stringdist backend", {
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
            "CASSLAPGATNEKLFF", "CASSLDRGEVFF")

  # Force stringdist path
  result_sd <- immGLIPH:::.find_motifs_stringdist(seqs, q = 2:3)
  # immApex path (findMotifs uses immApex when available)
  result_apex <- findMotifs(seqs, q = 2:3)

  # Both should have the same columns

  expect_true(all(c("motif", "V1") %in% colnames(result_sd)))
  expect_true(all(c("motif", "V1") %in% colnames(result_apex)))

  # Same motifs should be found, same counts
  merged <- merge(result_sd, result_apex, by = "motif",
                  suffixes = c(".sd", ".apex"))
  expect_equal(merged$V1.sd, merged$V1.apex)
})

test_that("immApex backend handles kmer_mindepth correctly", {
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")

  result_sd   <- immGLIPH:::.find_motifs_stringdist(seqs, q = 2, kmer_mindepth = 3)
  result_apex <- findMotifs(seqs, q = 2, kmer_mindepth = 3)

  # Both should only return motifs with V1 >= 3
  expect_true(all(result_sd$V1 >= 3))
  expect_true(all(result_apex$V1 >= 3))

  merged <- merge(result_sd, result_apex, by = "motif",
                  suffixes = c(".sd", ".apex"))
  expect_equal(merged$V1.sd, merged$V1.apex)
})

test_that("immApex backend handles discontinuous motifs correctly", {
  skip_if_not_installed("immApex")
  skip_if(!exists("calculateMotif", asNamespace("immApex")),
          "immApex::calculateMotif not available")

  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")

  result_sd   <- immGLIPH:::.find_motifs_stringdist(seqs, q = 2, discontinuous = TRUE)
  result_apex <- findMotifs(seqs, q = 2, discontinuous = TRUE)

  # Both should contain discontinuous motifs (with dots)
  disc_sd   <- result_sd[grep("\\.", result_sd$motif), ]
  disc_apex <- result_apex[grep("\\.", result_apex$motif), ]
  expect_true(nrow(disc_sd) > 0)
  expect_true(nrow(disc_apex) > 0)

  # Merge on motif and compare counts
  merged <- merge(disc_sd, disc_apex, by = "motif",
                  suffixes = c(".sd", ".apex"))
  expect_equal(merged$V1.sd, merged$V1.apex)
})
