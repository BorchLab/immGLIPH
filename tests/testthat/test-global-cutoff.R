test_that(".global_cutoff_stringdist returns correct structure", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLAPRQTNEKLFF", "CASSLDRGEVFF")
  motif_region <- substr(seqs, 4, nchar(seqs) - 3)
  sequences <- data.frame(
    CDR3b = seqs,
    TRBV  = c("TRBV5-1", "TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  # Use a high cutoff so at least some edges are found
  result <- .global_cutoff_stringdist(
    seqs         = seqs,
    motif_region = motif_region,
    sequences    = sequences,
    gccutoff     = 5,
    global_vgene = FALSE,
    no_cores     = 1,
    verbose      = FALSE
  )

  expect_type(result, "list")
  expect_true("edges" %in% names(result))
  expect_true("not_in_global_ids" %in% names(result))
  expect_s3_class(result$edges, "data.frame")
  if (nrow(result$edges) > 0) {
    expect_true(all(c("V1", "V2", "type") %in% colnames(result$edges)))
    expect_true(all(result$edges$type == "global"))
  }
})

test_that(".global_cutoff_stringdist returns empty for zero cutoff on different seqs", {
  seqs <- c("CASSLAPGATNEKLFF", "XYZABCDEFGHIJKLM")
  motif_region <- substr(seqs, 4, nchar(seqs) - 3)
  sequences <- data.frame(
    CDR3b = seqs,
    TRBV  = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  result <- .global_cutoff_stringdist(
    seqs         = seqs,
    motif_region = motif_region,
    sequences    = sequences,
    gccutoff     = 0,
    global_vgene = FALSE,
    no_cores     = 1,
    verbose      = FALSE
  )

  expect_equal(nrow(result$edges), 0)
})

test_that(".global_cutoff dispatches correctly", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF")
  motif_region <- substr(seqs, 4, nchar(seqs) - 3)
  sequences <- data.frame(
    CDR3b = seqs,
    TRBV  = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  # This should work regardless of whether immApex is installed
  result <- .global_cutoff(
    seqs         = seqs,
    motif_region = motif_region,
    sequences    = sequences,
    gccutoff     = 5,
    global_vgene = FALSE,
    no_cores     = 1,
    verbose      = FALSE
  )

  expect_type(result, "list")
  expect_s3_class(result$edges, "data.frame")
})

## ---- immApex backend equivalence ----------------------------------------

test_that("immApex buildNetwork backend matches stringdist backend", {
  skip_if_not_installed("immApex")

  seqs <- c("CASSLAPGATNEKLFF", "CASSLAPRQTNEKLFF",
            "CASSLDRGEVFF", "CASSLDRGQVFF")
  motif_region <- substr(seqs, 4, nchar(seqs) - 3)
  sequences <- data.frame(
    CDR3b = seqs,
    TRBV  = c("TRBV5-1", "TRBV5-1", "TRBV6-2", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  result_sd <- .global_cutoff_stringdist(
    seqs         = seqs,
    motif_region = motif_region,
    sequences    = sequences,
    gccutoff     = 3,
    global_vgene = FALSE,
    no_cores     = 1,
    verbose      = FALSE
  )

  result_apex <- .global_cutoff_immapex(
    seqs         = seqs,
    motif_region = motif_region,
    sequences    = sequences,
    gccutoff     = 3,
    global_vgene = FALSE,
    verbose      = FALSE
  )

  # Same number of edges
  expect_equal(nrow(result_sd$edges), nrow(result_apex$edges))

  # If edges found, same pairs (order may differ)
  if (nrow(result_sd$edges) > 0) {
    sd_pairs   <- paste(result_sd$edges$V1, result_sd$edges$V2)
    apex_pairs <- paste(result_apex$edges$V1, result_apex$edges$V2)
    expect_setequal(sd_pairs, apex_pairs)
  }
})
