# Tests for getRandomSubsample()

test_that("getRandomSubsample returns correct length without stratification", {
  ref_motifs <- paste0("MOTIF", seq_len(1000))
  sample_motifs <- paste0("MOTIF", seq_len(50))

  result <- getRandomSubsample(
    cdr3_len_stratify = FALSE,
    vgene_stratify = FALSE,
    refseqs_motif_region = ref_motifs,
    motif_region = sample_motifs,
    motif_lengths_list = list(),
    ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = list(),
    ref_motif_vgenes_id_list = list(),
    ref_lengths_vgenes_list = list(),
    lengths_vgenes_list = list()
  )

  expect_equal(length(result), length(sample_motifs))
  expect_true(all(result %in% ref_motifs))
})

test_that("getRandomSubsample with CDR3 length stratification returns correct length", {
  set.seed(42)
  ref_motifs <- c(rep("AABB", 200), rep("AABBCC", 300), rep("AABBCCDD", 500))
  sample_motifs <- c(rep("XXYY", 5), rep("XXYYCC", 10), rep("XXYYCCDD", 15))

  motif_lengths_list <- list("4" = 5, "6" = 10, "8" = 15)
  ref_motif_lengths_id_list <- list(
    "4" = which(nchar(ref_motifs) == 4),
    "6" = which(nchar(ref_motifs) == 6),
    "8" = which(nchar(ref_motifs) == 8)
  )

  result <- getRandomSubsample(
    cdr3_len_stratify = TRUE,
    vgene_stratify = FALSE,
    refseqs_motif_region = ref_motifs,
    motif_region = sample_motifs,
    motif_lengths_list = motif_lengths_list,
    ref_motif_lengths_id_list = ref_motif_lengths_id_list,
    motif_region_vgenes_list = list(),
    ref_motif_vgenes_id_list = list(),
    ref_lengths_vgenes_list = list(),
    lengths_vgenes_list = list()
  )

  expect_equal(length(result), length(sample_motifs))
  expect_true(all(result %in% ref_motifs))
})

test_that("getRandomSubsample with V-gene stratification returns correct length", {
  set.seed(42)
  ref_motifs <- paste0("REF", seq_len(500))
  sample_motifs <- paste0("SAM", seq_len(20))

  motif_region_vgenes_list <- list("TRBV5-1" = 12, "TRBV6-2" = 8)
  ref_motif_vgenes_id_list <- list(
    "TRBV5-1" = seq_len(300),
    "TRBV6-2" = 301:500
  )

  result <- getRandomSubsample(
    cdr3_len_stratify = FALSE,
    vgene_stratify = TRUE,
    refseqs_motif_region = ref_motifs,
    motif_region = sample_motifs,
    motif_lengths_list = list(),
    ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = motif_region_vgenes_list,
    ref_motif_vgenes_id_list = ref_motif_vgenes_id_list,
    ref_lengths_vgenes_list = list(),
    lengths_vgenes_list = list()
  )

  expect_equal(length(result), length(sample_motifs))
  expect_true(all(result %in% ref_motifs))
})

test_that("getRandomSubsample with both CDR3 and V-gene stratification", {
  set.seed(42)
  ref_motifs <- paste0("REF", seq_len(500))
  sample_motifs <- paste0("SAM", seq_len(20))

  motif_lengths_list <- list("10" = 12, "12" = 8)
  ref_motif_lengths_id_list <- list(
    "10" = seq_len(250),
    "12" = 251:500
  )
  motif_region_vgenes_list <- list("TRBV5-1" = 12, "TRBV6-2" = 8)
  ref_motif_vgenes_id_list <- list(
    "TRBV5-1" = seq_len(300),
    "TRBV6-2" = 301:500
  )
  lengths_vgenes_list <- list(
    "10" = list("TRBV5-1" = 7, "TRBV6-2" = 5),
    "12" = list("TRBV5-1" = 5, "TRBV6-2" = 3)
  )
  ref_lengths_vgenes_list <- list(
    "10" = list("TRBV5-1" = seq_len(100), "TRBV6-2" = 101:200),
    "12" = list("TRBV5-1" = 201:350, "TRBV6-2" = 351:500)
  )

  result <- getRandomSubsample(
    cdr3_len_stratify = TRUE,
    vgene_stratify = TRUE,
    refseqs_motif_region = ref_motifs,
    motif_region = sample_motifs,
    motif_lengths_list = motif_lengths_list,
    ref_motif_lengths_id_list = ref_motif_lengths_id_list,
    motif_region_vgenes_list = motif_region_vgenes_list,
    ref_motif_vgenes_id_list = ref_motif_vgenes_id_list,
    ref_lengths_vgenes_list = ref_lengths_vgenes_list,
    lengths_vgenes_list = lengths_vgenes_list
  )

  expect_equal(length(result), length(sample_motifs))
  expect_true(all(result %in% ref_motifs))
})

test_that("getRandomSubsample CDR3 length stratification handles overflow", {
  set.seed(42)
  # More sample motifs of length 4 than ref has
  ref_motifs <- c(rep("AABB", 3), rep("AABBCC", 300))
  sample_motifs <- c(rep("XXYY", 10), rep("XXYYCC", 5))

  motif_lengths_list <- list("4" = 10, "6" = 5)
  ref_motif_lengths_id_list <- list(
    "4" = which(nchar(ref_motifs) == 4),
    "6" = which(nchar(ref_motifs) == 6)
  )

  result <- getRandomSubsample(
    cdr3_len_stratify = TRUE,
    vgene_stratify = FALSE,
    refseqs_motif_region = ref_motifs,
    motif_region = sample_motifs,
    motif_lengths_list = motif_lengths_list,
    ref_motif_lengths_id_list = ref_motif_lengths_id_list,
    motif_region_vgenes_list = list(),
    ref_motif_vgenes_id_list = list(),
    ref_lengths_vgenes_list = list(),
    lengths_vgenes_list = list()
  )

  # Should still return correct total length even if some strata overflow
  expect_equal(length(result), length(sample_motifs))
  expect_true(all(result %in% ref_motifs))
})

test_that("getRandomSubsample V-gene stratification handles overflow", {
  set.seed(42)
  ref_motifs <- paste0("REF", seq_len(500))
  sample_motifs <- paste0("SAM", seq_len(20))

  # More sample TRBV5-1 than ref has
  motif_region_vgenes_list <- list("TRBV5-1" = 15, "TRBV6-2" = 5)
  ref_motif_vgenes_id_list <- list(
    "TRBV5-1" = seq_len(10),  # only 10 available, need 15
    "TRBV6-2" = 11:500
  )

  result <- getRandomSubsample(
    cdr3_len_stratify = FALSE,
    vgene_stratify = TRUE,
    refseqs_motif_region = ref_motifs,
    motif_region = sample_motifs,
    motif_lengths_list = list(),
    ref_motif_lengths_id_list = list(),
    motif_region_vgenes_list = motif_region_vgenes_list,
    ref_motif_vgenes_id_list = ref_motif_vgenes_id_list,
    ref_lengths_vgenes_list = list(),
    lengths_vgenes_list = list()
  )

  expect_equal(length(result), length(sample_motifs))
  expect_true(all(result %in% ref_motifs))
})
