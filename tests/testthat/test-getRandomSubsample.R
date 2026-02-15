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
  # All values should come from the reference
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
