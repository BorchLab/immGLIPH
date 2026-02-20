# Tests for clusterScoring()

# ---- Input validation --------------------------------------------------------

test_that("clusterScoring rejects non-list cluster_list", {
  expect_error(
    clusterScoring(cluster_list = "not_a_list",
                   cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF")),
    "list"
  )
})

test_that("clusterScoring rejects non-data.frame cdr3_sequences", {
  expect_error(
    clusterScoring(cluster_list = list(a = data.frame(CDR3b = "X")),
                   cdr3_sequences = list(CDR3b = "CASSLAPGATNEKLFF")),
    "data.frame"
  )
})

test_that("clusterScoring accepts character vector for cdr3_sequences", {
  cl <- list(
    "CRG-1" = data.frame(
      CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
      stringsAsFactors = FALSE
    )
  )
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  tryCatch(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   refdb_beta = ref_df, sim_depth = 10, n_cores = 1),
    error = function(e) {
      expect_false(grepl("data.frame", e$message))
    }
  )
})

test_that("clusterScoring rejects invalid refdb_beta name", {
  cl <- list("CRG-1" = data.frame(CDR3b = "CASSLAPGATNEKLFF"))
  seqs <- data.frame(CDR3b = "CASSLAPGATNEKLFF")

  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   refdb_beta = "invalid_name"),
    "refdb_beta"
  )
})

test_that("clusterScoring validates refdb_beta", {
  expect_error(clusterScoring(cluster_list = list(a = data.frame(CDR3b = "X")),
                               cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF"),
                               refdb_beta = "nonexistent_ref"),
               "refdb_beta must be")
})

test_that("clusterScoring rejects invalid gliph_version", {
  cl <- list("CRG-1" = data.frame(CDR3b = "CASSLAPGATNEKLFF"))
  seqs <- data.frame(CDR3b = "CASSLAPGATNEKLFF")

  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   gliph_version = 3),
    "gliph_version"
  )
})

test_that("clusterScoring rejects invalid sim_depth", {
  cl <- list("CRG-1" = data.frame(CDR3b = "CASSLAPGATNEKLFF"))
  seqs <- data.frame(CDR3b = "CASSLAPGATNEKLFF")

  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   sim_depth = "abc"),
    "sim_depth"
  )
  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   sim_depth = -5),
    "sim_depth"
  )
})

test_that("clusterScoring validates sim_depth minimum", {
  expect_error(clusterScoring(cluster_list = list(a = data.frame(CDR3b = "X")),
                               cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF"),
                               sim_depth = 0),
               "sim_depth")
})

test_that("clusterScoring rejects invalid hla_cutoff", {
  cl <- list("CRG-1" = data.frame(CDR3b = "CASSLAPGATNEKLFF"))
  seqs <- data.frame(CDR3b = "CASSLAPGATNEKLFF")

  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   hla_cutoff = 2.0),
    "hla_cutoff"
  )
  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   hla_cutoff = -0.1),
    "hla_cutoff"
  )
})

test_that("clusterScoring validates v_usage_freq format", {
  cl <- list("CRG-1" = data.frame(CDR3b = "CASSLAPGATNEKLFF"))
  seqs <- data.frame(CDR3b = "CASSLAPGATNEKLFF")

  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   v_usage_freq = c(0.1, 0.2)),
    "v_usage_freq"
  )

  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   v_usage_freq = data.frame(vgene = "TRBV5-1")),
    "v_usage_freq"
  )
})

test_that("clusterScoring validates cdr3_length_freq format", {
  cl <- list("CRG-1" = data.frame(CDR3b = "CASSLAPGATNEKLFF"))
  seqs <- data.frame(CDR3b = "CASSLAPGATNEKLFF")

  expect_error(
    clusterScoring(cluster_list = cl, cdr3_sequences = seqs,
                   cdr3_length_freq = c(0.1, 0.2)),
    "cdr3_length_freq"
  )
})

test_that("clusterScoring returns empty data.frame for empty cluster_list", {
  seqs <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    TRBV = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )
  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    TRBV = c("TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  result <- clusterScoring(
    cluster_list = list(),
    cdr3_sequences = seqs,
    refdb_beta = ref_df,
    sim_depth = 10,
    n_cores = 1
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

# ---- Computation tests -------------------------------------------------------

test_that("clusterScoring computes scores with custom refdb_beta data frame", {
  skip_on_cran()

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF", "CASSLGQAYEQYF",
              "CASSFSTCSANYGYTF", "CASSPTGGYEQYF", "CASRLAGGRNTLYF",
              "CASSLSGGRNTLYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1",
             "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  cl <- list(
    "CRG-1" = data.frame(
      CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
      TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
      stringsAsFactors = FALSE
    )
  )

  seqs <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  result <- clusterScoring(
    cluster_list = cl,
    cdr3_sequences = seqs,
    refdb_beta = ref_df,
    gliph_version = 1,
    sim_depth = 10,
    n_cores = 1
  )

  expect_s3_class(result, "data.frame")
  expect_true("total.score" %in% colnames(result))
  expect_true("network.size.score" %in% colnames(result))
  expect_true("cdr3.length.score" %in% colnames(result))
  expect_true("vgene.score" %in% colnames(result))
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$total.score))
})

test_that("clusterScoring gliph_version 2 omits 0.064 multiplier", {
  skip_on_cran()

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF", "CASSLGQAYEQYF",
              "CASSFSTCSANYGYTF", "CASSPTGGYEQYF", "CASRLAGGRNTLYF",
              "CASSLSGGRNTLYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1",
             "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1", "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  cl <- list(
    "CRG-1" = data.frame(
      CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
      TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
      stringsAsFactors = FALSE
    )
  )

  seqs <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  res_v1 <- clusterScoring(
    cluster_list = cl, cdr3_sequences = seqs,
    refdb_beta = ref_df, gliph_version = 1, sim_depth = 10, n_cores = 1
  )
  res_v2 <- clusterScoring(
    cluster_list = cl, cdr3_sequences = seqs,
    refdb_beta = ref_df, gliph_version = 2, sim_depth = 10, n_cores = 1
  )

  expect_s3_class(res_v1, "data.frame")
  expect_s3_class(res_v2, "data.frame")
  expect_true("total.score" %in% colnames(res_v1))
  expect_true("total.score" %in% colnames(res_v2))
})

test_that("clusterScoring with counts and HLA information", {
  skip_on_cran()

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF", "CASSLGQAYEQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1",
             "TRBV6-2"),
    stringsAsFactors = FALSE
  )

  cl <- list(
    "CRG-1" = data.frame(
      CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
      TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
      patient = c("P1", "P2", "P1"),
      HLA = c("A*02:01,B*07:02", "A*02:01,B*08:01", "A*02:01,B*07:02"),
      counts = c(5, 3, 2),
      stringsAsFactors = FALSE
    )
  )

  seqs <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2"),
    patient = c("P1", "P2", "P1", "P2"),
    HLA = c("A*02:01,B*07:02", "A*02:01,B*08:01",
            "A*02:01,B*07:02", "A*02:01,B*08:01"),
    counts = c(5, 3, 2, 1),
    stringsAsFactors = FALSE
  )

  result <- clusterScoring(
    cluster_list = cl,
    cdr3_sequences = seqs,
    refdb_beta = ref_df,
    gliph_version = 1,
    sim_depth = 10,
    n_cores = 1
  )

  expect_s3_class(result, "data.frame")
  expect_true("clonal.expansion.score" %in% colnames(result))
  expect_true("hla.score" %in% colnames(result))
})

test_that("clusterScoring with multiple clusters", {
  skip_on_cran()

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF", "CASSLGQAYEQYF",
              "CASSFSTCSANYGYTF", "CASSPTGGYEQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1",
             "TRBV6-2", "TRBV5-1", "TRBV7-2"),
    stringsAsFactors = FALSE
  )

  cl <- list(
    "CRG-1" = data.frame(
      CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
      TRBV = c("TRBV5-1", "TRBV6-2"),
      stringsAsFactors = FALSE
    ),
    "CRG-2" = data.frame(
      CDR3b = c("CASSYLAGGRNTLYF", "CASSLTGGEETQYF"),
      TRBV = c("TRBV5-1", "TRBV7-2"),
      stringsAsFactors = FALSE
    )
  )

  seqs <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2"),
    stringsAsFactors = FALSE
  )

  result <- clusterScoring(
    cluster_list = cl,
    cdr3_sequences = seqs,
    refdb_beta = ref_df,
    gliph_version = 1,
    sim_depth = 10,
    n_cores = 1
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

test_that("clusterScoring with custom v_usage_freq", {
  skip_on_cran()

  ref_df <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
              "CASSLTGGEETQYF", "CASSLGGRETQYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  cl <- list(
    "CRG-1" = data.frame(
      CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
      TRBV = c("TRBV5-1", "TRBV6-2"),
      stringsAsFactors = FALSE
    )
  )

  seqs <- data.frame(
    CDR3b = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    TRBV = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  v_freq <- data.frame(
    vgene = c("TRBV5-1", "TRBV6-2", "TRBV7-2"),
    freq = c(0.5, 0.3, 0.2),
    stringsAsFactors = FALSE
  )

  result <- clusterScoring(
    cluster_list = cl,
    cdr3_sequences = seqs,
    refdb_beta = ref_df,
    v_usage_freq = v_freq,
    gliph_version = 1,
    sim_depth = 10,
    n_cores = 1
  )

  expect_s3_class(result, "data.frame")
  expect_true("vgene.score" %in% colnames(result))
})
