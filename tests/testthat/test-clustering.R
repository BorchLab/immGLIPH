# Tests for .cluster_gliph1() and .cluster_gliph2()

# Register sequential backend for %dopar%
foreach::registerDoSEQ()

# ---- Helpers ----------------------------------------------------------------

.make_cluster_data <- function() {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF",
            "CASSLTGGEETQYF", "CASSLGGRETQYF")
  sequences <- data.frame(
    seq_ID  = seq_along(seqs),
    CDR3b   = seqs,
    TRBV    = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2", "TRBV5-1"),
    patient = c("P1", "P1", "P2", "P2", "P1"),
    stringsAsFactors = FALSE
  )
  list(seqs = seqs, sequences = sequences)
}

# ===========================================================================
# .cluster_gliph1() tests
# ===========================================================================

# ---- Output structure -------------------------------------------------------

test_that(".cluster_gliph1 returns list with expected elements", {
  d <- .make_cluster_data()

  clone_network <- data.frame(
    V1   = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF"),
    V2   = c("CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    type = c("local", "local"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph1(
    clone_network   = clone_network,
    sequences       = d$sequences,
    not_in_global_ids = integer(0),
    seqs            = d$seqs,
    vgene.info      = TRUE,
    patient.info    = TRUE,
    global_vgene    = FALSE,
    public_tcrs     = TRUE,
    cluster_min_size = 1,
    verbose         = FALSE
  )

  expect_type(result, "list")
  expect_true("cluster_properties" %in% names(result))
  expect_true("cluster_list" %in% names(result))
  expect_true("clone_network" %in% names(result))
  expect_true("save_cluster_list_df" %in% names(result))
})

# ---- Clustering logic -------------------------------------------------------

test_that(".cluster_gliph1 forms connected components correctly", {
  d <- .make_cluster_data()

  # Two connected components: {seq1, seq2, seq3} and {seq4, seq5}
  clone_network <- data.frame(
    V1   = c("CASSLAPGATNEKLFF", "CASSLDRGEVFF",
             "CASSLTGGEETQYF"),
    V2   = c("CASSLDRGEVFF", "CASSYLAGGRNTLYF",
             "CASSLGGRETQYF"),
    type = c("local", "local", "local"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = clone_network,
    sequences        = d$sequences,
    not_in_global_ids = integer(0),
    seqs             = d$seqs,
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = FALSE,
    public_tcrs      = TRUE,
    cluster_min_size = 1,
    verbose          = FALSE
  )

  expect_s3_class(result$cluster_properties, "data.frame")
  expect_true(nrow(result$cluster_properties) >= 2)
  expect_true(all(c("cluster_size", "tag", "members") %in%
                    colnames(result$cluster_properties)))
})

# ---- cluster_min_size -------------------------------------------------------

test_that(".cluster_gliph1 filters by cluster_min_size", {
  d <- .make_cluster_data()

  # One pair (size 2) and one pair (size 2) and one isolated (size 1)
  clone_network <- data.frame(
    V1   = c("CASSLAPGATNEKLFF", "CASSLTGGEETQYF"),
    V2   = c("CASSLDRGEVFF", "CASSLGGRETQYF"),
    type = c("local", "local"),
    stringsAsFactors = FALSE
  )

  # With min_size = 3, both pairs should be filtered out
  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = clone_network,
    sequences        = d$sequences,
    not_in_global_ids = integer(0),
    seqs             = d$seqs,
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = FALSE,
    public_tcrs      = TRUE,
    cluster_min_size = 3,
    verbose          = FALSE
  )

  if (!is.null(result$cluster_properties)) {
    expect_true(all(result$cluster_properties$cluster_size >= 3))
  }
})

# ---- Empty edge list --------------------------------------------------------

test_that(".cluster_gliph1 returns NULL for empty edge list", {
  d <- .make_cluster_data()

  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = NULL,
    sequences        = d$sequences,
    not_in_global_ids = integer(0),
    seqs             = d$seqs,
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = FALSE,
    public_tcrs      = TRUE,
    cluster_min_size = 2,
    verbose          = FALSE
  )

  expect_null(result$cluster_properties)
  expect_equal(length(result$cluster_list), 0)
})

# ---- Singleton handling -----------------------------------------------------

test_that(".cluster_gliph1 adds singletons from not_in_global_ids", {
  d <- .make_cluster_data()

  clone_network <- data.frame(
    V1   = "CASSLAPGATNEKLFF",
    V2   = "CASSLDRGEVFF",
    type = "local",
    stringsAsFactors = FALSE
  )

  # Sequence at index 5 is not in any edge
  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = clone_network,
    sequences        = d$sequences,
    not_in_global_ids = c(4L, 5L),  # indices not in global edges
    seqs             = d$seqs,
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = FALSE,
    public_tcrs      = TRUE,
    cluster_min_size = 1,
    verbose          = FALSE
  )

  # The clone_network should contain singleton rows
  expect_true(any(result$clone_network$type == "singleton"))
})

# ===========================================================================
# .cluster_gliph2() tests
# ===========================================================================

test_that(".cluster_gliph2 returns list with expected elements", {
  d <- .make_cluster_data()

  local_res <- data.frame(
    motif         = c("SL"),
    start         = c(1),
    stop          = c(2),
    num_in_sample = c(3),
    num_in_ref    = c(10),
    num_fold      = c(5.0),
    fisher.score  = c(0.001),
    members       = paste(d$seqs[1:3], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = local_res,
    global_res             = NULL,
    sequences              = d$sequences,
    local_similarities     = TRUE,
    global_similarities    = FALSE,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 1,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  expect_type(result, "list")
  expect_true("merged_clusters" %in% names(result))
  expect_true("cluster_list" %in% names(result))
  expect_true("clone_network" %in% names(result))
  expect_true("save_cluster_list_df" %in% names(result))
})

test_that(".cluster_gliph2 clone_network has expected edge types", {
  d <- .make_cluster_data()

  local_res <- data.frame(
    motif         = c("SL"),
    start         = c(1),
    stop          = c(2),
    num_in_sample = c(3),
    num_in_ref    = c(10),
    num_fold      = c(5.0),
    fisher.score  = c(0.001),
    members       = paste(d$seqs[1:3], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = local_res,
    global_res             = NULL,
    sequences              = d$sequences,
    local_similarities     = TRUE,
    global_similarities    = FALSE,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 10,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  if (!is.null(result$clone_network)) {
    expect_true("type" %in% colnames(result$clone_network))
    # Non-NA edge types should be local, global, or singleton
    type_col <- result$clone_network$type[!is.na(result$clone_network$type)]
    expect_true(all(type_col %in% c("local", "global", "singleton")))
    # At least some edges should be present
    expect_true(length(type_col) > 0)
  }
})

test_that(".cluster_gliph2 cluster_min_size filters small clusters", {
  d <- .make_cluster_data()

  local_res <- data.frame(
    motif         = c("SL"),
    start         = c(1),
    stop          = c(2),
    num_in_sample = c(2),
    num_in_ref    = c(10),
    num_fold      = c(5.0),
    fisher.score  = c(0.001),
    members       = paste(d$seqs[1:2], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = local_res,
    global_res             = NULL,
    sequences              = d$sequences,
    local_similarities     = TRUE,
    global_similarities    = FALSE,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 10,
    cluster_min_size       = 5,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  # Cluster of size 2 should be eliminated by min_size = 5
  if (!is.null(result$merged_clusters)) {
    expect_true(all(result$merged_clusters$cluster_size >= 5))
  } else {
    expect_null(result$merged_clusters)
  }
})

test_that(".cluster_gliph2 handles no local and no global similarities", {
  d <- .make_cluster_data()

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = NULL,
    global_res             = NULL,
    sequences              = d$sequences,
    local_similarities     = FALSE,
    global_similarities    = FALSE,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 1,
    cluster_min_size       = 2,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  expect_null(result$merged_clusters)
  expect_equal(length(result$cluster_list), 0)
})
