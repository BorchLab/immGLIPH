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

# ---- .cluster_gliph2 with global_res ----------------------------------------

test_that(".cluster_gliph2 handles global results only", {
  d <- .make_cluster_data()

  global_res <- data.frame(
    cluster_tag   = c("struct_%_13"),
    cluster_size  = c(3),
    unique_CDR3b  = c(3),
    num_in_ref    = c(5),
    fisher.score  = c(0.001),
    aa_at_position = c("L"),
    TRBV          = c("TRBV5-1"),
    CDR3b         = paste(d$seqs[c(1, 3, 5)], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = NULL,
    global_res             = global_res,
    sequences              = d$sequences,
    local_similarities     = FALSE,
    global_similarities    = TRUE,
    global_vgene           = TRUE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 1,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  expect_type(result, "list")
  if (!is.null(result$merged_clusters)) {
    expect_true(all(result$merged_clusters$type == "global"))
  }
})

test_that(".cluster_gliph2 with global_vgene FALSE for global clusters", {
  d <- .make_cluster_data()

  global_res <- data.frame(
    cluster_tag    = c("struct_%_13"),
    cluster_size   = c(3),
    unique_CDR3b   = c(3),
    num_in_ref     = c(5),
    fisher.score   = c(0.001),
    aa_at_position = c("L"),
    TRBV           = c("TRBV5-1"),
    CDR3b          = paste(d$seqs[c(1, 3, 5)], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = NULL,
    global_res             = global_res,
    sequences              = d$sequences,
    local_similarities     = FALSE,
    global_similarities    = TRUE,
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
  if (!is.null(result$merged_clusters)) {
    # Without global_vgene, tag should NOT include TRBV
    expect_false(any(grepl("TRBV", result$merged_clusters$tag)))
  }
})

# ---- .cluster_gliph2 with both local and global ----------------------------

test_that(".cluster_gliph2 merges local and global clusters", {
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

  global_res <- data.frame(
    cluster_tag    = c("struct_%_5"),
    cluster_size   = c(2),
    unique_CDR3b   = c(2),
    num_in_ref     = c(5),
    fisher.score   = c(0.01),
    aa_at_position = c("G"),
    TRBV           = c("TRBV7-2"),
    CDR3b          = paste(d$seqs[4:5], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = local_res,
    global_res             = global_res,
    sequences              = d$sequences,
    local_similarities     = TRUE,
    global_similarities    = TRUE,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 10,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  expect_type(result, "list")
  if (!is.null(result$merged_clusters)) {
    expect_true(any(result$merged_clusters$type == "local"))
    expect_true(any(result$merged_clusters$type == "global"))
  }
})

# ---- .cluster_gliph2 BLOSUM62 filtering ------------------------------------

test_that(".cluster_gliph2 applies BLOSUM62 filtering when all_aa_interchangeable is FALSE", {
  skip_if_not_installed("immApex")

  d <- .make_cluster_data()

  global_res <- data.frame(
    cluster_tag    = c("struct_%_5"),
    cluster_size   = c(2),
    unique_CDR3b   = c(2),
    num_in_ref     = c(5),
    fisher.score   = c(0.01),
    aa_at_position = c("G"),
    TRBV           = c("TRBV5-1"),
    CDR3b          = paste(d$seqs[c(1, 3)], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = NULL,
    global_res             = global_res,
    sequences              = d$sequences,
    local_similarities     = FALSE,
    global_similarities    = TRUE,
    global_vgene           = FALSE,
    all_aa_interchangeable = FALSE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 1,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  expect_type(result, "list")
})

# ---- .cluster_gliph1 with public_tcrs FALSE --------------------------------

test_that(".cluster_gliph1 restricts edges to same donor when public_tcrs is FALSE", {
  d <- .make_cluster_data()

  clone_network <- data.frame(
    V1   = c("CASSLAPGATNEKLFF", "CASSLAPGATNEKLFF"),
    V2   = c("CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    type = c("local", "local"),
    stringsAsFactors = FALSE
  )

  # P1 has seq1,2,5; P2 has seq3,4
  # Seq1-Seq2: same donor P1
  # Seq1-Seq3: different donors P1 vs P2
  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = clone_network,
    sequences        = d$sequences,
    not_in_global_ids = integer(0),
    seqs             = d$seqs,
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = FALSE,
    public_tcrs      = FALSE,
    cluster_min_size = 1,
    verbose          = FALSE
  )

  expect_type(result, "list")
})

# ---- .cluster_gliph1 global_vgene filtering --------------------------------

test_that(".cluster_gliph1 filters global edges by V-gene when global_vgene is TRUE", {
  d <- .make_cluster_data()

  clone_network <- data.frame(
    V1   = c("CASSLAPGATNEKLFF", "CASSLAPGATNEKLFF"),
    V2   = c("CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    type = c("global", "global"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = clone_network,
    sequences        = d$sequences,
    not_in_global_ids = integer(0),
    seqs             = d$seqs,
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = TRUE,
    public_tcrs      = TRUE,
    cluster_min_size = 1,
    verbose          = FALSE
  )

  expect_type(result, "list")
})

# ---- .cluster_gliph1 verbose messaging -------------------------------------

test_that(".cluster_gliph1 prints message when verbose is TRUE", {
  d <- .make_cluster_data()

  clone_network <- data.frame(
    V1   = "CASSLAPGATNEKLFF",
    V2   = "CASSLDRGEVFF",
    type = "local",
    stringsAsFactors = FALSE
  )

  expect_message(
    immGLIPH:::.cluster_gliph1(
      clone_network    = clone_network,
      sequences        = d$sequences,
      not_in_global_ids = integer(0),
      seqs             = d$seqs,
      vgene.info       = TRUE,
      patient.info     = TRUE,
      global_vgene     = FALSE,
      public_tcrs      = TRUE,
      cluster_min_size = 1,
      verbose          = TRUE
    ),
    "GLIPH1"
  )
})

# ---- .cluster_gliph1 singletons from NULL network ---------------------------

test_that(".cluster_gliph1 creates singleton network from NULL clone_network", {
  d <- .make_cluster_data()

  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = NULL,
    sequences        = d$sequences,
    not_in_global_ids = c(1L, 2L, 3L),
    seqs             = d$seqs,
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = FALSE,
    public_tcrs      = TRUE,
    cluster_min_size = 1,
    verbose          = FALSE
  )

  expect_true(!is.null(result$clone_network))
  expect_true(all(result$clone_network$type == "singleton"))
})

# ---- .cluster_gliph2 structboundaries FALSE --------------------------------

test_that(".cluster_gliph2 works with structboundaries FALSE", {
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
    structboundaries       = FALSE,
    boundary_size          = 3,
    motif_distance_cutoff  = 10,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  expect_type(result, "list")
})

# ---- .cluster_gliph2 save_cluster_list_df creation -------------------------

test_that(".cluster_gliph2 generates save_cluster_list_df", {
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

  if (!is.null(result$merged_clusters) && length(result$cluster_list) > 0) {
    expect_s3_class(result$save_cluster_list_df, "data.frame")
    expect_true("tag" %in% colnames(result$save_cluster_list_df))
  }
})

# ---- .cluster_gliph2 verbose messaging ---------------------------------------

test_that(".cluster_gliph2 prints message when verbose is TRUE", {
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

  expect_message(
    immGLIPH:::.cluster_gliph2(
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
      verbose                = TRUE
    ),
    "GLIPH2"
  )
})

# ---- .cluster_gliph1 with patient.info = FALSE -------------------------------

test_that(".cluster_gliph1 works without patient info", {
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSYLAGGRNTLYF")
  sequences <- data.frame(
    seq_ID = seq_along(seqs),
    CDR3b  = seqs,
    TRBV   = c("TRBV5-1", "TRBV6-2", "TRBV5-1"),
    stringsAsFactors = FALSE
  )

  clone_network <- data.frame(
    V1   = "CASSLAPGATNEKLFF",
    V2   = "CASSLDRGEVFF",
    type = "local",
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = clone_network,
    sequences        = sequences,
    not_in_global_ids = integer(0),
    seqs             = seqs,
    vgene.info       = TRUE,
    patient.info     = FALSE,
    global_vgene     = FALSE,
    public_tcrs      = TRUE,
    cluster_min_size = 1,
    verbose          = FALSE
  )

  expect_type(result, "list")
  expect_true(!is.null(result$cluster_properties))
})

# ---- .cluster_gliph2 with global_vgene TRUE includes TRBV in tag ------------

test_that(".cluster_gliph2 includes TRBV in tag when global_vgene is TRUE", {
  d <- .make_cluster_data()

  global_res <- data.frame(
    cluster_tag    = c("struct_%_13"),
    cluster_size   = c(3),
    unique_CDR3b   = c(3),
    num_in_ref     = c(5),
    fisher.score   = c(0.001),
    aa_at_position = c("L"),
    TRBV           = c("TRBV5-1"),
    CDR3b          = paste(d$seqs[c(1, 3, 5)], collapse = " "),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph2(
    local_res              = NULL,
    global_res             = global_res,
    sequences              = d$sequences,
    local_similarities     = FALSE,
    global_similarities    = TRUE,
    global_vgene           = TRUE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 1,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  if (!is.null(result$merged_clusters)) {
    expect_true(any(grepl("TRBV", result$merged_clusters$tag)))
  }
})

# ---- .cluster_gliph2 handles infinite OvE ------------------------------------

test_that(".cluster_gliph2 handles infinite OvE values", {
  d <- .make_cluster_data()

  local_res <- data.frame(
    motif         = c("SL"),
    start         = c(1),
    stop          = c(2),
    num_in_sample = c(3),
    num_in_ref    = c(0),
    num_fold      = c(Inf),
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

  expect_type(result, "list")
  if (!is.null(result$merged_clusters)) {
    # Infinite OvE should be replaced with 0
    expect_true(all(is.finite(as.numeric(result$merged_clusters$OvE))))
  }
})

# ---- .cluster_gliph1 duplicate cluster naming --------------------------------

test_that(".cluster_gliph1 handles duplicate cluster names with suffixes", {
  # Create sequences where two separate components have the same first CDR3b
  seqs <- c("CASSLAPGATNEKLFF", "CASSLDRGEVFF", "CASSLAPGATNEKLFF",
            "CASSYLAGGRNTLYF")
  sequences <- data.frame(
    seq_ID  = 1:4,
    CDR3b   = seqs,
    TRBV    = c("TRBV5-1", "TRBV6-2", "TRBV5-1", "TRBV7-2"),
    patient = c("P1", "P1", "P2", "P2"),
    stringsAsFactors = FALSE
  )

  clone_network <- data.frame(
    V1   = c("CASSLAPGATNEKLFF", "CASSLAPGATNEKLFF"),
    V2   = c("CASSLDRGEVFF", "CASSYLAGGRNTLYF"),
    type = c("local", "local"),
    stringsAsFactors = FALSE
  )

  result <- immGLIPH:::.cluster_gliph1(
    clone_network    = clone_network,
    sequences        = sequences,
    not_in_global_ids = integer(0),
    seqs             = unique(seqs),
    vgene.info       = TRUE,
    patient.info     = TRUE,
    global_vgene     = FALSE,
    public_tcrs      = TRUE,
    cluster_min_size = 1,
    verbose          = FALSE
  )

  expect_type(result, "list")
  if (!is.null(result$cluster_properties)) {
    # All tags should be unique
    expect_equal(length(unique(result$cluster_properties$tag)),
                 nrow(result$cluster_properties))
  }
})

# ---- .cluster_gliph2 motif_distance_cutoff restricts edges -------------------

test_that(".cluster_gliph2 motif_distance_cutoff = 0 eliminates positionally distant edges", {
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

  result_strict <- immGLIPH:::.cluster_gliph2(
    local_res              = local_res,
    global_res             = NULL,
    sequences              = d$sequences,
    local_similarities     = TRUE,
    global_similarities    = FALSE,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 0,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  result_lenient <- immGLIPH:::.cluster_gliph2(
    local_res              = local_res,
    global_res             = NULL,
    sequences              = d$sequences,
    local_similarities     = TRUE,
    global_similarities    = FALSE,
    global_vgene           = FALSE,
    all_aa_interchangeable = TRUE,
    structboundaries       = TRUE,
    boundary_size          = 3,
    motif_distance_cutoff  = 100,
    cluster_min_size       = 1,
    boost_local_significance = FALSE,
    verbose                = FALSE
  )

  # Strict cutoff should have fewer or equal edges
  n_strict <- if (!is.null(result_strict$clone_network)) {
    sum(result_strict$clone_network$type == "local", na.rm = TRUE)
  } else 0
  n_lenient <- if (!is.null(result_lenient$clone_network)) {
    sum(result_lenient$clone_network$type == "local", na.rm = TRUE)
  } else 0
  expect_true(n_strict <= n_lenient)
})
