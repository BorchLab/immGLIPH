test_that("clusterScoring validates input types", {
  expect_error(clusterScoring(cluster_list = "not_a_list",
                               cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF")),
               "class 'list'")
  expect_error(clusterScoring(cluster_list = list(),
                               cdr3_sequences = 42),
               "class 'data.frame'")
})

test_that("clusterScoring validates refdb_beta", {
  expect_error(clusterScoring(cluster_list = list(),
                               cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF"),
                               refdb_beta = "nonexistent_ref"),
               "refdb_beta must be")
})

test_that("clusterScoring validates gliph_version", {
  expect_error(clusterScoring(cluster_list = list(),
                               cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF"),
                               gliph_version = 3),
               "gliph_version")
})

test_that("clusterScoring validates sim_depth", {
  expect_error(clusterScoring(cluster_list = list(),
                               cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF"),
                               sim_depth = 0),
               "sim_depth")
})

test_that("clusterScoring validates hla_cutoff", {
  expect_error(clusterScoring(cluster_list = list(),
                               cdr3_sequences = data.frame(CDR3b = "CASSLAPGATNEKLFF"),
                               hla_cutoff = 2),
               "hla_cutoff")
})
