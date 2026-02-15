#' GLIPH1-style clustering via igraph connected components
#'
#' Builds a clone network from local and global edges, optionally filters
#' by V-gene and donor constraints, constructs an igraph graph, and extracts
#' connected components as convergence groups. Each component becomes one
#' cluster named \code{CRG-<first_CDR3b>}.
#'
#' @param clone_network Data frame. Edge list with columns \code{V1},
#'   \code{V2}, \code{type} (and optionally \code{tag}).
#' @param sequences Data frame. Full sample sequences with at least
#'   \code{CDR3b} and optionally \code{TRBV}, \code{patient}.
#' @param not_in_global_ids Integer vector. Indices of sequences without
#'   global neighbours (used to add singletons).
#' @param seqs Character vector. Unique CDR3b sequences.
#' @param vgene.info Logical. Whether V-gene info is available.
#' @param patient.info Logical. Whether patient info is available.
#' @param global_vgene Logical. Whether global edges require V-gene match.
#' @param public_tcrs Logical. If \code{FALSE}, restrict edges to same donor.
#' @param cluster_min_size Integer. Minimum cluster size to retain.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with:
#' \describe{
#'   \item{cluster_properties}{Data frame with cluster_size, tag, members.}
#'   \item{cluster_list}{Named list of data frames (member details).}
#'   \item{clone_network}{Data frame of the final edge list.}
#'   \item{save_cluster_list_df}{Data frame for saving cluster members.}
#' }
#'
#' @import foreach
#' @keywords internal
.cluster_gliph1 <- function(clone_network,
                             sequences,
                             not_in_global_ids,
                             seqs,
                             vgene.info,
                             patient.info,
                             global_vgene,
                             public_tcrs,
                             cluster_min_size,
                             verbose) {

  if (verbose) message("Clustering sequences (GLIPH1.0 method).")

  ## ---- Filter edges by donor constraints ----
  if (is.logical(public_tcrs) && !public_tcrs && patient.info &&
      !is.null(clone_network)) {
    temp <- apply(clone_network, 1, function(x) {
      any(sequences$patient[sequences$CDR3b == x[1]] %in%
            sequences$patient[sequences$CDR3b == x[2]])
    })
    clone_network <- clone_network[temp, ]
    if (nrow(clone_network) == 0) clone_network <- NULL
  }

  ## ---- Add singletons ----
  in_local_ids <- integer(0)
  if (!is.null(clone_network)) {
    in_network <- unique(c(clone_network$V1, clone_network$V2))
    # Track which sequence indices are in any edge
    in_local_ids <- which(seqs %in% in_network)
  }

  singleton_ids <- setdiff(not_in_global_ids, in_local_ids)
  if (length(singleton_ids) > 0 && !is.null(clone_network)) {
    clone_network <- rbind(
      clone_network,
      data.frame(
        V1   = seqs[singleton_ids],
        V2   = seqs[singleton_ids],
        type = rep("singleton", length(singleton_ids)),
        stringsAsFactors = FALSE
      )
    )
  } else if (is.null(clone_network) && length(singleton_ids) > 0) {
    clone_network <- data.frame(
      V1   = seqs[singleton_ids],
      V2   = seqs[singleton_ids],
      type = rep("singleton", length(singleton_ids)),
      stringsAsFactors = FALSE
    )
  }

  ## Empty results
  part_4_res <- NULL
  part_4_res_list <- list()
  save_cluster_list_df <- NULL

  if (is.null(clone_network) || nrow(clone_network) == 0) {
    return(list(
      cluster_properties  = part_4_res,
      cluster_list        = part_4_res_list,
      clone_network       = clone_network,
      save_cluster_list_df = save_cluster_list_df
    ))
  }

  clone_network[] <- lapply(clone_network, as.character)

  ## ---- Build tuple network ----
  x <- NULL
  temp_clone_network <- foreach::foreach(
    x = seq_len(nrow(clone_network))
  ) %dopar% {
    act_ids1 <- which(sequences$CDR3b == clone_network[x, 1])
    act_ids2 <- which(sequences$CDR3b == clone_network[x, 2])

    comb_ids <- expand.grid(act_ids1, act_ids2)
    comb_ids <- unique(comb_ids)

    act_infos1 <- sequences[comb_ids[, 1], ]
    act_infos2 <- sequences[comb_ids[, 2], ]

    ## Filter for same donor if required
    if (is.logical(public_tcrs) && !public_tcrs && patient.info) {
      exclude_rows <- act_infos1$patient != act_infos2$patient
      act_infos1 <- act_infos1[!exclude_rows, ]
      act_infos2 <- act_infos2[!exclude_rows, ]
    }

    ## Filter global edges for matching V-gene
    if (global_vgene && vgene.info &&
        clone_network[x, 3] == "global" && nrow(act_infos1) > 0) {
      exclude_rows <- act_infos1$TRBV != act_infos2$TRBV
      act_infos1 <- act_infos1[!exclude_rows, ]
      act_infos2 <- act_infos2[!exclude_rows, ]
    }

    if (nrow(act_infos1) > 0) {
      act_infos1 <- do.call(paste, c(act_infos1, sep = "$#$#$"))
      act_infos2 <- do.call(paste, c(act_infos2, sep = "$#$#$"))
      var_ret <- data.frame(act_infos1, act_infos2)
      t(var_ret)
    } else {
      NULL
    }
  }

  temp_clone_network <- data.frame(
    matrix(unlist(temp_clone_network), ncol = 2, byrow = TRUE),
    stringsAsFactors = FALSE
  )
  temp_clone_network[] <- lapply(temp_clone_network, as.character)

  if (nrow(temp_clone_network) == 0) {
    return(list(
      cluster_properties  = NULL,
      cluster_list        = list(),
      clone_network       = clone_network,
      save_cluster_list_df = NULL
    ))
  }

  ## ---- igraph clustering ----
  gr <- igraph::graph_from_edgelist(
    as.matrix(temp_clone_network[, c(1, 2)]),
    directed = FALSE
  )
  cm <- igraph::components(gr)

  all_leaders <- character(0)
  for (i in seq_len(cm$no)) {
    members_id <- names(cm$membership)[which(cm$membership == i)]
    members_id <- unique(members_id)
    members_id <- data.frame(
      matrix(unlist(strsplit(members_id, split = "$#$#$", fixed = TRUE)),
             ncol = ncol(sequences), byrow = TRUE),
      stringsAsFactors = FALSE
    )
    colnames(members_id) <- colnames(sequences)

    members <- sort(unique(members_id$CDR3b))
    csize   <- length(members)

    leader <- paste("CRG", members[1], sep = "-")
    while_counter <- 1
    temp_leader <- leader
    while (temp_leader %in% all_leaders) {
      temp_leader <- paste0(leader, "-", while_counter)
      while_counter <- while_counter + 1
    }
    leader <- temp_leader
    all_leaders <- c(all_leaders, leader)

    row_df <- data.frame(
      cluster_size = csize,
      tag          = leader,
      members      = paste(members, collapse = " "),
      stringsAsFactors = FALSE
    )

    if (i == 1) {
      part_4_res <- row_df
    } else {
      part_4_res <- rbind(part_4_res, row_df)
    }
    part_4_res_list[[leader]] <- members_id
  }

  ## ---- Filter by minimum cluster size ----
  eliminate_ids <- which(part_4_res$cluster_size < cluster_min_size)
  if (length(eliminate_ids) > 0) {
    part_4_res <- part_4_res[-eliminate_ids, ]
    part_4_res_list <- part_4_res_list[-eliminate_ids]
  }

  if (nrow(part_4_res) > 0) {
    save_cluster_list_df <- foreach::foreach(
      i = seq_along(part_4_res_list),
      .combine = "rbind"
    ) %dopar% {
      temp <- part_4_res_list[[i]]
      cbind(
        data.frame(tag = rep(names(part_4_res_list)[i], nrow(temp))),
        temp
      )
    }
  } else {
    part_4_res <- NULL
    part_4_res_list <- list()
    save_cluster_list_df <- NULL
  }

  list(
    cluster_properties   = part_4_res,
    cluster_list         = part_4_res_list,
    clone_network        = clone_network,
    save_cluster_list_df = save_cluster_list_df
  )
}


#' GLIPH2-style clustering by individual motif/struct tags
#'
#' Each local motif or global struct tag defines its own cluster. Edges are
#' restricted by motif distance for local connections and by BLOSUM62 for
#' global connections. Clusters are built as igraph components within each
#' tag, and enrichment is reassessed per cluster.
#'
#' @param local_res Data frame. Selected motifs from local enrichment
#'   (must have columns \code{motif}, \code{start}, \code{stop},
#'   \code{num_in_sample}, \code{num_in_ref}, \code{num_fold},
#'   \code{fisher.score}, \code{members}).
#' @param global_res Data frame. Global cluster list from
#'   \code{.global_fisher()} (columns: \code{cluster_tag},
#'   \code{cluster_size}, \code{unique_CDR3b}, \code{num_in_ref},
#'   \code{fisher.score}, \code{aa_at_position}, \code{TRBV},
#'   \code{CDR3b}).
#' @param sequences Data frame. Full sample data.
#' @param local_similarities Logical. Whether local similarities were run.
#' @param global_similarities Logical. Whether global similarities were found.
#' @param global_vgene Logical. Restrict global edges to matching V-gene.
#' @param all_aa_interchangeable Logical. BLOSUM62 filtering.
#' @param structboundaries Logical. Whether boundary trimming is active.
#' @param boundary_size Integer. Boundary trim size.
#' @param motif_distance_cutoff Integer. Max positional distance for
#'   local motifs.
#' @param cluster_min_size Integer. Minimum cluster size.
#' @param boost_local_significance Logical. Whether to boost local p-values
#'   using germline N-nucleotide information.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with:
#' \describe{
#'   \item{merged_clusters}{Data frame of cluster properties (type, tag,
#'     cluster_size, unique_cdr3_sample, unique_cdr3_ref, OvE, fisher.score,
#'     members).}
#'   \item{cluster_list}{Named list of data frames with member details.}
#'   \item{clone_network}{Data frame of edges (V1, V2, type, cluster_tag).}
#'   \item{save_cluster_list_df}{Data frame for saving cluster member details.}
#' }
#'
#' @import foreach
#' @keywords internal
.cluster_gliph2 <- function(local_res,
                             global_res,
                             sequences,
                             local_similarities,
                             global_similarities,
                             global_vgene,
                             all_aa_interchangeable,
                             structboundaries,
                             boundary_size,
                             motif_distance_cutoff,
                             cluster_min_size,
                             boost_local_significance,
                             verbose) {

  if (verbose) message("Clustering sequences (GLIPH2.0 method).")

  ## ---- N-nucleotide identification for germline boosting ----
  range_df <- NULL
  if (boost_local_significance) {
    gTRB <- NULL
    utils::data("gTRB", package = "immGLIPH", envir = environment())

    range_df <- data.frame(
      start_1 = rep(0, nrow(sequences)),
      stop_1  = rep(0, nrow(sequences)),
      start_2 = rep(0, nrow(sequences)),
      stop_2  = rep(0, nrow(sequences))
    )

    ## V gene fragments
    remaining_ids <- seq_len(nrow(sequences))
    for (i in max(gTRB$gTRBV[, 2]):min(gTRB$gTRBV[, 2])) {
      temp_bool <- substr(sequences$CDR3b[remaining_ids], 1, i) %in%
        gTRB$gTRBV[, 1][gTRB$gTRBV[, 2] == i]
      range_df$start_1[remaining_ids][temp_bool] <- i + 1
      remaining_ids <- remaining_ids[!temp_bool]
      if (length(remaining_ids) == 0) break
    }

    ## J gene fragments
    remaining_ids <- seq_len(nrow(sequences))
    for (i in max(gTRB$gTRBJ[, 2]):min(gTRB$gTRBJ[, 2])) {
      temp_bool <- substr(
        sequences$CDR3b[remaining_ids],
        nchar(sequences$CDR3b[remaining_ids]) - i + 1,
        nchar(sequences$CDR3b[remaining_ids])
      ) %in% gTRB$gTRBJ[, 1][gTRB$gTRBJ[, 2] == i]
      range_df$stop_2[remaining_ids][temp_bool] <- i
      remaining_ids <- remaining_ids[!temp_bool]
      if (length(remaining_ids) == 0) break
    }
    range_df$stop_2 <- nchar(sequences$CDR3b) - range_df$stop_2

    ## D gene fragments
    remaining_ids <- seq_len(nrow(sequences))
    if (length(remaining_ids) > 0) {
      act_subseqs <- substr(
        sequences$CDR3b[remaining_ids],
        range_df$start_1[remaining_ids],
        range_df$stop_2[remaining_ids]
      )
      for (i in max(gTRB$gTRBD[, 2]):min(gTRB$gTRBD[, 2])) {
        max_len <- max(nchar(act_subseqs))
        if ((max_len - i + 1) >= 1) {
          act_subseqs_matrix <- NULL
          for (j in seq_len(max_len - i + 1)) {
            col_data <- substr(act_subseqs, j, j + i - 1)
            if (j == 1) {
              act_subseqs_matrix <- matrix(col_data, ncol = 1)
            } else {
              act_subseqs_matrix <- cbind(act_subseqs_matrix, col_data)
            }
          }
          for (j in gTRB$gTRBD[gTRB$gTRBD[, 2] == i, 1]) {
            act_values <- do.call(
              pmax,
              data.frame(t(t(act_subseqs_matrix == j) *
                             (ncol(act_subseqs_matrix):1)))
            )
            act_values[act_values > 0] <-
              ncol(act_subseqs_matrix) - act_values[act_values > 0] + 1

            range_df$start_1[remaining_ids][act_values == 1] <-
              range_df$start_1[remaining_ids][act_values == 1] + i
            range_df$stop_1[remaining_ids][act_values > 1] <-
              range_df$start_1[remaining_ids][act_values > 1] +
              act_values[act_values > 1] - 2
            range_df$start_2[remaining_ids][act_values > 1] <-
              range_df$stop_1[remaining_ids][act_values > 1] + i + 1
            temp_bool <- (act_values > 0)

            remaining_ids <- remaining_ids[!temp_bool]
            act_subseqs <- act_subseqs[!temp_bool]
            if (is.matrix(act_subseqs_matrix)) {
              act_subseqs_matrix <- act_subseqs_matrix[!temp_bool, ,
                                                        drop = FALSE]
            }
            if (length(remaining_ids) == 0) break
          }
          if (length(remaining_ids) == 0) break
        }
      }
    }
  }

  ## ---- Mark non-germline positions in sequence ----
  sequences$ultCDR3b <- rep("", nrow(sequences))
  if (boost_local_significance) {
    sequences$ultCDR3b <- toupper(sequences$CDR3b)
    substr(sequences$ultCDR3b, range_df$start_1, range_df$stop_2) <-
      tolower(substr(sequences$ultCDR3b, range_df$start_1, range_df$stop_2))
    substr(
      sequences$ultCDR3b,
      range_df$stop_1 + 1,
      range_df$start_2 - 1
    ) <- toupper(substr(
      sequences$ultCDR3b,
      range_df$stop_1 + 1,
      range_df$start_2 - 1
    ))
  }

  ## ---- Combine local and global cluster information ----
  merged_clusters <- NULL

  if (local_similarities && !is.null(local_res) && nrow(local_res) > 0) {
    local_clusters <- data.frame(
      type              = rep("local", nrow(local_res)),
      tag               = paste(local_res$motif, local_res$start,
                                local_res$stop, sep = "_"),
      cluster_size      = rep(0, nrow(local_res)),
      unique_cdr3_sample = local_res$num_in_sample,
      unique_cdr3_ref   = local_res$num_in_ref,
      OvE               = local_res$num_fold,
      fisher.score      = local_res$fisher.score,
      members           = local_res$members,
      stringsAsFactors  = FALSE
    )
    merged_clusters <- local_clusters
  }

  if (global_similarities && !is.null(global_res) && nrow(global_res) > 0) {
    if (global_vgene) {
      g_tag <- paste(global_res$cluster_tag, global_res$TRBV,
                     global_res$aa_at_position, sep = "_")
    } else {
      g_tag <- paste(global_res$cluster_tag,
                     global_res$aa_at_position, sep = "_")
    }

    global_clusters <- data.frame(
      type              = rep("global", nrow(global_res)),
      tag               = g_tag,
      cluster_size      = global_res$cluster_size,
      unique_cdr3_sample = global_res$unique_CDR3b,
      unique_cdr3_ref   = global_res$num_in_ref,
      OvE               = rep(0, nrow(global_res)),
      fisher.score      = global_res$fisher.score,
      members           = global_res$CDR3b,
      stringsAsFactors  = FALSE
    )

    if (!is.null(merged_clusters)) {
      merged_clusters <- rbind(merged_clusters, global_clusters)
    } else {
      merged_clusters <- global_clusters
    }
  }

  ## ---- Build cluster list ----
  cluster_list <- list()
  clone_network <- NULL

  if (!is.null(merged_clusters) && nrow(merged_clusters) > 0) {
    ## Load BLOSUM vector for global filtering
    BlosumVec <- .get_blosum_vec()

    i <- NULL
    cluster_list <- foreach::foreach(
      i = seq_len(nrow(merged_clusters))
    ) %dopar% {
      if (merged_clusters$type[i] == "local") {
        act_seqs <- unlist(strsplit(merged_clusters$members[i], split = " "))
        return(sequences[sequences$CDR3b %in% act_seqs, ])
      }
      if (merged_clusters$type[i] == "global") {
        act_seqs <- unlist(strsplit(merged_clusters$members[i], split = " "))
        act_details <- unlist(strsplit(merged_clusters$tag[i], split = "_"))
        if (!global_vgene) {
          return(data.frame(
            sequences[sequences$CDR3b %in% act_seqs, ],
            stringsAsFactors = FALSE
          ))
        } else {
          return(data.frame(
            sequences[sequences$CDR3b %in% act_seqs &
                        sequences$TRBV == act_details[2], ],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    names(cluster_list) <- merged_clusters$tag

    ## Update local cluster sizes
    if (local_similarities && !is.null(local_res) && nrow(local_res) > 0) {
      for (i in which(merged_clusters$type == "local")) {
        merged_clusters$cluster_size[i] <-
          length(unique(cluster_list[[i]]$CDR3b))

        if (boost_local_significance && !is.null(range_df)) {
          motif_name <- strsplit(names(cluster_list)[i], split = "_")[[1]][[1]]
          start_motif <- as.numeric(
            strsplit(names(cluster_list)[i], split = "_")[[1]][[2]]
          )
          stop_motif <- start_motif + nchar(motif_name) - 1
          act_ids <- cluster_list[[i]]$seq_ID
          motif_in_range <-
            (range_df$start_1[act_ids] <= stop_motif &
               range_df$stop_1[act_ids] >= start_motif) |
            (range_df$start_2[act_ids] <= stop_motif &
               range_df$stop_2[act_ids] >= start_motif)
          merged_clusters$fisher.score[i] <-
            merged_clusters$fisher.score[i] / (2^sum(motif_in_range))
        }
      }
    }

    ## Minimum size filter
    eliminate_ids <- which(merged_clusters$cluster_size < cluster_min_size)
    if (length(eliminate_ids) > 0) {
      merged_clusters <- merged_clusters[-eliminate_ids, ]
      cluster_list <- cluster_list[-eliminate_ids]
    }
    if (nrow(merged_clusters) == 0) {
      merged_clusters <- NULL
      cluster_list <- list()
    }
  }

  ## ---- Build clone network edges ----
  if (!is.null(merged_clusters) && length(cluster_list) > 0) {
    ## Local edges
    if (local_similarities && any(merged_clusters$type == "local")) {
      local_clone_network <- foreach::foreach(
        i = which(merged_clusters$type == "local")
      ) %dopar% {
        temp_members <- cluster_list[[i]]$CDR3b
        if (structboundaries) {
          temp_members_frags <- substr(
            temp_members,
            boundary_size + 1,
            nchar(temp_members) - boundary_size
          )
        } else {
          temp_members_frags <- temp_members
        }

        motif_name <- strsplit(names(cluster_list)[i], split = "_")[[1]][[1]]
        temp_pos <- stringr::str_locate(temp_members_frags, motif_name)

        if (length(temp_members) >= 2) {
          combn_ids <- t(utils::combn(seq_along(temp_members), m = 2))
        } else {
          combn_ids <- t(utils::combn(rep(1, 2), m = 2))
        }
        temp_df <- data.frame(
          V1          = temp_members[combn_ids[, 1]],
          V2          = temp_members[combn_ids[, 2]],
          type        = rep("local", nrow(combn_ids)),
          cluster_tag = rep(names(cluster_list)[i], nrow(combn_ids)),
          stringsAsFactors = FALSE
        )

        ## Restrict by positional distance
        temp_df <- unique(temp_df[
          abs(temp_pos[combn_ids[, 1], 1] -
                temp_pos[combn_ids[, 2], 1]) < motif_distance_cutoff, ])

        t(temp_df)
      }
      local_clone_network <- data.frame(
        matrix(unlist(local_clone_network), ncol = 4, byrow = TRUE),
        stringsAsFactors = FALSE
      )
      colnames(local_clone_network) <- c("V1", "V2", "type", "cluster_tag")
      clone_network <- local_clone_network
    }

    ## Global edges
    if (global_similarities && any(merged_clusters$type == "global")) {
      global_clone_network <- foreach::foreach(
        i = which(merged_clusters$type == "global")
      ) %dopar% {
        temp_members <- cluster_list[[i]]$CDR3b

        if (length(temp_members) >= 2) {
          combn_ids <- t(utils::combn(seq_along(temp_members), m = 2))
        } else {
          combn_ids <- t(utils::combn(rep(1, 2), m = 2))
        }
        temp_df <- data.frame(
          V1          = temp_members[combn_ids[, 1]],
          V2          = temp_members[combn_ids[, 2]],
          type        = rep("global", nrow(combn_ids)),
          cluster_tag = rep(names(cluster_list)[i], nrow(combn_ids)),
          stringsAsFactors = FALSE
        )

        ## BLOSUM62 filtering at variable position
        if (!all_aa_interchangeable) {
          tag_name <- names(cluster_list)[i]
          temp_pos <- stringr::str_locate(tag_name, "%")
          if (structboundaries) temp_pos <- temp_pos + boundary_size
          temp_df <- unique(temp_df[
            paste0(
              substr(temp_df$V1, temp_pos[1], temp_pos[1]),
              substr(temp_df$V2, temp_pos[1], temp_pos[1])
            ) %in% BlosumVec, ])
        }

        t(temp_df)
      }
      global_clone_network <- data.frame(
        matrix(unlist(global_clone_network), ncol = 4, byrow = TRUE),
        stringsAsFactors = FALSE
      )
      colnames(global_clone_network) <- c("V1", "V2", "type", "cluster_tag")

      if (is.null(clone_network)) {
        clone_network <- global_clone_network
      } else {
        clone_network <- rbind(clone_network, global_clone_network)
      }
    }

    clone_network[] <- lapply(clone_network, as.character)

    ## Add singletons
    not_in_network <- sequences$CDR3b[
      !(sequences$CDR3b %in% c(clone_network$V1, clone_network$V2))
    ]
    if (length(not_in_network) > 0) {
      clone_network <- rbind(
        clone_network,
        data.frame(
          V1          = not_in_network,
          V2          = not_in_network,
          type        = rep("singleton", length(not_in_network)),
          cluster_tag = paste0("singleton_", seq_along(not_in_network)),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  ## ---- Save cluster list as data frame ----
  save_cluster_list_df <- NULL
  if (!is.null(merged_clusters) && length(cluster_list) > 0) {
    merged_clusters$fisher.score <- formatC(
      merged_clusters$fisher.score, digits = 1, format = "e"
    )
    merged_clusters$OvE[is.infinite(merged_clusters$OvE)] <- 0

    save_cluster_list_df <- foreach::foreach(
      i = seq_along(cluster_list),
      .combine = "rbind"
    ) %dopar% {
      temp <- cluster_list[[i]]
      cbind(
        data.frame(tag = rep(names(cluster_list)[i], nrow(temp))),
        temp
      )
    }
  }

  list(
    merged_clusters      = merged_clusters,
    cluster_list         = cluster_list,
    clone_network        = clone_network,
    save_cluster_list_df = save_cluster_list_df
  )
}
