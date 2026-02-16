#' Visualize TCR convergence group network
#'
#' Uses the visNetwork package to build an interactive network graph from the
#' clustering results produced by \code{\link{runGLIPH}}. Nodes represent 
#' individual CDR3b sequences and edges encode local or global sequence 
#' similarities. The resulting visualization is fully interactive: scroll to 
#' zoom, hover over a node for details, and click a node to highlight its 
#' direct neighbors.
#'
#' @param clustering_output Output list returned by \code{\link{runGLIPH}}.
#' **Default:** `NULL`
#' @param result_folder Path to the folder containing saved GLIPH output
#' files. When a non-empty path is supplied the results are loaded from disk
#' and \code{clustering_output} is ignored. **Default:** `""`
#' @param show_additional_columns Character vector of extra column names whose
#' values should be displayed in the node tooltips. Column names from the
#' original \code{cdr3_sequences} data frame and from
#' \code{clustering_output$cluster_properties} are accepted.
#' **Default:** `NULL`
#' @param color_info Column name used to colour the nodes. Accepts any column
#' from the input \code{cdr3_sequences} or
#' \code{clustering_output$cluster_properties}. Set to \code{"none"} to
#' colour all nodes grey, or \code{"color"} to use pre-assigned colour
#' values stored in that column. For numeric columns the viridis palette is
#' applied automatically (purple = low, yellow = high).
#' **Default:** `"total.score"`
#' @param color_palette A function that accepts a single integer \code{n} and
#' returns \code{n} colour values. **Default:** `viridis::viridis`
#' @param local_edge_color Colour applied to edges representing local
#' similarities. **Default:** `"orange"`
#' @param global_edge_color Colour applied to edges representing global
#' similarities. **Default:** `"#68bceb"`
#' @param size_info Column name whose numeric values determine node sizes.
#' Accepts columns from \code{cdr3_sequences} or
#' \code{clustering_output$cluster_properties}. **Default:** `NULL`
#' @param absolute_size If \code{TRUE} the raw values from the
#' \code{size_info} column are used as node sizes; otherwise the values are
#' linearly scaled to the range 12--20. **Default:** `FALSE`
#' @param cluster_min_size Minimum number of members a cluster must contain to
#' be included in the plot. **Default:** `3`
#' @param n_cores Number of cores for parallel processing. When \code{NULL}
#' the number of available cores minus two is used. **Default:** `1`
#'
#' @return A `visNetwork` object containing the interactive network graph.
#'
#' @examples
#' \dontrun{
#' utils::data("gliph_input_data")
#' res <- runGLIPH(cdr3_sequences = gliph_input_data[seq_len(200),],
#'                 method = "gliph1",
#'                 sim_depth = 100,
#'                 n_cores = 1)
#'
#' plotNetwork(clustering_output = res,
#'             n_cores = 1)
#' }
#'
#' @import viridis foreach grDevices
#' @export

plotNetwork <- function(clustering_output = NULL,
                         result_folder = "",
                         show_additional_columns = NULL,
                         color_info = "total.score",
                         color_palette = viridis::viridis,
                         local_edge_color = "orange",
                         global_edge_color = "#68bceb",
                         size_info = NULL,
                         absolute_size = FALSE,
                         cluster_min_size = 3,
                         n_cores = 1) {
  clustering_output <- clustering_output

  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################

  ### result_folder clustering_output
  if(!is.character(result_folder)) stop("result_folder has to be a character object")
  if(length(result_folder) > 1) stop("result_folder has to be a single path")
  if(result_folder != ""){
    clustering_output <- loadGLIPH(result_folder = result_folder)
    message()
  }

  ### show_additional_columns
  if(!is.character(show_additional_columns) && !(is.null(show_additional_columns))) stop("show_additional_columns must be a character vector representing column names.")

  ### color_info
  if(!is.character(color_info)) stop("color_info must be a character.")
  if(length(color_info) > 1) stop("color_info has to be a single column name")

  ### color_palette
  if(!is.function(color_palette)) stop("color_palette must be a function expecting a number as input and returning color values.")
  if(!plotfunctions::isColor(color_palette(1))) stop("color_palette must be a function expecting a number as input and returning color values.")

  ### local_edge_color
  if(!plotfunctions::isColor(local_edge_color)) stop("local_edge_color has to be a color.")

  ### global_edge_color
  if(!plotfunctions::isColor(global_edge_color)) stop("global_edge_color has to be a color.")

  ### size_info
  if(!is.character(size_info) && !(is.null(size_info))) stop("size_info must be either NULL or a character.")
  if(length(size_info) > 1  && !(is.null(size_info))) stop("size_info has to be a single column name")

  ### absolute_size
  if(!is.logical(absolute_size)) stop("absolute_size has to be logical")

  ### cluster_min_size
  if(!is.numeric(cluster_min_size)) stop("cluster_min_size has to be numeric")
  if(length(cluster_min_size) > 1) stop("cluster_min_size has to be a single number")
  if(cluster_min_size < 1) stop("cluster_min_size must be at least 1")
  cluster_min_size <- round(cluster_min_size)

  ### n_cores
  if(is.null(n_cores)) n_cores <- max(1, parallel::detectCores()-2) else {
    if(!is.numeric(n_cores)) stop("n_cores has to be numeric")
    if(length(n_cores) > 1) stop("n_cores has to be a single number")
    if(n_cores < 1) stop("n_cores must be at least 1")

    n_cores <- max(1, min(n_cores, parallel::detectCores()-2))
  }

  #################################################################
  ##                      Input preparation                      ##
  #################################################################

  ### Initiate parallelization
  .setup_parallel(n_cores)

  ### Get cluster_list and cluster_properties with at least cluster_min_size members
  # cluster_list:       contains all members of the cluster and the sequence specific additional information (e.g. patient, count, etc.)
  # cluster_properties: contains cluster specific information like all scores
  parameters <- clustering_output$parameters
  cluster_list <- clustering_output$cluster_list
  if(is.null(cluster_list)) stop("The specified clustering_output does not contain any clusters.")
  if(length(cluster_list) == 0) stop("The specified clustering_output does not contain any clusters.")
  cluster_properties <- clustering_output$cluster_properties

  hold_ids <- which(as.numeric(cluster_properties$cluster_size) >= cluster_min_size)
  if(length(hold_ids) == 0) stop("The specified clustering_output does not contain any clusters with a minimal cluster size of ",
                                            cluster_min_size,
                                            ".")
  cluster_list <- cluster_list[hold_ids]
  cluster_properties <- cluster_properties[hold_ids,]

  ### Pool sequence and cluster specific information (identical sequences in differenct clusters are treated as different entities)
  cluster_data_frame <- data.frame(foreach::foreach(i = seq_along(cluster_list), .combine = rbind) %dopar% {
      temp_df <- cluster_list[[i]]
      temp_adds <- unlist(cluster_properties[i,])
      temp_adds <- data.frame(matrix(rep(temp_adds, times = nrow(temp_df)), nrow = nrow(temp_df), byrow = TRUE), stringsAsFactors = FALSE)
      colnames(temp_adds) <- colnames(cluster_properties)
      temp_df <- cbind(temp_adds, temp_df)

      return(temp_df)
  }, stringsAsFactors = FALSE)
  cluster_data_frame[] <- lapply(cluster_data_frame, as.character)
  for(i in seq_len(ncol(cluster_data_frame))){
    if(suppressWarnings(any(is.na(as.numeric(cluster_data_frame[,i])))) == FALSE) cluster_data_frame[,i] <- as.numeric(cluster_data_frame[,i])
  }
  cluster_data_frame$ID <- seq_len(nrow(cluster_data_frame))

  # For better visualization insert a line break (<br>) between significant HLA alleles in the same cluster
  if("lowest.hlas" %in% colnames(cluster_data_frame)){
    cluster_data_frame$lowest.hlas <- gsub(", ", "<br>", cluster_data_frame$lowest.hlas)
    cluster_data_frame$lowest.hlas[cluster_data_frame$lowest.hlas != ""] <- paste0("<br>",
                                                                                   cluster_data_frame$lowest.hlas[cluster_data_frame$lowest.hlas != ""])
  }


  #################################################################
  ##                     Prepare graph edges                     ##
  #################################################################

  message("Preparing graph nodes and edges.")
  clone_network <- c()
  ### Different edge identification procedure in distinct GLIPH versions
  if(!("gliph_version" %in% names(parameters))){
    ### Merged function for GLIPH1.0 and GLIPH2.0 options

    ### Clustering method as in GLIPH1.0
    if(parameters$clustering_method == "GLIPH1.0"){
      # Get connections between sequences calculated by GLIPH; exclude all singletons
      ori_clone_net <- clustering_output$connections
      ori_clone_net <- ori_clone_net[ori_clone_net[,3] != "singleton",]
      if(nrow(ori_clone_net) == 0) stop("The specified clustering_output does not contain any clusters.")

      # Are patient and v gene dependent restrictions possible?
      patient.info <- FALSE
      vgene.info <- FALSE
      if("patient" %in% colnames(cluster_data_frame)) patient.info <- TRUE
      if("TRBV" %in% colnames(cluster_data_frame)) vgene.info <- TRUE

      ### Get connections between different tuples consisting of CDR3b sequence, v gene and donor
      x <- NULL
      clone_network <- foreach::foreach(x = seq_len(nrow(ori_clone_net))) %dopar% {
        # Identify all tuples with correpsonding CDR3b sequence
        act_ids1 <- cluster_data_frame$ID[cluster_data_frame$CDR3b == ori_clone_net[x,1]]
        act_ids2 <- cluster_data_frame$ID[cluster_data_frame$CDR3b == ori_clone_net[x,2]]

        # First assume connection between all tuples
        comb_ids <- expand.grid(act_ids1, act_ids2)
        comb_ids <- unique(comb_ids)

        # If requested (public_tcrs != "all"), exclude connections between tuples including different donors
        if((parameters$public_tcrs == "none" || parameters$public_tcrs == ori_clone_net[x,3]) && patient.info == TRUE && nrow(comb_ids) > 0){
          comb_ids <- comb_ids[cluster_data_frame$patient[comb_ids[,1]] == cluster_data_frame$patient[comb_ids[,2]],]
        }

        # If requested (vgene_match != "none"), exclude connections between tuples including different v genes
        if((parameters$vgene_match == "all" || parameters$vgene_match == ori_clone_net[x,3]) && vgene.info == TRUE && nrow(comb_ids) > 0){
          comb_ids <- comb_ids[cluster_data_frame$TRBV[comb_ids[,1]] == cluster_data_frame$TRBV[comb_ids[,2]],]
        }

        # add type of connection and return connections
        if(nrow(comb_ids) > 0){
          var_ret <- data.frame(comb_ids, stringsAsFactors = FALSE)
          var_ret <- cbind(var_ret, rep(ori_clone_net[x,3], nrow(var_ret)))
          return(unlist(t(var_ret)))
          t(var_ret)
        } else return(NULL)
      }
      clone_network <- data.frame(matrix(unlist(clone_network), ncol = 3, byrow = TRUE), stringsAsFactors = FALSE)
      colnames(clone_network) <- c("from", "to", "label")
      clone_network$from <- as.numeric(clone_network$from)
      clone_network$to <- as.numeric(clone_network$to)
    }
    ### Clustering method as in GLIPH2.0
    if(parameters$clustering_method == "GLIPH2.0"){
      ### local connections
      if(parameters$local_similarities == TRUE){
        ### Get local connections between different tuples consisting of CDR3b sequence, v gene and donor
        local_clone_network <- foreach::foreach(i = which(cluster_properties$type == "local")) %dopar% {

          # Get IDs, members and details of current cluster
          temp_ids <- cluster_data_frame$ID[cluster_data_frame$tag == cluster_properties$tag[i]]
          temp_members <- cluster_data_frame$CDR3b[temp_ids]

          # Extract current motif and starting position range for the motif
          if(parameters$structboundaries) temp_members_frags <- substr(x = temp_members,start = parameters$boundary_size + 1 ,stop = nchar(temp_members) - parameters$boundary_size) else temp_members_frags <- temp_members
          pattern <- strsplit(cluster_properties$tag[i], split = "_")[[1]][[1]]

          # Search all positions of the motif in the cluster members
          temp_members_frags <- rep(temp_members_frags, each = max(nchar(temp_members_frags)-nchar(pattern)+1))
          temp_subs <- substr(temp_members_frags,
                                    start = rep(seq_len(max(nchar(temp_members_frags)-nchar(pattern)+1)), times = length(temp_members)),
                                    stop = rep(nchar(pattern):max(nchar(temp_members_frags)), times = length(temp_members)))
          temp_which <- which(temp_subs == pattern)
          details <- data.frame(id = temp_ids[floor((temp_which-1)/max(nchar(temp_members_frags)-nchar(pattern)+1))+1],
                                      pos = ((temp_which-1) %% max(nchar(temp_members_frags)-nchar(pattern)+1))+1,
                                      stringsAsFactors = FALSE)

          # First assume connections between all cluster members
          if(length(temp_members) >= 2){
            combn_ids <- t(utils::combn(seq_along(temp_members), m = 2))
          } else {
            combn_ids <- t(utils::combn(rep(1, 2), m = 2))
          }
          combn_ids <- combn_ids[details$id[combn_ids[,1]] != details$id[combn_ids[,2]],]
          if(!is.matrix(combn_ids)) combn_ids <- matrix(combn_ids, ncol = 2, byrow = TRUE)

          temp_df <- data.frame(from = details$id[combn_ids[,1]], to = details$id[combn_ids[,2]],
                                      label = rep("local", nrow(combn_ids)), stringsAsFactors = FALSE)

          # Restrict local connections to sequences in which the starting position varies in a certain range of positions
          temp_df <- unique(temp_df[abs(details$pos[combn_ids[,1]]-details$pos[combn_ids[,2]]) < parameters$motif_distance_cutoff,])
          temp_df <- temp_df[cluster_data_frame$tag[temp_df$from] == cluster_data_frame$tag[temp_df$to],]
          temp_df <- t(temp_df)
          return(unlist(temp_df))
        }

        ### If available, set clone network to the local network
        if(length(local_clone_network) == 0) parameters$local_similarities == FALSE else {
          local_clone_network <- data.frame(matrix(unlist(local_clone_network), ncol = 3, byrow = TRUE),
                                                  stringsAsFactors = FALSE)
          colnames(local_clone_network) <- c("from", "to", "label")

          clone_network <- local_clone_network
        }

        message("Restrict local similarities to sequences with motif position varying between ", parameters$motif_distance_cutoff," positions.")
      }

      ### global connections
      if(parameters$global_similarities == TRUE){
        # load BlosumVec to determine interchangeable amino acids defined by BLOSUM62
        BlosumVec <- .get_blosum_vec()

        ### Get local connections between different tuples consisting of CDR3b sequence, v gene and donor
        global_clone_network <- foreach::foreach(i = which(cluster_properties$type == "global")) %dopar% {

          # Get IDs and members of current cluster
          temp_ids <- cluster_data_frame$ID[cluster_data_frame$tag == cluster_properties$tag[i]]
          temp_members <- cluster_data_frame$CDR3b[temp_ids]

          # Assume connections between all cluster members
          if(length(temp_members) >= 2){
            combn_ids <- t(utils::combn(seq_along(temp_members), m = 2))
          } else {
            combn_ids <- t(utils::combn(rep(1, 2), m = 2))
          }
          temp_df <- data.frame(from = temp_ids[combn_ids[,1]], to = temp_ids[combn_ids[,2]],
                                      label = rep("global", nrow(combn_ids)), stringsAsFactors = FALSE)

          # If requested, restrict global connections to sequences with a Hamming distance of 1 and a differing amino acid pair allowed by the BLOSUM62 matrix
          if(parameters$all_aa_interchangeable == FALSE){
            temp_pos <- stringr::str_locate(string = cluster_properties$tag[i], pattern = "%")
            if(parameters$structboundaries == TRUE) temp_pos <- temp_pos+parameters$boundary_size
            temp_df <- unique(temp_df[paste0(substr(cluster_data_frame$CDR3b[temp_df$from], temp_pos, temp_pos),
                                                         substr(cluster_data_frame$CDR3b[temp_df$to], temp_pos, temp_pos)) %in% BlosumVec,])
          }
          temp_df <- t(temp_df)
          return(unlist(temp_df))
        }
        if(parameters$all_aa_interchangeable == FALSE) message("Restrict global connections to sequences with a BLOSUM62 value for the amino acid substitution greater or equal to zero.")

        ### Append global network to the clone network
        if(length(global_clone_network) == 0) parameters$local_similarities == FALSE else {
          global_clone_network <- data.frame(matrix(unlist(global_clone_network), ncol = 3, byrow = TRUE),
                                                   stringsAsFactors = FALSE)
          colnames(global_clone_network) <- c("from", "to", "label")

          if(parameters$local_similarities == FALSE){
            clone_network <- global_clone_network
          } else {
            clone_network <- rbind(clone_network, global_clone_network)
          }
        }
      }
      clone_network$from <- as.numeric(clone_network$from)
      clone_network$to <- as.numeric(clone_network$to)
    }
  } else {
    if(parameters$gliph_version == 1){
      ### Original GLIPH
      # Get connections between sequences calculated by GLIPH; exclude all singletons
      ori_clone_net <- clustering_output$connections
      ori_clone_net <- ori_clone_net[ori_clone_net[,3] != "singleton",]
      if(nrow(ori_clone_net) == 0) stop("The specified clustering_output does not contain any clusters.")

      # Are patient and v gene dependent restrictions possible?
      patient.info <- FALSE
      vgene.info <- FALSE
      if("patient" %in% colnames(cluster_data_frame)) patient.info <- TRUE
      if("TRBV" %in% colnames(cluster_data_frame)) vgene.info <- TRUE

      ### Get connections between different tuples consisting of CDR3b sequence, v gene and donor
      x <- NULL
      clone_network <- foreach::foreach(x = seq_len(nrow(ori_clone_net))) %dopar% {
        # Identify all tuples with correpsonding CDR3b sequence
        act_ids1 <- cluster_data_frame$ID[cluster_data_frame$CDR3b == ori_clone_net[x,1]]
        act_ids2 <- cluster_data_frame$ID[cluster_data_frame$CDR3b == ori_clone_net[x,2]]

        # First assume connection between all tuples
        comb_ids <- expand.grid(act_ids1, act_ids2)
        comb_ids <- unique(comb_ids)

        # If requested (public_tcrs == FALSE), exclude connections between tuples including different donors
        if(parameters$public_tcrs == FALSE && patient.info == TRUE && nrow(comb_ids) > 0){
          comb_ids <- comb_ids[cluster_data_frame$patient[comb_ids[,1]] == cluster_data_frame$patient[comb_ids[,2]],]
        }

        # If requested (global_vgene == TRUE), exclude global connections between tuples including different v genes
        if(parameters$global_vgene == TRUE && vgene.info == TRUE && ori_clone_net[x,3] == "global" && nrow(comb_ids) > 0){
          comb_ids <- comb_ids[cluster_data_frame$TRBV[comb_ids[,1]] == cluster_data_frame$TRBV[comb_ids[,2]],]
        }

        # add type of connection and return connections
        if(nrow(comb_ids) > 0){
          var_ret <- data.frame(comb_ids, stringsAsFactors = FALSE)
          var_ret <- cbind(var_ret, rep(ori_clone_net[x,3], nrow(var_ret)))
          return(unlist(t(var_ret)))
          t(var_ret)
        } else return(NULL)
      }
      if(parameters$positional_motifs == TRUE) message("Restrict local similarities to sequences with identical N-terminal motif position.")
      if(parameters$public_tcrs == FALSE && patient.info == TRUE) message("Restrict similarities to sequences obtained from identical donor.")
      if(parameters$global_vgene == TRUE && vgene.info == TRUE) message("Restrict global similarities to sequences with identical v gene.")
      clone_network <- data.frame(matrix(unlist(clone_network), ncol = 3, byrow = TRUE), stringsAsFactors = FALSE)
      colnames(clone_network) <- c("from", "to", "label")
      clone_network$from <- as.numeric(clone_network$from)
      clone_network$to <- as.numeric(clone_network$to)
    }
    if(parameters$gliph_version == 2){
      ### GLIPH2

      ### local connections
      if(parameters$local_similarities == TRUE){
        ### Get local connections between different tuples consisting of CDR3b sequence, v gene and donor
        local_clone_network <- foreach::foreach(i = which(cluster_properties$type == "local")) %dopar% {

          # Get IDs, members and details of current cluster
          temp_ids <- cluster_data_frame$ID[cluster_data_frame$tag == cluster_properties$tag[i]]
          temp_members <- cluster_data_frame$CDR3b[temp_ids]

          # Extract current motif and starting position range for the motif
          if(parameters$structboundaries) temp_members_frags <- substr(x = temp_members,start = parameters$boundary_size + 1 ,stop = nchar(temp_members) - parameters$boundary_size) else temp_members_frags <- temp_members
          pattern <- strsplit(cluster_properties$tag[i], split = "_")[[1]][[1]]

          # Search all positions of the motif in the cluster members
          temp_members_frags <- rep(temp_members_frags, each = max(nchar(temp_members_frags)-nchar(pattern)+1))
          temp_subs <- substr(temp_members_frags,
                                    start = rep(seq_len(max(nchar(temp_members_frags)-nchar(pattern)+1)), times = length(temp_members)),
                                    stop = rep(nchar(pattern):max(nchar(temp_members_frags)), times = length(temp_members)))
          temp_which <- which(temp_subs == pattern)
          details <- data.frame(id = temp_ids[floor((temp_which-1)/max(nchar(temp_members_frags)-nchar(pattern)+1))+1],
                                      pos = ((temp_which-1) %% max(nchar(temp_members_frags)-nchar(pattern)+1))+1,
                                      stringsAsFactors = FALSE)

          # First assume connections between all cluster members
          if(length(temp_members) >= 2){
            combn_ids <- t(utils::combn(seq_along(temp_members), m = 2))
          } else {
            combn_ids <- t(utils::combn(rep(1, 2), m = 2))
          }
          combn_ids <- combn_ids[details$id[combn_ids[,1]] != details$id[combn_ids[,2]],]
          if(!is.matrix(combn_ids)) combn_ids <- matrix(combn_ids, ncol = 2, byrow = TRUE)

          temp_df <- data.frame(from = details$id[combn_ids[,1]], to = details$id[combn_ids[,2]],
                                      label = rep("local", nrow(combn_ids)), stringsAsFactors = FALSE)

          # Restrict local connections to sequences in which the starting position varies in a certain range of positions
          temp_df <- unique(temp_df[abs(details$pos[combn_ids[,1]]-details$pos[combn_ids[,2]]) < parameters$motif_distance_cutoff,])
          temp_df <- temp_df[cluster_data_frame$tag[temp_df$from] == cluster_data_frame$tag[temp_df$to],]
          temp_df <- t(temp_df)
          return(unlist(temp_df))
        }

        ### If available, set clone network to the local network
        if(length(local_clone_network) == 0) parameters$local_similarities == FALSE else {
          local_clone_network <- data.frame(matrix(unlist(local_clone_network), ncol = 3, byrow = TRUE),
                                                  stringsAsFactors = FALSE)
          colnames(local_clone_network) <- c("from", "to", "label")

          clone_network <- local_clone_network
        }

        message("Restrict local similarities to sequences with motif position varying between ", parameters$motif_distance_cutoff," positions.")
      }

      ### global connections
      if(parameters$global_similarities == TRUE){
        # load BlosumVec to determine interchangeable amino acids defined by BLOSUM62
        BlosumVec <- .get_blosum_vec()

        ### Get local connections between different tuples consisting of CDR3b sequence, v gene and donor
        global_clone_network <- foreach::foreach(i = which(cluster_properties$type == "global")) %dopar% {

          # Get IDs and members of current cluster
          temp_ids <- cluster_data_frame$ID[cluster_data_frame$tag == cluster_properties$tag[i]]
          temp_members <- cluster_data_frame$CDR3b[temp_ids]

          # Assume connections between all cluster members
          if(length(temp_members) >= 2){
            combn_ids <- t(utils::combn(seq_along(temp_members), m = 2))
          } else {
            combn_ids <- t(utils::combn(rep(1, 2), m = 2))
          }
          temp_df <- data.frame(from = temp_ids[combn_ids[,1]], to = temp_ids[combn_ids[,2]],
                                      label = rep("global", nrow(combn_ids)), stringsAsFactors = FALSE)

          # If requested, restrict global connections to sequences with a Hamming distance of 1 and a differing amino acid pair allowed by the BLOSUM62 matrix
          if(parameters$all_aa_interchangeable == FALSE){
            temp_pos <- stringr::str_locate(string = cluster_properties$tag[i], pattern = "%")
            if(parameters$structboundaries == TRUE) temp_pos <- temp_pos+parameters$boundary_size
            temp_df <- unique(temp_df[paste0(substr(cluster_data_frame$CDR3b[temp_df$from], temp_pos, temp_pos),
                                                         substr(cluster_data_frame$CDR3b[temp_df$to], temp_pos, temp_pos)) %in% BlosumVec,])
          }
          temp_df <- t(temp_df)
          return(unlist(temp_df))
        }
        if(parameters$all_aa_interchangeable == FALSE) message("Restrict global connections to sequences with a BLOSUM62 value for the amino acid substitution greater or equal to zero.")

        ### Append global network to the clone network
        if(length(global_clone_network) == 0) parameters$local_similarities == FALSE else {
          global_clone_network <- data.frame(matrix(unlist(global_clone_network), ncol = 3, byrow = TRUE),
                                                   stringsAsFactors = FALSE)
          colnames(global_clone_network) <- c("from", "to", "label")

          if(parameters$local_similarities == FALSE){
            clone_network <- global_clone_network
          } else {
            clone_network <- rbind(clone_network, global_clone_network)
          }
        }
      }
      clone_network$from <- as.numeric(clone_network$from)
      clone_network$to <- as.numeric(clone_network$to)
    }
  }
  clone_network <- unique(clone_network)

  ### To ensure that global connections are preferred over local connections for drawing
  clone_network <- clone_network[order(clone_network$label, decreasing = TRUE),]

  ### Set color of edges
  # Create data frame inlcuding the colors for local and global connections
  leg.col <- data.frame(color=c(local_edge_color,global_edge_color))
  leg.col[] <- lapply(leg.col,as.character)
  row.names(leg.col) <- c("local", "global")

  # Assign the corresponding color to each connection in the clone network
  temp_idx <- match(clone_network$label, row.names(leg.col))
  matched <- !is.na(temp_idx)
  clone_network$color[matched] <- leg.col[temp_idx[matched], 1]

  #################################################################
  ##                     Prepare graph nodes                     ##
  #################################################################

  ### Create the data frame framework for all node information
  vert.info <- data.frame(id = unique(c(clone_network$from,clone_network$to)),
                                stringsAsFactors = FALSE)
  vert.info$size <- rep(20, nrow(vert.info))
  vert.info$color <- rep("gray", nrow(vert.info))
  vert.info$label <- cluster_data_frame$CDR3b[vert.info$id]
  vert.info$title <- paste0("<p><b>",vert.info$label, "</b>")

  ### If requested, determine the node sizes
  if(!is.null(size_info)){
    if(!(size_info %in% colnames(cluster_data_frame))) stop("Column named ", size_info, " determining node size is not found.")
    if(suppressWarnings(any(is.na(as.numeric(cluster_data_frame[,size_info])))) == TRUE) stop("Column named ", size_info, " determining node size contains non-numeric values.")
  }

  if(is.null(size_info)){
    vert.info$size <- rep(20, nrow(vert.info))
  } else{
    vert.info$size <- cluster_data_frame[vert.info$id, size_info]

    if(absolute_size == FALSE){
      vert.info$size <- (vert.info$size-min(vert.info$size))/(max(vert.info$size)-min(vert.info$size))*8+12
    }
  }

  ### Determine the node colors
  if(!(color_info %in% c("none", "total.score", colnames(cluster_data_frame)))){
    stop("Column named ", color_info, " determining node color is not found.")
  }
  if(color_info == "color" && !(any(plotfunctions::isColor(cluster_data_frame[, color_info]) == FALSE))) stop("Column ", color_info, " determining node color has to contain only values that represent colors.")

  color.scale = ""
  if("color" %in% colnames(cluster_data_frame)){
    # Use the user specified colors
    vert.info$color <- cluster_data_frame$color[vert.info$id]
  } else if(color_info == "total.score"){
    # Color accoring to the total cluster score with a logarithmic scale
    min_score <- floor(log10(min(cluster_data_frame$total.score)))
    max_score <- ceiling(log10(max(cluster_data_frame$total.score)))

    node.col <- data.frame(logvalue = (seq(from = min_score, to = max_score, length.out = 50)), stringsAsFactors = FALSE)
    node.col$value <- 10^node.col$logvalue
    node.col$color <- color_palette(nrow(node.col))

    vert.info$color <- vapply(X=cluster_data_frame$total.score[vert.info$id], FUN = function(x){
      return(node.col$color[which.min(abs(node.col$value - x))])
    }, FUN.VALUE = c("#00ff00"))

    color.scale <- "log"
  } else if(color_info == "none"){} else {
    ### Automatic coloring based on the desired property

    # Coloring for numeric values by determining and using a linear or logarthmic scale
    if(is.numeric(cluster_data_frame[, color_info])){
      vals <- sort(unique(cluster_data_frame[, color_info]))
      lin_model <- stats::lm(vals ~ seq_along(vals))
      log_model <- stats::lm(log(vals) ~ seq_along(vals))

      if(summary(lin_model)$sigma < summary(log_model)$sigma) color.scale <- "linear" else color.scale <- "log"

      if(color.scale == "log"){
        min_score <- floor(log10(min(cluster_data_frame$total.score)))
        max_score <- ceiling(log10(max(cluster_data_frame$total.score)))
      } else {
        min_score <- floor(min(cluster_data_frame$total.score))
        max_score <- ceiling(max(cluster_data_frame$total.score))
      }


      node.col <- data.frame(logvalue = (seq(from = min_score, to = max_score, length.out = 50)), stringsAsFactors = FALSE)
      if(color.scale == "log") node.col$value <- 10^node.col$logvalue else node.col$value <- node.col$logvalue
      node.col$color <- color_palette(nrow(node.col))

      vert.info$color <- vapply(X=cluster_data_frame$total.score[vert.info$id], FUN = function(x){
        return(node.col$color[which.min(abs(node.col$value - x))])
      }, FUN.VALUE = c("#00ff00"))
    } else {
      # Coloring for categorial values
      vals <- sort(as.character(unique(cluster_data_frame[, color_info])))
      cols <- color_palette(length(vals))
      temp_idx <- match(as.character(cluster_data_frame[vert.info$id, color_info]), vals)
      matched <- !is.na(temp_idx)
      vert.info$color[matched] <- cols[temp_idx[matched]]

      if(!is.null(vals) && !is.null(cols)){
        # Create legend only if number of values is managable (<=6)
        if(length(vals) <= 6) vert.info$collab <- cluster_data_frame[vert.info$id, color_info]
      }
    }
  }

  ### Write the description for every node with the desired information (if possible)
  this_column_names <- c()
  if(!("gliph_version" %in% names(parameters))){
    this_column_names <- c("tag", "total.score", show_additional_columns)
  } else {
    if(parameters$gliph_version == 1) this_column_names <- c("tag", "total.score", show_additional_columns) else this_column_names <- c("tag", "fisher.score", "total.score", show_additional_columns)
  }

  all_column_names <- c(colnames(cluster_list[[1]]), colnames(cluster_properties))
  this_column_names <- this_column_names[this_column_names %in% all_column_names]
  for(actCol in this_column_names){
    vert.info$title <- paste0(vert.info$title, "<br>", actCol, ": ",cluster_data_frame[vert.info$id,actCol])
  }
  vert.info$title <- paste0(vert.info$title, "</p>")
  vert.info$group <- cluster_data_frame$tag[vert.info$id]


  ##################################################################
  ##       Convert edge and node information for visNetwork       ##
  ##################################################################
  vertex.info <- vert.info
  edge.info <- clone_network
  layout <- "layout_components"

  ### Node properties
  vertexes <- vertex.info$id
  num.v <- length(vertexes)
  v.size <- 8
  v.label <- vertexes
  v.group <- NA
  v.shape <- rep("dot")
  v.title <- vertexes
  v.color <- rep("gray", num.v)
  v.shadow <- rep(FALSE, num.v)
  if ("size" %in% names(vertex.info)) v.size = as.numeric(vertex.info$size)
  if ("label" %in% names(vertex.info)) v.label = as.character(vertex.info$label)
  if ("group" %in% names(vertex.info)) v.group = as.character(vertex.info$group)
  if ("shape" %in% names(vertex.info)) v.shape = as.character(vertex.info$shape)
  if ("title" %in% names(vertex.info)) v.title = as.character(vertex.info$title)
  if ("color" %in% names(vertex.info)) v.color = as.character(vertex.info$color)
  if ("shadow" %in% names(vertex.info)) v.shadow = as.logical(vertex.info$shadow)
  nodes <- data.frame(id = vertexes,
                            color = list(background = v.color, border = "black", highlight = "red"),
                            size=v.size,
                            label=v.label,
                            title=v.title,
                            group=v.group,
                            shape=v.shape,
                            shadow=v.shadow)

  ### Edge properties
  eds <- edge.info
  num.e <- dim(eds)[1]
  e.label <- NA
  e.length <- rep(150,num.e)
  e.width <- rep(10,num.e)
  e.color <- rep("gray", num.e)
  e.arrows <- rep("",num.e)
  e.dashes <- rep(FALSE,num.e)
  e.title <- rep("",num.e)
  e.smooth <- rep(FALSE,num.e)
  e.shadow <- rep(FALSE,num.e)
  if ("length" %in% names(edge.info)) e.length = as.numeric(edge.info$length)
  if ("label" %in% names(edge.info)) e.label = as.character(edge.info$label)
  if ("width" %in% names(edge.info)) e.width = as.numeric(edge.info$width)
  if ("color" %in% names(edge.info)) e.color = as.character(edge.info$color)
  if ("arrows" %in% names(edge.info)) e.arrows = as.character(edge.info$arrows)
  if ("dashes" %in% names(edge.info)) e.dashes = as.logical(edge.info$dashes)
  if ("title" %in% names(edge.info)) e.title = as.character(edge.info$title)
  if ("smooth" %in% names(edge.info)) e.smooth = as.logical(edge.info$smooth)
  if ("shadow" %in% names(edge.info)) e.shadow = as.logical(edge.info$shadow)
  edges <- data.frame(from = eds$from, to = eds$to,
                            length = e.length,
                            width = e.width,
                            color = e.color,
                            arrows = e.arrows,
                            dashes = e.dashes,
                            title = e.title,
                            smooth = e.smooth,
                            shadow = e.shadow )

  ### Legend properties
  ledges <- NULL
  if ("label" %in% names(edge.info) && "color" %in% names(edge.info)) {
    leg.info <- unique(cbind(as.character(edge.info$label),as.character(edge.info$color)))
    ledges <- data.frame(color = leg.info[,2], label = leg.info[,1],
                               arrows=rep("",length(leg.info[,1])), width=rep(4,length(leg.info[,1])))
  }
  lenodes <- NULL
  if ("color" %in% names(vertex.info) && "collab" %in% names(vertex.info)) {
    leg.info <- unique(cbind(as.character(vertex.info$collab),as.character(vertex.info$color)))
    lenodes <- data.frame(label = leg.info[,1], color = leg.info[,2], shape="dot", size=10)
  }
  if(color.scale != ""){
    leg.info <- node.col[seq(from = 1, to = nrow(node.col), by = (nrow(node.col)-1)/5), c(1, 3)]
    if(color.scale == "log") leg.info[,1] <- paste0("10^(", round(leg.info[,1], digits = 3), ")") else leg.info[,1] <- as.character(round(leg.info[,1], digits = 5))
    lenodes <- data.frame(label = leg.info[,1], color = leg.info[,2], shape="dot", size=10)
  }

  .stop_parallel()

  message("Drawing the graph.")

  ret <- visNetwork::visLegend(visNetwork::visOptions(visNetwork::visIgraphLayout(visNetwork::visNetwork(nodes = nodes, edges = edges), layout = "layout_components"),
                                                                   highlightNearest = list(enabled = TRUE, degree = 1,algorithm = "hierarchical"),
                                                                   selectedBy = list(variable = "group", multiple = TRUE), manipulation = TRUE),
                                            addEdges = ledges, addNodes = lenodes, useGroups = FALSE, position = "right",
                                            width=0.15, zoom = FALSE)

  ### Closing time!
  return(ret)
}
