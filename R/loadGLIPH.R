#' Load saved GLIPH results from disk
#'
#' Reads the tab-delimited output files produced by \code{\link{runGLIPH}} (when
#' \code{result_folder} was specified) and reconstructs the same list structure
#' that \code{runGLIPH()} returns.
#'
#' @param result_folder Path to the folder containing the saved GLIPH output
#'   files.
#'
#' @return A \code{list} with the same structure as the return value of
#'   \code{\link{runGLIPH}}, including elements such as \code{cluster_list},
#'   \code{cluster_properties}, \code{motif_enrichment}, \code{connections}, and
#'   \code{parameters}.
#'
#' @export
loadGLIPH <- function(result_folder = ""){

  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################

  ### result_folder
  if(!is.character(result_folder)) stop("result_folder has to be a character object")
  if(length(result_folder) > 1) stop("result_folder has to be a single path")
  if(result_folder == "") stop("To load output files, the path must be specified unter the parameter 'result_folder'.")
  if(substr(result_folder,nchar(result_folder),nchar(result_folder)) != "/") result_folder <- paste0(result_folder,"/")
  message("Files are loaded from the following folder: ", result_folder)
  if(!dir.exists(result_folder)) stop("Specified path under the parameter 'result folder' does not exist.")

  ##################################################################
  ##                          Load files                          ##
  ##################################################################

  ### load parameter file: parameters
  fname <- paste0(result_folder, "parameter.txt")
  if(!file.exists(fname)) stop("File named 'parameters.txt' is missing in the specified folder.")
  para_df <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
  parameters <- lapply(seq_len(nrow(para_df)), function(x){
    if(suppressWarnings(any(is.na(as.numeric(para_df[x,2])))) == FALSE) return(as.numeric(para_df[x,2])) else return(para_df[x,2])
    })
  names(parameters) <- para_df[,1]
  parameters$motif_length <- as.numeric(unlist(strsplit(parameters$motif_length, split = ",")))
  parameters$lcminove <- as.numeric(unlist(strsplit(parameters$lcminove, split = ",")))
  message("Loaded 'parameters' from parameter.txt")

  if(!("gliph_version" %in% names(parameters))){
    message("Output of gliph_combined function is loaded.")

    ### load sample log
    fname <- paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_log.txt")
    if(file.exists(fname)){
      sample_log <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
      message("Loaded 'sample_log' from ", fname)
    } else sample_log <- NULL

    ### load selected_motifs
    fname <- paste0(result_folder, "local_similarities_minp_",parameters$lcminp, "_minove_", paste(parameters$lcminove, collapse = "_"), "_kmer_mindepth_", parameters$lckmer_mindepth, ".txt")
    if(file.exists(fname)){
      selected_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      message("Loaded 'selected_motifs' from ", fname)
    } else selected_motifs <- NULL

    ### load all_motifs
    fname <- paste0(result_folder,"all_motifs.txt")
    if(file.exists(fname)){
      all_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      message("Loaded 'all_motifs' from ", fname)
    } else all_motifs <- NULL

    ### build motif_enrichment
    motif_enrichment <- list(selected_motifs = selected_motifs, all_motifs = all_motifs)

    ### load selected_structs
    fname <- paste0(result_folder, "global_similarities_minp_",parameters$gcminp, "_minove_", paste(parameters$gcminove, collapse = "_"), "_kmer_mindepth_", parameters$gckmer_mindepth, ".txt")
    if(file.exists(fname)){
      selected_structs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      message("Loaded 'selected_structs' from ", fname)
    } else selected_structs <- NULL

    ### load all_structs
    fname <- paste0(result_folder, "all_global_similarities.txt")
    if(file.exists(fname)){
      all_structs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      message("Loaded 'all_structs' from ", fname)
    } else all_structs <- NULL

    ### build global_enrichment
    global_enrichment <- list(selected_structs = selected_structs, all_structs = all_structs)

    ### load connections
    fname <- paste0(result_folder,"clone_network.txt")
    if(file.exists(fname)){
      connections <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
      message("Loaded 'connections' from ", fname)
    } else connections <- NULL

    ### load cluster_properties
    fname <- paste0(result_folder,"convergence_groups.txt")
    if(file.exists(fname)){
      cluster_properties <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      message("Loaded 'cluster_properties' from ", fname)
    } else cluster_properties <- NULL

    ### load cluster_list
    fname <- paste0(result_folder, "cluster_member_details.txt")
    if(file.exists(fname)){
      cluster_list <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      for(i in seq_len(ncol(cluster_list))){
        if(suppressWarnings(any(is.na(as.numeric(cluster_list[,i])))) == FALSE) cluster_list[,i] <- as.numeric(cluster_list[,i])
      }
      tag_names <- unique(cluster_list$tag)
      cluster_list <- lapply(tag_names, function(x){return(cluster_list[cluster_list$tag == x,-1])})
      names(cluster_list) <- tag_names

      message("Loaded 'cluster_list' from ", fname)
    } else cluster_list <- NULL

    output <- list(motif_enrichment = motif_enrichment,
                         global_enrichment = global_enrichment,
                         connections = connections,
                         cluster_properties = cluster_properties,
                         cluster_list = cluster_list,
                         parameters = parameters)

  } else {
    if(parameters$gliph_version == 1){
      message("Output of turboGliph function is loaded.")

      ### load sample log
      fname <- paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_log.txt")
      if(file.exists(fname)){
        sample_log <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
        message("Loaded 'sample_log' from ", fname)
      } else sample_log <- NULL


      ### load selected_motifs
      fname <- paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_minp",parameters$lcminp,"_ove", paste(parameters$lcminove, collapse = "_"),".txt")
      if(file.exists(fname)){
        selected_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        message("Loaded 'selected_motifs' from ", fname)
      } else selected_motifs <- NULL

      ### load all_motifs
      fname <- paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_all_motifs.txt")
      if(file.exists(fname)){
        all_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        message("Loaded 'all_motifs' from ", fname)
      } else all_motifs <- NULL

      ### build motif_enrichment
      motif_enrichment <- list(selected_motifs = selected_motifs, all_motifs = all_motifs)

      ### load connections
      fname <- paste0(result_folder,"clone_network.txt")
      if(file.exists(fname)){
        connections <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
        message("Loaded 'connections' from ", fname)
      } else connections <- NULL

      ### load cluster_properties
      fname <- paste0(result_folder,"convergence_groups.txt")
      if(file.exists(fname)){
        cluster_properties <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        message("Loaded 'cluster_properties' from ", fname)
      } else cluster_properties <- NULL

      ### load cluster_list
      fname <- paste0(result_folder, "cluster_member_details.txt")
      if(file.exists(fname)){
        cluster_list <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        for(i in seq_len(ncol(cluster_list))){
          if(suppressWarnings(any(is.na(as.numeric(cluster_list[,i])))) == FALSE) cluster_list[,i] <- as.numeric(cluster_list[,i])
        }
        tag_names <- unique(cluster_list$tag)
        cluster_list <- lapply(tag_names, function(x){return(cluster_list[cluster_list$tag == x,-1])})
        names(cluster_list) <- tag_names

        message("Loaded 'cluster_list' from ", fname)
      } else cluster_list <- NULL

      output <- list(sample_log = sample_log,
                           motif_enrichment = motif_enrichment,
                           connections = connections,
                           cluster_properties = cluster_properties,
                           cluster_list = cluster_list,
                           parameters = parameters)
    }
    if(parameters$gliph_version == 2){
      message("Output of gliph2 function is loaded.")

      ### load selected_motifs
      fname <- paste0(result_folder, "local_similarities_minp_",parameters$lcminp, "_minove_", paste(parameters$lcminove, collapse = "_"), "_kmer_mindepth_", parameters$kmer_mindepth, ".txt")
      if(file.exists(fname)){
        selected_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        message("Loaded 'selected_motifs' from ", fname)
      } else selected_motifs <- NULL

      ### load all_motifs
      fname <- paste0(result_folder, "all_motifs.txt")
      if(file.exists(fname)){
        all_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        message("Loaded 'all_motifs' from ", fname)
      } else all_motifs <- NULL

      ### build motif_enrichment
      motif_enrichment <- list(selected_motifs = selected_motifs, all_motifs = all_motifs)

      ### load global_enrichment
      fname <- paste0(result_folder, "global_similarities.txt")
      if(file.exists(fname)){
        global_enrichment <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        message("Loaded 'global_enrichment' from ", fname)
      } else global_enrichment <- NULL

      ### load connections
      fname <- paste0(result_folder,"clone_network.txt")
      if(file.exists(fname)){
        connections <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
        message("Loaded 'connections' from ", fname)
      } else connections <- NULL

      ### load cluster_properties
      fname <- paste0(result_folder,"convergence_groups.txt")
      if(file.exists(fname)){
        cluster_properties <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        message("Loaded 'cluster_properties' from ", fname)
      } else cluster_properties <- NULL

      ### load cluster_list
      fname <- paste0(result_folder, "cluster_member_details.txt")
      if(file.exists(fname)){
        cluster_list <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        for(i in seq_len(ncol(cluster_list))){
          if(suppressWarnings(any(is.na(as.numeric(cluster_list[,i])))) == FALSE) cluster_list[,i] <- as.numeric(cluster_list[,i])
        }
        tag_names <- unique(cluster_list$tag)
        cluster_list <- lapply(tag_names, function(x){return(cluster_list[cluster_list$tag == x,-1])})
        names(cluster_list) <- tag_names

        message("Loaded 'cluster_list' from ", fname)
      } else cluster_list <- NULL

      output <- list(motif_enrichment = motif_enrichment,
                           global_enrichment = global_enrichment,
                           connections = connections,
                           cluster_properties = cluster_properties,
                           cluster_list = cluster_list,
                           parameters = parameters)
    }
  }

  return(output)
}

