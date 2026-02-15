#' Score CDR3 clusters using the GLIPH or GLIPH2 algorithm
#'
#' Calculates scores for CDR3 clusters following the GLIPH and GLIPH2
#' scoring procedures. Depending on the information provided, a final
#' score is computed from up to five cluster properties: cluster size,
#' enrichment of CDR3 lengths, enrichment of V genes, enrichment of
#' clonal expansions, and enrichment of common HLA alleles.
#'
#' @param cluster_list A \code{list} where each element contains a
#' \code{data.frame} of CDR3b sequences and additional information
#' needed for scoring. Corresponds to the \code{$cluster_list} element
#' returned by \code{\link{runGLIPH}}.
#' @param cdr3_sequences A \code{vector} or \code{data.frame} of CDR3
#' sequences and optional metadata. The columns must be named as
#' specified below in arbitrary order:
#' \itemize{
#'   \item{\code{"CDR3b"}:}{ CDR3 sequences of beta chains.}
#'   \item{\code{"TRBV"}:}{ Optional. V-genes of beta chains.}
#'   \item{\code{"patient"}:}{ Optional. Donor index for the
#'     corresponding sequence, composed of a donor identifier and an
#'     optional experimental condition separated by a colon
#'     (e.g., \code{09/0410:MtbLys}). Only the identifier before the
#'     colon is used for HLA scoring.}
#'   \item{\code{"HLA"}:}{ Optional. Comma-separated HLA alleles for
#'     the corresponding donor in standard notation
#'     (e.g., \code{DPA1*01:03}). Information after the colon in each
#'     allele is ignored during HLA scoring.}
#'   \item{\code{"counts"}:}{ Optional. Clone frequency.}
#' }
#' @param refdb_beta A \code{character} string or \code{data.frame}
#' specifying the reference database. When a \code{data.frame} is
#' supplied, CDR3b sequences must be in the first column and V-gene
#' information (if available) in the second column. Built-in databases
#' include \code{"human_v1.0_CD4"}, \code{"human_v1.0_CD8"},
#' \code{"human_v1.0_CD48"}, \code{"human_v2.0_CD4"},
#' \code{"human_v2.0_CD8"}, \code{"human_v2.0_CD48"},
#' \code{"mouse_v1.0_CD4"}, \code{"mouse_v1.0_CD8"},
#' \code{"mouse_v1.0_CD48"}, and the legacy alias
#' \code{"gliph_reference"} (= \code{"human_v1.0_CD48"}).
#' See \code{\link{reference_list}} for details.
#' \strong{Default:} \code{"human_v2.0_CD48"}
#' @param v_usage_freq A \code{data.frame} with V-gene alleles in the
#' first column and their naive-repertoire frequencies in the second
#' column. \strong{Default:} \code{NULL}
#' @param cdr3_length_freq A \code{data.frame} with CDR3 lengths in the
#' first column and their naive-repertoire frequencies in the second
#' column. \strong{Default:} \code{NULL}
#' @param ref_cluster_size A \code{character} string defining which
#' cluster-size probabilities to use for scoring.
#' \itemize{
#'   \item{\code{"original"}:}{ Standard probabilities from the
#'     original algorithm, constant across sample sizes.}
#'   \item{\code{"simulated"}:}{ Probabilities estimated for different
#'     sample sizes via a 500-step simulation using random sequences
#'     from the reference database.}
#' }
#' \strong{Default:} \code{"original"}
#' @param gliph_version A \code{numeric} value indicating the algorithm
#' version.
#' \itemize{
#'   \item{\code{1}:}{ GLIPH scoring (product of individual scores
#'     multiplied by 0.064).}
#'   \item{\code{2}:}{ GLIPH2 scoring (product of individual scores
#'     only).}
#' }
#' \strong{Default:} \code{1}
#' @param sim_depth A \code{numeric} value for simulated resampling
#' depth in non-parametric convergence significance tests. Higher
#' values increase runtime but improve reproducibility.
#' \strong{Default:} \code{1000}
#' @param hla_cutoff A \code{numeric} threshold below which HLA
#' probability scores are considered significant.
#' \strong{Default:} \code{0.1}
#' @param n_cores A \code{numeric} value for the number of cores to
#' use. When \code{NULL}, it is set to the number of available cores
#' minus two. \strong{Default:} \code{1}
#'
#' @return A \code{data.frame} of cluster scoring results. The first
#' column contains the total score and additional columns contain up to
#' five individual scores (cluster size, CDR3 length enrichment, V-gene
#' enrichment, clonal expansion enrichment, and common HLA enrichment).
#'
#' @examples
#' \dontrun{
#' utils::data("gliph_input_data")
#'
#' res <- runGLIPH(cdr3_sequences = gliph_input_data[seq_len(200), ],
#'                 sim_depth = 100,
#'                 n_cores = 1)
#'
#' scoring_results <- clusterScoring(
#'     cluster_list = res$cluster_list,
#'     cdr3_sequences = gliph_input_data[seq_len(200), ],
#'     refdb_beta = "human_v2.0_CD48",
#'     gliph_version = 1,
#'     sim_depth = 100,
#'     n_cores = 1)
#' }
#'
#' @references Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
#' @references https://github.com/immunoengineer/gliph
#' @import foreach
#' @export
clusterScoring <- function(cluster_list,
                            cdr3_sequences,
                            refdb_beta = "human_v2.0_CD48",
                            v_usage_freq = NULL,
                            cdr3_length_freq = NULL,
                            ref_cluster_size = "original",
                            gliph_version = 1,
                            sim_depth = 1000,
                            hla_cutoff = 0.1,
                            n_cores = 1){
  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################

  ### cluster_list
  if(!is.list(cluster_list)) stop("parameter 'cluster_list' has to be an object of class 'list'.")

  ### cdr3_sequences
  if(is.atomic(cdr3_sequences)) cdr3_sequences <- data.frame("CDR3b" = cdr3_sequences)
  if(!is.data.frame(cdr3_sequences)) stop("parameter 'cdr3_sequences' has to be an object of class 'data.frame'.")
  cdr3_sequences[] <- lapply(cdr3_sequences, as.character)

  ### refdb_beta
  if(!is.data.frame(refdb_beta)){
    valid_names <- .valid_reference_names()
    if(length(refdb_beta) != 1 || !is.character(refdb_beta)){
      stop("refdb_beta must be a data frame or one of: ",
           paste(sQuote(valid_names), collapse = ", "))
    } else if(!(refdb_beta %in% valid_names)){
      stop("refdb_beta must be a data frame or one of: ",
           paste(sQuote(valid_names), collapse = ", "))
    }
  }

  ### v_usage_freq
  if(!is.null(v_usage_freq)){
    if(is.data.frame(v_usage_freq)){
      if(ncol(v_usage_freq) < 2) stop("v_usage_freq has to be a data frame containing V-gene information in the first column and the corresponding frequency in a naive  T-cell repertoire in the second column.")
      if(nrow(v_usage_freq) < 1) stop("v_usage_freq has to contain at least one row.")
      if(suppressWarnings(any(is.na(as.numeric(v_usage_freq[,2])))) == TRUE){
        stop("The second column of v_usage_freq must contain the frequency of the corresponding V-gene in the first column in a naive T-cell repertoire.")
      } else v_usage_freq[,2] <- as.numeric(v_usage_freq[,2])

    } else {stop("v_usage_freq has to be a data frame containing V-gene information in the first column and the corresponding frequency in a naive T-cell repertoire in the second column.")}
  }

  ### cdr3_length_freq
  if(!is.null(cdr3_length_freq)){
    if(is.data.frame(cdr3_length_freq)){
      if(ncol(cdr3_length_freq) < 2) stop("cdr3_length_freq has to be a data frame containing CDR3 lengths in the first column and the corresponding frequency in a naive  T-cell repertoire in the second column.")
      if(nrow(cdr3_length_freq) < 1) stop("cdr3_length_freq has to contain at least one row.")
      if(suppressWarnings(any(is.na(as.numeric(cdr3_length_freq[,2])))) == TRUE){
        stop("The second column of cdr3_length_freq must contain the frequency of the corresponding CDR3 length in the first column in a naive T-cell repertoire.")
      } else cdr3_length_freq[,2] <- as.numeric(cdr3_length_freq[,2])

    } else {stop("cdr3_length_freq has to be a data frame containing CDR3 lengths in the first column and the corresponding frequency in a naive T-cell repertoire in the second column.")}
  }

  ### ref_cluster_size
  if(!(ref_cluster_size %in% c("original", "simulated") || !is.character(ref_cluster_size) || length(ref_cluster_size) > 1)){
    stop("ref_cluster_size has to be either 'original' or 'simulated'.")
  }

  ### gliph_version
  if(!(gliph_version %in% c(1,2))) stop("gliph_version has to be either 1 or 2.")

  ### sim_depth
  if(!is.numeric(sim_depth)) stop("sim_depth has to be numeric")
  if(length(sim_depth) > 1) stop("sim_depth has to be a single number")
  if(sim_depth < 1) stop("sim_depth must be at least 1")
  sim_depth <- round(sim_depth)

  ### hla_cutoff
  if(!is.numeric(hla_cutoff)) stop("hla_cutoff has to be numeric")
  if(length(hla_cutoff) > 1) stop("hla_cutoff has to be a single number")
  if(hla_cutoff > 1 || hla_cutoff < 0) stop("hla_cutoff must be between 0 and 1")

  ### n_cores
  if(is.null(n_cores)) n_cores <- max(1, parallel::detectCores()-2) else {
    if(!is.numeric(n_cores)) stop("n_cores has to be numeric")
    if(length(n_cores) > 1) stop("n_cores has to be a single number")
    if(n_cores < 1) stop("n_cores must be at least 1")

    n_cores <- min(n_cores, parallel::detectCores()-2)
  }

  #################################################################
  ##                         Preparation                         ##
  #################################################################

  ### Which scores can be calculated from the dataset?
  score_names <- c("network.size.score","cdr3.length.score")
  if("TRBV" %in% colnames(cdr3_sequences)){
    vgene_info <- TRUE
    score_names <- c(score_names, "vgene.score")
  } else vgene_info <- FALSE
  if("counts" %in% colnames(cdr3_sequences)) {
    counts_info <- TRUE
    cdr3_sequences$counts <- as.numeric(cdr3_sequences$counts)
    cdr3_sequences$counts[is.na(cdr3_sequences$counts)] <- 1
    score_names <- c(score_names, "clonal.expansion.score")
  } else counts_info <- FALSE
  if("patient" %in% colnames(cdr3_sequences)) patient_info <- TRUE else patient_info <- FALSE
  if("HLA" %in% colnames(cdr3_sequences)) hla_info <- TRUE else hla_info <- FALSE
  if(hla_info == TRUE && patient_info == TRUE){
    if("patient" %in% colnames(cdr3_sequences) && "HLA" %in% colnames(cdr3_sequences)){
      cdr3_sequences <- cdr3_sequences[cdr3_sequences$HLA != "" & !is.na(cdr3_sequences$HLA),]
      if(nrow(cdr3_sequences > 0)) score_names <- c(score_names, "hla.score", "lowest.hlas") else hla_info <- FALSE
    }
  }

  ### load or generate reference tables from reference database
  # ref_cluster_sizes:            data frame containing the cluster size in the first and the probability of observing a cluster with this size
  #                               in a sample from the reference database in the second column
  # vgene_ref_frequencies:        vector containing the (relative) frequencies of v gene usage
  # cdr3_length_ref_frequencies:  vector containing the (relative) frequencies of CDR3b lengths
  utils::data(ref_cluster_sizes, envir = environment(), package = "immGLIPH")
  ref_cluster_sizes <- ref_cluster_sizes[[ref_cluster_size]]

  if(is.character(refdb_beta) && refdb_beta %in% .valid_reference_names()){
    reference_list <- .get_reference_list()
    refdb_name <- refdb_beta
    vgene_ref_frequencies <- reference_list[[refdb_name]]$vgene_frequencies$freq
    cdr3_length_ref_frequencies <- reference_list[[refdb_name]]$cdr3_length_frequencies$freq
  } else {
    # User-provided data frame: compute frequencies from the data
    if("TRBV" %in% colnames(refdb_beta)){
      vgene_ref_frequencies <- as.data.frame(table(refdb_beta$TRBV))
      vgene_ref_frequencies <- vgene_ref_frequencies$Freq/sum(vgene_ref_frequencies$Freq)
    } else {
      reference_list <- .get_reference_list()
      vgene_ref_frequencies <- reference_list[["gliph_reference"]]$vgene_frequencies$freq
    }
    cdr3_length_ref_frequencies <- as.data.frame(table(nchar(refdb_beta$CDR3b)))
    cdr3_length_ref_frequencies <- cdr3_length_ref_frequencies$Freq/sum(cdr3_length_ref_frequencies$Freq)
  }

  if(!is.null(v_usage_freq)) vgene_ref_frequencies <- as.numeric(v_usage_freq[,2])
  if(!is.null(cdr3_length_freq)) cdr3_length_ref_frequencies <- as.numeric(cdr3_length_freq[,2])

  ### Obtain the distribution of all patients and HLA alleles in the sample
  # all_patients:     vector containing all unique patient indices
  # all_hlas:         vector containing all unique HLA alleles
  # all_patient_hlas: list whose elements are named after the patients and contain the patient's HLA alleles in a vector.
  all_patient_hlas <- c()
  if(hla_info == TRUE && patient_info == TRUE){
    cdr3_sequences$patient <- gsub(":.*", "",cdr3_sequences$patient)
    all_patients <- sort(unique(cdr3_sequences$patient))
    all_patients <- all_patients[!is.na(all_patients)]
    all_hlas <- unlist(strsplit(unique(cdr3_sequences$HLA), ","))
    all_hlas <- all_hlas[!is.na(all_hlas)]
    all_hlas <- sort(unique(gsub(":.*", "", all_hlas, perl = TRUE)))
    num_patients <- length(all_patients)
    num_HLAs <- length(all_hlas)

    all_patient_hlas <- lapply(all_patients, function(x){
      sort(unique(gsub(":.*", "", unlist(strsplit(cdr3_sequences$HLA[cdr3_sequences$patient == x][1], ",")), perl = TRUE)))
    })
    names(all_patient_hlas) <- all_patients

    all_hlas <- data.frame(HLA = all_hlas)
    all_hlas$counts <- apply(all_hlas, 1, function(x){
      val <- 0
      for(pat in all_patients){
        if(x %in% all_patient_hlas[[pat]]) val <- val+1
      }
      val
    })
  }

  #################################################################
  ##                           Scoring                           ##
  #################################################################

  doParallel::registerDoParallel(n_cores)

  actCluster <- NULL
  res <- foreach::foreach(actCluster = seq_along(cluster_list)) %dopar% {

    ### Get sequences and information of current cluster
    act_seq_infos <- cluster_list[[actCluster]]
    num_members <- nrow(act_seq_infos) # number of ALL members
    ori_num_members <- length(unique(act_seq_infos$CDR3b)) # number of all unique CDR3b sequences
    all_scores <- c()

    ### Get network size score from lookup table
    score_network_size <- 1
    nearest_sample_size <- order(abs(1-as.numeric(colnames(ref_cluster_sizes)[-1])/nrow(cdr3_sequences)))[1]
    if(ori_num_members > 100){
      score_network_size <-  ref_cluster_sizes[100,nearest_sample_size+1]
    } else {
      score_network_size <-  ref_cluster_sizes[ori_num_members,nearest_sample_size+1]
    }
    all_scores <- c(all_scores, score_network_size)

    ### Enrichment of CDR3 length (spectratype) within cluster
    score_cdr3_length <- c()
    # calculate score of sample (product of all frequencies)
    pick_freqs <- data.frame(table(nchar(unique(act_seq_infos$CDR3b))))
    colnames(pick_freqs) <- c("object", "probs")
    pick_freqs$probs <- pick_freqs$probs/ori_num_members
    sample_score <- round(prod(pick_freqs$probs), digits = 3)

    # generate random subsamples
    random_subsample <- list()
    for(i in seq_len(sim_depth)){
      random_subsample[[i]] <- sample.int(n = length(cdr3_length_ref_frequencies), size = ori_num_members,
                                                prob = cdr3_length_ref_frequencies, replace = TRUE)
    }

    # calculate score of subsamples (product of all frequencies)
    pick_freqs <- stringdist::seq_qgrams(.list = random_subsample)[,-1]/ori_num_members
    pick_freqs[pick_freqs == 0] <- 1
    pick_freqs <- round(exp(colSums(log(pick_freqs))), digits = 3) # vectorized way to calculate the product of each column
    if(gliph_version == 1){
      score_cdr3_length <- sum(pick_freqs >= sample_score)/sim_depth
    } else {
      score_cdr3_length <- sum(pick_freqs > sample_score)/sim_depth
    }
    if(score_cdr3_length == 0) score_cdr3_length <- 1/sim_depth # minimum score of 1/sim_depth

    all_scores <- c(all_scores, score_cdr3_length)

    ### Enrichment of v genes within cluster
    score_vgene <- c()
    if(vgene_info == TRUE){

      # calculate score of sample (product of all frequencies)
      pick_freqs <- data.frame(table(act_seq_infos$TRBV))
      colnames(pick_freqs) <- c("object", "probs")
      pick_freqs$probs <- pick_freqs$probs/num_members
      sample_score <- round(prod(pick_freqs$probs), digits = 3)

      # generate random subsamples
      random_subsample <- list()
      for(i in seq_len(sim_depth)){
        random_subsample[[i]] <- sample.int(n = length(vgene_ref_frequencies), size = num_members,
                                                  prob = vgene_ref_frequencies, replace = TRUE)
      }

      # calculate score of subsamples (product of all frequencies)
      pick_freqs <- stringdist::seq_qgrams(.list = random_subsample)[,-1]/num_members
      pick_freqs[pick_freqs == 0] <- 1
      pick_freqs <- round(exp(colSums(log(pick_freqs))), digits = 3) # vectorized way to calculate the product of each column
      if(gliph_version == 1){
        score_vgene <- sum(pick_freqs >= sample_score)/sim_depth
      } else {
        score_vgene <- sum(pick_freqs > sample_score)/sim_depth
      }

      if(score_vgene == 0) score_vgene <- 1/sim_depth # minimum score of 1/sim_depth
      all_scores <- c(all_scores, score_vgene)
    }

    ### Enrichment of clonal expansion within cluster
    score_clonal_expansion <- c()
    if(counts_info == TRUE){
      sample_score <- sum(as.numeric(act_seq_infos$counts))/num_members
      counter <- 0
      for(i in seq_len(sim_depth)){
        random_subsample <- sample(x = cdr3_sequences$counts, size = num_members, replace = FALSE)
        test_score <- sum(as.numeric(random_subsample))/num_members
        if(test_score>=sample_score) counter <- counter+1
      }
      if(counter == 0) score_clonal_expansion <- 1/sim_depth else score_clonal_expansion <- counter/sim_depth
      score_clonal_expansion <- round(score_clonal_expansion, digits = 3)

      all_scores <- c(all_scores, score_clonal_expansion)
    }

    ### Enrichment of common HLA among donor TCR contributors in cluster
    score_hla <- c()
    lowest_hla <- ""
    if(hla_info == TRUE && patient_info == TRUE){
      act_seq_infos <- act_seq_infos[act_seq_infos$HLA != "" & !is.na(act_seq_infos$HLA),]

      if(nrow(act_seq_infos) > 0){
        act_seq_infos$patient <- gsub(":.*", "",act_seq_infos$patient)

      score_hla <- 1
      for(act_hla in seq_len(num_HLAs)){
        crg_patient_count <- length(unique(act_seq_infos$patient))
        crg_patient_hla_count <- sum(unlist(lapply(all_patient_hlas[unique(act_seq_infos$patient)], function(x){
          if(all_hlas$HLA[act_hla] %in% x) 1 else 0
        })))
        if(crg_patient_hla_count > 1){
          act_Prob <- sum(choose(all_hlas$counts[act_hla], crg_patient_hla_count:crg_patient_count)*choose(num_patients-all_hlas$counts[act_hla], crg_patient_count-crg_patient_hla_count:crg_patient_count)/choose(num_patients, crg_patient_count))
          if(act_Prob<score_hla) score_hla <- act_Prob
          if(act_Prob < hla_cutoff){
            if(lowest_hla == ""){
              lowest_hla <- paste(all_hlas$HLA[act_hla],
                                        " [(", crg_patient_hla_count, "/", crg_patient_count, ") vs (",
                                        all_hlas$counts[act_hla], "/", num_patients,
                                        ") = ",
                                        formatC(act_Prob, digits = 1, format = "e"),
                                        "]",
                                        sep = "" )
            } else {
              lowest_hla <- paste(lowest_hla, ", ", all_hlas$HLA[act_hla],
                                        " [(", crg_patient_hla_count, "/", crg_patient_count, ") vs (",
                                        all_hlas$counts[act_hla], "/", num_patients,
                                        ") = ",
                                        formatC(act_Prob, digits = 1, format = "e"),
                                        "]",
                                        sep = "" )
            }
          }
        }

      }
      } else {
      score_hla <- 1
      }

      all_scores <- c(all_scores, score_hla)
    }

    ### Total score
    if(gliph_version == 1) score_final <- prod(all_scores)*0.001*64 else if(gliph_version == 2) score_final <- prod(all_scores)

    ### Output
    all_scores <- c(score_final, all_scores)
    all_scores <- formatC(all_scores, digits = 1, format = "e")
    output <- c(names(cluster_list)[actCluster], all_scores)
    if(hla_info == TRUE && patient_info == TRUE){
      output <- c(output, lowest_hla)
    }
    output
  }

  doParallel::stopImplicitCluster()

  res <- data.frame(matrix(unlist(res), ncol = 2+length(score_names), byrow = TRUE))
  colnames(res) <- c("leader.tag", "total.score", score_names)

  for(i in c("total.score", score_names)) if(i != "lowest.hlas") res[,i] <- as.numeric(res[,i])

  # set all theoretical numeric values (actually characters) to numeric values
  if(is.data.frame(res)){
    for(i in seq_len(ncol(res))){
      if(suppressWarnings(any(is.na(as.numeric(res[,i])))) == FALSE) res[,i] <- as.numeric(res[,i])
    }
  }

  # Closing time!
  return(res[,-1])
}
