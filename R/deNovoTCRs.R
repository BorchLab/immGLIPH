#' Generate de novo TCR sequences
#'
#' De novo generation of CDR3 sequences based on GLIPH or GLIPH2 clustering
#' results. Using the position-specific abundance of amino acids in the CDR3
#' region of sequences within a convergence group, artificial sequences are
#' simulated following the approach established in Glanville et al. The
#' generated sequences are scored by a positional weight matrix (PWM) derived
#' from the convergence group, and optionally normalized against a reference
#' database. The top-scoring sequences are returned.
#'
#' @param convergence_group_tag Character. Tag of the convergence group to
#' use for prediction.
#' @param result_folder Character. Path to the folder containing clustering
#' output files and where results will be saved. If the value is \code{""},
#' results are not saved to disk and the clustering output must be provided
#' via \code{clustering_output}. **Default:** \code{""}
#' @param clustering_output List. The output list from \code{turbo_gliph} or
#' \code{gliph2}. Required when \code{result_folder} is \code{""}.
#' **Default:** \code{NULL}
#' @param refdb_beta Character or data.frame. Specifies the reference
#' database to use. When a data.frame is provided, the first column should
#' contain CDR3b sequences and the second column (optional) should contain
#' V genes. The following keyword can be used to select a built-in
#' database:
#' \itemize{
#' \item \code{"gliph_reference"}: 162,165 CDR3b sequences of naive human
#' CD4+ or CD8+ T cells from two individuals (GLIPH paper).
#' }
#' **Default:** \code{"gliph_reference"}
#' @param normalization Logical. If \code{TRUE}, calculated scores are
#' normalized to the reference database. The returned value represents the
#' probability that a reference sequence has a score greater than or equal
#' to the sample sequence score. When V gene information is available, only
#' sequences with identical V genes are compared.
#' **Default:** \code{FALSE}
#' @param accept_sequences_with_C_F_start_end Logical. If \code{TRUE}, only
#' sequences beginning with cysteine (C) and ending with phenylalanine (F)
#' are accepted. **Default:** \code{TRUE}
#' @param sims Numeric. Number of de novo CDR3 sequences to generate.
#' **Default:** \code{100000}
#' @param num_tops Numeric. Number of top-scoring de novo sequences to
#' return. **Default:** \code{1000}
#' @param min_length Numeric. Minimum CDR3 sequence length; also determines
#' the number of N-terminal positions used for PWM scoring.
#' **Default:** \code{10}
#' @param make_figure Logical. Whether to plot the \code{num_tops}
#' best-scoring de novo sequences as a function of rank.
#' **Default:** \code{FALSE}
#' @param n_cores Numeric. Number of cores for parallel computation. If
#' \code{NULL}, the number of available cores minus two is used.
#' **Default:** \code{1}
#'
#' @return A \code{list} with the following elements:
#' \describe{
#'   \item{de_novo_sequences}{A \code{data.frame} of the \code{num_tops}
#'   best-scoring generated sequences and their corresponding scores.}
#'   \item{sample_sequences_scores}{A \code{data.frame} of the convergence
#'   group sequences and their corresponding scores.}
#'   \item{cdr3_length_probability}{A \code{data.frame} with each observed
#'   CDR3 length and its probability of occurrence in the convergence group.
#'   The length distribution of generated sequences mirrors this
#'   distribution.}
#'   \item{PWM_Scoring}{A \code{data.frame} containing the positional weight
#'   matrix used for scoring. Columns represent amino acids and rows
#'   represent positions relative to the N-terminus.}
#'   \item{PWM_Prediction}{A \code{list} of \code{data.frame}s containing
#'   the positional weight matrices used for sequence generation, one per
#'   observed CDR3 length. Columns represent amino acids and rows represent
#'   positions relative to the N-terminus.}
#' }
#' If \code{result_folder} is specified, a tab-delimited file named
#' \code{<convergence_group_tag>_de_novo.txt} is also written to disk.
#'
#' @examples
#' utils::data("gliph_input_data")
#' res <- turbo_gliph(cdr3_sequences = gliph_input_data[seq_len(200),],
#' sim_depth = 100,
#' n_cores = 1)
#'
#' new_seqs <- deNovoTCRs(convergence_group_tag = res$cluster_properties$tag[1],
#' clustering_output = res,
#' sims = 10000,
#' make_figure = TRUE,
#' n_cores = 1)
#'
#' @references Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
#' @references https://github.com/immunoengineer/gliph
#' @import foreach
#' @export
deNovoTCRs <- function(convergence_group_tag,
                         result_folder = "",
                         clustering_output = NULL,
                         refdb_beta = "gliph_reference",
                         normalization = FALSE,
                         accept_sequences_with_C_F_start_end = TRUE,
                         sims = 100000,
                         num_tops = 1000,
                         min_length = 10,
                         make_figure = FALSE,
                         n_cores = 1){
  t1 <- Sys.time()

  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################

  ### convergence_group_tag
  if(!is.character(convergence_group_tag)) stop("convergence_group_tag has to be a character object")
  if(length(convergence_group_tag) > 1) stop("convergence_group_tag has to be a single character string")

  ### result_folder and clustering_output
  if(!is.character(result_folder)) stop("result_folder has to be a character object")
  if(length(result_folder) > 1) stop("result_folder has to be a single path")
  save_results <- FALSE
  if(result_folder != ""){
    if(substr(result_folder,nchar(result_folder),nchar(result_folder)) != "/") result_folder <- paste0(result_folder,"/")
    if (!dir.exists(result_folder)) dir.create(result_folder)
    save_results <- TRUE
  } else {
    if(!is.list(clustering_output)) stop("If 'result_folder' = \"\" the output list of clustering must be given by 'clustering_output'.")
  }
  ### refdb_beta
  # if(!(refdb_beta %in% c("gliph_reference", "human_v1.0_CD4", "human_v1.0_CD8", "human_v1.0_CD48", "human_v2.0_CD4",
  #                        "human_v2.0_CD8", "human_v2.0_CD48", "mouse_v1.0_CD4", "mouse_v1.0_CD8", "mouse_v1.0_CD48")) &&
  #    !is.data.frame(refdb_beta)){
  if(!is.data.frame(refdb_beta)){
    if(length(refdb_beta) != 1 || !is.character(refdb_beta)){
      stop("refdb_beta has to be a data frame (containing CDR3b sequences in the first column and optional V-gene information in the second column) or the value 'gliph_reference'")
    } else if(!(refdb_beta %in% c("gliph_reference"))){
      stop("refdb_beta has to be a data frame (containing CDR3b sequences in the first column and optional V-gene information in the second column) or the value 'gliph_reference'")
    }
  }

  ### accept_sequences_with_C_F_start_end
  if(!is.logical(accept_sequences_with_C_F_start_end)) stop("accept_sequences_with_C_F_start_end has to be logical")

  ### normalization
  if(!is.logical(normalization)) stop("normalization has to be logical")

  ### sims
  if(!is.numeric(sims)) stop("sims has to be numeric")
  if(length(sims) > 1) stop("sims has to be a single number")
  if(sims < 1) stop("sims must be at least 1")
  sims <- round(sims)

  ### num_tops
  if(!is.numeric(num_tops)) stop("num_tops has to be numeric")
  if(length(num_tops) > 1) stop("num_tops has to be a single number")
  if(num_tops < 1) stop("num_tops must be at least 1")
  num_tops <- round(num_tops)

  ### min_length
  if(!is.numeric(min_length)) stop("min_length has to be numeric")
  if(length(min_length) > 1) stop("min_length has to be a single number")
  if(min_length < 1) stop("min_length must be at least 1")
  min_length <- round(min_length)

  ### make_figure
  if(!is.logical(make_figure)) stop("make_figure has to be logical")

  ### n_cores
  if(is.null(n_cores)) n_cores <- max(1, parallel::detectCores()-2) else {
    if(!is.numeric(n_cores)) stop("n_cores has to be numeric")
    if(length(n_cores) > 1) stop("n_cores has to be a single number")
    if(n_cores < 1) stop("n_cores must be at least 1")

    n_cores <- min(n_cores, parallel::detectCores()-2)
  }

  # Amino acids one letter code
  aa_code <- LETTERS[-c(2, 10, 15, 21, 24, 26)]

  #################################################################
  ##                      Input preparation                      ##
  #################################################################

  ### load convergence groups from file or from input and only use specified convergence group
  message("Loading convergence group with tag ", convergence_group_tag, ".")
  if(is.null(clustering_output)){
    clustering_output <- loadGLIPH(result_folder = result_folder)
    crg <- clustering_output$cluster_list
  } else {
    clustering_output <- clustering_output
    crg <- clustering_output$cluster_list
  }
  if(!(convergence_group_tag %in% names(crg))) stop("Could not find convergence group with tag ",convergence_group_tag, " in clustering results stored in", result_folder, ".")
  crg <- crg[[convergence_group_tag]]
  all_crg_cdr3_seqs <- crg$CDR3b

  ### Filter sequences for a minimun sequenc length
  excluded <- which(nchar(all_crg_cdr3_seqs) < min_length)
  if(length(excluded) > 0){
    all_crg_cdr3_seqs <- all_crg_cdr3_seqs[-excluded]
    if(length(all_crg_cdr3_seqs) > 0){
      message("Warning: ", length(excluded), " sequences of the convergence group were excluded from the further procedure due to falling below a minimum length of ", min_length, ".")
    } else {
      stop("No sequences of the convergence group are of minimum length of ", min_length, ". For further procedure, adjust the parameter 'min_length'")
    }
  }

  #################################################################
  ##                          Main part                          ##
  #################################################################

  ### Initiate parallelization
  doParallel::registerDoParallel(n_cores)

  # get all sequences in cluster
  crg_cdr3_seqs <- all_crg_cdr3_seqs

  ### load non-redundant reference repertoire and V-genes of sequences in convergence group
  v_genes <- c()
  ref_vgenes <- c()
  refseqs <- c()
  v_gene_norm <- normalization
  if(normalization == TRUE){

    if(is.data.frame(refdb_beta)) {
      refseqs <- refdb_beta
      refseqs[] <- lapply(refseqs, as.character)

      if(ncol(refseqs) > 1){
        message("Notification: First column of reference database is considered as cdr3 sequences.")
      }
      if(ncol(refseqs) > 1 && v_gene_norm == TRUE){
        message("Notification: Second column of reference database is considered as V-gene information.")
      } else if(v_gene_norm == TRUE){
        message("Warning: Beta sequence reference database is missing column containing V-genes. Without V-gene information normalization may be inaccurate.")
        v_gene_norm <- FALSE
      }
      if(ncol(refseqs) == 1) refseqs <- cbind(refseqs, rep("", nrow(refseqs)))

      refseqs <- refseqs[, c(1,2)]
      colnames(refseqs) <- c("CDR3b", "TRBV")
      refseqs <- unique(refseqs)
      if(accept_sequences_with_C_F_start_end) refseqs <- refseqs[grep(pattern = "^C.*F$",x = refseqs$CDR3b,perl = TRUE),]
      refseqs <- refseqs[which(nchar(refseqs$CDR3b) >= min_length),]
      refseqs <- refseqs[grep("^[ACDEFGHIKLMNOPQRSTUVWY]*$", refseqs$CDR3b),]

      if(nrow(refseqs) == 0){
        normalization <- FALSE
        v_gene_norm <- FALSE
        message("Warning: No reference sequences with a minimum length of ", min_length, " given. Normalization therefore not possible. Adjust min_length to enable normalization.")
      } else {
        ref_vgenes <- as.character(refseqs$TRBV)
        refseqs <- as.character(refseqs$CDR3b)
      }
    } else {
      reference_list <- NULL
      utils::data("reference_list",envir = environment(), package = "immGLIPH")
      refseqs <- as.data.frame(reference_list[[refdb_beta]]$refseqs)
      refseqs[] <- lapply(refseqs, as.character)

      if(ncol(refseqs) > 1){
        message("Notification: First column of reference database is considered as cdr3 sequences.")
      }
      if(ncol(refseqs) > 1 && v_gene_norm == TRUE){
        message("Notification: Second column of reference database is considered as V-gene information.")
      } else if(v_gene_norm == TRUE){
        message("Warning: Beta sequence reference database is missing column containing V-genes. Without V-gene information normalization may be inaccurate.")
        v_gene_norm <- FALSE
      }
      if(ncol(refseqs) == 1) refseqs <- cbind(refseqs, rep("", nrow(refseqs)))

      refseqs <- refseqs[, c(1,2)]
      colnames(refseqs) <- c("CDR3b", "TRBV")
      refseqs <- unique(refseqs)
      if(accept_sequences_with_C_F_start_end) refseqs <- refseqs[grep(pattern = "^C.*F$",x = refseqs$CDR3b,perl = TRUE),]
      refseqs <- refseqs[which(nchar(refseqs$CDR3b) >= min_length),]
      refseqs <- refseqs[grep("^[ACDEFGHIKLMNPQRSTVWY]*$", refseqs$CDR3b),]

      if(nrow(refseqs) == 0){
        normalization <- FALSE
        v_gene_norm <- FALSE
        message("Warning: No reference sequences with a minimum length of ", min_length, " given. Normalization therefore not possible. Adjust min_length to enable normalization.")
      } else {
        ref_vgenes <- as.character(refseqs$TRBV)
        refseqs <- as.character(refseqs$CDR3b)
      }
    }

    # load V genes of cluster members
    if("TRBV" %in% colnames(crg) && v_gene_norm == TRUE){
      v_genes <- crg$TRBV
    } else if(v_gene_norm == TRUE){
      v_gene_norm <- FALSE
      message("Warning: Without V-gene information of sample sequences normalization may be inaccurate.")
    } else message("Warning: Without V-gene restriction normalization may be inaccurate.")
  }

  ### Initialization
  crg_num_seqs <- length(crg_cdr3_seqs)
  crg_cdr3_scores <- rep(1, crg_num_seqs)
  crg_cdr3_norm_scores <- rep(0, crg_num_seqs)
  crg_cdr3_lens <- nchar(crg_cdr3_seqs)
  max_cdr3_length <- max(crg_cdr3_lens)

  ### Build Positional Weight Matrix (PWM) of convergence group for min_length N-terminal positions
  message("Calculating positional weight matrix of convergence group.")

  # initialize the matrix
  crg_pwm_scoring <- as.data.frame(matrix(rep(0, min_length*length(aa_code)), ncol = length(aa_code)))
  colnames(crg_pwm_scoring) <- aa_code

  # for every position determine the amino acid frequency
  for(i in seq_len(nrow(crg_pwm_scoring))){
    aa_freqs <- rep(0, length(aa_code))
    act_letters <- substr(crg_cdr3_seqs, i, i)
    for(j in seq_along(aa_freqs)){
      aa_freqs[j] <- sum(act_letters == aa_code[j])
    }

    # Pseudocounts of 0.5% per aa
    zeroFreqs <- which(aa_freqs == 0)
    aa_freqs <- aa_freqs/sum(aa_freqs)*(1-length(zeroFreqs)*0.005)
    aa_freqs[zeroFreqs] <- 0.005
    crg_pwm_scoring[i,] <- aa_freqs
  }

  ### Calculate scores of convergence group members, only use first min_length N-terminal positions
  # score = product of amino acid frequencies from min_length N-terminal positions
  message("Calculating scores of convergence group members.")

  # is parallelization necessary?
  if(crg_num_seqs > 10000){
    # distribute sequences equally to all cores
    distribute <- lapply(seq_len(n_cores), function(x){return(((x-1)*floor(length(crg_cdr3_seqs)/n_cores)+1):(x*floor(length(crg_cdr3_seqs)/n_cores)))})

    # calculate the score
    crg_cdr3_scores <- foreach::foreach(i = seq_along(distribute), .combine = c) %dopar% {
      temp_scores <- rep(1, length(distribute[[i]]))
      temp_seqs <- crg_cdr3_seqs[distribute[[i]]]
      for(i in seq_len(nrow(crg_pwm_scoring))){
        temp_scores <- temp_scores*unlist(crg_pwm_scoring[i, substr(temp_seqs, i, i)])
      }

      return(temp_scores)
    }
  } else {
    # calculate the score
    for(i in seq_len(nrow(crg_pwm_scoring))){
      crg_cdr3_scores <- crg_cdr3_scores*unlist(crg_pwm_scoring[i, substr(crg_cdr3_seqs, i, i)])
    }
  }

  ### Normalize scores of convergence group members, only use first 10 N-terminal positions
  # calculate the probability that a score at least this high occurs in the reference database
  if(normalization == TRUE){
    message("Normalizing scores of convergence group members.")

    # Calculate scores of reference database (analogue to sample sequences)
    refseq_scores <- rep(1, length(refseqs))
    if(length(refseqs) > 10000){
      distribute <- lapply(seq_len(n_cores), function(x){return(((x-1)*floor(length(refseqs)/n_cores)+1):(x*floor(length(refseqs)/n_cores)))})

      refseq_scores <- foreach::foreach(i = seq_along(distribute), .combine = c) %dopar% {
        temp_scores <- rep(1, length(distribute[[i]]))
        temp_seqs <- refseqs[distribute[[i]]]
        for(i in seq_len(nrow(crg_pwm_scoring))){
          temp_scores <- temp_scores*unlist(crg_pwm_scoring[i, substr(temp_seqs, i, i)])
        }

        return(temp_scores)
      }
    } else {
      for(i in seq_len(nrow(crg_pwm_scoring))){
        refseq_scores <- refseq_scores*unlist(crg_pwm_scoring[i, substr(refseqs, i, i)])
      }
    }

    # Calculate normalized scores, if requested restrict score comparison to identical V genes
    crg_cdr3_norm_scores <- foreach::foreach(i = seq_len(crg_num_seqs), .combine = c) %dopar% {
      v_gene_penalty <- rep(0, length(refseqs))
      if(v_gene_norm == TRUE){
        v_gene_penalty[ref_vgenes != v_genes[i]] <- -2
      }

      return(sum((refseq_scores + v_gene_penalty) >= crg_cdr3_scores[i])/length(refseqs))
    }
  }

  ### Create global PWM of convergence group
  message("Calculating positional weight matrix for de novo CDR3b-sequence prediction.")

  # create a PWM for every CDR3b length in the cluster
  crg_pwm_predicting_list <- list()

  # get the probability of this CDR3b length to occur
  crg_len_prob <- data.frame(length = unique(crg_cdr3_lens),probability = rep(0, length(unique(crg_cdr3_lens))))
  for(i in seq_len(nrow(crg_len_prob))){
    crg_len_prob$probability[i] <- sum(crg_cdr3_lens == crg_len_prob$length[i])/crg_num_seqs
  }

  # receive the positional dependent amino acid frequencies for every CDR3b length
  for(len in unique(crg_cdr3_lens)){
    crg_pwm_predicting <- as.data.frame(matrix(rep(0,len*length(aa_code)), ncol = length(aa_code)))
    colnames(crg_pwm_predicting) <- aa_code

    seqs <- crg_cdr3_seqs[crg_cdr3_lens == len]
    for(i in seq_len(len)){
      aa_freqs <- rep(0, length(aa_code))
      act_letters <- substr(seqs, i, i)
      for(j in seq_along(aa_freqs)){
        aa_freqs[j] <- sum(act_letters == aa_code[j])
      }

      # Pseudocounts of 0.5% per aa
      zeroFreqs <- which(aa_freqs == 0)
      aa_freqs <- aa_freqs/sum(aa_freqs)*(1-length(zeroFreqs)*0.005)
      aa_freqs[zeroFreqs] <- 0.005
      crg_pwm_predicting[i,] <- aa_freqs
    }

    # if de novo sequences should only start with C and end with F, correct the PWM for this positions
    if(accept_sequences_with_C_F_start_end == TRUE){
      crg_pwm_predicting[1, ] <- rep(0, ncol(crg_pwm_predicting))
      crg_pwm_predicting[nrow(crg_pwm_predicting), ] <- rep(0, ncol(crg_pwm_predicting))
      crg_pwm_predicting$'C'[1] <- 1
      crg_pwm_predicting$'F'[nrow(crg_pwm_predicting)] <- 1
    }

    crg_pwm_predicting_list[[paste("Length", len)]] <- crg_pwm_predicting
  }


  ### Create a number of sims de novo TCR sequences
  message("Creating ", sims, " de novo sequences.")

  # randomly select the length of the sequences
  de_novo_lens <- sample(x = c(crg_len_prob$length, 0), size = sims, prob = c(crg_len_prob$probability,0), replace = TRUE)
  de_novo_seqs <- rep("", sims)

  # based on the PWM randomly create new sequences
  for(len in crg_len_prob$length){
    inds <- which(de_novo_lens == len)
    for(i in seq_len(len)){
      rands <- aa_code[sample.int(n = length(aa_code), size = length(inds), prob = crg_pwm_predicting_list[[paste("Length", len)]][i,], replace = TRUE)]
      de_novo_seqs[inds] <- paste(de_novo_seqs[inds], rands, sep = "")
    }
  }
  de_novo_seqs <- unique(de_novo_seqs)
  if(accept_sequences_with_C_F_start_end) de_novo_seqs <- grep(pattern = "^C.*F$",x = de_novo_seqs ,perl = TRUE,value = TRUE)
  de_novo_lens <- nchar(de_novo_seqs)

  # Score de_novo seqs as performed above
  message("Calculating scores of de novo sequences.")
  de_novo_seqs_scores <- rep(1, length(de_novo_seqs))

  if(length(de_novo_seqs) > 10000){
    distribute <- lapply(seq_len(n_cores), function(x){return(((x-1)*floor(length(de_novo_seqs)/n_cores)+1):(x*floor(length(de_novo_seqs)/n_cores)))})

    de_novo_seqs_scores <- foreach::foreach(i = seq_along(distribute), .combine = c) %dopar% {
      temp_scores <- rep(1, length(distribute[[i]]))
      temp_seqs <- de_novo_seqs[distribute[[i]]]
      for(i in seq_len(nrow(crg_pwm_scoring))){
        temp_scores <- temp_scores*unlist(crg_pwm_scoring[i, substr(temp_seqs, i, i)])
      }

      return(temp_scores)
    }
  } else {
    for(i in seq_len(nrow(crg_pwm_scoring))){
      de_novo_seqs_scores <- de_novo_seqs_scores*unlist(crg_pwm_scoring[i, substr(de_novo_seqs, i, i)])
    }
  }
  de_novo_seqs_scores <- as.numeric(formatC(de_novo_seqs_scores, digits = 1, format = "e"))

  # Normalize de_novo seqs scores as performed above
  if(normalization == TRUE){
    message("Normalizing scores of de novo sequences.")

    de_novo_norm_scores <- foreach::foreach(i = seq_along(de_novo_seqs), .combine = c) %dopar% {
      return(sum(refseq_scores >= de_novo_seqs_scores[i])/length(refseqs))
    }
    de_novo_norm_scores <- as.numeric(formatC(de_novo_norm_scores, digits = 1, format = "e"))
  }

  # sort sequences based on score (normalized scores are prioritized)
  if(normalization == TRUE){
    order_ids <- order(de_novo_norm_scores, decreasing=FALSE)
    de_novo_seqs <- de_novo_seqs[order_ids]
    de_novo_lens <- de_novo_lens[order_ids]
    de_novo_seqs_scores <- de_novo_seqs_scores[order_ids]
    de_novo_norm_scores <- de_novo_norm_scores[order_ids]
  } else{
    order_ids <- order(de_novo_seqs_scores, decreasing=TRUE)
    de_novo_seqs <- de_novo_seqs[order_ids]
    de_novo_lens <- de_novo_lens[order_ids]
    de_novo_seqs_scores <- de_novo_seqs_scores[order_ids]
  }

  # get num_tops heighest scored sequences
  if(length(de_novo_seqs) < num_tops) num_tops <- length(de_novo_seqs)
  de_novo_seqs <- de_novo_seqs[seq_len(num_tops)]
  de_novo_lens <- de_novo_lens[seq_len(num_tops)]
  de_novo_seqs_scores <- de_novo_seqs_scores[seq_len(num_tops)]
  if(normalization == TRUE){
    de_novo_norm_scores <- de_novo_norm_scores[seq_len(num_tops)]
    de_novo <- data.frame(length = de_novo_lens, seqs = de_novo_seqs, norm_score = de_novo_norm_scores, score = de_novo_seqs_scores)
  } else {
    de_novo <- data.frame(length = de_novo_lens, seqs = de_novo_seqs, score = de_novo_seqs_scores)
  }

  ### sort cluster sequences based on their score
  connected_inds <- seq_along(crg_cdr3_seqs)
  if(normalization == FALSE){
    crg_cdr3_seqs <- crg_cdr3_seqs[order(crg_cdr3_scores, decreasing = TRUE)]
    connected_inds <- connected_inds[order(crg_cdr3_scores, decreasing = TRUE)]
    crg_cdr3_scores <- crg_cdr3_scores[order(crg_cdr3_scores, decreasing = TRUE)]
  } else{
    crg_cdr3_seqs <- crg_cdr3_seqs[order(crg_cdr3_norm_scores, decreasing = FALSE)]
    connected_inds <- connected_inds[order(crg_cdr3_norm_scores, decreasing = FALSE)]
    crg_cdr3_scores <- crg_cdr3_scores[order(crg_cdr3_norm_scores, decreasing = FALSE)]
    crg_cdr3_norm_scores <- crg_cdr3_norm_scores[order(crg_cdr3_norm_scores, decreasing = FALSE)]
  }

  ### set all theoretical numeric values (actually characters) to numeric values
  if(is.data.frame(de_novo)){
    for(i in seq_len(ncol(de_novo))){
      if(suppressWarnings(any(is.na(as.numeric(de_novo[,i])))) == FALSE) de_novo[,i] <- as.numeric(de_novo[,i])
    }
  }
  if(is.data.frame(crg_cdr3_scores)){
    for(i in seq_len(ncol(crg_cdr3_scores))){
      if(suppressWarnings(any(is.na(as.numeric(crg_cdr3_scores[,i])))) == FALSE) crg_cdr3_scores[,i] <- as.numeric(crg_cdr3_scores[,i])
    }
  }
  if(normalization == TRUE){
    if(is.data.frame(crg_cdr3_norm_scores)){
      for(i in seq_len(ncol(crg_cdr3_norm_scores))){
        if(suppressWarnings(any(is.na(as.numeric(crg_cdr3_norm_scores[,i])))) == FALSE) crg_cdr3_norm_scores[,i] <- as.numeric(crg_cdr3_norm_scores[,i])
      }
    }
  }

  ### output
  if(normalization == FALSE){
    output <- list(de_novo_sequences = de_novo, sample_sequences_scores = data.frame(seqs = all_crg_cdr3_seqs[connected_inds], scores = crg_cdr3_scores),
                         cdr3_length_probability = crg_len_prob, PWM_Scoring = crg_pwm_scoring, PWM_Prediction = crg_pwm_predicting_list)
  } else {
    output <- list(de_novo_sequences = de_novo, sample_sequences_scores = data.frame(seqs = all_crg_cdr3_seqs[connected_inds], norm_scores = crg_cdr3_norm_scores, scores = crg_cdr3_scores),
                         cdr3_length_probability = crg_len_prob, PWM_Scoring = crg_pwm_scoring, PWM_Prediction = crg_pwm_predicting_list)
  }

  ### save
  fname <- paste0(result_folder, convergence_group_tag, "_de_novo.txt")
  if(save_results == TRUE) utils::write.table(x = de_novo,file = fname,quote = FALSE,sep = "\t",row.names = FALSE, col.names = TRUE)
  if(save_results == TRUE) message("Output: results are stored in ", fname)

  ### Print a graph with num_tops best scoring de novo sequences
  if(make_figure == TRUE){
    if(normalization == FALSE){
      graphics::plot(x = seq_len(num_tops), y = de_novo$score*100, xlab = paste("Top", num_tops, "predicted TCRs"), ylab = "TCR probability based on PWM in %", type = "p", log = "xy", col = "grey", pch = 19)
      graphics::points(x = seq_len(10), y = de_novo$score[seq_len(10)]*100, col = "red", pch = 19)

      # which cluster sequences are present in de novo created sequences?
      common_ids <- which(de_novo_seqs %in% crg_cdr3_seqs)
      graphics::points(x = common_ids, y = de_novo_seqs_scores[common_ids]*100, col = "yellow", pch = 20)
    } else{
      graphics::plot(x = seq_len(num_tops), y = de_novo$norm_score *100, xlab = paste("Top", num_tops, "predicted TCRs"), ylab = "Normalized TCR probability based on PWM in %",
                    type = "p", log = "xy", col = "grey", pch = 19)
      graphics::points(x = seq_len(10), y = de_novo$norm_score[seq_len(10)]*100, col = "red", pch = 19)

      # which cluster sequences are present in de novo created sequences?
      common_ids <- which(de_novo_seqs %in% crg_cdr3_seqs)
      graphics::points(x = common_ids, y = de_novo_norm_scores[common_ids]*100, col = "yellow", pch = 20)
    }

    graphics::legend("bottomleft", legend=c(paste0("Top ", num_tops," de novo scores"), "Top 10 de novo scores", "Convergence group members\nin de novo sequences"),
           col=c("grey", "red", "yellow"), pch = c(19,19,20),title="Legend", cex = 0.75)
  }

  t2 <- Sys.time()
  dt <- (t2-t1)
  message("Total time = ", dt, " ", units(dt))

  doParallel::stopImplicitCluster()


  ### Closing time!
  return(output)
}
