## ------------------------------------------------------------------
## Build reference_list.RData from the GLIPH reference files
##
## Raw reference data is downloaded from:
##   http://50.255.35.37:8080/downloads/human_v1.0.zip
##   http://50.255.35.37:8080/downloads/human_v2.0.zip
##   http://50.255.35.37:8080/downloads/mouse_v1.0.zip
##
## Original publications:
##   human_v1.0 - Glanville et al. (2017) Nature 547:94-98
##   human_v2.0 - Huang et al. (2020) Nature Biotech 38:1194-1202
##   mouse_v1.0 - Huang et al. (2020) Nature Biotech 38:1194-1202
##
## Each zip contains files for each T-cell subset (CD4, CD8, CD48):
##   ref_<subset>_*.txt  - Tab-delimited CDR3b, TRBV, TRBJ
##   ref_V_<subset>_*.txt - V-gene usage frequencies (vgene, freq)
##   ref_L_<subset>_*.txt - CDR3 length frequencies   (len, freq)
##
## The resulting reference_list is a named list where each element
## (e.g. "human_v1.0_CD4") is itself a list with three components:
##   $refseqs              - data.frame(CDR3b, TRBV)
##   $vgene_frequencies    - data.frame(vgene, freq)
##   $cdr3_length_frequencies - data.frame(len, freq)
##
## The legacy name "gliph_reference" is retained as an alias for
## "human_v1.0_CD48" (the original GLIPH publication reference).
##
## USAGE:
##   Run from the package root directory:
##     Rscript data-raw/build_reference_list.R
##
##   The output is saved to data-raw/reference_list.RData.
##   Upload to GitHub releases:
##     gh release create reference_data \
##       data-raw/reference_list.RData \
##       --title "Reference data" \
##       --notes "GLIPH reference repertoires (human v1.0, v2.0; mouse v1.0)"
## ------------------------------------------------------------------

# --- Download and extract reference zip files into a temp directory --------

base_url <- "http://50.255.35.37:8080/downloads"
zip_names <- c("human_v1.0.zip", "human_v2.0.zip", "mouse_v1.0.zip")

tmp_dir <- tempdir()

for (zf in zip_names) {
  url  <- file.path(base_url, zf)
  dest <- file.path(tmp_dir, zf)
  message("Downloading ", url, " ...")
  download.file(url, dest, mode = "wb", quiet = FALSE)
  message("Extracting ", zf, " ...")
  unzip(dest, exdir = tmp_dir)
}

# --- Define mapping from reference name to file paths within tmp_dir ------

ref_specs <- list(
  human_v1.0_CD4  = list(
    seqs  = file.path(tmp_dir, "human_v1.0", "ref_CD4_v1.0.txt"),
    vgene = file.path(tmp_dir, "human_v1.0", "ref_V_CD4_v1.0.txt"),
    len   = file.path(tmp_dir, "human_v1.0", "ref_L_CD4_v1.0.txt")
  ),
  human_v1.0_CD8  = list(
    seqs  = file.path(tmp_dir, "human_v1.0", "ref_CD8_v1.0.txt"),
    vgene = file.path(tmp_dir, "human_v1.0", "ref_V_CD8_v1.0.txt"),
    len   = file.path(tmp_dir, "human_v1.0", "ref_L_CD8_v1.0.txt")
  ),
  human_v1.0_CD48 = list(
    seqs  = file.path(tmp_dir, "human_v1.0", "ref_CD48_v1.0.txt"),
    vgene = file.path(tmp_dir, "human_v1.0", "ref_V_CD48_v1.0.txt"),
    len   = file.path(tmp_dir, "human_v1.0", "ref_L_CD48_v1.0.txt")
  ),
  human_v2.0_CD4  = list(
    seqs  = file.path(tmp_dir, "human_v2.0", "ref_CD4_v2.0.txt"),
    vgene = file.path(tmp_dir, "human_v2.0", "ref_V_CD4_v2.0.txt"),
    len   = file.path(tmp_dir, "human_v2.0", "ref_L_CD4_v2.0.txt")
  ),
  human_v2.0_CD8  = list(
    seqs  = file.path(tmp_dir, "human_v2.0", "ref_CD8_v2.0.txt"),
    vgene = file.path(tmp_dir, "human_v2.0", "ref_V_CD8_v2.0.txt"),
    len   = file.path(tmp_dir, "human_v2.0", "ref_L_CD8_v2.0.txt")
  ),
  human_v2.0_CD48 = list(
    seqs  = file.path(tmp_dir, "human_v2.0", "ref_CD48_v2.0.txt"),
    vgene = file.path(tmp_dir, "human_v2.0", "ref_V_CD48_v2.0.txt"),
    len   = file.path(tmp_dir, "human_v2.0", "ref_L_CD48_v2.0.txt")
  ),
  mouse_v1.0_CD4  = list(
    seqs  = file.path(tmp_dir, "mouse_v1.0", "ref_CD4_ms.txt"),
    vgene = file.path(tmp_dir, "mouse_v1.0", "ref_V_CD4_ms.txt"),
    len   = file.path(tmp_dir, "mouse_v1.0", "ref_L_CD4_ms.txt")
  ),
  mouse_v1.0_CD8  = list(
    seqs  = file.path(tmp_dir, "mouse_v1.0", "ref_CD8_ms.txt"),
    vgene = file.path(tmp_dir, "mouse_v1.0", "ref_V_CD8_ms.txt"),
    len   = file.path(tmp_dir, "mouse_v1.0", "ref_L_CD8_ms.txt")
  ),
  mouse_v1.0_CD48 = list(
    seqs  = file.path(tmp_dir, "mouse_v1.0", "ref_CD48_ms.txt"),
    vgene = file.path(tmp_dir, "mouse_v1.0", "ref_V_CD48_ms.txt"),
    len   = file.path(tmp_dir, "mouse_v1.0", "ref_L_CD48_ms.txt")
  )
)

# --- Build each reference entry -------------------------------------------

reference_list <- list()

for (ref_name in names(ref_specs)) {
  spec <- ref_specs[[ref_name]]
  message("Building ", ref_name, " ...")

  # --- CDR3b + TRBV sequences ---
  refseqs <- read.delim(
    spec$seqs,
    header           = FALSE,
    stringsAsFactors = FALSE
  )
  refseqs <- refseqs[, 1:2, drop = FALSE]
  colnames(refseqs) <- c("CDR3b", "TRBV")

  # --- V-gene frequencies ---
  vgene_freq <- read.delim(
    spec$vgene,
    header           = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(vgene_freq) <- c("vgene", "freq")
  vgene_freq$freq <- as.numeric(vgene_freq$freq)

  # --- CDR3 length frequencies ---
  len_freq <- read.delim(
    spec$len,
    header           = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(len_freq) <- c("len", "freq")
  len_freq$len  <- as.integer(len_freq$len)
  len_freq$freq <- as.numeric(len_freq$freq)

  reference_list[[ref_name]] <- list(
    refseqs                 = refseqs,
    vgene_frequencies       = vgene_freq,
    cdr3_length_frequencies = len_freq
  )
}

# Legacy alias: "gliph_reference" points to the original GLIPH v1.0 CD4+CD8
reference_list[["gliph_reference"]] <- reference_list[["human_v1.0_CD48"]]

# --- Save -----------------------------------------------------------------

message("Saving reference_list.RData ...")
save(reference_list, file = "data-raw/reference_list.RData", compress = "xz")
message("Done. reference_list contains: ",
        paste(names(reference_list), collapse = ", "))
message("\nTo upload to GitHub releases:")
message('  gh release create reference_data data-raw/reference_list.RData \\')
message('    --title "Reference data" \\')
message('    --notes "GLIPH reference repertoires (human v1.0, v2.0; mouse v1.0)"')
