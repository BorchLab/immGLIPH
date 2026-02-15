## Build example data objects for immGLIPH
##
## Creates:
##   data/gliph_input_data.rdata  - data.frame of TRB CDR3b, TRBV, patient
##   data/gliph_sce.rdata         - SingleCellExperiment with TCR in colData
##
## Source: scRepertoire example data (Yost et al. 2021, PMID 33622974)
##   4 patients with acute respiratory distress, lung + peripheral blood
##
## Requirements: scRepertoire, Seurat, immApex, SingleCellExperiment

library(scRepertoire)
library(Seurat)

data("contig_list", package = "scRepertoire")
data("scRep_example", package = "scRepertoire")

# --- Build combined TCR data from 8 samples ---------------------------------
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))

# --- Create SingleCellExperiment with TCR in colData -------------------------
sce_seurat <- combineExpression(combined, scRep_example, cloneCall = "aa")
gliph_sce  <- as.SingleCellExperiment(sce_seurat)

save(gliph_sce, file = "data/gliph_sce.rdata", compress = "xz")
message("Saved data/gliph_sce.rdata")

# --- Build gliph_input_data via immApex::getIR --------------------------------
ir  <- immApex::getIR(gliph_sce, chains = "TRB", sequence.type = "aa")
trb <- ir$TRB

gliph_input_data <- data.frame(
    CDR3b   = trb$cdr3_aa,
    TRBV    = trb$v,
    patient = gsub("_.*", "", trb$barcode),
    stringsAsFactors = FALSE
)

# Remove cells without TCR data
gliph_input_data <- gliph_input_data[!is.na(gliph_input_data$CDR3b), ]
rownames(gliph_input_data) <- NULL

save(gliph_input_data, file = "data/gliph_input_data.rdata", compress = "xz")
message("Saved data/gliph_input_data.rdata")

message(sprintf("gliph_input_data: %d rows, %d cols",
                nrow(gliph_input_data), ncol(gliph_input_data)))
message(sprintf("gliph_sce: %d genes x %d cells",
                nrow(gliph_sce), ncol(gliph_sce)))
