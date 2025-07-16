# ========================================
# IMC Data Analysis Script
# ========================================

# ------------------------------
# Load Required Libraries
# ------------------------------
library(imcRtools)
library(cytomapper)
library(ggsci)
library(paletteer)
library(dittoSeq)

# ------------------------------
# Set Working Directory and Source Utilities
# ------------------------------

# Define directories for saving results and Steinbock outputs
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"
saveDir <- file.path(workDir, "results")
figDir <- file.path(workDir, "figures")
steinResultDir <- file.path(workDir, "steinbock")

source(file.path(workDir,"code","utils.R"))

# Create the results directory if it doesn't exist
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}

if (!dir.exists(figDir)) {
  dir.create(figDir, recursive = TRUE)
}

# ------------------------------
# Load and Process Data
# ------------------------------

# Load data using Steinbock
spe <- read_steinbock(steinResultDir)

# Clean and rename feature names
extracted_protein_names <- sub("^[0-9]+[A-Za-z]+_", "", rowData(spe)$name)
rowData(spe)$Clean_Target <- extracted_protein_names
rownames(spe) <- rowData(spe)$Clean_Target

# ------------------------------
# Extract Patient IDs and Annotate Metadata
# ------------------------------
head(colData(spe))

sample_id_ <- unname(sapply(spe$sample_id,function(x){ ## ROI name
  name_ <- strsplit(x,"_")[[1]]
  return(paste0(name_[2],"_",name_[3]))
})) 
patient_id_ <- unname(sapply(sample_id_, function(x) { ## Patient ID
  return(strsplit(x[1],"_")[[1]][1])
}))
cell_id_ <- paste0(sample_id_,"_",spe$ObjectNumber)

spe$sample_id <- sample_id_
spe$patient_id <- patient_id_
rownames(colData(spe)) <- cell_id_

# ------------------------------
# Cell Quality Control
# ------------------------------
summary(spe$area)
spe <- spe[, spe$area >= as.numeric(quantile(spe$area,.01))]
spe <- spe[, spe$area <= as.numeric(quantile(spe$area,.99))]

# ------------------------------
# Normalize Data
# ------------------------------
# Flag channels to be used for normalization (exclude DNA and Histone channels)
rowData(spe)$use_channel <- !grepl("DNA|Histone", rownames(spe))

# Rename assay to 'exprs' for consistency
names(assays(spe)) <- "exprs"

# Perform normalization
normalized_spe <- spe
assay(normalized_spe, "exprs") <- asinh(assay(normalized_spe)/1)

# normalized_spe <- normData(
#   spe_ = spe,
#   censor_val = 0.99,
#   arcsinh = TRUE,
#   norm_method = "0-1",
#   savePath = figDir
# )

# ------------------------------
# Set color
# ------------------------------
sampleid_temp <- sapply(normalized_spe$sample_id,function(x){
  return(extract_sample_info(x))
})

cat("The number of patients was",length(names(table(normalized_spe$patient_id))),"\n")

## set color bar
colors_dittoColors <-  dittoColors()
colors_20 <- paletteer_dynamic("cartography::pastel.pal", 20)
colors_10 <- pal_npg("nrc")(10)

## combine to object
color_vectors <- list()
color_vectors$color_multi <- colors_dittoColors
color_vectors$color_20 <- colors_20
color_vectors$color_10 <- colors_10

metadata(normalized_spe)$color_vectors <- color_vectors

# ------------------------------
# Combine Data
# ------------------------------
# (read meta data and combine here)
sce <- readRDS("/mnt/public/lyx/IMC_HE_Merge/CRLM/all_sce_2025.rds")
colnames(colData(sce))

## Patient Information
features <- c(
  ## Basic information
  "Batch",
  "Chemical_therapy","taget_therapy",
  
  ## Survival
  "RFS_status","RFS_time",
  
  ## Contineous
  "fong_score","Age","TBS","CEA","CA199","CRLM_number","CRLM_size","T_stage",
  
  ## Category
  "Gender","Pathology","Recurrence_site",
  "Differential_grade","Lymph_positive",
  "KRAS_mutation","NAS_mutation","BRAF_mutation"
)

metadata_df <- as.data.frame(colData(sce))

for(feature_ in features){
  colData(normalized_spe)[,feature_] <- metadata_df[match(normalized_spe$patient_id,sce$PID),feature_]
}

## ROI information
roi_df <- read.csv(file.path(steinResultDir,"region type.csv"))
normalized_spe$Tissue <- roi_df[match(normalized_spe$sample_id,roi_df$Sample_id),]$tissue

# ------------------------------
# Save Normalized Data
# ------------------------------
saveRDS(normalized_spe, file.path(saveDir, "norm_spe.rds"))
