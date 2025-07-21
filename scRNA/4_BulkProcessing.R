library(GEOquery)
library(dplyr)

dataPath <- "/mnt/NAS_21T/ProjectData/IMC_CRLM/bulkRNA"

## GSE41258
GEO_id <- "GSE41258"

gse <- getGEO(filename = file.path(dataPath, paste0(GEO_id,"_series_matrix.txt.gz")), getGPL = TRUE, destdir = dataPath)

## Extract the expression matrix and sample information
exprSet <- exprs(gse)
sampleInfo <- pData(gse)


# Save the cleaned expression matrix
write.csv(final_expression, "GSE41258_final_expression.csv")

# Save the clinical/sample information
write.csv(sampleInfo, "GSE41258_clinical_data.csv")

## BCGSC
study_name <- "BCGSC"
exp_data <- read.csv(file.path(dataPath, study_name, "data_mrna_seq_rpkm_zscores_ref_all_samples.txt"),header = TRUE, sep = "\t")
clinical_data <- read.table(file.path(dataPath, study_name, "pog570_bcgsc_2020_clinical_data.tsv"), header = TRUE, sep = "\t")

# Process the expression data
gene_list <- read.table(file.path(dataPath, study_name, "protein-coding_gene.txt"), header = TRUE, sep = "\t") ## Load gene reference

exp_data <- exp_data[exp_data$Hugo_Symbol %in% gene_list$symbol, ]  # Filter for protein-coding genes
exp_data_aggregated <- exp_data %>% ## Aggregate expression data by gene symbol
  group_by(Hugo_Symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()
rownames(exp_data_aggregated) <- exp_data_aggregated$Hugo_Symbol
final_expression <- exp_data_aggregated[, -1]

# Match the clinical data with the expression data
head(clinical_data)
colnames(clinical_data) <- gsub("\\.", "_", colnames(clinical_data))

clinical_data$Sample_ID <- paste0("X", clinical_data$Sample_ID)
final_expression <- final_expression[, colnames(final_expression) %in% clinical_data$Sample_ID]

# Save the cleaned expression matrix
write.csv(final_expression, file.path(dataPath, study_name, "exp_data.csv"))

# Save the clinical/sample information
write.csv(clinical_data, file.path(dataPath, study_name, "clinical_data.csv"))