library(GEOquery)
library(dplyr)

dataPath <- "/mnt/NAS_21T/ProjectData/IMC_CRLM/bulkRNA"

# =============================================================================
## GSE41258
# =============================================================================

GEO_id <- "GSE41258"

gse <- getGEO(filename = file.path(dataPath, paste0(GEO_id,"_series_matrix.txt.gz")), getGPL = TRUE, destdir = dataPath)

## Extract the expression matrix and sample information
exprSet <- exprs(gse)
sampleInfo <- pData(gse)

# Save the cleaned expression matrix
write.csv(final_expression, "GSE41258_final_expression.csv")

# Save the clinical/sample information
write.csv(sampleInfo, "GSE41258_clinical_data.csv")

# =============================================================================
## BCGSC
# =============================================================================

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

# =============================================================================
## E-TABM-1112
# =============================================================================
library(tibble)
library(biomaRt)
library(dplyr)

# Set data path
dataPath <- "/mnt/NAS_21T/ProjectData/IMC_CRLM/bulkRNA/E-TABM-1112"
setwd(dataPath)

# Load the SDRF file
sdrf_file <- "E-TABM-1112.sdrf.txt"
sdrf <- read.delim(sdrf_file, stringsAsFactors = FALSE)
colnames(sdrf) <- gsub( "\\.", "_", colnames(sdrf))  # Replace dots with underscores in column names

# --- Create the Clinical Table ---

# Filter out the universal reference sample and select relevant columns
# We only want rows corresponding to the actual patient samples
clinical_table <- sdrf %>%
  filter(Factor_Value__DISEASE_STATE_ != "Reference") %>%
  dplyr::select(
    SampleName = Source_Name, 
    DiseaseState = Characteristics__disease_state_,
    OrganismPart = Characteristics__OrganismPart_
  ) %>%
  distinct() %>% # Remove duplicate rows for the same sample
  column_to_rownames(var = "SampleName") # Set sample names as row names

# View the first few rows of the clinical table
print(head(clinical_table))

# --- Create the Expression Matrix ---

# Get a list of unique hybridizations (i.e., unique microarray slides)
unique_hybridizations <- unique(sdrf$Hybridization_Name)

# Create a list to store the expression data for each sample
expression_list <- list()

# Loop through each unique slide
for (hyb in unique_hybridizations) {
  
  # Get the SDRF info for the current slide
  slide_info <- sdrf %>% filter(Hybridization_Name == hyb)
  
  # Find the row for the patient sample and the reference sample
  sample_row <- slide_info %>% filter(Factor_Value__DISEASE_STATE_ != "Reference")
  ref_row <- slide_info %>% filter(Factor_Value__DISEASE_STATE_ == "Reference")
  
  # Get the sample name and the data file name
  sample_name <- sample_row$Source_Name
  data_file <- sample_row$Derived_Array_Data_File
  
  # Read the derived data file for the slide
  # Using 'check.names = FALSE' to handle special characters if any
  expression_data <- read.delim(data_file, check.names = FALSE)
  
  # Identify which channel (Cy3 or Cy5) corresponds to the sample and reference
  sample_channel <- ifelse(grepl("Cy3", sample_row$Label), "Signal Norm_Cy3", "Signal Norm_Cy5")
  ref_channel <- ifelse(grepl("Cy3", ref_row$Label), "Signal Norm_Cy3", "Signal Norm_Cy5")
  
  # Calculate the log2 ratio (Sample / Reference)
  # Adding a small offset (1) to avoid log2(0)
  log_ratios <- log2((expression_data[[sample_channel]] + 1) / (expression_data[[ref_channel]] + 1))
  
  # Add the result to our list, with probe IDs as names
  names(log_ratios) <- expression_data$`Reporter Identifier`
  expression_list[[sample_name]] <- log_ratios
}

# Combine the list of vectors into a single matrix
probe_expression_matrix <- do.call(cbind, expression_list)

# View the dimensions and first few rows/columns of the matrix
print(dim(probe_expression_matrix))
print(probe_expression_matrix[1:5, 1:5])

# --- Annotate Probes to Gene Symbols ---

# --- 1. Load the local annotation file (ADF) ---

# Define the path to your annotation file
adf_file_path <- "A-UMCU-3.adf.txt" # Make sure this file is in your working directory

# Read all lines to find where the main data table starts
adf_lines <- readLines(adf_file_path)
start_line <- which(grepl("\\[main\\]", adf_lines))

# Read the data table, skipping the header metadata
probe_annotation_data <- read.delim(
  adf_file_path,
  skip = start_line,
  stringsAsFactors = FALSE
)

# Create a clean map from probe ID to Ensembl gene ID
# We only need two columns and will remove rows where the Ensembl ID is missing.
probe_to_ensembl_map <- probe_annotation_data %>%
  dplyr::select(
    probe_id = Reporter.Name,
    ensembl_gene_id = Reporter.Database.Entry.ensembl.
  ) %>%
  filter(ensembl_gene_id != "") %>%
  distinct() # Ensure unique mappings

print("Probe to Ensembl ID Map (first 6 rows):")
probe_to_ensembl_map$probe_id <- toupper(probe_to_ensembl_map$probe_id)
print(head(probe_to_ensembl_map))

# --- 2. Convert Ensembl IDs to Gene Symbols using biomaRt ---

# Get the unique list of Ensembl IDs from our map
unique_ensembl_ids <- unique(probe_to_ensembl_map$ensembl_gene_id)

# Connect to Ensembl
ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Create a map from Ensembl ID to HGNC gene symbol
ensembl_to_symbol_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique_ensembl_ids,
  mart = ensembl_mart
)

print("Ensembl to Gene Symbol Map (first 6 rows):")
print(head(ensembl_to_symbol_map))

# --- 3. Merge with Expression Data and Aggregate ---

# Convert the probe expression matrix to a data frame for joining
probe_exp_df <- as.data.frame(probe_expression_matrix) %>%
  rownames_to_column(var = "probe_id")

# Join everything together: Expression Data -> Probe Map -> Symbol Map
# inner_join automatically keeps only the rows that match in all tables.
final_expression_matrix <- probe_exp_df %>%
  inner_join(probe_to_ensembl_map, by = "probe_id") %>%
  inner_join(ensembl_to_symbol_map, by = "ensembl_gene_id") %>%
  filter(hgnc_symbol != "") %>% # Remove any mappings that resulted in an empty gene symbol
  dplyr::select(-probe_id, -ensembl_gene_id) %>% # Remove ID columns
  group_by(hgnc_symbol) %>% # Group by the final gene symbol
  summarise_all(mean) %>% # Take the mean for duplicate genes
  column_to_rownames(var = "hgnc_symbol") # Set gene symbols as row names

# View the final expression matrix with gene symbols as rows
print("Final Aggregated Expression Matrix (first 5 rows/cols):")
print(dim(final_expression_matrix))
print(final_expression_matrix[1:5, 1:5])

# --- Save the Final Expression Matrix and Clinical Data ---
write.csv(final_expression_matrix, "exp_data.csv")
write.csv(clinical_table, "clinical_data.csv")

# =============================================================================
## E-MTAB-1951
# =============================================================================
# Set data path
dataPath <- "/mnt/NAS_21T/ProjectData/IMC_CRLM/bulkRNA/E-MTAB-1951"
setwd(dataPath)

# Load required libraries
library(dplyr)
library(readr)
library(stringr)

# Set file paths (adjust these to your actual file paths)
expression_file <- "HI_2008_96samples_Sample_Probe_Profile_for_EMBL-EBI_update_101213.txt"
sdrf_file <- "E-MTAB-1951.sdrf.txt"
annotation_file <- "A-MEXP-1171.adf.txt"  # Array annotation file

## Process the expression data
lines <- readLines(expression_file) # Read the file and find where the actual data starts

# Find the line that starts with "PROBE_ID" (header of data section)
data_start <- which(str_detect(lines, "^PROBE_ID"))

# Read the data starting from the header line
expr_data <- read_tsv(expression_file, skip = data_start - 1, show_col_types = FALSE)

# Extract sample IDs from column names
sample_cols <- colnames(expr_data)[str_detect(colnames(expr_data), "\\d+\\.(Signal|AVG_Signal)$")]

# Get unique sample IDs
sample_ids <- unique(str_extract(sample_cols, "\\d+(?=\\.(Signal|AVG_Signal))"))

# Create expression matrix using Signal values (main expression values)
signal_cols <- paste0(sample_ids, ".Signal")

# Check which signal columns actually exist
existing_signal_cols <- intersect(signal_cols, colnames(expr_data))

if(length(existing_signal_cols) == 0) {
  # If no .Signal columns, try AVG_Signal
  signal_cols <- paste0(sample_ids, ".AVG_Signal")
  existing_signal_cols <- intersect(signal_cols, colnames(expr_data))
}

# Select probe ID and expression columns
expr_matrix <- expr_data %>%
  dplyr::select(PROBE_ID, all_of(existing_signal_cols))

# Rename columns to just sample IDs
colnames(expr_matrix)[-1] <- str_extract(colnames(expr_matrix)[-1], "\\d+")

## Process platform annotation data
  lines <- readLines(annotation_file)   # Read annotation file
  
  # Find the header line for probe annotation
  header_line <- which(str_detect(lines, "Block Column.*Reporter Name"))
  
  if(length(header_line) == 0) {
    # Alternative header pattern
    header_line <- which(str_detect(lines, "Reporter Name"))
  }
  
  # Read annotation data
  annotation <- read_tsv(annotation_file, skip = header_line - 1, show_col_types = FALSE)
  
  # Extract probe ID to gene symbol mapping
  probe_to_gene <- annotation %>%
    dplyr::select(
      probe_id = `Reporter Name`,
      gene_symbol = `Reporter Database Entry[hugo]`
    ) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "")

##  Merge expression data with gene symbols
  expr_with_genes <- expr_matrix %>%
    left_join(probe_to_gene, by = c("PROBE_ID" = "probe_id"))
  
  # Handle multiple probes per gene by taking the mean
  expr_matrix <- expr_with_genes %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    dplyr::select(-PROBE_ID) %>%
    group_by(gene_symbol) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = 'drop') %>%
    column_to_rownames("gene_symbol")
  
  expr_matrix <- as.matrix(expr_matrix)

## process_clinical_data
  clinical_raw <- read_tsv(sdrf_file, show_col_types = FALSE)
  
  # Extract relevant clinical information
  clinical_data <- clinical_raw %>%
    dplyr::select(
      sample_name = `Source Name`,
      assay_name = `Assay Name`,
      organism = `Characteristics[organism]`,
      specimen_type = `Characteristics[specimen with known storage state]`,
      disease = `Characteristics[disease]`,
      clinical_risk_score = `Characteristics[clinical risk score]`,
      tissue = `Characteristics[organism part]`,
      individual_id = `Characteristics[individual ID]`,
      risk_score_factor = `Factor Value[clinical risk score]`
    ) %>%
    # Use assay_name as the primary sample identifier (matches expression data)
    mutate(sample_id = assay_name) %>%
    # Remove duplicates if any
    distinct(sample_id, .keep_all = TRUE) %>%
    # Set sample_id as rownames
    column_to_rownames("sample_id")

  ## Ensure sample names match between expression and clinical data
  common_samples <- intersect(colnames(expr_matrix), rownames(clinical_data))
  
  if(length(common_samples) == 0) {
    warning("No common samples found between expression and clinical data!")
    cat("Expression samples:", head(colnames(expr_matrix)), "\n")
    cat("Clinical samples:", head(rownames(clinical_data)), "\n")
  }
  
  # Subset to common samples
  expr_matrix_final <- expr_matrix[, common_samples, drop = FALSE]
  clinical_data_final <- clinical_data[common_samples, , drop = FALSE]
  
  cat("Final data dimensions:\n")
  cat("Expression matrix:", dim(expr_matrix_final), "(genes x samples)\n")
  cat("Clinical data:", dim(clinical_data_final), "(samples x features)\n")
  cat("Number of matching samples:", length(common_samples), "\n")

# Save results
write.csv(expr_matrix_final, "exp_data.csv")
write.csv(clinical_data_final, "clinical_data.csv", row.names = TRUE)


# ==========================================
# QUALITY CONTROL CHECKS
# ==========================================

# Check for missing values
cat("\n=== QUALITY CONTROL ===\n")
cat("Missing values in expression matrix:", sum(is.na(expression_matrix)), "\n")
cat("Missing values in clinical data:", sum(is.na(clinical_data)), "\n")

# Check expression distribution
if(nrow(expression_matrix) > 0 && ncol(expression_matrix) > 0) {
  # Log2 transform if data appears to be on linear scale
  expr_sample <- expression_matrix[1:min(100, nrow(expression_matrix)), 1:min(5, ncol(expression_matrix))]
  if(max(expr_sample, na.rm = TRUE) > 100) {
    cat("Expression values appear to be on linear scale, consider log2 transformation\n")
    # Example: log2_expr <- log2(expression_matrix + 1)
  }
}

# Visualizations (optional)
# Uncomment if you want to create plots
# library(ggplot2)
# library(pheatmap)
# 
# # Sample correlation heatmap
# if(ncol(expression_matrix) > 1) {
#   sample_cor <- cor(expression_matrix, use = "pairwise.complete.obs")
#   pheatmap(sample_cor, main = "Sample Correlation")
# }
# 
# # Clinical data summary plot
# if("clinical_risk_score" %in% colnames(clinical_data)) {
#   ggplot(clinical_data, aes(x = clinical_risk_score)) +
#     geom_bar() +
#     theme_minimal() +
#     labs(title = "Distribution of Clinical Risk Scores")
# }