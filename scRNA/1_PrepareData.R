## scRNA-seq for CRLM

library(Seurat)
library(DoubletFinder)
library(harmony)
library(Matrix)

library(dplyr)
library(tidyr)

library(ggplot2)

library(RColorBrewer)
library(cowplot)
library(patchwork)
library(viridis)
library(dittoSeq)

dataPath <- "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
figurePath <- "/mnt/NAS_21T/ProjectResult/IMC-CRLM/scRNA"

source("./PrepareData_functions.R")

if (!dir.exists(figurePath)) {
    dir.create(figurePath, recursive = TRUE)
}

# =============================================================================
## 1. Load and Data pre-processing
# =============================================================================

## Che - Cell Discovery - 2021
# =============================================================================
Che_dataPath <- file.path(dataPath, "Che_2021_CellDiscovery")
cell_barcodes <- read.table(file.path(Che_dataPath, "GSE178318_barcodes.tsv.gz"))
genes <- read.table(file.path(Che_dataPath, "GSE178318_genes.tsv.gz"))
gene_matrix <- readMM(file.path(Che_dataPath, "GSE178318_matrix.mtx.gz"))

clean_gene_matrix <- resolve_duplicate_features(
    expression_matrix = gene_matrix,
    feature_names = genes$V2,
    is_count_data = TRUE
)
colnames(clean_gene_matrix) <- cell_barcodes$V1
Che_seu <- CreateSeuratObject(counts = clean_gene_matrix, project = "Che") ### Create Seurat object

### Process meta data
PaitentID <- sapply(cell_barcodes$V1, function(x) {
    temp <- strsplit(x, "_")[[1]]
    return(temp[2])
})
Tissue <- sapply(cell_barcodes$V1, function(x) {
    temp <- strsplit(x, "_")[[1]]
    return(temp[3])
})

Che_seu$patient <- PaitentID
Che_seu$tissue <- Tissue
Che_seu$patient_tissue <- paste0(PaitentID, "_", Tissue)

clinical_info <- read.csv(file.path(dataPath, "Che_2021_CellDiscovery", "clinical_info.csv"))
Che_seu@meta.data <- left_join(Che_seu@meta.data, clinical_info, by = "patient")
rownames(Che_seu@meta.data) <- cell_barcodes$V1
Che_seu$orig.ident <- Che_seu$patient_tissue

### Sample Filter
Che_seu <- Che_seu[, Che_seu$tissue == "LM"] ## Subset to Liver Metastasis (LM) samples
Che_seu <- Che_seu[, !Che_seu$Preoperative_Chemotherapy] ## Remove preoperative chemotherapy samples

### QC
Che_seu <- perform_qc(
    seurat_obj = Che_seu,
    nFeature_RNA_min = 500, nFeature_RNA_max = 6000, percent.mt = 15,
    n_hvgs = 2000, use_sct = FALSE, pca_dims = 1:50, verbose = TRUE
)
Che_seu_doublet <- perform_doublet_detect(seurat_obj = Che_seu)
Che_seu$doublet <- Che_seu_doublet@meta.data[, ncol(Che_seu_doublet@meta.data)]

### rename tissue
tissue_ <- as.character(Che_seu$tissue)
table(tissue_)
tissue_ <- ifelse(tissue_ == "LM", "TC", tissue_)
Che_seu$tissue <- tissue_

saveRDS(Che_seu, file.path(dataPath, "Che_seurat.rds"))

## Wu - Cancer Discovery - 2022
# =============================================================================

load(file.path(dataPath, "Wu_2022_CancerDiscovery", "exprmatrix.rda"))
load(file.path(dataPath, "Wu_2022_CancerDiscovery", "metadata.rda"))

Wu_seu <- CreateSeuratObject(counts = exprmatrix, meta.data = metadata, project = "Wu") ## Creat seurat object

### subset CRLM samples
Wu_seu <- Wu_seu[, Wu_seu$tissue %in% c("Liver_P", "Liver_T", "Liver_T1", "Liver_T2")]

### Rename Tissue
tissue_ <- as.character(Wu_seu$tissue)
tissue_ <- ifelse(tissue_ == "Liver_P", "PT", "TC")
Wu_seu$tissue <- tissue_

saveRDS(Wu_seu, file.path(dataPath, "Wu_seurat.rds"))

## Liu - Cancer Cell - 2022
# =============================================================================
liu_dataPath <- file.path(dataPath, "GSE164522_Liu_CancerCell_2022")

meta_data <- read.csv(file.path(liu_dataPath, "GSE164522_CRLM_metadata.csv.gz"), row.names = 1)
MN_matrix <- read.table(file.path(liu_dataPath, "GSE164522_CRLM_MN_expression.csv.gz"),row.names = 1, sep = ",", header = TRUE)
MT_matrix <- read.table(file.path(liu_dataPath, "GSE164522_CRLM_MT_expression.csv.gz"),row.names = 1, sep = ",", header = TRUE)

gene_matrix <- cbind(MN_matrix, MT_matrix)

# Check for duplicate gene names
print(table(duplicated(rownames(gene_matrix))))
clean_gene_matrix <- gene_matrix

# Check for cell barcodes
colnames(clean_gene_matrix) <- gsub("\\.", "-", colnames(clean_gene_matrix))
print(table(colnames(clean_gene_matrix)%in%rownames(meta_data)))

Liu_seu <- CreateSeuratObject(counts = clean_gene_matrix, project = "Liu") ### Create Seurat object
rm(gene_matrix, clean_gene_matrix, MN_matrix, MT_matrix)
gc()

### Process meta data
PaitentID <- meta_data$patient
Tissue <- meta_data$tissue

## Add metadata to Seurat object
cell_idx <- match(colnames(Liu_seu), rownames(meta_data))

Liu_seu$patient <- PaitentID[cell_idx]
Liu_seu$tissue <- Tissue[cell_idx]
# Liu_seu$tissue <- ifelse(Liu_seu$tissue == "metastasis normal", "PT", "TC") ## rename tissue
Liu_seu$patient_tissue <- paste0(Liu_seu$patient, "_", Liu_seu$tissue)
Liu_seu$orig.ident <- Liu_seu$patient_tissue

### QC
Liu_seu <- perform_qc(
    seurat_obj = Liu_seu,
    nFeature_RNA_min = 800, nFeature_RNA_max = 6000, percent.mt = 15,
    n_hvgs = 2000, use_sct = FALSE, pca_dims = 1:50, verbose = TRUE
)
Liu_seu_doublet <- perform_doublet_detect(seurat_obj = Liu_seu)
Liu_seu$doublet <- Liu_seu_doublet@meta.data[, ncol(Liu_seu_doublet@meta.data)]

saveRDS(Liu_seu, file.path(dataPath, "Liu_seurat.rds"))

## Wang - ScienceAdvances - 2023
# =============================================================================
Wang_dataPath <- file.path(dataPath, "GSE225857_Wang_ScienceAdvances_2023")

immune_meta <- read.table(file.path(Wang_dataPath, "GSM7058754_immune_meta.txt.gz"), row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
immune_exp <- read.table(file.path(Wang_dataPath, "GSM7058754_immune_counts.txt.gz"), row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nonimmune_meta <- read.table(file.path(Wang_dataPath, "GSM7058755_non_immune_meta.txt.gz"), row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nonimmune_exp <- read.table(file.path(Wang_dataPath, "GSM7058755_non_immune_counts.txt.gz"), row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Intersect genes
common_genes <- intersect(rownames(immune_exp), rownames(nonimmune_exp))
immune_exp <- immune_exp[common_genes, ]
nonimmune_exp <- nonimmune_exp[common_genes, ]

all_exp <- cbind(immune_exp, nonimmune_exp)

# Intersect meta data
common_colnames <- intersect(colnames(immune_meta), colnames(nonimmune_meta))
immune_meta <- immune_meta[, common_colnames]
nonimmune_meta <- nonimmune_meta[, common_colnames]

all_meta <- rbind(immune_meta, nonimmune_meta)

# Check for duplicate gene names
print(table(duplicated(rownames(all_exp))))
clean_gene_matrix <- all_exp
Wang_seu <- CreateSeuratObject(counts = clean_gene_matrix, project = "Wang") ### Create Seurat object

rm(all_exp, clean_gene_matrix, immune_exp, nonimmune_exp, immune_meta, nonimmune_meta)
gc()

### Process meta data
PaitentID <- all_meta$patients
Tissue <- all_meta$organs

# Rename tissue
Tissue <- ifelse(startsWith(Tissue, prefix = "L"), Tissue, "Other")
Tissue <- ifelse(startsWith(Tissue, prefix = "LC"), "TC", Tissue)
Tissue <- ifelse(startsWith(Tissue, prefix = "LN"), "PT", Tissue)

Wang_seu$patient <- PaitentID
Wang_seu$tissue <- Tissue
Wang_seu$patient_tissue <- paste0(PaitentID, "_", Tissue)
Wang_seu$doublet <- all_meta$doublet

Wang_seu$orig.ident <- Wang_seu$patient_tissue

### Sample Filter
Wang_seu <- Wang_seu[, Wang_seu$tissue != "Other"] ## Subset to Liver Metastasis (LM) samples

### QC
Wang_seu <- perform_qc(
    seurat_obj = Wang_seu,
    nFeature_RNA_min = 200, nFeature_RNA_max = 6000, percent.mt = 40,
    n_hvgs = 2500, use_sct = FALSE, pca_dims = 1:20, verbose = TRUE
)

saveRDS(Wang_seu, file.path(dataPath, "Wang_seu.rds"))

## FDZH - 6 Samples from invasive margin (IM)
# =============================================================================

### Load data
FDZS_dataPath <- file.path(dataPath, "FDZS")

# Get all sample directories
sample_dirs <- list.dirs(FDZS_dataPath, recursive = FALSE, full.names = FALSE)

cat("Found samples:", paste(sample_dirs, collapse = ", "), "\n")

# Initialize list to store individual Seurat objects
FDZS_seurat_list <- list()

### Load each sample individually
for (i in seq_along(sample_dirs)) {
    sample_id <- sample_dirs[i]
    sample_path <- file.path(FDZS_dataPath, sample_id)

    cat("Processing sample:", sample_id, "\n")

    # Load data
    cell_barcodes <- read.table(file.path(sample_path, "barcodes.tsv.gz"))
    genes <- read.table(file.path(sample_path, "features.tsv.gz"))
    gene_matrix <- readMM(file.path(sample_path, "matrix.mtx.gz"))

    clean_gene_matrix <- resolve_duplicate_features(
        expression_matrix = gene_matrix,
        feature_names = genes$V1,
        is_count_data = TRUE
    )
    colnames(clean_gene_matrix) <- cell_barcodes$V1

    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = clean_gene_matrix, project = "FDZS") ### Create Seurat object

    # Add sample metadata
    seurat_obj$patient <- sample_id
    seurat_obj$tissue <- "IM"
    seurat_obj$orig.ident <- paste0("FDZS_", sample_id)

    # Apply QC to each sample separately
    seurat_obj <- perform_qc(
        seurat_obj = seurat_obj,
        nFeature_RNA_min = 500,
        nFeature_RNA_max = 6000,
        percent.mt = 15,
        n_hvgs = 2000,
        use_sct = FALSE, # Start with traditional workflow
        verbose = TRUE # Reduce verbosity for multiple samples
    )

    # Apply doublet detection
    processed_obj <- perform_doublet_detect(seurat_obj = processed_obj)
    seurat_obj$doublet <- processed_obj@meta.data[, ncol(processed_obj@meta.data)]

    # Store in list
    FDZS_seurat_list[[sample_id]] <- seurat_obj

    cat("  Sample", sample_id, ":", nrow(seurat_obj), "features x", ncol(seurat_obj), "cells\n")
}

### Merge data and Integrated analysis
FDZS_merged <- merge(
    x = FDZS_seurat_list[[1]],
    y = FDZS_seurat_list[2:length(FDZS_seurat_list)],
    add.cell.ids = names(FDZS_seurat_list),
    project = "FDZS"
)

saveRDS(FDZS_merged, file.path(dataPath, "FDZS_seurat.rds"))

### Simply Visualization
if (F) {
    # Re-process merged data
    FDZS_merged <- NormalizeData(FDZS_merged, verbose = FALSE)
    FDZS_merged <- FindVariableFeatures(FDZS_merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    FDZS_merged <- ScaleData(FDZS_merged, verbose = FALSE)
    FDZS_merged <- RunPCA(FDZS_merged, npcs = 50, verbose = FALSE)

    FDZS_merged <- RunUMAP(FDZS_merged, reduction = "pca", dims = 1:30, verbose = FALSE)
    FDZS_merged <- FindNeighbors(FDZS_merged, reduction = "pca", dims = 1:30, verbose = FALSE)
    FDZS_merged <- FindClusters(FDZS_merged, resolution = 0.5, verbose = FALSE)

    ### Visualization
    cat("\n=== Creating visualizations ===\n")

    # Sample distribution plot
    p1 <- DimPlot(FDZS_merged, group.by = "patient", pt.size = 0.5) +
        ggtitle("FDZS Samples - UMAP") +
        theme_bw()

    # Clustering plot
    p2 <- DimPlot(FDZS_merged, label = TRUE, pt.size = 0.5) +
        ggtitle("FDZS Clusters - UMAP") +
        theme_bw()

    # Combined plot
    p_combined <- plot_grid(p1, p2, ncol = 2)

    # Save plots
    pdf(file.path(figurePath, "FDZS_UMAP_samples_and_clusters.pdf"), width = 16, height = 6)
    print(p_combined)
    dev.off()
}

# =============================================================================
## 2. Data integrate
# =============================================================================
# Step 0: Load objects first
Che_seu <- readRDS(file.path(dataPath, "Che_seurat.rds"))
Wu_seu <- readRDS(file.path(dataPath, "Wu_seurat.rds"))
FDZS_seu <- readRDS(file.path(dataPath, "FDZS_seurat.rds"))

# Step 1: Merge the objects
list_to_merge <- list(
    "Che" = Che_seu,
    "Wu" = Wu_seu,
    "FDZS" = FDZS_seu
)
merged <- merge(
    x = list_to_merge[[1]],
    y = list_to_merge[2:length(list_to_merge)],
    add.cell.ids = names(list_to_merge),
    project = "CRLM"
)
table(merged$orig.ident)
merged <- merged[, merged$doublet == "Singlet"]

rm(list_to_merge, Che_seu, Wu_seu, FDZS_seu)
gc()

# Step 2: Standard preprocessing
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 50)

# Step 3: Run Harmony integration
merged <- RunHarmony(merged, "orig.ident")

# Step 4: Use harmony reduction for downstream analysis
merged <- RunUMAP(merged, reduction = "pca", dims = 1:15)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:15)
merged <- FindClusters(merged, resolution = 0.5)
p1 <- DimPlot(merged, reduction = "umap", group.by = "orig.ident", pt.size = .5)
p2 <- DimPlot(merged, reduction = "umap", label = TRUE, pt.size = .5)
p_pca <- plot_grid(p1, p2)

merged <- RunUMAP(merged, reduction = "harmony", dims = 1:15)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:15)
merged <- FindClusters(merged, resolution = 0.5)
p1 <- DimPlot(merged, reduction = "umap", group.by = "orig.ident", pt.size = .5)
p2 <- DimPlot(merged, reduction = "umap", label = TRUE, pt.size = .5)
p_harmony <- plot_grid(p1, p2)

pdf(file.path(figurePath, paste0("UMAP of pre-harmony.pdf")), width = 18, height = 6)
print(p_pca)
dev.off()

pdf(file.path(figurePath, paste0("UMAP of post-harmony.pdf")), width = 18, height = 6)
print(p_harmony)
dev.off()

saveRDS(merged, file.path(dataPath, "merge_seurat.rds"))

## 3. Merge Data analysis
# =============================================================================
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:15)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:15)
merged <- FindClusters(merged, resolution = 0.5)

### Major annoatation
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.1)

features <- c("CD3D", "KLRF1", "CD79A", "COL1A1", "CLDN5", "LILRA4", "LYZ", "FCGR3B", "EPCAM")

### Features plot of markers
set.seed(619)
idx <- sample(1:ncol(merged), size = 10000)
merged_subset <- merged[, idx] ## Downsample

p1 <- DimPlot(merged_subset,
    reduction = "umap",
    label = TRUE,
    label.size = 4,
    label.color = "white",
    label.box = TRUE,
    pt.size = 0.8,
    raster = FALSE,
    shuffle = TRUE
) +
    theme_bw() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "UMAP Clustering") +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

pdf(file.path(figurePath, paste0("UMAP of seurat cluster.pdf")), width = 8, height = 6)
print(p1)
dev.off()

p <- FeaturePlot(merged_subset,
    features = features,
    pt.size = 0.5,
    raster = FALSE,
    ncol = 3,
    cols = c("lightgrey", "red")
) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold")
    )

pdf(file.path(figurePath, paste0("UMAP of major features.pdf")), width = 12, height = 10)
print(p)
dev.off()

### Assign major clusters
new.cluster.ids <- c(
    "T", "T", "T", "T", "NK",
    "Myeloid", "Myeloid", "Neutrophils", "Epithelial", "B",
    "B", "T", "Epithelial", "Stromal", "pDC"
)
names(new.cluster.ids) <- levels(merged_subset)
merged_subset <- RenameIdents(merged_subset, new.cluster.ids)
merged <- RenameIdents(merged, new.cluster.ids)
merged$major_type <- as.character(Idents(merged))

## Self defined color
celltypes <- unique(new.cluster.ids)
majortype_colors <- setNames(ggsci::pal_nejm("default")(length(celltypes)), celltypes)

p1 <- DimPlot(merged_subset,
    reduction = "umap", label = FALSE, pt.size = 0.5,
    cols = merged@misc$colors$majorcell_types
) +
    theme_bw() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "UMAP Clustering") +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

pdf(file.path(figurePath, paste0("UMAP of major celltype.pdf")), width = 8, height = 6)
print(p1)
dev.off()

### Save
saveRDS(merged, file.path(dataPath, "merge_seurat.rds"))


## 4. Only peritumor region Data analysis
# =============================================================================
Wu_seu <- readRDS(file.path(dataPath, "Wu_seurat.rds"))

Wu_seu_LM_P <- subset(Wu_seu, tissue == "Liver_P")
Wu_seu_LM_P <- perform_qc(Wu_seu_LM_P)
Idents(Wu_seu_LM_P) <- Wu_seu_LM_P$sub_cell_type

## Self defined color
celltypes <- unique(Wu_seu_LM_P$sub_cell_type)
colors <- setNames(dittoColors()[1:length(celltypes)], celltypes)

set.seed(619)
idx <- sample(1:ncol(Wu_seu_LM_P), size = 10000)
Wu_seu_LM_P_subset <- Wu_seu_LM_P[, idx] ## Downsample

p1 <- DimPlot(Wu_seu_LM_P_subset,
    reduction = "umap",
    cols = colors,
    label = FALSE,
    label.size = 4,
    label.color = "white",
    label.box = TRUE,
    pt.size = 0.5,
    raster = FALSE,
    shuffle = TRUE
) +
    theme_bw() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "UMAP Clustering") +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

pdf(file.path(figurePath, paste0("UMAP of Wu data peritumor annotation.pdf")), width = 14, height = 8)
print(p1)
dev.off()
