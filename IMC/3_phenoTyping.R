## clustering and Annotation

library(imcRtools)
library(cytomapper)
library(SpatialExperiment)
library(scater)
library(CATALYST)
library(dplyr)
library(batchelor)
library(bluster) 
library(BiocParallel)
library(scran)

## plot
library(RColorBrewer)
library(dittoSeq)
library(viridis)
library(patchwork)
library(cowplot)
library(scuttle)
library(ggplot2)
library(ggrepel)

set.seed(619)

# ------------------------------
# Set Working Directory and Source Utilities
# ------------------------------
workDir <- "/mnt/public/lyx/IMC_HE_Merge/CRLM"

# Define directories for saving results and Steinbock outputs
setwd(workDir)
source("./code/utils.R")
source("./code/visualize.R")

saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","1_Anno")

# Create the results directory if it doesn't exist
if(!dir.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

# Load data
spe <- readRDS(file.path(saveDir,"spe_correct.rds"))

# Load marker table for annotation
# write.csv(as.data.frame(rowData(spe)),file.path(saveDir,"panel.csv"))
anno_df <- read.csv(file.path(saveDir,"panel.csv"),header = T,row.names = 1)
rowData(spe) <- anno_df

## Cell phenotyping
### Step 1: Major clustering
if(T){
  
  ## Shared nearest neighbor
  if (F){
    # Select the corrected cell embeddings for clustering 
    mat <- reducedDim(spe, "fastMNN")
    
    # Perform the cluster sweep
    combinations <- clusterSweep(mat, BLUSPARAM=SNNGraphParam(),
                                 k=c(10L, 20L, 30L),
                                 type = c("rank", "jaccard"),
                                 cluster.fun = "louvain",
                                 BPPARAM = SerialParam(RNGseed = 230214))
    
    # Compute the average silhouette width per parameter combination
    sil <- vapply(as.list(combinations$clusters), function(x) mean(approxSilhouette(mat, x)$width), 0)
    
    # Visualize the average silhouette width per parameter
    # combination
    ggplot(data.frame(method = names(sil), sil = sil)) +
      geom_point(aes(method, sil), size = 3) +
      theme_classic(base_size = 15) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Cluster parameter combination") + ylab("Average silhouette width")
    
    clusters <- clusterCells(spe[rowData(spe)$Major,],
                             use.dimred = "fastMNN",
                             BLUSPARAM = SNNGraphParam(k = 20, cluster.fun = "louvain", type = "jaccard"))
    spe$nn_clusters <- clusters
  }
    
  ## Rphenograph
  if(F){
    library(Rphenograph)
    library(igraph)
    
    mat <- reducedDim(spe[rowData(spe)$Major,], "fastMNN")
    out <- Rphenograph(mat, k = 25) 
  }
    
  ## FlowSOM + ConsensusClusterPlus
  if(T){
    library(kohonen)
    library(ConsensusClusterPlus)
    
    # Select integrated cells
    mat <- reducedDim(spe[rowData(spe)$Major,], "fastMNN")
    
    # Perform SOM clustering
    set.seed(220410)
    som.out <- clusterRows(mat, SomParam(100), full = TRUE)
    
    # Cluster the 100 SOM codes into larger clusters
    maxK_ <- 20
    ccp <- ConsensusClusterPlus(t(som.out$objects$som$codes[[1]]),
                                maxK = maxK_,
                                reps = 100, 
                                distance = "euclidean", 
                                seed = 220410, 
                                plot = NULL)
    
    clusters <- ccp[[maxK_]][["consensusClass"]][som.out$clusters]
    spe$nn_clusters <- as.factor(clusters)
    
    cur_cells <- sample(seq_len(ncol(spe)), 10000)
    spe_subset <- spe[,cur_cells]
    spe_subset <- runUMAP(spe_subset, dimred= "fastMNN", name = "UMAP_mnnCorrected") 

    p1 <- dittoDimPlot(spe_subset, 
                 var = "nn_clusters", 
                 reduction.use = "UMAP_mnnCorrected", 
                 size = 0.2,do.label = TRUE) +
      ggtitle("SOM clusters on UMAP, integrated cells")
    
    pdf(file.path(figureDir,paste0("UMAP of Major cell type clustering via SOM.pdf")),width = 8,height = 6)
    print(p1)
    dev.off()
  }

  p <- plot_cluster_Heatmap(spe_ = spe, marker_use_col = "Major", clustercol = "nn_clusters", cutoff = FALSE)
  
  pdf(file.path(figureDir,paste0("Heatmap of Major cell type clustering via SOMt.pdf")),width = 20,height = 7.5)
  print(p)
  dev.off()
  
  celltype_mean_temp <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),
                                             ids = spe$nn_clusters, 
                                             statistics = "mean",
                                             use.assay.type = "exprs")
  frac_table <- as.data.frame(assay(celltype_mean_temp))
  write.csv(frac_table,file = file.path(figureDir, "major cluster mean value.csv") ,row.names = TRUE)
}

### Step 2: Mannual Annotate Major types
if(T){ ## Stromal, Myeloid, Lymphocyte, Epithelial, Unlabelled
    major_celltype <- recode(as.character(spe$nn_clusters),
                             "1" = "Stromal","2" = "Stromal","3" = "Lymphocyte","4" = "Lymphocyte","5" = "Myeloid",
                             "6" = "Epithelial","7" = "Epithelial","8" = "Epithelial","9" = "Myeloid","10" = "Stromal",
                             "11" = "Myeloid","12" = "Epithelial","13" = "Epithelial","14" = "Epithelial","15" = "Epithelial",
                             "16" = "Myeloid","17" = "Myeloid","18" = "Myeloid","19" = "Lymphocyte","20" = "Myeloid"
                             )

  cat("The major annotation results is:","\n")
  table(major_celltype)
  
  spe$Major_type <- major_celltype  
  saveRDS(spe,file.path(figureDir,"Major_SOM_cluster_0714.rds"))
  
  ## view the expression of major types
  if(T){
    celltype_mean_temp <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),
                                          ids = spe$Major_type, 
                                          statistics = "mean",
                                          use.assay.type = "exprs", 
                                          subset.row = rownames(spe)[rowData(spe)$Major])
    
    ## Expression Heatmap
    p <- dittoHeatmap(celltype_mean_temp,
                      assay = "exprs", 
                      cluster_cols = TRUE, show_colnames = TRUE,
                      annot.by = c("Major_type", "ncells"),
                      scale = "column",
                      annotation_colors = list(ncells = plasma(100)))
    pdf(file.path(figureDir,"Manual annotation of major celltypes.pdf"),width = 6,height = 7.5)
    print(p)
    dev.off()
    
    ## UMAP
    set.seed(619)
    cur_cells <- sample(seq_len(ncol(spe)), 10000)
    spe_subset <- spe[,cur_cells]
    spe_subset <- runUMAP(spe_subset, dimred= "fastMNN", name = "UMAP_mnnCorrected") 

    p1 <- dittoDimPlot(spe_subset, var = "Major_type", reduction.use = "UMAP_mnnCorrected", size = 0.4,do.label = TRUE) +
      ggtitle("SOM clusters on UMAP, integrated cells")
    
    pdf(file.path(figureDir,"UMAP of Manual annotation of major celltypes.pdf"),width = 8,height = 6)
    print(p1)
    dev.off()
    
    rm(celltype_mean_temp,p)
    gc()
  }
  
}

### Step 3: Clustering for subtype
majortypes <- unique(spe$Major_type)
majortypes <- majortypes[-match("Unlabelled",majortypes)] ## omit unlabel type

# collectDict <- list(
#   "Immune" = c("T","B"),
#   "Myeloid" = c("Myeloid"),
#   "Stromal" = c("Endothelial","SMC","Fibroblast"),
#   "Neuron" = c("Neuron","Astrocyte","Oligodendrocytes")
# )
collectDict <- list(
  "Lymphoid" = c("Lymphocyte"),
  "Myeloid" = c("Myeloid"),
  "Stromal" = c("Stromal"),
  "EC" = c("Epithelial")
)

if(T){
  maxK_ <- 20
  
  for(selective_types_ in names(collectDict)) {
    set.seed(619)
    
    # Subset data for the current major type
    types_ <- collectDict[[selective_types_]]
    spe_subset <- spe[, spe$Major_type %in% types_]
    
    # Perform unsupervised clustering on the subset
    mat <- reducedDim(spe_subset[rowData(spe_subset)$selective_types_,], "fastMNN")
    
    # Cluster the 100 SOM codes into larger clusters
    som.out <- clusterRows(mat, SomParam(100), full = TRUE)
    ccp <- ConsensusClusterPlus(t(som.out$objects$som$codes[[1]]),maxK = maxK_,reps = 100, distance = "euclidean",plot = NULL)
    
    clusters <- ccp[[maxK_]][["consensusClass"]][som.out$clusters]
    spe_subset$minor_som_clusters <- as.factor(clusters)
    
    # plot the heatmap of clustering results
    p <- plot_cluster_Heatmap(spe_ = spe_subset, marker_use_col = selective_types_, clustercol = "minor_som_clusters", cutoff = FALSE)
    pdf(file.path(figureDir,paste0("Heatmap of ",selective_types_," cell type clustering via SOM(k=",maxK_,").pdf")),width = 15,height = 6)
    print(p)
    dev.off()
    
    # Store the clustering results
    saveRDS(spe_subset,file.path(figureDir,paste0(selective_types_,"_SOM_cluster_0714.rds")))
  }
}

### Step 4: Annotation Subtypes
if(T){
  collectDict <- list(
    "Lymphoid" = c("Lymphocyte"),
    "Myeloid" = c("Myeloid"),
    "Stromal" = c("Stromal"),
    "EC" = c("Epithelial")
  )
  spe$sub_celltype <- "unlabelled"

  ## Lymphocyte
  if(T){ ## B, CD8T, CD4T, Treg, NK, unlabelled
    subtype_ <- "Lymphoid"
    spe_temp <- readRDS(file.path(figureDir,paste0(subtype_,"_SOM_cluster_0714.rds")))
    anno_cellsubtype <- recode(as.character(spe_temp$minor_som_clusters),
                               "1" = "unlabelled","2" = "B","3" = "B","4" = "unlabelled","5" = "B",
                               "6" = "unlabelled","7" = "B","8" = "B","9" = "unlabelled","10" = "unlabelled",
                               "11" = "unlabelled","12" = "unlabelled","13" = "CD8T","14" = "unlabelled","15" = "CD4T",
                               "16" = "unlabelled","17" = "Treg","18" = "unlabelled","19" = "NK","20" = "NK"
    )
    # Assign "unlabelled" to clusters that are not assigned  
    anno_cellsubtype[anno_cellsubtype%in%""] <- "unlabelled" 
    
    # Merge to all 
    spe_temp$sub_celltype <- anno_cellsubtype
    print(table(spe_temp$sub_celltype))
    
    idx <- intersect(rownames(colData(spe_temp)),rownames(colData(spe)))
    spe$sub_celltype[match(idx,rownames(colData(spe)))] <- anno_cellsubtype[match(idx,rownames(colData(spe_temp)))]
    
    print(table(spe$sub_celltype))
  }
  
  ## Myeloid
  if(T){ ## Mono_CD14, Macro_CD11b, Macro_CD163, Macro_HLADR, Macro_Other, DC, unlabelled
    subtype_ <- "Myeloid"
    spe_temp <- readRDS(file.path(figureDir,paste0(subtype_,"_SOM_cluster_0714.rds")))
    anno_cellsubtype <- recode(as.character(spe_temp$minor_som_clusters),
                               "1" = "unlabelled","2" = "Macro_HLADR","3" = "unlabelled","4" = "Macro_HLADR","5" = "unlabelled",
                               "6" = "Macro_HLADR","7" = "unlabelled","8" = "Macro_HLADR","9" = "unlabelled","10" = "Mono_CD14",
                               "11" = "Macro_Other","12" = "Macro_CD11b","13" = "Macro_CD11b","14" = "Macro_HLADR","15" = "Macro_CD11b",
                               "16" = "Mono_CD14","17" = "Macro_CD11b","18" = "Macro_CD11b","19" = "Mono_CD14","20" = "Mono_CD14"
    )
    # Assign "unlabelled" to clusters that are not assigned  
    anno_cellsubtype[anno_cellsubtype%in%""] <- "unlabelled" 
    
    # Merge to all 
    spe_temp$sub_celltype <- anno_cellsubtype

    idx <- intersect(rownames(colData(spe_temp)),rownames(colData(spe)))
    spe$sub_celltype[match(idx,rownames(colData(spe)))] <- anno_cellsubtype[match(idx,rownames(colData(spe_temp)))]
    print(table(spe_temp$sub_celltype))
    
    print(table(spe$sub_celltype))

  }
  
  ## Stromal
  if(T){ ## CAF (all high), Fibroblast (Collagen-I+ Vimentin+)
    subtype_ <- "Stromal"
    spe_temp <- readRDS(file.path(figureDir,paste0(subtype_,"_SOM_cluster_0714.rds")))
    anno_cellsubtype <- recode(as.character(spe_temp$minor_som_clusters),
                               "1" = "SC_aSMA_Vimentin","2" = "SC_aSMA_Vimentin","3" = "SC_Vimentin","4" = "SC_Vimentin","5" = "CAF",
                               "6" = "SC_aSMA","7" = "SC_aSMA_Collagen","8" = "SC_aSMA_Collagen","9" = "CAF","10" = "CAF",
                               "11" = "SC_Collagen_Vimentin","12" = "SC_aSMA_Collagen","13" = "SC_Collagen_Vimentin","14" = "SC_aSMA_Collagen","15" = "SC_aSMA_Collagen",
                               "16" = "SC_Collagen","17" = "SC_Collagen","18" = "SC_Collagen","19" = "SC_Collagen","20" = "unlabelled"
    )
    # Assign "unlabelled" to clusters that are not assigned  
    anno_cellsubtype[anno_cellsubtype%in%""] <- "unlabelled" 
    
    # Merge to all 
    spe_temp$sub_celltype <- anno_cellsubtype
    print(table(spe_temp$sub_celltype))
    
    idx <- intersect(rownames(colData(spe_temp)),rownames(colData(spe)))
    spe$sub_celltype[match(idx,rownames(colData(spe)))] <- anno_cellsubtype[match(idx,rownames(colData(spe_temp)))]
    
    print(table(spe$sub_celltype))

  }
  
  ## EC
  if(T){ ## EC_CAIX, EC_Vimentin, EC_Ki67, EC_EpCAM, EC_GLUT1
    subtype_ <- "EC"
    spe_temp <- readRDS(file.path(figureDir,paste0(subtype_,"_SOM_cluster_0714.rds")))
    
    anno_cellsubtype <- recode(as.character(spe_temp$minor_som_clusters),
                               "1" = "EC_CAIX","2" = "EC_CAIX","3" = "EC_CAIX","4" = "EC_Ki67","5" = "EC_EpCAM",
                               "6" = "EC_Vimentin","7" = "EC_CAIX","8" = "EC_EpCAM","9" = "EC_EpCAM","10" = "EC_EpCAM",
                               "11" = "EC_EpCAM","12" = "EC_Vimentin","13" = "EC_GLUT1","14" = "EC_Vimentin","15" = "EC_GLUT1",
                               "16" = "EC_Vimentin","17" = "EC_GLUT1","18" = "EC_Ki67","19" = "unlabelled","20" = "unlabelled"
    )
    # Assign "unlabelled" to clusters that are not assigned  
    anno_cellsubtype[anno_cellsubtype%in%""] <- "unlabelled" 
    
    # Merge to all 
    spe_temp$sub_celltype <- anno_cellsubtype
    print(table(spe_temp$sub_celltype))
    
    idx <- intersect(rownames(colData(spe_temp)),rownames(colData(spe)))
    spe$sub_celltype[match(idx,rownames(colData(spe)))] <- anno_cellsubtype[match(idx,rownames(colData(spe_temp)))]
    
    print(table(spe$sub_celltype))
  } 
  
  ## Unlabeled
  if(T){
    unlabelled_idx <- (spe$Major_type == "unlabelled" | spe$sub_celltype == "unlabelled")
    spe$sub_celltype[unlabelled_idx] <- "unlabelled"
    
    # Subset unlabel
    set.seed(619)
    spe_subset <- spe[, spe$sub_celltype == "unlabelled"]
    
    # Perform unsupervised clustering on the subset
    mat <- reducedDim(spe_subset[rowData(spe_subset)$use_channel,], "fastMNN")
    
    # Cluster the 500 SOM codes into larger clusters
    som.out <- clusterRows(mat, SomParam(500), full = TRUE)
    maxK_ = 50
    ccp <- ConsensusClusterPlus(t(som.out$objects$som$codes[[1]]),maxK = maxK_,reps = 500, distance = "euclidean",plot = NULL)
    
    clusters <- ccp[[maxK_]][["consensusClass"]][som.out$clusters]
    spe_subset$minor_som_clusters <- as.factor(clusters)
    
    # plot the heatmap of clustering results
    p <- plot_cluster_Heatmap(spe_ = spe_subset, marker_use_col = "use_channel", clustercol = "minor_som_clusters", cutoff = FALSE)
    pdf(file.path(figureDir,paste0("Heatmap of unlabelled cell type clustering via SOM(k=",maxK_,").pdf")),width = 15,height = 6)
    print(p)
    dev.off()
    
    # Store the clustering results
    saveRDS(spe_subset,file.path(figureDir,paste0("unlabelled","_SOM_cluster_0714.rds")))
    
    # Annotation
    subtype_ <- "unlabelled"
    spe_temp <- readRDS(file.path(figureDir,paste0(subtype_,"_SOM_cluster_0714.rds")))
    anno_cellsubtype <- recode(as.character(spe_temp$minor_som_clusters),
                               "1" = "Mono_CD14","2" = "Macro_Other","3" = "Macro_CD163","4" = "Macro_CD163","5" = "Macro_Other",
                               "6" = "DC","7" = "Macro_Other","8" = "Mono_CD14","9" = "Macro_CD163","10" = "Mono_CD14",
                               "11" = "Other_Immune","12" = "unlabelled","13" = "Mono_CD14","14" = "Macro_HLADR","15" = "Macro_CD163",
                               "16" = "Mono_CD16","17" = "SC_Collagen","18" = "SC_Collagen","19" = "Macro_Other","20" = "Other_Immune",
                               "21" = "SC_Vimentin","22" = "Macro_HLADR","23" = "B","24" = "SC_Vimentin","25" = "NK",
                               "26" = "Other_Immune","27" = "SC_Collagen_Vimentin","28" = "EC_Vimentin","29" = "SC_Vimentin","30" = "B",
                               "31" = "Other_Immune","32" = "DC","33" = "Other_Immune","34" = "B","35" = "CD4T",
                               "36" = "Macro_CD11b","37" = "SC_aSMA_Vimentin","38" = "Other_Immune","39" = "Other_Immune","40" = "SC_Collagen",
                               "41" = "CD4T","42" = "Macro_HLADR","43" = "SC_Collagen_Vimentin","44" = "B","45" = "CD8T",
                               "46" = "B","47" = "Other_Immune","48" = "B","49" = "B","50" = "Other_Immune"
    )
    # Assign "unlabelled" to clusters that are not assigned  
    anno_cellsubtype[anno_cellsubtype%in%""] <- "unlabelled" 
    
    # Merge to all 
    spe_temp$sub_celltype <- anno_cellsubtype
    print(table(spe_temp$sub_celltype))
    
    idx <- intersect(rownames(colData(spe_temp)),rownames(colData(spe)))
    spe$sub_celltype[match(idx,rownames(colData(spe)))] <- anno_cellsubtype[match(idx,rownames(colData(spe_temp)))]
    
    print(table(spe$sub_celltype))
  }
  
  ### Set color for spe
  celltype_color <- setNames(dittoColors()[1:length(unique(spe$sub_celltype))],unique(spe$sub_celltype))
  celltype_color["unlabelled"] <- "grey75"
  metadata(spe)$color_vectors["sub_celltype"] <- list(celltype_color)
  
  
  ### Mannual merge
  if(F){
    ## Merge subtype
    sub_celltype <- spe$sub_celltype
    sub_celltype[sub_celltype == "Monocyte"] <- "Mono_CD14"
    sub_celltype[sub_celltype == "unlabelled"] <- "SC_Collagen"
    sub_celltype[sub_celltype == "SC_aSMA"] <- "CAF"
    
    ## Set Major type
    subtype_to_major <- c(
      "CD4T" = "T_cells",
      "CD8T" = "T_cells", 
      "Treg" = "T_cells",
      "B" = "B_cells",
      "NK" = "NK_cells",
      "Mono_CD14" = "Monocytes",
      "Mono_CD16" = "Monocytes",
      "Macro_CD11b" = "Macrophages",
      "Macro_CD163" = "Macrophages",
      "Macro_HLADR" = "Macrophages",
      "Macro_Other" = "Macrophages",
      "DC" = "DCs",
      "EC_EpCAM" = "Epithelial",
      "EC_GLUT1" = "Epithelial",
      "EC_Ki67" = "Epithelial",
      "EC_Ki67_CAIX" = "Epithelial",
      "EC_CAIX" = "Epithelial",
      "EC_Vimentin" = "Epithelial",
      "CAF" = "Stromal",
      "SC_aSMA_Collagen" = "Stromal",
      "SC_aSMA_Vimentin" = "Stromal",
      "SC_Collagen" = "Stromal",
      "SC_Collagen_Vimentin" = "Stromal",
      "SC_Vimentin" = "Stromal",
      "Other_Immune" = "Other"
    )
    
    major_celltype <- subtype_to_major[sub_celltype]
    
    ## Set color
    major_colors <- c(
      # Immune cells (Blue family)
      "T_cells" = "#1f77b4",      # Base blue
      "B_cells" = "#3f8fbf",      # Slightly different blue
      "NK_cells" = "#5fa7ca",     # Another blue variant
      
      # Myeloid cells (Purple family)  
      "Monocytes" = "#9467bd",    # Base purple
      "Macrophages" = "#8c5aa6", # Slightly different purple
      "DCs" = "#a373c4",         # Another purple variant
      
      # Single types (distinct colors)
      "Epithelial" = "#2ca02c",   # Green
      "Stromal" = "#ff7f0e",      # Orange
      "Other" = "#7f7f7f"         # Gray
    )
    
    subtype_colors <- c(
      # T_cells, B_cells, NK_cells (Blues - Immune)
      "CD4T" = "#1f77b4",        # Base blue
      "CD8T" = "#aec7e8",        # Light blue
      "Treg" = "#0d47a1",        # Dark blue
      "B" = "#42a5f5",           # Medium blue
      "NK" = "#90caf9",          # Very light blue
      
      # Monocytes, Macrophages, DCs (Purples - Myeloid)
      "Mono_CD14" = "#9467bd",   # Base purple
      "Mono_CD16" = "#ce93d8",   # Light purple
      "Macro_CD11b" = "#673ab7", # Dark purple
      "Macro_CD163" = "#ba68c8", # Medium purple
      "Macro_HLADR" = "#4a148c", # Very dark purple
      "Macro_Other" = "#ab47bc", # Medium light purple
      "DC" = "#e1bee7",          # Very light purple
      
      # Epithelial (Greens)
      "EC_EpCAM" = "#2ca02c",    # Base green
      "EC_GLUT1" = "#66bb6a",    # Light green
      "EC_Ki67" = "#1b5e20",     # Dark green
      "EC_Ki67_CAIX" = "#2e7d32", # Dark medium green
      "EC_CAIX" = "#388e3c",     # Medium green
      "EC_Vimentin" = "#81c784", # Light medium green
      
      # Stromal (Oranges)
      "CAF" = "#ff7f0e",         # Base orange
      "SC_aSMA_Collagen" = "#e65100",    # Dark orange
      "SC_aSMA_Vimentin" = "#ff8f00",    # Medium dark orange
      "SC_Collagen" = "#ffb74d",         # Light orange
      "SC_Collagen_Vimentin" = "#ffcc80", # Very light orange
      "SC_Vimentin" = "#ffe0b2",         # Pale orange
      
      # Other (Gray)
      "Other_Immune" = "#7f7f7f"  # Gray
    )
    
    spe$sub_celltype <- sub_celltype
    spe$major_celltype <- major_celltype
    
    metadata(spe)$color_vectors$major_celltype <- major_colors
    metadata(spe)$color_vectors$sub_celltype <- subtype_colors
  }
  
  ### Save
  saveRDS(spe,file.path(saveDir,"subanno_spe_0714.rds"))
}

## calculate and save the fraction of clusters within each image
spe <- readRDS(file.path(saveDir,"subanno_spe_0714.rds"))

### Plot heatmap
### aggregate by cell type
celltype_order <- c(
  # B cells  
  "B",
  
  # T cells
  "CD4T", "CD8T", "Treg", "NK","Other_Immune", 
  
  # Myeloid cells
  "Mono_CD14", "Mono_CD16", 
  "Macro_CD163", "Macro_CD11b", "Macro_HLADR", "Macro_Other","DC",
  

  # Stromal cells
  "CAF", "SC_aSMA_Collagen", "SC_aSMA_Vimentin", 
  "SC_Collagen", "SC_Collagen_Vimentin", "SC_Vimentin",
  
  # Epithelial cells (tumor)
  "EC_EpCAM", "EC_GLUT1", "EC_Ki67", "EC_Ki67_CAIX", "EC_CAIX", "EC_Vimentin"
)
displayed_genes <- c(
  # B cell markers  
  "CD45", "CD20",
  
  # T cell markers
  "CD3", "CD4", "CD8a", "FoxP3", "CD279", "Ki67", "CD57",
  
  # Myeloid markers
  "CD14", "CD16", "CD68", "CD163", "CD11b", "HLA_DR", "CD11c", "CLEC9A",
  
  # Stromal markers
  "FAP", "Alpha_SMA", "Collagen_I", "VEGF",
  
  # Epithelial/tumor markers
  "EpCAM", "GLUT1", "CA_IX", "Vimentin"
)

celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),
                                      ids = spe$sub_celltype, 
                                      statistics = "mean",
                                      use.assay.type = "exprs", 
                                      #subset.row = rowData(spe)$use_channel)
                                      subset.row = displayed_genes)

### Expression Heatmap (all subpopulations)
p <- dittoHeatmap(celltype_mean,
                  assay = "exprs", genes = displayed_genes,
                  cluster_cols = FALSE, cluster_rows = FALSE, 
                  show_colnames = TRUE, show_rownames = TRUE,
                  order.by = match(celltype_order,colnames(celltype_mean)),
                  annot.by = c("sub_celltype", "ncells"),
                  scale = "column",
                  annotation_colors = metadata(spe)$color_vectors["sub_celltype"],
                  # heatmap.colors = custom_palette,
                  # breaks = breaks
)
pdf(file.path(figureDir,"All_subpopulations_markers_heatmap_ditto(scale celltype)(0714).pdf"),width = 10,height = 7.5)
print(p)
dev.off()

### Plot the umap of each subpopulations
if(T){
  ## All cell subpopulations
  set.seed(619)
  sample_size <- 1e5
  spe_subset <- spe[,sample.int(n=ncol(spe),size = sample_size, replace = F)]
  spe_subset <- spe_subset[,!spe_subset$sub_celltype %in% c("unlabelled")]
  
  spe_subset <- runUMAP(spe_subset, subset_row = displayed_genes, exprs_values = "exprs")

  ## UMAP for major types
  set.seed(619)
  cur_cells <- sample(seq_len(ncol(spe)), 20000)
  spe_subset <- spe[,cur_cells]
  spe_subset <- runUMAP(spe_subset, dimred= "fastMNN", name = "UMAP_mnnCorrected") 
  
  p1 <- dittoDimPlot(spe_subset, var = "major_celltype", reduction.use = "UMAP_mnnCorrected", size = 0.4,do.label = TRUE) +
    ggtitle("Annotation on major cell types")
  
  p2 <- dittoDimPlot(spe_subset, var = "sub_celltype", reduction.use = "UMAP_mnnCorrected", size = 0.4,do.label = TRUE) +
    ggtitle("Annotation on sub cell types")
  
  # save the plot
  pdf(file.path(figureDir,"UMAP of Manual annotation of major and sub celltypes.pdf"),width = 18,height = 6)
  print(p1+p2)
  dev.off()
}

### Fraction of types in ROI
abundancedf <- as.data.frame(colData(spe))[,c("sample_id","sub_celltype")]

cluster_counts <- abundancedf %>%
  group_by(sample_id, sub_celltype) %>%
  summarise(count = n()) %>%
  group_by(sample_id) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(fraction = count / total)

fraction_wide <- cluster_counts %>%
  select(sample_id, sub_celltype, total, fraction) %>%
  pivot_wider(names_from = sub_celltype, values_from = fraction, values_fill = 0)

### Absolute number of types in ROI
total_counts <- abundancedf %>%
  group_by(sub_celltype) %>%
  summarise(total_count = n(), .groups = 'drop')

write.csv(as.data.frame(fraction_wide),file.path(figureDir,"Celltype_fraction_within_img.csv"),row.names = F)
write.csv(as.data.frame(total_counts),file.path(figureDir,"total_number_of_celltypes.csv"),row.names = F)

# "1" = "","2" = "","3" = "","4" = "","5" = "",
# "6" = "","7" = "","8" = "","9" = "","10" = "",
# "11" = "","12" = "","13" = "","14" = "","15" = "",
# "16" = "","17" = "","18" = "","19" = "","20" = "",
# "21" = "","22" = "","23" = "","24" = "","25" = "",
# "26" = "","27" = "","28" = "","29" = "","30" = "",
# "31" = "","32" = "","33" = "","34" = "","35" = "",
# "36" = "","37" = "","38" = "","39" = "","40" = "",
# "41" = "","42" = "","43" = "","44" = "","45" = "",
# "46" = "","47" = "","48" = "","49" = "","50" = ""
