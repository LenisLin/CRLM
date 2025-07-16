# For spatial analysis

# ------------------------------
# Set Working Directory and Source Utilities
# ------------------------------
library(SpatialExperiment)
library(imcRtools)
library(igraph)

library(ggplot2)
library(ggridges)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(ggraph)
library(viridis)
library(patchwork)

library(dplyr)
library(tidyr)

library(BiocParallel)
library(lisaClust)
library(scales)

workDir <- "~/IMCDataAnalysis/horizontalProject/IMC_xuanwu"

# Define directories for saving results and Steinbock outputs
setwd(workDir)
source("./code/utils.R")
source("./code/visualize.R")

saveDir <- file.path(workDir, "results")
steinResultDir <- file.path(workDir, "steinbock")
figureDir <- file.path(workDir, "figures","3_spatial")

if (!dir.exists(figureDir)) {
  dir.create(figureDir, recursive = T)
}

## prepare new spe data
if (F){
  ## Load clincal data
  clinical <- read.csv(file.path(workDir,"data","clinical.csv"),header = TRUE)
  
  ## subset
  spe_new <- spe[,spe$sample_id %in% clinical$ROI]
  
  ## combine clinical information
  clinical_features <- c('Pediatric', 'Recurrent_hemorrhage', 'Seizure', 'Brainstem')
  for(feature_ in clinical_features){
    colData(spe_new)[,feature_] <- clinical[match(spe_new$sample_id,clinical$ROI),feature_]
  }
  head(colData(spe_new))
  
  saveRDS(spe_new,file.path(saveDir, "spatial_spe_0518.rds"))
}

## Load data
spe <- readRDS(file.path(saveDir, "spatial_spe_0606.rds"))

compare_groups <- list(
  'Muta_GroupA'=list(c("1", "2"),c("1","3"),c("2","3")),
  'Muta_GroupB'=list(c("0", "1")),
  'Muta_GroupC'=list(c("0", "1")),
  'Bleeding'=list(c("0", "1")),
  'Pediatric'=list(c("1", "2")),
  'Recurrent_hemorrhage'=list(c("0", "1")),
  'Seizure'=list(c("0", "1")),
  'Brainstem'=list(c("0", "1"))
  )

# Cellular neighborhood analysis
img_id_ <- "sample_id"

##### Cellular neighborhood analysis
# Choose K clusters
k_clusters <- c(10, 15, 20)
Pairnames <- colPairNames(spe)

for(pairname_ in Pairnames){
  
  ## aggregate neighbor
  spe <- aggregateNeighbors(
    spe,
    colPairName = pairname_,
    aggregate_by = "metadata",
    count_by = "sub_celltype"
  )
  
  ## Iterative on number of clusters
  for(k_cluster_ in k_clusters){
    ## create folder
    figureDir_ <- file.path(figureDir,paste0(pairname_," ",k_cluster_))
    if(!file.exists(figureDir_)){
      dir.create(figureDir_,recursive = T)
    }
    
    cn <- kmeans(spe$aggregatedNeighbors, centers = k_cluster_)
    
    # 创建动态列名避免覆盖
    cluster_colname <- paste0(pairname_,"_cluster_", k_cluster_)
    colData(spe)[[cluster_colname]] <- as.factor(cn$cluster)
    spe$cn_celltypes <- as.factor(colData(spe)[[cluster_colname]])
    
    # 绘制热图
    for_plot <- prop.table(table(as.character(spe$cn_celltypes), spe$sub_celltype), margin = 1)
    p <- pheatmap(for_plot, 
                  color = colorRampPalette(c("dark blue", "white", "dark red"))(100),
                  scale = "column")
    
    pdf(file.path(figureDir_, paste0("CN analysis of ",pairname_," in ",k_cluster_," cluster.pdf")),width = 15,height = 10)
    print(p)
    dev.off()
    
    # 绘制组成条形图 
    if(T){
      # get the number of celltypes within different CN
      countdf_ <- Transform_CellCountMat(
        spe_ = spe,
        clinicalFeatures = "cn_celltypes",
        img_id = img_id_,
        count_by = "sub_celltype",  # 这里保持使用临时列
        is.fraction = FALSE
      )
      
      # Identify cell subpopulation columns by excluding "PID" and "cn_celltypes"
      cell_subpop_cols <- setdiff(names(countdf_), c("PID", "cn_celltypes"))
      
      # Sum cell subpopulation counts for each cn_celltypes
      summed_df <- countdf_ %>%
        group_by(cn_celltypes) %>%
        summarise(across(all_of(cell_subpop_cols), sum, na.rm = TRUE)) %>%
        ungroup()
      
      # Convert cn_celltypes to factor for categorical plotting
      summed_df$cn_celltypes <- factor(summed_df$cn_celltypes)
      
      # Transform to long format for ggplot2
      long_df <- summed_df %>%
        pivot_longer(cols = all_of(cell_subpop_cols), 
                     names_to = "cell_subpopulation", 
                     values_to = "count")
      
      long_df <- as.data.frame(long_df)
      long_df$cell_subpopulation <- as.factor(long_df$cell_subpopulation)
      
      # Create the stacked barplot with horizontal bars
      p <- ggplot(long_df, aes(x = cn_celltypes, y = count, fill = cell_subpopulation)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(x = "CN Cell Types", 
             y = "Total Cell Count", 
             fill = "Cell Subpopulation") +
        theme_minimal() +
        scale_fill_manual(values = metadata(spe)$color_vectors[["sub_celltype"]]) +
        coord_flip()
      
      pdf(file.path(figureDir_, paste0("CNAtype_stacked_barplot_k", k_cluster_, ".pdf")),
          width = 10,
          height = 7.5)
      print(p)
      dev.off()
    }
    
    # Save Matrix
    TRG_countdf <- Transform_CellCountMat(spe_ = spe,clinicalFeatures = c(names(compare_groups)),img_id = img_id_,count_by = "cn_celltypes",is.fraction = TRUE)
    TRG_countdf <- as.data.frame(TRG_countdf)
    
    write.csv(TRG_countdf,file.path(figureDir_, paste0("CN analysis of ",pairname_," in ",k_cluster_," cluster.csv")))
    
    # Compare within different clinical matrix
    for(compare_group_ in names(compare_groups)){
      plotdf <- TRG_countdf[,c(1:k_cluster_,match(compare_group_,colnames(TRG_countdf)))]
      
      # Data transform
      plotdf <- pivot_longer(data = plotdf,cols = 1:k_cluster_,names_to = "CN_type",values_to = "Fraction")
      plotdf <- as.data.frame(plotdf)
      plotdf <- na.omit(plotdf)
      colnames(plotdf)[1] <- "group"
      
      plotdf$Fraction <- as.numeric(plotdf$Fraction)
      plotdf$group <- as.factor(plotdf$group)
      
      # Plot boxplot
      p1 <- ggplot(plotdf, aes(x = group, y = Fraction, fill = group)) +
        geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
        geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), size = 0.1, alpha = 0.3) +
        scale_fill_manual(values = ggsci::pal_jco("default")(length(unique(plotdf$group)))) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 12),
          strip.background = element_blank()
        ) +
        ggtitle(paste0("CN analysis of ", pairname_, " in ", compare_group_)) +
        facet_wrap(~CN_type,scales = "free_y") +
        stat_compare_means(
          comparisons = compare_groups[[clinical.info]], 
          method = "wilcox.test",        # More robust than t-test
          p.adjust.method = "fdr",       # Less conservative than Bonferroni
          label = "p.adj",               # Shows adjusted p-values
          hide.ns = FALSE                 # Cleaner visualization
        )
      
      pdf(file.path(figureDir_, paste0("CN analysis of ",compare_group_," in ",pairname_," with ",k_cluster_," cluster.pdf")),width = 12,height = 15)
      print(p1)
      dev.off()
    }
  }
}

## Visualize the CN types on images
pairname_ <- "knn_20"
all_images <- unique(spe$sample_id)

num_figures <- 10
length_ <- length(all_images) %/% num_figures
for(i in 1:num_figures){
  sample_image <- all_images[c((length_*(i-1)):(length_*(i)))]
  
  #print(sample_image)
  spe_subset <- spe[,spe$sample_id %in% sample_image]
  
  ## Plot
  p <- plotSpatial(spe_subset,  # spe
                   node_color_by = "knn_20_cluster_10", img_id = "sample_id",node_size_fix = 0.15) + 
    scale_color_brewer(palette = "Set3")
  
  pdf(file.path(figureDir, paste0("CN Types on select images of figure id ",i,".pdf")),width = 12,height = 9)
  print(p)
  dev.off()
}

## Cell-cell pairwise interaction analysis
out <- readRDS(file = file.path(saveDir, "Interaction_analysis_out.rds"))
out <- out[out$group_by %in% spe$sample_id, ] ## subset tumor interaction

out$group <- spe$patient_id[match(out$group_by, spe$sample_id)]
head(out)

## overall interaction
df_ <- out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE))
p <- df_ %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  geom_text(aes(from_label, to_label, label = sum_sigval), size = 2) +
  scale_fill_gradient2(low = muted("#7aa6dcff"),
                       mid = "white",
                       high = muted("#cd534cff")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(figureDir, "celltype pair interaction heatmap in all images.pdf"),width = 12,height = 9)
print(p)
dev.off()

# Cellular interation between different clinical matrix
for(clinical_ in names(compare_groups)){
  clinical_groups <- unique(unlist(compare_groups[[clinical_]]) )
  
  for(clinical_group_ in clinical_groups){
    idx_ <- colData(spe)[,clinical_] %in% clinical_group_
    select_roi_ <- unique(spe[,idx_]$sample_id)
    
    out_ <- out[out$group_by %in% select_roi_, ] %>% 
      as_tibble() %>%
      group_by(from_label, to_label) %>%
      summarize(sum_sigval = sum(sigval, na.rm = TRUE))
    
    p <- out_ %>%
      ggplot() +
      geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
      geom_text(aes(from_label, to_label, label = sum_sigval), size = 2) +
      scale_fill_gradient2(low = muted("#7aa6dcff"),
                           mid = "white",
                           high = muted("#cd534cff")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    pdf(file.path(figureDir,paste0("Pair interaction in ",clinical_,"_group_",clinical_group_,".pdf")) ,width = 12,height = 9)
    print(p)
    dev.off()
  }
}

## Display specific cell-cell interaction
print(unique(out$from_label))
sourceType <- "YAP+ Fibro"
group_ <- "Group_A"

out_ <- out[out$from_label %in% sourceType,]
out_ <- as.data.frame(out_[out_$group %in% group_,]) 

if(T){
  # Step 1: Sum sigval for each to_label across all ROIs
  sum_sigval <- out_ %>%
    group_by(to_label) %>%
    summarise(total_sigval = sum(sigval, na.rm = TRUE)) %>%
    ungroup()
  
  # Step 2: Compute the fraction
  total_rois <- length(unique(out_$group_by))
  sum_sigval <- sum_sigval %>%
    mutate(fraction = total_sigval / total_rois)
  
  # Step 3: Sort by fraction in decreasing order
  sum_sigval <- sum_sigval %>%
    arrange(desc(fraction))
  
  # Step 4: Set factor levels and add sign column
  sum_sigval$to_label <- factor(sum_sigval$to_label, levels = sum_sigval$to_label)
  sum_sigval <- sum_sigval %>%
    mutate(sign = ifelse(fraction > 0, "Positive", "Negative"))
  
  # Step 5: Create the barplot
  p <- ggplot(sum_sigval, aes(x = to_label, y = fraction, fill = sign)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Positive" = "#cd534cff", "Negative" = "#7aa6dcff")) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "To Label", y = "Fraction", title = "Fraction of Significant Cell-Cell Interactions")
  
  pdf(file.path(figureDir,paste0("celltype pair interaction barplot in ",group_,".pdf")) ,width = 8,height = 6)
  print(p)
  dev.off()
}

# ------------------------------------------------------------------------------------------------
##### Patch Analysis
colnames(colData(spe)) ## view column of data
patch_types <- colnames(colData(spe))[endsWith(colnames(colData(spe)),suffix = "_patch")] ## get all patch type
patch_names <- sapply(patch_types, function(x){
  return(strsplit(x,"_")[[1]][1])
})
all_images <- unique(spe$sample_id)

figureDir <- file.path(workDir, "figures", "4_Patch")
if(!file.exists(figureDir)){
  dir.create(figureDir,recursive = T)
}

## Plot different types of patches
custom_palette <- c(
  "EndoMT_patch" = "#A65628",
  "Endothelial_patch" = "#4DAF4A"
)

## Plot Patch on spatial
for(patch_type_ in patch_types){
  
  ## Sample for plotting
  sample_rois_ <- sample(all_images,12)
  spe_temp <- spe[,spe$sample_id %in% sample_rois_]
  
  patch_type_vec_ <- colData(spe_temp)[,patch_type_] ## get correspond patch type
  colData(spe_temp)[,patch_type_] <- ifelse(is.na(patch_type_vec_),NA,patch_type_)
  
  p <- plotSpatial(
    spe_temp,
    node_color_by = patch_type_,
    img_id = "sample_id",
    node_size_fix = 0.5
  ) +
    theme(legend.position = "right") +
    scale_color_manual(values = custom_palette)
  
  rm(spe_temp,sample_rois_,patch_type_vec_)
  gc()
  
  pdf(file.path(figureDir,paste0(patch_type_," in random ROIs.pdf")), width = 24,height = 18)
  print(p)
  dev.off()
}


### Analysis configuration in different patch types
analysisDF <- as.data.frame(colData(spe)[,c("sample_id","sub_celltype",names(compare_groups),patch_types,unname(patch_names))])
analysisDF_Back <- analysisDF

for(patch_type_ in patch_types){
  
  ## Get patch name
  patch_name_ <- strsplit(patch_type_,"_")[[1]]
  patch_name_ <- paste(patch_name_[1:length(patch_name_)-1],collapse =  "_")
  
  ## Get metadata
  analysisDF <- analysisDF[,match(c("sample_id","sub_celltype",names(compare_groups),patch_type_,patch_name_),
                                  colnames(analysisDF))]
  colnames(analysisDF)[c(11,12)] <- c("patch_type","patch_name")
  
  # --- 1. Patch Frequency and Size Analysis ---
  if(T){
    # Step 1: Calculate total number of cells per sample_id
    total_cells <- analysisDF %>%
      group_by(sample_id) %>%
      summarise(total_cells = n(), .groups = 'drop')
    
    # Step 2: Calculate number of cells with non-NA patch_type per sample_id
    patch_cells <- analysisDF %>%
      filter(!is.na(patch_type)) %>%
      group_by(sample_id) %>%
      summarise(patch_cells = n(), .groups = 'drop')
    
    # Step 3: Combine and calculate the ratio (patch_cells / total_cells)
    cell_ratio <- total_cells %>%
      left_join(patch_cells, by = "sample_id") %>%
      mutate(
        patch_cells = ifelse(is.na(patch_cells), 0, patch_cells),
        ratio = patch_cells / total_cells
      ) %>%
      arrange(desc(ratio))  # Sort in decreasing order
    
    # Print summary statistics
    cat("Summary of cell counts and ratios:\n")
    print(head(cell_ratio, 10))
    cat("\nTotal samples:", nrow(cell_ratio), "\n")
    cat("Samples with patch_type data:", sum(cell_ratio$patch_cells > 0), "\n")
    
    # Step 4 & 5: Create decreasing bar plot with top-5 sample labels
    # Identify top 5 samples with non-zero ratios for labeling
    top5_samples <- cell_ratio %>%
      filter(ratio > 0) %>%
      slice_head(n = 10) %>%
      pull(sample_id)
    
    # Create the lollipop plot
    p <- ggplot(cell_ratio, aes(x = reorder(sample_id, ratio), y = ratio)) +
      # Add the lollipop stick
      geom_segment(aes(xend = reorder(sample_id, ratio), yend = 0), 
                   color = "grey85", linewidth = 0.3) +
      # Add the lollipop circle
      geom_point(color = "#ffccccff", size = 4) +
      # Add text labels inside circles for top 5 samples
      geom_text(
        data = cell_ratio %>% filter(sample_id %in% top5_samples),
        #data = cell_ratio,
        aes(label = round(ratio,digits = 4)),
        color = "grey30",
        size = 2.5,
        fontface = "bold"
      ) +
      coord_flip() +  # Flip coordinates for better readability
      labs(
        title = "Proportion of Cells with Patch Type Information by Sample",
        x = "Sample ID",
        y = "Proportion of Cells with Patch Type",
        caption = paste("Total samples:", nrow(cell_ratio), "| Total cells:", sum(cell_ratio$total_cells))
      ) +
      theme_classic2() +
      theme(
        # axis.text.y = element_blank(),  # Hide y-axis labels 
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
        panel.grid.minor.y = element_blank()
      )
    
    # Display the plot
    pdf(file = file.path(figureDir,paste0(patch_type_," Patch size distribution.pdf")),width = 6,height = 12)
    print(p)
    dev.off()
    
  }
  patch_freq_roi <- analysisDF %>%
    filter(!is.na(patch_type)) %>%                  # Filter cells in patches
    group_by(sample_id) %>%
    summarise(n_patches = n_distinct(patch_type)) %>%  # Count unique patches
    left_join(analysisDF %>% select(sample_id, names(compare_groups)) %>% distinct(), by = "sample_id")
  
  patch_freq_roi$ratio <- cell_ratio[match(patch_freq_roi$sample_id,cell_ratio$sample_id),]$ratio 
  write.csv(patch_freq_roi,file.path(figureDir,paste0(patch_type_," Patch number distribution(ROI).csv")))
  
  # Boxplot by clinical group
  for(clinical.info in names(compare_groups)){
    patch_freq_roiTemp <- patch_freq_roi[,c("sample_id","n_patches","ratio",clinical.info)]
    patch_freq_roiTemp <- as.data.frame(na.omit(patch_freq_roiTemp)) 
    patch_freq_roiTemp$group <- as.factor(patch_freq_roiTemp[,clinical.info])
    
    p1 <- ggplot(patch_freq_roiTemp, aes(x = group, y = n_patches, fill = group)) +
      geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
      geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), size = 0.1, alpha = 0.3) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        strip.background = element_blank()
      ) +
      labs(title = paste0("Number of ",patch_type_," by ",clinical.info), x = "Group", y = "Number of Patches")+
      scale_fill_manual(values = ggsci::pal_jco("default")(length(unique(patch_freq_roiTemp$group)))) +
      stat_compare_means(
        comparisons = compare_groups[[clinical.info]], 
        method = "wilcox.test",        # More robust than t-test
        p.adjust.method = "fdr",       # Less conservative than Bonferroni
        label = "p.adj",               # Shows adjusted p-values
        hide.ns = FALSE                 # Cleaner visualization
      )
    
    p2 <- ggplot(patch_freq_roiTemp, aes(x = group, y = ratio, fill = group)) +
      geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
      geom_jitter(color = "darkgrey", position = position_jitter(width = 0.2), size = 0.1, alpha = 0.3) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        strip.background = element_blank()
      ) +
      labs(title = paste0("Fraction of ",patch_type_," by ",clinical.info), x = "Group", y = "Fraction of Patches")+
      scale_fill_manual(values = ggsci::pal_jco("default")(length(unique(patch_freq_roiTemp$group)))) +
      stat_compare_means(
        comparisons = compare_groups[[clinical.info]], 
        method = "wilcox.test",        # More robust than t-test
        p.adjust.method = "fdr",       # Less conservative than Bonferroni
        label = "p.adj",               # Shows adjusted p-values
        hide.ns = FALSE                 # Cleaner visualization
      )
    
    pdf(file = file.path(figureDir,paste0(patch_type_," boxplot for Patch number and size comparison in ",clinical.info,".pdf")),width = 10,height = 5)
    print(p1 + p2)
    dev.off()
  }
  
  # --- 2. Patch Cellular Composition ---  
  # Microenvironment composition
  microenv_comp <- analysisDF %>%
    filter(!is.na(patch_type) & !patch_name) %>%  # Non-CD45_TC cells in patches
    group_by(sample_id, patch_type, sub_celltype) %>%
    summarise(n_cells = n()) %>%
    group_by(sample_id, patch_type) %>%
    mutate(prop = n_cells / sum(n_cells)) %>%     # Proportion per patch
    group_by(sample_id, sub_celltype) %>%
    summarise(avg_prop = mean(prop)) %>%          # Average across patches
    left_join(analysisDF %>% select(sample_id, names(compare_groups)) %>% distinct(), by = "sample_id")
  
  for(clinical.info in names(compare_groups)){
    # Aggregate by clinical group
    microenv_compTemp <- microenv_comp[,c("sample_id","sub_celltype","avg_prop",clinical.info)]
    colnames(microenv_compTemp)[4] <- "group"
    microenv_compTemp <- microenv_compTemp[!is.na(microenv_compTemp$group),]
    
    microenv_group <- microenv_compTemp %>%
      group_by(group, sub_celltype) %>%
      summarise(mean_prop = mean(avg_prop))
    
    microenv_group$group <- as.factor(microenv_group$group)
    
    # Stacked bar plot for microenvironment composition
    p <- ggplot(microenv_group, aes(x = group, y = mean_prop, fill = sub_celltype)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_classic2()+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
      labs(title = "Microenvironment Composition in ",patch_name_," Patches", x = "Group", y = "Mean Proportion")+
      scale_fill_manual(values = metadata(spe)$color_vectors[["sub_celltype"]])
    
    pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," microenvironment composition.pdf")),width = 5,height = 8)
    print(p)
    dev.off()
    
    #write.csv(microenv_comp,file.path(figureDir,paste0(patch_type_," microenvironment composition.csv")))
    
    ## Compare across groups of subpopulation
    groups <- unique(microenv_compTemp$group)
    
    # Loop through each clinical group
    diff_results <- data.frame()
    for (grp in groups) {
      # Create a binary column: 1 if in the group, 0 otherwise
      microenv_compTemp$in_group <- ifelse(microenv_compTemp$group == grp, 1, 0)
      
      # Get unique cell subpopulations
      cell_types <- unique(microenv_compTemp$sub_celltype)
      
      # Perform Wilcoxon test for each cell type
      for (cell_type in cell_types) {
        # Filter data for the current cell type
        cell_data <- microenv_compTemp %>% filter(sub_celltype == cell_type)
        
        if(!(grp %in% unique(cell_data$group))){
          next;
        }
        if(length(unique(cell_data$group))==1){
          next;
        }
        # Perform Wilcoxon test (non-parametric test for comparing two groups)
        test_result <- wilcox.test(avg_prop ~ in_group, data = cell_data)
        
        # Calculate fold change (mean in group / mean in others)
        mean_in_group <- mean(cell_data$avg_prop[cell_data$in_group == 1], na.rm = TRUE)
        mean_others <- mean(cell_data$avg_prop[cell_data$in_group == 0], na.rm = TRUE)
        fold_change <- mean_in_group / mean_others
        log2_fc <- log2(fold_change)  # Log2 fold change for direction and magnitude
        
        # Store results
        diff_results <- rbind(diff_results, data.frame(
          group = grp,
          major_celltype = cell_type,
          p_value = test_result$p.value,
          fold_change = fold_change,
          log2_fc = log2_fc
        ))
      }
    }
    
    # Adjust p-values for multiple testing using Benjamini-Hochberg correction
    diff_results <- diff_results %>%
      group_by(group) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      ungroup()
    
    # Determine significance and direction
    p_sig_cutoff <- 0.05
    logfc_cutoff <- 1
    
    # Create the bubble plot
    if(T){
      diff_results <- diff_results %>%
        mutate(
          significance = case_when(
            p_value <= p_sig_cutoff & log2_fc >= logfc_cutoff ~ "Up",
            p_value <= p_sig_cutoff & log2_fc <= logfc_cutoff ~ "Down",
            TRUE ~ "NS"  # Not significant
          )
        )
      
      p <- ggplot(diff_results, aes(x = log2_fc, y = -log10(p_adj))) +
        # Plot all points in grey as background
        geom_point(color = "grey", size = 0.8) +
        # Overlay significant points in color
        geom_point(data = subset(diff_results, significance %in% c("Up", "Down")), 
                   aes(color = significance), size = 0.8) +
        # Set color scale for significant points21
        scale_color_manual(values = c("Up" = "#cd534cff", "Down" = "#7ca6dcff", "NS" = "grey")) +
        # Facet by group, one row with multiple columns
        facet_grid(. ~ group) +
        # Add threshold lines: vertical for log2_fc, horizontal for p_adj = 0.05
        geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey50", size = 0.5) +
        geom_hline(yintercept = -log10(p_sig_cutoff), linetype = "dashed", color = "grey50", size = 0.5) +
        # Use a clean theme similar to the new code
        theme_bw() +
        theme(
          panel.grid = element_blank(),           # Remove grid lines
          axis.text = element_text(size = 10),    # Axis text size from new code
          strip.text.x = element_text(size = 10, face = "bold"),  # Bold facet titles
          legend.position = "bottom",             # Legend at bottom, inspired by your code
          axis.title = element_text(face = "bold", size = 14),    # Bold titles from your code
          text = element_text(size = 12)          # General text size from your code
        ) +
        # Add labels and title
        labs(
          title = paste0("Differential Abundance of Cell Subpopulations by ",clinical.info),
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted P-value",
          color = "Significance"
        ) +
        # Add labels for significant points, matching your original approach
        geom_text_repel(
          data = subset(diff_results, significance %in% c("Up", "Down")),
          aes(label = major_celltype),
          size = 2,
          box.padding = 0.5,
          max.overlaps = Inf
        )
      
      pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," cell subpopulation fraction difference.pdf")),width = 8,height = 6)
      print(p)
      dev.off()
      
      # write.csv(diff_results,file.path(figureDir,paste0(patch_type_," cell subpopulation fraction difference.csv") ))
    }
    
  }
    # --- 3. Patch Cellular Expression ---  
  if(T){
    # Step 1: Filter cells by patch_type
    spe_filtered <- spe[,rownames(analysisDF)]
    
    # Step 2: Get expression matrix and metadata
    expr_matrix_all <- assay(spe_filtered)
    
    # Step 3: Calculate mean expression per celltype per sample
    
    for(clinical.info in names(compare_groups)){
      
      ## Get subset expression matrix and metadata
      cell_meta <- analysisDF[,c("sample_id","sub_celltype",clinical.info,"patch_type","patch_name")]
      cell_meta <- cell_meta[!is.na(cell_meta$patch_type),]
      cell_meta <- cell_meta[!is.na(cell_meta[,clinical.info]),]
      
      expr_matrix <- expr_matrix_all[,match(rownames(cell_meta),colnames(expr_matrix_all))]
      
      # Calculate mean expression for each marker in each group
      mean_expr_list <- list()
      
      for (i in 1:nrow(expr_matrix)) {
        marker_data <- data.frame(
          expression = expr_matrix[i, ],
          sample_id = cell_meta$sample_id,
          sub_celltype = cell_meta$sub_celltype,
          clinical_group = cell_meta[[clinical.info]],
          stringsAsFactors = FALSE
        )
        
        # Calculate mean expression per celltype per sample per clinical group
        marker_summary <- marker_data %>%
          group_by(sample_id, sub_celltype, clinical_group) %>%
          summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = 'drop') %>%
          mutate(marker = rownames(expr_matrix)[i])
        
        mean_expr_list[[i]] <- marker_summary
      }
      
      # Combine all markers
      all_mean_expr <- do.call(rbind, mean_expr_list)
      group_vals <- unique(all_mean_expr$clinical_group)
      n_groups <- length(group_vals)
      
      # Step 4: Perform statistical testing and calculate fold change
      results_list <- list()
      
      celltypes <- unique(all_mean_expr$sub_celltype)
      markers <- unique(all_mean_expr$marker)
      markers <- markers[!startsWith(markers,prefix = "DNA")]
      
      for (celltype in celltypes) {
        for (marker in markers) {
          # Get data for this celltype and marker
          subset_data <- all_mean_expr %>%
            filter(sub_celltype == celltype, marker == !!marker)
          
          if (nrow(subset_data) < n_groups * 3) next  # Need at least 3 samples per group
          
          # Check if we have data for all groups
          groups_present <- unique(subset_data$clinical_group)
          if (length(groups_present) < 2) next
          
          if (n_groups == 2) {
            # Two groups: use t-test (original logic)
            group1_data <- subset_data %>% filter(clinical_group == group_vals[1]) %>% pull(mean_expr)
            group2_data <- subset_data %>% filter(clinical_group == group_vals[2]) %>% pull(mean_expr)
            
            if (length(group1_data) == 0 || length(group2_data) == 0) next
            
            # Calculate fold change (group1 vs group2)
            mean_group1 <- mean(group1_data, na.rm = TRUE)
            mean_group2 <- mean(group2_data, na.rm = TRUE)
            fold_change <- log2((mean_group1 + 0.001) / (mean_group2 + 0.001))
            
            # Perform t-test
            if (length(group1_data) > 1 && length(group2_data) > 1) {
              test_result <- t.test(group1_data, group2_data)
              p_value <- test_result$p.value
            } else {
              p_value <- 1
            }
            
            results_list[[paste(celltype, marker, sep = "_")]] <- data.frame(
              sub_celltype = celltype,
              marker = marker,
              fold_change = fold_change,
              p_value = p_value,
              mean_group1 = mean_group1,
              mean_group2 = mean_group2,
              n_group1 = length(group1_data),
              n_group2 = length(group2_data),
              comparison = paste(group_vals[1], "vs", group_vals[2])
            )
            
          } else {
            # Multiple groups: use "1 vs. remaining" approach
            
            # Perform each group vs. all others
            for (target_group in group_vals) {
              
              # Split data into target group vs. all others
              target_data <- subset_data %>% filter(clinical_group == target_group) %>% pull(mean_expr)
              other_data <- subset_data %>% filter(clinical_group != target_group) %>% pull(mean_expr)
              
              if (length(target_data) == 0 || length(other_data) == 0) next
              
              # Calculate fold change (target vs others)
              mean_target <- mean(target_data, na.rm = TRUE)
              mean_others <- mean(other_data, na.rm = TRUE)
              fold_change <- log2((mean_target + 0.001) / (mean_others + 0.001))
              
              # Perform t-test: target group vs. all others
              if (length(target_data) > 1 && length(other_data) > 1) {
                test_result <- t.test(target_data, other_data)
                p_value <- test_result$p.value
              } else {
                p_value <- 1        
              }
              
              comparison_name <- paste(celltype, marker, target_group, "vs_others", sep = "_")
              results_list[[comparison_name]] <- data.frame(
                sub_celltype = celltype,
                marker = marker,
                fold_change = fold_change,
                p_value = p_value,
                mean_group1 = mean_target,
                mean_group2 = mean_others,
                n_group1 = length(target_data),
                n_group2 = length(other_data),
                comparison = paste(target_group, "vs Others")
              )
            }
          }
        }
      }
      
      # Combine results
      if (length(results_list) == 0) {
        cat("No valid comparisons found\n")
        next
      }
      
      results_df <- do.call(rbind, results_list)
      results_df$p_adj <- p.adjust(results_df$p_value, method = "BH") # Adjust p-values
      
      # Classify changes        
      results_df$change_type <- ifelse(results_df$p_adj <= 0.05 & results_df$fold_change >= logfc_cutoff, "up",
                                       ifelse(results_df$p_adj <= 0.05 & results_df$fold_change <= -logfc_cutoff, "down", 
                                              "n.s."))
      
      # Calculate transparency (-log10(p_adj))
      results_df$neg_log_p <- -log10(results_df$p_adj + 1e-10)  # Add small value to avoid infinite
      
      # Step 5: Create bubble plot
      if (n_groups == 2) {
        # Simple plot for two groups
        p <- ggplot(results_df, aes(x = sub_celltype, y = marker, 
                                    size = abs(fold_change), 
                                    color = change_type,
                                    alpha = neg_log_p)) +
          geom_point() +
          scale_color_manual(values = c("up" = "#cd534cff", "down" = "#7ca6dcff", "n.s." = "gray"),name = "Change") +
          scale_size_continuous(name = "Abs(log2FC)", range = c(1, 8)) +
          scale_alpha_continuous(name = "-log10(p.adj)", range = c(0.3, 1)) +
          theme_classic2() +
          theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "lightgray"),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5)
          ) +
          labs(title = paste("Differential Expression:",patch_name_ ),x = "Cell Type",y = "Marker") +
          guides(
            color = guide_legend(override.aes = list(size = 5, alpha = 1)),
            size = guide_legend(override.aes = list(alpha = 1)),
            alpha = guide_legend(override.aes = list(size = 5))
          )
        
        pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," cellular expression change.pdf")),width = 10,height = 7.5)
        print(p)
        dev.off()
      } 
      else {
        # Multiple comparison plot - create faceted plot by comparison
        p <- ggplot(results_df, aes(x = sub_celltype, y = marker, 
                                    size = abs(fold_change), 
                                    color = change_type,
                                    alpha = neg_log_p)) +
          geom_point() +
          facet_wrap(~comparison, scales = "free") +
          scale_color_manual(values = c("up" = "#cd534cff", "down" = "#7ca6dcff", "n.s." = "gray"),name = "Change") +
          scale_size_continuous(name = "Abs(log2FC)", range = c(1, 6)) +
          scale_alpha_continuous(name = "-log10(p.adj)", range = c(0.3, 1)) +
          theme_classic2() +
          theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "lightgray"),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
            axis.text.y = element_text(size = 7),
            plot.title = element_text(hjust = 0.5),
            strip.text = element_text(size = 8)
          ) +
          labs(title = paste("One vs Others Differential Expression:",patch_name_ ),x = "Cell Type",y = "Marker") +
          guides(
            color = guide_legend(override.aes = list(size = 5, alpha = 1)),
            size = guide_legend(override.aes = list(alpha = 1)),
            alpha = guide_legend(override.aes = list(size = 5)))
        
        pdf(file = file.path(figureDir,paste0(patch_type_," of ",clinical.info," one_vs_others cellular expression change.pdf")),width = 15,height = 7.5)
        print(p)
        dev.off()
      }
    }
  }

  analysisDF <- analysisDF_Back
}

# --- 3. Overlap and Network Analysis of Patches ---
# Assuming analysisDF is your data frame
# Identify patch columns
patch_cols <- grep("_patch$", names(analysisDF), value = TRUE)

# Define short labels for patch types
patch_labels <- c(
  "CD45_TC_patch" = "CD45",
  "YAP_TC_DTC_patch" = "YAP_DTC",
  "PD1_T_patch" = "PD1_T",
  "Treg_patch" = "Treg",
  "PDL1_Macro_patch" = "PDL1_M",
  "Endo_patch" = "Endo",
  "MC_2_patch" = "MC_2",
  "YAP_Fibro_patch" = "YAP_Fibro"
)
all_vertices <- unique(patch_labels)

# Get unique clinical groups
groups <- unique(analysisDF$group)

# Function to create a graph for a given group
create_group_graph <- function(df_grp, patch_cols, patch_labels) {
  # Get all pairwise combinations of patch columns
  pairs <- combn(patch_cols, 2, simplify = FALSE)
  
  # Calculate intersection sizes for each pair
  edge_list <- lapply(pairs, function(pair) {
    patch1 <- pair[1]
    patch2 <- pair[2]
    # Intersection size: count where both patches are non-NA
    weight <- sum(!is.na(df_grp[[patch1]]) & !is.na(df_grp[[patch2]]))
    c(patch_labels[patch1], patch_labels[patch2], weight)
  })
  
  # Convert to data frame
  edge_df <- do.call(rbind, edge_list)
  edge_df <- as.data.frame(edge_df)
  names(edge_df) <- c("from", "to", "weight")
  edge_df$weight <- as.numeric(edge_df$weight)
  
  # Filter out edges with zero weight
  edge_df <- edge_df %>% filter(weight > 0)
  
  # Create graph, ensuring all vertices are included (even isolates)
  g <- graph_from_data_frame(edge_df, directed = FALSE, 
                             vertices = data.frame(name = all_vertices))
  return(g)
}

# Create a list of graphs, one for each group
group_graphs <- lapply(groups, function(grp) {
  df_grp <- analysisDF %>% filter(group == grp)
  create_group_graph(df_grp, patch_cols, patch_labels)
})
names(group_graphs) <- groups
saveRDS(groups,file.path(figureDir,paste0("Network Analysis of Patches plot data.rds")))

# Visualize each graph using ggraph
plot_list <- lapply(groups, function(grp) {
  g <- group_graphs[[grp]]
  
  # Create the ggraph plot
  ggraph(g, layout = "circle") +
    # Edges with width proportional to intersection size
    geom_edge_link(aes(edge_width = weight), edge_colour = "grey50", alpha = 0.7) +
    # Nodes as points
    geom_node_point(size = 12, color = "lightblue", alpha = 0.9) +
    # Node labels
    geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
    # Scale edge width for better visualization
    scale_edge_width(range = c(0.5, 5), name = "Intersection Size") +
    # Clean theme
    theme_void() +
    # Titles and legend
    labs(title = paste("Group:", grp)) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )})

p <- wrap_plots(plot_list) + plot_layout(ncol = 2)

pdf(file = file.path(figureDir,"Patch Interaction.pdf"),width = 20,height = 16)
print(p)
dev.off()