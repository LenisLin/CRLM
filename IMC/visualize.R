library(SingleCellExperiment)

library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(pheatmap)
library(ggDoubleHeat)
library(ggrepel)
library(corrplot)

library(RColorBrewer)
library(ggsci)
library(circlize)

library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)

## Plot Heatmap for visualizing clustering
plot_cluster_Heatmap <- function(spe_, marker_use_col = "used_for_clustering", clustercol = "lcell_flowsom100pheno15", cutoff = FALSE) {
  if(class(rowData(spe_)[1,marker_use_col])  == class(1)){
    markes_use <- rowData(spe_)[,marker_use_col] == 1
  }
  if(class(rowData(spe_)[1,marker_use_col]) == class(TRUE)){
    markes_use <- rowData(spe_)[,marker_use_col]
  }
  
  marker_total <- rowData(spe_)[markes_use,]$Clean_Target
  
  transdata <-  data.frame(t(assay(spe_, "exprs")[marker_total,]), check.names = TRUE)
  transdata["metacluster"] <- colData(spe_)[clustercol][,1]
  colnames(transdata)
  
  tmp_count <- data.frame((table(transdata$metacluster)))

  heatmap_data <- transdata %>%
    group_by_at(c("metacluster")) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    data.frame(check.names = FALSE)
  
  if (cutoff != FALSE) {
    heatmap_data[, marker_total][(heatmap_data[, marker_total] > cutoff)] <- cutoff
  }
  
  htdf <- data.frame(t(heatmap_data[, marker_total]), check.names = FALSE)
  colname_ <- as.character(heatmap_data[, "metacluster"]) 
  for(i in 1:length(colname_)){
    colname_[i] <- paste0(colname_[i]," (n=",tmp_count[i,"Freq"],")")
  }
  colnames(htdf) <- colname_
  
  
  p1 <- pheatmap(data.frame(htdf, check.names = FALSE),
                 display_numbers = TRUE, cluster_rows = TRUE,
                 clustering_method = "average", angle_col = "90"
  )
  return(p1)
}

## Plot Heatmap for visualizing subtypes
plot_subtype_Heatmap <- function(spe_, markers) {
  # Get unique minor and major types
  MinorType <- unique(spe_$sub_celltype)
  MajorType <- sapply(MinorType, function(x) spe_$major_celltype[spe_$sub_celltype == x][1])
  
  # Create a dataframe with MajorType and SubType
  typedf <- data.frame("MajorType" = MajorType, "SubType" = MinorType)
  typedf <- typedf[order(typedf$MajorType), ]
  
  majortype <- unique(typedf$MajorType)
  
  # Assign color
  ann_colors <- list(
    setNames(MajorType <- c(brewer.pal(6, "YlOrRd")[3:5],
               brewer.pal(6, "PuBu")[3:6],
               brewer.pal(6, "YlGn")[3:5],
               brewer.pal(6, "BuPu")[3:6])[1:length(majortype)],
             majortype)
  )
  
  
  # Filter data to include only the specified markers
  data <- data.frame(t(assay(spe_, "exprs")[markers,]), check.names = TRUE)
  data$SubType <- spe_$sub_celltype
  
  # Compute mean expression for each marker within each subtype
  mean_expression <- data %>%
    pivot_longer(cols = markers, names_to = "Marker", values_to = "Expression") %>%
    group_by(SubType, Marker) %>%
    summarize(MeanExpression = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Marker, values_from = MeanExpression, values_fill = NA)
  
  # Add fraction information to the subtype names
  mean_expression$fraction <- sapply(mean_expression$SubType, function(st) {
    fraction <- nrow(data[data$SubType == st, ]) / nrow(data)
    paste0(st, " (", round(fraction, 3) * 100, "%)")
  })
  
  # Create a matrix for heatmap plotting
  mean_expression_matrix <- mean_expression %>%
    select(-SubType) %>%
    column_to_rownames(var = "fraction") %>%
    as.matrix()
  
  # Add MajorType column for annotation
  mean_expression$MajorType <- typedf$MajorType[match(mean_expression$SubType, typedf$SubType)]
  
  # Prepare row annotations
  annotationRow <- as.data.frame(mean_expression$MajorType)
  colnames(annotationRow) <- "MajorType"
  rownames(annotationRow) <- rownames(mean_expression_matrix)
  
  # Generate heatmap
  p <- pheatmap(
    mean_expression_matrix[, markers],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(colors = c("white", "red"))(100),
    annotation_row = annotationRow,
    annotation_colors = ann_colors
  )
  
  return(p)
}

plot_rich_heatmap <- function(spe_, celltype_column = "celltype", marker_class_column = "marker_class",
                              metadata_columns = c("indication", "patient_id"), spatial_columns = c("area", "n_neighbors")) {
  
  set.seed(22)
  
  ### 1. Heatmap bodies ###
  
  # Heatmap body color 
  col_exprs <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#440154FF", "#3B518BFF", "#20938CFF", "#6ACD5AFF", "#FDE725FF"))
  
  # Create Heatmap objects
  celltype_mean <- aggregateAcrossCells(as(spe_, "SingleCellExperiment"),  
                                        ids = spe_[[celltype_column]], 
                                        statistics = "mean",
                                        use.assay.type = "exprs", 
                                        subset.row = rownames(spe_)[rowData(spe_)[[marker_class_column]] == "type"])
  
  h_type <- Heatmap(t(assay(celltype_mean, "exprs")),
                    column_title = "type_markers",
                    col = col_exprs,
                    name = "mean exprs",
                    show_row_names = TRUE, 
                    show_column_names = TRUE)
  
  cellstate_mean <- aggregateAcrossCells(as(spe_, "SingleCellExperiment"),  
                                         ids = spe_[[celltype_column]], 
                                         statistics = "mean",
                                         use.assay.type = "exprs", 
                                         subset.row = rownames(spe_)[rowData(spe_)[[marker_class_column]] == "state"])
  
  h_state <- Heatmap(t(assay(cellstate_mean, "exprs")),
                     column_title = "state_markers",
                     col = col_exprs,
                     name = "mean exprs",
                     show_row_names = TRUE,
                     show_column_names = TRUE)
  
  ### 2. Heatmap annotation ###
  
  ### 2.1 Metadata features ###
  
  anno <- colData(celltype_mean) %>% as.data.frame %>% select(celltype = !!sym(celltype_column), ncells)
  
  # Proportion of indication per celltype
  indication <- unclass(prop.table(table(spe_[[celltype_column]], spe_[[metadata_columns[1]]]), margin = 1))
  
  # Number of contributing patients per celltype
  cluster_PID <- colData(spe_) %>% 
    as.data.frame() %>% 
    select(celltype = !!sym(celltype_column), patient_id = !!sym(metadata_columns[2])) %>% 
    group_by(celltype) %>% table() %>% 
    as.data.frame()
  
  n_PID <- cluster_PID %>% 
    filter(Freq > 0) %>% 
    group_by(celltype) %>% 
    count(name = "n_PID") %>% 
    column_to_rownames("celltype")
  
  # Create HeatmapAnnotation objects
  ha_anno <- HeatmapAnnotation(celltype = anno$celltype,
                               border = TRUE, 
                               gap = unit(1, "mm"),
                               col = list(celltype = metadata(spe_)$color_vectors[[celltype_column]]),
                               which = "row")
  
  ha_meta <- HeatmapAnnotation(n_cells = anno_barplot(anno$ncells, width = unit(10, "mm")),
                               n_PID = anno_barplot(n_PID, width = unit(10, "mm")),
                               indication = anno_barplot(indication, width = unit(10, "mm"),
                                                         gp = gpar(fill = metadata(spe_)$color_vectors[[metadata_columns[1]]])
                                                         ),
                               border = TRUE, 
                               annotation_name_rot = 90,
                               gap = unit(1, "mm"),
                               which = "row")
  
  ### 2.2 Spatial features ###
  
  # Add number of neighbors to spe object (saved in colPair)
  spe_$n_neighbors <- countLnodeHits(colPair(spe_, "neighborhood"))
  
  # Select spatial features and average over celltypes
  spatial <- colData(spe_) %>% 
    as.data.frame() %>% 
    select(all_of(spatial_columns), celltype = !!sym(celltype_column))
  
  spatial <- spatial %>% 
    select(-celltype) %>% 
    aggregate(by = list(celltype = spatial$celltype), FUN = mean) %>% 
    column_to_rownames("celltype")
  
  # Create HeatmapAnnotation object
  ha_spatial <- HeatmapAnnotation(area = spatial$area,
                                  n_neighbors = spatial$n_neighbors,
                                  border = TRUE,
                                  gap = unit(1, "mm"),
                                  which = "row")
  
  ### 3. Plot rich heatmap ###
  
  # Create HeatmapList object
  h_list <- h_type + h_state + ha_anno + ha_spatial + ha_meta
  
  # Add customized legend for anno_barplot()
  lgd <- Legend(title = "indication", 
                at = colnames(indication), 
                legend_gp = gpar(fill = metadata(spe_)$color_vectors[[metadata_columns[1]]]))
  
  # Return
  p <- draw(h_list, annotation_legend_list = list(lgd))
  return(p)
}

abudance_boxplot <- function(countdf, expCols,clinicalGroupCol = "TRG", test_ = "t.test", mean.by.patient = TRUE, return.mat = FALSE) {
  # Convert data to long format
  countdf$PID <- as.character(countdf$PID)
  
  if(mean.by.patient){
    countdf <- countdf %>%
      group_by(PID,group) %>%
      summarise(across(c(1:(ncol(countdf) - 2)), mean, na.rm = TRUE))
      }

  countdf <- as.data.frame(countdf)
  
  AbunBoxDF <- pivot_longer(countdf, cols = expCols, values_to = "Abundance", names_to = "Celltype")
  AbunBoxDF <- as.data.frame(AbunBoxDF)
  
  # Rename clinical group column and map clinical group values
  factor_vector <- factor(AbunBoxDF[,clinicalGroupCol])
  # Get the levels of the factor
  levels_vector <- levels(factor_vector)
  
  # Print the correspondence
  correspondence_df <- data.frame(
    Original_String = levels_vector,
    Factor_Level = seq_along(levels_vector)
  )
  print(correspondence_df)
  AbunBoxDF$ClinicalGroup <- factor_vector
  
  if(return.mat){
    return(AbunBoxDF)
  }
  
  # Create boxplot
  p <- ggplot(AbunBoxDF, aes(x = Celltype, y = Abundance, fill = ClinicalGroup)) +
    geom_boxplot(alpha = 0.7, color = "black", outlier.shape = NA) +
    scale_fill_manual(values = ggsci::pal_jco("default")(length(levels_vector))) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_blank()
    ) +
    stat_compare_means(aes(group = ClinicalGroup),
                       method = test_,
                       hide.ns = FALSE,
                       label = "..p..",
                       label.y.npc = "middle")
  
  return(p)
}

# Heatmap to visualize of cell interaction data
DoubleHeat <- function(data1, label1, group1, data2, label2, group2, plot = "circle", savePath) {
  # Multiply the data frames element-wise
  df1 <- data1 * label1
  df2 <- data2 * label2
  
  # Initialize an empty data frame
  plotdf <- matrix(data = 0, nrow = nrow(df1) * ncol(df1), ncol = 4)
  plotdf <- as.data.frame(plotdf)
  
  # Fill the data frame with the calculated values
  plotdf[, 1] <- rep(rownames(df1), times = ncol(df1))
  plotdf[, 2] <- rep(colnames(df1), each = nrow(df1))
  plotdf[, 3] <- as.numeric(as.matrix(df1))
  plotdf[, 4] <- as.numeric(as.matrix(df2))
  
  # Set the column names of the data frame
  colnames(plotdf) <- c("Celltype1", "Celltype2", "Interaction1", "Interaction2")
  
  # Check the specified plot type and create the corresponding plot
  if (plot == "heatmap") {
    p <- ggplot(plotdf, aes(x = Celltype1, y = Celltype2)) +
      geom_heat_tri(
        upper = Interaction1, lower = Interaction2,
        upper_name = c(group1), lower_name = c(group2),
        lower_colors = c("#075fd5", "white", "#fd6a78"),
        upper_colors = c("#075fd5", "white", "#fd6a78")
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  if (plot == "circle") {
    p <- ggplot(plotdf, aes(x = Celltype1, y = Celltype2)) +
      geom_heat_circle(
        outside = Interaction2,
        inside = Interaction1,
        outside_colors = c("#075fd5", "white", "#fd6a78"),
        inside_colors = c("#075fd5", "white", "#fd6a78")
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Save the plot as a PDF
  pdf(savePath, height = 6, width = 8)
  print(p)
  dev.off()
  
  return(NULL)
}

## 
corHeatmap <- function(M){
  col <- colorRampPalette(c("#BB4444", "#EE9988", 
                            "#FFFFFF", "#77AADD",
                            "#4477AA"))
  
  p <- corrplot(M, method = "color", col = col(200),
           #type = "upper", order = "hclust", 
           # addCoef.col = "black", # Add coefficient of correlation
           # tl.col="black", tl.srt = 45, # Text label color and rotation
           
           # Combine with significance
           #p.mat = p.mat, sig.level = 0.01, insig = "blank", 
           
           # hide correlation coefficient
           # on the principal diagonal
           diag = TRUE 
  )
  
  return(p)
}

## Define the Volcano Plot function
VolcanoPlot <- function(df, pthreshold = 0.05, fcthreshold = 1.4, clinicalFeatures, Qvalue = FALSE) {
  ## Fold change
  df$Foldchange <- as.numeric(df$Foldchange)
  df$P.value <- as.numeric(df$P.value)
  
  ## p value
  if (Qvalue) {
    df$Q.value <- p.adjust(df$P.value, method = "BH")
    df$change <- as.factor(ifelse(df$Q.value < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold),
                                  ifelse(log2(df$Foldchange) > log2(fcthreshold), "Up-regulate", "Down-regulate"), "Non-significant"
    ))
    
    ## label
    df$label <- ifelse(df[, "Q.value"] < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold), as.character(df[, 1]), "")
    
    ## plot
    p.vol <- ggplot(
      data = df,
      aes(x = log2(Foldchange), y = -log10(Q.value), colour = change, fill = change)
    ) +
      scale_color_manual(values = c("Down-regulate" = "blue", "Non-significant" = "grey", "Up-regulate" = "red")) +
      geom_point(alpha = 0.4, size = 3.5) +
      geom_text_repel(aes(x = log2(Foldchange), y = -log10(Q.value), label = label),
                      size = 3,
                      box.padding = unit(0.6, "lines"), point.padding = unit(0.7, "lines"),
                      segment.color = "black", show.legend = FALSE
      ) +
      geom_vline(xintercept = c(-(log2(fcthreshold)), (log2(fcthreshold))), lty = 4, col = "black", lwd = 0.8) +
      geom_hline(yintercept = -log10(pthreshold), lty = 4, col = "black", lwd = 0.8) +
      theme_bw() +
      labs(x = "log2(Fold Change)", y = "-log10(Q value)", title = paste0("Volcano Plot of Different Expression Markers in ", clinicalFeatures)) +
      theme(
        axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      )
  } else {
    df$change <- as.factor(ifelse(df$P.value < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold),
                                  ifelse(log2(df$Foldchange) > log2(fcthreshold), "Up-regulate", "Down-regulate"), "Non-significant"
    ))
    
    ## label
    df$label <- ifelse(df[, 3] < pthreshold & abs(log2(df$Foldchange)) > log2(fcthreshold), as.character(df[, 1]), "")
    
    ## plot
    p.vol <- ggplot(
      data = df,
      aes(x = log2(Foldchange), y = -log10(P.value), colour = change, fill = change)
    ) +
      scale_color_manual(values = c("Down-regulate" = "blue", "Non-significant" = "grey", "Up-regulate" = "red")) +
      geom_point(alpha = 0.4, size = 3.5) +
      geom_text_repel(aes(x = log2(Foldchange), y = -log10(P.value), label = label),
                      size = 3,
                      box.padding = unit(0.6, "lines"), point.padding = unit(0.7, "lines"),
                      segment.color = "black", show.legend = FALSE
      ) +
      geom_vline(xintercept = c(-(log2(fcthreshold)), (log2(fcthreshold))), lty = 4, col = "black", lwd = 0.8) +
      geom_hline(yintercept = -log10(pthreshold), lty = 4, col = "black", lwd = 0.8) +
      theme_bw() +
      labs(x = "log2(Fold Change)", y = "-log10(P value)", title = paste0("Volcano Plot of Different Expression Markers in ", clinicalFeatures)) +
      theme(
        axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      )
  }
  return(p.vol)
}

## Modify
HeatmapForDiff <- function(PvalueMat, FoldChangeMat, Sig = 0.05, savepath) {
  ## Set the FC less than cutoff as n.s.
  TempDF <- ifelse(FoldChangeMat == 0, 0, 1)
  PvalueMat <- PvalueMat * TempDF
  PvalueMat[is.na(PvalueMat)] <- 0
  FoldChangeMat[is.na(FoldChangeMat)] <- 0
  
  if (Sig == 0.05) {
    sig_mat <- matrix(sapply(PvalueMat, getSig.05), nrow = nrow(PvalueMat))
  }
  if (Sig == 0.01) {
    sig_mat <- matrix(sapply(PvalueMat, getSig.01), nrow = nrow(PvalueMat))
  }
  if (Sig == 0.1) {
    sig_mat <- matrix(sapply(PvalueMat, getSig.1), nrow = nrow(PvalueMat))
  }
  plotdf <- FoldChangeMat
  
  # Determine the range of fold changes in your data
  min_value <- min(plotdf)
  max_value <- max(plotdf)
  
  # Determine the absolute maximum value, either positive or negative
  abs_max <- max(abs(min_value), abs(max_value))
  
  # Create a symmetric sequence of breaks with 0 in the middle
  breaks <- seq(from = -abs_max, to = abs_max, length.out = 101)
  
  # Create custom color palette with blue for negative values, white for zero, and red for positive values
  color_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)
  
  # Create custom labels for the color legend
  legend_labels <- c("Response", "Non-Response")
  
  # Create the heatmap with custom breaks, color palette, and legend labels
  p <- pheatmap(
    plotdf,
    cellwidth = 16, cellheight = 12,
    cluster_row = F, cluster_col = F,
    angle_col = "90", display_numbers = sig_mat, fontsize_number = 15,
    breaks = breaks,
    color = color_palette,
    legend = TRUE,
    legend_labels = legend_labels,
    legend_breaks = c(min(breaks), max(breaks))
  )
  
  pdf(savepath, width = 10, height = 8)
  print(p)
  dev.off()
  
  return(NULL)
}

getSig.05 <- function(dc) {
  sc <- ""
  
  if (dc >= 3) {
    sc <- "***"
  } else if (dc >= 2) {
    sc <- "**"
  } else if (dc >= 1.3) {
    sc <- "*"
  }
  return(sc)
}

## display differential genes and neighbors
DEGs_pointPlot <- function(mat_foldchangeMat, pthreshold = 0.05, fcthreshold = 1.2 ){
  # Convert p-values to -log10(p-value) for better visualization
  mat_foldchangeMat$P.value <- as.numeric(mat_foldchangeMat$P.value)
  mat_foldchangeMat$Foldchange <- as.numeric(mat_foldchangeMat$Foldchange)
  mat_foldchangeMat$Celltype <- as.character(mat_foldchangeMat$Celltype)

  mat_foldchangeMat$neg_log10_pvalue <- -log10(mat_foldchangeMat$P.value)
  
  # Define colors for the plot
  mycol <- c("Up-Regulated" = "#E64B35FF", "Down-Regulated" = "#4DBBD5FF", "n.s." = "#808180FF")
  
  # Create a subset for significant genes to label
  mat_foldchangeMat$Label <- ifelse(
    mat_foldchangeMat$P.value <= pthreshold,
    ifelse(mat_foldchangeMat$Foldchange >= fcthreshold, "Up-Regulated",
           ifelse(mat_foldchangeMat$Foldchange <= 1 / fcthreshold, "Down-Regulated", "n.s.")),
    "n.s."
  )
  
  significant_genes <- mat_foldchangeMat[mat_foldchangeMat$Label != "n.s.", ]
  
  # Create the plot
  p1 <- ggplot(mat_foldchangeMat, aes(x = "Marker", y = Foldchange, color = Label)) +
    geom_jitter(aes(size = neg_log10_pvalue), width = 0.2, height = 0.2, alpha = 0.7) + # Use geom_jitter to avoid overplotting
    scale_color_manual(values = mycol) +
    xlab("Cell Type") +
    ylab("Fold Change") +
    theme_classic() +
    theme(
      legend.position = "right",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.text.x = element_text(colour = "black", size = 12),
      axis.text.y = element_text(colour = "black", size = 12),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold")
    ) +
    labs(size = "-log10(P-value)") +
    # Add labels for significant genes
    geom_text_repel(data = significant_genes, aes(label = Celltype), 
                    size = 4, max.overlaps = Inf)
  
  return(p1)
}
