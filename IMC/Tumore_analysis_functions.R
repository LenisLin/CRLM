library(SpatialExperiment)
library(imcRtools)
library(igraph)

library(ggplot2)
library(ggridges)
library(ggalluvial)
library(RColorBrewer)
library(ggpubr)
library(ggraph)
library(viridis)
library(patchwork)

library(dplyr)
library(tidyr)

# Function to compute subtype fraction for each ROI
compute_subtype_fractions <- function(subtype) {
  tc_im_mali_cells <- tc_im_cells[startsWith(tc_im_cells$sub_celltype,prefix = "EC"),]
  
  subtype_summary <- tc_im_mali_cells %>%
    group_by(sample_id, patient_id, Tissue, Treatment, RFS_status) %>%
    summarise(
      total_cells = n(),
      subtype_cells = sum(sub_celltype == subtype),
      subtype_fraction = subtype_cells / total_cells,
      .groups = "drop"
    )
  
  # Add early relapse labels
  subtype_summary$RFS_group <- ifelse(subtype_summary$RFS_status == 0, "No Early Relapse", "Early Relapse")
  
  return(subtype_summary)
}

# Function to create plots for each subtype
create_subtype_plots <- function(subtype, subtype_data, figureDir) {
  
  # Clean subtype name for file names
  clean_name <- gsub("[^A-Za-z0-9]", "_", subtype)
  
  cat(paste("Analyzing subtype:", subtype, "\n"))
  
  # Step 4: Stacked barplot
  p1_sub <- ggplot(subtype_data, aes(x = sample_id, y = subtype_fraction)) +
    geom_col(aes(fill = Tissue), alpha = 0.7) +
    facet_wrap(~paste(Treatment, "-", RFS_group), scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
    labs(title = paste(subtype, "Fraction by ROI"),
         x = "ROI (sample_id)", 
         y = paste(subtype, "Fraction"),
         fill = "Tissue Region") +
    scale_fill_manual(values = metadata(spe)$color_vectors[["tissue"]])
  
  print(p1_sub)
  ggsave(file.path(figureDir, paste0(clean_name, "_barplot.pdf")), p1_sub, width = 14, height = 8)
  
  # Step 5: Boxplot by RFS status
  comparison_data_sub <- subtype_data %>%
    select(sample_id, patient_id, Tissue, Treatment, RFS_status, RFS_group, subtype_fraction)
  
  p2_sub <- ggplot(comparison_data_sub, aes(x = RFS_group, y = subtype_fraction)) +
    geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_wrap(~Tissue) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_minimal() +
    labs(title = paste(subtype, "Fraction by Recurrence Status"),
         x = "Recurrence Status",
         y = paste(subtype, "Fraction"),
         fill = "RFS Status") +
    scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))
  
  print(p2_sub)
  ggsave(file.path(figureDir, paste0(clean_name, "_boxplot_rfs.pdf")), p2_sub, width = 10, height = 6)
  
  # Step 6: Boxplot by treatment and RFS status
  p3_sub <- ggplot(comparison_data_sub, aes(x = RFS_group, y = subtype_fraction)) +
    geom_boxplot(aes(fill = RFS_group), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_grid(Treatment ~ Tissue) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    theme_minimal() +
    labs(title = paste(subtype, "Fraction by Treatment and Recurrence Status"),
         x = "Recurrence Status",
         y = paste(subtype, "Fraction"),
         fill = "RFS Status") +
    scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))
  
  print(p3_sub)
  ggsave(file.path(figureDir, paste0(clean_name, "_boxplot_treatment.pdf")), p3_sub, width = 10, height = 8)
  
  # Step 7: Sankey diagram
  # Get paired TC-IM data for this subtype
  patient_summary_sub <- subtype_data %>%
    group_by(patient_id, Tissue, Treatment, RFS_status, RFS_group) %>%
    summarise(mean_subtype_fraction = mean(subtype_fraction, na.rm = TRUE), .groups = "drop")
  
  paired_data_sub <- patient_summary_sub %>%
    select(patient_id, Tissue, Treatment, RFS_status, RFS_group, mean_subtype_fraction) %>%
    pivot_wider(names_from = Tissue, values_from = mean_subtype_fraction, 
                names_prefix = "fraction_") %>%
    filter(!is.na(fraction_TC) & !is.na(fraction_IM))
  
  if(nrow(paired_data_sub) > 0) {
    paired_data_sub$TC_bin <- create_fraction_bins(paired_data_sub$fraction_TC)
    paired_data_sub$IM_bin <- create_fraction_bins(paired_data_sub$fraction_IM)
    
    # Prepare data for alluvial plot
    sankey_data_sub <- paired_data_sub %>%
      count(RFS_group, TC_bin, IM_bin) %>%
      rename(freq = n)
    
    # Create Sankey/Alluvial plot
    p4_sub <- ggplot(sankey_data_sub, aes(axis1 = TC_bin, axis2 = IM_bin, y = freq)) +
      geom_alluvium(aes(fill = RFS_group), alpha = 0.7) +
      geom_stratum() +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_x_discrete(limits = c("TC", "IM"), expand = c(0.1, 0.1)) +
      facet_wrap(~RFS_group) +
      theme_minimal() +
      labs(title = paste(subtype, "Fraction Changes from TC to IM"),
           subtitle = "Flow shows transition patterns by recurrence status",
           x = "Tissue Region",
           y = "Number of Patients",
           fill = "RFS Status") +
      scale_fill_manual(values = c("No Early Relapse" = "#0073c2b2", "Early Relapse" = "#efc000b2"))
    
    print(p4_sub)
    ggsave(file.path(figureDir, paste0(clean_name, "_sankey.pdf")), p4_sub, width = 12, height = 8)
  } else {
    cat(paste("No paired TC-IM data for", subtype, "- skipping Sankey plot\n"))
  }
  
  return(list(
    barplot = p1_sub,
    boxplot_rfs = p2_sub,
    boxplot_treatment = p3_sub,
    sankey = if(exists("p4_sub")) p4_sub else NULL,
    data = subtype_data,
    paired_data = if(exists("paired_data_sub")) paired_data_sub else NULL
  ))
}
