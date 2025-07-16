## functions for preprocessing data
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(arrow)

# Function to add clinical information from a dataframe to a SpatialExperiment object
add_clinical_info <- function(spe_, df_, add_columns_, match_col_1, match_col_2) {
  
  # Extract metadata in spe
  spe_meta <- colData(spe_)
  
  # Ensure the matching columns and additional columns exist in respective dataframes
  if (!match_col_1 %in% colnames(spe_meta)) stop("Index column not found in expression profile object")
  if (!match_col_2 %in% colnames(df_)) stop("Index column not found in clinical data")
  if (!all(add_columns_ %in% colnames(df_))) stop("Some add features not found in clinical data")
  
  # Subset df_ to include only the necessary columns for merging
  cat("Add following features:",add_columns_,"\n")
  df_subset <- df_[, c(match_col_2, add_columns_), drop = FALSE]
  
  # Merge the clinical information with the SpatialExperiment object
  for(feature_ in add_columns_){
    colData(spe_)[feature_] <- df_subset[match(spe_meta[,match_col_1],df_subset[,match_col_2]),
                                         match(feature_,colnames(df_subset))]
  }
  
  return(spe_)
}

# Function to normalize
normData <- function(spe_, censor_val = NULL, arcsinh = FALSE, norm_method = "0-1", savePath = NULL) {
  #' censor_val=0.999'
  exp_data <- data.frame(t(assay(spe_, "exprs")), check.names = TRUE)
  
  if (arcsinh) {
    exp_data <- asinh(exp_data)
  }
  

  exp_data$sample_id <- spe_$sample_id
  
  
    dat <- data.table(exp_data %>% pivot_longer(-sample_id, names_to = "channel", values_to = "mc_counts"))
    dat$mc_counts <- as.numeric(dat$mc_counts)
    
    if (!is.null(censor_val)) {
      dat[, c_counts := censor_dat(mc_counts, censor_val), by = channel]
      if (norm_method == "0-1") {
        dat[, c_counts_scaled := ((c_counts - min(c_counts)) / (max(c_counts) - min(c_counts))), by = channel]
        dat[c_counts_scaled < 0, c_counts_scaled := 0, by = channel]
      } else if (norm_method == "znorm") {
        dat[, c_counts_scaled := ((c_counts - mean(c_counts)) / sd(c_counts)), by = channel]
      } else if (norm_method == "null") {
        dat[, c_counts_scaled := c_counts, by = channel]
      }
    } else {
      if (norm_method == "0-1") {
        dat[, c_counts_scaled := ((mc_counts - min(mc_counts)) / (max(mc_counts) - min(mc_counts))), by = channel]
        dat[c_counts_scaled < 0, c_counts_scaled := 0, by = channel]
      } else if (norm_method == "znorm") {
        dat[, c_counts_scaled := ((mc_counts - mean(mc_counts)) / sd(mc_counts)), by = channel]
      }
    }
  
  
  a <- data.frame(dat$channel, dat$c_counts_scaled) %>%
    group_by(dat.channel) %>%
    dplyr::mutate(index = row_number()) %>%
    pivot_wider(
      names_from = dat.channel,
      values_from = dat.c_counts_scaled
    ) %>%
    dplyr::select(-index)
  
  histdata <- reshape2::melt(a, variable.name = "marker", value.name = "expression")
  p1 <- ggplot(data = histdata, aes(x = expression)) +
    geom_histogram(bins = 60, colour = "black", fill = "blue", alpha = 0.5) +
    facet_wrap(~marker, scale = "free")
  
  ggsave(filename = file.path(savePath,"Marker Expression level.pdf"), height = 8, width = 12)
  
  ## re-assign
  assay(spe_, "raw") <- assay(spe_, "exprs")
  norm_data <- t(a)
  colnames(norm_data) <- colnames(assay(spe_, "exprs"))
  assay(spe_, "exprs") <- norm_data
  return(spe_)
}

censor_dat <- function(x, quant = 0.999, ignorezero = TRUE, symmetric = F) {
  if (symmetric) {
    lower_quant <- (1 - quant) / 2
    quant <- quant + lower_quant
  }
  
  if (ignorezero) {
    q <- stats::quantile(x[x > 0], quant)
  } else {
    q <- stats::quantile(x, quant)
  }
  x[x > q] <- q
  
  if (symmetric) {
    q <- stats::quantile(x, lower_quant)
    x[x < q] <- q
  }
  return(x)
}

## calculate the fraction of cell types in each image
Transform_CellCountMat <- function(spe_, clinicalFeatures, img_id = "roi", count_by = "sub_celltype", is.fraction = TRUE) {
  
  # Extract cell metadata
  cellMeta <- colData(spe_)
  
  # Get unique ROI, subtypes, and other types
  ROIs <- unique(cellMeta[,img_id])
  SubTypes <- unique(cellMeta[[count_by]])
  alltypes <- unique(c(SubTypes))
  
  # Initialize the cell count matrix
  CellCountMat <- matrix(0, nrow = length(ROIs), ncol = length(alltypes))
  CellCountMat <- as.data.frame(CellCountMat)
  rownames(CellCountMat) <- ROIs
  colnames(CellCountMat) <- alltypes
  
  # Loop through each ROI and count cells
  for (ROI in ROIs) {
    coldataTemp <- cellMeta[as.character(cellMeta[,img_id]) %in% ROI, ]
    cellnum <- nrow(coldataTemp)
    
    # Count cells by subtype
    SubTem <- as.data.frame(table(coldataTemp[[count_by]]))
    
    # Update the cell count matrix
    for (i in 1:nrow(SubTem)) {
      subtype <- SubTem$Var1[i]
      subtype <- as.character(subtype)
      count <- SubTem$Freq[i]
      if (is.fraction) {
        CellCountMat[ROI, subtype] <- count / cellnum
      } else {
        CellCountMat[ROI, subtype] <- count
      }
    }
  }
  
  # Add clinical features to the cell count matrix
  CellCountMat$PID <- sapply(rownames(CellCountMat), function(x) strsplit(x, "_")[[1]][1])
  
  for (feature in clinicalFeatures) {
    CellCountMat[,feature] <- cellMeta[match(rownames(CellCountMat), cellMeta[,img_id]), feature]
  }
  
  return(CellCountMat)
}

# Function to split and extract SampleID and SampleID_ROIID  
extract_sample_info <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  
  ## Condition
  len <- length(parts)
  if(len==1){
    cat("Sample_id was error.","\n")
    return(NA)
  }
  if(len==2){
    sample_id <- paste0(parts,collapse = "_")
  }
  if(len>=3){
    sample_id <- paste0(parts[-1],collapse = "_")
  }
  
  ## return
  return(sample_id)  
}  

## save the expression matrix and meta data of spe
save_exp_meta <- function(spe, saveDir){
  exp_ <- as.data.frame(assay(spe))
  meta_ <- as.data.frame(colData(spe))  
  
  ## expression profile
  exp_ <- as_arrow_table(exp_)
  write_feather(x = exp_, sink = file.path(saveDir,"unanno_exp.feather"),version = 1)
  
  ## meta information
  meta_ <- as_arrow_table(meta_)
  write_feather(x = meta_, sink = file.path(saveDir,"unanno_meta.feather"))
  
  ## gene name
  write.table(x = rowData(spe), file = file.path(saveDir,"unanno_genename.csv"),sep = ',')
  
  ## cell name
  write.table(x = colnames(spe), file = file.path(saveDir,"unanno_cellname.csv"),sep = ',')
  
  rm(exp_,meta_)
  return(NULL)
}

