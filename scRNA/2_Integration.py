import os
import gc

import scanpy as sc

import pandas as pd
import numpy as np

## Function to load and process scRNA-seq datasets
from integration_functions import *

#%% Main here
## Set workDir
DataPath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA")

if not os.path.exists(figurePath):
    os.makedirs(figurePath)

## View all data
dataFiles = [x for x in os.listdir(DataPath) if os.path.isdir(os.path.join(DataPath, x))] 

preprocessed_datasets = {}

## Preprocess each dataset
for study_name in dataFiles:
    print(f"Study: {study_name}")

    ## Load and process each dataset
    data_path = os.path.join(DataPath, study_name)
    save_path = os.path.join(DataPath, f"preprocessed_{study_name}.h5ad")

    ## Check if preprocessed data already exists
    if os.path.exists(save_path):
        print(f"Preprocessed data for {study_name} already exists. Skipping preprocessing.")
        
        ## Read and combine preprocessed data
        preprocessed_adata = sc.read(save_path)
        preprocessed_datasets[study_name] = preprocessed_adata
        
        del preprocessed_adata
        gc.collect()
        continue
    
    else:
        ## Read preprocessed data
        adata = load_seurat_data(data_dir = data_path, prefix="processed")

        ## Start Preprocessing
        preprocessed_adata = preprocess_single_dataset(adata.copy(), study_name)
        preprocessed_adata.obs['study'] = study_name  # Add study name to metadata

        ## Combine preprocessed data
        preprocessed_datasets[study_name] = preprocessed_adata
            
        # Print summary
        print(f"\n📊 {study_name} Summary:")
        print(f"  Final shape: {preprocessed_adata.shape}")

        ## Save data
        preprocessed_adata.write(save_path)
        del adata, preprocessed_adata
        gc.collect()

## Perform integration
#%%
merged_adata = merge_datasets_for_integration(preprocessed_datasets)
n_top_genes=2000

# Phase 1: Quick method comparison (20K cells)
quick_results, quick_metrics, quick_fig = quick_method_comparison(
    merged_adata, 
    batch_key='batch',
    n_sample=10000,  # Fast comparison
    n_top_genes=n_top_genes
)

print("🔍 Quick comparison metrics:")
print(quick_metrics)

print(quick_fig)
quick_fig.savefig(
    os.path.join(figurePath, "quick_method_comparison.pdf"), 
    dpi=300, bbox_inches='tight')

# Phase 2: Integration without downsampling
print(f"\n🎯 PHASE 2: Detailed Integration")
print("="*60)

# Prepare data with larger sample
final_integrated = prepare_for_integration(
    merged_adata, 
    n_top_genes=n_top_genes,
    downsample_n=None,
    batch_key='batch'
)

# Modify the class in the final integrated data
# Columns typically representing continuous numerical values
numerical_cols = ['nCount_RNA', 'nFeature_RNA', 'percent.mt', 'complexity', 'Age',]

# Columns typically representing discrete categories or groups
categorical_cols = [
    'study', 'patient', 'tissue', 'sample_id', "doublet", 'patient_tissue',
    'tissueunique',  'Gender', 'TNM', 'AJCC_Stage', 'Tumor_Grade', 'MSI', 'batch'
]

# Assign data classes for obs
final_integrated = assign_data_classes_for_obs(final_integrated,
    numerical_cols=numerical_cols, 
    categorical_cols=categorical_cols
)

sc.external.pp.harmony_integrate(
    final_integrated, 
    key='batch',
    basis='X_pca',
    adjusted_basis='X_pca_harmony'
)

# Save final result
final_integrated.write(f"{DataPath}/final_integrated.h5ad")