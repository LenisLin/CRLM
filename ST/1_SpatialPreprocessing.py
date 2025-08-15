## Script for processing 10x Visium spatial data

import os
import json
import gc
import random
import warnings
import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
import numpy as np

warnings.filterwarnings('ignore')

from Preprocessing_functions import *

## Set up paths
data_dir = '/mnt/NAS_21T/ProjectData/IMC_CRLM/ST'
output_dir = '/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/raw'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#%% GSE225857_10xVisium_ScienceAdvances_2023
sample_dir = os.path.join(data_dir, 'GSE225857_10xVisium_ScienceAdvances_2023')
samples = os.listdir(sample_dir)

dataset_name = 'GSE225857'

## Load the 10x Visium data
for sample in samples:
    sample_path = os.path.join(sample_dir, sample)
    if not os.path.isdir(sample_path):
        continue
    
    print(f'Processing sample: {sample}')
    
    # Load the data
    # Read the MTX format files
    adata = sc.read_10x_mtx(
        os.path.join(sample_path, 'filtered_feature_bc_matrix'),
        var_names='gene_symbols',
        cache=True
    )

    # Read spatial information
    spatial_path = os.path.join(sample_path, 'spatial')

    # Read tissue positions
    positions = pd.read_csv(
        os.path.join(spatial_path, 'tissue_positions_list.csv'), 
        header=None, 
        index_col=0
    )
    positions.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']

    # Add spatial info to adata.obs
    adata.obs = adata.obs.join(positions, how='left')

    # Add spatial coordinates to obsm
    adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values

    # Read scale factors
    with open(os.path.join(spatial_path, 'scalefactors_json.json'), 'r') as f:
        scalefactors = json.load(f)

    # Load H&E images
    # Load high and low resolution images
    hires_img = np.array(Image.open(os.path.join(spatial_path, 'tissue_hires_image.png')))
    lowres_img = np.array(Image.open(os.path.join(spatial_path, 'tissue_lowres_image.png')))

    # Store spatial data in uns (following scanpy convention)
    adata.uns['spatial'] = {
        sample: {
            'scalefactors': scalefactors,
            'images': {
                'hires': hires_img,
                'lowres': lowres_img
            },
            'metadata': {
                'chemistry_description': "Visium",
                'software_version': "spaceranger-1.1.0"
            }
        }
    }

    # Make variable names unique (in case of duplicates)
    adata.var_names_make_unique()

    # Basic information about the loaded data
    print(f"Number of spots: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    print(f"Spatial coordinates available: {'spatial' in adata.obsm}")
    
    # Save the processed data
    adata.write_h5ad(os.path.join(output_dir, f'{dataset_name}_{sample}_raw.h5ad'))
    
    print(f'Sample {sample} processed and saved.')

#%% OEP001756_10xVisium_CancerDiscovery_2022
sample_dir = os.path.join(data_dir, 'OEP001756_10xVisium_CancerDiscovery_2022')
samples = os.listdir(sample_dir)

dataset_name = 'OEP001756'

## Load the 10x Visium data
for sample in samples:
    sample_path = os.path.join(sample_dir, sample)
    if not os.path.isdir(sample_path):
        continue
    
    print(f'Processing sample: {sample}')
    
    # Load the data
    # Read the MTX format files
    adata = sc.read_visium(os.path.join(sample_path))

    # Make variable names unique (in case of duplicates)
    adata.var_names_make_unique()

    # Basic information about the loaded data
    print(f"Number of spots: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    print(f"Spatial coordinates available: {'spatial' in adata.obsm}")
    
    # Save the processed data
    adata.write_h5ad(os.path.join(output_dir, f'{dataset_name}_{sample}_raw.h5ad'))
    
    print(f'Sample {sample} processed and saved.')


#%% FDZS Stereo-seq
# Stereo-seq data processing
import stereo as st
from stereo.plots import violin_distribution

import scanpy as sc
import squidpy as sq

# Load Stereo-seq data
sample_dir = os.path.join(data_dir, 'FDZS_stereoseq') 
dataset_name = 'FDZS'
sample = 'A04932G3_bin100'

## Start from the processed data
# # data = st.io.read_h5ad(os.path.join(sample_dir,A04932G3.cellbin_1.0.adjusted.h5ad")) ## Cellbin data
# data = st.io.read_h5ad(os.path.join(sample_dir,"A04932G3.bin50_1.0.h5ad")) ## Bin50 data
# adata = st.io.stereo_to_anndata(data)

# # # Find the correct attribute
# # all_attrs = [attr for attr in dir(data) if not attr.startswith('_')]
# # all_attrs

# adata = combine_stereo_process_info(adata, data)

## Start from the raw data
data = st.io.read_gef(file_path=os.path.join(sample_dir,"A04932G3.tissue.gef"), bin_size=100)
adata = st.io.stereo_to_anndata(data, output=os.path.join(sample_dir,'A04932G3_anndata.h5ad'))

## Adjust Gene Names
real_gene_name = [str(x) for x in adata.var['real_gene_name']]

adata.var_names = real_gene_name
adata.var_names_make_unique()

## Identify the TMA
identify_tma_cores(adata, eps=0.1, min_samples=50, plot=True,saveDir=sample_dir)

## Add the tme order 
core_mapping = {
            '1': 'C6', '2': 'C5', '3': 'C7', '4': 'A5', '5': 'A6', 
            '6': 'A7', '7': 'A1', '8': 'A2', '9': 'A3'}

# Create new core names column
adata.obs['core_name'] = adata.obs['tma_core'].map(core_mapping)

# Display new core distribution
current_counts = adata.obs['core_name'].value_counts().sort_index()
for core, count in current_counts.items():
    print(f"  Core {core}: {count:,} cells")

adata.obs['core_name'] = adata.obs['core_name'].astype('category') # Convert to categorical
adata.obs['RFS_status'] = adata.obs['core_name'].str.startswith('A').astype(int)

## Save the processed data
adata.write_h5ad(os.path.join(output_dir, f'{dataset_name}_{sample}_raw.h5ad'))