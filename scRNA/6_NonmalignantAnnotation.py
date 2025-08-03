#%% Annotate the non-malignant cells in the scRNA dataset
import os
import gc
import pickle
from re import T
import select
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

from annotation_functions import *

print("🧬 Starting non-malignant Cell Annotation")
print("="*70)

# Set paths
DataPath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","nonMalignant_Analysis")
result_path = figurePath

# create directories if they don't exist
os.makedirs(figurePath, exist_ok=True)

#%% Task 1: Extract Non-Malignant Cells with Original Expression
print("\n🎯 Task 1: Extracting Non-Malignant Cells")
print("="*60)

# Load integrated data if not already loaded
if 'final_integrated' not in locals():
    final_integrated = sc.read_h5ad(f"{DataPath}/major_anno_all.h5ad")

# Check major cell types
print("Major cell type distribution:")
print(final_integrated.obs['Major_type'].value_counts())
print(final_integrated.obs['tissue'].value_counts())

# Load the annotated malignant cells
if 'malignant_cells' not in locals():
    result_path_ = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","Malignant_Analysis")
    malignant_cells = sc.read_h5ad(os.path.join(result_path_, "malignant_epithelial_cells_annotated.h5ad"))

# Initialize the Sub_type column
final_integrated.obs["cell_id"] = final_integrated.obs.index
malignant_cells.obs["cell_id"] = malignant_cells.obs.index

final_integrated.obs["Sub_type"] = final_integrated.obs["Major_type"].copy()
mapping = malignant_cells.obs[['cell_id', 'Malignant_type']].set_index('cell_id')['Malignant_type'].to_dict()
matched_cells = final_integrated.obs['cell_id'].isin(mapping.keys())
final_integrated.obs.loc[matched_cells, 'Sub_type'] = final_integrated.obs.loc[matched_cells, 'cell_id'].map(mapping)
final_integrated.obs['Sub_type'] = final_integrated.obs['Sub_type'].astype(str)

# Filter for non-maligant cells only
nonmalignant_mask = (final_integrated.obs['Sub_type'] == 'Unknown')
adata = final_integrated[nonmalignant_mask].copy()
adata.write_h5ad(os.path.join(result_path, "nonmalignant_cells_raw.h5ad"))

del final_integrated, malignant_cells, mapping, nonmalignant_mask
gc.collect()

#%% Task 2: Downsample and Visualize Batch Effects
# Downsample for visualization (stratified by batch)
target_cells = 20000
sampled_indices = np.random.choice(adata.shape[0], size=target_cells, replace=False)
downsampled_adata = adata[sampled_indices].copy()

# Basic preprocessing for visualization
sc.pp.normalize_total(downsampled_adata, target_sum=1e4)
sc.pp.log1p(downsampled_adata)
sc.pp.highly_variable_genes(downsampled_adata, n_top_genes=2000)
sc.pp.scale(downsampled_adata)
sc.tl.pca(downsampled_adata)
sc.pp.neighbors(downsampled_adata, n_neighbors=15)
sc.tl.umap(downsampled_adata)

# Create batch visualization
fig, axes = plt.subplots(2, 2, figsize=(18, 12))

# Batch effect
sc.pl.umap(downsampled_adata, color='Major_type', ax=axes[0,0], show=False, 
          legend_loc='right margin', title='Major Cell Types')
sc.pl.umap(downsampled_adata, color='batch', ax=axes[1,0], show=False,
          legend_loc='right margin', title='Batch Effect')
sc.pl.umap(downsampled_adata, color='tissue', ax=axes[0,1], show=False,
          legend_loc='right margin', title='Tissue Type')
sc.pl.umap(downsampled_adata, color='percent.mt', ax=axes[1,1], show=False,
          title='MT Percentage', color_map='Reds')

plt.tight_layout()
plt.savefig(os.path.join(figurePath, "All nonmalignant cell UMAP.pdf"), dpi=300, bbox_inches='tight')
plt.show()

#%% Task 3: Clustering for each subtype
print("\n🎯 Task 3: Annotating Non-Malignant Cells")
print("="*60)

# Load the Major annotation data
if 'adata' not in locals():
    adata = sc.read_h5ad(os.path.join(result_path, "nonmalignant_cells_raw.h5ad"))

# Check the distribution of Major types
print("Major cell type distribution in non-malignant cells:")
print(adata.obs['Major_type'].value_counts())

marker_genes_dict = {
    'T': {
        # Pan T cell
        'Pan_T':["CD3E","CD8A","CD8B","CD4"],

        # CD4+ T cell subtypes
        'CD4_Naive': ['SELL', 'CCR7', 'TCF7'],
        'CD4_Memory_CCL5': ['CCL5', 'CCL4', 'IFNG'],  # Effector memory
        'CD4_Exhausted': ['TOX2', 'TIGIT', 'PDCD1'],  # Exhausted CD4+
        'CD4_Th17': ['KLRB1', 'CCR6', 'IL17A'],  # Th17 cells
        'CD4_LTB': ['LTB'],  # LTB+ CD4+ T cells
        'CD4_HSP': ['HSPA1A', 'HSP90AA1'],  # Heat shock response
        'CD4_CXCL13': ['CXCL13'],  # Follicular helper-like
        'Treg': ['FOXP3', 'TIGIT', 'CTLA4', 'IL2RA'],  # Regulatory T cells
        
        # CD8+ T cell subtypes
        'CD8_Naive': ['SELL', 'CCR7', 'TCF7'],  # Naive CD8+
        'CD8_GZMK': ['GZMK'],  # Memory CD8+
        'CD8_CD52': ['CD52'],  # CD52+ CD8+
        'CD8_Exhausted': ['CTLA4', 'TNFRSF9', 'HAVCR2', 'TIGIT'],  # Exhausted
        'CD8_GIMAP7': ['GIMAP7', 'GIMAP4', 'GIMAP1'],  # GIMAP7+ CD8+
        'CD8_Effector': ['FGFBP2', 'GZMB', 'PRF1'],  # Effector CD8+
        'CD8_XCL1': ['XCL1'],  # XCL1+ CD8+
        'CD8_HSP': ['HSPA6', 'HSPA1A', 'HSPD1'],  # Heat shock CD8+
        'CD8_CXCL13': ['CXCL13', 'PDCD1'],  # Tissue-resident memory
        
        # MAIT cells
        'MAIT': ['SLC4A10', 'DPP4', 'KLRB1'],
        
        # γδ T cells
        'GammaDelta_T': ['TRGC1', 'TRGC2', 'TRDV1']
    },
    
    'Myeloid': {
        # Myeloid cell subtypes
        'Macrophage': ['CD68','CD14'],  # General macrophage marker

        # Macrophage subtypes
        'Mac_M1': ['CXCL10', 'TNF', 'IL1B', 'CCL3'],  # M1-like pro-inflammatory
        'Mac_M2': ['MRC1', 'CCL18', 'MARCO'],  # M2-like anti-inflammatory
        'Mac_SPP1': ['SPP1', 'MARCO', 'LIPA'],  # SPP1+ tumor-associated
        'Mac_Proliferating': ['MKI67'],  # Proliferating macrophages
        
        # Monocyte subtypes
        'Mono_Classical': ['CD14', 'S100A8', 'S100A9', 'S100A12'],  # Classical monocytes
        'Mono_HSP': ['HSPA1B', 'HSPH1', 'HSPA1A', 'HSPB1'],  # Heat shock monocytes
        
        # Dendritic cell subtypes
        'cDC1': ['CLEC9A', 'XCR1', 'BATF3'],  # Conventional DC1
        'cDC2': ['CLEC10A', 'CD1C', 'FCER1A'],  # Conventional DC2
        'pDC': ['LILRA4', 'IRF7', 'CLEC4C'],  # Plasmacytoid DC
        
        # Neutrophils
        'Neutrophil': ['FCGR3B', 'CXCR2', 'CEACAM8']
    },
    
    'NK': {
        # NK cell subtypes
        'NK_GZMK': ['GZMK', 'XCL1'],  # GZMK+ NK cells
        'NK_GZMK_IL7R': ['GZMK', 'IL7R', 'GPR183'],  # GZMK+ IL7R+ NK
        'NK_GZMB': ['GZMB', 'GNLY'],  # GZMB+ cytotoxic NK
        'NK_GZMB_MYOM2': ['GZMB', 'MYOM2'],  # GZMB+ MYOM2+ NK
        'NK_Adaptive': ['FCGR3A', 'CD57']  # Adaptive NK cells
    },
    
    'Fibroblast': {
        # Fibroblast subtypes
        'Fib_F2_MCAM': ['MCAM', 'JAG1', 'NOTCH3'],  # F2_MCAM (LM-enriched)
        'Fib_F4_F3': ['F3', 'NRG1'],  # F4_F3 (CC-enriched)
        'Fib_F5_CCL11': ['CCL11'],  # F5_CCL11 (CC-enriched)
        'CAF_FAP': ['FAP', 'PDPN'],  # FAP+ cancer-associated fibroblasts
        'CAF_SMA': ['ACTA2'],  # α-SMA+ myofibroblasts
        'Fib_Matrix': ['COL3A1', 'COL4A1', 'LUM', 'VCAN'],  # Matrix-producing
        'Fib_Proliferating': ['MKI67']  # Proliferating fibroblasts
    },
    
    'B_Plasma': {
        # B cell subtypes
        'Pan_B': ['CD19', 'MS4A1'],  # General B cell markers
        'B_Naive': ['TCL1A', 'SELL'],  # Naive B cells
        'B_Memory': ['AIM2'],  # Memory B cells
        'B_Activated': ['EGR1'],  # Activated B cells
        'B_Proliferating': ['MKI67'],  # Proliferating B cells
        'B_HSP': ['HSPA1A', 'HSP90AA1', 'HSPE1'],  # Heat shock B cells
        'B_GC': ['BCL6', 'AICDA'],  # Germinal center B cells

        # Plasma cell subtypes
        'Plasma_IgA': ['IGHA1', 'IGHA2'],  # IgA+ plasma cells
        'Plasma_IgG': ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4'],  # IgG+ plasma cells
        'Plasma_IgM': ['IGHM'],  # IgM+ plasma cells
        'Plasma_Proliferating': ['MKI67']  # Proliferating plasma cells
    },
    
    'Endothelial': {
        # Pan-endothelial markers
        'Endo_Pan': ['PECAM1', 'CD34', 'VWF'],  # General endothelial markers
        # Endothelial subtypes
        'Endo_Arterial': ['EFNB2', 'HEY1'],  # Arterial endothelial
        'ADD':["SELP","DLL4","NOTCH3","CD36","CLEC4G"],
        'Endo_Venous': ['NR2F2', 'ACKR1'],  # Venous endothelial
        'Endo_Capillary': ['CA4', 'RGCC'],  # Capillary endothelial
        'Endo_Lymphatic': ['LYVE1', 'PROX1', 'FLT4'],  # Lymphatic endothelial
        'Endo_Angiogenic': ['VEGFA', 'ANGPT2', 'ESM1'],  # Angiogenic endothelial
        'Endo_Proliferating': ['MKI67']  # Proliferating endothelial
    }
}

# Cell types to annotate (excluding Epithelial, Mast, Hepatocyte as requested)
cell_types_to_annotate = ['T'] # ,'Myeloid', 'NK',  'Fibroblast', 'B_Plasma', 'Endothelial'
# 
# Different resolutions for different cell types (adjust as needed)
resolutions = {
    'T': 1.0,           # T cells are diverse, need higher resolution
    'Myeloid': 0.6,     # Myeloid cells have clear subtypes
    'NK': 0.3,          # NK cells are more homogeneous
    'Fibroblast': 0.6,  # Moderate diversity
    'B_Plasma': 0.6,           # Moderate diversity  
    'Endothelial': 0.6  # Less diverse
}

# Annotate each cell type
for cell_type in cell_types_to_annotate:
    if cell_type in marker_genes_dict:
        preprocessing_for_subtype_anno(
            adata.copy(), 
            cell_type, 
            marker_genes = marker_genes_dict[cell_type], 
            figurePath = figurePath,
            resolution=resolutions.get(cell_type, 0.5)
        )

#%% Task 4: Annotation T cells

## Load Anndata of Subtype
cell_type = "T"
sub_anndata = sc.read_h5ad(os.path.join(figurePath, f"{cell_type}_analysis", f"{cell_type}_processed.h5ad"))

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))
sc.pl.dotplot(sub_anndata, var_names=marker_genes_dict[cell_type], groupby="leiden",show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1)
plt.xticks(rotation=45, fontsize=10) # Add styling after the plot is created
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type}_specific_marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()

new_cluster_names = [
    'CD8T_GZMK','CD4T_','MAIT_SLC4A10','_IL7R','CD8T_GZMB',
    'Treg_FOXP3','CD8T_GZMK','CD8T_MKI67','CD4T_CXCL13','9',
    '10','Unknown','NaiveT_TCF7'
]

## Map cluster names to cell types
cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}
sub_anndata.obs['Sub_type'] = sub_anndata.obs['leiden'].map(cluster_to_celltype)

## Mapping to major anndata
adata.obs.loc[:,"Major_type"] = adata.obs['Major_type'].astype(str)  # Ensure Sub_type is string type
adata.obs.loc[:,"Sub_type"] = adata.obs['Sub_type'].astype(str)  # Ensure Sub_type is string type

idx = sub_anndata.obs.index
adata.obs.loc[idx,"Sub_type"] = sub_anndata.obs['Sub_type']

## Mannual adjustment for Sub_type
adata.obs["Sub_type"].value_counts()
adata.obs.loc[adata.obs["Sub_type"] == "Unknown","Major_type"] = "Unknown"
sub_anndata = sub_anndata[sub_anndata.obs["Sub_type"].str.startswith(("B", "Plasma"))]

# Create the UMAP plot with enhanced styling
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
sc.pl.umap(
    sub_anndata, color="Sub_type",
    ax=ax, show=False, frameon=True, size=1.5, alpha=0.8,
    #palette=malignant_subtype_colors,
    legend_loc='right margin', legend_fontsize=7, legend_fontweight='normal'
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')
plt.tight_layout()

# Save with high quality settings for publication
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type} subtypes UMAP.pdf'), 
    dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf'
)
plt.show()
plt.close()

# Create a more appealing dot plot
selected_marker_genes = {
    # B cell subtypes
    'B_MS4A1': ['CD19', 'MS4A1','SELL'],  # General B cell markers

    # Plasma cell subtypes
    'Plasma_IgA': ['IGHA1', 'IGHA2'],  # IgA+ plasma cells
    'Plasma_IgG': ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4'],  # IgG+ plasma cells
    'Plasma_CD59': ['CD59',"IGLL5"],  # CD59+ plasma cells
    'Plasma_Proliferating': ['MKI67']  # Proliferating plasma cells
}

plt.figure(figsize=(15, 8))
sc.pl.dotplot(
    sub_anndata, var_names=selected_marker_genes, groupby="Sub_type",
    show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath,f"{cell_type}_analysis", f'{cell_type} select marker genes dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()
plt.close()

#%% Task 5: Annotation Myeloid cells

## Load Anndata of Subtype
cell_type = "Myeloid"
sub_anndata = sc.read_h5ad(os.path.join(figurePath, f"{cell_type}_analysis", f"{cell_type}_processed.h5ad"))

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))
sc.pl.dotplot(sub_anndata, var_names=marker_genes_dict[cell_type], groupby="leiden",show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1)
plt.xticks(rotation=45, fontsize=10) # Add styling after the plot is created
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type}_specific_marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()

new_cluster_names = [
    'Macro_VCAN', 'Macro_APOE', 'cDC_CD1C', 'Macro_APOE', 'Neutrophil', 
    'Macro_SPP1', 'cDC_CLEC9A', 'Macro_S100A12', 'Macro_FOLR2', 'pDC_LILRA4',
    'Macro_HSPH1', 'Macro_APOE', 'Macro_MKI67', 'Macro_APOE'
]

## Map cluster names to cell types
cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}
sub_anndata.obs['Sub_type'] = sub_anndata.obs['leiden'].map(cluster_to_celltype)

## Mapping to major anndata
adata.obs.loc[:,"Major_type"] = adata.obs['Major_type'].astype(str)  # Ensure Sub_type is string type
adata.obs.loc[:,"Sub_type"] = adata.obs['Sub_type'].astype(str)  # Ensure Sub_type is string type

idx = sub_anndata.obs.index
adata.obs.loc[idx,"Sub_type"] = sub_anndata.obs['Sub_type']

# Create the UMAP plot with enhanced styling
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
sc.pl.umap(
    sub_anndata, color="Sub_type",
    ax=ax, show=False, frameon=True, size=1.5, alpha=0.8,
    #palette=malignant_subtype_colors,
    legend_loc='right margin', legend_fontsize=7, legend_fontweight='normal'
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')
plt.tight_layout()

# Save with high quality settings for publication
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type} subtypes UMAP.pdf'), 
    dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf'
)
plt.show()
plt.close()

# Create a more appealing dot plot
selected_marker_genes = {
    "Myeloid": ['CD68','CD14'],

    # Macrophage subtypes
    'Macro_SPP1': ['SPP1', 'MARCO'],  # SPP1+ tumor-associated
    'Macro_MKI67': ['MKI67',"STMN1"],  # Proliferating macrophages
    'Macro_VCAN': ["VCAN"],
    'Macro_APOE': ["APOE"],
    'Macro_S100A12': ["S100A12"],
    'Macro_FOLR2': ["FOLR2"],
    'Macro_HSPH1': ['HSPA1B', 'HSPA1A', 'HSPB1'],  # Heat shock macrophages

    # Dendritic cell subtypes
    'cDC_CLEC9A': ['CLEC9A', 'XCR1', 'BATF3'],  # Conventional DC1
    'cDC_CD1C': ['CLEC10A', 'CD1C', 'FCER1A'],  # Conventional DC2
    'pDC_LILRA4': ['LILRA4', 'CLEC4C'],  # Plasmacytoid DC

    # Neutrophils
    'Neutrophil': ['FCGR3B', 'CXCR2']
}

plt.figure(figsize=(15, 8))
sc.pl.dotplot(
    sub_anndata, var_names=selected_marker_genes, groupby="Sub_type",
    show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath,f"{cell_type}_analysis", f'{cell_type} select marker genes dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()
plt.close()

#%% Task 6: Annotation NK cells

## Load Anndata of Subtype
cell_type = "NK"
sub_anndata = sc.read_h5ad(os.path.join(figurePath, f"{cell_type}_analysis", f"{cell_type}_processed.h5ad"))

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))
sc.pl.dotplot(sub_anndata, var_names=marker_genes_dict[cell_type], groupby="leiden",show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1)
plt.xticks(rotation=45, fontsize=10) # Add styling after the plot is created
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type}_specific_marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()

new_cluster_names = ['NK_GZMK', 'NK_FCGR3A', 'Unknown', 'Unknown']

## Map cluster names to cell types
cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}
sub_anndata.obs['Sub_type'] = sub_anndata.obs['leiden'].map(cluster_to_celltype)

## Mapping to major anndata
idx = sub_anndata.obs.index
adata.obs.loc[idx,"Sub_type"] = sub_anndata.obs['Sub_type'].tolist()

## Mannual adjustment for Sub_type
adata.obs["Sub_type"].value_counts()
adata.obs.loc[adata.obs["Sub_type"] == "Unknown","Major_type"] = "Unknown"
sub_anndata = sub_anndata[sub_anndata.obs["Sub_type"].str.startswith("NK")] ## Subset

# Create the UMAP plot with enhanced styling
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
sc.pl.umap(
    sub_anndata, color="Sub_type",
    ax=ax, show=False, frameon=True, size=1.5, alpha=0.8,
    #palette=malignant_subtype_colors,
    legend_loc='right margin', legend_fontsize=7, legend_fontweight='normal'
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')
plt.tight_layout()

# Save with high quality settings for publication
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type} subtypes UMAP.pdf'), 
    dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf'
)
plt.show()
plt.close()

# Create a more appealing dot plot
selected_marker_genes = {
    'NK': ['NKG7',"KLRD1"],
    'NK_FCGR3A': ["FCGR3A","GZMB","GNLY"],  # FCGR3A+ NK cells
    'NK_GZMK': ["GZMK","XCL1","CD160"]  # GZMK+ NK cells
}

plt.figure(figsize=(15, 8))
sc.pl.dotplot(
    sub_anndata, var_names=selected_marker_genes, groupby="Sub_type",
    show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath,f"{cell_type}_analysis", f'{cell_type} select marker genes dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()
plt.close()

#%% Task 7: Annotation Fibroblast

## Load Anndata of Subtype
cell_type = "Fibroblast"
sub_anndata = sc.read_h5ad(os.path.join(figurePath, f"{cell_type}_analysis", f"{cell_type}_processed.h5ad"))

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))
sc.pl.dotplot(sub_anndata, var_names=marker_genes_dict[cell_type], groupby="leiden",show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1)
plt.xticks(rotation=45, fontsize=10) # Add styling after the plot is created
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type}_specific_marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()

new_cluster_names = [
    'CAF_POSTN', 'Fibro_RGS5', 'CAF_POSTN', 'Fibro_RGS5', 'CAF_PDPN',
    'Fibro_COL3A1', 'Fibro_RGS5', 'CAF_POSTN', 'CAF_CXCL12', 'Fibro_MKI67',
    'CAF_CXCL12'
]

## Map cluster names to cell types
cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}
sub_anndata.obs['Sub_type'] = sub_anndata.obs['leiden'].map(cluster_to_celltype)

## Mapping to major anndata
adata.obs.loc[:,"Major_type"] = adata.obs['Major_type'].astype(str)  # Ensure Sub_type is string type
adata.obs.loc[:,"Sub_type"] = adata.obs['Sub_type'].astype(str)  # Ensure Sub_type is string type

idx = sub_anndata.obs.index
adata.obs.loc[idx,"Sub_type"] = sub_anndata.obs['Sub_type'].tolist()

# ## Mannual adjustment for Sub_type
# adata.obs["Sub_type"].value_counts()
# adata.obs.loc[adata.obs["Sub_type"] == "Mast","Major_type"] = "Mast"
# adata.obs.loc[adata.obs["Sub_type"] == "T","Major_type"] = "T"
# sub_anndata = sub_anndata[sub_anndata.obs["Sub_type"].str.startswith("NK")] ## Subset

# Create the UMAP plot with enhanced styling
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
sc.pl.umap(
    sub_anndata, color="Sub_type",
    ax=ax, show=False, frameon=True, size=1.5, alpha=0.8,
    #palette=malignant_subtype_colors,
    legend_loc='right margin', legend_fontsize=7, legend_fontweight='normal'
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')
plt.tight_layout()

# Save with high quality settings for publication
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type} subtypes UMAP.pdf'), 
    dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf'
)
plt.show()
plt.close()

# Create a more appealing dot plot
selected_marker_genes = {
    'Stromal': ['LUM', "DCN",'FAP', 'ACTA2'],
    
    # Fibroblast subtypes
    'Fibro_RGS5': ['RGS5', 'MCAM', 'NOTCH3'],  # F2_MCAM (LM-enriched)
    'Fibro_COL3A1': ['COL3A1',"COL1A2"],  # F4_F3 (CC-enriched)
    'Fibro_MKI67': ['MKI67'],  # F5_CCL11 (CC-enriched)
    'CAF_POSTN': ['POSTN'],
    'CAF_PDPN': ['PDPN',"C3"],
    'CAF_CXCL12': ['CXCL12',"CCL11"],
}


plt.figure(figsize=(15, 8))
sc.pl.dotplot(
    sub_anndata, var_names=selected_marker_genes, groupby="Sub_type",
    show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath,f"{cell_type}_analysis", f'{cell_type} select marker genes dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()
plt.close()

#%% Task 8: Annotation B and Plasma cells

## Load Anndata of Subtype
cell_type = "B_Plasma"
sub_anndata = sc.read_h5ad(os.path.join(figurePath, f"{cell_type}_analysis", f"{cell_type}_processed.h5ad"))

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))
sc.pl.dotplot(sub_anndata, var_names=marker_genes_dict[cell_type], groupby="leiden",show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1)
plt.xticks(rotation=45, fontsize=10) # Add styling after the plot is created
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type}_specific_marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()

new_cluster_names = [
    'B_MS4A1', 'Plasma_IGHG', 'B_CD37', 'B_MS4A1', 'Plasma_IGHA',
    'Plasma_MKI67', 'Plasma_CD59', 'Unknown', 'B_CD37', 'Unknown',
    'Unknown', 'B_MS4A1', 'Plasma_IGHA', 'B_MS4A1'
]

## Map cluster names to cell types
cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}
sub_anndata.obs['Sub_type'] = sub_anndata.obs['leiden'].map(cluster_to_celltype)

## Mapping to major anndata
adata.obs.loc[:,"Major_type"] = adata.obs['Major_type'].astype(str)  # Ensure Sub_type is string type
adata.obs.loc[:,"Sub_type"] = adata.obs['Sub_type'].astype(str)  # Ensure Sub_type is string type

idx = sub_anndata.obs.index
adata.obs.loc[idx,"Sub_type"] = sub_anndata.obs['Sub_type']

## Mannual adjustment for Sub_type
adata.obs["Sub_type"].value_counts()
adata.obs.loc[adata.obs["Sub_type"] == "Unknown","Major_type"] = "Unknown"
sub_anndata = sub_anndata[sub_anndata.obs["Sub_type"].str.startswith(("B", "Plasma"))]

# Create the UMAP plot with enhanced styling
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
sc.pl.umap(
    sub_anndata, color="Sub_type",
    ax=ax, show=False, frameon=True, size=1.5, alpha=0.8,
    #palette=malignant_subtype_colors,
    legend_loc='right margin', legend_fontsize=7, legend_fontweight='normal'
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')
plt.tight_layout()

# Save with high quality settings for publication
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type} subtypes UMAP.pdf'), 
    dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf'
)
plt.show()
plt.close()

# Create a more appealing dot plot
selected_marker_genes = {
    # B cell subtypes
    'B_MS4A1': ['CD19', 'MS4A1','SELL'],  # General B cell markers

    # Plasma cell subtypes
    'Plasma_IgA': ['IGHA1', 'IGHA2'],  # IgA+ plasma cells
    'Plasma_IgG': ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4'],  # IgG+ plasma cells
    'Plasma_CD59': ['CD59',"IGLL5"],  # CD59+ plasma cells
    'Plasma_Proliferating': ['MKI67']  # Proliferating plasma cells
}

plt.figure(figsize=(15, 8))
sc.pl.dotplot(
    sub_anndata, var_names=selected_marker_genes, groupby="Sub_type",
    show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath,f"{cell_type}_analysis", f'{cell_type} select marker genes dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()
plt.close()

#%% Task 9: Annotation Endothelial cells

## Load Anndata of Subtype
cell_type = "Endothelial"
sub_anndata = sc.read_h5ad(os.path.join(figurePath, f"{cell_type}_analysis", f"{cell_type}_processed.h5ad"))

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))
sc.pl.dotplot(sub_anndata, var_names=marker_genes_dict[cell_type], groupby="leiden",show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1)
plt.xticks(rotation=45, fontsize=10) # Add styling after the plot is created
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type}_specific_marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()

new_cluster_names = [
    'EC_EFNB2', 'EC_ACKR1', 'EC_EFNB2', 'EC_ANGPT2', 'EC_LYVE1',
    'EC_ACKR1', 'EC_EFNB2', 'Myeloid', 'EC_EFNB2',
    'EC_MKI67'
]

## Map cluster names to cell types
adata.obs.loc[:,"Major_type"] = adata.obs['Major_type'].astype(str)  # Ensure Sub_type is string type
adata.obs.loc[:,"Sub_type"] = adata.obs['Sub_type'].astype(str)  # Ensure Sub_type is string type

cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}
sub_anndata.obs['Sub_type'] = sub_anndata.obs['leiden'].map(cluster_to_celltype)

## Mapping to major anndata
idx = sub_anndata.obs.index
adata.obs.loc[idx,"Sub_type"] = sub_anndata.obs['Sub_type']

## Mannual adjustment for Sub_type
adata.obs["Sub_type"].value_counts()
adata.obs.loc[adata.obs["Sub_type"] == "Myeloid","Major_type"] = "Myeloid"
sub_anndata = sub_anndata[sub_anndata.obs["Sub_type"].str.startswith("EC")] ## Subset

# Create the UMAP plot with enhanced styling
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)
sc.pl.umap(
    sub_anndata, color="Sub_type",
    ax=ax, show=False, frameon=True, size=1.5, alpha=0.8,
    #palette=malignant_subtype_colors,
    legend_loc='right margin', legend_fontsize=7, legend_fontweight='normal'
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')
plt.tight_layout()

# Save with high quality settings for publication
plt.savefig(
    os.path.join(figurePath, f"{cell_type}_analysis", f'{cell_type} subtypes UMAP.pdf'), 
    dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf'
)
plt.show()
plt.close()

# Create a more appealing dot plot
selected_marker_genes = {
    # Pan-endothelial markers
    'Endo_Pan': ['PECAM1', 'CD34', 'VWF'],  # General endothelial markers
    # Endothelial subtypes
    'EC_EFNB2': ['EFNB2', 'HEY1'],  # Arterial endothelial
    'EC_ACKR1':["ACKR1"],
    'EC_ANGPT2': ['ANGPT2', 'NOTCH3', "VEGFA"],  # Venous endothelial
    'EC_LYVE1': ['LYVE1', 'SELP'],  # Capillary endothelial
    'EC_MKI67': ['MKI67']  # Proliferating endothelial
}

plt.figure(figsize=(15, 8))
sc.pl.dotplot(
    sub_anndata, var_names=selected_marker_genes, groupby="Sub_type",
    show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath,f"{cell_type}_analysis", f'{cell_type} select marker genes dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()
plt.close()

# new_cluster_names = [
#     '0', '1', '2', '3', '4', 
#     '5', '6', '7', '8', '9',
#     '10', '11', '12', '13', '14', 
#     '15', '16', '17', '18', '19',
#     '20', '21', '22', '23', '24', 
#     '25', '26', '27', '28', '29'
# ]