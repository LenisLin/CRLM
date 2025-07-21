import os
import gc

import scanpy as sc

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.utils import resample

import warnings
warnings.filterwarnings('ignore')

## Function to load and process scRNA-seq datasets
from annotation_functions import *

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

## Load integrated data
## Set workDir
DataPath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA")

final_integrated = sc.read_h5ad(f"{DataPath}/final_integrated.h5ad")

#%% Step 0: Remove doublets
print("🔍 Step 0: Removing doublets...")
print("="*60)

## subsample doublet prediction
target_cells = 20000
sampled_indices = np.random.choice(final_integrated.shape[0], size=target_cells, replace=False)
downsampled_adata = final_integrated[sampled_indices].copy()

sc.pp.neighbors(downsampled_adata, n_neighbors=10, n_pcs=40, use_rep='X_pca_harmony')
sc.tl.umap(downsampled_adata)

## Doublet prediction visualize
downsampled_adata.obs['doublet'].value_counts(dropna=False)
final_integrated.obs.doublet = final_integrated.obs.doublet.replace({"singlet": "Singlet"})
downsampled_adata.obs.doublet = downsampled_adata.obs.doublet.replace({"singlet": "Singlet"})

doublet_colors = {'Singlet': '#2E86AB', 'Doublet': '#F24236'}  # Blue for singlet, Red for doublet
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(downsampled_adata, color='doublet', 
           palette=doublet_colors, ax=ax, show=False, frameon=False, size=30)
ax.set_title('Predicted Doublets', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(figurePath, 'doublet_prediction_umap.pdf'), dpi=300, bbox_inches='tight')
plt.show()

## Print doublet statistics
n_doublets = (final_integrated.obs['doublet'] == "Doublet").sum()
n_total = final_integrated.shape[0]
doublet_rate = n_doublets / n_total * 100

print(f"📊 Doublet Statistics:")
print(f"   Total cells: {n_total:,}")
print(f"   Predicted doublets: {n_doublets:,}")
print(f"   Doublet rate: {doublet_rate:.2f}%")

## Remove doublet in integrated dataset
final_integrated = final_integrated[final_integrated.obs['doublet'] != "Doublet"].copy()

#%% Step 1: Leiden clustering on all cells
print("🔬 Step 1: Performing Leiden clustering on all cells...")
print("="*60)

# Build neighborhood graph if not already done
if 'neighbors' not in final_integrated.uns:
    print("Building neighborhood graph...")
    sc.pp.neighbors(final_integrated, n_neighbors=15, n_pcs=50, use_rep='X_pca_harmony')

# Perform leiden clustering with multiple resolutions
resolutions = [0.5] # [0.3, 0.5, 0.8, 1.0, 1.2]
for res in resolutions:
    print(f"Leiden clustering with resolution {res}...")
    sc.tl.leiden(final_integrated, resolution=res, key_added=f'leiden_res_{res}')

# Use resolution 0.5 as default (you can adjust based on results)
final_integrated.obs['leiden'] = final_integrated.obs['leiden_res_0.5']

print(f"✅ Clustering completed!")
print(f"Number of clusters (res=0.5): {len(final_integrated.obs['leiden'].unique())}")
print(f"Cluster distribution:")
print(final_integrated.obs['leiden'].value_counts().sort_index())

#%% Step 2: Downsample cells
print("\n🎯 Step 2: Downsampling...")
print("="*60)

# Check if we have enough cells
n_cells = final_integrated.shape[0]
target_cells = 30000

if n_cells <= target_cells:
    print(f"Dataset has {n_cells} cells, which is ≤ {target_cells}. Using all cells.")
    downsampled_adata = final_integrated.copy()
else:
    print(f"Downsampling from {n_cells} to {target_cells} cells...")
    
    # Stratified sampling to maintain cluster proportions
    cluster_counts = final_integrated.obs['leiden'].value_counts()
    cluster_proportions = cluster_counts / cluster_counts.sum()
    
    # Calculate target cells per cluster
    target_per_cluster = (cluster_proportions * target_cells).round().astype(int)
    
    # Ensure we don't exceed available cells per cluster
    target_per_cluster = np.minimum(target_per_cluster, cluster_counts)
    
    # Adjust if total is not exactly 40000
    diff = target_cells - target_per_cluster.sum()
    if diff > 0:
        # Add extra cells to largest clusters
        largest_clusters = cluster_counts.nlargest(diff).index
        target_per_cluster[largest_clusters] += 1
    elif diff < 0:
        # Remove cells from largest clusters
        largest_clusters = cluster_counts.nlargest(-diff).index
        target_per_cluster[largest_clusters] -= 1
    
    print("Target cells per cluster:")
    for cluster, target in target_per_cluster.items():
        print(f"  Cluster {cluster}: {target} cells")
    
    # Perform stratified sampling
    sampled_indices = []
    for cluster in final_integrated.obs['leiden'].unique():
        cluster_mask = final_integrated.obs['leiden'] == cluster
        cluster_indices = np.where(cluster_mask)[0]
        
        n_target = target_per_cluster[cluster]
        if len(cluster_indices) >= n_target:
            sampled_cluster_indices = np.random.choice(
                cluster_indices, size=n_target, replace=False
            )
        else:
            sampled_cluster_indices = cluster_indices
        
        sampled_indices.extend(sampled_cluster_indices)
    
    # Create downsampled dataset
    downsampled_adata = final_integrated[sampled_indices].copy()
    
    print(f"✅ Downsampling completed!")
    print(f"Final sample size: {downsampled_adata.shape[0]} cells")
    print(f"Maintained cluster proportions:")
    print(downsampled_adata.obs['leiden'].value_counts().sort_index())

#%% Step 3: Calculate UMAP
print("\n🗺️ Step 3: Calculating UMAP...")
print("="*60)

# Calculate UMAP embedding
sc.pp.neighbors(downsampled_adata, n_neighbors=10, n_pcs=40, use_rep='X_pca_harmony')
sc.tl.umap(downsampled_adata)
print("✅ UMAP calculation completed in downsample!")

#%% Step 4: Comprehensive visualization
print("\n🎨 Step 4: Creating visualizations...")
print("="*60)

# Set up the plotting parameters
plt.rcParams.update({'font.size': 12})

# Define color palettes
cluster_colors = plt.cm.tab20(np.linspace(0, 1, len(downsampled_adata.obs['leiden'].unique())))

# Create a comprehensive figure
fig = plt.figure(figsize=(25, 7))

# 1. Cluster visualization
ax1 = plt.subplot(1, 3, 1)
sc.pl.umap(downsampled_adata, color='leiden', 
           legend_loc='on data', legend_fontsize=8, 
           ax=ax1, show=False, frameon=False)
ax1.set_title('Leiden Clusters', fontsize=14, fontweight='bold')

# 2. Tissue_Type
ax2 = plt.subplot(1, 3, 2)
sc.pl.umap(downsampled_adata, color='tissue', 
           ax=ax2, show=False, frameon=False)
ax2.set_title('Tissue_Type', fontsize=14, fontweight='bold')

# 3. Batch effect
ax3 = plt.subplot(1, 3, 3)
sc.pl.umap(downsampled_adata, color='batch', 
           ax=ax3, show=False, frameon=False)
ax3.set_title('Batch', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(figurePath, "comprehensive_umap_visualization.pdf"), 
            dpi=300, bbox_inches='tight')
plt.show()

#%% Additional summary statistics
print("\n📊 Summary Statistics:")
print("="*40)
print(f"Total cells after downsampling: {downsampled_adata.shape[0]}")
print(f"Total features: {downsampled_adata.shape[1]}")
print(f"Number of clusters: {len(downsampled_adata.obs['leiden'].unique())}")

print("\nCluster sizes:")
cluster_sizes = downsampled_adata.obs['leiden'].value_counts().sort_index()
for cluster, size in cluster_sizes.items():
    percentage = size / len(downsampled_adata) * 100
    print(f"  Cluster {cluster}: {size} cells ({percentage:.1f}%)")

print("\nBatch distribution:")
batch_dist = downsampled_adata.obs['batch'].value_counts()
for batch, count in batch_dist.items():
    percentage = count / len(downsampled_adata) * 100
    print(f"  {batch}: {count} cells ({percentage:.1f}%)")

#%% Save the downsampled dataset
print(f"\n💾 Saving downsampled dataset...")
downsampled_adata.write(os.path.join(DataPath, "downsampled_integrated_30k.h5ad"))
print("✅ Dataset saved!")

print(f"\n🎉 Analysis completed successfully!")
print(f"📁 Visualizations saved to: {figurePath}")
print(f"📁 Data saved to: {DataPath}")

#%% Step 5: Annotation
print("🧬 Part 5: Cell Type Annotation using Marker Genes")
print("="*60)

# Define comprehensive marker genes for major cell types
marker_genes = {
    # Epithelial cells
    'Epithelial': ['EPCAM','KRT19','CDH1','ALB'],
    # 'Malignant_Epithelial': ['EPCAM', 'KRT19', 'KRT8', 'KRT18', 'CDH1']
    # 'Hepatocytes': ['ALB', 'APOE', 'CYP3A4', 'G6PC', 'PCK1']
    # T cells

    'T': ['CD3E'],
    'NKT': ['GNLY', 'NKG7'],

    # B cells
    'B': ['CD19', 'MS4A1'],  # MS4A1 is CD20
    'Plasma': ['IGKC'],  # MS4A1 is CD20

    # Myeloid cells
    'Myeloid': ['CD14', 'CD68', 'LYZ'],
    'cDC': ['ITGAX','CLEC9A'],  # ITGAX is CD11c
    
    # Neutrophils (important in cancer)
    'Neutrophils': ['FCGR3B', 'S100A9'],

    # Stromal cells
    'Stromal': ['COL1A2', 'FAP','POSTN'],
    
    # Endothelial cells
    'Endothelial': ['PECAM1', 'VWF'],  # PECAM1 is CD31
    
    # Mast cells
    'Mast_cells': ['KIT']  # KIT is CD117
}

# Calculate marker gene scores for each cell type
print("Calculating marker gene scores...")
for cell_type, genes in marker_genes.items():
    # Filter genes that exist in the dataset
    available_genes = [g for g in genes if g in downsampled_adata.var_names]
    
    if available_genes:
        print(f"  {cell_type}: {available_genes}")
        sc.tl.score_genes(downsampled_adata, available_genes, 
                          score_name=f'{cell_type}_score', use_raw=False)
    else:
        print(f"  {cell_type}: No genes found in dataset")

# 2.2 Create comprehensive cell type marker visualization
print("Creating marker gene expression plots...")

# Get available marker genes
all_marker_genes = []
for genes in marker_genes.values():
    all_marker_genes.extend(genes)

print(f"Selected marker genes: {all_marker_genes}")

# 2.3 Exam the marker of cluster 4
sc.tl.rank_genes_groups(downsampled_adata, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    downsampled_adata, groupby="leiden", standard_scale="var", n_genes=5
)
plt.tight_layout()
plt.savefig(os.path.join(figurePath, "markers of all clusters.pdf"), 
            dpi=300, bbox_inches='tight')
plt.show()

# Create marker gene expression plot
if all_marker_genes:
    fig = plt.figure(figsize=(24, 16))
    n_genes = len(all_marker_genes)
    n_cols = 6
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    for i, gene in enumerate(all_marker_genes):
        ax = plt.subplot(n_rows, n_cols, i+1)
        sc.pl.umap(downsampled_adata, color=gene, use_raw=True, ax=ax, show=False, 
                  frameon=False, size=8, color_map='Reds')
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig(os.path.join(figurePath, 'marker_genes_expression.pdf'), 
                dpi=300, bbox_inches='tight')
    plt.show()


#%% 3. Annotation
new_cluster_names = [
    'T', 'T', 'Myeloid', 'T', 'Epithelial', 
    'NK', 'NK', 'T', 'Fibroblast', 'B', 
    'Plasma', 'T', 'Endothelial', 'Myeloid', 'Myeloid', 
    'Mast', 'Hepatocyte', 'T', 
    ]
cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}

downsampled_adata.obs['Major_type'] = downsampled_adata.obs['leiden'].map(cluster_to_celltype)
final_integrated.obs['Major_type'] = final_integrated.obs['leiden'].map(cluster_to_celltype)

# Set up Nature/Cell journal style
plt.rcParams.update({
    'font.size': 8,
    'axes.linewidth': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'xtick.minor.size': 1,
    'ytick.minor.size': 1,
    'legend.frameon': False,
    'legend.fontsize': 7,
    'pdf.fonttype': 42,  # Important for Nature journals
    'ps.fonttype': 42
})

cell_type_colors = {
    'B': '#FFD700',          # Gold/appealing yellow
    'Plasma': '#FFFF99',     # Light yellow
    'Epithelial': '#33A02C',  # Green
    'Hepatocyte': '#B2DF8A',  # Light green
    'Mast': '#B15928',       # Brown
    'Fibroblast': '#FF7F00',  # Orange
    'Endothelial': '#E31A1C', # Red - represents blood vessels
    'Myeloid': '#6A3D9A',    # Purple - strong immune color
    'NK': '#A6CEE3',         # Light blue
    'T': '#1F78B4',          # Blue
}

downsampled_adata.uns['Major_type_colors'] = cell_type_colors
final_integrated.uns['Major_type_colors'] = cell_type_colors

# Set up the figure with appropriate size for Nature/Cell (typically 85mm or 180mm width)
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)

# Create the UMAP plot with enhanced styling
sc.pl.umap(
    downsampled_adata, 
    color="Major_type",
    ax=ax,
    show=False,
    frameon=True,
    size=1.5,  # Smaller point size for 30k cells to avoid overcrowding
    alpha=0.8,  # Slight transparency for better visualization
    palette=cell_type_colors,
    legend_loc='right margin',
    legend_fontsize=7,
    legend_fontweight='normal'
)

# Add a subtle border
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color('black')

# Adjust layout
plt.tight_layout()

# Save with high quality settings for publication
plt.savefig(
    os.path.join(figurePath, 'major_type_annotation_30k_nature_style.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    facecolor='white',
    edgecolor='none',
    format='pdf'
)

print(f"\n💾 Saving the annotated dataset...")
downsampled_adata.write(os.path.join(DataPath, "major_anno_30k.h5ad"))
final_integrated.write(os.path.join(DataPath, "major_anno_all.h5ad"))
print("✅ Dataset saved!")

# Create Marker Gene Expression Dot Plot
print("🎯 Creating marker gene dot plot...")

# Flatten marker genes list for easier plotting
all_marker_genes = []
for cell_type, genes in marker_genes.items():
    all_marker_genes.extend(genes)

# Remove duplicates while preserving order
all_marker_genes = list(dict.fromkeys(all_marker_genes))

# Define the order of cell types to match marker gene groupings
cell_type_order = [
    'Epithelial',    # EPCAM, KRT19, CDH1, ALB
    'Hepatocyte',    # ALB (shared with Epithelial)
    'T',             # CD3E
    'NK',            # GNLY, NKG7
    'B',             # CD19, MS4A1
    'Plasma',        # IGKC
    'Myeloid',       # CD14, CD68, LYZ, ITGAX, CLEC9A, FCGR3B, S100A9
    'Fibroblast',    # COL1A2, FAP, POSTN
    'Endothelial',   # PECAM1, VWF
    'Mast'           # KIT
]

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))

sc.pl.dotplot(
    downsampled_adata, 
    var_names=all_marker_genes, 
    groupby="Major_type",
    show=False,
    categories_order=cell_type_order,
    standard_scale='var',
    color_map='Reds',
    dot_max=0.8,
    dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)

plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, 'marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()
print("✅ Enhanced marker gene dot plot created!")

#%% 4.1 Create some data descriptions
# Explore all clinical variables
key_vars = ['Major_type', 'Tissue_Type', 'Treatment_Strategy', 'Microsatellite_Status', 
           'TNM', 'Tumor_stage', 'Gender', 'Age', 'Response', 'Treatment_Stage']

results = explore_clinical_data(final_integrated, key_vars)
patient_stats = explore_patient_level(final_integrated)

# Set Cancer Cell journal style
plt.rcParams.update({
    'font.size': 7,
    'font.family': 'Arial',
    'axes.linewidth': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'legend.frameon': False,
    'legend.fontsize': 6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 300
})

# Save all color mappings
cancer_cell_colors = color_setting()
save_color_maps(final_integrated, cancer_cell_colors)
print("🎨 Creating Cancer Cell style descriptive figures...") # Execute all plots

# Generate all fancy figures
print("📊 Figure 1: Fancy Cell Type Composition (Treemap + Donut)...")
plot_fancy_cell_type_composition(final_integrated, 
                                save_path=os.path.join(figurePath, 'Fig1_fancy_cell_composition.pdf'))

print("📊 Figure 2: Clinical Overview (Violin + Bubble + Ridgeline)...")
plot_fancy_clinical_overview(final_integrated,
                            save_path=os.path.join(figurePath, 'Fig2_fancy_clinical_overview.pdf'))

print("📊 Figure 3: Alluvial Flow Diagram...")
plot_alluvial_celltype_tissue(final_integrated,
                             save_path=os.path.join(figurePath, 'Fig3_alluvial_celltype_tissue.pdf'))

print("📊 Figure 4: Network Treatment Response (Waterfall + Sunburst)...")
plot_network_treatment_response(final_integrated,
                               save_path=os.path.join(figurePath, 'Fig4_network_treatment_response.pdf'))

print("📊 Figure 5: Patient Journey Overview (Timeline + Matrix)...")
plot_patient_journey_overview(final_integrated,
                             save_path=os.path.join(figurePath, 'Fig5_patient_journey_overview.pdf'))

print("✅ All fancy figures generated and saved!")
print(f"📁 Figures saved in: {figurePath}/")
print("🎨 Color mappings saved in AnnData.uns")

# Reset matplotlib to defaults
plt.rcdefaults()