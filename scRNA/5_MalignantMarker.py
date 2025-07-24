#%% Comprehensive Malignant Cells Analysis Pipeline
import os
import gc
import pickle
from re import T
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ranksums # For Wilcoxon rank-sum test (Mann-Whitney U test)

import warnings
warnings.filterwarnings('ignore')

from malignant_functions import *

print("🧬 Starting Comprehensive Malignant Cells Analysis Pipeline")
print("="*70)

# Set paths
DataPath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","NMF_results")
result_path = figurePath

# create directories if they don't exist
os.makedirs(result_path, exist_ok=True)
os.makedirs(figurePath, exist_ok=True)

#%% Task 1: Extract Malignant-Only Cells with Original Expression
print("\n🎯 Task 1: Extracting Malignant-Only Cells")
print("="*60)

# Load integrated data if not already loaded
if 'final_integrated' not in locals():
    final_integrated = sc.read_h5ad(f"{DataPath}/major_anno_all.h5ad")

# Check major cell types
print("Major cell type distribution:")
print(final_integrated.obs['Major_type'].value_counts())
print(final_integrated.obs['tissue'].value_counts())

# Filter for tumor-origin epithelial cells only (malignant)
malignant_mask = (final_integrated.obs['Major_type'] == 'Epithelial') & \
                 (final_integrated.obs['tissue'] != 'PT')

print(f"Tumor epithelial cells: {malignant_mask.sum()}")

# Create malignant-only object with original expression matrix
malignant_cells = sc.AnnData(
    X=final_integrated[malignant_mask].raw.X,
    obs=final_integrated[malignant_mask].obs.copy(),
    var=final_integrated.raw.var.copy()
)
malignant_cells.obs['Major_type'] = 'Malignant'

print(f"\n✅ Extracted malignant cells:")
print(f"   Total malignant cells: {malignant_cells.n_obs:,}")
print(f"   Total features: {malignant_cells.n_vars:,}")
print(f"   Batch distribution:")
print(malignant_cells.obs['study'].value_counts())

# Save malignant-only object
malignant_cells.write(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))
print(f"✅ Saved malignant cells to: {result_path}/malignant_epithelial_cells.h5ad")

#%% Task 2: Downsample and Visualize Batch Effects
print("\n🎨 Task 2: Downsampling and Batch Visualization")
print("="*60)

if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))

# Downsample for visualization (stratified by batch)
malignant_subset = malignant_cells.copy()

# Basic preprocessing for visualization
sc.pp.normalize_total(malignant_subset, target_sum=1e4)
sc.pp.log1p(malignant_subset)
sc.pp.highly_variable_genes(malignant_subset, n_top_genes=2000)
sc.pp.scale(malignant_subset)
sc.tl.pca(malignant_subset)
sc.pp.neighbors(malignant_subset, n_neighbors=15)
sc.tl.umap(malignant_subset)

# Create batch visualization
fig, axes = plt.subplots(2, 2, figsize=(18, 12))

# Batch effect
sc.pl.umap(malignant_subset, color='batch', ax=axes[0,0], show=False, 
          legend_loc='right margin', title='Batch Effect')

# Patient ID (to check patient mixing)
sc.pl.umap(malignant_subset, color='sample_id', ax=axes[1,0], show=False,
          legend_loc=None, title='Patient Distribution')

# Quality metrics
sc.pl.umap(malignant_subset, color='nCount_RNA', ax=axes[0,1], show=False,
          title='Gene Counts', color_map='viridis')

sc.pl.umap(malignant_subset, color='percent.mt', ax=axes[1,1], show=False,
          title='MT Percentage', color_map='Reds')

plt.tight_layout()
plt.savefig(os.path.join(figurePath, "malignant_batch_effects.pdf"), dpi=300, bbox_inches='tight')
plt.show()

print("✅ Batch visualization saved")

#%% Task 3: Malignant Cell Batch correction using Harmony
print("\n🔄 Task 3: Batch Correction using Harmony")

# Prepare data with larger sample
if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))

malignant_cells.obs["sample_id"] = malignant_cells.obs["batch"]

# Basic preprocessing for visualization
sc.pp.normalize_total(malignant_cells, target_sum=1e4)
sc.pp.log1p(malignant_cells)
sc.pp.highly_variable_genes(malignant_cells, n_top_genes=2000)
sc.pp.scale(malignant_cells)
sc.tl.pca(malignant_cells)

sc.external.pp.harmony_integrate(
    malignant_cells, 
    key='batch',
    basis='X_pca',
    adjusted_basis='X_pca_harmony'
)

# Calculate UMAP embedding
sc.pp.neighbors(malignant_cells, n_neighbors=10, n_pcs=50, use_rep='X_pca_harmony')
resolution=0.8
key_added_ = f'leiden_res_{resolution}'
sc.tl.leiden(malignant_cells, resolution=resolution, key_added=key_added_)
malignant_cells.obs['leiden'] = malignant_cells.obs[key_added_]
sc.tl.umap(malignant_cells)
print("✅ UMAP calculation completed in downsample!")

# Comprehensive visualization
plt.rcParams.update({'font.size': 12}) # Set up the plotting parameters
cluster_colors = plt.cm.tab20(np.linspace(0, 1, len(malignant_cells.obs['leiden'].unique()))) # Define color palettes

# Create a comprehensive figure
fig = plt.figure(figsize=(15, 7))

# 1. Cluster visualization
ax1 = plt.subplot(1,2,1)
sc.pl.umap(malignant_cells, color='leiden', 
           legend_loc='on data', legend_fontsize=8, 
           ax=ax1, show=False, frameon=False)
ax1.set_title('Leiden Clusters', fontsize=14, fontweight='bold')

# 3. Batch effect
ax3 = plt.subplot(1,2,2)
sc.pl.umap(malignant_cells, color='batch', 
           ax=ax3, show=False, frameon=False)
ax3.set_title('Batch', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(figurePath, "umap of malignant cluster (harmony).pdf"), 
            dpi=300, bbox_inches='tight')
plt.show()


#%% Task 4: Malignant Cell Annotation
# Define comprehensive marker genes for major cell types
marker_genes = {
    # Epithelial cells
    'Malignant_Epithelial': ['EPCAM','KRT19','CDH1'],
    'Hypoxic': ['CA9','HIF1A'],  # CA9 encodes CA-IX protein
    'Glycolytic': ['SLC2A1', 'HK2', 'FASN'],  # SLC2A1 encodes GLUT1, metabolic genes
    'Proliferation': ['MKI67'],  # MKI67 encodes Ki-67, proliferation genes
    'EMT': ['VIM', 'FN1','VEGFA']  # EMT-like malignant cells
}

# Calculate marker gene scores for each cell type
print("Calculating marker gene scores...")
for cell_type, genes in marker_genes.items():
    # Filter genes that exist in the dataset
    available_genes = [g for g in genes if g in malignant_cells.var_names]
    
    if available_genes:
        print(f"  {cell_type}: {available_genes}")
        sc.tl.score_genes(malignant_cells, available_genes, 
                          score_name=f'{cell_type}_score', use_raw=False)
    else:
        print(f"  {cell_type}: No genes found in dataset")

# Get available marker genes
all_marker_genes = []
for genes in marker_genes.values():
    all_marker_genes.extend(genes)

print(f"Selected marker genes: {all_marker_genes}")

# Get all markers
sc.tl.rank_genes_groups(malignant_cells, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    malignant_cells, groupby="leiden", standard_scale="var", n_genes=8, show = False
)
plt.tight_layout()
plt.savefig(os.path.join(figurePath, "markers of all malignant cells.pdf"), 
            dpi=300, bbox_inches='tight')
plt.show()

# Create comprehensive cell type marker visualization
print("Creating marker gene expression plots...")
if all_marker_genes:
    fig = plt.figure(figsize=(24, 16))
    n_genes = len(all_marker_genes)
    n_cols = 6
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    for i, gene in enumerate(all_marker_genes):
        ax = plt.subplot(n_rows, n_cols, i+1)
        sc.pl.umap(malignant_cells, color=gene, use_raw=True, ax=ax, show=False, 
                  frameon=False, size=8, color_map='Reds')
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig(os.path.join(figurePath, 'malignant marker gene expression.pdf'), 
                dpi=300, bbox_inches='tight')
    plt.show()

# Annotation
new_cluster_names = [
    'EC_Quiescent', 'EC_Quiescent', 'EC_Quiescent', 'EC_EMT', 'EC_EMT', 
    'EC_Glycolysis', 'EC_Hypoxia', 'EC_EMT', 'EC_Proliferation', 'EC_Glycolysis', 
    'EC_EMT', 'EC_Hypoxia', 'EC_Quiescent', 'EC_Hypoxia', 'EC_EMT', 
    'Unknown', 'EC_EMT', 'Unknown', 'Unknown', 'EC_Glycolysis',
    'EC_Glycolysis'
]
cluster_to_celltype = {str(i): celltype for i, celltype in enumerate(new_cluster_names)}

malignant_cells.obs['Malignant_type'] = malignant_cells.obs['leiden'].map(cluster_to_celltype)

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

malignant_subtype_colors = {
    'EC_Quiescent': '#2E86AB',      # Blue - Epithelial marker
    'EC_EMT': '#E63946',        # Red - EMT/Mesenchymal marker   
    'EC_Hypoxia': '#8E44AD',        # Purple - Hypoxia marker
    'EC_Proliferation': '#F1C40F',      # Gold - Signaling/Immune marker
    'EC_Glycolysis': '#229954',        # Dark green - Metabolic marker (similar to FASN)
    'Unknown': '#BDC3C7'        # Light gray
}

malignant_cells.uns['Malignant_type_colors'] = malignant_subtype_colors

# Set up the figure with appropriate size for Nature/Cell (typically 85mm or 180mm width)
fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)

# Create the UMAP plot with enhanced styling
sc.pl.umap(
    malignant_cells, 
    color="Malignant_type",
    ax=ax,
    show=False,
    frameon=True,
    size=1.5,  # Smaller point size for 30k cells to avoid overcrowding
    alpha=0.8,  # Slight transparency for better visualization
    palette=malignant_subtype_colors,
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
    os.path.join(figurePath, 'malignant subtypes.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    facecolor='white',
    edgecolor='none',
    format='pdf'
)
plt.show()
plt.close()

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))

sc.pl.dotplot(
    malignant_cells, 
    var_names=all_marker_genes, groupby="Malignant_type",
    show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1
)

# Add styling after the plot is created
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=10)
plt.title('Expression Level', fontsize=12, fontweight='bold', pad=5)
plt.tight_layout()
plt.savefig(
    os.path.join(figurePath, 'malignant_marker_genes_dotplot.pdf'), 
    dpi=300, 
    bbox_inches='tight',
    format='pdf'
)
plt.show()

## Save the final annotated malignant cells object
malignant_cells.write(os.path.join(result_path, "malignant_epithelial_cells_annotated.h5ad"))

#%% Task 5: Malignant Cell Marker pathway analysis
# Scanpy-based Pathway Enrichment and Bulk Survival Analysis
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

from malignant_functions_marker import *

save_path="./malignant_subtype_analysis"
os.makedirs(save_path, exist_ok=True)

print("="*60)
print("MALIGNANT SUBTYPE COMPREHENSIVE ANALYSIS")
print("="*60)

# Step 1: Identify marker genes
if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells_annotated.h5ad"))

malignant_cells_back = malignant_cells.copy()
malignant_cells = malignant_cells[malignant_cells.obs['Malignant_type'] != 'Unknown']
marker_dict = identify_subtype_markers(malignant_cells,min_pct=0.25,logfc_threshold=1)

# Step 2: Get HALLMARK pathways
hallmark_pathways = get_hallmark_pathways(use_IMC_pathway=True)

# Step 3: Pathway enrichment analysis
enrichment_results = perform_pathway_enrichment(marker_dict, hallmark_pathways, 
                                                gene_size = malignant_cells.shape[1],
                                                n_top_gene = 50)

# Step 4: Create pathway enrichment figure
pathway_fig, enrichment_results = create_pathway_enrichment_figure(
    enrichment_results, save_path=save_path
)

# Step 5: Bulk RNA-seq data loading and survival analysis
bulk_datapath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/bulkRNA"
dataIDs = ["BCGSC","E-TABM-1112","E-MTAB-1951"]

bulk_datasets = {}

for dataID in dataIDs:
    bulk_exp = pd.read_csv(f"{bulk_datapath}/{dataID}/exp_data.csv", index_col=0)
    bulk_clinical = pd.read_csv(f"{bulk_datapath}/{dataID}/clinical_data.csv", index_col=0)

    if dataID == "E-TABM-1112":
        bulk_clinical["Sample_ID"] = bulk_clinical.index.tolist()

    elif dataID == "E-MTAB-1951":
        bulk_clinical['DiseaseState'] = bulk_clinical['risk_score_factor']
        bulk_clinical["Sample_ID"] = bulk_clinical.index.tolist()

    bulk_datasets[dataID] = (bulk_exp, bulk_clinical)

# Step 5: Bulk survival analysis (if bulk data provided)
if bulk_datasets:
    print("\n" + "="*40)
    print("BULK SURVIVAL ANALYSIS")
    print("="*40)
    
    for dataset_name, (bulk_exp, bulk_clinical) in bulk_datasets.items():
        print(f"\nAnalyzing dataset: {dataset_name}")
        
        # Calculate subtype scores
        scores_df = calculate_subtype_scores_ssgsea(bulk_exp, marker_dict)
        # scores_df = calculate_subtype_scores_zscore(bulk_exp, marker_dict,score_method="zscore_mean")
        
        # Check if survival data is available
        if 'OS_time' in bulk_clinical.columns and 'OS_status' in bulk_clinical.columns:
            survival_results = perform_survival_analysis(
                scores_df, bulk_clinical, save_path=f"{save_path}/{dataset_name}"
            )
            
            print(f"Survival analysis results for {dataset_name}:")
            for subtype, results in survival_results.items():
                significance = "***" if results['p_value'] < 0.001 else \
                                "**" if results['p_value'] < 0.01 else \
                                "*" if results['p_value'] < 0.05 else "ns"
                print(f"  {subtype}: p = {results['p_value']:.3f} {significance}")
        
        elif 'DiseaseState' in bulk_clinical.columns:
            # Perform disease state comparison
            comparison_results = compare_bulk_datasets(scores_df, bulk_clinical,save_path=f"{save_path}/{dataset_name}")
            
            print(f"Disease state comparison for {dataset_name}:")
            for subtype, results in comparison_results.items():
                significance = "***" if results['p_value'] < 0.001 else \
                                "**" if results['p_value'] < 0.01 else \
                                "*" if results['p_value'] < 0.05 else "ns"
                print(f"  {subtype}: p = {results['p_value']:.3f} {significance}")

print("\n" + "="*60)
print("ANALYSIS COMPLETE!")
print(f"Results saved to: {save_path}")
print("="*60)

print("Scanpy pathway and survival analysis functions ready!")
print("Main function: run_complete_analysis(malignant_cells, bulk_datasets, save_path)")
print("This will perform:")
print("1. Marker gene identification for each Malignant_type")
print("2. HALLMARK pathway enrichment analysis") 
print("3. Bulk survival analysis (if data provided)")
print("4. Comprehensive visualizations and statistical tests")