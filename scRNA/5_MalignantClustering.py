#%% Comprehensive Malignant Cells Analysis Pipeline
import os
import gc
import pickle
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ranksums # For Wilcoxon rank-sum test (Mann-Whitney U test)

import warnings
warnings.filterwarnings('ignore')

from malignant_functions_marker import *

print("🧬 Starting Comprehensive Malignant Cells Analysis Pipeline")
print("="*70)

# Set paths
DataPath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","Malignant_Analysis")
result_path = figurePath

# create directories if they don't exist
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
malignant_cells.obs["sample_id"] = malignant_cells.obs["batch"]

malignant_cells.raw = malignant_cells

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
sc.pl.umap(malignant_subset, color='study', ax=axes[1,0], show=False,
          legend_loc='right margin', title='Cohort Distribution')

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

# Manual HVG selection based on your IMC findings:
metabolism_gene_sets = get_imc_aligned_hallmark_pathways()
imc_guided_genes = ['EPCAM','CDH1','KRT19','CA9']
for pathway_name, genes in metabolism_gene_sets.items():
    imc_guided_genes.extend([g for g in genes if g in malignant_cells.var_names])

# Remove duplicates
imc_guided_genes = list(set(imc_guided_genes))

# Basic preprocessing for visualization
sc.pp.normalize_total(malignant_cells, target_sum=1e4)
sc.pp.log1p(malignant_cells)
# sc.pp.highly_variable_genes(malignant_cells, n_top_genes=2000)
malignant_cells = manual_hvg(malignant_cells.copy(),imc_guided_genes)
sc.pp.scale(malignant_cells)
sc.tl.pca(malignant_cells)
sc.external.pp.harmony_integrate(
    malignant_cells, 
    key='batch',
    basis='X_pca',
    adjusted_basis='X_pca_harmony'
)

# Calculate UMAP embedding
sc.pp.neighbors(malignant_cells, n_neighbors=10, n_pcs=50,  use_rep='X_pca_harmony')
resolution=0.6
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
    "Epithelial": ['EPCAM','CDH1'],
    'Malignant': ['CDX2','KRT20','KRT19'],
    'Other': ['Hep-Par1','KRT7'],
    'Hypoxic': ['CA9', 'HIF1A', 'BNIP3'],  # CA9 encodes CA-IX protein
    'Lipids_metabolism': ['FASN', 'ACLY', 'SREBF1'],
    'Glycolytic': ['SLC2A1', 'HK2', 'PFKP', 'ALDOA'],  # SLC2A1 encodes GLUT1, metabolic genes
    'Proliferation': ['MKI67', 'PCNA', 'CCNB1', 'CDC20', 'TOP2A'],  # MKI67 encodes Ki-67, proliferation genes
    'EMT': ['VIM', 'FN1']  # EMT-like malignant cells
}

# Get available marker genes
all_marker_genes = []
for genes in marker_genes.values():
    all_marker_genes.extend(genes)

print(f"Selected marker genes: {all_marker_genes}")

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

# Get all markers
sc.tl.rank_genes_groups(malignant_cells, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    malignant_cells, groupby="leiden", standard_scale="var", n_genes=8, show = False
)
plt.tight_layout()
plt.savefig(os.path.join(figurePath, "markers of all malignant cells.pdf"), 
            dpi=300, bbox_inches='tight')
plt.show()

# Create UMAP for cell Marker visualization
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

# Create a more appealing dot plot
plt.figure(figsize=(15, 8))
sc.pl.dotplot(malignant_cells, var_names=all_marker_genes, groupby=key_added_,show=False,standard_scale='var',color_map='Reds',dot_max=0.8,dot_min=0.1)
plt.xticks(rotation=45, fontsize=10) # Add styling after the plot is created
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

# Annotation
## create 
new_cluster_names = [
    'EC_LipidMeta', 'EC_Glycolysis', 'EC_Quiescent', 'EC_Quiescent', 'EC_Quiescent',
    'EC_Glycolysis', 'EC_Proliferation', 'EC_LipidMeta', 'EC_EMT', 'Unknown'
    # 'EC_Quiescent', 'EC_EMT', 'EC_Glycolysis', 'EC_Hypoxia', 'EC_EMT', 'EC_Proliferation', 'Unknown'
]
# new_cluster_names = ["EC"+"_" + str(x) for x in range(0, 20)]
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
    'EC_LipidMeta': '#8E44AD',        # Purple - Hypoxia marker
    'EC_Proliferation': '#F1C40F',      # Gold - Signaling/Immune marker
    'EC_Glycolysis': '#229954',        # Dark green - Metabolic marker (similar to FASN)
    'Unknown': '#BDC3C7'        # Light gray
}
malignant_cells.uns['Malignant_type_colors'] = malignant_subtype_colors

## remove the Unknown type
malignant_cells = malignant_cells[malignant_cells.obs['Malignant_type'] != 'Unknown'] 

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
    #palette=malignant_subtype_colors,
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
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","Malignant_Analysis")
result_path = figurePath
save_path = os.path.join(figurePath,"pathway_analysis")
os.makedirs(save_path, exist_ok=True)

print("="*60)
print("MALIGNANT SUBTYPE COMPREHENSIVE ANALYSIS")
print("="*60)

# Step 1: Identify marker genes
if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells_annotated.h5ad"))

malignant_cells_back = malignant_cells.copy()
# malignant_cells = malignant_cells[malignant_cells.obs['Malignant_type'] != 'Unknown']
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

#%% Step 6: Bulk RNA-seq data loading and survival analysis
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

#%% Step 7: Transcriptomic Trajectory Analysis
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

from malignant_functions_marker import *
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","Malignant_Analysis")
result_path = figurePath
save_path = os.path.join(figurePath,"trajectory_analysis")
os.makedirs(save_path, exist_ok=True)

## Load malignant cells
if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells_annotated.h5ad"))

# Step 1: Preprocess data for trajectory analysis
sc.tl.draw_graph(malignant_cells)
# sc.tl.diffmap(malignant_cells)
# sc.pp.neighbors(malignant_cells, n_neighbors=10, use_rep="X_diffmap")
sc.pl.draw_graph(malignant_cells, color="Malignant_type", legend_loc="on data")

# Step 2: Clustering and PAGA
## Already annotation
sc.tl.paga(malignant_cells, groups="Malignant_type")
sc.pl.paga(malignant_cells, color=["Malignant_type"])

sc.tl.draw_graph(malignant_cells, init_pos="paga")
sc.pl.draw_graph(malignant_cells, color="Malignant_type", legend_loc="on data")

# Step 3: Reconstructing gene changes along PAGA paths for a given set of genes
malignant_cells.uns["iroot"] = np.flatnonzero(malignant_cells.obs["Malignant_type"] == "EC_Quiescent")[0]
sc.tl.dpt(malignant_cells)
sc.pl.draw_graph(malignant_cells, color=["Malignant_type", "dpt_pseudotime"], legend_loc="on data")

#%% Step 8: GRN
import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose
from collections import Counter
from IPython.display import HTML, display

import scipy.io
import math
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss
from pyscenic.cli.utils import load_signatures

sc.settings.njobs = 12 # Set maximum number of jobs for Scanpy.

figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","Malignant_Analysis")
result_path = figurePath
save_path = os.path.join(figurePath,"GRN_analysis")
os.makedirs(save_path, exist_ok=True)

# Define Folder structure.
RESOURCES_FOLDERNAME = "/home/lenislin/Experiment/projects/zsz_CRC/pySCENIC/resources"
RESULTS_FOLDERNAME = save_path
FIGURES_FOLDERNAME = save_path
if not os.path.exists(RESULTS_FOLDERNAME):
    os.makedirs(RESULTS_FOLDERNAME,exist_ok=True)

# Auxilliary functions.
BASE_URL = save_path
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"

# Download Auxilliary data sets.
# Downloaded fromm pySCENIC github repo: https://github.com/aertslab/pySCENIC/tree/master/resources
HUMAN_TFS_FNAME = os.path.join(RESOURCES_FOLDERNAME, 'allTFs_hg38.txt')
# Ranking databases. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
RANKING_DBS_FNAMES = list(map(lambda fn: os.path.join(RESOURCES_FOLDERNAME, fn),
                       ['hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
                        'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather']))
# Motif annotations. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDERNAME, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

################### pyscenic ###################

if 'adata' not in locals():
    adata = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells_annotated.h5ad"))

## Results save
DATASET_ID="CRLM_Malignant"

METADATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.metadata.csv'.format(DATASET_ID))
EXP_MTX_QC_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.qc.tpm.csv'.format(DATASET_ID))
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(DATASET_ID))
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.motifs.csv'.format(DATASET_ID))
REGULONS_DAT_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.regulons.dat'.format(DATASET_ID))
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.auc.csv'.format(DATASET_ID))
BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.bin.csv'.format(DATASET_ID))
THR_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.thresholds.csv'.format(DATASET_ID))
ANNDATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.h5ad'.format(DATASET_ID))
LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.loom'.format(DATASET_ID))

adata.write_h5ad(ANNDATA_FNAME) # Categorical dtypes are created.
adata.to_df().to_csv(EXP_MTX_QC_FNAME)

## Load functions
def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', -1)
    display(HTML(df.head().to_html(escape=False)))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)

def mapping_arr_to_matrix(exp_,colnames,meta_):
    # Create a mapping from labels to indices
    label_to_index = {label: idx for idx, label in enumerate(colnames)}

    # Find indices of the labels in characteristic_labels
    characteristic_labels = meta_.index
    col_indices = [label_to_index[label] for label in characteristic_labels if label in label_to_index]

    # Select the corresponding columns from the sparse matrix
    exp_ = exp_.tocsr()  # For efficient row slicing
    selected_columns = exp_[:, col_indices]

    # (Optional) Convert the result to a dense array for display
    # selected_columns_dense = selected_columns.toarray()
    
    return selected_columns

# REGULON CREATION
def derive_regulons(motifs, db_names=(#'hg19-tss-centered-10kb-10species', 'hg19-500bp-upstream-10species','hg19-tss-centered-5kb-10species'
        "hg38-limited-upstream500-tss-downstream100-full-transcript",
        "hg38-limited-upstream10000-tss-downstream10000-full-transcript"
                                      )):
    # motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]

    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]
    
    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0) 
                                                                      & ((motifs['Annotation'] == 'gene is directly annotated')
                                                                        | (motifs['Annotation'].str.startswith('gene is orthologous to')
                                                                           & motifs['Annotation'].str.endswith('which is directly annotated for motif')))
                                                                     ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))


def saveimg(fname: str, fig, folder: str=FIGURES_FOLDERNAME) -> None:
    """
    Save figure as vector-based PDF image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname))


def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

def clusterplot(meta, bin_mtx_, auc_mtx_, cellid_col = "cell_id", celltype_col = "sub_celltype", saveName = None):
        
    COLORS = [color['color'] for color in mpl.rcParams["axes.prop_cycle"]]

    cell_type_color_lut = dict(zip(meta[celltype_col].dtype.categories, COLORS))
    cell_id2cell_type_lut = meta.set_index(cellid_col)[celltype_col].to_dict()
    bw_palette = sns.xkcd_palette(["white", "black"])

    sns.set_style("whitegrid")
    fig = palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])
    saveimg(f'{saveName}-legend on_off.png', fig)

    sns.set(font_scale=0.8)
    fig = palplot(sns.color_palette(COLORS), meta[celltype_col].dtype.categories, size=1.0)
    saveimg(f'{saveName}-legend cell_type_colors.png', fig)

    sns.set(font_scale=1.0)
    sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})
    g = sns.clustermap(bin_mtx_.T, 
                col_colors=auc_mtx_.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
                cmap=bw_palette, figsize=(20,20))
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.set_xticks([])
    g.ax_heatmap.set_xlabel('Cells')
    g.ax_heatmap.set_ylabel('Regulons')
    g.ax_col_colors.set_yticks([0.5])
    g.ax_col_colors.set_yticklabels(['Cell Type'])
    g.cax.set_visible(False)
    g.savefig(os.path.join(FIGURES_FOLDERNAME, f'{saveName}-GRN clustermap.pdf'), format='pdf')

    return None

## STEP 1: Network inference based on GRNBoost2 from CLI
os.system(f"pyscenic grn {EXP_MTX_QC_FNAME} {HUMAN_TFS_FNAME} -o {ADJACENCIES_FNAME} --num_workers 8")

## STEP 2-3: Regulon prediction aka cisTarget from CLI
DBS_PARAM = ' '.join(RANKING_DBS_FNAMES)
os.system(f"pyscenic ctx {ADJACENCIES_FNAME} {DBS_PARAM} --annotations_fname {MOTIF_ANNOTATIONS_FNAME} --expression_mtx_fname {EXP_MTX_QC_FNAME} --output {MOTIFS_FNAME} --num_workers 12")

df_motifs = load_motifs(MOTIFS_FNAME)
df_motifs.head()

## STEP 4: Cellular enrichment aka AUCell (skip this step: https://github.com/aertslab/pySCENIC/issues/199)
# regulons = derive_regulons(df_motifs)
regulons = load_signatures(MOTIFS_FNAME)
with open(REGULONS_DAT_FNAME, 'wb') as f:
    pickle.dump(regulons, f)

with open(REGULONS_DAT_FNAME, 'rb') as f:
    regulons = pickle.load(f)

# AUCELL
df_tpm = pd.read_csv(EXP_MTX_QC_FNAME, index_col=0)
auc_mtx = aucell(df_tpm, regulons, num_workers=12) 
auc_mtx.to_csv(AUCELL_MTX_FNAME)

## OPTIONAL STEP 5 - Regulon activity binarization¶
auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0)

bin_mtx, thresholds = binarize(auc_mtx) 
bin_mtx.to_csv(BIN_MTX_FNAME) 
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME)

bin_mtx = pd.read_csv(BIN_MTX_FNAME, index_col=0)
thresholds = pd.read_csv(THR_FNAME, index_col=0).threshold

sort_auc_mtx = pd.DataFrame(data=np.sum(auc_mtx,axis=0))
sort_auc_mtx = sort_auc_mtx.sort_values(by=0, ascending=False)

sort_auc_mtx.to_csv(os.path.join(FIGURES_FOLDERNAME,'Sorted AUC of GRN in Malignant cells.csv'))
top_grn = [x for x in sort_auc_mtx.head(8).index]

fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(8, 4), dpi=100)
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

for i, grn in enumerate(top_grn):
    # Select the corresponding axis (ax1 to ax8)
    ax = axes[i]
    plot_binarization(auc_mtx, grn, thresholds[grn], ax=ax)

saveimg(os.path.join(FIGURES_FOLDERNAME,'GRN of Malignant cells.pdf'), fig)

## Create heatmap with binarized regulon activity
# adata = sc.read_h5ad(ANNDATA_FNAME)
# df_metadata = adata.obs
# df_metadata.head()
# df_metadata["cell_id"] = df_metadata.index

# meta = df_metadata[df_metadata['sub_celltype'].str.startswith(target_type, na=False)]
# meta['sub_celltype'] = meta['sub_celltype'].cat.remove_unused_categories()

# bin_mtx_ = bin_mtx.loc[meta.index,]
# auc_mtx_ = auc_mtx.loc[meta.index,]
# 
# clusterplot(meta=meta,bin_mtx_=bin_mtx_,auc_mtx_=auc_mtx_)

## STEP 6: Non-linear projection and clustering
adata = sc.read_h5ad(ANNDATA_FNAME)
auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0)
with open(REGULONS_DAT_FNAME, 'rb') as f:
    regulons = pickle.load(f)
f.close()

add_scenic_metadata(adata, auc_mtx, regulons)

rss = regulon_specificity_scores(auc_mtx, adata.obs.Malignant_type)
rss.head()
rss.to_csv(os.path.join(FIGURES_FOLDERNAME, "rss.csv"))

## CELL TYPE SPECIFIC REGULATORS - RSS
# List of features to plot
# features = [
#     'Mph_APOE', 'Mph_CCL20', 'Mph_S100A8', 'Mph_SPP1',
#     'Fibro_CCL11', 'Fibro_RGS5', 'Fibro_CXCL8', 'Fibro_MYH11'
# ]
features = [x for x in set(adata.obs["Malignant_type"])]
len_feat = math.ceil(len(features) / 2)

# Flatten the axes for easier indexing
sns.set(style='whitegrid', font_scale=1)
fig, axes = plt.subplots(2, len_feat, figsize=((len_feat + 0.5)*2, 8), dpi=300) # Create a 2x4 grid for 8 plots

axes = axes.flatten()

# Plot each feature
for ax, feature in zip(axes, features):
    plot_rss(rss, feature, ax=ax)
    ax.set_xlabel('')
    ax.set_ylabel('')

plt.tight_layout()
saveimg(os.path.join(FIGURES_FOLDERNAME,'GRN-RSS of Malignant cells.pdf'), fig)

## CELL TYPE SPECIFIC REGULATORS - Z-SCORE
df_obs = adata.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['Malignant_type']]
df_results = ((df_scores.groupby(by='Malignant_type').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results[(df_results.Z >= 2.0)].sort_values('Z', ascending=False)

# Create pivot table with only top regulons
top_5_regulons = df_results.groupby('Malignant_type').apply(lambda x: x.nlargest(5, 'Z')).reset_index(drop=True)
df_heatmap_top5 = pd.pivot_table(data=top_5_regulons,index='Malignant_type', columns='regulon', values='Z')

# Plot heatmap
fig, ax1 = plt.subplots(1, 1, figsize=(12, 8))
sns.heatmap(df_heatmap_top5, ax=ax1, annot=True, fmt=".1f", linewidths=.7, 
            cbar=False, square=True, linecolor='gray', 
            cmap="YlGnBu", annot_kws={"size": 6})
ax1.set_ylabel('')
saveimg('Top5_Malignant_cell_regulons_heatmap.pdf', fig)
