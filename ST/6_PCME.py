"""
PCME analysis on stereo-seq data
"""
import os
import pickle
import pandas as pd
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_ind, mannwhitneyu
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from PCME_functions import *

# Set up scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

print("=== Step 1.1: Data Loading and Quality Control ===")
print("Focus: FDZS A/C samples only")

# Set data path
DATA_PATH = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/analysis/cell2location_map"
SAVE_PATH = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/PCME"

if not os.path.exists(SAVE_PATH):
    os.makedirs(SAVE_PATH)

#%% Step 1.1: Load and Check data

# # Define FDZS A/C samples only
# fdzs_samples = {
#     # A group (Early Recurrence)
#     'FDZS_A04932G3_bin50_A1': {'group': 'A', 'recurrence': 'Early_Recurrence'},
#     'FDZS_A04932G3_bin50_A2': {'group': 'A', 'recurrence': 'Early_Recurrence'},
#     'FDZS_A04932G3_bin50_A3': {'group': 'A', 'recurrence': 'Early_Recurrence'},
#     'FDZS_A04932G3_bin50_A5': {'group': 'A', 'recurrence': 'Early_Recurrence'},
#     'FDZS_A04932G3_bin50_A6': {'group': 'A', 'recurrence': 'Early_Recurrence'},
#     'FDZS_A04932G3_bin50_A7': {'group': 'A', 'recurrence': 'Early_Recurrence'},
    
#     # C group (Non-Early Recurrence)
#     'FDZS_A04932G3_bin50_C5': {'group': 'C', 'recurrence': 'Non_Early_Recurrence'},
#     'FDZS_A04932G3_bin50_C6': {'group': 'C', 'recurrence': 'Non_Early_Recurrence'},
#     'FDZS_A04932G3_bin50_C7': {'group': 'C', 'recurrence': 'Non_Early_Recurrence'},
# }

# print(f"\nTarget samples: {len(fdzs_samples)} FDZS slides")
# print(f"A group (Early Recurrence): {len([s for s in fdzs_samples.values() if s['group'] == 'A'])} samples")
# print(f"C group (Non-Early Recurrence): {len([s for s in fdzs_samples.values() if s['group'] == 'C'])} samples")

# # Load each sample
# data_path = Path(DATA_PATH)
# adata_list = []

# print(f"\nLoading samples from: {data_path}")

# for sample_name, sample_info in fdzs_samples.items():
#     print(f"\nProcessing {sample_name}...")
    
#     # Find h5ad file in sample folder
#     sample_path = data_path / sample_name
#     h5ad_files = list(sample_path.glob("*.h5ad"))
    
#     if len(h5ad_files) == 0:
#         print(f"  ❌ No h5ad file found in {sample_path}")
#         continue
    
#     # Load the data
#     adata = sc.read_h5ad(h5ad_files[0])

#     # Add meta information
#     adata.obs["n_spots"] = adata.n_obs
#     adata.obs["n_genes"] = adata.n_vars
#     adata.obs["sample_id"] = adata.obs['sample'].astype(str).values + "_" + adata.obs['core_name'].astype(str).values

#     # Add to list
#     adata_list.append(adata)
#     print(f"  ✅ Loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")

# # Check loading results
# print(f"Successfully loaded: {len(adata_list)} / {len(fdzs_samples)} samples")

# # Combine all samples
# print(f"\nCombining {len(adata_list)} samples...")
# adata_combined = sc.concat(
#     adata_list, 
#     join='outer', 
#     label='batch_sample',
#     keys=[adata.obs['core_name'].iloc[0] for adata in adata_list]
# )

adata_combined = sc.read_h5ad(f"{DATA_PATH}/FDZS_A04932G3_bin100/sp.h5ad")
print(f"Combined dataset: {adata_combined.n_obs:,} spots, {adata_combined.n_vars:,} genes")

# Create metadata DataFrame
metadata_df = pd.DataFrame(adata_combined.obs)

# Display sample distribution
print(f"\n📊 Sample Distribution:")
print(metadata_df.groupby('RFS_status').size())

# Extract cell abundance data from cell2location results
print(f"\n🔬 Extracting cell abundance data...")

# Check available abundance matrices
abundance_keys = [key for key in adata_combined.obsm.keys() if 'abundance' in key]
print(f"Available abundance matrices: {abundance_keys}")

# Use the main results (means)
if 'q05_cell_abundance_w_sf' in adata_combined.obsm.keys():
    abundance_data = adata_combined.obsm['q05_cell_abundance_w_sf'].copy()
    print(f"Using 'q05_cell_abundance_w_sf' with shape: {abundance_data.shape}")

# Clean cell type names
if abundance_data is not None:
    # Remove prefix from column names
    new_columns = []
    for col in abundance_data.columns:
        if col.startswith('q05_cell_abundance_w_sf'):
            clean_name = col.replace('q05_cell_abundance_w_sf', '')
            new_columns.append(clean_name)
        else:
            new_columns.append(col)
    
    abundance_data.columns = new_columns
    
    # Add sample information
    abundance_df = abundance_data.copy()
    abundance_df['sample_id'] = adata_combined.obs['sample'].astype(str).values + "_" + adata_combined.obs['core_name'].astype(str).values
    abundance_df['RFS_status'] = adata_combined.obs['RFS_status'].values
    
    # Get cell type columns
    cell_type_cols = [col for col in abundance_df.columns 
                     if col not in ['sample_id', 'RFS_status']]
    
    print(f"Extracted {len(cell_type_cols)} cell types")

# Quality Control Summary
print(f"\n📈 Quality Control Summary:")
print(f"Total spots: {adata_combined.n_obs:,}")
print(f"Total genes: {adata_combined.n_vars:,}")
print(f"Total samples: {len(metadata_df)}")

# Save results
print(f"\n💾 Saving results...")
adata_combined.write(os.path.join(SAVE_PATH,"step1_fdzs_combined_data.h5ad")) # Save main data
metadata_df.to_csv(os.path.join(SAVE_PATH,"step1_fdzs_metadata.csv"), index=False) # Save metadata
abundance_df.to_csv(os.path.join(SAVE_PATH,"step1_fdzs_abundance.csv"), index=False) # Save abundance data

#%% Step 1.2: Overall description and visualization data
from matplotlib.patches import Patch

adata = sc.read_h5ad(f"{SAVE_PATH}/step1_fdzs_combined_data.h5ad")
metadata_df = pd.read_csv(f"{SAVE_PATH}/step1_fdzs_metadata.csv")
abundance_df = pd.read_csv(f"{SAVE_PATH}/step1_fdzs_abundance.csv")

# Create visualizations
print("Creating sample overview plots...")
fig1= plot_sample_overview(metadata_df, figsize=(15, 9))
fig1.savefig(os.path.join(SAVE_PATH,'step1_sample_overview.pdf'), dpi=300, bbox_inches='tight')
plt.show()

print("Creating cell type overview plots...")
fig2, _ = plot_cell_type_overview(abundance_df, figsize=(24, 20))
fig2.savefig(os.path.join(SAVE_PATH,'step1_cell_type_overview.pdf'), dpi=300, bbox_inches='tight')
plt.show()

#%% Step 2.1: Cholangiocyte Spot Identification
# Load the combined data
SAVE_PATH = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/PCME"
adata = sc.read_h5ad(f"{SAVE_PATH}/step1_fdzs_combined_data.h5ad")
metadata_df = pd.read_csv(f"{SAVE_PATH}/step1_fdzs_metadata.csv")
abundance_df = pd.read_csv(f"{SAVE_PATH}/step1_fdzs_abundance.csv")

print(f"Loaded data: {adata.n_obs:,} spots, {adata.n_vars:,} genes")
print(f"RFS status distribution: {adata.obs['RFS_status'].value_counts().to_dict()}")

# Define cholangiocyte signature genes based on literature
cholangiocyte_genes = [
    'KRT7',      # Cytokeratin 7 - classic cholangiocyte marker
    'KRT19',     # Cytokeratin 19 - classic cholangiocyte marker
    'EPCAM',     # Epithelial cell adhesion molecule
    'SOX9',      # SRY-box transcription factor 9
    'AQP1',      # Aquaporin 1
]

print(f"\n🔍 Searching for cholangiocyte signature genes...")
print(f"Target genes: {cholangiocyte_genes}")

# Check which genes are available in the data
available_genes = []

for gene in cholangiocyte_genes:
    if gene in adata.var_names:
        available_genes.append(gene)
        print(f"  ✅ {gene} - Found")
    else:
        print(f"  ❌ {gene} - Not found")

# Calculate cholangiocyte signature score using available genes
print(f"\n📊 Calculating cholangiocyte signature score...")

# Get gene expression matrix for available genes
gene_subset = adata[:, available_genes].X.toarray()
gene_names = adata[:, available_genes].var_names.tolist()

# Calculate mean expression across available genes
cholangiocyte_signature = np.mean(gene_subset, axis=1)

# Add signature score to adata
adata.obs['cholangiocyte_signature'] = cholangiocyte_signature

print(f"Signature score range: {cholangiocyte_signature.min():.3f} - {cholangiocyte_signature.max():.3f}")
print(f"Signature score mean: {cholangiocyte_signature.mean():.3f}")

# Get epithelial cell abundance from cell2location results
epithelial_abundance = abundance_df['q05cell_abundance_w_sf_Epithelial'].values
adata.obs['epithelial_abundance'] = epithelial_abundance

print(f"\n🧬 Epithelial cell abundance:")
print(f"Range: {epithelial_abundance.min():.3f} - {epithelial_abundance.max():.3f}")
print(f"Mean: {epithelial_abundance.mean():.3f}")

# Define cholangiocyte-enriched spots using combined criteria
print(f"\n🎯 Defining cholangiocyte-enriched spots...")

## cholangiocyte signature basaed
signature_threshold = np.percentile(cholangiocyte_signature, 95) # 1.0 # Top 99%
high_signature_spots = cholangiocyte_signature > signature_threshold

## epithelial fraction signature basaed
# fraction_threshold = np.percentile(epithelial_abundance, 99) # 1.0 # Top 99%
# high_signature_spots = epithelial_abundance > fraction_threshold

adata.obs['cholangiocyte_enriched'] = high_signature_spots
print(f"  High signature spots: {high_signature_spots.sum():,} ({high_signature_spots.mean()*100:.1f}%)")

# Analyze distribution by RFS status
print(f"\n📈 Cholangiocyte spots by RFS status:")
cholangiocyte_by_rfs = pd.crosstab(
    adata.obs['RFS_status'], 
    adata.obs['cholangiocyte_enriched'], 
    normalize='index'
) * 100
print(cholangiocyte_by_rfs.round(1))

# Sample-level analysis
adata.obs["sample_id"] = [str(x) for x in adata.obs["core_name"]]
sample_cholangiocyte_stats = []

for sample_id in adata.obs['sample_id'].unique():
    sample_mask = adata.obs['sample_id'] == sample_id
    sample_data = adata.obs[sample_mask]
    
    total_spots = len(sample_data)
    cholangiocyte_spots = sample_data['cholangiocyte_enriched'].sum()
    cholangiocyte_percentage = (cholangiocyte_spots / total_spots) * 100
    rfs_status = sample_data['RFS_status'].iloc[0]
    
    sample_stats = {
        'sample_id': sample_id,
        'RFS_status': rfs_status,
        'total_spots': total_spots,
        'cholangiocyte_spots': cholangiocyte_spots,
        'cholangiocyte_percentage': cholangiocyte_percentage,
        'mean_signature': sample_data['cholangiocyte_signature'].mean(),
        'mean_epithelial': sample_data['epithelial_abundance'].mean()
    }
    sample_cholangiocyte_stats.append(sample_stats)

sample_stats_df = pd.DataFrame(sample_cholangiocyte_stats)
print(sample_stats_df.round(2))

# Statistical comparison between RFS groups at sample level
rfs0_percentages = sample_stats_df[sample_stats_df['RFS_status'] == 0]['cholangiocyte_percentage']
rfs1_percentages = sample_stats_df[sample_stats_df['RFS_status'] == 1]['cholangiocyte_percentage']

if len(rfs0_percentages) > 0 and len(rfs1_percentages) > 0:
    stat, p_value = stats.mannwhitneyu(rfs1_percentages, rfs0_percentages, alternative='two-sided')
    print(f"\nStatistical comparison (sample-level):")
    print(f"RFS_status 0: {rfs0_percentages.mean():.1f}% ± {rfs0_percentages.std():.1f}%")
    print(f"RFS_status 1: {rfs1_percentages.mean():.1f}% ± {rfs1_percentages.std():.1f}%") 
    print(f"Mann-Whitney U test p-value: {p_value:.4f}")

# Create comprehensive visualization
fig = plt.figure(figsize=(15, 8))

# 1. RFS comparison - cholangiocyte percentage
ax1 = plt.subplot(1, 3, 1)
rfs_groups = [
    sample_stats_df[sample_stats_df['RFS_status'] == 0]['cholangiocyte_percentage'],
    sample_stats_df[sample_stats_df['RFS_status'] == 1]['cholangiocyte_percentage']
]
bp = ax1.boxplot(rfs_groups, labels=['RFS 0', 'RFS 1'], patch_artist=True)
bp['boxes'][0].set_facecolor('lightblue')
bp['boxes'][1].set_facecolor('lightcoral')
ax1.set_ylabel('Cholangiocyte Spots (%)')
ax1.set_title('Cholangiocyte % by RFS Status')

# Add individual points
for i, group_data in enumerate(rfs_groups):
    x = np.random.normal(i+1, 0.04, size=len(group_data))
    ax1.scatter(x, group_data, alpha=0.7, s=50)

# 2. Cholangiocyte spots per sample
ax2 = plt.subplot(1, 3, 2)
x_pos = range(len(sample_stats_df))
bars = ax2.bar(x_pos, sample_stats_df['cholangiocyte_percentage'], 
               color=['red' if rfs == 1 else 'blue' for rfs in sample_stats_df['RFS_status']])
ax2.set_xlabel('Sample')
ax2.set_ylabel('Cholangiocyte Spots (%)')
ax2.set_title('Cholangiocyte % per Sample')
ax2.set_xticks(x_pos)
ax2.set_xticklabels(sample_stats_df['sample_id'], rotation=45)

# Add value labels
for bar, val in zip(bars, sample_stats_df['cholangiocyte_percentage']):
    ax2.text(bar.get_x() + bar.get_width()/2., val,
             f'{val:.1f}%', ha='center', va='bottom', fontsize=8)

# 3. Validation: Gene expression in cholangiocyte spots
ax3 = plt.subplot(1, 3, 3)

if len(available_genes) > 0:
    # Compare gene expression in cholangiocyte vs non-cholangiocyte spots
    chol_spots = adata.obs['cholangiocyte_enriched']
    
    expression_comparison = []
    for gene in available_genes:  # Show top 4 genes
        chol_expr = adata[chol_spots, gene].X.toarray().flatten()
        non_chol_expr = adata[~chol_spots, gene].X.toarray().flatten()
        
        expression_comparison.append({
            'Gene': gene,
            'Cholangiocyte_mean': chol_expr.mean(),
            'Non_cholangiocyte_mean': non_chol_expr.mean(),
            'Fold_change': chol_expr.mean() / (non_chol_expr.mean() + 1e-10)
        })
    
    expr_df = pd.DataFrame(expression_comparison)
    
    x = range(len(expr_df))
    width = 0.35
    
    ax3.bar([i - width/2 for i in x], expr_df['Cholangiocyte_mean'], 
             width, label='Cholangiocyte spots', alpha=0.8)
    ax3.bar([i + width/2 for i in x], expr_df['Non_cholangiocyte_mean'], 
             width, label='Other spots', alpha=0.8)
    
    ax3.set_xlabel('Genes')
    ax3.set_ylabel('Mean Expression')
    ax3.set_title('Gene Expression Validation')
    ax3.set_xticks(x)
    ax3.set_xticklabels(expr_df['Gene'], rotation=45)
    ax3.legend()
    
    # Print fold changes
    print(f"\nGene expression validation:")
    for _, row in expr_df.iterrows():
        print(f"  {row['Gene']}: {row['Fold_change']:.2f}x higher in cholangiocyte spots")

plt.tight_layout()
plt.savefig(f'{SAVE_PATH}/step2_cholangiocyte_identification.pdf', dpi=300, bbox_inches='tight')
plt.show()

# 4. Draw spatial distribution of Score
fig = plt.figure(figsize=(20, 20))

# Get sample data
sample_ids = adata.obs['sample_id'].unique()
for i, sample_id_ in enumerate(sample_ids):
    ax = plt.subplot(3, 3, i + 1)
    sample_mask = adata.obs['sample_id'] == sample_id_
    sample_adata = adata[sample_mask]
    
    # Get spatial coordinates
    if 'spatial' in sample_adata.obsm.keys():
        coords = sample_adata.obsm['spatial']
    elif 'X_spatial' in sample_adata.obsm.keys():
        coords = sample_adata.obsm['X_spatial']  
    else:
        coords = np.column_stack([sample_adata.obs['x'].values, sample_adata.obs['y'].values])
    
    # Plot background (non-enriched) spots first
    non_enriched_mask = sample_adata.obs['cholangiocyte_enriched'] == 0
    enriched_mask = sample_adata.obs['cholangiocyte_enriched'] == 1
    
    # Background spots in light gray
    ax.scatter(coords[non_enriched_mask, 0], coords[non_enriched_mask, 1], 
              c='lightgray', s=1, alpha=0.5, label='Non-enriched')
    
    # Highlight enriched spots in bright red
    if enriched_mask.sum() > 0:
        ax.scatter(coords[enriched_mask, 0], coords[enriched_mask, 1], 
                  c='red', s=10, alpha=0.9, marker='o', 
                  edgecolors='darkred', linewidth=0.5, label='Cholangiocyte enriched')
    
    rfs_status = sample_adata.obs['RFS_status'].iloc[0]
    cholangiocyte_pct = (sample_adata.obs['cholangiocyte_enriched'].sum() / len(sample_adata)) * 100
    cholangiocyte_count = sample_adata.obs['cholangiocyte_enriched'].sum()
    
    ax.set_title(f'{sample_id_}\nRFS: {rfs_status}, Chol: {cholangiocyte_count} spots ({cholangiocyte_pct:.1f}%)')
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    
    # Add legend only to first subplot
    if i == 0:
        ax.legend(markerscale=2, fontsize=8)

plt.tight_layout()
plt.savefig(f'{SAVE_PATH}/step2_cholangiocyte_spatial_distribution.pdf', dpi=300, bbox_inches='tight')
plt.show()

#%% Step 2.2: CAIX Expression Analysis in Cholangiocytes

# Check for CA9 (CAIX) gene availability
caix_gene = 'CA9'  # Carbonic anhydrase IX gene symbol
print(f"\n🔍 Searching for CAIX gene (CA9)...")

## Use raw CAIX gene
# if caix_gene in adata.var_names:
#     print(f"  ✅ {caix_gene} found in the data")
#     caix_expression = adata[:, caix_gene].X.toarray().flatten()
#     adata.obs['CA9_expression'] = caix_expression

## Use hypoxic signature
hypoxia_genes = ['CA9', 'HIF1A',] #  'VEGFA', 'LDHA', 'PKM', 'SLC2A1', 'BNIP3'
hypoxia_expr = adata[:, hypoxia_genes].X.toarray()
caix_expression = np.mean(hypoxia_expr, axis=1)
adata.obs['CA9_expression'] = caix_expression

# Focus analysis on cholangiocyte-enriched spots only
cholangiocyte_spots = adata.obs['cholangiocyte_enriched'] == True
cholangiocyte_adata = adata[cholangiocyte_spots].copy()

print(f"\n🎯 Focusing on cholangiocyte-enriched spots...")
print(f"Cholangiocyte spots for analysis: {cholangiocyte_adata.n_obs:,}")
print(f"RFS distribution in cholangiocytes: {cholangiocyte_adata.obs['RFS_status'].value_counts().to_dict()}")

# CAIX expression analysis in cholangiocytes
cholangiocyte_caix = cholangiocyte_adata.obs['CA9_expression'].values

print(f"\nCAIX in cholangiocytes:")
print(f"  Range: {cholangiocyte_caix.min():.3f} - {cholangiocyte_caix.max():.3f}")
print(f"  Mean: {cholangiocyte_caix.mean():.3f}")

# Define CAIX+ vs CAIX- cholangiocytes
# caix_threshold = np.percentile(cholangiocyte_caix, 90)  # Top 10%
caix_threshold = 0.3

cholangiocyte_adata.obs['CAIX_positive'] = cholangiocyte_adata.obs['CA9_expression'] > caix_threshold

caix_positive_count = cholangiocyte_adata.obs['CAIX_positive'].sum()
print(f"\nCAIX+ cholangiocytes: {caix_positive_count} ({caix_positive_count/len(cholangiocyte_adata)*100:.1f}%)")
print(f"CAIX threshold: {caix_threshold:.3f}")

# Compare CAIX expression between RFS groups in cholangiocytes
print(f"\n📊 CAIX Expression: RFS Group Comparison (Cholangiocytes Only)")

rfs0_cholangiocytes = cholangiocyte_adata[cholangiocyte_adata.obs['RFS_status'] == 0]
rfs1_cholangiocytes = cholangiocyte_adata[cholangiocyte_adata.obs['RFS_status'] == 1]

if len(rfs0_cholangiocytes) > 0 and len(rfs1_cholangiocytes) > 0:
    rfs0_caix = rfs0_cholangiocytes.obs['CA9_expression'].values
    rfs1_caix = rfs1_cholangiocytes.obs['CA9_expression'].values
    
    print(f"RFS 0 cholangiocytes: n={len(rfs0_caix)}, CAIX = {rfs0_caix.mean():.3f} ± {rfs0_caix.std():.3f}")
    print(f"RFS 1 cholangiocytes: n={len(rfs1_caix)}, CAIX = {rfs1_caix.mean():.3f} ± {rfs1_caix.std():.3f}")
    
    # Statistical test
    stat, p_value_caix = stats.mannwhitneyu(rfs1_caix, rfs0_caix, alternative='two-sided')
    print(f"Mann-Whitney U test p-value: {p_value_caix:.4f}")
    
    # Effect size (fold change)
    fold_change = rfs1_caix.mean() / (rfs0_caix.mean() + 1e-10)
    print(f"Fold change (RFS 1 vs RFS 0): {fold_change:.2f}")
    
else:
    print("❌ Insufficient data for RFS group comparison")
    p_value_caix = 1.0
    fold_change = 1.0

# Sample-level CAIX analysis in cholangiocytes
print(f"\n📈 Sample-level CAIX Analysis (Cholangiocytes)")

sample_caix_stats = []
for sample_id in cholangiocyte_adata.obs['sample_id'].unique():
    sample_cholangiocytes = cholangiocyte_adata[cholangiocyte_adata.obs['sample_id'] == sample_id]
    
    if len(sample_cholangiocytes) > 0:
        sample_stats = {
            'sample_id': sample_id,
            'RFS_status': sample_cholangiocytes.obs['RFS_status'].iloc[0],
            'n_cholangiocytes': len(sample_cholangiocytes),
            'mean_CAIX': sample_cholangiocytes.obs['CA9_expression'].mean(),
            'median_CAIX': sample_cholangiocytes.obs['CA9_expression'].median(),
            'CAIX_positive_count': sample_cholangiocytes.obs['CAIX_positive'].sum(),
            'CAIX_positive_percentage': (sample_cholangiocytes.obs['CAIX_positive'].sum() / len(sample_cholangiocytes)) * 100
        }
        sample_caix_stats.append(sample_stats)

sample_caix_df = pd.DataFrame(sample_caix_stats)
print("\nSample-level CAIX statistics:")
print(sample_caix_df.round(3))

# Statistical test at sample level
if len(sample_caix_df) > 0:
    rfs0_samples = sample_caix_df[sample_caix_df['RFS_status'] == 0]
    rfs1_samples = sample_caix_df[sample_caix_df['RFS_status'] == 1]
    
    if len(rfs0_samples) > 0 and len(rfs1_samples) > 0:
        # Test mean CAIX
        stat_sample, p_value_sample = stats.mannwhitneyu(
            rfs1_samples['mean_CAIX'], rfs0_samples['mean_CAIX'], alternative='two-sided'
        )
        
        # Test CAIX+ percentage
        stat_pct, p_value_pct = stats.mannwhitneyu(
            rfs1_samples['CAIX_positive_percentage'], rfs0_samples['CAIX_positive_percentage'], alternative='two-sided'
        )
        
        print(f"\nSample-level statistical tests:")
        print(f"Mean CAIX (RFS 0 vs 1): p = {p_value_sample:.4f}")
        print(f"CAIX+ percentage (RFS 0 vs 1): p = {p_value_pct:.4f}")
        
        print(f"\nSample-level summary:")
        print(f"RFS 0 samples: {rfs0_samples['mean_CAIX'].mean():.3f} ± {rfs0_samples['mean_CAIX'].std():.3f}")
        print(f"RFS 1 samples: {rfs1_samples['mean_CAIX'].mean():.3f} ± {rfs1_samples['mean_CAIX'].std():.3f}")

# CAIX+ cholangiocyte distribution analysis
print(f"\n🔍 CAIX+ Cholangiocyte Distribution by RFS Status")

caix_rfs_crosstab = pd.crosstab(
    cholangiocyte_adata.obs['RFS_status'], 
    cholangiocyte_adata.obs['CAIX_positive'], 
    normalize='index'
) * 100

print("Percentage of CAIX+ cholangiocytes per RFS group:")
print(caix_rfs_crosstab.round(1))

# Create comprehensive visualization
fig = plt.figure(figsize=(16, 6))

# 1. CAIX expression by RFS status (cholangiocytes only)
ax1 = plt.subplot(1, 4, 1)
if len(rfs0_cholangiocytes) > 0 and len(rfs1_cholangiocytes) > 0:
    bp = ax1.boxplot([rfs0_caix, rfs1_caix], labels=['RFS 0', 'RFS 1'], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    
    # Add individual points
    for i, data in enumerate([rfs0_caix, rfs1_caix]):
        x = np.random.normal(i+1, 0.04, size=len(data))
        ax1.scatter(x, data, alpha=0.6, s=10)
    
    ax1.set_ylabel('CAIX Expression')
    ax2.set_title(f'CAIX in Cholangiocytes\n(p = {p_value_caix:.4f})')
    
    # Add statistics text
    ax1.text(0.02, 0.98, f'Fold change: {fold_change:.2f}', transform=ax1.transAxes, 
             verticalalignment='top', fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# 2. CAIX+ percentage by sample
ax2 = plt.subplot(1, 4, 2)
if len(sample_caix_df) > 0:
    x_pos = range(len(sample_caix_df))
    bars = ax2.bar(x_pos, sample_caix_df['CAIX_positive_percentage'],
                   color=['red' if rfs == 1 else 'blue' for rfs in sample_caix_df['RFS_status']])
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('CAIX+ Cholangiocytes (%)')
    ax2.set_title('CAIX+ Percentage per Sample')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(sample_caix_df['sample_id'], rotation=45)
    
    # Add value labels
    for bar, val in zip(bars, sample_caix_df['CAIX_positive_percentage']):
        ax2.text(bar.get_x() + bar.get_width()/2., val,
                 f'{val:.1f}%', ha='center', va='bottom', fontsize=8)

# 3. CAIX+ vs CAIX- cholangiocyte comparison across samples
ax3 = plt.subplot(1, 4, 3)
caix_comparison_data = []

for _, row in sample_caix_df.iterrows():
    caix_comparison_data.append({
        'Sample': row['sample_id'],
        'RFS_Status': f"RFS_{row['RFS_status']}",
        'CAIX_Negative': 100 - row['CAIX_positive_percentage'],
        'CAIX_Positive': row['CAIX_positive_percentage']
    })

if caix_comparison_data:
    comp_df = pd.DataFrame(caix_comparison_data)
    
    x = range(len(comp_df))
    width = 0.8
    
    ax3.bar(x, comp_df['CAIX_Negative'], width, label='CAIX-', color='lightblue', alpha=0.8)
    ax3.bar(x, comp_df['CAIX_Positive'], width, bottom=comp_df['CAIX_Negative'], 
            label='CAIX+', color='red', alpha=0.8)
    
    ax3.set_xlabel('Sample')
    ax3.set_ylabel('Cholangiocyte Percentage')
    ax3.set_title('CAIX+ vs CAIX- Distribution')
    ax3.set_xticks(x)
    ax3.set_xticklabels(comp_df['Sample'], rotation=45)
    ax3.legend()

# 4. Sample-level CAIX statistics by RFS
ax4 = plt.subplot(1, 4, 4)
if len(rfs0_samples) > 0 and len(rfs1_samples) > 0:
    rfs_groups_sample = [rfs0_samples['mean_CAIX'], rfs1_samples['mean_CAIX']]
    bp = ax4.boxplot(rfs_groups_sample, labels=['RFS 0', 'RFS 1'], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    
    # Add individual points
    for i, data in enumerate(rfs_groups_sample):
        x = np.random.normal(i+1, 0.04, size=len(data))
        ax4.scatter(x, data, alpha=0.7, s=50)
    
    ax4.set_ylabel('Mean CAIX (Sample Level)')
    ax4.set_title(f'Sample-Level CAIX\n(p = {p_value_sample:.4f})')


plt.tight_layout()
plt.savefig(f'{SAVE_PATH}/step2_2_CAIX_analysis.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Spatial distribution of CAIX in cholangiocytes (show all samples)
sample_ids = cholangiocyte_adata.obs['sample_id'].unique()

fig = plt.figure(figsize=(20, 20))

for i, sample_id in enumerate(sample_ids):
    ax = plt.subplot(3, 3, 1 + i)
    
    # Get all spots for this sample (not just cholangiocytes)
    sample_all_spots = adata[adata.obs['sample_id'] == sample_id]
    
    # Get coordinates
    if 'spatial' in sample_all_spots.obsm.keys():
        coords = sample_all_spots.obsm['spatial']
    else:
        coords = np.column_stack([sample_all_spots.obs['x'].values, sample_all_spots.obs['y'].values])
    
    # Plot background spots
    ax.scatter(coords[:, 0], coords[:, 1], c='lightgray', s=1, alpha=0.3)
    
    # Overlay cholangiocytes colored by CAIX expression
    sample_cholangiocytes = cholangiocyte_adata[cholangiocyte_adata.obs['sample_id'] == sample_id]
    if len(sample_cholangiocytes) > 0:
        # Get coordinates for cholangiocytes
        chol_indices = sample_cholangiocytes.obs.index
        chol_coords_mask = np.isin(sample_all_spots.obs.index, chol_indices)
        chol_coords = coords[chol_coords_mask]
        
        scatter = ax.scatter(chol_coords[:, 0], chol_coords[:, 1], 
                           c=sample_cholangiocytes.obs['CA9_expression'], 
                           s=30, alpha=0.8, cmap='Reds', 
                           edgecolors='black', linewidth=0.5)
        
        # Add colorbar
        if i == 0:  # Only add colorbar to first plot
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label('CAIX Expression')
    
    rfs_status = sample_all_spots.obs['RFS_status'].iloc[0]
    n_cholangiocytes = len(sample_cholangiocytes)
    mean_caix = sample_cholangiocytes.obs['CA9_expression'].mean() if len(sample_cholangiocytes) > 0 else 0
    
    ax.set_title(f'{sample_id}\nRFS: {rfs_status}, n={n_cholangiocytes}\nCAIX: {mean_caix:.3f}')
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')

plt.tight_layout()
plt.savefig(f'{SAVE_PATH}/step2_2_CAIX_spatial_distribution.pdf', dpi=300, bbox_inches='tight')
plt.show()

#%% Step 2.3: Cholangiocyte Differential Gene Expression Analysis

# Focus on cholangiocyte-enriched spots
cholangiocyte_spots = adata.obs['cholangiocyte_enriched'] == True
cholangiocyte_adata = adata[cholangiocyte_spots].copy()

print(f"\n🎯 Focusing on cholangiocyte-enriched spots for DEG analysis...")
print(f"Cholangiocyte spots: {cholangiocyte_adata.n_obs:,}")
print(f"RFS distribution: {cholangiocyte_adata.obs['RFS_status'].value_counts().to_dict()}")

# Check sample distribution
sample_distribution = cholangiocyte_adata.obs.groupby(['sample_id', 'RFS_status']).size().reset_index(name='count')
print(f"\nSample distribution in cholangiocytes:")
print(sample_distribution)

# Prepare data for differential expression analysis
print(f"\n📊 Preparing data for differential expression analysis...")

# Set up comparison groups
rfs0_mask = cholangiocyte_adata.obs['RFS_status'] == 0
rfs1_mask = cholangiocyte_adata.obs['RFS_status'] == 1

n_rfs0 = rfs0_mask.sum()
n_rfs1 = rfs1_mask.sum()

print(f"RFS 0 (Non-early recurrence): {n_rfs0:,} cholangiocyte spots")
print(f"RFS 1 (Early recurrence): {n_rfs1:,} cholangiocyte spots")

if n_rfs0 < 10 or n_rfs1 < 10:
    print("⚠️  Warning: Low sample size for differential expression analysis")

# Perform differential expression analysis using scanpy
print(f"\n🧬 Performing differential expression analysis...")

# Add group labels for scanpy
cholangiocyte_adata.obs['RFS_group'] = cholangiocyte_adata.obs['RFS_status'].astype(str)

# Use scanpy's rank_genes_groups for differential expression
sc.tl.rank_genes_groups(
    cholangiocyte_adata, 
    'RFS_group', 
    groups=['1'],  # Compare RFS 1 vs rest (RFS 0)
    reference='0',
    method='wilcoxon',
    use_raw=False,
    n_genes=None  # Get all genes
)

# Extract results
deg_results = sc.get.rank_genes_groups_df(cholangiocyte_adata, group='1')
deg_results = deg_results.sort_values('pvals_adj')

print(f"Differential expression analysis complete:")
print(f"  Total genes tested: {len(deg_results)}")
print(f"  Significant genes (padj < 0.05): {(deg_results['pvals_adj'] <= 0.05).sum()}")
print(f"  Highly significant (padj < 0.01): {(deg_results['pvals_adj'] <= 0.01).sum()}")

# Manual calculation for additional statistics
print(f"\n📈 Calculating additional statistics...")
deg_stats = []

for gene in adata.var_names:
    if gene in cholangiocyte_adata.var_names:
        rfs0_expr = cholangiocyte_adata[rfs0_mask, gene].X.toarray().flatten()
        rfs1_expr = cholangiocyte_adata[rfs1_mask, gene].X.toarray().flatten()
        
        # Basic statistics
        mean_rfs0 = np.mean(rfs0_expr)
        mean_rfs1 = np.mean(rfs1_expr)
        
        # Calculate fold change
        if mean_rfs0 > 0:
            fold_change = mean_rfs1 / mean_rfs0
        else:
            fold_change = np.inf if mean_rfs1 > 0 else 1.0
        
        log2_fc = np.log2(fold_change) if fold_change > 0 and fold_change != np.inf else 0
        
        # Expression percentages
        pct_rfs0 = (rfs0_expr > 0).mean() * 100
        pct_rfs1 = (rfs1_expr > 0).mean() * 100
        
        deg_stats.append({
            'gene': gene,
            'mean_rfs0': mean_rfs0,
            'mean_rfs1': mean_rfs1,
            'fold_change': fold_change,
            'log2_fc': log2_fc,
            'pct_rfs0': pct_rfs0,
            'pct_rfs1': pct_rfs1
        })

deg_stats_df = pd.DataFrame(deg_stats)

# Merge scanpy results with manual calculations
deg_combined = deg_results.merge(deg_stats_df, left_on='names', right_on='gene', how='left')
deg_combined = deg_combined.sort_values('pvals_adj')

print(f"Combined DEG analysis complete")
deg_combined.to_csv(f"{SAVE_PATH}/step2_3_cholangiocyte_DEG_results.csv", index=False) # Save DEG results

# Filter for significant genes
significant_genes = deg_combined[deg_combined['pvals'] <= 0.05].copy()
print(f"\nSignificant genes (padj < 0.05): {len(significant_genes)}")

if len(significant_genes) > 0:
    print(f"\nTop 10 most significant upregulated genes (RFS 1 > RFS 0):")
    upregulated = significant_genes[significant_genes['log2_fc'] > 0].head(10)
    for _, gene in upregulated.iterrows():
        print(f"  {gene['names']}: log2FC={gene['log2_fc']:.2f}, padj={gene['pvals_adj']:.2e}")
    
    print(f"\nTop 10 most significant downregulated genes (RFS 1 < RFS 0):")
    downregulated = significant_genes[significant_genes['log2_fc'] < 0].head(10)
    for _, gene in downregulated.iterrows():
        print(f"  {gene['names']}: log2FC={gene['log2_fc']:.2f}, padj={gene['pvals_adj']:.2e}")

# Focus on PCME-relevant gene categories
print(f"\n🔍 Analyzing PCME-relevant gene categories...")

# Define gene categories of interest
gene_categories = load_pathway_databases(database_name = "KEGG")

# Get significant upregulated and downregulated genes separately
significant_up = significant_genes[significant_genes['log2_fc'] > 0]['names'].tolist()
significant_down = significant_genes[significant_genes['log2_fc'] < 0]['names'].tolist()
all_tested_genes = deg_combined['names'].tolist()

print(f"Performing enrichment analysis:")
print(f"  Upregulated genes: {len(significant_up)}")
print(f"  Downregulated genes: {len(significant_down)}")
print(f"  Background genes: {len(all_tested_genes)}")

# Perform enrichment for upregulated genes
if len(significant_up) > 0:
    enrichment_up = perform_pathway_enrichment(
        significant_up, gene_categories, all_tested_genes, p_threshold=0.05
    )
    print(f"  Significantly enriched pathways (upregulated): {len(enrichment_up)}")
else:
    enrichment_up = pd.DataFrame()

# Perform enrichment for downregulated genes
if len(significant_down) > 0:
    enrichment_down = perform_pathway_enrichment(
        significant_down, gene_categories, all_tested_genes, p_threshold=0.05
    )
    print(f"  Significantly enriched pathways (downregulated): {len(enrichment_down)}")
else:
    enrichment_down = pd.DataFrame()

# Print top enriched pathways
if len(enrichment_up) > 0:
    print(f"\nTop upregulated pathways:")
    for _, pathway in enrichment_up.head(5).iterrows():
        print(f"  {pathway['pathway']}: {pathway['overlap_size']}/{pathway['pathway_size']} genes, "
              f"p_adj={pathway['p_adjusted']:.2e}, ratio={pathway['enrichment_ratio']:.2f}")

if len(enrichment_down) > 0:
    print(f"\nTop downregulated pathways:")
    for _, pathway in enrichment_down.head(5).iterrows():
        print(f"  {pathway['pathway']}: {pathway['overlap_size']}/{pathway['pathway_size']} genes, "
              f"p_adj={pathway['p_adjusted']:.2e}, ratio={pathway['enrichment_ratio']:.2f}")

# Create comprehensive visualization
fig = plt.figure(figsize=(16, 8))

# 1. Volcano plot
ax1 = plt.subplot(1, 4, 1)
x_vals = deg_combined['log2_fc'].fillna(0)
y_vals = -np.log10(deg_combined['pvals_adj'].fillna(1))

# Color points based on significance and fold change
colors = []
for _, row in deg_combined.iterrows():
    if row['pvals_adj'] < 0.01 and abs(row['log2_fc']) > 1:
        colors.append('red')
    elif row['pvals_adj'] < 0.05 and abs(row['log2_fc']) > 0.5:
        colors.append('orange')
    elif row['pvals_adj'] < 0.05:
        colors.append('lightcoral')
    else:
        colors.append('gray')

ax1.scatter(x_vals, y_vals, c=colors, alpha=0.6, s=20)
ax1.set_xlabel('Log2 Fold Change (RFS 1 vs 0)')
ax1.set_ylabel('-Log10 Adjusted P-value')
ax1.set_title('Volcano Plot: Cholangiocyte DEGs')

# Add significance lines
ax1.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p=0.05')
ax1.axhline(y=-np.log10(0.01), color='red', linestyle='--', alpha=0.5, label='p=0.01')
ax1.axvline(x=1, color='black', linestyle='--', alpha=0.5)
ax1.axvline(x=-1, color='black', linestyle='--', alpha=0.5)
ax1.legend()

# Label top genes
if len(significant_genes) > 0:
    top_genes = significant_genes.head(10)
    for _, gene in top_genes.iterrows():
        if abs(gene['log2_fc']) > 0.5:  # Only label genes with substantial fold change
            ax1.annotate(gene['names'], 
                        (gene['log2_fc'], -np.log10(gene['pvals'])),
                        xytext=(5, 5), textcoords='offset points', 
                        fontsize=8, alpha=0.8)

# 2. Top upregulated genes
ax2 = plt.subplot(1, 4, 2)
if len(significant_genes) > 0:
    # Sort by log2_fc descending (most upregulated first)
    top_up = significant_genes[significant_genes['log2_fc'] > 0].sort_values(
        'log2_fc', ascending=False).head(15)
    
    if len(top_up) > 0:
        # Create gradient colors based on fold change
        norm = plt.Normalize(vmin=top_up['log2_fc'].min(), vmax=top_up['log2_fc'].max())
        colors = plt.cm.Reds(norm(top_up['log2_fc']))
        
        y_pos = range(len(top_up))
        bars = ax2.barh(y_pos, top_up['log2_fc'], color=colors, alpha=0.8, 
                        edgecolor='darkred', linewidth=0.5)
        
        # Reverse y-axis to show highest values at top
        ax2.invert_yaxis()
        ax2.set_yticks(y_pos)
        
        # Clean gene names (remove prefix if needed)
        clean_names = [name.replace('_', ' ').title() for name in top_up['names']]
        ax2.set_yticklabels(clean_names, fontsize=8, fontweight='bold')
        
        ax2.set_xlabel('Log₂ Fold Change', fontsize=10, fontweight='bold')
        ax2.set_title('Top Upregulated Genes\n(RFS 1 > RFS 0)', 
                        fontsize=11, fontweight='bold', color='darkred')
        
        # Add p-values with better positioning
        for i, (_, gene) in enumerate(top_up.iterrows()):
            p_val = gene["pvals"]
            if p_val < 0.001:
                p_text = '***'
            elif p_val < 0.01:
                p_text = '**'
            elif p_val < 0.05:
                p_text = '*'
            else:
                p_text = f'{p_val:.1e}'
            
            ax2.text(gene['log2_fc'] + max(top_up['log2_fc']) * 0.02, i, p_text, 
                    va='center', fontsize=7, fontweight='bold', 
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Add grid for better readability
        ax2.grid(axis='x', alpha=0.3, linestyle='--')
        ax2.set_axisbelow(True)

# 3. Top downregulated genes (Enhanced)
ax3 = plt.subplot(1, 4, 3)
if len(significant_genes) > 0:
    # Sort by log2_fc ascending (most downregulated first)
    top_down = significant_genes[significant_genes['log2_fc'] < 0].sort_values(
        'log2_fc', ascending=True).head(15)
    
    if len(top_down) > 0:
        # Create gradient colors based on fold change (blues for downregulated)
        norm = plt.Normalize(vmin=top_down['log2_fc'].min(), vmax=top_down['log2_fc'].max())
        colors = plt.cm.Blues_r(norm(top_down['log2_fc']))
        
        y_pos = range(len(top_down))
        bars = ax3.barh(y_pos, top_down['log2_fc'], color=colors, alpha=0.8,
                        edgecolor='darkblue', linewidth=0.5)
        
        # Reverse y-axis to show most downregulated at top
        ax3.invert_yaxis()
        ax3.set_yticks(y_pos)
        
        # Clean gene names
        clean_names = [name.replace('_', ' ').title() for name in top_down['names']]
        ax3.set_yticklabels(clean_names, fontsize=8, fontweight='bold')
        
        ax3.set_xlabel('Log₂ Fold Change', fontsize=10, fontweight='bold')
        ax3.set_title('Top Downregulated Genes\n(RFS 1 < RFS 0)', 
                        fontsize=11, fontweight='bold', color='darkblue')
        
        # Add p-values
        for i, (_, gene) in enumerate(top_down.iterrows()):
            p_val = gene["pvals"]
            if p_val < 0.001:
                p_text = '***'
            elif p_val < 0.01:
                p_text = '**'
            elif p_val < 0.05:
                p_text = '*'
            else:
                p_text = f'{p_val:.1e}'
            
            ax3.text(gene['log2_fc'] - abs(min(top_down['log2_fc'])) * 0.02, i, p_text, 
                    va='center', ha='right', fontsize=7, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
        
        # Add grid
        ax3.grid(axis='x', alpha=0.3, linestyle='--')
        ax3.set_axisbelow(True)

# 4. Enhanced GSEA Pathway Enrichment
ax4 = plt.subplot(1, 4, 4)
# Combine and sort enrichment results
combined_enrichment = []

if len(enrichment_up) > 0:
    enrichment_up_plot = enrichment_up.copy()
    enrichment_up_plot['direction'] = 'Upregulated'
    enrichment_up_plot['neg_log_p'] = -np.log10(enrichment_up_plot['p_adjusted'])
    combined_enrichment.append(enrichment_up_plot)

if len(enrichment_down) > 0:
    enrichment_down_plot = enrichment_down.copy()
    enrichment_down_plot['direction'] = 'Downregulated'
    enrichment_down_plot['neg_log_p'] = -np.log10(enrichment_down_plot['p_adjusted'])
    # Make enrichment ratio negative for visual separation
    enrichment_down_plot['enrichment_ratio'] = -enrichment_down_plot['enrichment_ratio']
    combined_enrichment.append(enrichment_down_plot)

if combined_enrichment:
    plot_enrichment = pd.concat(combined_enrichment, ignore_index=True)
    
    # Sort by significance (p_adjusted) then by enrichment ratio
    plot_enrichment = plot_enrichment.sort_values(
        ['p_adjusted', 'enrichment_ratio'], 
        ascending=[True, False]
    ).head(20)
    
    # Clean pathway names for better display
    def clean_pathway_name(pathway):
        # Remove common prefixes
        clean_name = pathway.replace('HALLMARK_', '').replace('KEGG_', '').replace('REACTOME_', '')
        clean_name = clean_name.replace('_', ' ').title()
        
        # Truncate long names smartly
        if len(clean_name) > 35:
            words = clean_name.split()
            truncated = []
            current_length = 0
            for word in words:
                if current_length + len(word) + 1 <= 35:
                    truncated.append(word)
                    current_length += len(word) + 1
                else:
                    break
            clean_name = ' '.join(truncated) + '...'
        
        return clean_name
    
    plot_enrichment['clean_pathway'] = plot_enrichment['pathway'].apply(clean_pathway_name)
    
    # Create the plot
    y_pos = range(len(plot_enrichment))
    
    # Enhanced color scheme
    colors = ['#E31A1C' if d == 'Upregulated' else '#1F78B4' for d in plot_enrichment['direction']]
    
    # Size based on overlap_size with better scaling
    sizes = np.clip(plot_enrichment['overlap_size'] * 8, 20, 200)
    
    # Create scatter plot
    scatter = ax4.scatter(plot_enrichment['enrichment_ratio'], y_pos, 
                            c=colors, s=sizes, alpha=0.8, 
                            edgecolors='black', linewidth=0.8)
    
    # Reverse y-axis to show most significant at top
    ax4.invert_yaxis()
    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(plot_enrichment['clean_pathway'], fontsize=8, fontweight='bold')
    
    ax4.set_xlabel('Enrichment Ratio', fontsize=10, fontweight='bold')
    ax4.set_title('GSEA Pathway Enrichment\n(Top 20 Significant Pathways)', 
                    fontsize=11, fontweight='bold')
    
    # Add reference lines
    ax4.axvline(x=0, color='gray', linestyle='-', alpha=0.5, linewidth=1)
    ax4.axvline(x=1, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
    ax4.axvline(x=-1, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
    
    # Enhanced legends
    # Direction legend
    legend_elements = [
        Patch(facecolor='#E31A1C', alpha=0.8, label='Upregulated in RFS 1'),
        Patch(facecolor='#1F78B4', alpha=0.8, label='Downregulated in RFS 1')
    ]
    direction_legend = ax4.legend(handles=legend_elements, loc='lower right', 
                                fontsize=8, framealpha=0.9)
    
    # Size legend with better positioning
    size_values = [5, 15, 30]
    size_legend_elements = [plt.scatter([], [], s=s*8, c='gray', alpha=0.8, 
                                        edgecolors='black', linewidth=0.8) 
                            for s in size_values]
    size_legend = ax4.legend(size_legend_elements, [str(s) for s in size_values], 
                            title='Gene Count', loc='lower left', fontsize=7, 
                            framealpha=0.9, title_fontsize=8)
    ax4.add_artist(size_legend)
    
    # Enhanced statistics box
    up_count = len(enrichment_up) if len(enrichment_up) > 0 else 0
    down_count = len(enrichment_down) if len(enrichment_down) > 0 else 0
    total_shown = len(plot_enrichment)
    
    stats_text = f"Showing: {total_shown}/{up_count + down_count} pathways\n" \
                f"Up: {up_count} | Down: {down_count}\n" \
                f"p < 0.05 (FDR adjusted)"
    
    ax4.text(0.02, 0.98, stats_text, transform=ax4.transAxes, 
            verticalalignment='top', fontsize=8, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', 
                        alpha=0.9, edgecolor='black', linewidth=0.5))
    
    # Add grid
    ax4.grid(axis='x', alpha=0.3, linestyle='--')
    ax4.set_axisbelow(True)
    
plt.tight_layout()
plt.subplots_adjust(wspace=0.3)
plt.savefig(f'{SAVE_PATH}/step2_3_cholangiocyte_DEG_analysis.pdf', dpi=300, bbox_inches='tight')
plt.show()

#%% Step 3.1: Spatial Distance Analysis Framework
# Load data and reconstruct cholangiocyte annotations
SAVE_PATH = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/PCME"
adata = sc.read_h5ad(f"{SAVE_PATH}/step1_fdzs_combined_data.h5ad")
abundance_df = pd.read_csv(f"{SAVE_PATH}/step1_fdzs_abundance.csv")

adata.obs['sample_id'] = adata.obs['core_name'].astype(str)
abundance_df['sample_id'] = adata.obs['sample_id'].tolist()

# Reconstruct cholangiocyte annotations
cholangiocyte_genes = ['KRT7', 'KRT19', 'EPCAM', 'SOX9', 'AQP1']
available_genes = [gene for gene in cholangiocyte_genes if gene in adata.var_names]

if len(available_genes) > 0:
    gene_subset = adata[:, available_genes].X.toarray()
    cholangiocyte_signature = np.mean(gene_subset, axis=1)
    adata.obs['cholangiocyte_signature'] = cholangiocyte_signature
    
    signature_threshold = np.percentile(cholangiocyte_signature, 95) # 1.0 # Top 99%
    high_signature_spots = cholangiocyte_signature > signature_threshold
    adata.obs['cholangiocyte_enriched'] = high_signature_spots

    print(f"Reconstructed cholangiocyte annotations using {len(available_genes)} genes")
    print(f"Cholangiocyte-enriched spots: {adata.obs['cholangiocyte_enriched'].sum():,}")
else:
    print("❌ Cannot reconstruct cholangiocyte annotations")
    exit()

print(f"Data loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")
print(f"RFS distribution: {adata.obs['RFS_status'].value_counts().to_dict()}")

# Get cell type columns from abundance data
cell_type_cols = [col for col in abundance_df.columns 
                 if col not in ['sample_id', 'RFS_status']]

print(f"Cell types for analysis: {len(cell_type_cols)}")

# Define distance bins for microenvironment analysis (matching IMC approach)
distance_bins = [
    (50, 60),      # Immediate vicinity (0-25μm)
    (60, 80),    # Close proximity (25-50μm)  
    (80, 120),   # Intermediate distance (50-100μm)
    (120,150),
]

print(f"\n📏 Distance bins for microenvironment analysis:")
for i, (start, end) in enumerate(distance_bins):
    print(f"  Bin {i+1}: {start}-{end} μm")

# Analyze all samples using graph-based cholangiocyte zone detection
print(f"\n🔬 Starting graph-based cholangiocyte zone analysis...")
sample_results = {}
sample_ids = adata.obs['sample_id'].unique()

for sample_id in sample_ids:
    result = analyze_sample_cholangiocyte_zones(
        sample_id, adata, abundance_df, 
        max_distance_um=20,  # IMC-aligned threshold
        distance_bins=distance_bins
    )
    if result is not None:
        sample_results[sample_id] = result

print(f"\nCompleted analysis for {len(sample_results)} samples")

# Create zone-level data matrix for PCME classification
print(f"\n📈 Creating zone-level data matrix...")
zone_data = []

for sample_id, sample_result in sample_results.items():
    for zone_id, zone_info in sample_result['zones'].items():
        # Create one row per zone with overall microenvironment signatures
        zone_row = {
            'sample_id': sample_id,
            'zone_id': f"{sample_id}_{zone_id}",
            'rfs_status': zone_info['rfs_status'],
            'n_cholangiocytes': zone_info['n_cholangiocytes'],
            'total_microenv_spots': zone_info['total_microenv_spots'],
            'immune_signature': zone_info['overall_immune_signature'],
            'stromal_signature': zone_info['overall_stromal_signature'],
            'immune_stromal_ratio': zone_info['overall_immune_stromal_ratio'],
            'center_x': zone_info['center'][0],
            'center_y': zone_info['center'][1],
            **zone_info['overall_cell_abundances']
        }
        
        # Add distance-specific data including cholangiocytes
        for bin_name, bin_data in zone_info['distance_analysis'].items():
            zone_row[f'{bin_name}_immune_signature'] = bin_data['immune_signature']
            zone_row[f'{bin_name}_stromal_signature'] = bin_data['stromal_signature']
            zone_row[f'{bin_name}_immune_stromal_ratio'] = bin_data['immune_stromal_ratio']
            zone_row[f'{bin_name}_n_spots'] = bin_data['n_spots']
            
            # ADDED: Include individual cell type abundances for each distance
            for cell_type, abundance in bin_data['cell_abundances'].items():
                zone_row[f'{bin_name}_{cell_type}'] = abundance
        
        zone_data.append(zone_row)

zone_df = pd.DataFrame(zone_data)
print(f"Created zone matrix: {len(zone_df)} cholangiocyte zones across {len(sample_results)} samples")

# Summary statistics
print(f"\n📊 Zone summary:")
print(f"  Total zones: {len(zone_df)}")
print(f"  RFS distribution: {zone_df['rfs_status'].value_counts().to_dict()}")
print(f"  Zones per sample: {zone_df.groupby('sample_id').size().describe()}")
print(f"  Zone size distribution: {zone_df['n_cholangiocytes'].describe()}")

# Compare zones between RFS groups
print(f"\n📈 Comparing zones between RFS groups...")

rfs0_zones = zone_df[zone_df['rfs_status'] == 0]
rfs1_zones = zone_df[zone_df['rfs_status'] == 1]

if len(rfs0_zones) > 0 and len(rfs1_zones) > 0:
    print(f"RFS 0: {len(rfs0_zones)} zones from {rfs0_zones['sample_id'].nunique()} samples")
    print(f"RFS 1: {len(rfs1_zones)} zones from {rfs1_zones['sample_id'].nunique()} samples")
    
    # Test overall microenvironment differences
    stat_immune, p_immune = stats.mannwhitneyu(
        rfs1_zones['immune_signature'], rfs0_zones['immune_signature'], alternative='two-sided'
    )
    stat_stromal, p_stromal = stats.mannwhitneyu(
        rfs1_zones['stromal_signature'], rfs0_zones['stromal_signature'], alternative='two-sided'
    )
    stat_ratio, p_ratio = stats.mannwhitneyu(
        rfs1_zones['immune_stromal_ratio'], rfs0_zones['immune_stromal_ratio'], alternative='two-sided'
    )
    
    print(f"\nZone-level differences (overall microenvironment):")
    print(f"  Immune signature: RFS 0 = {rfs0_zones['immune_signature'].mean():.3f}, RFS 1 = {rfs1_zones['immune_signature'].mean():.3f}, p = {p_immune:.4f}")
    print(f"  Stromal signature: RFS 0 = {rfs0_zones['stromal_signature'].mean():.3f}, RFS 1 = {rfs1_zones['stromal_signature'].mean():.3f}, p = {p_stromal:.4f}")
    print(f"  I/S ratio: RFS 0 = {rfs0_zones['immune_stromal_ratio'].mean():.3f}, RFS 1 = {rfs1_zones['immune_stromal_ratio'].mean():.3f}, p = {p_ratio:.4f}")
    
    # Test distance-specific differences
    print(f"\nDistance-specific differences:")
    for bin_idx, (min_dist, max_dist) in enumerate(distance_bins):
        bin_col = f'bin_{bin_idx}_{min_dist}-{max_dist}um_immune_stromal_ratio'
        if bin_col in zone_df.columns:
            rfs0_bin = rfs0_zones[bin_col].dropna()
            rfs1_bin = rfs1_zones[bin_col].dropna()
            
            if len(rfs0_bin) > 0 and len(rfs1_bin) > 0:
                stat_bin, p_bin = stats.mannwhitneyu(rfs1_bin, rfs0_bin, alternative='two-sided')
                print(f"  {min_dist}-{max_dist}μm I/S ratio: RFS 0 = {rfs0_bin.mean():.3f}, RFS 1 = {rfs1_bin.mean():.3f}, p = {p_bin:.4f}")

# Create enhanced visualization
fig = plt.figure(figsize=(15, 12))

# 1. Zone count distribution per sample
ax1 = plt.subplot(2, 3, 1)
sample_zone_counts = zone_df.groupby(['sample_id', 'rfs_status']).size().reset_index(name='count')
rfs0_counts = sample_zone_counts[sample_zone_counts['rfs_status'] == 0]['count']
rfs1_counts = sample_zone_counts[sample_zone_counts['rfs_status'] == 1]['count']

bp = ax1.boxplot([rfs0_counts, rfs1_counts], labels=['RFS 0', 'RFS 1'], patch_artist=True)
bp['boxes'][0].set_facecolor('lightblue')
bp['boxes'][1].set_facecolor('lightcoral')
ax1.set_ylabel('Zones per Sample')
ax1.set_title('Zone Count Distribution\n(Graph-Based, 30μm threshold)')

for i, data in enumerate([rfs0_counts, rfs1_counts]):
    x = np.random.normal(i+1, 0.04, size=len(data))
    ax1.scatter(x, data, alpha=0.7, s=50)

# 3. Overall immune vs stromal signatures
ax3 = plt.subplot(2, 3, 2)
scatter = ax3.scatter(zone_df['immune_signature'], zone_df['stromal_signature'],
                     c=zone_df['rfs_status'], cmap='RdYlBu_r', alpha=0.7, s=30)
ax3.set_xlabel('Immune Signature')
ax3.set_ylabel('Stromal Signature')
ax3.set_title('Zone Immune vs Stromal')
cbar = plt.colorbar(scatter, ax=ax3)
cbar.set_label('RFS Status')

# 4. Overall I/S ratio distribution
ax4 = plt.subplot(2, 3, 3)
ax4.hist(rfs0_zones['immune_stromal_ratio'], bins=20, alpha=0.6, 
         label='RFS 0', color='lightblue', density=True)
ax4.hist(rfs1_zones['immune_stromal_ratio'], bins=20, alpha=0.6, 
         label='RFS 1', color='lightcoral', density=True)
ax4.set_xlabel('Overall I/S Ratio')
ax4.set_ylabel('Density')
ax4.set_title(f'I/S Ratio Distribution\n(p = {p_ratio:.4f})')
ax4.legend()

# 5-7. Distance-specific I/S ratios
for i, (bin_idx, (min_dist, max_dist)) in enumerate([(0, distance_bins[0]), (1, distance_bins[1]), (2, distance_bins[2])]):
    if i < 3:  # Only plot first 3 bins
        ax = plt.subplot(2, 3, 4 + i)
        bin_col = f'bin_{bin_idx}_{min_dist}-{max_dist}um_immune_stromal_ratio'
        
        if bin_col in zone_df.columns:
            rfs0_bin = rfs0_zones[bin_col].dropna()
            rfs1_bin = rfs1_zones[bin_col].dropna()
            
            if len(rfs0_bin) > 0 and len(rfs1_bin) > 0:
                stat_bin, p_bin = stats.mannwhitneyu(rfs1_bin, rfs0_bin, alternative='two-sided')
                
                bp = ax.boxplot([rfs0_bin, rfs1_bin], labels=['RFS 0', 'RFS 1'], patch_artist=True)
                bp['boxes'][0].set_facecolor('lightblue')
                bp['boxes'][1].set_facecolor('lightcoral')
                
                ax.set_ylabel('I/S Ratio')
                ax.set_title(f'{min_dist}-{max_dist}μm I/S Ratio\n(p = {p_bin:.4f})')

plt.tight_layout()
plt.savefig(f'{SAVE_PATH}/step3_1_graph_based_zone_analysis.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Spatial maps showing zones
fig = plt.figure(figsize=(20, 16))
selected_samples = list(sample_results.keys())
spatial_axes = [plt.subplot(3, 3, 1 + i) for i in range(len(selected_samples))]

# Define colors for each assignment type
color_mapping = {
    'cholangiocyte': '#FF0000',        # Bright red
    'bin_0_50-60um': '#0066CC',        # Deep blue  
    'bin_1_60-80um': '#00AA44',        # Green
    'bin_2_80-120um': '#FF8800',       # Orange (more visible than yellow)
    'bin_3_120-150um': '#8800CC',      # Purple
    'beyond_range': '#E0E0E0',         # Light gray
    'unassigned': '#808080'            # Medium gray
}

# Size mapping - larger for sparse important categories
size_mapping = {
    'cholangiocyte': 20,               # Largest
    'bin_0_50-60um': 12,               # Large for distance bins
    'bin_1_60-80um': 12,
    'bin_2_80-120um': 12,
    'bin_3_120-150um': 12,
    'beyond_range': 1,                 # Tiny for background
    'unassigned': 5                    # Medium
}

# Alpha mapping - less transparent for important categories
alpha_mapping = {
    'cholangiocyte': 1.0,              # Fully opaque
    'bin_0_50-60um': 0.9,              # High opacity
    'bin_1_60-80um': 0.9,
    'bin_2_80-120um': 0.9,
    'bin_3_120-150um': 0.9,
    'beyond_range': 0.6,              # Very transparent
    'unassigned': 0.6                  # Medium transparency
}

for i, (sample_id, ax) in enumerate(zip(selected_samples, spatial_axes)):
    result = sample_results[sample_id]
    all_assignments = pd.concat([result['sample_adata'].obs['distance_assignment'] 
                            for result in sample_results.values()])
    overall_counts = all_assignments.value_counts()

    coords = result['sample_adata'].obs.loc[:,['x','y']]
    zone_assignments = result['sample_adata'].obs['distance_assignment']
    
    # Plot in layers - background first, important categories last
    plot_order = ['beyond_range', 'unassigned', 'bin_3_120-150um', 
                  'bin_2_80-120um', 'bin_1_60-80um', 'bin_0_50-60um', 'cholangiocyte']
    
    for assignment in plot_order:
        if assignment in zone_assignments.values:
            mask = zone_assignments == assignment
            if mask.sum() > 0:
                coords_subset = coords[mask]
                color = color_mapping.get(assignment, 'gray')
                size = size_mapping.get(assignment, 3)
                alpha = alpha_mapping.get(assignment, 0.7)
                
                # Add edge color for cholangiocytes to make them stand out more
                edgecolor = 'black' if assignment == 'cholangiocyte' else None
                linewidth = 1.5 if assignment == 'cholangiocyte' else 0
                
                scatter = ax.scatter(coords_subset.iloc[:, 0], coords_subset.iloc[:, 1], 
                                   c=color, s=size, alpha=alpha,
                                   edgecolors=edgecolor, linewidths=linewidth)
    
    # Enhanced title with more statistics
    rfs_status = result['rfs_status']
    chol_count = result['total_cholangiocytes']
    sample_counts = zone_assignments.value_counts()
    bin_spots = sum([sample_counts.get(f'bin_{j}_{ranges}', 0) 
                    for j, ranges in enumerate(['50-60um', '60-80um', '80-120um', '120-150um'])])
    
    ax.set_title(f'{sample_id}\nRFS: {rfs_status} | Chol: {chol_count} | Bins: {bin_spots}',
                fontweight='bold', fontsize=11)
    ax.set_xlabel('X coordinate', fontweight='bold')
    ax.set_ylabel('Y coordinate', fontweight='bold')
    
    # Add subtle grid for better spatial reference
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_aspect('equal', adjustable='box')

# Enhanced legend with counts and better formatting
fig = plt.gcf()
legend_handles = []

# Create legend in logical order with counts
legend_order = ['cholangiocyte', 'bin_0_50-60um', 'bin_1_60-80um', 
               'bin_2_80-120um', 'bin_3_120-150um', 'unassigned', 'beyond_range']

for assignment in legend_order:
    if assignment in overall_counts.index:
        color = color_mapping.get(assignment, 'gray')
        count = overall_counts[assignment]
        percentage = (count / len(all_assignments)) * 100
        
        # Create clean display names with counts
        if assignment == 'cholangiocyte':
            display_name = f'Cholangiocyte (n={count})'
            # Special marker for cholangiocytes with border
            patch = plt.Line2D([0], [0], marker='o', color='w', 
                             markerfacecolor=color, markersize=10,
                             markeredgecolor='black', markeredgewidth=1.5,
                             label=display_name, linestyle='None')
        elif assignment.startswith('bin_'):
            distance_range = assignment.split('_')[1] + assignment.split('_')[2]
            distance_range = distance_range.replace('um', 'μm')
            display_name = f'Zone {distance_range} (n={count})'
            patch = mpatches.Patch(color=color, label=display_name)
        elif assignment == 'beyond_range':
            display_name = f'Beyond Range (n={count}, {percentage:.1f}%)'
            patch = mpatches.Patch(color=color, label=display_name)
        elif assignment == 'unassigned':
            display_name = f'Unassigned (n={count})'
            patch = mpatches.Patch(color=color, label=display_name)
        else:
            display_name = f'{assignment} (n={count})'
            patch = mpatches.Patch(color=color, label=display_name)
        
        legend_handles.append(patch)

# Enhanced figure-level legend
legend = fig.legend(handles=legend_handles,
                   bbox_to_anchor=(0.98, 0.85),
                   loc='upper left',
                   frameon=True,
                   framealpha=0.95,
                   edgecolor='black',
                   title='Distance Assignment',
                   title_fontsize=12,
                   fontsize=10,
                   handlelength=1.5,
                   handletextpad=0.5,
                   borderpad=0.8)

# Enhanced figure title with overall statistics
total_samples = len(selected_samples)
total_spots = len(all_assignments)
total_chol = overall_counts.get('cholangiocyte', 0)
total_bins = sum([overall_counts.get(f'bin_{j}_{ranges}', 0) 
                 for j, ranges in enumerate(['50-60um', '60-80um', '80-120um', '120-150um'])])

fig.suptitle(f'Spatial Distance Analysis - {total_samples} Samples\n'
            f'Total: {total_spots:,} spots | Cholangiocytes: {total_chol} | Zone spots: {total_bins}',
            fontsize=16, fontweight='bold', y=0.96)

plt.tight_layout()
plt.subplots_adjust(right=0.85, top=0.90)  # Make room for legend and title
plt.savefig(f'{SAVE_PATH}/step3_1_spatial_distance_analysis.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Save zone-level results
print(f"\n💾 Saving graph-based zone analysis results...")
zone_df.to_csv(f"{SAVE_PATH}/step3_1_zone_data_graph_based.csv", index=False)

# Save sample results for next step
with open(f"{SAVE_PATH}/step3_1_sample_results_graph_based.pkl", 'wb') as f:
    pickle.dump(sample_results, f)

print(f"✅ Step 3.1 complete: Graph-based zone analysis saved")

#%% Step 3.2: PCME-I vs PCME-S Classification

# Load data from Step 3.1
SAVE_PATH = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/PCME"
zone_df = pd.read_csv(f"{SAVE_PATH}/step3_1_zone_data_graph_based.csv")
with open(f"{SAVE_PATH}/step3_1_sample_results_graph_based.pkl", 'rb') as f:
    sample_results = pickle.load(f)

print(f"Loaded zone data: {len(zone_df)} cholangiocyte zones")
print(f"RFS distribution: {zone_df['rfs_status'].value_counts().to_dict()}")

# Try different methods
classification_results = {}
for method in ['percentile', 'kmeans', 'adaptive']:
    classifications, threshold = classify_zones_pcme(zone_df, method=method, 
                                                     distence_prefix='cholangiocyte_distance_0', use_overall=False)
    classification_results[method] = {
        'classifications': classifications,
        'threshold': threshold
    }

# Analyze relationship with RFS status
print(f"\n📈 Analyzing PCME classification vs RFS status...")
zone_df['pcme_classification'] = classification_results['kmeans']['classifications']

# Zone-level analysis
pcme_rfs_crosstab = pd.crosstab(
    zone_df['pcme_classification'], 
    zone_df['rfs_status'], 
    normalize='columns'
) * 100

print(f"PCME classification by RFS status (column percentages):")
print(pcme_rfs_crosstab.round(1))

# Statistical test
from scipy.stats import chi2_contingency
pcme_rfs_contingency = pd.crosstab(zone_df['pcme_classification'], zone_df['rfs_status'])
chi2, p_value, dof, expected = chi2_contingency(pcme_rfs_contingency)
print(f"\nChi-square test: χ²={chi2:.3f}, p={p_value:.4f}")

# Sample-level aggregation
print(f"\n📊 Sample-level PCME aggregation...")
sample_pcme_df = calculate_sample_pcme_metrics_zones(zone_df)
print(f"Sample-level PCME summary:")
print(sample_pcme_df.round(2))
# Test sample-level differences
if len(sample_pcme_df) > 0:
    rfs0_samples = sample_pcme_df[sample_pcme_df['rfs_status'] == 0]
    rfs1_samples = sample_pcme_df[sample_pcme_df['rfs_status'] == 1]
    
    if len(rfs0_samples) > 0 and len(rfs1_samples) > 0:
        # Test PCME-I percentage
        stat_i, p_val_i = mannwhitneyu(
            rfs1_samples['pcme_i_percentage'], rfs0_samples['pcme_i_percentage'], alternative='two-sided'
        )
        
        # Test PCME-S percentage  
        stat_s, p_val_s = mannwhitneyu(
            rfs1_samples['pcme_s_percentage'], rfs0_samples['pcme_s_percentage'], alternative='two-sided'
        )
        
        # Test weighted I/S ratio
        stat_ratio, p_val_ratio = mannwhitneyu(
            rfs1_samples['weighted_is_ratio'], rfs0_samples['weighted_is_ratio'], alternative='two-sided'
        )
        
        print(f"\nSample-level statistical tests:")
        print(f"PCME-I percentage: RFS 0 = {rfs0_samples['pcme_i_percentage'].mean():.1f}%, RFS 1 = {rfs1_samples['pcme_i_percentage'].mean():.1f}%, p = {p_val_i:.4f}")
        print(f"PCME-S percentage: RFS 0 = {rfs0_samples['pcme_s_percentage'].mean():.1f}%, RFS 1 = {rfs1_samples['pcme_s_percentage'].mean():.1f}%, p = {p_val_s:.4f}")
        print(f"Weighted I/S ratio: RFS 0 = {rfs0_samples['weighted_is_ratio'].mean():.3f}, RFS 1 = {rfs1_samples['weighted_is_ratio'].mean():.3f}, p = {p_val_ratio:.4f}")

# Create comprehensive visualization
fig = plt.figure(figsize=(20, 12))

# Define celltypes
cell_types = define_celltypes()
immune_cell_types = cell_types["immune_cell_types"]
stromal_cell_types = cell_types["stromal_cell_types"]

# Define consistent colors
pcme_colors = {'PCME-I': '#E31A1C', 'PCME-S': '#1F78B4', 'Intermediate': '#888888'}

# ================== TOP ROW: OVERVIEW & CLASSIFICATION ==================

# 1. PCME classification scatter plot (Main result)
ax1 = plt.subplot(2, 4, 1)
for pcme_type, color in pcme_colors.items():
    mask = zone_df['pcme_classification'] == pcme_type
    if np.any(mask):
        ax1.scatter(zone_df.loc[mask, 'immune_signature'], 
                   zone_df.loc[mask, 'stromal_signature'],
                   c=color, label=pcme_type, alpha=0.7, s=40, edgecolors='black', linewidth=0.5)

ax1.set_xlabel('Immune Signature Score', fontweight='bold')
ax1.set_ylabel('Stromal Signature Score', fontweight='bold')
ax1.set_title('PCME Classification Overview', fontweight='bold', fontsize=12)
ax1.legend(framealpha=0.9, edgecolor='black')
ax1.grid(True, alpha=0.3)

# 2. PCME distribution by RFS status (Clinical relevance)
ax2 = plt.subplot(2, 4, 2)
pcme_rfs_counts = pd.crosstab(zone_df['pcme_classification'], zone_df['rfs_status'])

# Create stacked bar plot with better colors
pcme_rfs_counts.plot(kind='bar', ax=ax2, 
                    color=['#87CEEB', '#FF6B6B'], alpha=0.8, 
                    edgecolor='black', linewidth=0.5)
ax2.set_xlabel('PCME Classification', fontweight='bold')
ax2.set_ylabel('Number of Zones', fontweight='bold')
ax2.set_title('PCME Distribution by RFS Status', fontweight='bold', fontsize=12)
ax2.legend(['RFS 0 (Better)', 'RFS 1 (Worse)'], framealpha=0.9)
ax2.tick_params(axis='x', rotation=45)

# Add count labels on bars
for container in ax2.containers:
    ax2.bar_label(container, label_type='center', fontweight='bold')

# 3. Immune/Stromal ratio distribution (Mechanism)
ax3 = plt.subplot(2, 4, 3)
for pcme_type, color in pcme_colors.items():
    pcme_data = zone_df[zone_df['pcme_classification'] == pcme_type]
    if len(pcme_data) > 0:
        ax3.hist(pcme_data['immune_stromal_ratio'], bins=15, alpha=0.6, 
                label=pcme_type, color=color, density=True, edgecolor='black', linewidth=0.5)

ax3.set_xlabel('Immune/Stromal Ratio', fontweight='bold')
ax3.set_ylabel('Density', fontweight='bold')
ax3.set_title('I/S Ratio Distribution by PCME Type', fontweight='bold', fontsize=12)
ax3.legend(framealpha=0.9)
ax3.grid(True, alpha=0.3)

# 4. Summary statistics table
ax4 = plt.subplot(2, 4, 4)
ax4.axis('off')

# Create summary statistics
summary_stats = []
for pcme_type in ['PCME-I', 'PCME-S', 'Intermediate']:
    pcme_data = zone_df[zone_df['pcme_classification'] == pcme_type]
    if len(pcme_data) > 0:
        stats = {
            'PCME Type': pcme_type,
            'Count': len(pcme_data),
            'Percentage': f"{len(pcme_data)/len(zone_df)*100:.1f}%",
            'Mean I/S Ratio': f"{pcme_data['immune_stromal_ratio'].mean():.2f}",
            'RFS 1 Rate': f"{(pcme_data['rfs_status'] == 1).mean()*100:.1f}%"
        }
        summary_stats.append(stats)

summary_df = pd.DataFrame(summary_stats)

# Create table
table = ax4.table(cellText=summary_df.values,
                 colLabels=summary_df.columns,
                 cellLoc='center',
                 loc='center',
                 bbox=[0, 0, 1, 1])

table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 2)

# Style the table
for i in range(len(summary_df.columns)):
    table[(0, i)].set_facecolor('#D3D3D3')
    table[(0, i)].set_text_props(weight='bold')

for i in range(1, len(summary_df) + 1):
    pcme_type = summary_df.iloc[i-1]['PCME Type']
    color = pcme_colors.get(pcme_type, '#FFFFFF')
    for j in range(len(summary_df.columns)):
        table[(i, j)].set_facecolor(color)
        table[(i, j)].set_alpha(0.3)

ax4.set_title('PCME Classification Summary', fontweight='bold', fontsize=12, pad=20)

# ================== BOTTOM ROW: DETAILED CELL TYPE ANALYSIS ==================

# 5. Key immune cell types in PCME-I vs PCME-S
ax5 = plt.subplot(2, 4, 5)
pcme_i_data = zone_df[zone_df['pcme_classification'] == 'PCME-I']
pcme_s_data = zone_df[zone_df['pcme_classification'] == 'PCME-S']

if len(pcme_i_data) > 0 and len(pcme_s_data) > 0 and len(immune_cell_types) > 0:
    x = range(len(immune_cell_types))
    width = 0.35
    
    pcme_i_means = [pcme_i_data[ct].mean() for ct in immune_cell_types]
    pcme_s_means = [pcme_s_data[ct].mean() for ct in immune_cell_types]
    
    bars1 = ax5.bar([i - width/2 for i in x], pcme_i_means, width, 
                   label='PCME-I', color=pcme_colors['PCME-I'], alpha=0.8,
                   edgecolor='black', linewidth=0.5)
    bars2 = ax5.bar([i + width/2 for i in x], pcme_s_means, width, 
                   label='PCME-S', color=pcme_colors['PCME-S'], alpha=0.8,
                   edgecolor='black', linewidth=0.5)
    
    ax5.set_xlabel('Immune Cell Type', fontweight='bold')
    ax5.set_ylabel('Mean Abundance', fontweight='bold')
    ax5.set_title('Immune Cell Types: PCME-I vs PCME-S', fontweight='bold', fontsize=12)
    ax5.set_xticks(x)
    
    # Clean cell type names
    clean_labels = []
    for ct in immune_cell_types:
        clean_name = ct.replace('q05cell_abundance_w_sf_', '').replace('_', ' ')
        # Wrap long names
        if len(clean_name) > 10:
            words = clean_name.split()
            if len(words) >= 2:
                clean_name = words[0] + '\n' + ' '.join(words[1:])
        clean_labels.append(clean_name)
    
    ax5.set_xticklabels(clean_labels, fontsize=8, rotation=0)
    ax5.legend(framealpha=0.9)
    ax5.grid(axis='y', alpha=0.3)

# 6. Key stromal cell types in PCME-I vs PCME-S
ax6 = plt.subplot(2, 4, 6)
if len(pcme_i_data) > 0 and len(pcme_s_data) > 0 and len(stromal_cell_types) > 0:
    x = range(len(stromal_cell_types))
    width = 0.35
    
    pcme_i_means = [pcme_i_data[ct].mean() for ct in stromal_cell_types]
    pcme_s_means = [pcme_s_data[ct].mean() for ct in stromal_cell_types]
    
    bars1 = ax6.bar([i - width/2 for i in x], pcme_i_means, width, 
                   label='PCME-I', color=pcme_colors['PCME-I'], alpha=0.8,
                   edgecolor='black', linewidth=0.5)
    bars2 = ax6.bar([i + width/2 for i in x], pcme_s_means, width, 
                   label='PCME-S', color=pcme_colors['PCME-S'], alpha=0.8,
                   edgecolor='black', linewidth=0.5)
    
    ax6.set_xlabel('Stromal Cell Type', fontweight='bold')
    ax6.set_ylabel('Mean Abundance', fontweight='bold')
    ax6.set_title('Stromal Cell Types: PCME-I vs PCME-S', fontweight='bold', fontsize=12)
    ax6.set_xticks(x)
    
    # Clean stromal cell type names
    clean_labels = []
    for ct in stromal_cell_types:
        clean_name = ct.replace('q05cell_abundance_w_sf_', '').replace('_', ' ')
        if len(clean_name) > 10:
            words = clean_name.split()
            if len(words) >= 2:
                clean_name = words[0] + '\n' + ' '.join(words[1:])
        clean_labels.append(clean_name)
    
    ax6.set_xticklabels(clean_labels, fontsize=8, rotation=0)
    ax6.legend(framealpha=0.9)
    ax6.grid(axis='y', alpha=0.3)

# 7. Statistical significance heatmap
ax7 = plt.subplot(2, 4, 7)

# Calculate p-values for all cell types between PCME-I and PCME-S
all_cell_types = immune_cell_types + stromal_cell_types
distence_prefix='cholangiocyte_distance_0'

if len(pcme_i_data) > 0 and len(pcme_s_data) > 0:
    pval_results = []
    fold_changes = []
    
    for ct in all_cell_types:
        if ct in pcme_i_data.columns and ct in pcme_s_data.columns:
            i_values = pcme_i_data[f"{distence_prefix}_{ct}"].values
            s_values = pcme_s_data[f"{distence_prefix}_{ct}"].values
            
            # Perform t-test
            try:
                _, pval = ttest_ind(i_values, s_values)
                fold_change = np.log2((np.mean(i_values) + 1e-8) / (np.mean(s_values) + 1e-8))
            except:
                pval = 1.0
                fold_change = 0.0
            
            pval_results.append(pval)
            fold_changes.append(fold_change)
    
    # Create significance matrix
    significance_data = pd.DataFrame({
        'Cell_Type': [ct.replace('q05cell_abundance_w_sf_', '') for ct in all_cell_types],
        'P_Value': pval_results,
        'Log2_FC': fold_changes,
        'Significant': [p < 0.05 for p in pval_results]
    })
    
    # Plot as horizontal bar chart of fold changes
    y_pos = range(len(significance_data))
    colors = ['red' if fc > 0 else 'blue' for fc in significance_data['Log2_FC']]
    bars = ax7.barh(y_pos, significance_data['Log2_FC'], color=colors, alpha=0.7,
                   edgecolor='black', linewidth=0.5)
    
    # Mark significant ones
    for i, (_, row) in enumerate(significance_data.iterrows()):
        if row['Significant']:
            ax7.text(row['Log2_FC'] + 0.01 if row['Log2_FC'] > 0 else row['Log2_FC'] - 0.01, 
                    i, '*', va='center', ha='left' if row['Log2_FC'] > 0 else 'right',
                    fontsize=12, fontweight='bold')
    
    ax7.set_yticks(y_pos)
    ax7.set_yticklabels([name.replace('_', ' ') for name in significance_data['Cell_Type']], 
                       fontsize=8)
    ax7.set_xlabel('Log₂ Fold Change (PCME-I vs PCME-S)', fontweight='bold')
    ax7.set_title('Cell Type Enrichment\n(* p < 0.05)', fontweight='bold', fontsize=12)
    ax7.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    ax7.grid(axis='x', alpha=0.3)

# 8. PCME composition pie chart
ax8 = plt.subplot(2, 4, 8)
pcme_counts = zone_df['pcme_classification'].value_counts()
colors_list = [pcme_colors[pcme_type] for pcme_type in pcme_counts.index]

wedges, texts, autotexts = ax8.pie(pcme_counts.values, labels=pcme_counts.index, 
                                  colors=colors_list, autopct='%1.1f%%',
                                  startangle=90, explode=[0.05]*len(pcme_counts),
                                  textprops={'fontweight': 'bold'})

ax8.set_title('PCME Type Distribution', fontweight='bold', fontsize=12)

# Enhance pie chart text
for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontsize(10)
    autotext.set_fontweight('bold')

# Overall figure adjustments
plt.suptitle('PCME Classification Analysis', fontsize=16, fontweight='bold', y=0.95)
plt.tight_layout()
plt.subplots_adjust(top=0.92, hspace=0.3, wspace=0.3)

# Save the figure
plt.savefig(f'{SAVE_PATH}/step3_2_PCME_classification.pdf', dpi=300, bbox_inches='tight')
plt.show()

zone_df.to_csv(f"{SAVE_PATH}/step3_2_zone_pcme_data_graph_based.csv", index=False)
sample_pcme_df.to_csv(f"{SAVE_PATH}/step3_2_sample_pcme_summary_graph_based.csv", index=False)


#%% Step 4.1: PCME-Specific Pathway Analysis

# Load PCME classification results from Step 3.2
SAVE_PATH = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/PCME"
zone_df = pd.read_csv(f"{SAVE_PATH}/step3_2_zone_pcme_data_graph_based.csv")
sample_pcme_df = pd.read_csv(f"{SAVE_PATH}/step3_2_sample_pcme_summary_graph_based.csv")

with open(f"{SAVE_PATH}/step3_1_sample_results_graph_based.pkl", 'rb') as f:
    sample_results = pickle.load(f)

# Load original spatial data
# adata = sc.read_h5ad(f"{SAVE_PATH}/step1_fdzs_combined_data.h5ad")
# adata.obs['sample_id'] = adata.obs['core_name'].astype(str)

print(f"Loaded data: {len(zone_df)} zones, {adata.n_obs:,} spots")
print(f"PCME distribution: {zone_df['pcme_classification'].value_counts().to_dict()}")

# 1. Extract PCME region expression data
pcme_expression_data = extract_pcme_region_expression(adata, zone_df, sample_results, max_radius = 50)

# 2. Perform differential expression analysis
deg_results = perform_pcme_differential_analysis(adata, pcme_expression_data)

# 3. Pathway enrichment analysis
enrichment_results = perform_pathway_enrichment_analysis(deg_results)
    
# 4. Save results
print(f"\n💾 Saving Step 4.1 results...")
deg_results.to_csv(f"{SAVE_PATH}/step4_1_pcme_deg_results.csv", index=False)

# Save enrichment results
for direction, enrichment_df in enrichment_results.items():
    enrichment_df.to_csv(f"{SAVE_PATH}/step4_1_pathway_enrichment_{direction}.csv", index=False)

# Set style for publication-quality figures
plt.style.use('default')
sns.set_palette("husl")

# 1. Create individual volcano plot
print("Creating volcano plot...")
fig_volcano, ax_volcano = create_volcano_plot(deg_results)
fig_volcano.savefig(f'{SAVE_PATH}/step4_1_volcano_plot.pdf', dpi=300, bbox_inches='tight')
plt.show()

# 2. Create individual pathway enrichment plots
print("Creating pathway enrichment plots...")
fig_pathways = create_pathway_enrichment_plots(enrichment_results)
fig_pathways.savefig(f'{SAVE_PATH}/step4_1_pathway_enrichment.pdf', dpi=300, bbox_inches='tight')
plt.show()

print("✅ All figures created and saved successfully!")
print(f"📁 Figures saved to: {SAVE_PATH}")
print("   - step4_1_volcano_plot.pdf")
print("   - step4_1_pathway_enrichment.pdf") 
print("   - step4_1_comprehensive_analysis.pdf")