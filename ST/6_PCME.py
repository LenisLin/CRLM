"""
PCME analysis on stereo-seq data
"""
import os
import pandas as pd
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

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

# Define FDZS A/C samples only
fdzs_samples = {
    # A group (Early Recurrence)
    'FDZS_A04932G3_bin50_A1': {'group': 'A', 'recurrence': 'Early_Recurrence'},
    'FDZS_A04932G3_bin50_A2': {'group': 'A', 'recurrence': 'Early_Recurrence'},
    'FDZS_A04932G3_bin50_A3': {'group': 'A', 'recurrence': 'Early_Recurrence'},
    'FDZS_A04932G3_bin50_A5': {'group': 'A', 'recurrence': 'Early_Recurrence'},
    'FDZS_A04932G3_bin50_A6': {'group': 'A', 'recurrence': 'Early_Recurrence'},
    'FDZS_A04932G3_bin50_A7': {'group': 'A', 'recurrence': 'Early_Recurrence'},
    
    # C group (Non-Early Recurrence)
    'FDZS_A04932G3_bin50_C5': {'group': 'C', 'recurrence': 'Non_Early_Recurrence'},
    'FDZS_A04932G3_bin50_C6': {'group': 'C', 'recurrence': 'Non_Early_Recurrence'},
    'FDZS_A04932G3_bin50_C7': {'group': 'C', 'recurrence': 'Non_Early_Recurrence'},
}

print(f"\nTarget samples: {len(fdzs_samples)} FDZS slides")
print(f"A group (Early Recurrence): {len([s for s in fdzs_samples.values() if s['group'] == 'A'])} samples")
print(f"C group (Non-Early Recurrence): {len([s for s in fdzs_samples.values() if s['group'] == 'C'])} samples")

# Load each sample
data_path = Path(DATA_PATH)
adata_list = []

print(f"\nLoading samples from: {data_path}")

for sample_name, sample_info in fdzs_samples.items():
    print(f"\nProcessing {sample_name}...")
    
    # Find h5ad file in sample folder
    sample_path = data_path / sample_name
    h5ad_files = list(sample_path.glob("*.h5ad"))
    
    if len(h5ad_files) == 0:
        print(f"  ❌ No h5ad file found in {sample_path}")
        continue
    
    # Load the data
    adata = sc.read_h5ad(h5ad_files[0])

    # Add meta information
    adata.obs["n_spots"] = adata.n_obs
    adata.obs["n_genes"] = adata.n_vars
    adata.obs["sample_id"] = adata.obs['sample'].astype(str).values + "_" + adata.obs['core_name'].astype(str).values

    # Add to list
    adata_list.append(adata)
    print(f"  ✅ Loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")

# Check loading results
print(f"Successfully loaded: {len(adata_list)} / {len(fdzs_samples)} samples")

# Combine all samples
print(f"\nCombining {len(adata_list)} samples...")
adata_combined = sc.concat(
    adata_list, 
    join='outer', 
    label='batch_sample',
    keys=[adata.obs['core_name'].iloc[0] for adata in adata_list]
)

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
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.patches import Patch

from PCME_functions import plot_sample_overview, plot_cell_type_overview

adata = sc.read_h5ad("step1_combined_data.h5ad")
metadata_df = pd.read_csv("step1_sample_metadata.csv")
abundance_df = pd.read_csv("step1_cell_abundance.csv")

# Create visualizations
print("Creating sample overview plots...")
fig1 = plot_sample_overview(metadata_df)
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
signature_threshold = 1.0 # np.percentile(cholangiocyte_signature, 99) # Top 99%
high_signature_spots = cholangiocyte_signature > signature_threshold

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
from scipy import stats

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
    for gene in available_genes[:4]:  # Show top 4 genes
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
hypoxia_genes = ['CA9', 'HIF1A', 'VEGFA',] #  'LDHA', 'SLC2A1', 'PKM', 'BNIP3'
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
caix_threshold = 0.5  # Top 10%

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

print(f"Data loaded: {adata.n_obs:,} spots, {adata.n_vars:,} genes")
print(f"RFS distribution: {adata.obs['RFS_status'].value_counts().to_dict()}")

# Get cell type columns from abundance data
cell_type_cols = [col for col in abundance_df.columns 
                 if col not in ['sample_id', 'RFS_status']]

print(f"Cell types for analysis: {len(cell_type_cols)}")

# Define distance bins for analysis (in micrometers)
distance_bins = [
    (0, 25),      # Immediate vicinity
    (25, 50),    # Close proximity  
    (50, 100),   # Intermediate distance
    (100, 150),   # Distant proximity
]

print(f"\n📏 Distance bins defined:")
for i, (start, end) in enumerate(distance_bins):
    print(f"  Bin {i+1}: {start}-{end} μm")

# Analyze all samples
print(f"\n🔬 Starting spatial analysis for all samples...")
sample_results = {}
sample_ids = adata.obs['sample_id'].unique()

for sample_id in sample_ids:
    result = analyze_sample_spatial_zones(sample_id, adata, abundance_df, distance_bins)
    if result is not None:
        sample_results[sample_id] = result

print(f"\nCompleted analysis for {len(sample_results)} samples")

# Aggregate results across samples
print(f"\n📈 Aggregating results across samples...")

# Create zone-level cell abundance matrix
zone_abundance_data = []
zone_info_data = []

for sample_id, result in sample_results.items():
    rfs_status = result['rfs_status']
    
    for zone_name, zone_data in result['zone_composition'].items():
        # Zone information
        zone_info = {
            'sample_id': sample_id,
            'rfs_status': rfs_status,
            'zone_name': zone_name,
            'spot_count': zone_data['spot_count'],
            'percentage': zone_data['percentage'],
            'mean_distance': zone_data['mean_distance'],
            'std_distance': zone_data['std_distance']
        }
        zone_info_data.append(zone_info)
        
        # Cell abundance data
        abundance_row = {
            'sample_id': sample_id,
            'rfs_status': rfs_status,
            'zone_name': zone_name,
            **zone_data['cell_abundances']
        }
        zone_abundance_data.append(abundance_row)
