import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

import scipy
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.spatial.distance import cdist
from scipy.stats import gaussian_kde, mannwhitneyu
from sklearn.preprocessing import MinMaxScaler

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%% COLOR DEFINITIONS - Consistent color scheme for all plots
def create_color_mapping(cell_types):
    """Create consistent color mapping for cell types"""
    # Define color categories
    color_categories = {
        'Tumor': '#E74C3C',      # Red
        'Endothelial': '#3498DB', # Blue  
        'Immune': '#2ECC71',      # Green
        'Stromal': '#F39C12',     # Orange
        'Epithelial': '#9B59B6',  # Purple
        'Hepatocyte': '#1ABC9C',  # Teal
        'Others': '#95A5A6'       # Gray
    }
    
    # Specific cell type colors
    cell_type_colors = {}
    
    for cell_type in cell_types:
        if cell_type.startswith('TC_'):
            if 'Glycolysis' in cell_type:
                cell_type_colors[cell_type] = '#C0392B'  # Dark red
            elif 'Quiescent' in cell_type:
                cell_type_colors[cell_type] = '#E74C3C'  # Red
            else:
                cell_type_colors[cell_type] = '#E74C3C'  # Red (tumor)
        elif cell_type.startswith('EC_'):
            cell_type_colors[cell_type] = '#3498DB'  # Blue (endothelial)
        elif any(x in cell_type for x in ['CD4T_', 'CD8T_', 'NK_', 'Macro_', 'B_']):
            cell_type_colors[cell_type] = '#2ECC71'  # Green (immune)
        elif any(x in cell_type for x in ['CAF_', 'Fibro_']):
            cell_type_colors[cell_type] = '#F39C12'  # Orange (stromal)
        elif 'Epithelial' in cell_type:
            cell_type_colors[cell_type] = '#9B59B6'  # Purple
        elif 'Hepatocyte' in cell_type:
            cell_type_colors[cell_type] = '#1ABC9C'  # Teal
        else:
            cell_type_colors[cell_type] = '#95A5A6'  # Gray (others)
    
    return cell_type_colors, color_categories

def perform_wilcox_test(data, group_col, value_col):
    """Perform Wilcoxon rank-sum test between groups"""
    groups = data[group_col].unique()
    if len(groups) != 2:
        return None, None
    
    group1_data = data[data[group_col] == groups[0]][value_col].dropna()
    group2_data = data[data[group_col] == groups[1]][value_col].dropna()
    
    if len(group1_data) == 0 or len(group2_data) == 0:
        return None, None
    
    statistic, p_value = mannwhitneyu(group1_data, group2_data, alternative='two-sided')
    return statistic, p_value

def add_significance_bar(ax, x1, x2, y, p_value, height_offset=0.05):
    """Add significance bar to plot"""
    if p_value < 0.001:
        sig_text = '***'
    elif p_value < 0.01:
        sig_text = '**'
    elif p_value < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'
    
    y_max = ax.get_ylim()[1]
    bar_height = y_max * (1 + height_offset)
    
    ax.plot([x1, x1, x2, x2], [y, bar_height, bar_height, y], 'k-', linewidth=1)
    ax.text((x1 + x2) * 0.5, bar_height, f'{sig_text}\np={p_value:.3f}', 
            ha='center', va='bottom', fontsize=10)

#%% PHASE I: Data Integration and Visualization

# Define paths
run_name = '/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/analysis/cell2location_map'
figure_path = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/multi_sample_analysis"

if not os.path.exists(figure_path):
    os.makedirs(figure_path, exist_ok=True)

# Get all slides and filter relevant ones
all_slides = os.listdir(run_name)
print("All available slides:")
for slide in all_slides:
    print(f"  {slide}")

# Step 1: Multi-Sample Data Integration with Treatment Annotation
print("\n" + "="*60)
print("STEP 1: MULTI-SAMPLE DATA INTEGRATION")
print("="*60)

# Define relevant slides and treatment status
relevant_slides = []
treatment_mapping = {}

for slide in all_slides:
    if slide.startswith('FDZS'):
        print(f"Excluding peritumor-focused slide: {slide}")
        continue
    elif 'Untreated' in slide:
        relevant_slides.append(slide)
        treatment_mapping[slide] = 'untreated'
    elif 'NACPR' in slide:
        relevant_slides.append(slide) 
        treatment_mapping[slide] = 'pre_chemotherapy'
    elif slide.startswith('GSE225857'):
        relevant_slides.append(slide)
        # Assume these are untreated based on context
        treatment_mapping[slide] = 'untreated'

print(f"\nSelected {len(relevant_slides)} slides for analysis:")
for slide in relevant_slides:
    print(f"  {slide} -> {treatment_mapping[slide]}")

# Load and integrate all relevant slides
adata_list = []
sample_info = []

for slide_name in relevant_slides:
    print(f"\nLoading slide: {slide_name}")
    slide_path = os.path.join(run_name, slide_name, "sp.h5ad")
    
    if os.path.exists(slide_path):
        adata_slide = sc.read_h5ad(slide_path)
        # Add cell type abundances from cell2location
        adata_slide.obs[adata_slide.uns['mod']['factor_names']] = adata_slide.obsm['q05_cell_abundance_w_sf']
        
        # Add metadata
        adata_slide.obs['sample_id'] = slide_name
        adata_slide.obs['treatment_status'] = treatment_mapping[slide_name]
        
        adata_list.append(adata_slide)
        sample_info.append({
            'sample_id': slide_name,
            'treatment': treatment_mapping[slide_name],
            'n_spots': adata_slide.n_obs
        })
        
        print(f"  Loaded {adata_slide.n_obs} spots")
    else:
        print(f"  Warning: File not found - {slide_path}")

# Concatenate all slides
if len(adata_list) > 0:
    adata_combined = sc.concat(adata_list, axis=0, join='outer', index_unique='-')
    print(f"\nCombined dataset: {adata_combined.n_obs} spots, {adata_combined.n_vars} genes")
    print(f"Treatment groups: {adata_combined.obs['treatment_status'].value_counts().to_dict()}")
else:
    raise ValueError("No valid slides found!")

# Create sample info DataFrame
sample_df = pd.DataFrame(sample_info)
print("\nSample Summary:")
print(sample_df)

#%% Step 2: Comprehensive Data Visualization
print("\n" + "="*60)
print("STEP 2: COMPREHENSIVE DATA VISUALIZATION")
print("="*60)

# Get all cell types for visualization
all_cell_types = adata_slide.obs[adata_slide.uns['mod']['factor_names']].columns

# Create consistent color mapping
cell_type_colors, color_categories = create_color_mapping(all_cell_types)

# Focus on key cell types for main analysis
tumor_types = [col for col in all_cell_types if col.startswith('TC_')]
ec_types = [col for col in all_cell_types if col.startswith('EC_')]
immune_types = [col for col in all_cell_types if any(x in col for x in ['CD4T_', 'CD8T_', 'NK_', 'Macro_', 'B_'])]
stromal_types = [col for col in all_cell_types if any(x in col for x in ['CAF_', 'Fibro_'])]

print(f"Tumor cell types ({len(tumor_types)}): {tumor_types}")
print(f"Endothelial cell types ({len(ec_types)}): {ec_types}")
print(f"Immune cell types ({len(immune_types)}): {immune_types[:5]}...")  # Show first 5
print(f"Stromal cell types ({len(stromal_types)}): {stromal_types}")

# A) Sample Overview Visualizations

# Calculate cell type composition per sample
composition_data = []
for sample in adata_combined.obs['sample_id'].unique():
    sample_mask = adata_combined.obs['sample_id'] == sample
    sample_data = adata_combined.obs[sample_mask]
    treatment = sample_data['treatment_status'].iloc[0]
    
    # Calculate mean abundance for each cell type
    for cell_type in all_cell_types:
        if cell_type in sample_data.columns:
            mean_abundance = sample_data[cell_type].mean()
            composition_data.append({
                'sample_id': sample,
                'treatment': treatment,
                'cell_type': cell_type,
                'mean_abundance': mean_abundance
            })

composition_df = pd.DataFrame(composition_data)

# Create pie charts for each sample with consistent colors
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for i, sample in enumerate(adata_combined.obs['sample_id'].unique()):
    if i < len(axes):
        sample_data = composition_df[composition_df['sample_id'] == sample]
        treatment = sample_data['treatment'].iloc[0]
        
        # Get top cell types (>1% abundance)
        major_types = sample_data[sample_data['mean_abundance'] > 0.01].copy()
        major_types = major_types.sort_values('mean_abundance', ascending=False)
        
        # Group small cell types as "Others"
        if len(major_types) > 8:  # Show top 8 + others
            top_types = major_types.head(8)
            others_abundance = major_types.tail(-8)['mean_abundance'].sum()
            if others_abundance > 0:
                others_row = pd.DataFrame({
                    'cell_type': ['Others'],
                    'mean_abundance': [others_abundance]
                })
                plot_data = pd.concat([top_types[['cell_type', 'mean_abundance']], others_row])
            else:
                plot_data = top_types[['cell_type', 'mean_abundance']]
        else:
            plot_data = major_types[['cell_type', 'mean_abundance']]
        
        # Use consistent colors
        colors = [cell_type_colors.get(ct, '#95A5A6') for ct in plot_data['cell_type']]
        
        # Create pie chart
        wedges, texts, autotexts = axes[i].pie(plot_data['mean_abundance'], 
                                              labels=plot_data['cell_type'],
                                              autopct='%1.1f%%',
                                              colors=colors)
        
        axes[i].set_title(f'{sample}\n({treatment})', fontsize=10, fontweight='bold')
        
        # Adjust text size
        for autotext in autotexts:
            autotext.set_fontsize(8)
        for text in texts:
            text.set_fontsize(8)

# Remove empty subplots
for i in range(len(adata_combined.obs['sample_id'].unique()), len(axes)):
    axes[i].remove()

plt.suptitle('Cell Type Composition by Sample', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(figure_path, "sample_composition_pie_charts.pdf"), 
           bbox_inches='tight', dpi=300)
plt.show()

# Stacked bar plot for treatment comparison with consistent colors
# Aggregate by treatment group and major cell categories
major_categories = {
    'Tumor': tumor_types,
    'Endothelial': ec_types, 
    'Immune': immune_types,  # Top 10 immune types
    'Stromal': stromal_types,
    'Epithelial': ['Epithelial'],
    'Hepatocyte': ['Hepatocyte']
}

treatment_composition = []
for treatment in adata_combined.obs['treatment_status'].unique():
    treatment_mask = adata_combined.obs['treatment_status'] == treatment
    treatment_data = adata_combined.obs[treatment_mask]
    
    for category, cell_types in major_categories.items():
        # Calculate mean abundance for this category
        available_types = [ct for ct in cell_types if ct in treatment_data.columns]
        if available_types:
            category_abundance = treatment_data[available_types].sum(axis=1).mean()
            treatment_composition.append({
                'treatment': treatment,
                'category': category,
                'abundance': category_abundance
            })

treatment_df = pd.DataFrame(treatment_composition)

# Create stacked bar plot with consistent colors
plt.figure(figsize=(10, 6))
treatment_pivot = treatment_df.pivot(index='treatment', columns='category', values='abundance').fillna(0)

# Use consistent category colors
category_color_list = [color_categories.get(cat, '#95A5A6') for cat in treatment_pivot.columns]

treatment_pivot.plot(kind='bar', stacked=True, ax=plt.gca(), 
                    color=category_color_list, figsize=(10, 6))
plt.title('Cell Type Composition by Treatment Group', fontsize=14, fontweight='bold')
plt.xlabel('Treatment Status', fontsize=12)
plt.ylabel('Mean Cell Abundance', fontsize=12)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(figure_path, "treatment_composition_stacked_bar.pdf"), 
           bbox_inches='tight', dpi=300)
plt.show()

# B) Tumor Subtype Focus Visualizations

# Bar plot: TC_Glycolysis vs TC_Quiescent by treatment with consistent colors
tumor_focus_data = []
for treatment in adata_combined.obs['treatment_status'].unique():
    treatment_mask = adata_combined.obs['treatment_status'] == treatment
    treatment_data = adata_combined.obs[treatment_mask]
    
    glycolysis_mean = treatment_data['TC_Glycolysis'].mean()
    quiescent_mean = treatment_data['TC_Quiescent'].mean()
    
    tumor_focus_data.extend([
        {'treatment': treatment, 'tumor_type': 'TC_Glycolysis', 'abundance': glycolysis_mean},
        {'treatment': treatment, 'tumor_type': 'TC_Quiescent', 'abundance': quiescent_mean}
    ])

tumor_focus_df = pd.DataFrame(tumor_focus_data)

# Use consistent colors for tumor types
tumor_type_colors = ['#C0392B', '#E74C3C']  # Glycolysis: dark red, Quiescent: red

plt.figure(figsize=(8, 6))
sns.barplot(data=tumor_focus_df, x='treatment', y='abundance', hue='tumor_type', 
           palette=tumor_type_colors)
plt.title('Key Tumor Subtypes by Treatment', fontsize=14, fontweight='bold')
plt.xlabel('Treatment Status', fontsize=12)
plt.ylabel('Mean Cell Abundance', fontsize=12)
plt.legend(title='Tumor Subtype')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(figure_path, "tumor_subtypes_by_treatment.pdf"), 
           bbox_inches='tight', dpi=300)
plt.show()

# Violin plots with statistical tests: Distribution differences
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# TC_Glycolysis
sns.violinplot(data=adata_combined.obs, x='treatment_status', y='TC_Glycolysis', 
              ax=axes[0], color='#C0392B')
axes[0].set_title('TC_Glycolysis Distribution', fontweight='bold')
axes[0].set_ylabel('Cell Abundance')

# Perform Wilcoxon test for TC_Glycolysis
stat, p_val = perform_wilcox_test(adata_combined.obs, 'treatment_status', 'TC_Glycolysis')
if p_val is not None:
    add_significance_bar(axes[0], 0, 1, axes[0].get_ylim()[1] * 0.9, p_val)

# TC_Quiescent
sns.violinplot(data=adata_combined.obs, x='treatment_status', y='TC_Quiescent', 
              ax=axes[1], color='#E74C3C')
axes[1].set_title('TC_Quiescent Distribution', fontweight='bold')
axes[1].set_ylabel('Cell Abundance')

# Perform Wilcoxon test for TC_Quiescent
stat, p_val = perform_wilcox_test(adata_combined.obs, 'treatment_status', 'TC_Quiescent')
if p_val is not None:
    add_significance_bar(axes[1], 0, 1, axes[1].get_ylim()[1] * 0.9, p_val)

plt.suptitle('Tumor Subtype Distributions by Treatment', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(figure_path, "tumor_subtype_distributions.pdf"), 
           bbox_inches='tight', dpi=300)
plt.show()

print("\nStep 1-2 Complete!")
print(f"Integrated {len(relevant_slides)} slides:")
print(f"  - Untreated: {sum(1 for t in treatment_mapping.values() if t == 'untreated')} samples")
print(f"  - Pre-chemotherapy: {sum(1 for t in treatment_mapping.values() if t == 'pre_chemotherapy')} samples")
print(f"Total spots: {adata_combined.n_obs}")
print(f"Available cell types: {len(all_cell_types)}")

# Save the integrated dataset for next steps
print("\nSaving integrated dataset...")
adata_combined.write_h5ad(os.path.join(figure_path, "integrated_dataset.h5ad"))
print("Integrated dataset saved as 'integrated_dataset.h5ad'")