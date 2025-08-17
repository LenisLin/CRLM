import os
import pickle
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from gseapy import Msigdb
import gseapy as gp

import scipy
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, pearsonr
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler

from Vascular_functions import *

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

#%% REFINED SPATIAL ANALYSIS: GLUT1+ METASTASIS-RESISTANCE HYPOTHESIS

# Define paths
run_name = '/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/analysis/cell2location_map'
figure_path = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/vessel_analysis"

if not os.path.exists(figure_path):
    os.makedirs(figure_path, exist_ok=True)

# Define slide information
all_slides = os.listdir(run_name)
relevant_slides = []
treatment_mapping = {}

for slide in all_slides:
    if slide.startswith('FDZS'):
        continue  # Exclude peritumor-focused slides
    elif 'Untreated' in slide:
        relevant_slides.append(slide)
        treatment_mapping[slide] = 'untreated'
    elif 'NACPR' in slide:
        relevant_slides.append(slide) 
        treatment_mapping[slide] = 'pre_chemotherapy'
    elif slide.startswith('GSE225857'):
        relevant_slides.append(slide)
        treatment_mapping[slide] = 'untreated'

print("GLUT1+ Metastasis-Resistance Analysis")
print("="*50)
print("Core Hypothesis: GLUT1+ cancer cells represent metastatic subtype")
print("that associates with vessels → chemotherapy resistance")
print(f"\nAnalyzing {len(relevant_slides)} slides:")
for slide in relevant_slides:
    print(f"  {slide} -> {treatment_mapping[slide]}")

#%% Main Analysis Execution

print("\n" + "="*60)
print("EXECUTING GLUT1+ METASTASIS-RESISTANCE ANALYSIS")
print("="*60)

# Store results for all slides
all_slide_results = {}
all_treatment_metrics = []
all_pathway_results = []

# Load pathway databases
vessel_pathways = load_pathway_databases()

# Run analysis for all slides
for slide_name in relevant_slides:
    print(f"\n{'='*20} {slide_name} {'='*20}")
    
    # Load slide data
    adata = load_slide_data(slide_name, run_name, treatment_mapping)
    
    # Create output directory
    slide_dir = create_slide_output_dir(figure_path, slide_name)
    
    # Step 2: Visualize the malignant distribution
    visualize_malignant_distribution(adata, slide_dir, slide_name)

    # Step 3: Setup tumor-focused vessel analysis
    vessel_success = setup_tumor_focused_vessel_analysis(adata, slide_name)
    
    ## Visualize the distribution of vessal and tumor region
    visualize_vascular_tumor_regions(adata, slide_dir, slide_name)

    # Step 4: Tumor-vessel spatial associations
    associations = analyze_tumor_vessel_spatial_associations(adata)
    distance_profiles = create_distance_profiles(adata, slide_name)
    visualize_tumor_vessel_associations(adata, slide_name, distance_profiles, associations, slide_dir)    

    # Step 5: Vascular hotspot co-occurrence (metastasis mechanism)
    cooccurrence_results = analyze_vascular_hotspot_cooccurrence(adata, slide_name)
    visualize_hotspot_cooccurrence(adata, slide_name, cooccurrence_results, slide_dir)
    
    # Step 6: Treatment metrics extraction
    treatment_metrics = analyze_treatment_effects_single_slide(
        adata, slide_name, associations, cooccurrence_results, distance_profiles
    )    
    all_treatment_metrics.append(treatment_metrics)
    
    # Step 8: Identify vessel-associated genes
    gene_correlations, vessel_associated_genes = identify_vessel_proximity_genes(adata, slide_name)
    
    # Step 9: Visualize vessel-gene associations
    visualize_vessel_gene_heatmap(adata, gene_correlations, vessel_associated_genes, slide_dir, slide_name)
        
    # Step 10: Pathway enrichment analysis
    enrichment_results, pathway_scores = perform_vessel_pathway_enrichment(adata, gene_correlations, vessel_pathways, slide_name)

    # Step 11: Visualize pathway enrichment
    visualize_pathway_enrichment(enrichment_results, pathway_scores, adata, slide_dir, slide_name)
    
    # Step 12: Cell type-pathway correlations
    correlation_results = analyze_cell_pathway_correlations(adata, pathway_scores, slide_name)
    
    visualize_cell_pathway_correlations(adata, correlation_results, slide_dir, slide_name)    

    # Store comprehensive results
    all_slide_results[slide_name] = {
        'slide_name': slide_name,
        'treatment': adata.obs['treatment_status'].iloc[0],
        'gene_correlations': gene_correlations,
        'vessel_associated_genes': vessel_associated_genes,
        'enrichment_results': enrichment_results,
        'pathway_scores': pathway_scores,
        'correlation_results': correlation_results
    }
    
    print(f"✓ Completed analysis for {slide_name}")

print(f"\n✅ Completed analysis for {len(all_slide_results)} slides")

#%% Treatment Comparison Analysis (Multi-slide level) - Modified for Target Analysis

print(f"\n--- TREATMENT COMPARISON ANALYSIS ---")

# Convert to DataFrame for analysis
treatment_df = pd.DataFrame(all_treatment_metrics)

print(f"Comparing {len(treatment_df)} slides:")
treatment_counts = treatment_df['treatment'].value_counts()
for treatment, count in treatment_counts.items():
    print(f"  {treatment}: {count} slides")

# Calculate tumor fraction for each slide
treatment_df['tumor_fraction'] = treatment_df['tumor_spots'] / treatment_df['n_spots']

#%% 1. TUMOR FRACTION CHANGE ANALYSIS
print(f"\n=== 1. TUMOR FRACTION CHANGE ANALYSIS ===")

# Compare tumor fractions between treatments
untreated_tumor_fraction = treatment_df[treatment_df['treatment'] == 'untreated']['tumor_fraction']
treated_tumor_fraction = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy']['tumor_fraction']

if len(untreated_tumor_fraction) > 0 and len(treated_tumor_fraction) > 0:
    untreated_mean = untreated_tumor_fraction.mean()
    treated_mean = treated_tumor_fraction.mean()
    
    fold_change = treated_mean / (untreated_mean + 1e-6)
    direction = "↑" if fold_change > 1.1 else "↓" if fold_change < 0.9 else "="
    
    print(f"Tumor Fraction Changes:")
    print(f"  Untreated: {untreated_mean:.3f} ± {untreated_tumor_fraction.std():.3f}")
    print(f"  Treated: {treated_mean:.3f} ± {treated_tumor_fraction.std():.3f}")
    print(f"  Change: {fold_change:.2f}x {direction}")
    
    # Statistical test
    stat, p_value = mannwhitneyu(untreated_tumor_fraction, treated_tumor_fraction, alternative='two-sided')
    significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
    print(f"  Statistical significance: p={p_value:.3f} {significance}")

#%% 2. RELATIVE TUMOR TYPES ANALYSIS
print(f"\n=== 2. RELATIVE TUMOR TYPES ANALYSIS ===")

# Get all tumor cell type abundance columns
tumor_type_cols = [col for col in treatment_df.columns if col.endswith('_mean_abundance')]
tumor_types = [col.replace('_mean_abundance', '').replace('TC_', '') for col in tumor_type_cols]

print(f"Analyzing tumor cell types: {tumor_types}")

# Calculate relative proportions for each slide
for idx, row in treatment_df.iterrows():
    total_abundance = sum([row[col] for col in tumor_type_cols])
    for col in tumor_type_cols:
        relative_col = col.replace('_mean_abundance', '_relative_proportion')
        treatment_df.loc[idx, relative_col] = row[col] / (total_abundance + 1e-6)

# Compare relative proportions between treatments
relative_cols = [col for col in treatment_df.columns if col.endswith('_relative_proportion')]

print(f"\nRelative Tumor Type Proportions:")
for col in relative_cols:
    tumor_type = col.replace('TC_', '').replace('_relative_proportion', '')
    
    untreated_values = treatment_df[treatment_df['treatment'] == 'untreated'][col]
    treated_values = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy'][col]
    
    if len(untreated_values) > 0 and len(treated_values) > 0:
        untreated_mean = untreated_values.mean()
        treated_mean = treated_values.mean()
        
        fold_change = treated_mean / (untreated_mean + 1e-6)
        direction = "↑" if fold_change > 1.2 else "↓" if fold_change < 0.8 else "="
        
        print(f"  {tumor_type}: {untreated_mean:.3f} -> {treated_mean:.3f} ({fold_change:.2f}x {direction})")

#%% 3. VESSEL PROXIMITY TYPES CHANGE ANALYSIS
print(f"\n=== 3. VESSEL PROXIMITY TYPES CHANGE ANALYSIS ===")

# Get vessel-related metrics
vessel_correlation_cols = [col for col in treatment_df.columns if col.endswith('_vessel_correlation')]
vessel_cooccurrence_cols = [col for col in treatment_df.columns if col.endswith('_cooccurrence_rate')]
vessel_enrichment_cols = [col for col in treatment_df.columns if col.endswith('_enrichment_ratio')]

print(f"Vessel Correlation Changes (negative = vessel-seeking):")
for col in vessel_correlation_cols:
    tumor_type = col.replace('TC_', '').replace('_vessel_correlation', '')
    
    untreated_values = treatment_df[treatment_df['treatment'] == 'untreated'][col]
    treated_values = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy'][col]
    
    if len(untreated_values) > 0 and len(treated_values) > 0:
        untreated_mean = untreated_values.mean()
        treated_mean = treated_values.mean()
        
        # For correlation, more negative = more vessel-seeking
        change = treated_mean - untreated_mean
        direction = "more vessel-seeking" if change < -0.1 else "less vessel-seeking" if change > 0.1 else "stable"
        
        print(f"  {tumor_type}: {untreated_mean:+.3f} -> {treated_mean:+.3f} ({direction})")

print(f"\nVessel Co-occurrence Rate Changes:")
for col in vessel_cooccurrence_cols:
    tumor_type = col.replace('TC_', '').replace('_cooccurrence_rate', '')
    
    untreated_values = treatment_df[treatment_df['treatment'] == 'untreated'][col]
    treated_values = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy'][col]
    
    if len(untreated_values) > 0 and len(treated_values) > 0:
        untreated_mean = untreated_values.mean()
        treated_mean = treated_values.mean()
        
        fold_change = treated_mean / (untreated_mean + 1e-6)
        direction = "↑" if fold_change > 1.2 else "↓" if fold_change < 0.8 else "="
        
        print(f"  {tumor_type}: {untreated_mean:.2f}% -> {treated_mean:.2f}% ({fold_change:.2f}x {direction})")

#%% COMPREHENSIVE VISUALIZATION
print(f"\n=== CREATING COMPREHENSIVE VISUALIZATIONS ===")

# Create large figure with custom layout
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.patch.set_facecolor('white')

# Define enhanced color palette
treatment_palette = {
    'untreated': '#2E8B57',      # Sea green
    'pre_chemotherapy': '#CD5C5C' # Indian red
}
accent_colors = ['#4A90E2', '#F5A623', '#D0021B', '#7ED321']

# Set consistent styling
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 10
})

# Panel A (Top Left): Enhanced Tumor Fraction Comparison
ax_a = axes[0, 0]
tumor_fraction_data = []
for _, row in treatment_df.iterrows():
    tumor_fraction_data.append({
        'slide': row['slide_name'],
        'treatment': row['treatment'],
        'tumor_fraction': row['tumor_fraction']
    })

tumor_fraction_df = pd.DataFrame(tumor_fraction_data)

# Enhanced violin plot with box plot overlay
parts = ax_a.violinplot([tumor_fraction_df[tumor_fraction_df['treatment'] == 'untreated']['tumor_fraction'],
                         tumor_fraction_df[tumor_fraction_df['treatment'] == 'pre_chemotherapy']['tumor_fraction']],
                        positions=[0, 1], widths=0.6, showmeans=True, showmedians=True)

# Customize violin plot colors
colors = [treatment_palette['untreated'], treatment_palette['pre_chemotherapy']]
for pc, color in zip(parts['bodies'], colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.7)
    pc.set_edgecolor('white')
    pc.set_linewidth(2)

# Add individual data points
for i, treatment in enumerate(['untreated', 'pre_chemotherapy']):
    data = tumor_fraction_df[tumor_fraction_df['treatment'] == treatment]['tumor_fraction']
    y = data.values
    x = np.random.normal(i, 0.04, size=len(y))  # Add some jitter
    ax_a.scatter(x, y, alpha=0.8, s=60, color='white', edgecolor=colors[i], linewidth=2, zorder=3)

# Statistical annotation
untreated_vals = tumor_fraction_df[tumor_fraction_df['treatment'] == 'untreated']['tumor_fraction']
treated_vals = tumor_fraction_df[tumor_fraction_df['treatment'] == 'pre_chemotherapy']['tumor_fraction']

if len(untreated_vals) > 0 and len(treated_vals) > 0:
    stat, p_val = mannwhitneyu(untreated_vals, treated_vals)
    significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
    
    y_max = tumor_fraction_df['tumor_fraction'].max()
    ax_a.plot([0, 1], [y_max*1.05, y_max*1.05], 'k-', linewidth=1.5)
    ax_a.text(0.5, y_max*1.08, f'{significance}', ha='center', va='bottom', 
             fontweight='bold', fontsize=12)
    ax_a.text(0.5, y_max*1.12, f'p={p_val:.3f}', ha='center', va='bottom', 
             fontsize=9, style='italic')

ax_a.set_title('A. Tumor Burden Response', fontweight='bold', pad=20)
ax_a.set_ylabel('Tumor Fraction', fontweight='bold')
ax_a.set_xticks([0, 1])
ax_a.set_xticklabels(['Untreated', 'Pre-Chemotherapy'], fontweight='bold')
ax_a.grid(True, alpha=0.3, axis='y', linestyle='--')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Panel B (Top Right): Enhanced Tumor Composition Heatmap
ax_b = axes[0, 1]
relative_data = treatment_df.groupby('treatment')[relative_cols].mean()
relative_data.columns = [col.replace('TC_', '').replace('_relative_proportion', '') 
                        for col in relative_data.columns]

# Calculate statistics for annotations
stats_annotations = []
for col in relative_cols:
    untreated_vals = treatment_df[treatment_df['treatment'] == 'untreated'][col]
    treated_vals = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy'][col]
    
    if len(untreated_vals) > 0 and len(treated_vals) > 0:
        stat, p_val = mannwhitneyu(untreated_vals, treated_vals)
        fold_change = treated_vals.mean() / (untreated_vals.mean() + 1e-6)
        significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        stats_annotations.append(significance)
    else:
        stats_annotations.append("")

# Create annotation matrix
annot_matrix = relative_data.copy()
for i, (idx, row) in enumerate(relative_data.iterrows()):
    for j, val in enumerate(row):
        if i == 1:  # Pre-chemotherapy row
            annot_matrix.iloc[i, j] = f'{val:.3f}\n{stats_annotations[j]}'
        else:
            annot_matrix.iloc[i, j] = f'{val:.3f}'

# Enhanced heatmap
sns.heatmap(relative_data.T, annot=annot_matrix.T, fmt='', cmap='RdYlBu_r', 
            center=relative_data.values.mean(), ax=ax_b, 
            cbar_kws={'label': 'Relative Proportion', 'shrink': 0.8},
            linewidths=1, linecolor='white')

ax_b.set_title('B. Tumor Composition Changes', fontweight='bold', pad=20)
ax_b.set_xlabel('Treatment Status', fontweight='bold')
ax_b.set_ylabel('Tumor Cell Type', fontweight='bold')

# Panel C (Bottom Left): Enhanced Co-occurrence Analysis
ax_c = axes[1, 0]
cooccur_data = []
for col in vessel_cooccurrence_cols[:4]:
    tumor_type = col.replace('TC_', '').replace('_cooccurrence_rate', '')
    for _, row in treatment_df.iterrows():
        cooccur_data.append({
            'tumor_type': tumor_type,
            'treatment': row['treatment'],
            'cooccurrence_rate': row[col]
        })

cooccur_df = pd.DataFrame(cooccur_data)

# Enhanced grouped bar plot with error bars
sns.barplot(data=cooccur_df, x='tumor_type', y='cooccurrence_rate', 
            hue='treatment', palette=treatment_palette, ax=ax_c,
            capsize=0.1, errwidth=1.5, alpha=0.8)

# Add individual data points
for i, tumor_type in enumerate(cooccur_df['tumor_type'].unique()):
    for j, treatment in enumerate(['untreated', 'pre_chemotherapy']):
        data = cooccur_df[(cooccur_df['tumor_type'] == tumor_type) & 
                         (cooccur_df['treatment'] == treatment)]['cooccurrence_rate']
        if len(data) > 0:
            x_pos = i + (j - 0.5) * 0.4  # Offset for grouped bars
            ax_c.scatter([x_pos] * len(data), data.values, 
                        color='white', edgecolor='black', s=30, alpha=0.8, zorder=3)

ax_c.set_title('C. Vascular Co-occurrence Patterns', fontweight='bold', pad=20)
ax_c.set_ylabel('Co-occurrence Rate (%)', fontweight='bold')
ax_c.set_xlabel('Tumor Cell Type', fontweight='bold')
ax_c.tick_params(axis='x', rotation=45)
ax_c.legend(title='Treatment', frameon=True, fancybox=True, shadow=True)
ax_c.grid(True, alpha=0.3, axis='y', linestyle='--')
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# Panel D (Bottom Right): Enhanced Treatment Effect Summary
ax_d = axes[1, 1]

# Calculate key metrics with effect sizes
summary_metrics = ['tumor_fraction'] + relative_cols[:3]
summary_data = []

for metric in summary_metrics:
    untreated_vals = treatment_df[treatment_df['treatment'] == 'untreated'][metric]
    treated_vals = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy'][metric]
    
    if len(untreated_vals) > 0 and len(treated_vals) > 0:
        fold_change = treated_vals.mean() / (untreated_vals.mean() + 1e-6)
        stat, p_val = mannwhitneyu(untreated_vals, treated_vals)
        
        # Calculate Cohen's d for effect size
        pooled_std = np.sqrt(((len(untreated_vals) - 1) * untreated_vals.std()**2 + 
                             (len(treated_vals) - 1) * treated_vals.std()**2) / 
                            (len(untreated_vals) + len(treated_vals) - 2))
        cohens_d = abs(treated_vals.mean() - untreated_vals.mean()) / pooled_std
        
        summary_data.append({
            'metric': metric.replace('TC_', '').replace('_relative_proportion', '').replace('_', ' '),
            'fold_change': fold_change,
            'p_value': p_val,
            'effect_size': cohens_d,
            'significant': p_val < 0.05
        })

if summary_data:
    summary_df = pd.DataFrame(summary_data)
    
    # Enhanced color scheme based on significance and effect size
    colors = []
    for _, row in summary_df.iterrows():
        if row['significant']:
            if row['effect_size'] > 0.8:  # Large effect
                color = '#D32F2F' if row['fold_change'] > 1 else '#1976D2'
            elif row['effect_size'] > 0.5:  # Medium effect
                color = '#F57C00' if row['fold_change'] > 1 else '#0288D1'
            else:  # Small effect
                color = '#FFA726' if row['fold_change'] > 1 else '#03A9F4'
        else:
            color = '#9E9E9E'  # Gray for non-significant
        colors.append(color)
    
    bars = ax_d.bar(range(len(summary_df)), summary_df['fold_change'], 
                    color=colors, alpha=0.8, edgecolor='white', linewidth=2)
    
    # Add significance and effect size annotations
    for i, (bar, row) in enumerate(zip(bars, summary_df.itertuples())):
        height = bar.get_height()
        
        # Significance stars
        significance = "***" if row.p_value < 0.001 else "**" if row.p_value < 0.01 else "*" if row.p_value < 0.05 else ""
        if significance:
            ax_d.text(i, height + 0.05, significance, ha='center', va='bottom', 
                     fontweight='bold', fontsize=12)
        
        # Effect size indicator
        effect_text = f'd={row.effect_size:.2f}'
        ax_d.text(i, height/2, effect_text, ha='center', va='center', 
                 fontweight='bold', fontsize=8, color='white',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.6))
    
    ax_d.axhline(y=1, color='black', linestyle='-', alpha=0.8, linewidth=2)
    ax_d.set_title('D. Treatment Effect Summary', fontweight='bold', pad=20)
    ax_d.set_ylabel('Fold Change (Treated/Untreated)', fontweight='bold')
    ax_d.set_xlabel('Metric', fontweight='bold')
    ax_d.set_xticks(range(len(summary_df)))
    ax_d.set_xticklabels(summary_df['metric'], rotation=45, ha='right')
    ax_d.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax_d.spines['top'].set_visible(False)
    ax_d.spines['right'].set_visible(False)

# Overall styling enhancements
plt.suptitle('Spatial Tumor Dynamics: Treatment Response Analysis', 
             fontsize=16, fontweight='bold', y=0.98)

# Add a subtle background color
fig.patch.set_facecolor('#FAFAFA')

# Adjust spacing
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Add panel borders
for ax in axes.flat:
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('#333333')

plt.savefig(os.path.join(figure_path, "enhanced_four_panel_treatment_analysis.pdf"), 
            bbox_inches='tight', dpi=300, facecolor='white', edgecolor='none')
plt.show()

# Save extended results with new metrics
final_results = {
    'slide_results': all_slide_results,
    'treatment_metrics': all_treatment_metrics,
    'treatment_summary': treatment_df.to_dict('records') if len(treatment_df) > 0 else [],
    'tumor_fraction_analysis': {
        'untreated_mean': untreated_tumor_fraction.mean() if len(untreated_tumor_fraction) > 0 else None,
        'treated_mean': treated_tumor_fraction.mean() if len(treated_tumor_fraction) > 0 else None,
        'statistical_test': {'statistic': stat, 'p_value': p_value} if 'stat' in locals() else None
    },
    'analysis_parameters': {
        'perivascular_threshold': 55,
        'intermediate_threshold': 165,
        'tumor_threshold_percentile': 25,
        'vascular_hotspot_percentile': 95
    }
}
# Save results
with open(os.path.join(figure_path, "glut1_metastasis_resistance_results.pkl"), 'wb') as f:
    pickle.dump(final_results, f)

if len(treatment_df) > 0:
    treatment_df.to_csv(os.path.join(figure_path, "treatment_comparison_summary.csv"), index=False)

final_pathway_results = {
    'individual_slides': all_pathway_results,
    'pathway_databases': vessel_pathways,
    'analysis_parameters': PATHWAY_ANALYSIS_PARAMS,
    'methods_summary': {
        'step1': 'Vessel-proximity gene identification via distance correlation',
        'step2': 'Pathway enrichment analysis using MSigDB Hallmark pathways',
        'step3': 'Cell type-pathway correlation analysis in perivascular regions',
        'clinical_relevance': 'Mechanistic insights for vessel-targeting therapeutics'
    }
}

with open(os.path.join(figure_path, "vessel_pathway_analysis_results.pkl"), 'wb') as f:
    pickle.dump(final_pathway_results, f)

print(f"\n✅ VESSEL PATHWAY ANALYSIS COMPLETE!")
print(f"📁 Results saved in: {figure_path}")
print(f"🧬 Pathway mechanisms identified for vessel-tumor interactions")
print(f"🎯 Ready for therapeutic target validation!")
