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
sc.settings.set_figure_params(dpi=80, facecolor='white')

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

# Define EMT-associated genes (replacing poor TC_EMT definition)
EMT_GENES = ['VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'ZEB2', 'TWIST1', 'CDH2', 'FN1']

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

    # Calculate EMT score from genes
    calculate_emt_score(adata, EMT_GENES)
    
    # Create output directory
    slide_dir = create_slide_output_dir(figure_path, slide_name)
    
    # Step 2: Visualize the malignant distribution
    visualize_malignant_distribution(adata, slide_dir, slide_name)

    # Step 3: Setup tumor-focused vessel analysis
    vessel_success = setup_tumor_focused_vessel_analysis(adata, slide_name)
    
    # Step 4: Tumor-vessel spatial associations
    associations = analyze_tumor_vessel_spatial_associations(adata, slide_name)
    distance_profiles = create_distance_profiles(adata, slide_name)
    visualize_tumor_vessel_associations(adata, slide_name, distance_profiles, slide_dir)
    
    # Step 5: Vascular hotspot co-occurrence (metastasis mechanism)
    cooccurrence_results = analyze_vascular_hotspot_cooccurrence(adata, slide_name)
    visualize_hotspot_cooccurrence(adata, slide_name, cooccurrence_results, slide_dir)
    
    # Step 6: Treatment metrics extraction
    treatment_metrics = analyze_treatment_effects_single_slide(adata, slide_name, associations, cooccurrence_results)
    all_treatment_metrics.append(treatment_metrics)
    
    # Step 7: Spatial competition and EMT analysis
    pairwise_results, emt_results = analyze_spatial_competition_emt(adata, slide_name)
    visualize_competition_emt(adata, slide_name, pairwise_results, emt_results, slide_dir)

    # Step 8: Identify vessel-associated genes
    gene_correlations, vessel_associated_genes = identify_vessel_proximity_genes(adata, slide_name)
    
    if len(gene_correlations) != 0:

        # Step 9: Visualize vessel-gene associations
        visualize_vessel_gene_heatmap(adata, gene_correlations, vessel_associated_genes, slide_dir, slide_name)
            
        # Step 10: Pathway enrichment analysis
        enrichment_results, pathway_scores = perform_vessel_pathway_enrichment(
            adata, gene_correlations, vessel_pathways, slide_name)
    
        # Step 11: Visualize pathway enrichment
        visualize_pathway_enrichment(enrichment_results, pathway_scores, adata, slide_dir, slide_name)
        
        # Step 12: Cell type-pathway correlations
        correlation_results = analyze_cell_pathway_correlations(adata, pathway_scores, slide_name)
        
        if correlation_results:
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
fig = plt.figure(figsize=(20, 14))
gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1], width_ratios=[1, 1, 1, 0.8])

# Define consistent colors
treatment_palette = {'untreated': '#2ECC71', 'pre_chemotherapy': '#E74C3C'}

# 1. Enhanced Tumor Fraction Comparison
ax1 = fig.add_subplot(gs[0, 0])
tumor_fraction_data = []
for _, row in treatment_df.iterrows():
    tumor_fraction_data.append({
        'slide': row['slide_name'],
        'treatment': row['treatment'],
        'tumor_fraction': row['tumor_fraction']
    })

tumor_fraction_df = pd.DataFrame(tumor_fraction_data)

# Enhanced boxplot with swarm overlay
sns.boxplot(data=tumor_fraction_df, x='treatment', y='tumor_fraction', 
            palette=treatment_palette, ax=ax1, width=0.6)
sns.swarmplot(data=tumor_fraction_df, x='treatment', y='tumor_fraction', 
                color='white', edgecolor='black', linewidth=1, size=8, ax=ax1)

# Statistical annotation
untreated_vals = tumor_fraction_df[tumor_fraction_df['treatment'] == 'untreated']['tumor_fraction']
treated_vals = tumor_fraction_df[tumor_fraction_df['treatment'] == 'pre_chemotherapy']['tumor_fraction']

if len(untreated_vals) > 0 and len(treated_vals) > 0:
    stat, p_val = mannwhitneyu(untreated_vals, treated_vals)
    significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
    
    # Add statistical annotation
    y_max = tumor_fraction_df['tumor_fraction'].max()
    ax1.plot([0, 1], [y_max*1.1, y_max*1.1], 'k-', linewidth=1)
    ax1.text(0.5, y_max*1.12, f'p={p_val:.3f} {significance}', 
            ha='center', va='bottom', fontweight='bold')

ax1.set_title('Tumor Burden Response', fontsize=14, fontweight='bold', pad=15)
ax1.set_ylabel('Tumor Fraction', fontsize=12, fontweight='bold')
ax1.set_xlabel('Treatment Status', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y')

# 2. Enhanced Relative Tumor Types Heatmap
ax2 = fig.add_subplot(gs[0, 1:3])
relative_data = treatment_df.groupby('treatment')[relative_cols].mean()
relative_data.columns = [col.replace('TC_', '').replace('_relative_proportion', '') 
                        for col in relative_data.columns]

# Calculate fold changes for annotation
fold_changes = relative_data.loc['pre_chemotherapy'] / relative_data.loc['untreated']

# Enhanced heatmap
sns.heatmap(relative_data.T, annot=True, fmt='.3f', cmap='RdYlBu_r', 
            center=relative_data.values.mean(), ax=ax2, 
            cbar_kws={'label': 'Relative Proportion'})

# Add fold change annotations
for i, tumor_type in enumerate(relative_data.columns):
    fc = fold_changes[tumor_type]
    if fc > 1.2:
        arrow = '↗'
        color = 'red'
    elif fc < 0.8:
        arrow = '↘'
        color = 'blue'
    else:
        arrow = '→'
        color = 'gray'
    
    ax2.text(1.5, i + 0.5, f'{arrow} {fc:.2f}x', 
            ha='center', va='center', fontweight='bold', 
            color=color, fontsize=11)

ax2.set_title('Tumor Composition Changes\n(Arrows: Treatment Effect)', 
                fontsize=14, fontweight='bold', pad=15)
ax2.set_xlabel('Treatment Status', fontsize=12, fontweight='bold')

# 3. Enhanced Vessel Correlation Matrix
ax3 = fig.add_subplot(gs[1, :2])
vessel_corr_data = []
for col in vessel_correlation_cols[:4]:
    tumor_type = col.replace('TC_', '').replace('_vessel_correlation', '')
    for _, row in treatment_df.iterrows():
        vessel_corr_data.append({
            'tumor_type': tumor_type,
            'treatment': row['treatment'],
            'vessel_correlation': row[col]
        })

vessel_corr_df = pd.DataFrame(vessel_corr_data)

# Create matrix for heatmap
aggregated_df = vessel_corr_df.groupby(['tumor_type', 'treatment']).mean().reset_index()
corr_matrix = aggregated_df.pivot(index='tumor_type', columns='treatment', values='vessel_correlation')

# Enhanced heatmap with diverging colormap
sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
            ax=ax3, cbar_kws={'label': 'Vessel Correlation\n(-: seeking, +: avoiding)'})

ax3.set_title('Vessel-Seeking Behavior by Treatment\n(Red: Vessel-Avoiding, Blue: Vessel-Seeking)', 
                fontsize=14, fontweight='bold', pad=15)
ax3.set_xlabel('Treatment Status', fontsize=12, fontweight='bold')
ax3.set_ylabel('Tumor Cell Type', fontsize=12, fontweight='bold')

# 4. Enhanced Co-occurrence Comparison
ax4 = fig.add_subplot(gs[1, 2:])
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

# Enhanced grouped bar plot
sns.barplot(data=cooccur_df, x='tumor_type', y='cooccurrence_rate', 
            hue='treatment', palette=treatment_palette, ax=ax4)

ax4.set_title('Metastatic Potential by Treatment\n(Co-occurrence Rate)', 
                fontsize=14, fontweight='bold', pad=15)
ax4.set_ylabel('Co-occurrence Rate (%)', fontsize=12, fontweight='bold')
ax4.set_xlabel('Tumor Cell Type', fontsize=12, fontweight='bold')
ax4.tick_params(axis='x', rotation=45)
ax4.legend(title='Treatment', fontsize=11, title_fontsize=12)
ax4.grid(True, alpha=0.3, axis='y')

# 5. Enhanced Summary Statistics
ax5 = fig.add_subplot(gs[2, :2])

# Calculate key metrics
summary_metrics = ['tumor_fraction'] + relative_cols[:3]
summary_data = []

for metric in summary_metrics:
    untreated_vals = treatment_df[treatment_df['treatment'] == 'untreated'][metric]
    treated_vals = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy'][metric]
    
    if len(untreated_vals) > 0 and len(treated_vals) > 0:
        fold_change = treated_vals.mean() / (untreated_vals.mean() + 1e-6)
        stat, p_val = mannwhitneyu(untreated_vals, treated_vals)
        
        summary_data.append({
            'metric': metric.replace('TC_', '').replace('_relative_proportion', '').replace('_', ' '),
            'fold_change': fold_change,
            'p_value': p_val,
            'significant': p_val < 0.05
        })

if summary_data:
    summary_df = pd.DataFrame(summary_data)
    
    # Color based on significance and direction
    colors = []
    for _, row in summary_df.iterrows():
        if row['significant']:
            if row['fold_change'] > 1.2:
                colors.append('#E74C3C')  # Red for significant increase
            elif row['fold_change'] < 0.8:
                colors.append('#3498DB')  # Blue for significant decrease
            else:
                colors.append('#F39C12')  # Orange for significant but small change
        else:
            colors.append('#95A5A6')  # Gray for non-significant
    
    bars = ax5.bar(summary_df['metric'], summary_df['fold_change'], color=colors,
                    edgecolor='white', linewidth=1.5, alpha=0.8)
    
    # Add significance annotations
    for bar, row in zip(bars, summary_df.itertuples()):
        height = bar.get_height()
        significance = "***" if row.p_value < 0.001 else "**" if row.p_value < 0.01 else "*" if row.p_value < 0.05 else "ns"
        ax5.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                significance, ha='center', va='bottom', fontweight='bold', fontsize=12)
    
    ax5.axhline(y=1, color='black', linestyle='-', alpha=0.8, linewidth=2)
    ax5.set_title('Treatment Effect Summary\n(Fold Change: Treated/Untreated)', 
                    fontsize=14, fontweight='bold', pad=15)
    ax5.set_ylabel('Fold Change', fontsize=12, fontweight='bold')
    ax5.set_xlabel('Metric', fontsize=12, fontweight='bold')
    ax5.tick_params(axis='x', rotation=45)
    ax5.grid(True, alpha=0.3, axis='y')

# 6. Enhanced Key Findings Panel
ax6 = fig.add_subplot(gs[2, 2:])
ax6.axis('off')

# Create enhanced summary text
summary_text = "KEY TREATMENT EFFECTS\n" + "="*25 + "\n\n"

# Tumor fraction effect
if len(untreated_vals) > 0 and len(treated_vals) > 0:
    tf_change = treated_vals.mean() / untreated_vals.mean()
    tf_direction = "↗ INCREASED" if tf_change > 1.1 else "↘ DECREASED" if tf_change < 0.9 else "→ STABLE"
    summary_text += f"🎯 TUMOR BURDEN: {tf_direction}\n"
    summary_text += f"   ({tf_change:.2f}x change)\n\n"

# Top changing tumor types
fold_changes = {}
for col in relative_cols:
    tumor_type = col.replace('TC_', '').replace('_relative_proportion', '')
    untreated_vals = treatment_df[treatment_df['treatment'] == 'untreated'][col]
    treated_vals = treatment_df[treatment_df['treatment'] == 'pre_chemotherapy'][col]
    
    if len(untreated_vals) > 0 and len(treated_vals) > 0:
        fold_changes[tumor_type] = treated_vals.mean() / (untreated_vals.mean() + 1e-6)

# Sort by fold change magnitude
sorted_changes = sorted(fold_changes.items(), key=lambda x: abs(x[1] - 1), reverse=True)

summary_text += "🔄 TUMOR COMPOSITION:\n"
for tumor_type, fc in sorted_changes[:3]:
    if fc > 1.2:
        direction = "↗ ENRICHED"
        color_indicator = "🔴"
    elif fc < 0.8:
        direction = "↘ DEPLETED"
        color_indicator = "🔵"
    else:
        direction = "→ STABLE"
        color_indicator = "⚪"
    
    summary_text += f"   {color_indicator} {tumor_type}: {direction} ({fc:.2f}x)\n"

summary_text += f"\n🩸 VESSEL ASSOCIATIONS:\n"
summary_text += f"   See correlation matrix →\n"
summary_text += f"   Red = Vessel-avoiding\n"
summary_text += f"   Blue = Vessel-seeking\n\n"

summary_text += f"📊 METASTATIC POTENTIAL:\n"
summary_text += f"   See co-occurrence plots →"

# Enhanced text styling
ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, 
        fontsize=11, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle="round,pad=0.5", facecolor='lightgray', alpha=0.3))

# Overall title
plt.suptitle('Comprehensive Treatment Response Analysis\nSpatial Tumor Dynamics & Metastatic Potential', 
                fontsize=18, fontweight='bold', y=0.98)

plt.tight_layout()
plt.savefig(os.path.join(figure_path, "enhanced_comprehensive_treatment_analysis.pdf"), 
            bbox_inches='tight', dpi=300, facecolor='white')
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
        'emt_genes': EMT_GENES,
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
