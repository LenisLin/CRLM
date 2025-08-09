import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib import patheffects
from gseapy import Msigdb
import gseapy as gp

import seaborn as sns
from pathlib import Path
import os

import scipy
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, pearsonr
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler

from scipy import stats

# Set publication-quality style
plt.style.use('default')
sns.set_palette("Set2")

# Define enhanced color palettes
VESSEL_SEEKING_COLOR = '#E74C3C'  # Red
VESSEL_AVOIDING_COLOR = '#3498DB'  # Blue  
NEUTRAL_COLOR = '#95A5A6'  # Gray
GLUT1_COLOR = '#E67E22'  # Orange
HIGH_RISK_COLOR = '#8E44AD'  # Purple
TREATMENT_COLORS = {'untreated': '#2ECC71', 'pre_chemotherapy': '#F39C12'}

# Define pathway analysis parameters
PATHWAY_ANALYSIS_PARAMS = {
    'distance_correlation_threshold': 0.25,  # Minimum correlation with vessel distance
    'perivascular_zone': 100,  # μm for pathway enrichment analysis
    'intermediate_zone': 200,  # μm
    'min_spots_for_analysis': 10,  # Minimum spots for reliable analysis
    'pathway_score_method': 'ssgsea',  # or 'ssgsea'
    'enrichment_pvalue_threshold': 0.05
}

#%% Highly flexible reusable functions

def load_slide_data(slide_name, run_name, treatment_mapping):
    """Load and prepare slide data"""
    slide_path = os.path.join(run_name, slide_name, "sp.h5ad")
    
    if not os.path.exists(slide_path):
        print(f"Warning: File not found - {slide_path}")
        return None
    
    adata = sc.read_h5ad(slide_path)
    adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
    adata.obs['sample_id'] = slide_name
    adata.obs['treatment_status'] = treatment_mapping[slide_name]
    
    print(f"Loaded {slide_name}: {adata.n_obs} spots, {adata.n_vars} genes")
    return adata

def calculate_emt_score(adata, emt_genes):
    """Calculate EMT score from gene expression"""
    available_genes = adata.var.index.tolist()
    available_emt_genes = [g for g in emt_genes if g in available_genes]
    
    if len(available_emt_genes) >= 2:
        emt_expr = adata[:, available_emt_genes].X
        if hasattr(emt_expr, 'toarray'):
            emt_expr = emt_expr.toarray()
        emt_score = np.mean(emt_expr, axis=1)
        adata.obs['EMT_Score'] = emt_score
        print(f"  EMT score calculated from {len(available_emt_genes)} genes: {available_emt_genes}")
        return True
    else:
        print(f"  Insufficient EMT genes ({len(available_emt_genes)} available)")
        adata.obs['EMT_Score'] = 0
        return False

def create_slide_output_dir(figure_path, slide_name):
    """Create slide-specific output directory"""
    slide_dir = os.path.join(figure_path, slide_name)
    if not os.path.exists(slide_dir):
        os.makedirs(slide_dir, exist_ok=True)
    return slide_dir

#%% Step 2: Visualize the malignant distribution
def visualize_malignant_distribution(adata, slide_dir, slide_name):
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]

    with mpl.rc_context({'axes.facecolor':  'black','figure.figsize': [4.5, 5]}):
        sc.pl.spatial(adata, cmap='magma',
                    # show first 8 cell types
                    color=tumor_types,
                    ncols=4, size=1.3,
                    img_key='hires',
                    # limit color scale at 99.2% quantile of cell abundance
                    vmin=0, vmax='p99.2', show = False
                    )
        
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"tumor_spatial_distribution_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()

    return None

#%% Step 3: Spatial Foundation Setup (Tumor-focused vessel analysis)

def setup_tumor_focused_vessel_analysis(adata, slide_name):
    """Setup vessel analysis focused on tumor regions only"""
    
    print(f"\n--- Step 3: Tumor-focused vessel analysis ({slide_name}) ---")
    
    # Get endothelial and tumor cell types
    ec_types = [col for col in adata.obs.columns if col.startswith('EC_')]
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    
    if len(ec_types) == 0 or len(tumor_types) == 0:
        print(f"  Missing cell types: EC={len(ec_types)}, TC={len(tumor_types)}")
        return False
    
    # Calculate vascular density
    adata.obs['Vascular_Sum'] = adata.obs[ec_types].sum(axis=1)
    
    # Calculate total tumor cell abundance to define tumor regions
    adata.obs['Tumor_Sum'] = adata.obs[tumor_types].sum(axis=1)
    
    # Define tumor regions (spots with significant tumor cell presence)
    # tumor_threshold = np.percentile(adata.obs['Tumor_Sum'], 95)
    tumor_threshold = 1.5
    adata.obs['Is_Tumor_Region'] = adata.obs['Tumor_Sum'] >= tumor_threshold

    n_tumor_spots = adata.obs['Is_Tumor_Region'].sum()
    tumor_percentage = n_tumor_spots / len(adata.obs) * 100
    
    print(f"  Tumor regions: {n_tumor_spots} spots ({tumor_percentage:.1f}%)")
    
    # Focus vessel analysis on tumor regions only
    tumor_mask = adata.obs['Is_Tumor_Region']
    tumor_vascular = adata.obs.loc[tumor_mask, 'Vascular_Sum']
    
    if len(tumor_vascular) == 0:
        print(f"  No tumor regions found")
        return False
    
    # Define vascular hotspots within tumor regions
    vascular_threshold = np.percentile(tumor_vascular, 95)
    adata.obs['Is_Vascular_Hotspot'] = False
    tumor_hotspot_mask = tumor_mask & (adata.obs['Vascular_Sum'] >= vascular_threshold)
    adata.obs.loc[tumor_hotspot_mask, 'Is_Vascular_Hotspot'] = True
    
    n_hotspots = adata.obs['Is_Vascular_Hotspot'].sum()
    hotspot_percentage = n_hotspots / n_tumor_spots * 100 if n_tumor_spots > 0 else 0
    
    print(f"  Vascular hotspots in tumor: {n_hotspots} spots ({hotspot_percentage:.1f}% of tumor)")
    
    # Calculate distances to vessels within tumor regions only
    coords = adata.obs.loc[:,['array_row','array_col']]
    hotspot_coords = coords[adata.obs['Is_Vascular_Hotspot']]
    
    if len(hotspot_coords) > 0:
        distance_matrix = cdist(coords, hotspot_coords)
        nearest_distances = distance_matrix.min(axis=1)
        adata.obs['Distance_to_Vessel_um'] = nearest_distances * 100 
        
        # Define spatial zones within tumor regions
        perivascular_max = 100    # 0-55μm: immediate perivascular
        intermediate_max = 200   # 55-165μm: intermediate zone
        
        conditions = [
            adata.obs['Distance_to_Vessel_um'] <= perivascular_max,
            (adata.obs['Distance_to_Vessel_um'] > perivascular_max) & 
            (adata.obs['Distance_to_Vessel_um'] <= intermediate_max),
            adata.obs['Distance_to_Vessel_um'] > intermediate_max
        ]
        
        choices = ['Perivascular', 'Intermediate', 'Distant']
        adata.obs['Spatial_Zone'] = np.select(conditions, choices, default='Non_Tumor')
        
        # Report zones within tumor regions only
        tumor_zones = adata.obs[tumor_mask]['Spatial_Zone'].value_counts()
        print(f"  Spatial zones in tumor:")
        for zone in ['Perivascular', 'Intermediate', 'Distant']:
            if zone in tumor_zones:
                count = tumor_zones[zone]
                percentage = count / n_tumor_spots * 100 if n_tumor_spots > 0 else 0
                print(f"    {zone}: {count} spots ({percentage:.1f}%)")
        
        return True
    else:
        print(f"  No vascular hotspots found in tumor regions")
        return False

#%% Step 4: Vascular-Tumor Subtype Spatial Associations

def analyze_tumor_vessel_spatial_associations(adata, slide_name):
    """Analyze spatial associations between tumor subtypes and vessels"""
    
    print(f"  Analyzing tumor-vessel spatial associations...")
    
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    tumor_mask = adata.obs['Is_Tumor_Region']
    
    if tumor_mask.sum() == 0:
        return {}
    
    associations = {}
    
    # Focus analysis on tumor regions only
    tumor_data = adata.obs[tumor_mask]
    
    print(f"    Distance correlations in tumor regions:")
    
    for tc_type in tumor_types:
        # Distance correlation (negative = prefers vessels)
        valid_mask = np.isfinite(tumor_data['Distance_to_Vessel_um'])
        
        if valid_mask.sum() > 10:
            distances = tumor_data.loc[valid_mask, 'Distance_to_Vessel_um']
            abundances = tumor_data.loc[valid_mask, tc_type]
            
            corr = np.corrcoef(distances, abundances)[0, 1]
            
            # Zone-specific abundances (tumor regions only)
            zone_abundances = {}
            for zone in ['Perivascular', 'Intermediate', 'Distant']:
                zone_mask = tumor_data['Spatial_Zone'] == zone
                if zone_mask.sum() > 0:
                    zone_abundances[zone] = tumor_data.loc[zone_mask, tc_type].mean()
            
            associations[tc_type] = {
                'distance_correlation': corr,
                'zone_abundances': zone_abundances
            }
            
            # Interpretation
            vessel_association = "vessel-seeking" if corr < -0.2 else \
                               "vessel-avoiding" if corr > 0.2 else "neutral"
            
            print(f"      {tc_type}: r={corr:+.3f} ({vessel_association})")
    
    return associations

def create_distance_profiles(adata, slide_name):
    """Create distance-dependent association curves for tumor subtypes"""
    
    print(f"  Creating distance profiles for {slide_name}...")
    
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    tumor_mask = adata.obs['Is_Tumor_Region']
    tumor_data = adata.obs[tumor_mask]
    
    # Create distance bins (55μm resolution)
    distance_bins = np.arange(0, 400, 55)  # Up to 400μm
    bin_centers = distance_bins[:-1] + 27.5
    
    distance_profiles = {}
    
    for tc_type in tumor_types:
        bin_profiles = []
        bin_distances = []
        
        for bin_center in bin_centers:
            bin_min = bin_center - 27.5
            bin_max = bin_center + 27.5
            
            bin_mask = (tumor_data['Distance_to_Vessel_um'] >= bin_min) & \
                      (tumor_data['Distance_to_Vessel_um'] < bin_max)
            
            if bin_mask.sum() >= 3:  # Minimum spots for reliable statistics
                mean_abundance = tumor_data.loc[bin_mask, tc_type].mean()
                bin_profiles.append(mean_abundance)
                bin_distances.append(bin_center)
        
        if len(bin_profiles) > 2:
            distance_profiles[tc_type] = {
                'distances': np.array(bin_distances),
                'abundances': np.array(bin_profiles)
            }
    
    return distance_profiles

def visualize_tumor_vessel_associations(adata, slide_name, distance_profiles, slide_dir):
    """Enhanced distance profiles with publication-quality styling"""
    
    treatment = adata.obs['treatment_status'].iloc[0]
    
    # Create figure with better proportions
    fig, (ax_main, ax_legend) = plt.subplots(1, 2, figsize=(16, 8), 
                                           gridspec_kw={'width_ratios': [4, 1]})
    
    # Enhanced color scheme for tumor types
    tumor_type_colors = {
        'TC_Glycolysis': '#E74C3C',    # Red - highlight as most important
        'TC_EMT': '#9B59B6',           # Purple 
        'TC_Proliferation': '#3498DB',  # Blue
        'TC_Quiescent': '#2ECC71',     # Green
        'TC_LipidMeta': '#F39C12'      # Orange
    }
    
    # Track lines for legend
    legend_elements = []
    
    for tc_type, profile in distance_profiles.items():
        distances = profile['distances']
        abundances = profile['abundances']
        
        if len(distances) > 2:
            # Get color and styling
            color = tumor_type_colors.get(tc_type, '#34495E')
            linewidth = 4 if 'Glycolysis' in tc_type else 3
            alpha = 1.0 if 'Glycolysis' in tc_type else 0.8
            
            # LOWESS smoothing with confidence intervals
            try:
                lowess_result = lowess(abundances, distances, frac=0.4, it=3, return_sorted=True)
                smooth_distances = lowess_result[:, 0]
                smooth_abundances = lowess_result[:, 1]
                
                # Plot main line with shadow effect for GLUT1
                if 'Glycolysis' in tc_type:
                    # Add shadow effect for emphasis
                    shadow_line = ax_main.plot(smooth_distances, smooth_abundances, 
                                             color='black', linewidth=linewidth+2, alpha=0.3, zorder=1)
                
                main_line = ax_main.plot(smooth_distances, smooth_abundances, 
                                       color=color, linewidth=linewidth, alpha=alpha, 
                                       label=tc_type.replace('TC_', ''), zorder=3)
                
                # Add confidence intervals (bootstrap-style)
                if len(distances) > 5:
                    # Simple confidence interval estimation
                    residuals = abundances - np.interp(distances, smooth_distances, smooth_abundances)
                    std_err = np.std(residuals)
                    ci_upper = smooth_abundances + 1.96 * std_err
                    ci_lower = smooth_abundances - 1.96 * std_err
                    
                    ax_main.fill_between(smooth_distances, ci_lower, ci_upper, 
                                       color=color, alpha=0.15, zorder=1)
                
                # Original data points with better styling
                ax_main.scatter(distances, abundances, color=color, s=40, 
                              alpha=0.7, edgecolors='white', linewidths=1, zorder=4)
                
                # Store for legend
                legend_elements.append(mpatches.Patch(color=color, label=tc_type.replace('TC_', '')))
                
            except Exception as e:
                # Fallback to simple plot
                ax_main.plot(distances, abundances, 'o-', color=color, 
                           label=tc_type.replace('TC_', ''), linewidth=linewidth, markersize=6)
    
    # Enhanced zone visualization
    zone_colors = {'Perivascular': '#E74C3C', 'Intermediate': '#F39C12', 'Distant': '#3498DB'}
    zone_alphas = {'Perivascular': 0.15, 'Intermediate': 0.10, 'Distant': 0.08}
    
    # Zone backgrounds
    ax_main.axvspan(0, 100, alpha=zone_alphas['Perivascular'], 
                   color=zone_colors['Perivascular'], zorder=0)
    ax_main.axvspan(100, 200, alpha=zone_alphas['Intermediate'], 
                   color=zone_colors['Intermediate'], zorder=0)
    ax_main.axvspan(200, 350, alpha=zone_alphas['Distant'], 
                   color=zone_colors['Distant'], zorder=0)
    
    # Zone boundary lines
    ax_main.axvline(x=100, color='gray', linestyle='--', alpha=0.7, linewidth=2)
    ax_main.axvline(x=200, color='gray', linestyle='--', alpha=0.7, linewidth=2)
    
    # Enhanced annotations
    ax_main.text(50, ax_main.get_ylim()[1]*0.95, 'Perivascular\n(0-100μm)', 
                ha='center', va='top', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    ax_main.text(150, ax_main.get_ylim()[1]*0.95, 'Intermediate\n(100-200μm)', 
                ha='center', va='top', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    ax_main.text(275, ax_main.get_ylim()[1]*0.95, 'Distant\n(>200μm)', 
                ha='center', va='top', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    
    # Enhanced main plot styling
    ax_main.set_xlabel('Distance to Nearest Vessel (μm)', fontsize=14, fontweight='bold')
    ax_main.set_ylabel('Tumor Cell Abundance', fontsize=14, fontweight='bold')
    ax_main.set_title(f'Vessel-Distance Profiles: {slide_name}\n({treatment.replace("_", " ").title()})', 
                     fontsize=16, fontweight='bold', pad=20)
    
    # Enhanced grid and spines
    ax_main.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.spines['left'].set_linewidth(1.5)
    ax_main.spines['bottom'].set_linewidth(1.5)
    
    # Enhanced legend in separate subplot
    ax_legend.axis('off')
    
    # Create custom legend
    legend_title = ax_legend.text(0.1, 0.9, 'Tumor Cell Types', fontsize=14, fontweight='bold',
                                 transform=ax_legend.transAxes)
    
    y_pos = 0.8
    for i, element in enumerate(legend_elements):
        # Color patch
        rect = Rectangle((0.1, y_pos - 0.02), 0.15, 0.04, 
                        facecolor=element.get_facecolor(), transform=ax_legend.transAxes)
        ax_legend.add_patch(rect)
        
        # Label with emphasis for GLUT1
        label_text = element.get_label()
        if 'Glycolysis' in label_text:
            label_text += ' ★'  # Star for emphasis
            fontweight = 'bold'
            fontsize = 12
        else:
            fontweight = 'normal'
            fontsize = 11
            
        ax_legend.text(0.3, y_pos, label_text, fontsize=fontsize, fontweight=fontweight,
                      va='center', transform=ax_legend.transAxes)
        y_pos -= 0.12
    
    # Add interpretation guide
    ax_legend.text(0.1, 0.4, 'Interpretation:', fontsize=12, fontweight='bold',
                  transform=ax_legend.transAxes)
    ax_legend.text(0.1, 0.32, '• Downward slope:\n  Vessel-seeking behavior', fontsize=10,
                  transform=ax_legend.transAxes, color=VESSEL_SEEKING_COLOR)
    ax_legend.text(0.1, 0.20, '• Upward slope:\n  Vessel-avoiding behavior', fontsize=10,
                  transform=ax_legend.transAxes, color=VESSEL_AVOIDING_COLOR)
    ax_legend.text(0.1, 0.08, '★ = Metastasis-associated', fontsize=10, fontweight='bold',
                  transform=ax_legend.transAxes, color=HIGH_RISK_COLOR)
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"enhanced_distance_profiles_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300, facecolor='white')
    plt.show()

#%% Step 5: Vascular Hotspot Co-occurrence Analysis (Metastasis Mechanism)

def analyze_vascular_hotspot_cooccurrence(adata, slide_name):
    """Analyze tumor-vessel co-occurrence in same spots (metastasis mechanism)"""
    
    print(f"  Analyzing vascular hotspot co-occurrence (metastasis potential)...")
    
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    
    # Focus on tumor regions
    tumor_mask = adata.obs['Is_Tumor_Region']
    hotspot_mask = adata.obs['Is_Vascular_Hotspot']
    
    # Co-occurrence: same spot has both vessels and tumor cells
    cooccurrence_results = {}
    
    print(f"    Co-occurrence analysis (same 100μm spot):")
    
    for tc_type in tumor_types:
        # Define tumor-positive spots (within tumor regions only)
        tumor_region_values = adata.obs.loc[tumor_mask, tc_type]
        tumor_positive_threshold = np.percentile(tumor_region_values, 90)  # Top 25% within tumor regions
        tumor_positive_mask = adata.obs[tc_type] >= tumor_positive_threshold
        
        # Co-occurrence: hotspot + tumor in same spot (within tumor regions)
        cooccurrence_mask = hotspot_mask & tumor_positive_mask & tumor_mask
        
        # Control groups (within tumor regions)
        hotspot_only = hotspot_mask & (~tumor_positive_mask) & tumor_mask
        tumor_only = tumor_positive_mask & (~hotspot_mask) & tumor_mask
        neither = (~hotspot_mask) & (~tumor_positive_mask) & tumor_mask
        
        n_cooccurrence = cooccurrence_mask.sum()
        n_hotspot_only = hotspot_only.sum()
        n_tumor_only = tumor_only.sum()
        n_neither = neither.sum()
        n_total_tumor = tumor_mask.sum()
        
        # Calculate metastasis potential metrics
        if n_total_tumor > 0:
            cooccurrence_rate = n_cooccurrence / n_total_tumor * 100
            
            # FIXED: Statistical enrichment within tumor regions only
            hotspot_rate_in_tumor = (hotspot_mask & tumor_mask).sum() / n_total_tumor
            tumor_rate_in_tumor = (tumor_positive_mask & tumor_mask).sum() / n_total_tumor
            expected_cooccurrence = hotspot_rate_in_tumor * tumor_rate_in_tumor * n_total_tumor
            
            enrichment_ratio = n_cooccurrence / (expected_cooccurrence + 1e-6)
            
            # ADDED: Fisher's exact test for significance
            contingency_table = [[n_cooccurrence, n_tumor_only], 
                               [n_hotspot_only, n_neither]]
            odds_ratio, p_value = fisher_exact(contingency_table)
            significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            
            cooccurrence_results[tc_type] = {
                'cooccurrence_spots': n_cooccurrence,
                'cooccurrence_rate': cooccurrence_rate,
                'enrichment_ratio': enrichment_ratio,
                'odds_ratio': odds_ratio,
                'p_value': p_value,
                'total_tumor_spots': n_total_tumor
            }
            
            print(f"      {tc_type}: {n_cooccurrence}/{n_total_tumor} spots ({cooccurrence_rate:.1f}%), "
                  f"enrichment={enrichment_ratio:.2f}x, OR={odds_ratio:.2f}, p={p_value:.3f} {significance}")
    
    return cooccurrence_results

def visualize_hotspot_cooccurrence(adata, slide_name, cooccurrence_results, slide_dir):
    """Enhanced metastasis potential with statistical annotations"""
    
    treatment = adata.obs['treatment_status'].iloc[0]
    
    if len(cooccurrence_results) > 0:
        # Create enhanced figure
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        tumor_types = list(cooccurrence_results.keys())
        tumor_types_clean = [tc.replace('TC_', '') for tc in tumor_types]
        
        # Data extraction
        cooccurrence_rates = [cooccurrence_results[tc]['cooccurrence_rate'] for tc in tumor_types]
        enrichment_ratios = [cooccurrence_results[tc]['enrichment_ratio'] for tc in tumor_types]
        p_values = [cooccurrence_results[tc]['p_value'] for tc in tumor_types]
        odds_ratios = [cooccurrence_results[tc]['odds_ratio'] for tc in tumor_types]
        
        # 1. Enhanced Co-occurrence Rates
        colors_rate = [HIGH_RISK_COLOR if 'Glycolysis' in tc else 
                      VESSEL_SEEKING_COLOR if rate > 2.0 else 
                      NEUTRAL_COLOR for tc, rate in zip(tumor_types, cooccurrence_rates)]
        
        bars1 = axes[0].bar(tumor_types_clean, cooccurrence_rates, color=colors_rate,
                           edgecolor='white', linewidth=1.5, alpha=0.8)
        
        # Add value labels on bars
        for i, (bar, rate, p_val) in enumerate(zip(bars1, cooccurrence_rates, p_values)):
            height = bar.get_height()
            significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
            axes[0].text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{rate:.1f}%\n{significance}', ha='center', va='bottom',
                        fontweight='bold', fontsize=10)
        
        axes[0].set_title('Vessel Co-occurrence Rate\n(Metastatic Potential)', 
                         fontsize=14, fontweight='bold', pad=20)
        axes[0].set_ylabel('Co-occurrence Rate (%)', fontsize=12, fontweight='bold')
        axes[0].set_xlabel('Tumor Cell Type', fontsize=12, fontweight='bold')
        axes[0].tick_params(axis='x', rotation=45, labelsize=11)
        axes[0].grid(True, alpha=0.3, axis='y')
        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)
        
        # 2. Enhanced Enrichment Ratios with significance zones
        colors_enrich = [HIGH_RISK_COLOR if ratio > 2.0 else 
                        VESSEL_SEEKING_COLOR if ratio > 1.5 else 
                        VESSEL_AVOIDING_COLOR if ratio < 0.67 else 
                        NEUTRAL_COLOR for ratio in enrichment_ratios]
        
        bars2 = axes[1].bar(tumor_types_clean, enrichment_ratios, color=colors_enrich,
                           edgecolor='white', linewidth=1.5, alpha=0.8)
        
        # Reference lines
        axes[1].axhline(y=1, color='black', linestyle='-', alpha=0.8, linewidth=2)
        axes[1].axhline(y=2.0, color=HIGH_RISK_COLOR, linestyle='--', alpha=0.7, 
                       label='High Risk (2x)')
        axes[1].axhline(y=1.5, color=VESSEL_SEEKING_COLOR, linestyle='--', alpha=0.7, 
                       label='Enriched (1.5x)')
        axes[1].axhline(y=0.67, color=VESSEL_AVOIDING_COLOR, linestyle='--', alpha=0.7, 
                       label='Depleted (0.67x)')
        
        # Add value labels
        for bar, ratio, p_val in zip(bars2, enrichment_ratios, p_values):
            height = bar.get_height()
            significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
            axes[1].text(bar.get_x() + bar.get_width()/2., height + 0.05,
                        f'{ratio:.2f}\n{significance}', ha='center', va='bottom',
                        fontweight='bold', fontsize=10)
        
        axes[1].set_title('Vessel Co-occurrence Enrichment\n(vs Expected)', 
                         fontsize=14, fontweight='bold', pad=20)
        axes[1].set_ylabel('Enrichment Ratio\n(Observed/Expected)', fontsize=12, fontweight='bold')
        axes[1].set_xlabel('Tumor Cell Type', fontsize=12, fontweight='bold')
        axes[1].tick_params(axis='x', rotation=45, labelsize=11)
        axes[1].legend(fontsize=10, loc='upper right')
        axes[1].grid(True, alpha=0.3, axis='y')
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['right'].set_visible(False)
        
        # 3. Enhanced Risk Assessment Heatmap
        # Create risk matrix
        risk_data = []
        for i, tc in enumerate(tumor_types_clean):
            risk_score = (cooccurrence_rates[i] / 5.0) * (enrichment_ratios[i])  # Normalized risk
            p_val = p_values[i]
            risk_data.append({
                'Tumor_Type': tc,
                'Co-occurrence_Rate': cooccurrence_rates[i],
                'Enrichment_Ratio': enrichment_ratios[i],
                'Risk_Score': risk_score,
                'P_Value': p_val,
                'Significant': p_val < 0.05
            })
        
        risk_df = pd.DataFrame(risk_data)
        
        # Create heatmap data
        heatmap_data = risk_df.set_index('Tumor_Type')[['Co-occurrence_Rate', 'Enrichment_Ratio', 'Risk_Score']]
        
        # Custom colormap
        cmap = sns.diverging_palette(240, 10, as_cmap=True)
        
        im = axes[2].imshow(heatmap_data.T, cmap='Reds', aspect='auto')
        
        # Set ticks and labels
        axes[2].set_xticks(range(len(tumor_types_clean)))
        axes[2].set_xticklabels(tumor_types_clean, rotation=45, ha='right')
        axes[2].set_yticks(range(len(heatmap_data.columns)))
        axes[2].set_yticklabels(['Co-occurrence\nRate (%)', 'Enrichment\nRatio', 'Risk\nScore'])
        
        # Add text annotations
        for i in range(len(tumor_types_clean)):
            for j, col in enumerate(heatmap_data.columns):
                value = heatmap_data.iloc[i, j]
                text = f'{value:.1f}' if j == 0 else f'{value:.2f}'
                
                # Add significance markers
                if risk_df.iloc[i]['Significant']:
                    text += '*'
                
                axes[2].text(i, j, text, ha='center', va='center', 
                           fontweight='bold', fontsize=10,
                           color='white' if value > heatmap_data.iloc[:, j].mean() else 'black')
        
        axes[2].set_title('Metastatic Risk Assessment\n(* = p<0.05)', 
                         fontsize=14, fontweight='bold', pad=20)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=axes[2], shrink=0.8)
        cbar.set_label('Relative Risk Level', rotation=270, labelpad=20, fontweight='bold')
        
        plt.suptitle(f'Metastasis Potential Analysis: {slide_name}\n({treatment.replace("_", " ").title()})', 
                     fontsize=16, fontweight='bold', y=1.02)
        
        plt.tight_layout()
        plt.savefig(os.path.join(slide_dir, f"enhanced_metastasis_potential_{slide_name}.pdf"), 
                   bbox_inches='tight', dpi=300, facecolor='white')
        plt.show()

#%% Step 6: Treatment Effect Analysis

def analyze_treatment_effects_single_slide(adata, slide_name, associations, cooccurrence_results):
    """Extract treatment-relevant metrics for cross-slide comparison"""
    
    treatment = adata.obs['treatment_status'].iloc[0]
    
    # Treatment metrics
    treatment_metrics = {
        'slide_name': slide_name,
        'treatment': treatment,
        'n_spots': len(adata.obs),
        'tumor_spots': adata.obs['Is_Tumor_Region'].sum(),
        'vascular_hotspots': adata.obs['Is_Vascular_Hotspot'].sum()
    }
    
    # Tumor cell abundances
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    tumor_mask = adata.obs['Is_Tumor_Region']
    
    for tc_type in tumor_types:
        if tumor_mask.sum() > 0:
            treatment_metrics[f'{tc_type}_mean_abundance'] = adata.obs.loc[tumor_mask, tc_type].mean()
    
    # Vessel associations
    for tc_type, assoc in associations.items():
        treatment_metrics[f'{tc_type}_vessel_correlation'] = assoc['distance_correlation']
    
    # Co-occurrence metrics
    for tc_type, cooc in cooccurrence_results.items():
        treatment_metrics[f'{tc_type}_cooccurrence_rate'] = cooc['cooccurrence_rate']
        treatment_metrics[f'{tc_type}_enrichment_ratio'] = cooc['enrichment_ratio']
    
    return treatment_metrics

#%% Step 7: Spatial Competition and EMT Transition Analysis

def analyze_spatial_competition_emt(adata, slide_name):
    """Analyze tumor subtype spatial competition and EMT transitions"""
    
    print(f"  Analyzing spatial competition and EMT transitions...")
    
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    tumor_mask = adata.obs['Is_Tumor_Region']
    tumor_data = adata.obs[tumor_mask]
    
    # Pairwise correlations
    tumor_pairs = [
        ('TC_Glycolysis', 'TC_Quiescent'),
        ('TC_Glycolysis', 'TC_EMT'),
        ('TC_Quiescent', 'TC_EMT'),
        ('TC_EMT', 'TC_Proliferation')
    ]
    
    pairwise_results = {}
    
    print(f"    Pairwise correlations in tumor regions:")
    
    for tc1, tc2 in tumor_pairs:
        if tc1 in tumor_types and tc2 in tumor_types:
            if tumor_data[tc1].var() > 0 and tumor_data[tc2].var() > 0:
                corr, p_value = spearmanr(tumor_data[tc1], tumor_data[tc2])
                
                relationship = "co-occurrence" if corr > 0.2 else \
                             "mutual exclusion" if corr < -0.2 else "independent"
                
                pairwise_results[f'{tc1}_vs_{tc2}'] = {
                    'correlation': corr,
                    'p_value': p_value,
                    'relationship': relationship
                }
                
                significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                print(f"      {tc1} vs {tc2}: r={corr:+.3f} ({relationship}) p={p_value:.3f} {significance}")
    
    # EMT region analysis using gene-based score
    emt_results = {}
    
    if 'EMT_Score' in adata.obs.columns:
        # Define high EMT regions using gene-based score
        emt_threshold = np.percentile(tumor_data['EMT_Score'], 75)  # Top 25%
        high_emt_mask = tumor_data['EMT_Score'] >= emt_threshold
        
        n_high_emt = high_emt_mask.sum()
        print(f"    High EMT regions: {n_high_emt} spots ({n_high_emt/len(tumor_data)*100:.1f}% of tumor)")
        
        if n_high_emt > 5:
            # Analyze other tumor types in EMT regions
            for tc_type in tumor_types:
                if tc_type != 'TC_EMT':  # Skip TC_EMT itself
                    high_emt_abundance = tumor_data.loc[high_emt_mask, tc_type].mean()
                    low_emt_abundance = tumor_data.loc[~high_emt_mask, tc_type].mean()
                    
                    fold_change = high_emt_abundance / (low_emt_abundance + 1e-6)
                    
                    emt_results[f'{tc_type}_in_high_emt'] = {
                        'fold_enrichment': fold_change,
                        'high_emt_abundance': high_emt_abundance,
                        'low_emt_abundance': low_emt_abundance
                    }
                    
                    print(f"      {tc_type} in high EMT: {fold_change:.2f}x enrichment")
    
    return pairwise_results, emt_results

def visualize_competition_emt(adata, slide_name, pairwise_results, emt_results, slide_dir):
    """Visualize spatial competition and EMT patterns"""
    
    treatment = adata.obs['treatment_status'].iloc[0]
    
    # 1. Pairwise correlation matrix
    if len(pairwise_results) > 0:
        tumor_pairs = list(pairwise_results.keys())
        all_tumor_types = list(set([pair.split('_vs_')[0] for pair in tumor_pairs] + 
                                  [pair.split('_vs_')[1] for pair in tumor_pairs]))
        
        # Create correlation matrix
        n_types = len(all_tumor_types)
        corr_matrix = np.eye(n_types)
        
        for pair, results in pairwise_results.items():
            tc1, tc2 = pair.split('_vs_')
            if tc1 in all_tumor_types and tc2 in all_tumor_types:
                i1 = all_tumor_types.index(tc1)
                i2 = all_tumor_types.index(tc2)
                corr_matrix[i1, i2] = results['correlation']
                corr_matrix[i2, i1] = results['correlation']
        
        plt.figure(figsize=(8, 6))
        sns.heatmap(corr_matrix, 
                   xticklabels=all_tumor_types,
                   yticklabels=all_tumor_types,
                   annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                   cbar_kws={'label': 'Correlation'})
        
        plt.title(f'Tumor Subtype Spatial Competition - {slide_name} ({treatment})', fontweight='bold')
        plt.xticks(rotation=45)
        plt.yticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(slide_dir, f"spatial_competition_{slide_name}.pdf"), 
                   bbox_inches='tight', dpi=300)
        plt.show()
    
    # 2. EMT region analysis
    if 'EMT_Score' in adata.obs.columns:
        tumor_mask = adata.obs['Is_Tumor_Region']
        adata_tumor = adata[tumor_mask].copy()
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # EMT score spatial distribution
        sc.pl.spatial(adata_tumor, color='EMT_Score', cmap='viridis',
                     ax=axes[0], size=1.2, img_key='hires', show=False,
                     title='EMT Score (Gene-based)')
        
        # EMT enrichment in other tumor types
        if len(emt_results) > 0:
            tumor_types_emt = []
            enrichments = []
            
            for key, results in emt_results.items():
                if key.endswith('_in_high_emt'):
                    tumor_type = key.replace('_in_high_emt', '')
                    tumor_types_emt.append(tumor_type)
                    enrichments.append(results['fold_enrichment'])
            
            if len(enrichments) > 0:
                colors = ['red' if e > 1.2 else 'blue' if e < 0.8 else 'gray' for e in enrichments]
                bars = axes[1].bar(tumor_types_emt, enrichments, color=colors)
                axes[1].axhline(y=1, color='black', linestyle='-', alpha=0.5)
                axes[1].set_title('Enrichment in High EMT Regions')
                axes[1].set_ylabel('Fold Enrichment')
                axes[1].tick_params(axis='x', rotation=45)
        
        plt.suptitle(f'EMT Transition Analysis - {slide_name} ({treatment})', 
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(slide_dir, f"emt_transitions_{slide_name}.pdf"), 
                   bbox_inches='tight', dpi=300)
        plt.show()

#%% STEP 1: IDENTIFY VESSEL-ASSOCIATED GENES

def identify_vessel_proximity_genes(adata, slide_name):
    """
    Step 1: Identify genes associated with vessel proximity
    Returns genes correlated with distance to vessels
    """
    print(f"\n--- Step 1: Identifying vessel-proximity genes ({slide_name}) ---")
    
    if 'Distance_to_Vessel_um' not in adata.obs.columns:
        print("  No vessel distance data found - run vessel analysis first")
        return None, None
    
    # Focus on tumor regions only
    tumor_mask = adata.obs.get('Is_Tumor_Region', True)
    
    if tumor_mask.sum() < PATHWAY_ANALYSIS_PARAMS['min_spots_for_analysis']:
        print(f"  Insufficient tumor spots for analysis: {tumor_mask.sum()}")
        return None, None
    
    # Get gene expression data
    tumor_data = adata[tumor_mask].copy()
    distances = tumor_data.obs['Distance_to_Vessel_um'].values
    
    # Calculate gene-distance correlations
    gene_correlations = {}
    vessel_associated_genes = {'vessel_seeking': [], 'vessel_avoiding': []}
    
    print(f"  Analyzing {tumor_data.n_vars} genes across {tumor_data.n_obs} tumor spots...")
    
    if hasattr(tumor_data.X, 'toarray'):
        gene_expr_all = tumor_data.X.toarray()
    else:
        gene_expr_all = tumor_data.X

    for gene_idx, gene_name in enumerate(tumor_data.var.index):
        if gene_idx % 1000 == 0:
            print(f"    Processed {gene_idx}/{tumor_data.n_vars} genes")
        
        # Get gene expression
        gene_expr = gene_expr_all[:, gene_idx].flatten()
        
        # Skip genes with low variance
        if np.var(gene_expr) < 1e-6:
            continue
        
        # Calculate correlation with distance
        if len(gene_expr) == len(distances):
            corr, p_value = pearsonr(distances, gene_expr)
            
            if abs(corr) >= PATHWAY_ANALYSIS_PARAMS['distance_correlation_threshold'] and p_value < 0.05:
                gene_correlations[gene_name] = {
                    'correlation': corr,
                    'p_value': p_value,
                    'abs_correlation': abs(corr)
                }
                
                # Classify genes
                if corr < -PATHWAY_ANALYSIS_PARAMS['distance_correlation_threshold']:
                    vessel_associated_genes['vessel_seeking'].append(gene_name)  # Higher near vessels
                elif corr > PATHWAY_ANALYSIS_PARAMS['distance_correlation_threshold']:
                    vessel_associated_genes['vessel_avoiding'].append(gene_name)  # Lower near vessels
    
    print(f"  Found {len(gene_correlations)} vessel-associated genes:")
    print(f"    Vessel-seeking (high near vessels): {len(vessel_associated_genes['vessel_seeking'])}")
    print(f"    Vessel-avoiding (low near vessels): {len(vessel_associated_genes['vessel_avoiding'])}")
    
    return gene_correlations, vessel_associated_genes

def visualize_vessel_gene_heatmap(adata, gene_correlations, vessel_associated_genes, slide_dir, slide_name):
    """
    Visualize vessel-associated genes with heatmap
    """
    if not gene_correlations:
        return
    
    print(f"  Creating vessel-gene association heatmap...")
    
    # Get top genes for visualization
    sorted_genes = sorted(gene_correlations.items(), key=lambda x: x[1]['abs_correlation'], reverse=True)
    top_genes = [gene for gene, data in sorted_genes[:50]]  # Top 50 genes
    
    if len(top_genes) == 0:
        return
    
    # Focus on tumor regions
    tumor_mask = adata.obs.get('Is_Tumor_Region', True)
    tumor_data = adata[tumor_mask].copy()
    
    # Create spatial zones based on vessel distance
    distances = tumor_data.obs['Distance_to_Vessel_um'].values
    
    # Define zones
    perivascular_mask = distances <= PATHWAY_ANALYSIS_PARAMS['perivascular_zone']
    intermediate_mask = (distances > PATHWAY_ANALYSIS_PARAMS['perivascular_zone']) & \
                       (distances <= PATHWAY_ANALYSIS_PARAMS['intermediate_zone'])
    distant_mask = distances > PATHWAY_ANALYSIS_PARAMS['intermediate_zone']
    
    zones = ['Perivascular', 'Intermediate', 'Distant']
    zone_masks = [perivascular_mask, intermediate_mask, distant_mask]
    
    # Calculate mean expression in each zone
    expression_matrix = []
    gene_names_clean = []
    
    for gene in top_genes:
        if gene in tumor_data.var.index:
            gene_idx = tumor_data.var.index.get_loc(gene)
            
            if hasattr(tumor_data.X, 'toarray'):
                gene_expr = tumor_data.X[:, gene_idx].toarray().flatten()
            else:
                gene_expr = tumor_data.X[:, gene_idx]
            
            zone_expressions = []
            for mask in zone_masks:
                if mask.sum() > 0:
                    zone_expr = gene_expr[mask].mean()
                else:
                    zone_expr = 0
                zone_expressions.append(zone_expr)
            
            expression_matrix.append(zone_expressions)
            gene_names_clean.append(gene.replace('_', '-'))  # Clean gene names
    
    if len(expression_matrix) == 0:
        return
    
    # Create heatmap
    fig, axes = plt.subplots(1, 2, figsize=(16, 12))
    
    # 1. Expression heatmap
    expression_df = pd.DataFrame(expression_matrix, columns=zones, index=gene_names_clean)
    
    # Normalize by row (gene) for better visualization
    expression_df_norm = expression_df.div(expression_df.max(axis=1), axis=0)
    
    sns.heatmap(expression_df_norm, ax=axes[0], cmap='magma', 
                cbar_kws={'label': 'Normalized Expression'})
    axes[0].set_title(f'Vessel-Associated Gene Expression\n{slide_name}\n(Top {len(top_genes)} genes)', 
                     fontweight='bold')
    axes[0].set_xlabel('Spatial Zone')
    axes[0].set_ylabel('Genes')
    
    # 2. Correlation with distance
    correlations = [gene_correlations[gene]['correlation'] for gene in top_genes if gene in gene_correlations]
    p_values = [gene_correlations[gene]['p_value'] for gene in top_genes if gene in gene_correlations]
    
    if len(correlations) > 0:
        # Color code by significance
        colors = ['red' if p < 0.001 else 'orange' if p < 0.01 else 'yellow' if p < 0.05 else 'gray' 
                 for p in p_values]
        
        bars = axes[1].barh(range(len(correlations)), correlations, color=colors, alpha=0.7)
        axes[1].set_yticks(range(len(gene_names_clean)))
        axes[1].set_yticklabels(gene_names_clean)
        axes[1].set_xlabel('Correlation with Vessel Distance')
        axes[1].set_title('Gene-Vessel Distance Correlations\n(Red: p<0.001, Orange: p<0.01, Yellow: p<0.05)')
        axes[1].axvline(x=0, color='black', linestyle='-', alpha=0.5)
        
        # Add correlation value labels
        for i, (bar, corr) in enumerate(zip(bars, correlations)):
            axes[1].text(corr + (0.01 if corr > 0 else -0.01), bar.get_y() + bar.get_height()/2, 
                        f'{corr:.2f}', ha='left' if corr > 0 else 'right', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"vessel_gene_heatmap_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()

#%% STEP 2: PATHWAY ENRICHMENT ANALYSIS

def load_pathway_databases():
    """
    Load pathway databases for enrichment analysis
    """
    print(f"\n--- Loading pathway databases ---")
    msig = Msigdb()

    try:
        print(f"  Fetching Hallmark gene sets...")
        category = 'h.all'  # Hallmark gene sets
        dbver = "2025.1.Hs"  # Specify the version
        print(f"Fetching Hallmark gene sets from {category} with version {dbver}...")

        hallmark = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
        
        # Define vessel-relevant pathway keywords
        vessel_keywords = ['angiogenesis', 'hypoxia', 'glycolysis', 'oxidative', 'metabolism', 
                          'emt', 'epithelial', 'migration', 'invasion', 'inflammation', 'immune']
        
        # Filter relevant pathways
        vessel_pathways = {}
        for pathway_name, genes in hallmark.items():
            pathway_lower = pathway_name.lower()
            if any(keyword in pathway_lower for keyword in vessel_keywords):
                clean_name = pathway_name.replace('HALLMARK_', '')
                vessel_pathways[clean_name] = genes
        
        print(f"  ✅ Loaded {len(vessel_pathways)} vessel-relevant Hallmark pathways")
        
        # Add custom vessel-related pathways
        # custom_pathways = {
        #     'VESSEL_COMMUNICATION': ['VEGFA', 'VEGFB', 'VEGFC', 'FLT1', 'KDR', 'FLT4', 'ANGPT1', 'ANGPT2', 'TEK'],
        #     'METABOLIC_ADAPTATION': ['GLUT1', 'GLUT3', 'HK1', 'HK2', 'PFKL', 'PFKM', 'PKM', 'LDHA', 'LDHB'],
        #     'INVASION_SIGNALS': ['MMP2', 'MMP9', 'MMP14', 'TIMP1', 'TIMP2', 'ITGA1', 'ITGA2', 'ITGB1', 'ITGB3'],
        #     'CHEMOTAXIS': ['CXCL12', 'CXCR4', 'CCL2', 'CCR2', 'CXCL8', 'CXCR1', 'CXCR2']
        # }
        
        # vessel_pathways.update(custom_pathways)
        # print(f"  ✅ Added {len(custom_pathways)} custom vessel pathways")
        
        return vessel_pathways
        
    except Exception as e:
        print(f"  ❌ Error loading pathways: {e}")
        print(f"  Using backup pathway definitions...")
        
        # Backup pathway definitions
        backup_pathways = {
            'GLYCOLYSIS': ['GLUT1', 'HK1', 'HK2', 'PFKL', 'ALDOA', 'TPI1', 'GAPDH', 'PGK1', 'PGAM1', 'ENO1', 'PKM', 'LDHA'],
            'HYPOXIA_RESPONSE': ['HIF1A', 'VEGFA', 'CA9', 'GLUT1', 'LDHA', 'PDK1', 'PGK1', 'ALDOA', 'ENO1'],
            'ANGIOGENESIS': ['VEGFA', 'VEGFB', 'VEGFC', 'ANGPT1', 'ANGPT2', 'FLT1', 'KDR', 'TEK', 'PDGFA', 'PDGFB'],
            'EMT': ['VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'ZEB2', 'TWIST1', 'CDH1', 'CDH2', 'FN1', 'FOXC2'],
            'CELL_MIGRATION': ['ROCK1', 'ROCK2', 'RHOA', 'CDC42', 'RAC1', 'ACTB', 'MYH9', 'VCL', 'PTK2', 'ITGB1']
        }
        
        return backup_pathways

def calculate_pathway_scores(adata, vessel_pathways, method='mean'):
    """
    Calculate pathway enrichment scores for each spot
    """
    print(f"  Calculating pathway scores using {method} method...")
    
    pathway_scores = {}
    available_genes = adata.var.index.tolist()

    if hasattr(adata.X, 'toarray'):
        all_expr = adata.X.toarray()
    else:
        all_expr = adata.X
    
    for pathway_name, pathway_genes in vessel_pathways.items():
        # Find available genes in the pathway
        available_pathway_genes = [gene for gene in pathway_genes if gene in available_genes]
        
        if len(available_pathway_genes) >= 3:  # Minimum genes for reliable score
            if method == 'mean':
                # Simple mean expression method
                pathway_indices = [adata.var.index.get_loc(gene) for gene in available_pathway_genes]
                pathway_expr = all_expr[:, pathway_indices]
                
                pathway_score = np.mean(pathway_expr, axis=1)
                pathway_scores[pathway_name] = pathway_score
                
                print(f"    {pathway_name}: {len(available_pathway_genes)}/{len(pathway_genes)} genes available")
    
    return pathway_scores

def perform_vessel_pathway_enrichment(adata, gene_correlations, vessel_pathways, slide_name):
    """
    Step 2: Perform pathway enrichment analysis on vessel-associated genes
    """
    print(f"\n--- Step 2: Pathway enrichment analysis ({slide_name}) ---")
    
    if not gene_correlations:
        return None
    
    # Get vessel-seeking and vessel-avoiding genes
    vessel_seeking_genes = [gene for gene, data in gene_correlations.items() 
                           if data['correlation'] < -PATHWAY_ANALYSIS_PARAMS['distance_correlation_threshold']]
    vessel_avoiding_genes = [gene for gene, data in gene_correlations.items() 
                            if data['correlation'] > PATHWAY_ANALYSIS_PARAMS['distance_correlation_threshold']]
    
    print(f"  Vessel-seeking genes: {len(vessel_seeking_genes)}")
    print(f"  Vessel-avoiding genes: {len(vessel_avoiding_genes)}")
    
    # Calculate pathway scores for all spots
    pathway_scores = calculate_pathway_scores(adata, vessel_pathways)
    
    # Store pathway scores in adata
    for pathway_name, scores in pathway_scores.items():
        adata.obs[f'pathway_{pathway_name}'] = scores
    
    # Perform enrichment analysis
    enrichment_results = {}
    
    for gene_set_name, genes in [('vessel_seeking', vessel_seeking_genes), 
                                ('vessel_avoiding', vessel_avoiding_genes)]:
        
        if len(genes) < 3:
            continue
        
        set_enrichments = {}
        
        for pathway_name, pathway_genes in vessel_pathways.items():
            # Calculate overlap
            overlap_genes = set(genes) & set(pathway_genes)
            overlap_count = len(overlap_genes)
            
            if overlap_count >= 2:  # Minimum overlap
                # Simple enrichment calculation (can be enhanced with Fisher's exact test)
                total_genes_in_pathway = len(pathway_genes)
                total_vessel_genes = len(genes)
                enrichment_ratio = (overlap_count / total_vessel_genes) / (total_genes_in_pathway / 20000)  # Assume ~20k total genes
                
                set_enrichments[pathway_name] = {
                    'overlap_count': overlap_count,
                    'overlap_genes': list(overlap_genes),
                    'enrichment_ratio': enrichment_ratio,
                    'pathway_size': total_genes_in_pathway
                }
        
        enrichment_results[gene_set_name] = set_enrichments
        
        print(f"  {gene_set_name} enrichment: {len(set_enrichments)} pathways enriched")
    
    return enrichment_results, pathway_scores

def visualize_pathway_enrichment(enrichment_results, pathway_scores, adata, slide_dir, slide_name):
    """
    Visualize pathway enrichment results in enhanced 1x3 layout
    """
    if not enrichment_results:
        return
    
    print(f"  Creating pathway enrichment visualizations...")
    
    # Create 1x3 subplot layout
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    
    # 1. Vessel-seeking gene enrichment
    create_enrichment_barplot(enrichment_results, 'vessel_seeking', axes[0], 
                             'Vessel-Seeking Gene Enrichment', 'Reds')
    
    # 2. Vessel-avoiding gene enrichment  
    create_enrichment_barplot(enrichment_results, 'vessel_avoiding', axes[1],
                             'Vessel-Avoiding Gene Enrichment', 'Blues')
    
    # 3. Enhanced pathway score comparison with statistics
    create_pathway_comparison_plot(pathway_scores, adata, axes[2])
    
    # Overall styling
    fig.suptitle(f'Pathway Enrichment Analysis - {slide_name}', 
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"pathway_enrichment_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()


def create_enrichment_barplot(enrichment_results, gene_set_type, ax, title, colormap):
    """
    Create enhanced enrichment bar plot
    """
    if gene_set_type not in enrichment_results:
        ax.text(0.5, 0.5, f'No {gene_set_type} data available', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontweight='bold', fontsize=12)
        return
    
    enrichment_data = enrichment_results[gene_set_type]
    
    if len(enrichment_data) == 0:
        ax.text(0.5, 0.5, f'No enriched pathways found', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontweight='bold', fontsize=12)
        return
    
    # Prepare data
    pathways = list(enrichment_data.keys())
    enrichments = [enrichment_data[p]['enrichment_ratio'] for p in pathways]
    overlaps = [enrichment_data[p]['overlap_count'] for p in pathways]
    
    # Sort by enrichment ratio
    sorted_data = sorted(zip(pathways, enrichments, overlaps), 
                        key=lambda x: x[1], reverse=True)
    pathways_sorted, enrichments_sorted, overlaps_sorted = zip(*sorted_data)
    
    # Take top 10 for better visualization
    n_show = min(10, len(pathways_sorted))
    pathways_show = pathways_sorted[:n_show]
    enrichments_show = enrichments_sorted[:n_show]
    overlaps_show = overlaps_sorted[:n_show]
    
    # Enhanced coloring
    colors = plt.cm.get_cmap(colormap)([o/max(overlaps_show) for o in overlaps_show])
    
    # Create bars
    bars = ax.barh(range(len(pathways_show)), enrichments_show, 
                   color=colors, alpha=0.8, edgecolor='white', linewidth=0.5)
    
    # Formatting
    ax.set_yticks(range(len(pathways_show)))
    ax.set_yticklabels([p.replace('_', ' ').title() for p in pathways_show], fontsize=10)
    ax.set_xlabel('Enrichment Ratio', fontweight='bold', fontsize=11)
    ax.set_title(title, fontweight='bold', fontsize=12, pad=15)
    ax.axvline(x=1, color='black', linestyle='--', alpha=0.6, linewidth=1)
    
    # Add value labels on bars
    for i, (bar, enrichment) in enumerate(zip(bars, enrichments_show)):
        ax.text(enrichment + 0.05, bar.get_y() + bar.get_height()/2, 
               f'{enrichment:.2f}', ha='left', va='center', fontsize=9, 
               fontweight='bold')
    
    # Style improvements
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.3, linestyle='-', linewidth=0.5)


def create_pathway_comparison_plot(pathway_scores, adata, ax):
    """
    Create enhanced pathway comparison plot with statistical tests
    """
    if not pathway_scores:
        ax.text(0.5, 0.5, 'No pathway scores available', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Pathway Score Comparison', fontweight='bold', fontsize=12)
        return
    
    # Get top pathways for visualization
    top_pathways = list(pathway_scores.keys())  # Show top 6 pathways
    
    # Prepare data for plotting
    plot_data = []
    pathway_region_data = {}
    
    for pathway in top_pathways:
        pathway_col = f'pathway_{pathway}'
        if pathway_col in adata.obs.columns:
            tumor_mask = adata.obs.get('Is_Tumor_Region', True)
            tumor_data = adata.obs[tumor_mask]
            
            distances = tumor_data['Distance_to_Vessel_um'].values
            scores = tumor_data[pathway_col].values
            
            # Define spatial zones
            perivascular_scores = []
            intermediate_scores = []
            distant_scores = []
            
            for distance, score in zip(distances, scores):
                if distance <= PATHWAY_ANALYSIS_PARAMS['perivascular_zone']:
                    perivascular_scores.append(score)
                    plot_data.append({
                        'pathway': pathway.replace('_', ' ').title(),
                        'region': 'Perivascular',
                        'score': score
                    })
                elif distance <= PATHWAY_ANALYSIS_PARAMS['intermediate_zone']:
                    intermediate_scores.append(score)
                    plot_data.append({
                        'pathway': pathway.replace('_', ' ').title(),
                        'region': 'Intermediate', 
                        'score': score
                    })
                else:
                    distant_scores.append(score)
                    plot_data.append({
                        'pathway': pathway.replace('_', ' ').title(),
                        'region': 'Distant',
                        'score': score
                    })
            
            # Store for statistical testing
            pathway_region_data[pathway] = {
                'Perivascular': perivascular_scores,
                'Intermediate': intermediate_scores,
                'Distant': distant_scores
            }
    
    if not plot_data:
        ax.text(0.5, 0.5, 'No pathway data found', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Pathway Score Comparison', fontweight='bold', fontsize=12)
        return
    
    # Create DataFrame
    plot_df = pd.DataFrame(plot_data)
    
    # Create enhanced box plot
    box_plot = sns.boxplot(data=plot_df, x='pathway', y='score', hue='region', 
                          ax=ax, palette=['#FF6B6B', '#4ECDC4', '#45B7D1'])
    
    # Add statistical annotations
    add_statistical_annotations(ax, pathway_region_data, plot_df)
    
    # Formatting
    ax.set_xlabel('Pathways', fontweight='bold', fontsize=11)
    ax.set_ylabel('Pathway Score', fontweight='bold', fontsize=11)
    ax.set_title('Pathway Scores by Spatial Region', fontweight='bold', fontsize=12, pad=15)
    
    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    
    # Enhanced legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='Spatial Region',
             frameon=True, fancybox=True, shadow=True, loc='upper right')
    
    # Style improvements
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)


def add_statistical_annotations(ax, pathway_region_data, plot_df):
    """
    Add statistical significance annotations using Wilcoxon test
    """
    pathways = list(pathway_region_data.keys())
    x_positions = range(len(pathways))
    
    for i, pathway in enumerate(pathways):
        region_data = pathway_region_data[pathway]
        
        # Compare Perivascular vs Distant (most biologically relevant)
        if (len(region_data['Perivascular']) > 0 and 
            len(region_data['Distant']) > 0):
            
            try:
                # Perform Mann-Whitney U test (non-parametric alternative to t-test)
                statistic, p_value = mannwhitneyu(
                    region_data['Perivascular'], 
                    region_data['Distant'],
                    alternative='two-sided'
                )
                
                # Add significance annotation
                y_max = plot_df[plot_df['pathway'] == pathway.replace('_', ' ').title()]['score'].max()
                y_annotation = y_max + (y_max * 0.1)
                
                # Significance levels
                if p_value < 0.001:
                    sig_text = '***'
                elif p_value < 0.01:
                    sig_text = '**'
                elif p_value < 0.05:
                    sig_text = '*'
                else:
                    sig_text = 'ns'
                
                # Add significance text
                ax.text(i, y_annotation, sig_text, ha='center', va='bottom',
                       fontweight='bold', fontsize=10)
                
            except Exception as e:
                continue  # Skip if statistical test fails


def get_pathway_analysis_params():
    """
    Helper function to get pathway analysis parameters
    """
    return {
        'perivascular_zone': 50,  # μm
        'intermediate_zone': 150,  # μm
        # Add other parameters as needed
    }

#%% STEP 3: CELL TYPE-PATHWAY CORRELATIONS

def analyze_cell_pathway_correlations(adata, pathway_scores, slide_name):
    """
    Step 3: Analyze correlation between cell types and pathway scores
    """
    print(f"\n--- Step 3: Cell type-pathway correlations ({slide_name}) ---")
    
    if not pathway_scores:
        return None
    
    # Focus on tumor regions and perivascular zones
    tumor_mask = adata.obs.get('Is_Tumor_Region', True)
    perivascular_mask = adata.obs['Distance_to_Vessel_um'] <= PATHWAY_ANALYSIS_PARAMS['perivascular_zone']
    
    analysis_mask = tumor_mask & perivascular_mask
    
    if analysis_mask.sum() < PATHWAY_ANALYSIS_PARAMS['min_spots_for_analysis']:
        print(f"  Insufficient perivascular tumor spots: {analysis_mask.sum()}")
        return None
    
    analysis_data = adata.obs[analysis_mask]
    
    # Get malignant cell types
    malignant_cell_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    
    # Also include other relevant cell types
    other_cell_types = [col for col in adata.obs.columns 
                       if any(cell_type in col for cell_type in ['CD8', 'CD4', 'DC_', 'Mac_', 'CAF_'])]
    
    all_cell_types = malignant_cell_types + other_cell_types
    
    # Calculate correlations
    correlation_results = {}
    
    print(f"  Analyzing correlations for {len(all_cell_types)} cell types and {len(pathway_scores)} pathways...")
    
    for cell_type in all_cell_types:
        if cell_type in analysis_data.columns:
            cell_abundances = analysis_data[cell_type].values
            
            if np.var(cell_abundances) < 1e-6:  # Skip cells with no variance
                continue
            
            cell_correlations = {}
            
            for pathway_name in pathway_scores.keys():
                pathway_col = f'pathway_{pathway_name}'
                
                if pathway_col in analysis_data.columns:
                    pathway_values = analysis_data[pathway_col].values
                    
                    if np.var(pathway_values) > 1e-6:
                        corr, p_value = pearsonr(cell_abundances, pathway_values)
                        
                        if abs(corr) > 0.2 and p_value < 0.05:  # Significant correlations
                            cell_correlations[pathway_name] = {
                                'correlation': corr,
                                'p_value': p_value,
                                'abs_correlation': abs(corr)
                            }
            
            if cell_correlations:
                correlation_results[cell_type] = cell_correlations
                
                # Print top correlations for malignant cells
                if cell_type.startswith('TC_'):
                    top_pathways = sorted(cell_correlations.items(), 
                                        key=lambda x: x[1]['abs_correlation'], reverse=True)[:3]
                    
                    print(f"  {cell_type} top correlations:")
                    for pathway, data in top_pathways:
                        corr = data['correlation']
                        p_val = data['p_value']
                        print(f"    {pathway}: r={corr:+.3f}, p={p_val:.3f}")
    
    return correlation_results

def visualize_cell_pathway_correlations(adata, correlation_results, slide_dir, slide_name):
    """
    Visualize cell type-pathway correlations
    """
    if not correlation_results:
        return
    
    print(f"  Creating cell-pathway correlation visualizations...")
    
    # Create correlation matrix
    all_pathways = set()
    all_cell_types = list(correlation_results.keys())
    
    for cell_correlations in correlation_results.values():
        all_pathways.update(cell_correlations.keys())
    
    all_pathways = sorted(list(all_pathways))
    
    # Build correlation matrix
    correlation_matrix = []
    
    for cell_type in all_cell_types:
        cell_row = []
        for pathway in all_pathways:
            if pathway in correlation_results[cell_type]:
                cell_row.append(correlation_results[cell_type][pathway]['correlation'])
            else:
                cell_row.append(0)
        correlation_matrix.append(cell_row)
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    
    # 1. Correlation heatmap
    correlation_df = pd.DataFrame(correlation_matrix, 
                                 index=[ct.replace('TC_', '').replace('_', ' ') for ct in all_cell_types],
                                 columns=[pw.replace('_', ' ') for pw in all_pathways])
    
    sns.heatmap(correlation_df, ax=axes[0], cmap='RdBu_r', center=0, 
               annot=True, fmt='.2f', cbar_kws={'label': 'Correlation Coefficient'})
    axes[0].set_title(f'Cell Type-Pathway Correlations\n{slide_name}\n(Perivascular Region)')
    axes[0].set_xlabel('Pathways')
    axes[0].set_ylabel('Cell Types')
    
    # 2. Focus on malignant cells
    malignant_correlations = {ct: corrs for ct, corrs in correlation_results.items() 
                            if ct.startswith('TC_')}
    
    if malignant_correlations:
        malignant_data = []
        
        for cell_type, pathway_corrs in malignant_correlations.items():
            for pathway, corr_data in pathway_corrs.items():
                malignant_data.append({
                    'cell_type': cell_type.replace('TC_', ''),
                    'pathway': pathway.replace('_', ' '),
                    'correlation': corr_data['correlation'],
                    'abs_correlation': corr_data['abs_correlation'],
                    'p_value': corr_data['p_value']
                })
        
        malignant_df = pd.DataFrame(malignant_data)
        
        # Create scatter plot
        scatter = axes[1].scatter(malignant_df['correlation'], malignant_df['abs_correlation'],
                                c=malignant_df['p_value'], cmap='viridis_r', 
                                s=100, alpha=0.7, edgecolors='black')
        
        # Add text annotations for significant correlations
        for _, row in malignant_df.iterrows():
            if row['p_value'] < 0.01 and row['abs_correlation'] > 0.4:
                axes[1].annotate(f"{row['cell_type']}\n{row['pathway']}", 
                               (row['correlation'], row['abs_correlation']),
                               xytext=(5, 5), textcoords='offset points', 
                               fontsize=8, alpha=0.8)
        
        axes[1].set_xlabel('Correlation Coefficient')
        axes[1].set_ylabel('Absolute Correlation')
        axes[1].set_title('Malignant Cell-Pathway Correlations\n(Color: p-value)')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=axes[1])
        cbar.set_label('p-value')
        
        # Add significance thresholds
        axes[1].axhline(y=0.3, color='red', linestyle='--', alpha=0.5, label='|r| = 0.3')
        axes[1].axvline(x=0, color='black', linestyle='-', alpha=0.3)
        axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"cell_pathway_correlations_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()