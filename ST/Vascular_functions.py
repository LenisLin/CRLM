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
                    ncols=3, size=1.3,
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

def visualize_vascular_tumor_regions(adata, slide_dir, slide_name):
    """
    Visualize vascular hotspots, tumor regions, and spatial zones
    Similar to visualize_malignant_distribution but for vessel analysis results
    """
    
    print(f"  Creating vascular-tumor region visualizations for {slide_name}...")
    
    # Check if required columns exist
    required_cols = ['Is_Tumor_Region', 'Is_Vascular_Hotspot', 'Spatial_Zone', 'Distance_to_Vessel_um']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"  Warning: Missing columns {missing_cols}. Skipping visualization.")
        return
    
    # Create figure with 2x2 subplots
    with mpl.rc_context({'axes.facecolor': 'black', 'figure.figsize': [12, 10]}):
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        
        # 1. Tumor Regions (boolean)
        tumor_values = adata.obs['Is_Tumor_Region'].astype(int)
        adata.obs['tumor_regions_viz'] = tumor_values
        
        sc.pl.spatial(adata, 
                     color='tumor_regions_viz',
                     cmap='Reds',
                     ax=axes[0],
                     size=1.3,
                     img_key='hires',
                     vmin=0, vmax=1,
                     show=False,
                     colorbar_loc=None)
        axes[0].set_title(f'Tumor Regions\n({tumor_values.sum()} spots)', 
                         color='white', fontsize=12, fontweight='bold')
        
        # 2. Vascular Hotspots (boolean)
        hotspot_values = adata.obs['Is_Vascular_Hotspot'].astype(int)
        adata.obs['vascular_hotspots_viz'] = hotspot_values
        
        sc.pl.spatial(adata,
                     color='vascular_hotspots_viz', 
                     cmap='Blues',
                     ax=axes[1],
                     size=1.3,
                     img_key='hires',
                     vmin=0, vmax=1,
                     show=False,
                     colorbar_loc=None)
        axes[1].set_title(f'Vascular Hotspots\n({hotspot_values.sum()} spots)', 
                         color='white', fontsize=12, fontweight='bold')
        
        # 3. Spatial Zones (categorical)
        # Create numeric mapping for spatial zones
        zone_mapping = {
            'Perivascular': 3,
            'Intermediate': 2, 
            'Distant': 1
        }
        
        zone_numeric = adata.obs['Spatial_Zone'].map(zone_mapping)
        adata.obs['spatial_zones_viz'] = zone_numeric
        
        # Custom colormap for zones
        zone_colors = ['#3498DB', '#F39C12', '#E74C3C']  # Non_Tumor, Distant, Intermediate, Perivascular
        zone_cmap = mpl.colors.ListedColormap(zone_colors)
        
        sc.pl.spatial(adata,
                     color='spatial_zones_viz',
                     cmap=zone_cmap,
                     ax=axes[2], 
                     size=1.3,
                     img_key='hires',
                     vmin=0, vmax=3,
                     show=False,
                     colorbar_loc=None)
        
        # Add custom legend for zones
        zone_labels = ['Distant', 'Intermediate', 'Perivascular']
        zone_handles = [mpl.patches.Patch(color=zone_colors[i], label=zone_labels[i]) 
                       for i in range(len(zone_labels))]
        axes[2].legend(handles=zone_handles, loc='upper right', fontsize=8)
        axes[2].set_title('Spatial Zones', color='white', fontsize=12, fontweight='bold')
        
        # 4. Distance to Vessels (continuous)
        distance_values = adata.obs['Distance_to_Vessel_um']
        # Cap extreme values for better visualization
        distance_cap = np.percentile(distance_values[distance_values > 0], 95)
        distance_capped = np.clip(distance_values, 0, distance_cap)
        adata.obs['distance_vessels_viz'] = distance_capped
        
        sc.pl.spatial(adata,
                     color='distance_vessels_viz',
                     cmap='viridis_r',  # Reversed so close = hot colors
                     ax=axes[3],
                     size=1.3, 
                     img_key='hires',
                     vmin=0, 
                     vmax=distance_cap,
                     show=False,
                     colorbar_loc=None)
        axes[3].set_title(f'Distance to Vessels\n(max: {distance_cap:.0f}μm)', 
                         color='white', fontsize=12, fontweight='bold')
        
        # Overall styling
        for ax in axes:
            ax.set_facecolor('black')
            ax.tick_params(colors='white')
            for spine in ax.spines.values():
                spine.set_color('white')
        
        # Add treatment info
        treatment = adata.obs['treatment_status'].iloc[0] if 'treatment_status' in adata.obs.columns else 'Unknown'
        plt.suptitle(f'Vascular-Tumor Analysis: {slide_name}\nTreatment: {treatment.replace("_", " ").title()}', 
                    color='white', fontsize=14, fontweight='bold', y=0.95)
        
        plt.tight_layout(rect=[0, 0, 1, 0.93])
        
        # Save figure
        plt.savefig(f"{slide_dir}/vascular_tumor_regions_{slide_name}.pdf", 
                   bbox_inches='tight', dpi=300, facecolor='black')
        plt.show()
    
    # Print summary statistics
    print(f"  Summary for {slide_name}:")
    print(f"    Total spots: {len(adata.obs)}")
    print(f"    Tumor regions: {tumor_values.sum()} ({tumor_values.mean()*100:.1f}%)")
    print(f"    Vascular hotspots: {hotspot_values.sum()} ({hotspot_values.mean()*100:.1f}%)")
    
    return None
#%% Step 4: Vascular-Tumor Subtype Spatial Associations

def analyze_tumor_vessel_spatial_associations(adata):
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

def create_distance_profiles(adata, slide_name, bin_size=55, max_distance=500, min_spots=3):
    """Create distance-dependent association curves for tumor subtypes"""
    
    print(f"  Creating distance profiles for {slide_name}...")
    
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    tumor_mask = adata.obs['Is_Tumor_Region']
    tumor_data = adata.obs[tumor_mask]
    
    # Create distance bins (55μm resolution)
    distance_bins = np.arange(0, max_distance + bin_size, bin_size)
    distance_profiles = {}
    
    for tc_type in tumor_types:
        bin_profiles = []
        bin_distances = []
        bin_ranges = []
        bin_counts = []
        bin_stds = []
        
        total_spots_processed = 0

        # Process each distance range
        for i in range(len(distance_bins) - 1):
            bin_min = distance_bins[i]
            bin_max = distance_bins[i + 1]
            
            # Cap at max_distance
            if bin_min >= max_distance:
                break
            if bin_max > max_distance:
                bin_max = max_distance
            
            # Create mask for current distance range
            bin_mask = (tumor_data['Distance_to_Vessel_um'] >= bin_min) & \
                      (tumor_data['Distance_to_Vessel_um'] < bin_max)
            
            spot_count = bin_mask.sum()
            total_spots_processed += spot_count
            
            if spot_count >= min_spots:
                abundance_values = tumor_data.loc[bin_mask, tc_type]
                mean_abundance = abundance_values.mean()
                std_abundance = abundance_values.std()
                bin_center = (bin_min + bin_max) / 2
                
                bin_profiles.append(mean_abundance)
                bin_distances.append(bin_center)
                bin_ranges.append(f"{bin_min}-{bin_max}")
                bin_counts.append(spot_count)
                bin_stds.append(std_abundance)
        
        if len(bin_profiles) > 2:
            distance_profiles[tc_type] = {
                'distances': np.array(bin_distances),
                'abundances': np.array(bin_profiles),
                'std_abundances': np.array(bin_stds),
                'ranges': bin_ranges,
                'spot_counts': np.array(bin_counts),
                'bin_edges': distance_bins[:len(bin_profiles)+1],
                'total_spots': total_spots_processed,
                'coverage_rate': total_spots_processed / len(tumor_data) * 100
            }
    
    return distance_profiles

def visualize_tumor_vessel_associations(adata, slide_name, distance_profiles, associations, slide_dir):
    """
    Enhanced distance profiles with error bars and publication-quality styling
    """
    
    if not distance_profiles:
        print(f"No distance profiles to visualize for {slide_name}")
        return
    
    treatment = adata.obs['treatment_status'].iloc[0] if 'treatment_status' in adata.obs.columns else 'Unknown'
    
    # Create figure with better proportions
    fig, (ax_main, ax_stats) = plt.subplots(1, 2, figsize=(16, 8), 
                                           gridspec_kw={'width_ratios': [3, 1]})
    
    # Enhanced color scheme for tumor types
    tumor_type_colors = {
        'TC_Glycolysis': '#E74C3C',    # Red - keep original color
        'TC_EMT': '#9B59B6',           # Purple 
        'TC_Proliferation': '#3498DB',  # Blue
        'TC_Quiescent': '#2ECC71',     # Green
        'TC_LipidMeta': '#F39C12'      # Orange
    }
    
    # Track data for statistics
    profile_stats = []
    legend_elements = []
    
    print(f"Plotting distance profiles for {slide_name}:")
    
    for tc_type, profile in distance_profiles.items():
        distances = profile['distances']
        abundances = profile['abundances']
        
        if len(distances) < 2:
            continue
            
        # Get color and styling (uniform for all types)
        color = tumor_type_colors.get(tc_type, '#34495E')
        linewidth = 2  # Same for all types
        markersize = 6  # Same for all types
        alpha = 0.8    # Same for all types
        
        # Plot main line with points
        main_line = ax_main.plot(distances, abundances, 'o-', 
                               color=color, linewidth=linewidth, 
                               markersize=markersize, alpha=alpha,
                               label=tc_type.replace('TC_', ''), 
                               markerfacecolor=color, 
                               markeredgecolor='white',
                               markeredgewidth=1,
                               zorder=3)
        
        # Add error bars if available
        if 'std_abundances' in profile and profile['std_abundances'] is not None:
            std_abundances = profile['std_abundances']
            ax_main.errorbar(distances, abundances, yerr=std_abundances,
                           color=color, alpha=0.6, capsize=4, capthick=2,
                           linestyle='none', zorder=2)
            
            print(f"  {tc_type}: {len(distances)} points, range {abundances.min():.3f}-{abundances.max():.3f}")
            
            # Use pre-calculated correlation from associations
            if tc_type in associations:
                correlation = associations[tc_type]['distance_correlation']
                vessel_behavior = "vessel-seeking" if correlation < -0.2 else \
                                "vessel-avoiding" if correlation > 0.2 else "neutral"
                
                print(f"    Correlation with distance: {correlation:+.3f} ({vessel_behavior})")
                
                profile_stats.append({
                    'tumor_type': tc_type.replace('TC_', ''),
                    'n_points': len(distances),
                    'mean_abundance': abundances.mean(),
                    'abundance_range': abundances.max() - abundances.min(),
                    'distance_correlation': correlation,
                    'vessel_behavior': vessel_behavior,
                    'color': color
                })
            else:
                print(f"    Warning: No association data found for {tc_type}")
        
        # Store for legend
        legend_elements.append(mpatches.Patch(color=color, 
                                            label=tc_type.replace('TC_', '')))
    
    # Enhanced zone visualization
    zone_colors = {'Perivascular': '#E74C3C', 'Intermediate': '#F39C12', 'Distant': '#3498DB'}
    zone_alphas = {'Perivascular': 0.15, 'Intermediate': 0.10, 'Distant': 0.08}
    
    # Zone backgrounds
    ax_main.axvspan(0, 100, alpha=zone_alphas['Perivascular'], 
                   color=zone_colors['Perivascular'], zorder=0, label='Perivascular')
    ax_main.axvspan(100, 200, alpha=zone_alphas['Intermediate'], 
                   color=zone_colors['Intermediate'], zorder=0, label='Intermediate')
    ax_main.axvspan(200, 500, alpha=zone_alphas['Distant'], 
                   color=zone_colors['Distant'], zorder=0, label='Distant')
    
    # Zone boundary lines
    ax_main.axvline(x=100, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    ax_main.axvline(x=200, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Dynamic y-limits for annotations
    y_max = ax_main.get_ylim()[1]
    
    # Enhanced annotations
    ax_main.text(50, y_max*0.95, 'Perivascular\n(0-100μm)', 
                ha='center', va='top', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    ax_main.text(150, y_max*0.95, 'Intermediate\n(100-200μm)', 
                ha='center', va='top', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    ax_main.text(350, y_max*0.95, 'Distant\n(>200μm)', 
                ha='center', va='top', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    
    # Enhanced main plot styling
    ax_main.set_xlabel('Distance to Nearest Vessel (μm)', fontsize=12, fontweight='bold')
    ax_main.set_ylabel('Tumor Cell Abundance', fontsize=12, fontweight='bold')
    ax_main.set_title(f'Vessel-Distance Profiles: {slide_name}\n({treatment.replace("_", " ").title()})', 
                     fontsize=14, fontweight='bold', pad=20)
    
    # Enhanced grid and spines
    ax_main.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.spines['left'].set_linewidth(1.5)
    ax_main.spines['bottom'].set_linewidth(1.5)
    
    # Add legend to main plot
    ax_main.legend(loc='upper right', fontsize=10, frameon=True, 
                  fancybox=True, shadow=True)
    
    # Statistics panel
    ax_stats.axis('off')
    
    if profile_stats:
        # Title
        ax_stats.text(0.05, 0.95, 'Vessel Association Analysis', 
                     fontsize=12, fontweight='bold', transform=ax_stats.transAxes)
        
        # Table header
        y_pos = 0.85
        ax_stats.text(0.05, y_pos, 'Tumor Type', fontsize=10, fontweight='bold', 
                     transform=ax_stats.transAxes)
        ax_stats.text(0.65, y_pos, 'Behavior', fontsize=10, fontweight='bold', 
                     transform=ax_stats.transAxes)
        
        # Add separator line
        y_pos -= 0.05
        ax_stats.plot([0.05, 0.95], [y_pos, y_pos], 'k-', alpha=0.3, 
                     transform=ax_stats.transAxes)
        
        # Statistics for each tumor type
        for stat in profile_stats:
            y_pos -= 0.08
            
            # Color indicator
            rect = Rectangle((0.02, y_pos-0.01), 0.02, 0.04, 
                           facecolor=stat['color'], transform=ax_stats.transAxes)
            ax_stats.add_patch(rect)
            
            # Tumor type name
            ax_stats.text(0.05, y_pos, stat['tumor_type'], fontsize=9, 
                         transform=ax_stats.transAxes)
            
            # Behavior with color coding
            behavior = stat['vessel_behavior']
            corr = stat['distance_correlation']
            
            if behavior == 'vessel-seeking':
                behavior_color = VESSEL_SEEKING_COLOR
                behavior_text = f'Seeking\n(r={corr:+.2f})'
            elif behavior == 'vessel-avoiding':
                behavior_color = VESSEL_AVOIDING_COLOR  
                behavior_text = f'Avoiding\n(r={corr:+.2f})'
            else:
                behavior_color = NEUTRAL_COLOR
                behavior_text = f'Neutral\n(r={corr:+.2f})'
            
            ax_stats.text(0.65, y_pos, behavior_text, fontsize=8, 
                         color=behavior_color, transform=ax_stats.transAxes)
        
        # Add interpretation guide
        y_pos -= 0.15
        ax_stats.text(0.05, y_pos, 'Interpretation:', fontsize=10, fontweight='bold',
                     transform=ax_stats.transAxes)
        
        y_pos -= 0.08
        ax_stats.text(0.05, y_pos, '• r < -0.2: Vessel-seeking', fontsize=9,
                     color=VESSEL_SEEKING_COLOR, transform=ax_stats.transAxes)
        
        y_pos -= 0.06
        ax_stats.text(0.05, y_pos, '• r > +0.2: Vessel-avoiding', fontsize=9,
                     color=VESSEL_AVOIDING_COLOR, transform=ax_stats.transAxes)
        
        y_pos -= 0.06
        ax_stats.text(0.05, y_pos, '• |r| ≤ 0.2: Neutral', fontsize=9,
                     color=NEUTRAL_COLOR, transform=ax_stats.transAxes)
        
        # Data quality info
        y_pos -= 0.12
        total_spots = list(distance_profiles.values())[0]['total_spots']
        coverage = list(distance_profiles.values())[0]['coverage_rate']
        n_points = len(list(distance_profiles.values())[0]['distances'])
        
        ax_stats.text(0.05, y_pos, f'Data Quality:', fontsize=10, fontweight='bold',
                     transform=ax_stats.transAxes)
        y_pos -= 0.06
        ax_stats.text(0.05, y_pos, f'Spots: {total_spots}', fontsize=9,
                     transform=ax_stats.transAxes)
        y_pos -= 0.05
        ax_stats.text(0.05, y_pos, f'Coverage: {coverage:.1f}%', fontsize=9,
                     transform=ax_stats.transAxes)
        y_pos -= 0.05
        ax_stats.text(0.05, y_pos, f'Distance bins: {n_points}', fontsize=9,
                     transform=ax_stats.transAxes)
    
    plt.tight_layout()
    
    # Save figure
    save_path = os.path.join(slide_dir, f"enhanced_distance_profiles_{slide_name}.pdf")
    plt.savefig(save_path, bbox_inches='tight', dpi=300, facecolor='white')
    print(f"Distance profiles saved to: {save_path}")
    
    plt.show()
    
    return None
#%% Step 5: Vascular Hotspot Co-occurrence Analysis (Metastasis Mechanism)
def analyze_vascular_hotspot_cooccurrence(adata, slide_name):
    """
    Analyze tumor-vessel co-occurrence in same spots (metastasis mechanism)
    
    Fixed version with corrected statistical logic
    """
    
    print(f"  Analyzing vascular hotspot co-occurrence (metastasis potential)...")
    
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    
    # Focus on tumor regions only
    tumor_mask = adata.obs['Is_Tumor_Region']
    hotspot_mask = adata.obs['Is_Vascular_Hotspot']
    
    # Restrict all analysis to tumor regions only
    tumor_data = adata.obs[tumor_mask].copy()
    
    if len(tumor_data) == 0:
        print(f"    Warning: No tumor regions found in {slide_name}")
        return {}
    
    # Get hotspot status within tumor regions
    hotspot_in_tumor = tumor_data['Is_Vascular_Hotspot']
    
    cooccurrence_results = {}
    
    print(f"    Co-occurrence analysis (same 100μm spot within tumor regions):")
    print(f"    Total tumor spots: {len(tumor_data)}")
    print(f"    Vascular hotspots in tumor: {hotspot_in_tumor.sum()} ({hotspot_in_tumor.mean()*100:.1f}%)")
    
    for tc_type in tumor_types:
        if tc_type not in tumor_data.columns:
            continue
            
        # Define tumor-positive spots (within tumor regions only)
        tumor_values = tumor_data[tc_type]
        tumor_positive_threshold = np.percentile(tumor_values, 90)  # Top 10% within tumor regions
        tumor_positive_in_tumor = tumor_values >= tumor_positive_threshold
        
        # Four categories (all within tumor regions)
        cooccurrence_mask = hotspot_in_tumor & tumor_positive_in_tumor    # Both high
        hotspot_only = hotspot_in_tumor & (~tumor_positive_in_tumor)      # Vessels only
        tumor_only = (~hotspot_in_tumor) & tumor_positive_in_tumor        # Tumor only  
        neither = (~hotspot_in_tumor) & (~tumor_positive_in_tumor)        # Neither
        
        # Counts
        n_cooccurrence = cooccurrence_mask.sum()
        n_hotspot_only = hotspot_only.sum()
        n_tumor_only = tumor_only.sum()
        n_neither = neither.sum()
        n_total_tumor = len(tumor_data)
        
        # Verify counts add up
        total_check = n_cooccurrence + n_hotspot_only + n_tumor_only + n_neither
        assert total_check == n_total_tumor, f"Count mismatch: {total_check} != {n_total_tumor}"
        
        # Calculate metastasis potential metrics
        if n_total_tumor > 0:
            cooccurrence_rate = n_cooccurrence / n_total_tumor * 100
            
            # Expected co-occurrence under independence (within tumor regions)
            hotspot_rate = hotspot_in_tumor.mean()
            tumor_positive_rate = tumor_positive_in_tumor.mean()
            expected_cooccurrence = hotspot_rate * tumor_positive_rate * n_total_tumor
            
            # Enrichment ratio
            enrichment_ratio = n_cooccurrence / (expected_cooccurrence + 1e-6)
            
            # Fisher's exact test for significance
            # Test: Are vascular hotspots and tumor-positive status independent?
            contingency_table = [[n_cooccurrence, n_tumor_only],      # Tumor-positive row
                               [n_hotspot_only, n_neither]]           # Tumor-negative row
            
            try:
                odds_ratio, p_value = fisher_exact(contingency_table)
            except ValueError as e:
                print(f"      Warning: Fisher's test failed for {tc_type}: {e}")
                odds_ratio, p_value = np.nan, 1.0
            
            significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                        
            cooccurrence_results[tc_type] = {
                'cooccurrence_spots': n_cooccurrence,
                'cooccurrence_rate': cooccurrence_rate,
                'enrichment_ratio': enrichment_ratio,
                'odds_ratio': odds_ratio,
                'p_value': p_value,
                'significance': significance,
                'total_tumor_spots': n_total_tumor,
                'expected_cooccurrence': expected_cooccurrence,
                'hotspot_rate': hotspot_rate * 100,
                'tumor_positive_rate': tumor_positive_rate * 100,
                'contingency_table': contingency_table,
                'counts': {
                    'cooccurrence': n_cooccurrence,
                    'hotspot_only': n_hotspot_only, 
                    'tumor_only': n_tumor_only,
                    'neither': n_neither
                }
            }
            
            print(f"      {tc_type}:")
            print(f"        Co-occurrence: {n_cooccurrence}/{n_total_tumor} spots ({cooccurrence_rate:.1f}%)")
            print(f"        Expected: {expected_cooccurrence:.1f}, Enrichment: {enrichment_ratio:.2f}x")
            print(f"        Odds Ratio: {odds_ratio:.2f}, p-value: {p_value:.3f} {significance}")
    
    return cooccurrence_results

def visualize_hotspot_cooccurrence(adata, slide_name, cooccurrence_results, slide_dir):
    """
    Enhanced co-occurrence visualization with meaningful biological interpretation
    """
    
    if not cooccurrence_results:
        print(f"No co-occurrence results to visualize for {slide_name}")
        return
    
    treatment = adata.obs['treatment_status'].iloc[0] if 'treatment_status' in adata.obs.columns else 'Unknown'
    
    # Create figure with 2x2 layout
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    # Extract data for visualization
    tumor_types = []
    enrichment_ratios = []
    cooccurrence_rates = []
    p_values = []
    significance_levels = []
    observed_counts = []
    expected_counts = []
    
    for tc_type, results in cooccurrence_results.items():
        tumor_types.append(tc_type.replace('TC_', ''))
        enrichment_ratios.append(results['enrichment_ratio'])
        cooccurrence_rates.append(results['cooccurrence_rate'])
        p_values.append(results['p_value'])
        significance_levels.append(results['significance'])
        observed_counts.append(results['cooccurrence_spots'])
        expected_counts.append(results['expected_cooccurrence'])
    
    # Color scheme based on vessel association behavior
    colors = []
    for ratio in enrichment_ratios:
        if ratio > 1.2:
            colors.append('#E74C3C')  # Red: Vessel-seeking
        elif ratio < 0.8:
            colors.append('#3498DB')  # Blue: Vessel-avoiding
        else:
            colors.append('#95A5A6')  # Gray: Neutral
    
    # 1. Enhanced Enrichment Ratio Plot
    bars1 = axes[0].bar(tumor_types, enrichment_ratios, color=colors, 
                       alpha=0.8, edgecolor='white', linewidth=1.5)
    
    # Add reference lines
    axes[0].axhline(y=1, color='black', linestyle='-', linewidth=2, alpha=0.8, label='Expected (1.0x)')
    axes[0].axhline(y=1.2, color='#E74C3C', linestyle='--', alpha=0.7, label='Vessel-seeking (>1.2x)')
    axes[0].axhline(y=0.8, color='#3498DB', linestyle='--', alpha=0.7, label='Vessel-avoiding (<0.8x)')
    
    # Add value labels with significance
    for i, (bar, ratio, sig) in enumerate(zip(bars1, enrichment_ratios, significance_levels)):
        height = bar.get_height()
        y_pos = height + 0.05 if height > 0 else height - 0.05
        va = 'bottom' if height > 0 else 'top'
        
        # Format significance
        sig_marker = {"***": "***", "**": "**", "*": "*", "ns": ""}[sig]
        label_text = f'{ratio:.2f}\n{sig_marker}'
        
        axes[0].text(bar.get_x() + bar.get_width()/2., y_pos, label_text,
                    ha='center', va=va, fontweight='bold', fontsize=10)
    
    axes[0].set_title('Vessel Association Pattern\n(Enrichment vs Expected)', 
                     fontsize=14, fontweight='bold', pad=20)
    axes[0].set_ylabel('Enrichment Ratio\n(Observed/Expected)', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Tumor Cell Type', fontsize=12, fontweight='bold')
    axes[0].tick_params(axis='x', rotation=45, labelsize=11)
    axes[0].legend(loc='upper right', fontsize=10)
    axes[0].grid(True, alpha=0.3, axis='y')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    
    # 2. Observed vs Expected Co-occurrence
    x_pos = np.arange(len(tumor_types))
    width = 0.35
    
    bars2a = axes[1].bar(x_pos - width/2, observed_counts, width, label='Observed', 
                        color='#E74C3C', alpha=0.8)
    bars2b = axes[1].bar(x_pos + width/2, expected_counts, width, label='Expected', 
                        color='#95A5A6', alpha=0.6)
    
    # Add value labels
    for bars in [bars2a, bars2b]:
        for bar in bars:
            height = bar.get_height()
            axes[1].text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{height:.1f}', ha='center', va='bottom', fontsize=9)
    
    axes[1].set_title('Co-occurrence Spot Counts\n(Observed vs Expected)', 
                     fontsize=14, fontweight='bold', pad=20)
    axes[1].set_ylabel('Number of Co-occurring Spots', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Tumor Cell Type', fontsize=12, fontweight='bold')
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(tumor_types, rotation=45)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3, axis='y')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    
    # 3. Statistical Significance Landscape
    # Create a heatmap-style visualization of p-values
    p_value_matrix = np.array(p_values).reshape(1, -1)
    
    # Create custom colormap for p-values
    from matplotlib.colors import ListedColormap
    p_colors = ['#E74C3C', '#F39C12', '#F1C40F', '#95A5A6']  # red, orange, yellow, gray
    p_cmap = ListedColormap(p_colors)
    
    # Convert p-values to categories
    p_categories = []
    for p in p_values:
        if p < 0.001:
            p_categories.append(0)  # Highly significant
        elif p < 0.01:
            p_categories.append(1)  # Significant
        elif p < 0.05:
            p_categories.append(2)  # Marginally significant
        else:
            p_categories.append(3)  # Not significant
    
    p_category_matrix = np.array(p_categories).reshape(1, -1)
    
    im = axes[2].imshow(p_category_matrix, cmap=p_cmap, aspect='auto', vmin=0, vmax=3)
    
    # Customize the heatmap
    axes[2].set_xticks(range(len(tumor_types)))
    axes[2].set_xticklabels(tumor_types, rotation=45)
    axes[2].set_yticks([0])
    axes[2].set_yticklabels(['P-value'])
    axes[2].set_title('Statistical Significance\n(Fisher\'s Exact Test)', 
                     fontsize=14, fontweight='bold', pad=20)
    
    # Add text annotations with actual p-values
    for i, (p_val, sig) in enumerate(zip(p_values, significance_levels)):
        color = 'white' if p_categories[i] < 2 else 'black'
        axes[2].text(i, 0, f'{p_val:.3f}\n({sig})', ha='center', va='center', 
                    fontweight='bold', fontsize=10, color=color)
    
    # Add colorbar legend
    p_labels = ['p<0.001\n(***)', 'p<0.01\n(**)', 'p<0.05\n(*)', 'p≥0.05\n(ns)']
    legend_patches = [mpatches.Patch(color=p_colors[i], label=p_labels[i]) for i in range(4)]
    axes[2].legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left')

    # 4. Precision vs Recall Analysis
    scatter1 = axes[3].scatter(cooccurrence_rates, enrichment_ratios, c=colors, s=100, 
                              edgecolors='black', linewidths=1, alpha=0.8)
    for i, tumor_type in enumerate(tumor_types):
        axes[3].annotate(tumor_type, (cooccurrence_rates[i], enrichment_ratios[i]), 
                        xytext=(5, 5), textcoords='offset points', fontsize=10)
    axes[3].axhline(y=1, color='black', linestyle='-', linewidth=2, alpha=0.8)
    axes[3].set_xlabel('Co-occurrence Rate (%)')
    axes[3].set_ylabel('Enrichment Ratio')
    axes[3].set_title('Option 1: Spatial Association Summary')
    axes[3].grid(True, alpha=0.3)

    # Overall styling
    plt.suptitle(f'Vessel-Tumor Co-occurrence Analysis: {slide_name}\n'
                f'Treatment: {treatment.replace("_", " ").title()}', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save figure
    save_path = os.path.join(slide_dir, f"enhanced_cooccurrence_analysis_{slide_name}.pdf")
    plt.savefig(save_path, bbox_inches='tight', dpi=300, facecolor='white')
    print(f"Co-occurrence analysis saved to: {save_path}")
    
    plt.show()
    
    return None

#%% Step 6: Treatment Effect Analysis

def analyze_treatment_effects_single_slide(adata, slide_name, associations, cooccurrence_results, distance_profiles=None):
    """
    Extract treatment-relevant metrics for cross-slide comparison
    Uses pre-calculated results from associations, cooccurrence_results, and distance_profiles
    """
    
    treatment = adata.obs['treatment_status'].iloc[0]
    
    print(f"  Extracting treatment metrics for {slide_name} ({treatment})...")
    
    # Basic treatment metrics
    treatment_metrics = {
        'slide_name': slide_name,
        'treatment': treatment,
        'n_spots': len(adata.obs),
        'tumor_spots': adata.obs['Is_Tumor_Region'].sum(),
        'vascular_hotspots': adata.obs['Is_Vascular_Hotspot'].sum()
    }
    
    # Calculate basic coverage metrics
    tumor_mask = adata.obs['Is_Tumor_Region']
    hotspot_mask = adata.obs['Is_Vascular_Hotspot']
    
    treatment_metrics['tumor_coverage_pct'] = (tumor_mask.sum() / len(adata.obs)) * 100
    treatment_metrics['vascular_density_pct'] = (hotspot_mask.sum() / len(adata.obs)) * 100
    treatment_metrics['vascular_in_tumor_pct'] = ((hotspot_mask & tumor_mask).sum() / tumor_mask.sum()) * 100 if tumor_mask.sum() > 0 else 0
    
    # 1. Overall tumor abundances (basic stats only)
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    for tc_type in tumor_types:
        if tumor_mask.sum() > 0:
            treatment_metrics[f'{tc_type}_mean_abundance'] = adata.obs.loc[tumor_mask, tc_type].mean()
            treatment_metrics[f'{tc_type}_std_abundance'] = adata.obs.loc[tumor_mask, tc_type].std()
            treatment_metrics[f'{tc_type}_max_abundance'] = adata.obs.loc[tumor_mask, tc_type].max()
    
    # 2. Zone-specific abundances from associations (DIRECT USE)
    print(f"    Using pre-calculated zone abundances from associations...")
    for tc_type, assoc in associations.items():
        # Distance correlation (already calculated)
        treatment_metrics[f'{tc_type}_vessel_correlation'] = assoc['distance_correlation']
        
        # Vessel behavior classification
        correlation = assoc['distance_correlation']
        if correlation < -0.2:
            behavior = 'vessel_seeking'
        elif correlation > 0.2:
            behavior = 'vessel_avoiding'
        else:
            behavior = 'neutral'
        treatment_metrics[f'{tc_type}_vessel_behavior'] = behavior
        
        # Zone abundances (directly from associations)
        if 'zone_abundances' in assoc:
            zone_abundances = assoc['zone_abundances']
            
            # Store each zone abundance directly
            for zone, abundance in zone_abundances.items():
                treatment_metrics[f'{tc_type}_{zone}_abundance'] = abundance
            
            # Calculate derived zone metrics
            if 'Perivascular' in zone_abundances and 'Distant' in zone_abundances:
                peri_abundance = zone_abundances['Perivascular']
                distant_abundance = zone_abundances['Distant']
                gradient = (distant_abundance - peri_abundance) / (peri_abundance + 1e-6)
                treatment_metrics[f'{tc_type}_peri_to_distant_gradient'] = gradient
                
            # Zone preference metrics
            if zone_abundances:
                max_zone = max(zone_abundances, key=zone_abundances.get)
                max_abundance = zone_abundances[max_zone]
                min_abundance = min(zone_abundances.values())
                zone_selectivity = (max_abundance - min_abundance) / (max_abundance + 1e-6)
                treatment_metrics[f'{tc_type}_zone_selectivity'] = zone_selectivity
                treatment_metrics[f'{tc_type}_preferred_zone'] = max_zone
    
    # 3. Distance bin abundances from distance profiles (DIRECT USE)
    if distance_profiles:
        print(f"    Using pre-calculated distance profiles...")
        for tc_type, profile in distance_profiles.items():
            # Distance profile data (already calculated)
            if 'distances' in profile and 'abundances' in profile:
                distances = profile['distances']
                abundances = profile['abundances']
                ranges = profile.get('ranges', [])
                
                # Store distance bin abundances using pre-calculated ranges
                for i, (dist, abund) in enumerate(zip(distances, abundances)):
                    if i < len(ranges):
                        range_label = ranges[i]  # Use pre-calculated range labels
                    else:
                        range_label = f"{int(dist-27.5)}-{int(dist+27.5)}"  # Fallback
                    
                    treatment_metrics[f'{tc_type}_dist_{range_label}_abundance'] = abund
                
                # Distance profile summary statistics (direct from arrays)
                treatment_metrics[f'{tc_type}_dist_profile_range'] = abundances.max() - abundances.min()
                treatment_metrics[f'{tc_type}_dist_profile_mean'] = abundances.mean()
                treatment_metrics[f'{tc_type}_dist_profile_std'] = abundances.std()
                
                # Peak distance (where abundance is highest)
                peak_idx = np.argmax(abundances)
                treatment_metrics[f'{tc_type}_peak_distance'] = distances[peak_idx]
                treatment_metrics[f'{tc_type}_peak_abundance'] = abundances[peak_idx]
                
                # Additional pre-calculated metrics from distance profiles
                if 'total_spots' in profile:
                    treatment_metrics[f'{tc_type}_dist_total_spots'] = profile['total_spots']
                if 'coverage_rate' in profile:
                    treatment_metrics[f'{tc_type}_dist_coverage_rate'] = profile['coverage_rate']
    
    # 4. Spatial zone distributions (minimal calculation)
    if 'Spatial_Zone' in adata.obs.columns:
        # Overall zone distribution
        zone_counts = adata.obs['Spatial_Zone'].value_counts()
        total_spots = len(adata.obs)
        
        for zone in ['Perivascular', 'Intermediate', 'Distant']:
            count = zone_counts.get(zone, 0)
            treatment_metrics[f'{zone}_zone_spots'] = count
            treatment_metrics[f'{zone}_zone_pct'] = (count / total_spots) * 100
        
        # Tumor-specific zone distributions
        if tumor_mask.sum() > 0:
            tumor_zones = adata.obs.loc[tumor_mask, 'Spatial_Zone'].value_counts()
            tumor_total = tumor_mask.sum()
            
            for zone in ['Perivascular', 'Intermediate', 'Distant']:
                count = tumor_zones.get(zone, 0)
                treatment_metrics[f'{zone}_tumor_spots'] = count
                treatment_metrics[f'{zone}_tumor_pct'] = (count / tumor_total) * 100
    
    # 5. Co-occurrence metrics (DIRECT USE)
    print(f"    Using pre-calculated co-occurrence results...")
    for tc_type, cooc in cooccurrence_results.items():
        # Core co-occurrence metrics (already calculated)
        treatment_metrics[f'{tc_type}_cooccurrence_rate'] = cooc['cooccurrence_rate']
        treatment_metrics[f'{tc_type}_enrichment_ratio'] = cooc['enrichment_ratio']
        treatment_metrics[f'{tc_type}_cooccurrence_pvalue'] = cooc['p_value']
        treatment_metrics[f'{tc_type}_cooccurrence_significant'] = cooc['significance'] != 'ns'
        treatment_metrics[f'{tc_type}_odds_ratio'] = cooc['odds_ratio']
        
        # Additional pre-calculated metrics
        if 'expected_cooccurrence' in cooc:
            treatment_metrics[f'{tc_type}_expected_cooccurrence'] = cooc['expected_cooccurrence']
        if 'hotspot_rate' in cooc:
            treatment_metrics[f'{tc_type}_hotspot_rate'] = cooc['hotspot_rate']
        if 'tumor_positive_rate' in cooc:
            treatment_metrics[f'{tc_type}_tumor_positive_rate'] = cooc['tumor_positive_rate']
        if 'precision' in cooc:
            treatment_metrics[f'{tc_type}_precision'] = cooc['precision']
        if 'recall' in cooc:
            treatment_metrics[f'{tc_type}_recall'] = cooc['recall']
        
        # Co-occurrence classification
        enrichment = cooc['enrichment_ratio']
        p_value = cooc['p_value']
        
        if p_value < 0.05:
            if enrichment > 1.2:
                cooc_class = 'significant_attraction'
            elif enrichment < 0.8:
                cooc_class = 'significant_avoidance'
            else:
                cooc_class = 'significant_neutral'
        else:
            cooc_class = 'not_significant'
        
        treatment_metrics[f'{tc_type}_cooccurrence_class'] = cooc_class
    
    # 6. Summary spatial metrics (using pre-calculated zone abundances)
    print(f"    Calculating summary spatial metrics from pre-calculated data...")
    
    # Spatial heterogeneity index using zone abundances from associations
    for tc_type in associations.keys():
        if 'zone_abundances' in associations[tc_type]:
            zone_abunds = list(associations[tc_type]['zone_abundances'].values())
            if len(zone_abunds) > 1:
                cv = np.std(zone_abunds) / (np.mean(zone_abunds) + 1e-6)
                treatment_metrics[f'{tc_type}_spatial_heterogeneity'] = cv
    
    # Dominant tumor type in each zone (using associations data)
    for zone in ['Perivascular', 'Intermediate', 'Distant']:
        zone_abundances = {}
        for tc_type in associations.keys():
            if 'zone_abundances' in associations[tc_type]:
                zone_abundances[tc_type] = associations[tc_type]['zone_abundances'].get(zone, 0)
        
        if zone_abundances:
            dominant_type = max(zone_abundances, key=zone_abundances.get)
            dominant_abundance = zone_abundances[dominant_type]
            treatment_metrics[f'{zone}_dominant_tumor_type'] = dominant_type
            treatment_metrics[f'{zone}_dominant_abundance'] = dominant_abundance
    
    print(f"    ✓ Extracted {len(treatment_metrics)} metrics for {slide_name}")
    
    return treatment_metrics


#%% STEP 8: IDENTIFY VESSEL-ASSOCIATED GENES

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

def calculate_pathway_scores(adata, vessel_pathways, method='ssgsea'):
    """
    Calculate pathway enrichment scores for each spot using ssGSEA-like algorithm
    """
    print(f"  Calculating pathway scores using {method} method...")
    
    pathway_scores = {}
    available_genes = adata.var.index.tolist()
    
    if hasattr(adata.X, 'toarray'):
        expr_matrix = adata.X.toarray()
    else:
        expr_matrix = adata.X
    
    n_spots, n_genes = expr_matrix.shape
    
    for pathway_name, pathway_genes in vessel_pathways.items():
        # Find available genes in the pathway
        available_pathway_genes = [gene for gene in pathway_genes if gene in available_genes]
        
        if len(available_pathway_genes) >= 3:  # Minimum genes for reliable score
            print(f"    {pathway_name}: {len(available_pathway_genes)}/{len(pathway_genes)} genes available")
            
            if method == 'ssgsea':
                # ssGSEA-like enrichment score calculation
                pathway_scores_list = []
                
                # Get indices of pathway genes
                pathway_indices = [adata.var.index.get_loc(gene) for gene in available_pathway_genes]
                pathway_set = set(pathway_indices)
                
                for spot_idx in range(n_spots):
                    # Get expression values for this spot
                    spot_expr = expr_matrix[spot_idx, :]
                    
                    # Rank genes by expression (descending order)
                    gene_ranks = np.argsort(-spot_expr)  # Negative for descending
                    
                    # Calculate enrichment score
                    enrichment_score = 0
                    max_enrichment = 0
                    min_enrichment = 0
                    
                    # Parameters for weighting
                    n_pathway_genes = len(pathway_indices)
                    n_total_genes = len(gene_ranks)
                    
                    # Walk through ranked gene list
                    for rank_pos, gene_idx in enumerate(gene_ranks):
                        if gene_idx in pathway_set:
                            # Gene is in pathway - increase score
                            enrichment_score += 1.0 / n_pathway_genes
                        else:
                            # Gene not in pathway - decrease score
                            enrichment_score -= 1.0 / (n_total_genes - n_pathway_genes)
                        
                        # Track maximum and minimum enrichment
                        if enrichment_score > max_enrichment:
                            max_enrichment = enrichment_score
                        if enrichment_score < min_enrichment:
                            min_enrichment = enrichment_score
                    
                    # Final enrichment score is the maximum absolute deviation
                    if abs(max_enrichment) > abs(min_enrichment):
                        final_score = max_enrichment
                    else:
                        final_score = min_enrichment
                    
                    pathway_scores_list.append(final_score)
                
                pathway_scores[pathway_name] = np.array(pathway_scores_list)
                
            elif method == 'mean':
                # Simple mean expression method (your original)
                pathway_indices = [adata.var.index.get_loc(gene) for gene in available_pathway_genes]
                pathway_expr = expr_matrix[:, pathway_indices]
                pathway_score = np.mean(pathway_expr, axis=1)
                pathway_scores[pathway_name] = pathway_score
            
            elif method == 'zscore_mean':
                # Z-score normalized mean
                pathway_indices = [adata.var.index.get_loc(gene) for gene in available_pathway_genes]
                pathway_expr = expr_matrix[:, pathway_indices]
                
                # Z-score normalize each gene across spots
                pathway_expr_zscore = (pathway_expr - np.mean(pathway_expr, axis=0)) / (np.std(pathway_expr, axis=0) + 1e-8)
                pathway_score = np.mean(pathway_expr_zscore, axis=1)
                pathway_scores[pathway_name] = pathway_score
    
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
                enrichment_ratio = (overlap_count / total_vessel_genes) / (total_genes_in_pathway / len(adata.var_names))  # Assume ~20k total genes
                
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
    Add statistical significance annotations using Wilcoxon test with BH correction
    """
    from statsmodels.stats.multitest import multipletests
    
    pathways = list(pathway_region_data.keys())
    x_positions = range(len(pathways))
    
    # First pass: collect all p-values
    p_values = []
    valid_comparisons = []
    
    for i, pathway in enumerate(pathways):
        region_data = pathway_region_data[pathway]
        
        # Compare Perivascular vs Distant (most biologically relevant)
        if (len(region_data['Perivascular']) > 0 and 
            len(region_data['Distant']) > 0):
            
            try:
                # Perform Mann-Whitney U test
                statistic, p_value = mannwhitneyu(
                    region_data['Perivascular'], 
                    region_data['Distant'],
                    alternative='two-sided'
                )
                p_values.append(p_value)
                valid_comparisons.append((i, pathway))
                
            except Exception as e:
                continue  # Skip if statistical test fails
    
    # Apply Benjamini-Hochberg correction
    if len(p_values) > 0:
        rejected, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
        
        # Second pass: add annotations with corrected p-values
        for idx, (i, pathway) in enumerate(valid_comparisons):
            corrected_p = corrected_p_values[idx]
            
            # Get y position for annotation
            y_max = plot_df[plot_df['pathway'] == pathway.replace('_', ' ').title()]['score'].max()
            y_annotation = y_max + (y_max * 0.1)
            
            # Significance levels based on corrected p-values
            if corrected_p < 0.001:
                sig_text = '***'
            elif corrected_p < 0.01:
                sig_text = '**'
            elif corrected_p < 0.05:
                sig_text = '*'
            else:
                sig_text = 'ns'
            
            # Add significance text
            ax.text(i, y_annotation, sig_text, ha='center', va='bottom',
                   fontweight='bold', fontsize=10)

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
    Enhanced visualization of cell type-pathway correlations
    """
    if not correlation_results:
        return
    
    print(f"  Creating enhanced cell-pathway correlation visualizations...")
    
    # Prepare data
    all_pathways = set()
    all_cell_types = list(correlation_results.keys())
    
    for cell_correlations in correlation_results.values():
        all_pathways.update(cell_correlations.keys())
    
    all_pathways = sorted(list(all_pathways))
    
    # Build correlation matrix for panel 1
    correlation_matrix = []
    for cell_type in all_cell_types:
        cell_row = []
        for pathway in all_pathways:
            if pathway in correlation_results[cell_type]:
                cell_row.append(correlation_results[cell_type][pathway]['correlation'])
            else:
                cell_row.append(0)
        correlation_matrix.append(cell_row)
    
    # Create enhanced visualization
    fig, axes = plt.subplots(2, 2, figsize=(24, 20))
    
    # Panel 1: Enhanced correlation heatmap
    correlation_df = pd.DataFrame(correlation_matrix, 
                                 index=[ct.replace('TC_', '').replace('_', ' ') for ct in all_cell_types],
                                 columns=[pw.replace('_', ' ') for pw in all_pathways])
    
    # Create mask for non-significant correlations (optional)
    significance_matrix = []
    for cell_type in all_cell_types:
        cell_row = []
        for pathway in all_pathways:
            if pathway in correlation_results[cell_type]:
                p_val = correlation_results[cell_type][pathway]['p_value']
                cell_row.append(p_val < 0.05)  # Significant if p < 0.05
            else:
                cell_row.append(False)
        significance_matrix.append(cell_row)
    
    significance_df = pd.DataFrame(significance_matrix, 
                                  index=correlation_df.index,
                                  columns=correlation_df.columns)
    
    sns.heatmap(correlation_df, ax=axes[0,0], cmap='RdBu_r', center=0, 
               annot=True, fmt='.2f', cbar_kws={'label': 'Correlation Coefficient'},
               mask=~significance_df, annot_kws={'size': 8})
    axes[0,0].set_title(f'Cell Type-Pathway Correlations\n{slide_name}\n(Only significant correlations shown, p<0.05)')
    axes[0,0].set_xlabel('Pathways')
    axes[0,0].set_ylabel('Cell Types')
    
    # Panel 2: Bubble plot - much more informative!
    bubble_data = []
    for cell_type, pathway_corrs in correlation_results.items():
        for pathway, corr_data in pathway_corrs.items():
            bubble_data.append({
                'cell_type': cell_type.replace('TC_', '').replace('_', ' '),
                'pathway': pathway.replace('_', ' '),
                'correlation': corr_data['correlation'],
                'abs_correlation': corr_data['abs_correlation'],
                'p_value': corr_data['p_value'],
                'neg_log_p': -np.log10(corr_data['p_value']),
                'is_significant': corr_data['p_value'] < 0.05
            })
    
    bubble_df = pd.DataFrame(bubble_data)
    
    # Create bubble plot
    scatter = axes[0,1].scatter(
        range(len(bubble_df)), 
        bubble_df['neg_log_p'],
        s=bubble_df['abs_correlation'] * 400,  # Size represents correlation strength
        c=bubble_df['correlation'],
        cmap='RdBu_r',
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5
    )
    
    # Add significance threshold line
    axes[0,1].axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
    axes[0,1].axhline(y=-np.log10(0.01), color='darkred', linestyle='--', alpha=0.7, label='p=0.01')
    
    # Annotate top correlations
    top_correlations = bubble_df.nlargest(8, 'abs_correlation')
    for idx, row in top_correlations.iterrows():
        if row['is_significant']:
            axes[0,1].annotate(f"{row['cell_type']}\n{row['pathway'][:15]}...", 
                              (bubble_df.index[bubble_df.index == idx][0], row['neg_log_p']),
                              xytext=(5, 5), textcoords='offset points', 
                              fontsize=7, alpha=0.8,
                              bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.3))
    
    axes[0,1].set_xlabel('Correlation Index')
    axes[0,1].set_ylabel('-log10(p-value)')
    axes[0,1].set_title('Correlation Significance vs Strength\n(Bubble size = |correlation|, Color = correlation direction)')
    axes[0,1].legend()
    
    # Add colorbar for correlation direction
    cbar1 = plt.colorbar(scatter, ax=axes[0,1])
    cbar1.set_label('Correlation Coefficient')
    
    # Panel 3: Volcano plot style - correlation vs significance
    significant = bubble_df['is_significant']
    
    # Plot non-significant points
    axes[1,0].scatter(bubble_df.loc[~significant, 'correlation'], 
                     bubble_df.loc[~significant, 'neg_log_p'],
                     c='lightgray', alpha=0.6, s=50, label='Non-significant')
    
    # Plot significant positive correlations
    pos_sig = significant & (bubble_df['correlation'] > 0)
    axes[1,0].scatter(bubble_df.loc[pos_sig, 'correlation'], 
                     bubble_df.loc[pos_sig, 'neg_log_p'],
                     c='red', alpha=0.8, s=80, label='Positive (p<0.05)')
    
    # Plot significant negative correlations
    neg_sig = significant & (bubble_df['correlation'] < 0)
    axes[1,0].scatter(bubble_df.loc[neg_sig, 'correlation'], 
                     bubble_df.loc[neg_sig, 'neg_log_p'],
                     c='blue', alpha=0.8, s=80, label='Negative (p<0.05)')
    
    # Add threshold lines
    axes[1,0].axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
    axes[1,0].axvline(x=0, color='black', linestyle='-', alpha=0.3)
    axes[1,0].axvline(x=0.3, color='orange', linestyle='--', alpha=0.5, label='|r|=0.3')
    axes[1,0].axvline(x=-0.3, color='orange', linestyle='--', alpha=0.5)
    
    # Annotate strongest significant correlations
    strong_sig = bubble_df[(bubble_df['is_significant']) & (bubble_df['abs_correlation'] > 0.4)]
    for _, row in strong_sig.iterrows():
        axes[1,0].annotate(f"{row['cell_type']}\n{row['pathway'][:15]}", 
                          (row['correlation'], row['neg_log_p']),
                          xytext=(5, 5), textcoords='offset points', 
                          fontsize=7, alpha=0.9,
                          bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow', alpha=0.7))
    
    axes[1,0].set_xlabel('Correlation Coefficient')
    axes[1,0].set_ylabel('-log10(p-value)')
    axes[1,0].set_title('Volcano Plot: Correlation Strength vs Significance')
    axes[1,0].legend()
    
    # Panel 4: Top correlations bar plot
    top_n = 15
    top_correlations_plot = bubble_df.nlargest(top_n, 'abs_correlation')
    
    # Create labels
    labels = [f"{row['cell_type']}\n{row['pathway'][:15]}" for _, row in top_correlations_plot.iterrows()]
    
    # Color bars by correlation direction and significance
    colors = []
    for _, row in top_correlations_plot.iterrows():
        if not row['is_significant']:
            colors.append('lightgray')
        elif row['correlation'] > 0:
            colors.append('red')
        else:
            colors.append('blue')
    
    bars = axes[1,1].barh(range(len(top_correlations_plot)), 
                         top_correlations_plot['abs_correlation'],
                         color=colors, alpha=0.7, edgecolor='black')
    
    # Add correlation values as text
    for i, (_, row) in enumerate(top_correlations_plot.iterrows()):
        axes[1,1].text(row['abs_correlation'] + 0.01, i, 
                      f"{row['correlation']:.3f}\n(p={row['p_value']:.1e})",
                      va='center', fontsize=8)
    
    axes[1,1].set_yticks(range(len(top_correlations_plot)))
    axes[1,1].set_yticklabels(labels, fontsize=9)
    axes[1,1].set_xlabel('Absolute Correlation')
    axes[1,1].set_title(f'Top {top_n} Strongest Correlations\n(Red=Positive, Blue=Negative, Gray=Non-sig)')
    axes[1,1].grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"enhanced_cell_pathway_correlations_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()