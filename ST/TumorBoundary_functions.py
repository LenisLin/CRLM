import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, deque
from gseapy import Msigdb
import gseapy as gp

import scipy
from scipy import stats
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, pearsonr
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.preprocessing import StandardScaler

import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# Set modern styling
plt.style.use('default')
sns.set_palette("husl")

# Define analysis parameters based on Nature Communications paper
BOUNDARY_ANALYSIS_PARAMS = {
    'boundary_width': 300,  # ±150μm from interface (Nature paper)
    'liver_distant_threshold': -150,  # μm
    'tumor_center_threshold': 150,   # μm
    'invasion_detection_threshold': 2,  # cells for budding detection
    'immune_exclusion_distance': 100,  # μm
    'vessel_boundary_radius': 50,  # μm for vessel-boundary analysis
}

# Define pathway analysis parameters specific to boundary
BOUNDARY_PATHWAY_PARAMS = {
    'position_correlation_threshold': 0.3,  # Minimum correlation with boundary position
    'boundary_zone_width': 150,  # μm for pathway enrichment analysis
    'min_spots_for_analysis': 10,  # Minimum spots for reliable analysis
    'pathway_score_method': 'mean',  # or 'ssgsea'
    'enrichment_pvalue_threshold': 0.05
}

# Define key cell types and pathways
BOUNDARY_CELL_TYPES = {
    'immune_cells': ['CD8_Teff', 'CD8_Tem', 'CD4_CXCL13', 'CD4_Treg', 'DC_LAMP3', 'Mac_M1', 'Mac_M2'],
    'stromal_cells': ['CAF_CXCL14', 'CAF_normal', 'Endothelial', 'Fibroblasts'],
    'tumor_cells': ['TC_Glycolysis', 'TC_EMT', 'TC_Proliferation', 'TC_Quiescent'],
    'liver_cells': ['Hepatocytes', 'Cholangiocytes', 'Stellate_cells']
}


def find_connected_components(interface_indices, coords, max_distance=1.5):
    """
    Find connected components of interface spots using graph connectivity
    Much simpler and more accurate than DBSCAN clustering
    
    Parameters:
    - interface_indices: list of interface spot indices
    - coords: coordinate array
    - max_distance: maximum distance to consider spots as connected
    
    Returns:
    - List of connected components (each component is a list of indices)
    """
    if len(interface_indices) == 0:
        return []
    
    # Build adjacency graph
    interface_coords = coords[interface_indices]
    distance_matrix = cdist(interface_coords, interface_coords)
    
    # Create adjacency list
    adjacency = defaultdict(set)
    for i in range(len(interface_indices)):
        for j in range(i + 1, len(interface_indices)):
            if distance_matrix[i, j] <= max_distance:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Find connected components using DFS
    visited = set()
    components = []
    
    for start_idx in range(len(interface_indices)):
        if start_idx not in visited:
            # DFS to find all connected nodes
            component = []
            stack = deque([start_idx])
            
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    component.append(interface_indices[current])  # Store original index
                    
                    # Add unvisited neighbors to stack
                    for neighbor in adjacency[current]:
                        if neighbor not in visited:
                            stack.append(neighbor)
            
            if component:
                components.append(component)
    
    return components

def define_tumor_boundary(adata, slide_name, 
                                          tumor_threshold=1.5, 
                                          adjacency_distance=1.5,
                                          min_interface_component_size=20,
                                          connectivity_distance=1.5,
                                          vessel_exclusion=True):
    """
    Enhanced tumor boundary definition with graph-based connectivity refinement
    
    Key improvements:
    1. Data-driven thresholds for tumor/liver classification
    2. Spatial adjacency for initial interface detection  
    3. Graph connectivity to ensure interface continuity
    4. Component size filtering to remove artifacts
    5. Vessel exclusion to prevent intratumoral vessels from being misclassified
    
    Parameters:
    - tumor_threshold: Threshold for tumor signature (default 1.5)
    - adjacency_distance: Max distance for adjacency detection (default 1.5)
    - min_interface_component_size: Minimum spots in interface component (default 10)
    - connectivity_distance: Max distance for interface connectivity (default 1.5)
    - vessel_exclusion: Whether to exclude vessel spots from liver classification
    """
    print(f"\n--- Enhanced Tumor Boundary Detection with Graph Connectivity ({slide_name}) ---")
    
    # Get cell type abundances
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    vessel_types = [col for col in adata.obs.columns if col.startswith('EC_')]  # Endothelial cells
    
    if len(tumor_types) == 0:
        print(f"  No tumor cell types found")
        return False
    
    # Calculate signatures
    adata.obs['Tumor_Signature'] = adata.obs[tumor_types].sum(axis=1)
    
    # Vessel signature for exclusion (if requested)
    if vessel_exclusion and len(vessel_types) > 0:
        adata.obs['Vessel_Signature'] = adata.obs[vessel_types].sum(axis=1)
        vessel_threshold = np.percentile(adata.obs['Vessel_Signature'], 90)  # Top 10% as vessels
        vessel_mask = adata.obs['Vessel_Signature'] >= vessel_threshold
        print(f"  Vessel exclusion enabled: {vessel_mask.sum()} vessel spots identified")
    else:
        vessel_mask = np.zeros(len(adata.obs), dtype=bool)
    
    print(f"  Thresholds: Tumor={tumor_threshold:.2f}, Adjacency={adjacency_distance:.1f}")
    
    # STEP 1: Define core regions (exclude vessels from liver classification)
    tumor_core_mask = adata.obs['Tumor_Signature'] >= tumor_threshold
    
    # Liver is everything that's not tumor and not vessel (if vessel_exclusion is True)
    if vessel_exclusion:
        liver_core_mask = (~tumor_core_mask) & (~vessel_mask)
    else:
        liver_core_mask = (~tumor_core_mask)
    
    print(f"  Core regions: Tumor={tumor_core_mask.sum()}, Liver={liver_core_mask.sum()}")
    if vessel_exclusion:
        print(f"  Excluded vessels: {vessel_mask.sum()}")
    
    if tumor_core_mask.sum() == 0 or liver_core_mask.sum() == 0:
        print(f"  Insufficient core regions")
        return False
    
    # STEP 2: Initial adjacency-based interface detection
    coords = adata.obs[['array_row', 'array_col']].values
    
    tumor_coords = coords[tumor_core_mask]
    liver_coords = coords[liver_core_mask]
    
    # Calculate minimum distances
    tumor_to_liver_dist = cdist(tumor_coords, liver_coords)
    liver_to_tumor_dist = cdist(liver_coords, tumor_coords)
    
    # Find adjacent spots
    tumor_adjacent = np.any(tumor_to_liver_dist <= adjacency_distance, axis=1)
    liver_adjacent = np.any(liver_to_tumor_dist <= adjacency_distance, axis=1)
    
    # Initial interface classification
    tissue_type = np.full(len(adata.obs), 'Other', dtype=object)
    tissue_type[tumor_core_mask] = 'Tumor_Core'
    tissue_type[liver_core_mask] = 'Normal_Liver'
    if vessel_exclusion:
        tissue_type[vessel_mask] = 'Vessel'
    
    # Collect initial interface indices
    tumor_indices = np.where(tumor_core_mask)[0]
    liver_indices = np.where(liver_core_mask)[0]
    
    initial_interface_indices = []
    initial_interface_indices.extend(tumor_indices[tumor_adjacent])
    initial_interface_indices.extend(liver_indices[liver_adjacent])
    
    print(f"  Initial interface spots: {len(initial_interface_indices)}")
    
    # STEP 3: Graph connectivity analysis for refinement
    if len(initial_interface_indices) >= min_interface_component_size:
        
        # Find connected components using graph connectivity
        components = find_connected_components(
            initial_interface_indices, coords, connectivity_distance
        )
        
        print(f"  Graph connectivity results:")
        print(f"    Connected components found: {len(components)}")
        
        # Filter components by size and analyze
        refined_interface_indices = set()
        component_info = []
        
        for i, component in enumerate(components):
            component_size = len(component)
            component_info.append((i, component_size))
            
            # Keep component if it meets size threshold
            if component_size >= min_interface_component_size:
                refined_interface_indices.update(component)
                print(f"    Component {i}: {component_size} spots (KEPT)")
            else:
                print(f"    Component {i}: {component_size} spots (REMOVED - too small)")
        
        # Update tissue classification with refined interface
        for idx in initial_interface_indices:
            if idx in refined_interface_indices:
                tissue_type[idx] = 'Interface_Region'
            else:
                # Revert small components back to their original classification
                if tumor_core_mask[idx]:
                    tissue_type[idx] = 'Tumor_Core'
                elif liver_core_mask[idx]:
                    tissue_type[idx] = 'Normal_Liver'
                elif vessel_exclusion and vessel_mask[idx]:
                    tissue_type[idx] = 'Vessel'
        
        print(f"  Refined interface spots: {len(refined_interface_indices)}")
        print(f"  Removed artifacts: {len(initial_interface_indices) - len(refined_interface_indices)}")
        
        # Store component information
        adata.obs['Interface_Component_ID'] = -1  # Default for non-interface
        for i, component in enumerate(components):
            if len(component) >= min_interface_component_size:
                for idx in component:
                    adata.obs.iloc[idx, adata.obs.columns.get_loc('Interface_Component_ID')] = i
        
    else:
        print(f"  Too few initial interface spots for connectivity analysis")
        return False
    
    # STEP 4: Final boundary distance calculations
    adata.obs['Tissue_Type'] = tissue_type
    interface_mask = adata.obs['Tissue_Type'] == 'Interface_Region'
    
    if interface_mask.sum() > 0:
        interface_coords = coords[interface_mask]
        distance_matrix = cdist(coords, interface_coords)
        boundary_distances = distance_matrix.min(axis=1)
        adata.obs['Boundary_Distance_um'] = boundary_distances * 100  # 100μm spacing
        
        # Boundary positions (negative for liver, positive for tumor)
        boundary_positions = np.zeros(len(adata.obs))
        tumor_final = adata.obs['Tissue_Type'] == 'Tumor_Core'
        liver_final = adata.obs['Tissue_Type'] == 'Normal_Liver'
        
        boundary_positions[tumor_final] = adata.obs.loc[tumor_final, 'Boundary_Distance_um']
        boundary_positions[liver_final] = -adata.obs.loc[liver_final, 'Boundary_Distance_um']
        boundary_positions[interface_mask] = 0
        
        adata.obs['Boundary_Position'] = boundary_positions
        
        print(f"  ✅ Enhanced boundary detection completed:")
        print(f"    Tumor core: {tumor_final.sum()} spots")
        print(f"    Normal liver: {liver_final.sum()} spots")
        print(f"    Interface: {interface_mask.sum()} spots")
        if vessel_exclusion:
            print(f"    Vessels: {(adata.obs['Tissue_Type'] == 'Vessel').sum()} spots")
        print(f"    Interface efficiency: {interface_mask.sum()/len(initial_interface_indices)*100:.1f}%")
        
        return True
    
    return False

def create_boundary_zones(adata, slide_name):
    """
    Step 2: Create spatial zones based on Nature Communications methodology
    """
    print(f"  Creating boundary zones...")
    
    # Define zones based on boundary position
    params = BOUNDARY_ANALYSIS_PARAMS
    
    conditions = [
        adata.obs['Boundary_Position'] < params['liver_distant_threshold'],
        (adata.obs['Boundary_Position'] >= params['liver_distant_threshold']) & 
        (adata.obs['Boundary_Position'] <= params['tumor_center_threshold']),
        adata.obs['Boundary_Position'] > params['tumor_center_threshold']
    ]
    
    choices = ['Distant_Liver', 'Boundary_Region', 'Tumor_Center']
    adata.obs['Boundary_Zone'] = np.select(conditions, choices, default='Unknown')
    
    # Report zone statistics
    zone_counts = adata.obs['Boundary_Zone'].value_counts()
    total_spots = len(adata.obs)
    
    print(f"  Spatial zones created:")
    for zone in ['Distant_Liver', 'Boundary_Region', 'Tumor_Center']:
        if zone in zone_counts:
            count = zone_counts[zone]
            percentage = count / total_spots * 100
            print(f"    {zone}: {count} spots ({percentage:.1f}%)")
    
    return True

def characterize_boundary_composition(adata, slide_name):
    """
    Step 3: Characterize cell composition across boundary zones
    """
    print(f"  Characterizing boundary cell composition...")
    
    # Get available cell types
    available_cell_types = {}
    for category, cell_list in BOUNDARY_CELL_TYPES.items():
        available = [col for col in adata.obs.columns 
                    if any(cell_type in col for cell_type in cell_list)]
        if available:
            available_cell_types[category] = available
    
    # Calculate cell densities by zone
    boundary_composition = {}
    
    for zone in ['Distant_Liver', 'Boundary_Region', 'Tumor_Center']:
        zone_mask = adata.obs['Boundary_Zone'] == zone
        
        if zone_mask.sum() > 0:
            zone_composition = {}
            
            for category, cell_types in available_cell_types.items():
                for cell_type in cell_types:
                    if cell_type in adata.obs.columns:
                        zone_composition[cell_type] = adata.obs.loc[zone_mask, cell_type].mean()
            
            boundary_composition[zone] = zone_composition
    
    return boundary_composition

#%% PHASE 2: FUNCTIONAL BOUNDARY ANALYSIS

def analyze_invasion_patterns(adata, slide_name):
    """
    Step 4: Analyze tumor invasion patterns at boundary
    """
    print(f"\n--- Step 4: Analyzing invasion patterns ---")
    
    # Focus on boundary region
    boundary_mask = adata.obs['Boundary_Zone'] == 'Boundary_Region'
    
    if boundary_mask.sum() == 0:
        return {}
    
    boundary_data = adata.obs[boundary_mask]
    
    # Detect tumor budding (isolated tumor cells in liver-adjacent areas)
    tumor_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    
    invasion_results = {}
    
    # Calculate invasion metrics
    for tc_type in tumor_types:
        tumor_abundance = boundary_data[tc_type]
        
        # Invasion intensity (tumor cells in boundary region)
        invasion_intensity = tumor_abundance.mean()
        
        # Invasion penetration (how far tumor cells extend)
        liver_side_mask = boundary_data['Boundary_Position'] < 0
        if liver_side_mask.sum() > 0:
            liver_invasion = boundary_data.loc[liver_side_mask, tc_type].mean()
        else:
            liver_invasion = 0
        
        # Invasion gradient (change across boundary)
        positions = boundary_data['Boundary_Position'].values
        abundances = boundary_data[tc_type].values
        
        if len(positions) > 5:
            gradient_corr, _ = pearsonr(positions, abundances)
        else:
            gradient_corr = 0
        
        invasion_results[tc_type] = {
            'invasion_intensity': invasion_intensity,
            'liver_penetration': liver_invasion,
            'invasion_gradient': gradient_corr,
            'boundary_abundance': invasion_intensity
        }
        
        print(f"    {tc_type}: intensity={invasion_intensity:.3f}, penetration={liver_invasion:.3f}")
    
    # Calculate EMT gradient across boundary
    if 'EMT_Score' in adata.obs.columns:
        emt_gradient = calculate_emt_boundary_gradient(adata, boundary_mask)
        invasion_results['EMT_gradient'] = emt_gradient
    
    return invasion_results

def calculate_emt_boundary_gradient(adata, boundary_mask):
    """Calculate EMT score gradient across boundary"""
    boundary_data = adata.obs[boundary_mask]
    
    if len(boundary_data) > 5:
        positions = boundary_data['Boundary_Position'].values
        emt_scores = boundary_data['EMT_Score'].values
        
        gradient_corr, p_value = pearsonr(positions, emt_scores)
        
        return {
            'emt_gradient_correlation': gradient_corr,
            'p_value': p_value,
            'mean_emt_boundary': emt_scores.mean()
        }
    
    return {'emt_gradient_correlation': 0, 'p_value': 1, 'mean_emt_boundary': 0}

def analyze_microenvironment_permissiveness(adata, slide_name):
    """
    Step 5: Analyze microenvironment that permits/restricts invasion
    """
    print(f"\n--- Step 5: Analyzing microenvironment permissiveness ---")
    
    boundary_mask = adata.obs['Boundary_Zone'] == 'Boundary_Region'
    
    if boundary_mask.sum() == 0:
        return {}
    
    boundary_data = adata.obs[boundary_mask]
    
    microenv_results = {}
    
    # 1. Immune suppression analysis
    immune_metrics = calculate_immune_suppression_metrics(boundary_data)
    microenv_results['immune_suppression'] = immune_metrics
    
    # 2. CAF activation analysis  
    caf_metrics = calculate_caf_activation_metrics(boundary_data)
    microenv_results['caf_activation'] = caf_metrics
    
    # 3. ECM organization analysis
    ecm_metrics = calculate_ecm_organization_metrics(boundary_data)
    microenv_results['ecm_organization'] = ecm_metrics
    
    return microenv_results

def calculate_immune_suppression_metrics(boundary_data):
    """Calculate immune suppression metrics in boundary"""
    immune_metrics = {}
    
    # CD8/Treg ratio (higher = less suppressed)
    cd8_cols = [col for col in boundary_data.columns if col.startswith('CD8')]
    treg_cols = [col for col in boundary_data.columns if col.startswith('Treg')]
    
    if cd8_cols and treg_cols:
        cd8_total = boundary_data[cd8_cols].sum(axis=1).mean()
        treg_total = boundary_data[treg_cols].sum(axis=1).mean()
        
        immune_metrics['cd8_treg_ratio'] = cd8_total / (treg_total + 1e-6)
        immune_metrics['cd8_abundance'] = cd8_total
        immune_metrics['treg_abundance'] = treg_total
    
    return immune_metrics

def calculate_caf_activation_metrics(boundary_data):
    """Calculate CAF activation metrics"""
    caf_metrics = {}
    
    # CAFs (from Nature paper - associated with immune exclusion)
    cxcl_cols = [col for col in boundary_data.columns if col.startswith('CAF_')]
    normal_fib_cols = [col for col in boundary_data.columns if col.startswith('Fibro_')]
    
    if cxcl_cols:
        caf_metrics['caf_abundance'] = boundary_data[cxcl_cols].sum(axis=1).mean()
    
    if normal_fib_cols:
        caf_metrics['normal_fibroblast_abundance'] = boundary_data[normal_fib_cols].sum(axis=1).mean()
    
    # CAF activation ratio
    if cxcl_cols and normal_fib_cols:
        caf_total = boundary_data[cxcl_cols].sum(axis=1).mean()
        normal_total = boundary_data[normal_fib_cols].sum(axis=1).mean()
        caf_metrics['caf_activation_ratio'] = caf_total / (normal_total + 1e-6)
    
    return caf_metrics

def calculate_ecm_organization_metrics(boundary_data):
    """Calculate ECM organization metrics"""
    ecm_metrics = {}
    
    # Look for ECM-related genes/signatures
    ecm_related = [col for col in boundary_data.columns 
                   if any(term in col.upper() for term in ['COL', 'MATRIX', 'ECM', 'FIBRO'])]
    
    if ecm_related:
        ecm_metrics['ecm_signature'] = boundary_data[ecm_related].sum(axis=1).mean()
        ecm_metrics['ecm_organization_score'] = np.std(boundary_data[ecm_related].values)
    
    return ecm_metrics

def integrate_boundary_vessel_analysis(adata, slide_name):
    """
    Step 6: Integrate boundary analysis with existing vessel analysis
    """
    print(f"\n--- Step 6: Integrating boundary-vessel analysis ---")
    
    if 'Is_Vascular_Hotspot' not in adata.obs.columns:
        print("  No vessel analysis found - run vessel analysis first")
        return {}
    
    boundary_mask = adata.obs['Boundary_Zone'] == 'Boundary_Region'
    
    # Vessel density at boundary
    boundary_vessel_density = adata.obs.loc[boundary_mask, 'Vascular_Sum'].mean() if 'Vascular_Sum' in adata.obs.columns else 0
    
    # Vessel hotspots in boundary
    boundary_hotspots = (adata.obs['Boundary_Zone'] == 'Boundary_Region') & adata.obs['Is_Vascular_Hotspot']
    boundary_hotspot_rate = boundary_hotspots.sum() / boundary_mask.sum() * 100 if boundary_mask.sum() > 0 else 0
    
    # GLUT1+ cells at boundary
    glut1_boundary_results = {}
    if 'TC_Glycolysis' in adata.obs.columns:
        glut1_boundary = adata.obs.loc[boundary_mask, 'TC_Glycolysis'].mean()
        glut1_boundary_results['glut1_boundary_abundance'] = glut1_boundary
        
        # GLUT1+ vessel association at boundary
        boundary_vessel_mask = boundary_mask & adata.obs['Is_Vascular_Hotspot']
        if boundary_vessel_mask.sum() > 0:
            glut1_vessel_boundary = adata.obs.loc[boundary_vessel_mask, 'TC_Glycolysis'].mean()
            glut1_boundary_results['glut1_vessel_cooccurrence'] = glut1_vessel_boundary
    
    boundary_vessel_results = {
        'vessel_density_boundary': boundary_vessel_density,
        'vessel_hotspot_rate': boundary_hotspot_rate,
        'glut1_boundary': glut1_boundary_results
    }
    
    return boundary_vessel_results

def analyze_treatment_effects_boundary(adata, slide_name, boundary_metrics):
    """
    Step 7: Analyze treatment effects on boundary characteristics
    """
    treatment = adata.obs['treatment_status'].iloc[0]
    
    treatment_effects = {
        'slide_name': slide_name,
        'treatment': treatment,
        'boundary_metrics': boundary_metrics
    }
    
    return treatment_effects

#%% PHASE 4: CROSS-SLIDE ANALYSIS AND VISUALIZATION

def visualize_boundary_analysis(adata, slide_name, boundary_composition, invasion_results, slide_dir):
    """Comprehensive boundary analysis visualization"""
    
    treatment = adata.obs['treatment_status'].iloc[0]
    
    # Create comprehensive figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. Boundary zones spatial map
    if 'Boundary_Zone' in adata.obs.columns:
        sc.pl.spatial(adata, color='Boundary_Zone',
                     palette={'Distant_Liver': 'lightblue', 'Boundary_Region': 'orange', 'Tumor_Center': 'red'},
                     ax=axes[0,0], size=1.2, img_key='hires', show=False,
                     title='Boundary Zones')
    
    # 2. Tissue type classification
    if 'Tissue_Type' in adata.obs.columns:
        sc.pl.spatial(adata, color='Tissue_Type',
                     palette={'Normal_Liver': 'green', 'Interface_Region': 'yellow', 'Tumor_Core': 'red', 'Vessel': 'blue', 'Other': 'grey'},
                     ax=axes[0,1], size=1.2, img_key='hires', show=False,
                     title='Tissue Classification')
    
    # 3. Boundary distance map
    if 'Boundary_Position' in adata.obs.columns:
        sc.pl.spatial(adata, color='Boundary_Position', cmap='RdBu_r',
                     ax=axes[0,2], size=1.2, img_key='hires', show=False,
                     title='Distance from Boundary (μm)')
    
    # 4. Cell composition across zones
    if boundary_composition:
        plot_boundary_composition(boundary_composition, axes[1,0])
    
    # 5. Invasion patterns
    if invasion_results:
        plot_invasion_patterns(invasion_results, axes[1,1])
    
    # 6. Risk assessment summary
    # This will be populated by risk scoring
    axes[1,2].text(0.1, 0.5, f'Boundary Analysis\n{slide_name}\n({treatment})', 
                   transform=axes[1,2].transAxes, fontsize=12, fontweight='bold')
    axes[1,2].axis('off')
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"boundary_analysis_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()

def plot_boundary_composition(boundary_composition, ax):
    """Plot cell composition across boundary zones"""
    
    # Prepare data for heatmap
    zones = list(boundary_composition.keys())
    all_cell_types = set()
    for zone_data in boundary_composition.values():
        all_cell_types.update(zone_data.keys())
    
    # Create matrix
    composition_matrix = []
    cell_type_names = []
    
    for cell_type in all_cell_types:
        if any('TC_' in cell_type or 'CD' in cell_type or 'CAF' in cell_type for cell_type in [cell_type]):
            row = [boundary_composition[zone].get(cell_type, 0) for zone in zones]
            composition_matrix.append(row)
            cell_type_names.append(cell_type.replace('TC_', '').replace('CD', 'CD'))
    
    if composition_matrix:
        sns.heatmap(composition_matrix, xticklabels=zones, yticklabels=cell_type_names,
                   annot=True, fmt='.2f', cmap='viridis', ax=ax)
        ax.set_title('Cell Composition Across Boundary')

def plot_invasion_patterns(invasion_results, ax):
    """Plot invasion patterns summary"""
    
    tumor_types = []
    invasion_scores = []
    
    for tc_type, results in invasion_results.items():
        if isinstance(results, dict) and 'invasion_intensity' in results:
            tumor_types.append(tc_type.replace('TC_', ''))
            invasion_scores.append(results['invasion_intensity'])
    
    if tumor_types:
        colors = ['red' if 'Glycolysis' in tc else 'blue' for tc in tumor_types]
        bars = ax.bar(tumor_types, invasion_scores, color=colors, alpha=0.7)
        ax.set_title('Invasion Intensity by Tumor Type')
        ax.set_ylabel('Invasion Score')
        ax.tick_params(axis='x', rotation=45)

def compare_boundary_across_treatments(all_boundary_results, figure_path):
    """
    Cross-slide analysis: Compare boundary characteristics across treatments
    """
    print(f"\n--- CROSS-SLIDE BOUNDARY COMPARISON ---")
    
    # Convert to DataFrame
    boundary_df = pd.DataFrame(all_boundary_results)
    
    if len(boundary_df) == 0:
        print("No boundary results to compare")
        return
    
    # Extract key metrics for comparison
    comparison_metrics = []
    
    for result in all_boundary_results:
        metrics = {'slide_name': result['slide_name'], 'treatment': result['treatment']}
        
        # Add boundary metrics
        if 'boundary_metrics' in result:
            boundary_data = result['boundary_metrics']
            
            # Risk scores
            if 'risk_score' in boundary_data:
                metrics.update(boundary_data['risk_score'])
            
            # Invasion metrics
            if 'invasion_results' in boundary_data:
                invasion_data = boundary_data['invasion_results']
                for tc_type, inv_results in invasion_data.items():
                    if isinstance(inv_results, dict):
                        for metric, value in inv_results.items():
                            metrics[f'{tc_type}_{metric}'] = value
        
        comparison_metrics.append(metrics)
    
    comparison_df = pd.DataFrame(comparison_metrics)
    
    # Statistical comparisons
    if len(comparison_df) > 0:
        create_boundary_treatment_comparison_plots(comparison_df, figure_path)
    
    return comparison_df

def create_boundary_treatment_comparison_plots(comparison_df, figure_path):
    """
    Create a comprehensive, publication-ready visualization of boundary analysis
    """
    
    # Modern color palette
    colors = {
        'untreated': '#3498db',      # Professional blue
        'pre_chemotherapy': '#e74c3c',  # Professional red
        'treated': '#2ecc71'         # Professional green (if needed)
    }
    
    # Create figure with custom layout
    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(4, 4, figure=fig, hspace=0.3, wspace=0.25)
    
    # =============================================================================
    # 1. MAIN HEATMAP: All tumor cell invasion metrics
    # =============================================================================
    ax_heatmap = fig.add_subplot(gs[0:2, 0:2])
    
    # Prepare heatmap data
    tumor_metrics = ['invasion_intensity', 'liver_penetration', 'invasion_gradient', 'boundary_abundance']
    tumor_types = ['TC_EMT', 'TC_Glycolysis', 'TC_LipidMeta', 'TC_Proliferation', 'TC_Quiescent']
    
    heatmap_data = []
    heatmap_labels = []
    
    for tumor_type in tumor_types:
        for metric in tumor_metrics:
            col_name = f"{tumor_type}_{metric}"
            if col_name in comparison_df.columns:
                # Calculate fold change (treated/untreated)
                untreated_mean = comparison_df[comparison_df['treatment'] == 'untreated'][col_name].mean()
                treated_mean = comparison_df[comparison_df['treatment'] == 'pre_chemotherapy'][col_name].mean()
                
                if untreated_mean > 0:
                    fold_change = treated_mean / untreated_mean
                    heatmap_data.append(fold_change)
                else:
                    heatmap_data.append(treated_mean)
                
                # Clean label
                clean_tumor = tumor_type.replace('TC_', '')
                clean_metric = metric.replace('_', ' ').title()
                heatmap_labels.append(f"{clean_tumor}\n{clean_metric}")
    
    # Create heatmap matrix
    heatmap_matrix = np.array(heatmap_data).reshape(len(tumor_types), len(tumor_metrics))
    
    # Plot heatmap
    im = ax_heatmap.imshow(heatmap_matrix, cmap='RdBu_r', aspect='auto', vmin=0.5, vmax=2.0)
    
    # Customize heatmap
    ax_heatmap.set_xticks(range(len(tumor_metrics)))
    ax_heatmap.set_xticklabels([m.replace('_', ' ').title() for m in tumor_metrics], fontsize=11)
    ax_heatmap.set_yticks(range(len(tumor_types)))
    ax_heatmap.set_yticklabels([t.replace('TC_', '') for t in tumor_types], fontsize=11)
    ax_heatmap.set_title('Treatment Effect on Invasion Metrics\n(Fold Change: Pre-chemo/Untreated)', 
                        fontsize=14, fontweight='bold', pad=20)
    
    # Add values to heatmap
    for i in range(len(tumor_types)):
        for j in range(len(tumor_metrics)):
            text = ax_heatmap.text(j, i, f'{heatmap_matrix[i, j]:.2f}', 
                                 ha="center", va="center", fontweight='bold',
                                 color='white' if abs(heatmap_matrix[i, j] - 1) > 0.3 else 'black')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax_heatmap, shrink=0.8)
    cbar.set_label('Fold Change', fontsize=12)
    
    # =============================================================================
    # 2. INVASION INTENSITY COMPARISON (Ridge Plot Style)
    # =============================================================================
    ax_ridge = fig.add_subplot(gs[0, 2:])
    
    invasion_cols = [col for col in comparison_df.columns if 'invasion_intensity' in col]
    
    y_positions = np.arange(len(invasion_cols))
    colors_ridge = plt.cm.Set3(np.linspace(0, 1, len(invasion_cols)))
    
    for i, col in enumerate(invasion_cols):
        tumor_type = col.split('_')[1]
        
        # Plot distributions for each treatment
        for j, treatment in enumerate(['untreated', 'pre_chemotherapy']):
            data = comparison_df[comparison_df['treatment'] == treatment][col]
            if len(data) > 0:
                # Create violin plot manually for ridge effect
                parts = ax_ridge.violinplot([data], positions=[i + j*0.4], widths=0.35, 
                                          showmeans=True, showmedians=True)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[treatment])
                    pc.set_alpha(0.7)
    
    ax_ridge.set_yticks(y_positions + 0.2)
    ax_ridge.set_yticklabels([col.split('_')[1] for col in invasion_cols], fontsize=11)
    ax_ridge.set_xlabel('Invasion Intensity', fontsize=12)
    ax_ridge.set_title('Invasion Intensity by Tumor Type', fontsize=14, fontweight='bold')
    
    # Add legend
    legend_elements = [mpatches.Patch(color=colors['untreated'], label='Untreated'),
                      mpatches.Patch(color=colors['pre_chemotherapy'], label='Pre-chemotherapy')]
    ax_ridge.legend(handles=legend_elements, loc='upper right')
    
    # =============================================================================
    # 3. EMT GRADIENT ANALYSIS
    # =============================================================================
    ax_emt = fig.add_subplot(gs[1, 2:])
    
    # Scatter plot with EMT metrics
    if 'EMT_gradient_emt_gradient_correlation' in comparison_df.columns:
        scatter = ax_emt.scatter(comparison_df['EMT_gradient_emt_gradient_correlation'],
                               comparison_df['EMT_gradient_mean_emt_boundary'],
                               c=[colors[t] for t in comparison_df['treatment']],
                               s=100, alpha=0.7, edgecolors='white', linewidth=2)
        
        # Add sample labels
        for idx, row in comparison_df.iterrows():
            ax_emt.annotate(row['slide_name'].split('_')[0], 
                          (row['EMT_gradient_emt_gradient_correlation'], 
                           row['EMT_gradient_mean_emt_boundary']),
                          xytext=(5, 5), textcoords='offset points', fontsize=9)
        
        ax_emt.set_xlabel('EMT Gradient Correlation', fontsize=12)
        ax_emt.set_ylabel('Mean EMT at Boundary', fontsize=12)
        ax_emt.set_title('EMT Spatial Gradient Analysis', fontsize=14, fontweight='bold')
        
        # Add trend line
        z = np.polyfit(comparison_df['EMT_gradient_emt_gradient_correlation'], 
                      comparison_df['EMT_gradient_mean_emt_boundary'], 1)
        p = np.poly1d(z)
        ax_emt.plot(comparison_df['EMT_gradient_emt_gradient_correlation'], 
                   p(comparison_df['EMT_gradient_emt_gradient_correlation']), 
                   "r--", alpha=0.5)
    
    # =============================================================================
    # 4. LIVER PENETRATION ANALYSIS
    # =============================================================================
    ax_liver = fig.add_subplot(gs[2, 0:2])
    
    penetration_cols = [col for col in comparison_df.columns if 'liver_penetration' in col]
    
    # Create grouped bar plot
    x_pos = np.arange(len(penetration_cols))
    width = 0.35
    
    untreated_means = []
    treated_means = []
    
    for col in penetration_cols:
        untreated_means.append(comparison_df[comparison_df['treatment'] == 'untreated'][col].mean())
        treated_means.append(comparison_df[comparison_df['treatment'] == 'pre_chemotherapy'][col].mean())
    
    bars1 = ax_liver.bar(x_pos - width/2, untreated_means, width, 
                        label='Untreated', color=colors['untreated'], alpha=0.8)
    bars2 = ax_liver.bar(x_pos + width/2, treated_means, width, 
                        label='Pre-chemotherapy', color=colors['pre_chemotherapy'], alpha=0.8)
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax_liver.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                     f'{height:.3f}', ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        ax_liver.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                     f'{height:.3f}', ha='center', va='bottom', fontsize=9)
    
    ax_liver.set_xlabel('Tumor Cell Type', fontsize=12)
    ax_liver.set_ylabel('Liver Penetration Score', fontsize=12)
    ax_liver.set_title('Liver Penetration by Treatment and Cell Type', fontsize=14, fontweight='bold')
    ax_liver.set_xticks(x_pos)
    ax_liver.set_xticklabels([col.split('_')[1] for col in penetration_cols], rotation=45)
    ax_liver.legend()
    ax_liver.grid(True, alpha=0.3)
    
    # =============================================================================
    # 5. BOUNDARY ABUNDANCE ANALYSIS
    # =============================================================================
    ax_boundary = fig.add_subplot(gs[2, 2:])
    
    abundance_cols = [col for col in comparison_df.columns if 'boundary_abundance' in col]
    
    # Violin plot with swarm overlay
    abundance_data = []
    abundance_labels = []
    abundance_treatments = []
    
    for col in abundance_cols:
        tumor_type = col.split('_')[1]
        for treatment in comparison_df['treatment'].unique():
            data = comparison_df[comparison_df['treatment'] == treatment][col].values
            abundance_data.extend(data)
            abundance_labels.extend([tumor_type] * len(data))
            abundance_treatments.extend([treatment] * len(data))
    
    plot_df = pd.DataFrame({
        'abundance': abundance_data,
        'tumor_type': abundance_labels,
        'treatment': abundance_treatments
    })
    
    sns.violinplot(data=plot_df, x='tumor_type', y='abundance', hue='treatment', 
                   palette=[colors['untreated'], colors['pre_chemotherapy']], 
                   ax=ax_boundary, inner='box')
    
    ax_boundary.set_xlabel('Tumor Cell Type', fontsize=12)
    ax_boundary.set_ylabel('Boundary Abundance', fontsize=12)
    ax_boundary.set_title('Boundary Abundance Distribution', fontsize=14, fontweight='bold')
    ax_boundary.tick_params(axis='x', rotation=45)
    
    # =============================================================================
    # FINAL STYLING
    # =============================================================================
    
    # Main title
    fig.suptitle('Comprehensive Boundary Analysis: Treatment Effects on Tumor Invasion Patterns', 
                fontsize=18, fontweight='bold', y=0.98)
    
    # Add subtle background
    fig.patch.set_facecolor('white')
    
    # Save high-quality figure
    plt.savefig(f"{figure_path}/enhanced_boundary_analysis.pdf", 
               bbox_inches='tight', dpi=300, facecolor='white')
    
    plt.show()
    
    return fig

def load_boundary_pathway_databases():
    """
    Load pathway databases for boundary enrichment analysis
    """
    print(f"\n--- Loading boundary pathway databases ---")
    
    try:
        # Load MSigDB Hallmark pathways
        msig = Msigdb()
        
        print(f"  Fetching Hallmark gene sets...")
        hallmark = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
        
        # Define boundary-relevant pathway keywords
        boundary_keywords = ['emt', 'epithelial', 'migration', 'invasion', 'inflammation', 'immune',
                           'matrix', 'adhesion', 'apoptosis', 'proliferation', 'interferon', 'complement']
        
        # Filter relevant pathways
        boundary_pathways = {}
        for pathway_name, genes in hallmark.items():
            pathway_lower = pathway_name.lower()
            if any(keyword in pathway_lower for keyword in boundary_keywords):
                clean_name = pathway_name.replace('HALLMARK_', '')
                boundary_pathways[clean_name] = genes
        
        print(f"  ✅ Loaded {len(boundary_pathways)} boundary-relevant Hallmark pathways")
        
        # # Add custom boundary-related pathways
        # custom_pathways = {
        #     'INVASION_CORE': ['MMP2', 'MMP9', 'MMP14', 'TIMP1', 'TIMP2', 'ITGA1', 'ITGA2', 'ITGB1', 'ITGB3'],
        #     'EMT_CORE': ['VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'ZEB2', 'TWIST1', 'CDH1', 'CDH2', 'FN1', 'FOXC2'],
        #     'IMMUNE_EXCLUSION': ['CXCL14', 'TGFB1', 'TGFB2', 'IL10', 'FOXP3', 'CTLA4', 'PDCD1', 'CD274'],
        #     'CAF_ACTIVATION': ['FAP', 'ACTA2', 'COL1A1', 'COL1A2', 'PDGFRA', 'PDGFRB', 'CXCL12'],
        #     'HEPATOCYTE_RESPONSE': ['ALB', 'AFP', 'CRP', 'FGB', 'APOB', 'CYP1A2', 'CYP2E1'],
        #     'LIVER_INJURY': ['HMGB1', 'S100A8', 'S100A9', 'IL1B', 'TNF', 'IL6', 'CXCL1', 'CXCL2']
        # }
        
        # boundary_pathways.update(custom_pathways)
        # print(f"  ✅ Added {len(custom_pathways)} custom boundary pathways")
        
        return boundary_pathways
        
    except Exception as e:
        print(f"  ❌ Error loading pathways: {e}")
        print(f"  Using backup pathway definitions...")
        
        # Backup pathway definitions
        backup_pathways = {
            'EMT': ['VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'ZEB2', 'TWIST1', 'CDH1', 'CDH2', 'FN1', 'FOXC2'],
            'INVASION': ['MMP2', 'MMP9', 'MMP14', 'TIMP1', 'TIMP2', 'ITGA1', 'ITGA2', 'ITGB1', 'ITGB3'],
            'IMMUNE_SUPPRESSION': ['FOXP3', 'CTLA4', 'PDCD1', 'CD274', 'TGFB1', 'IL10', 'CXCL14'],
            'INFLAMMATION': ['TNF', 'IL1B', 'IL6', 'CXCL1', 'CXCL2', 'CCL2', 'NFKB1'],
            'PROLIFERATION': ['MKI67', 'PCNA', 'CCND1', 'CCNE1', 'CDK1', 'CDK2', 'CDK4']
        }
        
        return backup_pathways
    

#%% STEP 1: IDENTIFY BOUNDARY-ASSOCIATED GENES

def identify_boundary_proximity_genes(adata, slide_name):
    """
    Step 1: Identify genes associated with boundary position
    Returns genes correlated with position relative to tumor-liver boundary
    """
    print(f"\n--- Step 1: Identifying boundary-proximity genes ({slide_name}) ---")
    
    if 'Boundary_Position' not in adata.obs.columns:
        print("  No boundary position data found - run boundary analysis first")
        return None, None
    
    # Focus on boundary region only (±150μm from interface)
    boundary_mask = adata.obs['Boundary_Zone'] == 'Boundary_Region'
    
    if boundary_mask.sum() < BOUNDARY_PATHWAY_PARAMS['min_spots_for_analysis']:
        print(f"  Insufficient boundary spots for analysis: {boundary_mask.sum()}")
        return None, None
    
    # Get gene expression data
    boundary_data = adata[boundary_mask].copy()
    boundary_positions = boundary_data.obs['Boundary_Position'].values
    
    # Calculate gene-position correlations
    gene_correlations = {}
    boundary_associated_genes = {'invasion_associated': [], 'liver_associated': []}
    
    print(f"  Analyzing {boundary_data.n_vars} genes across {boundary_data.n_obs} boundary spots...")
    
    if hasattr(boundary_data.X, 'toarray'):
        gene_expr_all = boundary_data.X.toarray()
    else:
        gene_expr_all = boundary_data.X
    
    for gene_idx, gene_name in enumerate(boundary_data.var.index):
        if gene_idx % 1000 == 0:
            print(f"    Processed {gene_idx}/{boundary_data.n_vars} genes")
        
        # Get gene expression
        gene_expr = gene_expr_all[:, gene_idx]
        
        # Skip genes with low variance
        if np.var(gene_expr) < 1e-6:
            continue
        
        # Calculate correlation with boundary position
        if len(gene_expr) == len(boundary_positions):
            corr, p_value = pearsonr(boundary_positions, gene_expr)
            
            if abs(corr) >= BOUNDARY_PATHWAY_PARAMS['position_correlation_threshold'] and p_value < 0.05:
                gene_correlations[gene_name] = {
                    'correlation': corr,
                    'p_value': p_value,
                    'abs_correlation': abs(corr)
                }
                
                # Classify genes
                if corr > BOUNDARY_PATHWAY_PARAMS['position_correlation_threshold']:
                    boundary_associated_genes['invasion_associated'].append(gene_name)  # Higher in tumor direction
                elif corr < -BOUNDARY_PATHWAY_PARAMS['position_correlation_threshold']:
                    boundary_associated_genes['liver_associated'].append(gene_name)  # Higher in liver direction
    
    print(f"  Found {len(gene_correlations)} boundary-associated genes:")
    print(f"    Invasion-associated (high tumor side): {len(boundary_associated_genes['invasion_associated'])}")
    print(f"    Liver-associated (high liver side): {len(boundary_associated_genes['liver_associated'])}")
    
    return gene_correlations, boundary_associated_genes

def visualize_boundary_gene_heatmap(adata, gene_correlations, boundary_associated_genes, slide_dir, slide_name):
    """
    Visualize boundary-associated genes with heatmap
    """
    if not gene_correlations:
        return
    
    print(f"  Creating boundary-gene association heatmap...")
    
    # Get top genes for visualization
    sorted_genes = sorted(gene_correlations.items(), key=lambda x: x[1]['abs_correlation'], reverse=True)
    top_genes = [gene for gene, data in sorted_genes[:50]]  # Top 50 genes
    
    if len(top_genes) == 0:
        return
    
    # Focus on boundary regions
    boundary_mask = adata.obs['Boundary_Zone'] == 'Boundary_Region'
    boundary_data = adata[boundary_mask].copy()
    
    # Create spatial zones based on boundary position
    positions = boundary_data.obs['Boundary_Position'].values
    
    # Define zones
    liver_side_mask = positions < -50  # Liver side
    interface_mask = (positions >= -50) & (positions <= 50)  # Interface
    tumor_side_mask = positions > 50  # Tumor side
    
    zones = ['Liver Side', 'Interface', 'Tumor Side']
    zone_masks = [liver_side_mask, interface_mask, tumor_side_mask]
    
    # Calculate mean expression in each zone
    expression_matrix = []
    gene_names_clean = []
    
    for gene in top_genes:
        if gene in boundary_data.var.index:
            gene_idx = boundary_data.var.index.get_loc(gene)
            
            if hasattr(boundary_data.X, 'toarray'):
                gene_expr = boundary_data.X[:, gene_idx].toarray().flatten()
            else:
                gene_expr = boundary_data.X[:, gene_idx]
            
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
    
    sns.heatmap(expression_df_norm, ax=axes[0], cmap='viridis', 
                cbar_kws={'label': 'Normalized Expression'})
    axes[0].set_title(f'Boundary-Associated Gene Expression\n{slide_name}\n(Top {len(top_genes)} genes)', 
                     fontweight='bold')
    axes[0].set_xlabel('Boundary Zone')
    axes[0].set_ylabel('Genes')
    
    # 2. Correlation with boundary position
    correlations = [gene_correlations[gene]['correlation'] for gene in top_genes if gene in gene_correlations]
    p_values = [gene_correlations[gene]['p_value'] for gene in top_genes if gene in gene_correlations]
    
    if len(correlations) > 0:
        # Color code by significance and direction
        colors = []
        for corr, p in zip(correlations, p_values):
            if p < 0.001:
                color = 'red' if corr > 0 else 'blue'
            elif p < 0.01:
                color = 'orange' if corr > 0 else 'cyan'
            elif p < 0.05:
                color = 'yellow' if corr > 0 else 'lightblue'
            else:
                color = 'gray'
            colors.append(color)
        
        bars = axes[1].barh(range(len(correlations)), correlations, color=colors, alpha=0.7)
        axes[1].set_yticks(range(len(gene_names_clean)))
        axes[1].set_yticklabels(gene_names_clean)
        axes[1].set_xlabel('Correlation with Boundary Position\n(+: Tumor side, -: Liver side)')
        axes[1].set_title('Gene-Boundary Position Correlations\n(Red/Orange/Yellow: Tumor-associated, Blue/Cyan/LightBlue: Liver-associated)')
        axes[1].axvline(x=0, color='black', linestyle='-', alpha=0.5)
        
        # Add correlation value labels
        for i, (bar, corr) in enumerate(zip(bars, correlations)):
            axes[1].text(corr + (0.01 if corr > 0 else -0.01), bar.get_y() + bar.get_height()/2, 
                        f'{corr:.2f}', ha='left' if corr > 0 else 'right', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"boundary_gene_heatmap_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()

def calculate_boundary_pathway_scores(adata, boundary_pathways, method='mean'):
    """
    Calculate pathway enrichment scores for each spot in boundary region
    """
    print(f"  Calculating boundary pathway scores using {method} method...")
    
    pathway_scores = {}
    available_genes = adata.var.index.tolist()

    if hasattr(adata.X, 'toarray'):
        pathway_expr_all = adata.X.toarray()
    else:
        pathway_expr_all = adata.X
    
    for pathway_name, pathway_genes in boundary_pathways.items():
        # Find available genes in the pathway
        available_pathway_genes = [gene for gene in pathway_genes if gene in available_genes]
        
        if len(available_pathway_genes) >= 3:  # Minimum genes for reliable score
            if method == 'mean':
                # Simple mean expression method
                pathway_indices = [adata.var.index.get_loc(gene) for gene in available_pathway_genes]
                
                pathway_expr = pathway_expr_all[:, pathway_indices]
                pathway_score = np.mean(pathway_expr, axis=1)
                pathway_scores[pathway_name] = pathway_score
                
                print(f"    {pathway_name}: {len(available_pathway_genes)}/{len(pathway_genes)} genes available")
    
    return pathway_scores

def perform_boundary_pathway_enrichment(adata, gene_correlations, boundary_pathways, slide_name):
    """
    Step 2: Perform pathway enrichment analysis on boundary-associated genes
    """
    print(f"\n--- Step 2: Boundary pathway enrichment analysis ({slide_name}) ---")
    
    if not gene_correlations:
        return None
    
    # Get invasion-associated and liver-associated genes
    invasion_genes = [gene for gene, data in gene_correlations.items() 
                     if data['correlation'] > BOUNDARY_PATHWAY_PARAMS['position_correlation_threshold']]
    liver_genes = [gene for gene, data in gene_correlations.items() 
                  if data['correlation'] < -BOUNDARY_PATHWAY_PARAMS['position_correlation_threshold']]
    
    print(f"  Invasion-associated genes: {len(invasion_genes)}")
    print(f"  Liver-associated genes: {len(liver_genes)}")
    
    # Calculate pathway scores for all spots
    pathway_scores = calculate_boundary_pathway_scores(adata, boundary_pathways)
    
    # Store pathway scores in adata
    for pathway_name, scores in pathway_scores.items():
        adata.obs[f'boundary_pathway_{pathway_name}'] = scores
    
    # Perform enrichment analysis
    enrichment_results = {}
    
    for gene_set_name, genes in [('invasion_associated', invasion_genes), 
                                ('liver_associated', liver_genes)]:
        
        if len(genes) < 3:
            continue
        
        set_enrichments = {}
        
        for pathway_name, pathway_genes in boundary_pathways.items():
            # Calculate overlap
            overlap_genes = set(genes) & set(pathway_genes)
            overlap_count = len(overlap_genes)
            
            if overlap_count >= 2:  # Minimum overlap
                # Simple enrichment calculation
                total_genes_in_pathway = len(pathway_genes)
                total_boundary_genes = len(genes)
                enrichment_ratio = (overlap_count / total_boundary_genes) / (total_genes_in_pathway / 20000)  # Assume ~20k total genes
                
                set_enrichments[pathway_name] = {
                    'overlap_count': overlap_count,
                    'overlap_genes': list(overlap_genes),
                    'enrichment_ratio': enrichment_ratio,
                    'pathway_size': total_genes_in_pathway
                }
        
        enrichment_results[gene_set_name] = set_enrichments
        
        print(f"  {gene_set_name} enrichment: {len(set_enrichments)} pathways enriched")
    
    return enrichment_results, pathway_scores

def visualize_boundary_pathway_enrichment(enrichment_results, pathway_scores, adata, slide_dir, slide_name):
    """
    Visualize boundary pathway enrichment results
    """
    if not enrichment_results:
        return
    
    print(f"  Creating boundary pathway enrichment visualizations...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Invasion-associated gene enrichment
    if 'invasion_associated' in enrichment_results:
        invasion_enriched = enrichment_results['invasion_associated']
        
        if len(invasion_enriched) > 0:
            pathways = list(invasion_enriched.keys())
            enrichments = [invasion_enriched[p]['enrichment_ratio'] for p in pathways]
            overlaps = [invasion_enriched[p]['overlap_count'] for p in pathways]
            
            # Sort by enrichment ratio
            sorted_data = sorted(zip(pathways, enrichments, overlaps), key=lambda x: x[1], reverse=True)
            pathways_sorted, enrichments_sorted, overlaps_sorted = zip(*sorted_data)
            
            # Color by overlap count
            colors = plt.cm.Reds([o/max(overlaps_sorted) for o in overlaps_sorted])
            
            bars = axes[0,0].barh(range(len(pathways_sorted)), enrichments_sorted, color=colors)
            axes[0,0].set_yticks(range(len(pathways_sorted)))
            axes[0,0].set_yticklabels([p.replace('_', ' ') for p in pathways_sorted])
            axes[0,0].set_xlabel('Enrichment Ratio')
            axes[0,0].set_title('Invasion-Associated Gene Enrichment\n(Tumor Side of Boundary)')
            axes[0,0].axvline(x=1, color='black', linestyle='--', alpha=0.5)
    
    # 2. Liver-associated gene enrichment
    if 'liver_associated' in enrichment_results:
        liver_enriched = enrichment_results['liver_associated']
        
        if len(liver_enriched) > 0:
            pathways = list(liver_enriched.keys())
            enrichments = [liver_enriched[p]['enrichment_ratio'] for p in pathways]
            overlaps = [liver_enriched[p]['overlap_count'] for p in pathways]
            
            # Sort by enrichment ratio
            sorted_data = sorted(zip(pathways, enrichments, overlaps), key=lambda x: x[1], reverse=True)
            pathways_sorted, enrichments_sorted, overlaps_sorted = zip(*sorted_data)
            
            # Color by overlap count
            colors = plt.cm.Blues([o/max(overlaps_sorted) for o in overlaps_sorted])
            
            bars = axes[0,1].barh(range(len(pathways_sorted)), enrichments_sorted, color=colors)
            axes[0,1].set_yticks(range(len(pathways_sorted)))
            axes[0,1].set_yticklabels([p.replace('_', ' ') for p in pathways_sorted])
            axes[0,1].set_xlabel('Enrichment Ratio')
            axes[0,1].set_title('Liver-Associated Gene Enrichment\n(Liver Side of Boundary)')
            axes[0,1].axvline(x=1, color='black', linestyle='--', alpha=0.5)
    
    # 3. Pathway score distributions across boundary position
    if pathway_scores:
        # Get top enriched pathways for visualization
        top_pathways = list(pathway_scores.keys())
        
        # Create violin plot of pathway scores by boundary position
        position_data = []
        
        for pathway in top_pathways:
            pathway_col = f'boundary_pathway_{pathway}'
            if pathway_col in adata.obs.columns:
                boundary_mask = adata.obs['Boundary_Zone'] == 'Boundary_Region'
                boundary_data = adata.obs[boundary_mask]
                
                positions = boundary_data['Boundary_Position'].values
                scores = boundary_data[pathway_col].values
                
                # Define position bins
                for position, score in zip(positions, scores):
                    if position < -50:
                        pos_bin = 'Liver Side'
                    elif position <= 50:
                        pos_bin = 'Interface'
                    else:
                        pos_bin = 'Tumor Side'
                    
                    position_data.append({
                        'pathway': pathway.replace('_', ' '),
                        'position_bin': pos_bin,
                        'score': score
                    })
        
        if position_data:
            position_df = pd.DataFrame(position_data)
            
            sns.boxplot(data=position_df, x='position_bin', y='score', hue='pathway', ax=axes[1,0])
            axes[1,0].set_title('Pathway Scores by Boundary Position')
            axes[1,0].set_xlabel('Boundary Position')
            axes[1,0].set_ylabel('Pathway Score')
            axes[1,0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 4. Summary enrichment heatmap
    if enrichment_results:
        # Create combined enrichment matrix
        all_pathways = set()
        for gene_set_enrichments in enrichment_results.values():
            all_pathways.update(gene_set_enrichments.keys())
        
        enrichment_matrix = []
        pathway_names = []
        
        for pathway in all_pathways:
            enrichment_row = []
            for gene_set in ['invasion_associated', 'liver_associated']:
                if gene_set in enrichment_results and pathway in enrichment_results[gene_set]:
                    enrichment_row.append(enrichment_results[gene_set][pathway]['enrichment_ratio'])
                else:
                    enrichment_row.append(0)
            
            enrichment_matrix.append(enrichment_row)
            pathway_names.append(pathway.replace('_', ' '))
        
        if enrichment_matrix:
            enrichment_df = pd.DataFrame(enrichment_matrix, 
                                       columns=['Invasion-Associated', 'Liver-Associated'], 
                                       index=pathway_names)
            
            sns.heatmap(enrichment_df, ax=axes[1,1], cmap='RdBu_r', center=1, 
                       annot=True, fmt='.2f', cbar_kws={'label': 'Enrichment Ratio'})
            axes[1,1].set_title('Boundary Pathway Enrichment Summary')
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"boundary_pathway_enrichment_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()

def analyze_boundary_cell_pathway_correlations(adata, pathway_scores, slide_name):
    """
    Step 3: Analyze correlation between cell types and pathway scores in boundary region
    """
    print(f"\n--- Step 3: Boundary cell type-pathway correlations ({slide_name}) ---")
    
    if not pathway_scores:
        return None
    
    # Focus on boundary region
    boundary_mask = adata.obs['Boundary_Zone'] == 'Boundary_Region'
    
    if boundary_mask.sum() < BOUNDARY_PATHWAY_PARAMS['min_spots_for_analysis']:
        print(f"  Insufficient boundary spots: {boundary_mask.sum()}")
        return None
    
    analysis_data = adata.obs[boundary_mask]
    
    # Get cell types of interest
    tumor_cell_types = [col for col in adata.obs.columns if col.startswith('TC_')]
    immune_cell_types = [col for col in adata.obs.columns 
                        if any(cell_type in col for cell_type in ['CD8', 'CD4', 'DC_', 'Mac_', 'Treg'])]
    stromal_cell_types = [col for col in adata.obs.columns 
                         if any(cell_type in col for cell_type in ['CAF_', 'Fibroblast', 'Endothelial'])]
    
    all_cell_types = tumor_cell_types + immune_cell_types + stromal_cell_types
    
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
                pathway_col = f'boundary_pathway_{pathway_name}'
                
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
                
                # Print top correlations for tumor cells
                if cell_type.startswith('TC_'):
                    top_pathways = sorted(cell_correlations.items(), 
                                        key=lambda x: x[1]['abs_correlation'], reverse=True)[:3]
                    
                    print(f"  {cell_type} top correlations:")
                    for pathway, data in top_pathways:
                        corr = data['correlation']
                        p_val = data['p_value']
                        print(f"    {pathway}: r={corr:+.3f}, p={p_val:.3f}")
    
    return correlation_results

def visualize_boundary_cell_pathway_correlations(adata, correlation_results, slide_dir, slide_name):
    """
    Visualize boundary cell type-pathway correlations
    """
    if not correlation_results:
        return
    
    print(f"  Creating boundary cell-pathway correlation visualizations...")
    
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
    axes[0].set_title(f'Cell Type-Pathway Correlations\n{slide_name}\n(Boundary Region)')
    axes[0].set_xlabel('Pathways')
    axes[0].set_ylabel('Cell Types')
    
    # 2. Focus on tumor cells
    tumor_correlations = {ct: corrs for ct, corrs in correlation_results.items() 
                         if ct.startswith('TC_')}
    
    if tumor_correlations:
        tumor_data = []
        
        for cell_type, pathway_corrs in tumor_correlations.items():
            for pathway, corr_data in pathway_corrs.items():
                tumor_data.append({
                    'cell_type': cell_type.replace('TC_', ''),
                    'pathway': pathway.replace('_', ' '),
                    'correlation': corr_data['correlation'],
                    'abs_correlation': corr_data['abs_correlation'],
                    'p_value': corr_data['p_value']
                })
        
        tumor_df = pd.DataFrame(tumor_data)
        
        # Create scatter plot
        scatter = axes[1].scatter(tumor_df['correlation'], tumor_df['abs_correlation'],
                                c=tumor_df['p_value'], cmap='viridis_r', 
                                s=100, alpha=0.7, edgecolors='black')
        
        # Add text annotations for significant correlations
        for _, row in tumor_df.iterrows():
            if row['p_value'] < 0.01 and row['abs_correlation'] > 0.4:
                axes[1].annotate(f"{row['cell_type']}\n{row['pathway']}", 
                               (row['correlation'], row['abs_correlation']),
                               xytext=(5, 5), textcoords='offset points', 
                               fontsize=8, alpha=0.8)
        
        axes[1].set_xlabel('Correlation Coefficient')
        axes[1].set_ylabel('Absolute Correlation')
        axes[1].set_title('Tumor Cell-Pathway Correlations\n(Color: p-value)')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=axes[1])
        cbar.set_label('p-value')
        
        # Add significance thresholds
        axes[1].axhline(y=0.3, color='red', linestyle='--', alpha=0.5, label='|r| = 0.3')
        axes[1].axvline(x=0, color='black', linestyle='-', alpha=0.3)
        axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(slide_dir, f"boundary_cell_pathway_correlations_{slide_name}.pdf"), 
               bbox_inches='tight', dpi=300)
    plt.show()