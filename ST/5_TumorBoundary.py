#!/usr/bin/env python3
"""
Comprehensive Tumor Boundary Analysis for CRLM Post-Surgery Treatment
Combines vessel analysis with boundary characterization for recurrence prediction
"""

import os
import pickle
import scanpy as sc
import pickle

from Vascular_functions import *
from TumorBoundary_functions import *

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%% TUMOR BOUNDARY ANALYSIS: CRLM RECURRENCE PREDICTION

# Define paths
run_name = '/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/analysis/cell2location_map'
figure_path = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/boundary_analysis"

if not os.path.exists(figure_path):
    os.makedirs(figure_path, exist_ok=True)

# Define analysis parameters based on Nature Communications paper
BOUNDARY_ANALYSIS_PARAMS = {
    'boundary_width': 300,  # ±150μm from interface (Nature paper)
    'liver_distant_threshold': -150,  # μm
    'tumor_center_threshold': 150,   # μm
    'invasion_detection_threshold': 2,  # cells for budding detection
    'immune_exclusion_distance': 100,  # μm
    'vessel_boundary_radius': 50,  # μm for vessel-boundary analysis
}

# Define key cell types and pathways
BOUNDARY_CELL_TYPES = {
    'immune_cells': ['CD8_Teff', 'CD8_Tem', 'CD4_CXCL13', 'CD4_Treg', 'DC_LAMP3', 'Mac_M1', 'Mac_M2'],
    'stromal_cells': ['CAF_CXCL14', 'CAF_normal', 'Endothelial', 'Fibroblasts'],
    'tumor_cells': ['TC_Glycolysis', 'TC_EMT', 'TC_Proliferation', 'TC_Quiescent'],
    'liver_cells': ['Hepatocytes', 'Cholangiocytes', 'Stellate_cells']
}

INVASION_PATHWAYS = ['EMT', 'Matrix_remodeling', 'Cell_migration', 'Angiogenesis']
IMMUNE_PATHWAYS = ['T_cell_activation', 'Immune_suppression', 'Chemotaxis', 'Inflammation']
METABOLISM_PATHWAYS = ['Glycolysis', 'Oxidative_phosphorylation', 'Lipid_metabolism', 'Hypoxia']
    
# Define slides (using your existing slide definition logic)
all_slides = os.listdir(run_name)
relevant_slides = []
treatment_mapping = {}

for slide in all_slides:
    if slide.startswith('FDZS'):
        continue
    elif 'Untreated' in slide:
        relevant_slides.append(slide)
        treatment_mapping[slide] = 'untreated'
    elif 'NACPR' in slide:
        relevant_slides.append(slide)
        treatment_mapping[slide] = 'pre_chemotherapy'
    elif slide.startswith('GSE225857'):
        relevant_slides.append(slide)
        treatment_mapping[slide] = 'untreated'

print("TUMOR BOUNDARY ANALYSIS FOR CRLM RECURRENCE PREDICTION")
print("="*60)
print(f"Analyzing {len(relevant_slides)} slides for boundary characteristics")

# Single slide analysis
all_boundary_results = []

# Load pathway databases
boundary_pathways = load_boundary_pathway_databases()
    
for slide_name in relevant_slides:
    print(f"\n{'='*60}")
    print(f"BOUNDARY ANALYSIS: {slide_name}")
    print(f"{'='*60}")
    
    # Load slide data (using your existing function)
    adata = load_slide_data(slide_name, run_name, treatment_mapping)  

    # Calculate EMT score (using your existing function)  
    calculate_emt_score(adata, ['VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'ZEB2', 'TWIST1', 'CDH2', 'FN1'])
    
    # Create output directory
    slide_dir = create_slide_output_dir(figure_path, slide_name)
    
    # Phase 1: Boundary definition
    boundary_success = define_tumor_boundary(adata, slide_name, min_interface_component_size=20, vessel_exclusion=False)
    if not boundary_success:
        print(f"{slide_name} was inable to define boundary")
        continue;
    create_boundary_zones(adata, slide_name)
    boundary_composition = characterize_boundary_composition(adata, slide_name)
    
    # Phase 2: Functional analysis
    invasion_results = analyze_invasion_patterns(adata, slide_name)
    microenv_results = analyze_microenvironment_permissiveness(adata, slide_name)

    treatment_effects = analyze_treatment_effects_boundary(adata, slide_name, {
        'boundary_composition': boundary_composition,
        'invasion_results': invasion_results,
        'microenv_results': microenv_results    
    })
    
    # Visualization
    visualize_boundary_analysis(adata, slide_name, boundary_composition, 
                               invasion_results, slide_dir)
    
    print(f"✓ Completed boundary analysis for {slide_name}")

    all_boundary_results.append(treatment_effects)

    # Phase 3: Pathway analysis
    # Step 1: Identify boundary-associated genes
    gene_correlations, boundary_associated_genes = identify_boundary_proximity_genes(adata, slide_name)
    
    if gene_correlations:
        # Visualize boundary-gene associations
        visualize_boundary_gene_heatmap(adata, gene_correlations, boundary_associated_genes, slide_dir, slide_name)
        
        # Step 2: Pathway enrichment analysis
        enrichment_results, pathway_scores = perform_boundary_pathway_enrichment(
            adata, gene_correlations, boundary_pathways, slide_name)
        
        if enrichment_results:
            # Visualize pathway enrichment
            visualize_boundary_pathway_enrichment(enrichment_results, pathway_scores, adata, slide_dir, slide_name)
            
            # Step 3: Cell type-pathway correlations
            correlation_results = analyze_boundary_cell_pathway_correlations(adata, pathway_scores, slide_name)
            
            if correlation_results:
                visualize_boundary_cell_pathway_correlations(adata, correlation_results, slide_dir, slide_name)
        
        final_boundary_pathway_results = {
            'slide_name': slide_name,
            'treatment': adata.obs['treatment_status'].iloc[0],
            'gene_correlations': gene_correlations,
            'boundary_associated_genes': boundary_associated_genes,
            'enrichment_results': enrichment_results,
            'pathway_scores': pathway_scores,
            'correlation_results': correlation_results,
            'boundary_composition': boundary_composition
        }


# Cross-slide analysis
print(f"\n{'='*60}")
print("CROSS-SLIDE BOUNDARY COMPARISON")
print(f"{'='*60}")

comparison_results = compare_boundary_across_treatments(all_boundary_results, figure_path)

# Save results

final_boundary_results = {
    'individual_slides': all_boundary_results,
    'treatment_comparison': comparison_results.to_dict('records') if comparison_results is not None else [],
    'analysis_parameters': BOUNDARY_ANALYSIS_PARAMS,
    'methods_summary': {
        'boundary_definition': 'Interface-based ±150μm zones (Nature Communications method)',
        'invasion_analysis': 'Tumor budding and EMT gradient analysis',
        'microenvironment': 'Immune suppression and CAF activation metrics',
        'risk_scoring': 'Composite invasion + immune + vessel + CAF score',
        'clinical_relevance': 'Post-surgery recurrence prediction for CRLM'
    }
}

with open(os.path.join(figure_path, "boundary_analysis_results.pkl"), 'wb') as f:
    pickle.dump(final_boundary_results, f)

if comparison_results is not None:
    comparison_results.to_csv(os.path.join(figure_path, "boundary_comparison_summary.csv"), index=False)

print(f"\n✅ BOUNDARY ANALYSIS COMPLETE!")
print(f"📁 Results saved in: {figure_path}")
print(f"🎯 Clinical insights: Boundary signatures for CRLM recurrence prediction")
print(f"🚀 Ready for clinical validation and biomarker development!")