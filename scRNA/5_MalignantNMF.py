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

from malignant_functions import *

print("🧬 Starting Comprehensive Malignant Cells Analysis Pipeline")
print("="*70)

# Set paths
DataPath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","NMF_results")
result_path = figurePath

# create directories if they don't exist
os.makedirs(result_path, exist_ok=True)
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
sc.pl.umap(malignant_subset, color='sample_id', ax=axes[1,0], show=False,
          legend_loc=None, title='Patient Distribution')

# Quality metrics
sc.pl.umap(malignant_subset, color='nCount_RNA', ax=axes[0,1], show=False,
          title='Gene Counts', color_map='viridis')

sc.pl.umap(malignant_subset, color='percent.mt', ax=axes[1,1], show=False,
          title='MT Percentage', color_map='Reds')

plt.tight_layout()
plt.savefig(os.path.join(figurePath, "malignant_batch_effects.pdf"), dpi=300, bbox_inches='tight')
plt.show()

print("✅ Batch visualization saved")

#%% Task 3: Extract Metabolism-Associated Genes using MSigDB API
print("\n🧬 Task 3: Extracting Metabolism-Associated Genes from MSigDB")
print("="*60)

if 'malignant_cells' not in locals():
    malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))

# Get available gene sets from MSigDB
# metabolism_gene_sets = get_metabolism_pathways()
# stress_gene_sets = get_stress_resistance_pathways()
# oncogenic_gene_sets = get_oncogenic_pathways()
# metabolism_gene_sets.update(stress_gene_sets)
# metabolism_gene_sets.update(oncogenic_gene_sets)

metabolism_gene_sets = get_imc_aligned_hallmark_pathways()
coverage_report = validate_imc_pathway_coverage(metabolism_gene_sets, malignant_cells)

print(f"\nTotal gene sets collected: {len(metabolism_gene_sets)}")

# Create pathway coverage report
pathway_coverage_report = []
for pathway_name, genes in metabolism_gene_sets.items():
    available_genes = [g for g in genes if g in malignant_cells.var_names]
    coverage = len(available_genes) / len(genes) * 100 if genes else 0
    
    pathway_coverage_report.append({
        'pathway': pathway_name,
        'total_genes': len(genes),
        'available_genes': len(available_genes),
        'coverage_percent': coverage,
        'available_gene_list': available_genes
    })

# Convert to DataFrame for easy viewing
coverage_df = pd.DataFrame(pathway_coverage_report)
coverage_df = coverage_df.sort_values('coverage_percent', ascending=False)

print(f"\nTop 10 best covered pathways:")
print(coverage_df[['pathway', 'total_genes', 'available_genes', 'coverage_percent']].head(10))

# Filter pathways with good coverage (>30% and >5 genes)
well_covered_pathways = coverage_df[
    (coverage_df['coverage_percent'] > 30) & 
    (coverage_df['available_genes'] > 5)
]['pathway'].tolist()


# Extract all unique metabolism genes for NMF
all_metabolism_genes = set()
for pathway, genes in metabolism_gene_sets.items():
    all_metabolism_genes.update(genes)

all_metabolism_genes = list(all_metabolism_genes)
print(f"Total unique metabolism genes: {len(all_metabolism_genes)}")

# Check availability in malignant cells dataset
available_metabolism_genes = [gene for gene in all_metabolism_genes if gene in malignant_cells.var_names]
print(f"Available in dataset: {len(available_metabolism_genes)}")
print(f"Coverage: {len(available_metabolism_genes)/len(all_metabolism_genes)*100:.1f}%")
print(f"\n🎯 Ready for NMF with {len(available_metabolism_genes)} metabolism genes")

#%% Task 4: Robust Per-Patient NMF Analysis (Following CSCC Strategy)
print("\n🔬 Task 4: Robust Per-Patient Metabolism Program Discovery")
print("="*60)

## Use HVG genes for NMF    
# if 'malignant_cells' not in locals():
#     malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells.h5ad"))

# adata = malignant_cells.copy()
# if adata.X.max() > 50:  # Check if data looks like raw counts
#     sc.pp.normalize_total(adata, target_sum=1e4)
#     sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = 2000)
# available_metabolism_genes = adata.var_names[adata.var['highly_variable']].tolist()

# Run the robust NMF analysis
malignant_cells.obs["sample_id"] = malignant_cells.obs["batch"]
robust_results = run_robust_nmf_analysis(malignant_cells, available_metabolism_genes, result_path)

# Manually merge
# updated_results, updated_cells = merge_meta_programs(
#     malignant_cells=malignant_cells,
#     robust_results=robust_results,
#     source_mp=3,  # MP to merge (eliminate) # 7
#     target_mp=2,  # MP to merge into (expand) # 1
#     save_path=result_path,
#     update_scoring=True
# )

# Manually split
# updated_results, updated_cells = split_meta_program(
#     malignant_cells=malignant_cells,
#     robust_results=robust_results,
#     mp_to_split=4,
#     n_splits=2,
#     split_method='program_clustering',
#     save_path=result_path
# )

# Validate the merger
# validation_stats = validate_merger(updated_results, updated_cells, save_path=result_path)

# malignant_cells = updated_cells.copy()
# robust_results = updated_results.copy()

# Load the robust results if not already loaded
if not 'robust_results' in locals():
    with open(os.path.join(result_path, "robust_metabolism_nmf_results.pkl"), 'rb') as f:
        robust_results = pickle.load(f)

# Save robust NMF results
print("\n💾 Exporting Results for R Statistical Analysis")
print("="*60)

export_dir = export_results_for_r_analysis(malignant_cells, robust_results, result_path)
print(f"\n🎉 Ready for R analysis!")
print(f"📁 Data exported to: {export_dir}")

# Save final results
malignant_cells.write(os.path.join(result_path, "malignant_cells_with_robust_mps.h5ad"))
# malignant_cells = sc.read_h5ad(os.path.join(result_path, "malignant_cells_with_robust_mps.h5ad"))

print(f"\n🎉 Robust Meta-Program Analysis Completed!")
print(f"📁 Results saved to: {result_path}")

#%% Downstream analysis and visualizations
## Plot UMAP with Meta-Programs
print("\n📊 Plotting UMAP with Meta-Programs")
available_metabolism_genes = pd.read_csv(os.path.join(result_path, "NMF_results_for_R", "gene_names.csv"))["gene_name"].tolist()
create_mp_umap_visualization(malignant_cells, available_metabolism_genes, figurePath)

# Enhanced Clinical Analysis with Meta-Programs (move to R)
print("\n🔍 Enhanced Clinical Analysis with Meta-Programs")
print("="*70)

adata = malignant_cells.copy()

# Save all MP results
sample_mp_df = calculate_mp_fractions_per_sample(adata, mp_column='MP_assignment')
sample_mp_df.to_csv(f"{DataPath}/all_sample_mp_fraction_data.csv", index=False)

print(f"\n✅ Complete analysis pipeline finished!")

# Other figures
# Figure 1: NMF Robustness Analysis
create_nmf_robustness_plots(robust_results, figurePath)

# Figure 2: MP Assignment Quality
create_mp_assignment_quality_plots(malignant_cells, figurePath)

# Figure 3: Hierarchical Clustering Dendrogram
create_clustering_dendrogram(robust_results, figurePath)

# Figure 4: Create all MP visualization heatmaps
print("🎨 Creating comprehensive MP heatmap visualizations...")
print("="*60)

# Create comprehensive overview heatmap
fig1 = create_mp_clustering_heatmap(malignant_cells, robust_results, figurePath)

# Create focused signature comparison heatmap  
fig2 = create_mp_signature_comparison_heatmap(malignant_cells, robust_results, figurePath)

print("\n🎉 All MP heatmaps created successfully!")

# Execute the analysis
fig, enrichment_results = create_mp_pathway_analysis(malignant_cells, robust_results, figurePath)


#%% Included bulk transcript for analysis
# Load bulk transcriptomics data
bulk_datapath = "/mnt/NAS_21T/ProjectData/IMC_CRLM/bulkRNA"
dataIDs = ["BCGSC","E-TABM-1112","E-MTAB-1951"]

for dataID in dataIDs:
    bulk_exp = pd.read_csv(f"{bulk_datapath}/{dataID}/exp_data.csv", index_col=0)
    bulk_clinical = pd.read_csv(f"{bulk_datapath}/{dataID}/clinical_data.csv", index_col=0)

    # 1. Prepare the meta-program signatures dictionary
    mp_signatures = {f"MP_{key}": robust_results['meta_programs'][key]['signature_genes'] 
                    for key in robust_results['meta_programs']}

    # 2. Calculate enrichment scores for all MPs
    mp_enrichment_scores = calculate_mp_scores(bulk_exp, mp_signatures)

    # 3. Create the full output directory path
    output_base_dir = os.path.join(figurePath, "Bulk_analysis", dataID)

    # Ensure the directory exists
    os.makedirs(output_base_dir, exist_ok=True)
    print(f"Figures will be saved in: {os.path.abspath(output_base_dir)}")

    # 3. Loop through each MP and perform survival analysis
    if dataID == "BCGSC":
        survival_p_values = {}
        MPs = mp_enrichment_scores["MP"].unique()

        for mp_id in MPs:
            print("\n" + "="*40)
            print(f"🚀 Running Survival Analysis for {mp_id}")
            print("="*40)
            
            p_val = perform_survival_analysis(
                mp_scores=mp_enrichment_scores,
                clinical_data=bulk_clinical,
                mp_id=mp_id,
                time_col='OS_time',
                status_col='OS_status',
                save_path = output_base_dir
            )
            survival_p_values[mp_id] = p_val

        # 4. Display a summary of the results
        print("\n\n--- Survival Analysis Summary ---")
        summary_df = pd.DataFrame.from_dict(
            survival_p_values, 
            orient='index', 
            columns=['LogRank_PValue']
        )
        summary_df['Significant (p < 0.05)'] = summary_df['LogRank_PValue'] < 0.05
        print(summary_df.sort_values('LogRank_PValue'))
    
    else :
        # --- 2. Merge the two DataFrames ---
        bulk_clinical['Sample_ID'] = [str(x) for x in bulk_clinical.index] 
        mp_enrichment_scores['Sample_ID'] = mp_enrichment_scores['Sample_ID'].astype(str)
        merged_df = pd.merge(mp_enrichment_scores, bulk_clinical, left_on='Sample_ID', right_on='Sample_ID')

        unique_mps = merged_df['MP'].unique()

        for mp in unique_mps:
            # Filter data for the current metabolic pathway
            mp_df = merged_df[merged_df['MP'] == mp]

            # Check if there are at least two disease states to compare
            if dataID == "E-TABM-1112":
                disease_states = mp_df['DiseaseState'].unique()
            elif dataID == "E-MTAB-1951":
                disease_states = mp_df['risk_score_factor'].unique()

            if len(disease_states) < 2:
                print(f"Skipping '{mp}': Not enough disease states for comparison.")
                continue

            # Extract NES scores for each disease state
            # Assuming there are exactly two unique disease states for comparison (e.g., 'Control' and 'Diseased')
            # If you have more than two, you'll need to specify which pairs to compare.
            group1_name = disease_states[0]
            group2_name = disease_states[1]

            if dataID == "E-TABM-1112":
                group1_nes = mp_df[mp_df['DiseaseState'] == group1_name]['NES']
                group2_nes = mp_df[mp_df['DiseaseState'] == group2_name]['NES']
            elif dataID == "E-MTAB-1951":
                group1_nes = mp_df[mp_df['risk_score_factor'] == group1_name]['NES']
                group2_nes = mp_df[mp_df['risk_score_factor'] == group2_name]['NES']

            # Perform Wilcoxon rank-sum test (Mann-Whitney U test)
            # This test is for independent samples.
            # It returns (statistic, pvalue)
            if not group1_nes.empty and not group2_nes.empty:
                statistic, p_value = ranksums(group1_nes, group2_nes)
                p_value_str = f"p = {p_value:.3f}"
                if p_value < 0.001:
                    p_value_str = "p < 0.001"
            else:
                p_value_str = "Insufficient data for test"
                p_value = 1.0 # Assign a non-significant p-value if data is missing

            print(f"--- {mp} ---, p-value: {p_value_str}")

            # --- 3. Create an independent Boxplot for the current MP ---
            plt.figure(figsize=(8, 6)) # Adjust figure size for each plot

            # Ensure DiseaseState is categorical and NES is float for seaborn
            mp_df.loc[:, 'NES'] = mp_df['NES'].astype(float)

            if dataID == "E-TABM-1112":
                mp_df.loc[:, 'DiseaseState'] = mp_df['DiseaseState'].astype('category')
                sns.boxplot(data=mp_df, x='DiseaseState', y='NES', hue='DiseaseState', palette='viridis')

            elif dataID == "E-MTAB-1951":
                mp_df.loc[:, 'risk_score_factor'] = mp_df['risk_score_factor'].astype('category')
                sns.boxplot(data=mp_df, x='risk_score_factor', y='NES', hue='risk_score_factor', palette='viridis')

            # --- 4. Customize and show the plot ---
            plt.title(f'NES Score for {mp} by Disease State\n({p_value_str})', fontsize=14)
            plt.xlabel('Disease State', fontsize=12)
            plt.ylabel('Normalized Enrichment Score (NES)', fontsize=12)
            plt.grid(axis='y', linestyle='--', alpha=0.7)

            # Add p-value annotation directly on the plot
            # Find the max NES value to place the annotation slightly above the boxplots
            max_nes = mp_df['NES'].max()
            min_nes = mp_df['NES'].min()
            y_pos = max_nes + (max_nes - min_nes) * 0.1 # Position above the highest boxplot

            # Add a line connecting the two groups for visual clarity of comparison
            # Get x-coordinates for the two disease states (0 and 1 for the first two categories on the x-axis)
            x1, x2 = 0, 1
            plt.plot([x1, x2], [y_pos, y_pos], color='black', lw=1.5) # Horizontal line
            plt.text((x1 + x2) / 2, y_pos + (max_nes - min_nes) * 0.03, p_value_str,
                    ha='center', va='bottom', color='black', fontsize=10,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'))

            plt.tight_layout()

            # --- Save the figure ---
            # Sanitize the MP name for use in a filename (replace spaces and special characters)
            mp_filename = "".join(c if c.isalnum() else "_" for c in mp)
            output_path = os.path.join(output_base_dir, f"boxplot_{mp_filename}.pdf")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Saved figure: {output_path}")

            # --- Show and close the plot ---
            # plt.show() # Uncomment this line if you want to display each plot as it's generated
            plt.close() # Close the plot to free up memory, especially when generating many figures