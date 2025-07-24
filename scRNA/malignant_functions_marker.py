# Scanpy-based Pathway Enrichment and Bulk Survival Analysis
from math import log
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import gseapy as gp
from gseapy import Msigdb
from scipy.stats import hypergeom, ranksums, zscore
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

# ================================================================================
# STEP 1: IDENTIFY MARKER GENES FOR EACH MALIGNANT SUBTYPE
# ================================================================================

def identify_subtype_markers(malignant_cells, groupby='Malignant_type', 
                           method='wilcoxon', min_pct=0.25, logfc_threshold=0.25):
    """
    Identify marker genes for each malignant subtype using scanpy
    """
    print("=== IDENTIFYING SUBTYPE MARKER GENES ===")
    
    # Check if the groupby column exists
    if groupby not in malignant_cells.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in malignant_cells.obs")
    
    # Print subtype information
    subtype_counts = malignant_cells.obs[groupby].value_counts()
    print(f"Subtypes found: {len(subtype_counts)}")
    for subtype, count in subtype_counts.items():
        print(f"  {subtype}: {count} cells")
    
    # Find marker genes using scanpy
    print(f"Finding marker genes using {method} test...")
    sc.tl.rank_genes_groups(
        malignant_cells, 
        groupby=groupby, 
        method=method,
        use_raw=True if malignant_cells.raw is not None else False,
        min_pct=min_pct,
        logfc_threshold=logfc_threshold,
        key_added='subtype_markers'
    )
    
    # Extract marker genes into a dictionary
    marker_dict = {}
    subtypes = malignant_cells.obs[groupby].unique()
    
    for subtype in subtypes:
        # Get top marker genes for this subtype
        markers_df = sc.get.rank_genes_groups_df(
            malignant_cells, 
            group=subtype, 
            key='subtype_markers'
        )
        
        # Filter for significant markers
        significant_markers = markers_df[
            (markers_df['pvals_adj'] < 0.05) & 
            (markers_df['logfoldchanges'] > logfc_threshold)
        ]
        
        marker_dict[subtype] = {
            'signature_genes': significant_markers['names'].head(50).tolist(),
            'top_genes': significant_markers['names'].head(20).tolist(),  # Top 20 for visualization
            'all_markers': significant_markers
        }
        
        print(f"  {subtype}: {len(marker_dict[subtype]['signature_genes'])} significant markers")
    
    return marker_dict

# ================================================================================
# STEP 2: HALLMARK PATHWAY ANALYSIS
# ================================================================================

def get_hallmark_pathways(use_IMC_pathway=False):
    """
    Get HALLMARK pathways from MSigDB or use local definitions
    """
    print("=== LOADING HALLMARK PATHWAYS ===")
    msig = Msigdb()
    hallmark_pathways = msig.get_gmt(category='h.all', dbver="2025.1.Hs")
    select_pathways = hallmark_pathways

    if not use_IMC_pathway:
        print(f"Loaded {len(select_pathways)} HALLMARK pathways from MSigDB")

    else:
        select_pathways = {}
        # Local HALLMARK pathway definitions (subset)
        imc_relevant_hallmarks = {
            # Metabolism pathways (align with GLUT1, HK2, FASN)
            'HALLMARK_GLYCOLYSIS': 'Genes involved in glycolysis (→ GLUT1, HK2)',
            'HALLMARK_OXIDATIVE_PHOSPHORYLATION': 'Genes encoding oxidative phosphorylation machinery',
            'HALLMARK_FATTY_ACID_METABOLISM': 'Genes regulating fatty acid metabolism (→ FASN)',
            
            # Proliferation pathways (align with Ki67)
            'HALLMARK_E2F_TARGETS': 'Genes encoding E2F transcription factor targets (→ Ki67)',
            'HALLMARK_G2M_CHECKPOINT': 'Genes involved in G2/M checkpoint (→ Ki67)',
            'HALLMARK_MYC_TARGETS_V1': 'Genes encoding MYC targets (proliferation)',
            'HALLMARK_MYC_TARGETS_V2': 'Genes encoding MYC targets (metabolic)',
            
            # Stress/Hypoxia pathways (align with CA-IX, VEGF)
            'HALLMARK_HYPOXIA': 'Genes up-regulated in hypoxia (→ CA-IX, VEGF)',
            'HALLMARK_ANGIOGENESIS': 'Genes involved in angiogenesis (→ VEGF)',
            
            # EMT/Invasion pathways (align with Vimentin, EpCAM)
            'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION': 'Genes defining EMT (→ VIM, EPCAM)',
            
            # Growth/Anabolic pathways
            'HALLMARK_MTORC1_SIGNALING': 'Genes up-regulated by mTORC1 activation',
            
            # DNA damage/Cell cycle control
            'HALLMARK_P53_PATHWAY': 'Genes involved in p53 pathway and DNA damage response',
            'HALLMARK_DNA_REPAIR': 'Genes involved in DNA repair',
            
            # Apoptosis (complement to proliferation)
            'HALLMARK_APOPTOSIS': 'Genes mediating programmed cell death'
        }

        # Extract only the IMC-relevant pathways
        for pathway_name, description in imc_relevant_hallmarks.items():
            if pathway_name in hallmark_pathways:
                if pathway_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION":
                    select_pathways[pathway_name] = hallmark_pathways[pathway_name]
                    select_pathways[pathway_name].extend(["H3-3A","H3-3B","NME1-NME2"])
                else:
                    select_pathways[pathway_name] = hallmark_pathways[pathway_name]


                print(f"✅ {pathway_name}: {len(hallmark_pathways[pathway_name])} genes")
            else:
                print(f"⚠️ {pathway_name}: Not found in HALLMARK database")

        print(f"Using {len(select_pathways)} local HALLMARK pathways")

    return select_pathways

def perform_pathway_enrichment(marker_dict, hallmark_pathways, gene_size=20000, n_top_gene = None):
    """
    Perform pathway enrichment analysis for each subtype
    """
    print("=== PERFORMING PATHWAY ENRICHMENT ANALYSIS ===")
    
    enrichment_results = {}
    
    for subtype, marker_info in marker_dict.items():
        print(f"Analyzing {subtype}...")
        
        subtype_enrichment = {}

        if n_top_gene:
            signature_genes = set(marker_info['all_markers']['names'].head(n_top_gene).tolist())
        else:
            signature_genes = set(marker_info['all_markers']['names'].tolist())
        
        for pathway_name, pathway_genes in hallmark_pathways.items():
            pathway_genes_set = set(pathway_genes)
            
            # Calculate overlap
            overlap = signature_genes & pathway_genes_set
            overlap_count = len(overlap)
            
            if overlap_count > 0:
                # Hypergeometric test
                # Population size (total genes considered)
                population_size = gene_size  # Approximate human genome protein-coding genes
                
                # Sample size (signature genes)
                sample_size = len(signature_genes)
                
                # Successes in population (pathway genes)
                success_population = len(pathway_genes_set)
                
                # Observed successes (overlap)
                observed_successes = overlap_count
                
                # Calculate p-value using hypergeometric test
                p_value = hypergeom.sf(observed_successes - 1, population_size, 
                                     success_population, sample_size)
                
                # Calculate enrichment score
                expected = (sample_size * success_population) / population_size
                enrichment_score = observed_successes / expected if expected > 0 else 0
                
                subtype_enrichment[pathway_name] = {
                    'overlap_genes': list(overlap),
                    'overlap_count': overlap_count,
                    'pathway_size': len(pathway_genes_set),
                    'signature_size': len(signature_genes),
                    'p_value': p_value,
                    'enrichment_score': enrichment_score,
                    'expected': expected
                }
        
        # Filter for significant enrichments
        significant_enrichments = {
            k: v for k, v in subtype_enrichment.items() 
            if v['p_value'] < 0.05 and v['enrichment_score'] > 1.0
        }
        
        enrichment_results[subtype] = significant_enrichments
        print(f"  {subtype}: {len(significant_enrichments)} significant pathways")
    
    return enrichment_results

def create_pathway_enrichment_figure(enrichment_results, save_path=None, figsize=(16, 12)):
    """
    Create comprehensive pathway enrichment visualization
    """
    print("=== CREATING PATHWAY ENRICHMENT FIGURE ===")
    
    if not enrichment_results:
        print("No enrichment results to visualize")
        return None
    
    # Prepare data for visualization
    all_pathways = set()
    for subtype_results in enrichment_results.values():
        all_pathways.update(subtype_results.keys())
    
    all_pathways = sorted(list(all_pathways))
    all_subtypes = sorted(enrichment_results.keys())
    
    if not all_pathways:
        print("No significant pathways found")
        return None
    
    # Create enrichment matrix
    enrichment_matrix = np.zeros((len(all_pathways), len(all_subtypes)))
    pvalue_matrix = np.ones((len(all_pathways), len(all_subtypes)))
    
    for j, subtype in enumerate(all_subtypes):
        for i, pathway in enumerate(all_pathways):
            if pathway in enrichment_results[subtype]:
                enrichment_matrix[i, j] = enrichment_results[subtype][pathway]['enrichment_score']
                pvalue_matrix[i, j] = enrichment_results[subtype][pathway]['p_value']
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=figsize, 
                            gridspec_kw={'height_ratios': [3, 1], 'width_ratios': [4, 1]})
    
    # Main heatmap
    ax_main = axes[0, 0]
    
    # Log transform for better visualization
    enrichment_matrix_log = np.log2(enrichment_matrix + 1)
    
    im = ax_main.imshow(enrichment_matrix_log, aspect='auto', cmap='Reds', 
                       vmin=0, vmax=np.percentile(enrichment_matrix_log[enrichment_matrix_log > 0], 95))
    
    # Add significance markers
    for i in range(len(all_pathways)):
        for j in range(len(all_subtypes)):
            if pvalue_matrix[i, j] < 0.001:
                marker = '***'
            elif pvalue_matrix[i, j] < 0.01:
                marker = '**'
            elif pvalue_matrix[i, j] < 0.05:
                marker = '*'
            else:
                marker = ''
            
            if marker and enrichment_matrix[i, j] > 1:
                ax_main.text(j, i, marker, ha='center', va='center', 
                           fontsize=8, fontweight='bold', color='white')
    
    # Format main heatmap
    ax_main.set_xticks(range(len(all_subtypes)))
    ax_main.set_xticklabels(all_subtypes, rotation=45, ha='right', fontsize=10)
    ax_main.set_yticks(range(len(all_pathways)))
    
    clean_pathway_names = [p.replace('HALLMARK_', '').replace('_', ' ').title() 
                          for p in all_pathways]
    ax_main.set_yticklabels(clean_pathway_names, fontsize=9)
    
    ax_main.set_title('Malignant Subtype Pathway Enrichment', 
                     fontsize=14, fontweight='bold')
    ax_main.set_xlabel('Malignant Subtypes', fontsize=12)
    ax_main.set_ylabel('HALLMARK Pathways', fontsize=12)
    
    # Colorbar
    cbar_ax = axes[0, 1]
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label('log₂(Enrichment Score + 1)', rotation=270, labelpad=20, fontsize=11)
    
    # Network visualization
    ax_network = axes[1, 0]
    
    # Create network for top enrichments
    G = nx.Graph()
    
    # Add subtype nodes
    subtype_colors = plt.cm.Set2(np.linspace(0, 1, len(all_subtypes)))
    pos = {}
    
    # Position subtypes in a circle
    subtype_angles = np.linspace(0, 2*np.pi, len(all_subtypes), endpoint=False)
    subtype_radius = 1.5
    
    for i, subtype in enumerate(all_subtypes):
        x = subtype_radius * np.cos(subtype_angles[i])
        y = subtype_radius * np.sin(subtype_angles[i])
        pos[subtype] = (x, y)
        G.add_node(subtype, node_type='subtype', color=subtype_colors[i])
    
    # Add top pathway nodes and edges
    top_pathways = []
    for subtype in all_subtypes:
        if subtype in enrichment_results:
            subtype_enrichments = enrichment_results[subtype]
            top_subtype_pathways = sorted(subtype_enrichments.items(), 
                                        key=lambda x: x[1]['enrichment_score'], reverse=True)[:2]
            for pathway, data in top_subtype_pathways:
                if data['p_value'] < 0.05 and data['enrichment_score'] > 1.5:
                    top_pathways.append(pathway)
    
    top_pathways = list(set(top_pathways))  # Remove duplicates
    
    if top_pathways:
        pathway_angles = np.linspace(0, 2*np.pi, len(top_pathways), endpoint=False)
        pathway_radius = 0.8
        
        for i, pathway in enumerate(top_pathways):
            clean_name = pathway.replace('HALLMARK_', '').replace('_', '\n')
            if len(clean_name) > 15:
                clean_name = clean_name[:12] + '...'
            
            x = pathway_radius * np.cos(pathway_angles[i])
            y = pathway_radius * np.sin(pathway_angles[i])
            pos[clean_name] = (x, y)
            G.add_node(clean_name, node_type='pathway', color='lightblue')
            
            # Add edges for significant enrichments
            for subtype in all_subtypes:
                if subtype in enrichment_results and pathway in enrichment_results[subtype]:
                    enrichment_score = enrichment_results[subtype][pathway]['enrichment_score']
                    p_value = enrichment_results[subtype][pathway]['p_value']
                    
                    if p_value < 0.05 and enrichment_score > 1.5:
                        G.add_edge(subtype, clean_name, weight=enrichment_score)
        
        # Draw network
        subtype_nodes = [node for node in G.nodes() if G.nodes[node]['node_type'] == 'subtype']
        pathway_nodes = [node for node in G.nodes() if G.nodes[node]['node_type'] == 'pathway']
        
        if subtype_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=subtype_nodes, 
                                 node_color=[G.nodes[node]['color'] for node in subtype_nodes],
                                 node_size=800, ax=ax_network, alpha=0.8)
        
        if pathway_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=pathway_nodes,
                                 node_color='lightblue', node_size=400, ax=ax_network, alpha=0.6)
        
        # Draw edges
        edges = G.edges()
        if edges:
            weights = [G[u][v]['weight'] for u, v in edges]
            nx.draw_networkx_edges(G, pos, edgelist=edges, 
                                 width=[max(1, w/2) for w in weights], alpha=0.6, ax=ax_network)
        
        # Add labels
        nx.draw_networkx_labels(G, pos, labels={node: node for node in subtype_nodes}, 
                               font_size=8, font_weight='bold', ax=ax_network)
        nx.draw_networkx_labels(G, pos, labels={node: node for node in pathway_nodes}, 
                               font_size=6, ax=ax_network)
    
    ax_network.set_title('Top Enrichments Network', fontsize=12, fontweight='bold')
    ax_network.axis('off')
    
    # Summary table
    ax_table = axes[1, 1]
    ax_table.axis('off')
    
    # Create summary text
    summary_text = f"Analysis Summary:\n"
    summary_text += f"• {len(all_subtypes)} Malignant Subtypes\n"
    summary_text += f"• {len(all_pathways)} Enriched Pathways\n"
    summary_text += f"• Significance: * p<0.05, ** p<0.01, *** p<0.001\n"
    summary_text += f"• Enrichment Score = Observed/Expected\n"
    
    ax_table.text(0.1, 0.5, summary_text, transform=ax_table.transAxes, 
                 fontsize=10, verticalalignment='center',
                 bbox=dict(boxstyle="round,pad=0.5", facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(f"{save_path}/malignant_subtype_pathway_enrichment.pdf", 
                   dpi=300, bbox_inches='tight')
        print(f"Saved pathway enrichment figure to: {save_path}")
    
    plt.show()
    return fig, enrichment_results

# ================================================================================
# STEP 3: BULK SURVIVAL ANALYSIS
# ================================================================================

def calculate_subtype_scores_ssgsea(bulk_exp: pd.DataFrame, marker_dict: dict, min_genes: int = 5) -> pd.DataFrame:
    """
    Calculates single-sample Gene Set Enrichment Analysis (ssGSEA) scores for
    each subtype in each sample of a bulk RNA-seq dataset, using marker gene lists.

    Args:
        bulk_exp (pd.DataFrame): DataFrame of bulk expression data, with genes
                                 as rows and samples as columns. Typically log-normalized counts.
        marker_dict (dict): A dictionary where keys are subtype names (str) and values
                            are dictionaries containing 'signature_genes' (list of str).
                            Example: {'SubtypeA': {'signature_genes': ['GENE1', 'GENE2', ...]}}
        min_genes (int): Minimum number of available signature genes required for ssGSEA.
                         Gene sets with fewer genes in the expression data will be skipped
                         by gseapy, or result in NaN scores.

    Returns:
        pd.DataFrame: A DataFrame with 'Sample_ID', 'Subtype', 'Score' (ssGSEA NES),
                      'Available_genes', and 'Total_genes' columns.
    """
    print("=== CALCULATING SUBTYPE ENRICHMENT SCORES (ssGSEA) IN BULK DATA ===")

    # --- Input Validation ---
    if not isinstance(bulk_exp, pd.DataFrame):
        raise TypeError("bulk_exp must be a pandas DataFrame.")
    if bulk_exp.empty:
        print("Warning: bulk_exp DataFrame is empty. Returning empty scores.")
        return pd.DataFrame(columns=['Sample_ID', 'Subtype', 'Score', 'Available_genes', 'Total_genes'])

    if not isinstance(marker_dict, dict) or not marker_dict:
        print("Warning: marker_dict is empty or not a dictionary. Returning empty scores.")
        return pd.DataFrame(columns=['Sample_ID', 'Subtype', 'Score', 'Available_genes', 'Total_genes'])

    # --- Prepare gene_sets for gseapy ---
    # gseapy expects gene_sets as {'set_name': ['gene1', 'gene2', ...]}
    gseapy_gene_sets = {}
    total_genes_map = {} # To store total genes for reporting
    for subtype, marker_info in marker_dict.items():
        if 'signature_genes' in marker_info and isinstance(marker_info['signature_genes'], list):
            # Filter marker genes to only those present in bulk_exp index
            # This is important as gseapy expects genes to be in the expression data
            # and it will filter internally, but we want to report 'available_genes' accurately
            valid_genes = [g for g in marker_info['signature_genes'] if g in bulk_exp.index]
            if len(valid_genes) >= min_genes:
                gseapy_gene_sets[subtype] = valid_genes
                total_genes_map[subtype] = len(marker_info['signature_genes'])
            else:
                print(f"  Warning: Subtype '{subtype}' has only {len(valid_genes)} genes (out of {len(marker_info['signature_genes'])}) present in bulk_exp. Skipping for ssGSEA as it's below min_genes ({min_genes}).")
        else:
            print(f"  Warning: Marker information for subtype '{subtype}' is missing 'signature_genes' or it's not a list. Skipping.")

    if not gseapy_gene_sets:
        print("No valid gene sets found for ssGSEA after filtering for min_genes and presence in bulk_exp. Returning empty scores.")
        return pd.DataFrame(columns=['Sample_ID', 'Subtype', 'Score', 'Available_genes', 'Total_genes'])

    # --- Run ssGSEA using gseapy ---
    print(f"Running ssGSEA for {len(gseapy_gene_sets)} subtypes across {bulk_exp.shape[1]} samples...")
    try:
        ssgsea_results = gp.ssgsea(
            data=bulk_exp,
            gene_sets=gseapy_gene_sets,
            min_size=min_genes,  # Use min_genes for gseapy's internal filtering
            scale=True,          # Normalizes scores by default (Good for comparison)
            verbose=False,       # Suppress gseapy verbose output
            outdir=None          # Do not write results to disk by default
        )
    except Exception as e:
        print(f"Error during ssGSEA calculation: {e}")
        return pd.DataFrame(columns=['Sample_ID', 'Subtype', 'Score', 'Available_genes', 'Total_genes'])

    # --- Process gseapy output ---
    ssgsea_scores_df = ssgsea_results.res2d.copy()

    # Let's assume ssgsea_results.res2d has gene_sets as index and samples as columns.
    if isinstance(ssgsea_scores_df.index, pd.MultiIndex):
        # Handle case where index might be multi-level (e.g. from gsea)
        # For ssGSEA, it's usually just gene_set name as index
        ssgsea_scores_df.reset_index(inplace=True)
        ssgsea_scores_df.rename(columns={'index': 'Subtype'}, inplace=True)
    elif ssgsea_scores_df.index.name != 'Subtype':
        ssgsea_scores_df.index.name = 'Subtype'
        ssgsea_scores_df.reset_index(inplace=True)

    # Now, pivot to get samples as index, Subtype as column
    # Then melt to get Sample_ID, Subtype, Score (NES)
    scores_list = []
    for subtype in gseapy_gene_sets.keys():
        # Extract scores for this subtype across all samples
        subtype_scores_row = ssgsea_scores_df[ssgsea_scores_df['Term'] == subtype]
        # Get the NES column for all samples (assuming NES is always returned)
        # The column for sample_id is 'Name' in res2d, and score is 'NES'
        if 'Name' in subtype_scores_row.columns and 'NES' in subtype_scores_row.columns:
            for idx, row in subtype_scores_row.iterrows():
                sample_id = row['Name']
                es_score = row['ES']
                nes_score = row['NES']
                scores_list.append({
                    'Sample_ID': sample_id,
                    'Subtype': subtype,
                    'ES_Score': es_score,
                    'NES_Score': nes_score,
                    'Score': nes_score,
                    'Available_genes': len(gseapy_gene_sets[subtype]), # Genes actually used by ssGSEA
                    'Total_genes': total_genes_map.get(subtype, 0) # Original total genes
                })

    scores_df = pd.DataFrame(scores_list)

    if not scores_df.empty:
        print(f"Calculated ssGSEA scores for {len(scores_df['Sample_ID'].unique())} samples.")
        print(f"Subtypes scored: {', '.join(scores_df['Subtype'].unique())}.")
    else:
        print("No scores were calculated. Check input data, marker dictionary, and min_genes threshold.")

    print("✅ ssGSEA scores calculated successfully.")
    return scores_df


def calculate_subtype_scores_zscore(bulk_exp, marker_dict, score_method='mean', min_genes=5):
    """
    Calculate subtype signature scores in bulk RNA-seq data based on scRNA-seq marker lists.

    Parameters:
    - bulk_exp (pd.DataFrame): Bulk RNA-seq expression data. Genes should be
                                in the index and samples in columns.
                                (e.g., typically log-normalized counts or TPM).
    - marker_dict (dict): A dictionary where keys are subtype names (str) and values
                          are dictionaries containing 'signature_genes' (list of str).
                          Example: {'SubtypeA': {'signature_genes': ['GENE1', 'GENE2', ...]}}
    - score_method (str): Method to calculate the enrichment score.
                          'mean': Simple mean of available signature genes.
                          'zscore_mean': Calculate mean of Z-scored expression for available genes.
                                         (Z-scoring is done per gene across all samples in bulk_exp
                                         to standardize expression before averaging).
                          'gsva_like': A simplified, conceptual GSVA-like score (more advanced,
                                       would typically require dedicated libraries like `GSVA` in R
                                       or a Python equivalent for full implementation).
                                       For this function, 'mean' or 'zscore_mean' are more direct.
                                       Let's focus on 'mean' and 'zscore_mean' as readily implementable.
    - min_genes (int): Minimum number of available signature genes required to calculate a score.
                       If fewer genes are found, the score for that subtype will be NaN.

    Returns:
    - pd.DataFrame: A DataFrame with 'Sample_ID', 'Subtype', 'Score', 'Available_genes',
                    and 'Total_genes' columns.
    """
    print("=== CALCULATING SUBTYPE SCORES IN BULK DATA ===")

    # Ensure bulk_exp has genes as index and samples as columns
    if not isinstance(bulk_exp, pd.DataFrame):
        raise TypeError("bulk_exp must be a pandas DataFrame.")
    if bulk_exp.empty:
        print("Warning: bulk_exp DataFrame is empty. Returning empty scores.")
        return pd.DataFrame(columns=['Sample_ID', 'Subtype', 'Score', 'Available_genes', 'Total_genes'])

    # Optional: Z-score transformation of bulk expression data if requested
    processed_bulk_exp = bulk_exp.copy()
    if score_method == 'zscore_mean':
        print("Applying Z-score normalization per gene across all samples...")
        # Z-score each gene across all samples (row-wise operation)
        processed_bulk_exp = processed_bulk_exp.apply(zscore, axis=1)
        # Handle cases where zscore might return NaNs for constant genes
        processed_bulk_exp.fillna(0, inplace=True) # Replace NaNs (from std dev = 0) with 0

    scores_list = []

    # Iterate through each sample
    for sample_id in processed_bulk_exp.columns:
        sample_exp = processed_bulk_exp[sample_id]

        # Iterate through each subtype and its marker genes
        for subtype, marker_info in marker_dict.items():
            if 'signature_genes' not in marker_info or not isinstance(marker_info['signature_genes'], list):
                print(f"Warning: Marker information for subtype '{subtype}' is missing 'signature_genes' or it's not a list. Skipping.")
                continue

            signature_genes = marker_info['signature_genes']

            # Identify genes present in both marker list and bulk expression data
            available_genes = [g for g in signature_genes if g in sample_exp.index]

            if len(available_genes) >= min_genes:
                # Calculate enrichment score based on the chosen method
                if score_method == 'mean' or score_method == 'zscore_mean':
                    subtype_score = sample_exp[available_genes].mean()
                # You could add other scoring methods here if needed, e.g.,
                # elif score_method == 'sum':
                #     subtype_score = sample_exp[available_genes].sum()
                else:
                    print(f"Warning: Unknown score_method '{score_method}'. Using 'mean' by default for {subtype}.")
                    subtype_score = sample_exp[available_genes].mean()

                scores_list.append({
                    'Sample_ID': sample_id,
                    'Subtype': subtype,
                    'Score': subtype_score,
                    'Available_genes': len(available_genes),
                    'Total_genes': len(signature_genes)
                })
            else:
                # If not enough genes, append with NaN score
                scores_list.append({
                    'Sample_ID': sample_id,
                    'Subtype': subtype,
                    'Score': np.nan,  # Use NaN for scores that couldn't be calculated
                    'Available_genes': len(available_genes),
                    'Total_genes': len(signature_genes)
                })
                print(f"  Warning: Not enough available genes ({len(available_genes)} < {min_genes}) for Subtype '{subtype}' in Sample '{sample_id}'. Score set to NaN.")

    scores_df = pd.DataFrame(scores_list)

    if not scores_df.empty:
        print(f"Calculated scores for {len(scores_df['Sample_ID'].unique())} samples.")
        print(f"Subtypes scored: {', '.join(scores_df['Subtype'].unique())}.")
    else:
        print("No scores were calculated. Check input data and marker dictionary.")

    return scores_df

def perform_survival_analysis(scores_df, clinical_data, time_col='OS_time', 
                            status_col='OS_status', save_path=None):
    """
    Perform survival analysis for each subtype
    """
    print("=== PERFORMING SURVIVAL ANALYSIS ===")
    
    # Merge scores with clinical data
    clinical_data['Sample_ID'] = [f'X{x}' for x in clinical_data['Patient_ID']]
    scores_df['Sample_ID'] = scores_df['Sample_ID'].astype(str)
    
    merged_df = pd.merge(scores_df, clinical_data, on='Sample_ID')
    
    if time_col not in merged_df.columns or status_col not in merged_df.columns:
        print(f"Error: Required columns {time_col} or {status_col} not found")
        return None
    
    survival_results = {}
    subtypes = merged_df['Subtype'].unique()
    
    # Create figure for all survival curves
    n_subtypes = len(subtypes)
    n_cols = min(3, n_subtypes)
    n_rows = (n_subtypes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_subtypes == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for idx, subtype in enumerate(subtypes):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]
        
        print(f"Analyzing {subtype}...")
        
        # Get data for this subtype
        subtype_data = merged_df[merged_df['Subtype'] == subtype].copy()
        
        # Split into high/low based on median
        median_score = subtype_data['Score'].median()
        subtype_data['Score_group'] = subtype_data['Score'].apply(
            lambda x: 'High' if x > median_score else 'Low'
        )
        
        # Prepare survival data
        high_group = subtype_data[subtype_data['Score_group'] == 'High']
        low_group = subtype_data[subtype_data['Score_group'] == 'Low']
        
        if len(high_group) < 5 or len(low_group) < 5:
            print(f"  Insufficient samples for {subtype} (High: {len(high_group)}, Low: {len(low_group)})")
            continue
        
        # Kaplan-Meier analysis
        kmf_high = KaplanMeierFitter()
        kmf_low = KaplanMeierFitter()
        
        kmf_high.fit(high_group[time_col], high_group[status_col], label=f'{subtype} High')
        kmf_low.fit(low_group[time_col], low_group[status_col], label=f'{subtype} Low')
        
        # Log-rank test
        try:
            logrank_result = logrank_test(
                high_group[time_col], low_group[time_col],
                high_group[status_col], low_group[status_col]
            )
            p_value = logrank_result.p_value
        except:
            p_value = 1.0
        
        # Plot survival curves
        kmf_high.plot_survival_function(ax=ax, color='red', alpha=0.7)
        kmf_low.plot_survival_function(ax=ax, color='blue', alpha=0.7)
        
        ax.set_title(f'{subtype} Signature\n(p = {p_value:.3f})', fontsize=12)
        ax.set_xlabel('Time')
        ax.set_ylabel('Survival Probability')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Store results
        survival_results[subtype] = {
            'p_value': p_value,
            'high_group_size': len(high_group),
            'low_group_size': len(low_group),
            'median_score': median_score
        }
        
        print(f"  {subtype}: p = {p_value:.3f} (High: {len(high_group)}, Low: {len(low_group)})")
    
    # Hide unused subplots
    for idx in range(n_subtypes, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]
        ax.set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(f"{save_path}_malignant_subtype_survival_analysis.pdf", 
                   dpi=300, bbox_inches='tight')
        print(f"Saved survival analysis to: {save_path}")
    
    plt.show()
    
    return survival_results

def compare_bulk_datasets(scores_df, bulk_clinical, save_path="./bulk_comparison"):
    """
    Compare subtype scores across a bulk dataset with boxplot visualization and
    statistical testing (Wilcoxon rank-sum test).

    Parameters:
    - scores_df (pd.DataFrame): DataFrame containing 'Sample_ID', 'Subtype', and 'Score'.
    - bulk_clinical (pd.DataFrame): DataFrame with 'Sample_ID' (or index as Sample_ID)
                                    and 'DiseaseState'.
    - save_path (str): Base path to save the output plot. A directory will be
                       created if it doesn't exist.

    Returns:
    - dict: A dictionary containing statistical results for each subtype comparison.
    """
    # Ensure save directory exists
    save_dir = os.path.dirname(save_path)
    if save_dir and not os.path.exists(save_dir):
        os.makedirs(save_dir)
        print(f"Created directory: {save_dir}")

    # Merge with clinical data
    # Ensure Sample_ID in bulk_clinical is a column, not just the index, for merging
    if 'Sample_ID' not in bulk_clinical.columns and bulk_clinical.index.name != 'Sample_ID':
        bulk_clinical = bulk_clinical.reset_index()
        bulk_clinical = bulk_clinical.rename(columns={bulk_clinical.index.name: 'Sample_ID'})
    elif bulk_clinical.index.name == 'Sample_ID':
        bulk_clinical = bulk_clinical.reset_index()

    bulk_clinical['Sample_ID'] = bulk_clinical['Sample_ID'].astype(str)
    scores_df['Sample_ID'] = scores_df['Sample_ID'].astype(str)

    merged_df = pd.merge(scores_df, bulk_clinical, on='Sample_ID')

    # Ensure 'Score' column is numeric
    merged_df['Score'] = pd.to_numeric(merged_df['Score'], errors='coerce')
    merged_df.dropna(subset=['Score'], inplace=True)

    # Perform statistical tests and create plots
    dataset_results = {}
    subtypes = merged_df['Subtype'].unique() # Use merged_df to ensure only existing subtypes are processed

    # Create subplots for all subtypes in this dataset
    n_subtypes = len(subtypes)
    if n_subtypes == 0:
        print("No subtypes found in the merged data to plot.")
        return {}

    n_cols = min(3, n_subtypes) # Changed to max 3 columns for better visualization
    n_rows = (n_subtypes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))

    # Flatten the axes array for easier iteration if it's 2D
    if n_rows * n_cols > 1:
        axes = axes.flatten()
    else: # If only one subplot, axes is not an array
        axes = [axes]

    plot_idx = 0

    print("Starting comparison and plotting for subtypes:")
    for i, subtype in enumerate(subtypes):
        subtype_data = merged_df[merged_df['Subtype'] == subtype]

        disease_states = subtype_data['DiseaseState'].unique()
        if len(disease_states) >= 2:
            # We are comparing only the first two unique disease states found
            # If you need to compare all pairs, the logic needs to be expanded.
            group1_name = disease_states[0]
            group2_name = disease_states[1]

            group1 = subtype_data[subtype_data['DiseaseState'] == group1_name]['Score']
            group2 = subtype_data[subtype_data['DiseaseState'] == group2_name]['Score']

            if len(group1) > 0 and len(group2) > 0:
                # Statistical test: Wilcoxon rank-sum test for non-parametric data
                statistic, p_value = ranksums(group1, group2)

                # Store results
                dataset_results[subtype] = {
                    'statistic': statistic,
                    'p_value': p_value,
                    'group1_mean': group1.mean(),
                    'group2_mean': group2.mean(),
                    'group1_name': group1_name,
                    'group2_name': group2_name,
                    'group1_n': len(group1),
                    'group2_n': len(group2)
                }

                ax = axes[plot_idx]

                # Create boxplot using seaborn
                sns.boxplot(data=subtype_data, x='DiseaseState', y='Score', ax=ax, palette='viridis',
                            order=[group1_name, group2_name]) # Ensure order for boxplot
                sns.stripplot(data=subtype_data, x='DiseaseState', y='Score', ax=ax,
                              color='black', alpha=0.6, size=3, order=[group1_name, group2_name])

                # Format p-value
                if p_value < 0.001:
                    p_text = "p < 0.001"
                else:
                    p_text = f"p = {p_value:.3f}"

                # Add significance annotation
                y_max = subtype_data['Score'].max()
                y_min = subtype_data['Score'].min()
                y_range = y_max - y_min
                # Adjust y_pos to be slightly above the max score, but within reasonable limits
                y_pos = y_max + y_range * 0.1
                if y_range == 0: # Handle case where all scores are the same
                    y_pos = y_max + 0.1 # Small arbitrary offset

                # Draw significance line (between the two groups)
                ax.plot([0, 1], [y_pos, y_pos], color='black', linewidth=1)
                ax.plot([0, 0], [y_pos, y_pos - y_range * 0.03], color='black', linewidth=1)
                ax.plot([1, 1], [y_pos, y_pos - y_range * 0.03], color='black', linewidth=1)

                # Add p-value text
                ax.text(0.5, y_pos + y_range * 0.05, p_text, ha='center', va='bottom',
                        fontsize=10, fontweight='bold')

                # Format plot
                ax.set_title(f'{subtype}', fontsize=12, fontweight='bold')
                ax.set_xlabel('Disease State', fontsize=10)
                ax.set_ylabel('Signature Score', fontsize=10)
                ax.grid(True, alpha=0.3)

                # Add sample sizes to x-axis labels
                current_labels = ax.get_xticklabels()
                new_labels = []
                for label_obj in current_labels:
                    state_name = label_obj.get_text()
                    n_samples = len(subtype_data[subtype_data['DiseaseState'] == state_name])
                    new_labels.append(f'{state_name}\n(n={n_samples})')
                ax.set_xticklabels(new_labels)

                plot_idx += 1
                print(f"  - Subtype '{subtype}': {p_text} (statistic: {statistic:.2f})")
            else:
                print(f"  - Skipping subtype '{subtype}': Not enough samples in one or both groups ({group1_name}, {group2_name}).")
        else:
            print(f"  - Skipping subtype '{subtype}': Less than two unique disease states found.")

    # Hide unused subplots
    for idx in range(plot_idx, n_rows * n_cols):
        axes[idx].set_visible(False)

    plt.tight_layout()

    # Save plot
    plot_filename = f"{save_path}_subtype_comparisons.pdf"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"\nComparison plot saved to: {plot_filename}")
    print("Statistical results:")
    for subtype, res in dataset_results.items():
        print(f"  {subtype}:")
        print(f"    - {res['group1_name']} (n={res['group1_n']}) Mean Score: {res['group1_mean']:.3f}")
        print(f"    - {res['group2_name']} (n={res['group2_n']}) Mean Score: {res['group2_mean']:.3f}")
        print(f"    - Wilcoxon Statistic: {res['statistic']:.3f}, p-value: {res['p_value']:.4f}")
    
    return dataset_results