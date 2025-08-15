import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import matplotlib.patches as mpatches

import warnings
warnings.filterwarnings('ignore')

from scipy import stats
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from gseapy import Msigdb
import gseapy as gp
import networkx as nx

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

#%% Visualize Data
def plot_sample_overview(metadata_df, figsize=(18, 12)):
    """
    Create sample overview plots based on actual data structure
    
    Parameters:
    -----------
    metadata_df : DataFrame
        Sample metadata with columns: sample, RFS_status, tma_core, core_name, x, y
    figsize : tuple
        Figure size
    """
    # Prepare sample-level summary data
    sample_summary = metadata_df.groupby(['sample', 'RFS_status', 'tma_core', 'core_name']).agg({
        'x': 'count',  # Count spots per sample
    }).reset_index()
    
    # Rename columns for clarity
    sample_summary = sample_summary.rename(columns={
        'x': 'n_spots',
        'sample': 'sample_id'
    })
    
    sample_summary = sample_summary.loc[sample_summary["n_spots"] > 0]

    # Set up the figure with better layout
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    fig.suptitle('Sample Overview - Quality Control & Distribution Analysis', 
                 fontsize=16, fontweight='bold')
    
    # Define consistent colors
    rfs_colors = {0: '#2E86AB', 1: '#F24236'}  # Blue for RFS=0, Red for RFS=1
    
    # 1. RFS Status Distribution (Pie Chart)
    ax1 = axes[0]
    rfs_counts = metadata_df['RFS_status'].value_counts().sort_index()
    colors = [rfs_colors[status] for status in rfs_counts.index]
    
    wedges, texts, autotexts = ax1.pie(rfs_counts.values, 
                                      labels=[f'RFS {status}' for status in rfs_counts.index],
                                      autopct='%1.1f%%', 
                                      colors=colors, 
                                      startangle=90,
                                      explode=[0.05, 0.05])
    
    ax1.set_title('RFS Status Distribution\n(Overall Spots)', fontweight='bold')
    
    # Enhance pie chart text
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(11)
    
    # 2. TMA Core Distribution
    ax2 = axes[1]
    
    sns.barplot(data=sample_summary, x='core_name', y='n_spots', 
                hue='RFS_status', palette=rfs_colors, ax=ax2)
    ax2.set_title('Spots per TMA Core', fontweight='bold')
    ax2.set_xlabel('TMA Core', fontweight='bold')
    ax2.set_ylabel('Number of Spots', fontweight='bold')
    ax2.legend(title='RFS Status', title_fontsize=10, fontsize=9)
    ax2.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for container in ax2.containers:
        ax2.bar_label(container, fontsize=9, padding=3)
    
    # 4. Box plot comparison
    ax4 = axes[2]
    
    # Create comparison data
    comparison_data = []
    for _, row in sample_summary.iterrows():
        comparison_data.append({
            'RFS_Status': f"RFS {row['RFS_status']}",
            'Spots_per_Sample': row['n_spots'],
            'TMA_Core': row['tma_core'],
            'Core_Name': row['core_name']
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    sns.boxplot(data=comparison_df, x='RFS_Status', y='Spots_per_Sample', 
                palette=[rfs_colors[0], rfs_colors[1]], ax=ax4)
    sns.stripplot(data=comparison_df, x='RFS_Status', y='Spots_per_Sample', 
                 color='black', alpha=0.6, size=6, ax=ax4)
    
    ax4.set_title('Spots per Sample\nby RFS Status', fontweight='bold')
    ax4.set_xlabel('RFS Status', fontweight='bold')
    ax4.set_ylabel('Spots per Sample', fontweight='bold')
    ax4.grid(axis='y', alpha=0.3)
    
    # Add statistical annotation
    rfs0_spots = comparison_df[comparison_df['RFS_Status'] == 'RFS 0']['Spots_per_Sample']
    rfs1_spots = comparison_df[comparison_df['RFS_Status'] == 'RFS 1']['Spots_per_Sample']
    
    if len(rfs0_spots) > 0 and len(rfs1_spots) > 0:
        stat, pval = stats.mannwhitneyu(rfs0_spots, rfs1_spots, alternative='two-sided')
        
        y_max = comparison_df['Spots_per_Sample'].max()
        ax4.text(0.5, y_max * 1, f'p = {pval:.3f}', 
                 ha='center', va='bottom', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout()
    
    return fig

def clean_cell_type_names(cell_type_cols):
    """Clean cell type names by removing prefix"""
    clean_names = []
    for col in cell_type_cols:
        clean_name = col.replace('q05cell_abundance_w_sf_', '')
        clean_names.append(clean_name)
    return clean_names

def plot_cell_type_overview(abundance_df, figsize=(24, 16)):
    """
    Create comprehensive cell type abundance analysis
    
    Parameters:
    -----------
    abundance_df : DataFrame
        Cell abundance data with sample_id and RFS_status columns
    figsize : tuple
        Figure size
    """
    # Get cell type columns
    cell_type_cols = [col for col in abundance_df.columns 
                     if col not in ['sample_id', 'RFS_status']]
    
    # Clean cell type names
    clean_names = clean_cell_type_names(cell_type_cols)
    
    # Aggregate to slide level for graphs 2-4
    slide_df = aggregate_to_slide_level(abundance_df, cell_type_cols)
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    fig.suptitle('Cell Type Abundance Analysis (Slide-Level)', fontsize=16, fontweight='bold')
    
    # 1. Barplot - Overall fraction values (keep spot-level for better representation)
    ax1 = fig.add_subplot(gs[0, :2])
    
    # Calculate mean fractions from all spots
    mean_fractions = abundance_df[cell_type_cols].mean().sort_values(ascending=False)
    clean_names_sorted = [name.replace('q05cell_abundance_w_sf_', '') for name in mean_fractions.index]
    
    bars = ax1.bar(range(len(mean_fractions)), mean_fractions.values)
    ax1.set_title('Mean Fraction of Each Cell Subpopulation (All Spots)')
    ax1.set_xlabel('Cell Types')
    ax1.set_ylabel('Mean Fraction')
    ax1.set_xticks(range(len(mean_fractions)))
    ax1.set_xticklabels(clean_names_sorted, rotation=45, ha='right')
    
    # Add value labels on bars
    for i, v in enumerate(mean_fractions.values):
        ax1.text(i, v + 0.001, f'{v:.3f}', ha='center', va='bottom', fontsize=8)
    
    # 2. Heatmap - Cell types vs Slides (slide-level)
    ax2 = fig.add_subplot(gs[1, :])
    
    # Prepare data for heatmap using slide-level data
    heatmap_data = slide_df.set_index('sample_id')[cell_type_cols].T
    
    # Use clean names for y-axis
    heatmap_data.index = clean_names
    
    # Create heatmap
    sns.heatmap(heatmap_data, cmap='viridis', ax=ax2, cbar_kws={'label': 'Mean Fraction per Slide'})
    ax2.set_title('Cell Type Mean Fractions Across Slides')
    ax2.set_xlabel('Sample ID (Slides)')
    ax2.set_ylabel('Cell Types')
    
    # Rotate x-axis labels for better readability
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
    
    # 3. Volcano plot (slide-level)
    ax3 = fig.add_subplot(gs[2, 0])
    
    # Perform statistical analysis at slide level
    stats_results = perform_wilcox_test(slide_df, cell_type_cols)
    
    if len(stats_results) > 0:
        # Create volcano plot
        x_vals = stats_results['log2_fc']
        y_vals = -np.log10(stats_results['p_adj'])
        
        # Color points based on significance and fold change
        colors = []
        for _, row in stats_results.iterrows():
            if row['p_adj'] < 0.05 and abs(row['log2_fc']) > 0.5:
                colors.append('red')
            elif row['p_adj'] < 0.05:
                colors.append('orange')
            else:
                colors.append('gray')
        
        scatter = ax3.scatter(x_vals, y_vals, c=colors, alpha=0.7)
        ax3.set_xlabel('Log2 Fold Change (RFS_status 1 vs 0)')
        ax3.set_ylabel('-Log10 Adjusted P-value')
        ax3.set_title('Volcano Plot: RFS Status Comparison')
        
        # Add significance lines
        ax3.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        ax3.axvline(x=0.5, color='black', linestyle='--', alpha=0.5)
        ax3.axvline(x=-0.5, color='black', linestyle='--', alpha=0.5)
        
        # Add legend
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
                                    markersize=8, label='Significant & |FC|>1.4'),
                         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', 
                                    markersize=8, label='Significant'),
                         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
                                    markersize=8, label='Not significant')]
        ax3.legend(handles=legend_elements, loc='upper right')
    
    # 4. Boxplots for key cell types (slide-level)
    ax4 = fig.add_subplot(gs[2, 1:])
    
    # Select key cell types (top 6 most variable or significant)
    if len(stats_results) > 0:
        # Get top significant cell types
        significant_cells = stats_results[stats_results['p_adj'] < 0.1].sort_values('p_adj')
        significant_cells = significant_cells.sort_values(by = "log2_fc",ascending = False)
        key_cell_types = significant_cells.head(6)['cell_type'].tolist()
    else:
        # Fall back to most abundant
        key_cell_types = slide_df[cell_type_cols].mean().nlargest(8).index.tolist()
    
    # Prepare data for boxplot using slide-level data
    plot_data = []
    for cell_type in key_cell_types:
        clean_name = cell_type.replace('q05cell_abundance_w_sf_', '')
        for _, row in slide_df.iterrows():
            plot_data.append({
                'Cell_Type': clean_name,
                'Fraction': row[cell_type],
                'RFS_Status': f"Group_{row['RFS_status']}"
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create boxplot
    sns.boxplot(data=plot_df, x='Cell_Type', y='Fraction', hue='RFS_Status', ax=ax4)
    ax4.set_title('Key Cell Types: RFS Status Comparison (Slide-Level)')
    ax4.set_xlabel('Cell Types')
    ax4.set_ylabel('Mean Fraction per Slide')
    ax4.tick_params(axis='x', rotation=45)
    
    # Add p-values to boxplot
    y_max = plot_df['Fraction'].max()
    for i, cell_type in enumerate(key_cell_types):
        # Get p-value for this cell type
        if len(stats_results) > 0:
            p_val = stats_results[stats_results['cell_type'] == cell_type]['p_adj'].iloc[0]
            if p_val < 0.001:
                p_text = '***'
            elif p_val < 0.01:
                p_text = '**'
            elif p_val < 0.05:
                p_text = '*'
            else:
                p_text = 'ns'
            
            ax4.text(i, y_max * 0.85, p_text, ha='center', va='bottom', fontweight='bold')
    
    plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9,wspace=0.2,hspace=0.2)
    #plt.tight_layout()
    plt.show()
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"Total spots: {len(abundance_df)}")
    print(f"Total slides: {len(slide_df)}")
    print(f"Spots per slide (mean): {len(abundance_df) / len(slide_df):.1f}")
    print(f"RFS_status distribution (slides): {slide_df['RFS_status'].value_counts().to_dict()}")
    print(f"Total cell types analyzed: {len(cell_type_cols)}")
    
    if len(stats_results) > 0:
        sig_count = len(stats_results[stats_results['p_adj'] < 0.05])
        print(f"Significantly different cell types (p_adj < 0.05): {sig_count}")
        
        print("\nTop 5 most significant cell types (slide-level analysis):")
        top_sig = stats_results.nsmallest(5, 'p_adj')[['cell_type', 'log2_fc', 'p_adj', 'n_slides_group_0', 'n_slides_group_1']]
        top_sig['clean_name'] = top_sig['cell_type'].str.replace('q05cell_abundance_w_sf_', '')
        print(top_sig[['clean_name', 'log2_fc', 'p_adj', 'n_slides_group_0', 'n_slides_group_1']].to_string(index=False))
    
    return fig, stats_results

#%% Sample level analysis

def aggregate_to_slide_level(abundance_df, cell_type_cols):
    """Aggregate spot-level data to slide-level means"""
    # Group by sample_id and calculate mean for each slide
    slide_level_df = abundance_df.groupby(['sample_id', 'RFS_status'])[cell_type_cols].mean().reset_index()
    return slide_level_df

def perform_wilcox_test(slide_df, cell_type_cols):
    """Perform wilcoxon test between RFS_status groups at slide level"""
    results = []
    
    # Get groups at slide level
    group_0 = slide_df[slide_df['RFS_status'] == 0]
    group_1 = slide_df[slide_df['RFS_status'] == 1]
    
    for cell_type in cell_type_cols:
        values_0 = group_0[cell_type].values
        values_1 = group_1[cell_type].values
        
        # Calculate fold change (ratio of means)
        mean_0 = np.mean(values_0)
        mean_1 = np.mean(values_1)
        
        # Avoid division by zero
        if mean_0 > 0:
            fold_change = mean_1 / mean_0
        else:
            fold_change = np.inf if mean_1 > 0 else 1.0
        
        # Perform wilcoxon test
        try:
            stat, p_value = stats.mannwhitneyu(values_1, values_0, alternative='two-sided')
        except:
            p_value = 1.0
        
        results.append({
            'cell_type': cell_type,
            'fold_change': fold_change,
            'log2_fc': np.log2(fold_change) if fold_change > 0 and fold_change != np.inf else 0,
            'p_value': p_value,
            'mean_group_0': mean_0,
            'mean_group_1': mean_1,
            'n_slides_group_0': len(values_0),
            'n_slides_group_1': len(values_1)
        })
    
    results_df = pd.DataFrame(results)
    
    # BH adjustment
    if len(results_df) > 0:
        _, results_df['p_adj'], _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    
    return results_df

#%% Pathway analysis

def filter_pathways_by_keywords(pathways, keywords=None, exclude_keywords=None):
    """Filter pathways by including/excluding keywords"""
    if keywords is None and exclude_keywords is None:
        return pathways
    
    filtered = {}
    for k, v in pathways.items():
        include = True
        
        # Check inclusion keywords
        if keywords:
            include = any(keyword.lower() in k.lower() for keyword in keywords)
        
        # Check exclusion keywords
        if exclude_keywords and include:
            include = not any(keyword.lower() in k.lower() for keyword in exclude_keywords)
        
        if include:
            filtered[k] = v
    
    return filtered

def load_msigdb_pathways(category, dbver="2025.1.Hs", keywords=None, exclude_keywords=None):
    """Load MSigDB pathways for specific category"""
    try:
        msig = Msigdb()
        print(f"  Fetching {category} gene sets...")
        pathways = msig.get_gmt(category=category, dbver=dbver)
        
        # Apply filtering if specified
        if keywords or exclude_keywords:
            pathways = filter_pathways_by_keywords(pathways, keywords, exclude_keywords)
        
        return pathways
    except Exception as e:
        print(f"  ⚠️ Failed to load {category}: {e}")
        return {}

def load_gseapy_pathways(library_name, keywords=None, exclude_keywords=None, max_pathways=None):
    """Load pathways from gseapy libraries"""
    try:
        print(f"  Fetching {library_name}...")
        pathways = gp.get_library(name=library_name, organism='Human')
        
        # Apply filtering if specified
        if keywords or exclude_keywords:
            pathways = filter_pathways_by_keywords(pathways, keywords, exclude_keywords)
        
        # Limit number of pathways if specified
        if max_pathways and len(pathways) > max_pathways:
            pathways = dict(list(pathways.items())[:max_pathways])
        
        return pathways
    except Exception as e:
        print(f"  ⚠️ Failed to load {library_name}: {e}")
        return {}

def get_predefined_keywords(category_type):
    """Get predefined keywords for common pathway categories"""
    keyword_sets = {
        'metabolism': [
            'metabolism', 'metabolic', 'glycolysis', 'glycolytic',
            'oxidative phosphorylation', 'electron transport', 'atp synthesis',
            'fatty acid', 'lipid', 'cholesterol', 'steroid',
            'amino acid', 'protein', 'nucleotide', 'purine', 'pyrimidine',
            'tca cycle', 'citrate cycle', 'krebs cycle',
            'pentose phosphate', 'hexosamine', 'fructose',
            'gluconeogenesis', 'glucose', 'insulin signaling'
        ],
        'immune': [
            'immune', 'inflammation', 'interferon', 'interleukin',
            'tnf', 'nfkb', 't cell', 'b cell', 'macrophage',
            'cytokine', 'chemokine', 'antigen', 'antibody',
            'complement', 'toll like', 'innate immunity'
        ],
        'cancer': [
            'cancer', 'tumor', 'oncogene', 'tumor suppressor',
            'apoptosis', 'cell cycle', 'dna repair', 'p53',
            'proliferation', 'metastasis', 'angiogenesis',
            'epithelial mesenchymal', 'stem cell'
        ],
        'signaling': [
            'signaling', 'pathway', 'cascade', 'receptor',
            'kinase', 'phosphorylation', 'transcription',
            'wnt', 'notch', 'hedgehog', 'tgf', 'mapk',
            'pi3k', 'akt', 'mtor', 'jak stat'
        ],
        'development': [
            'development', 'differentiation', 'morphogenesis',
            'embryonic', 'organogenesis', 'cell fate',
            'specification', 'patterning', 'migration'
        ],
        'stress': [
            'stress', 'response', 'heat shock', 'oxidative',
            'hypoxia', 'unfolded protein', 'autophagy',
            'dna damage', 'repair', 'survival'
        ]
    }
    return keyword_sets.get(category_type, [])

def load_pathway_databases(database_name="HALLMARK", category_filter=None, 
                         custom_keywords=None, exclude_keywords=None, 
                         max_pathways=None):
    """
    Enhanced pathway database loader supporting multiple databases and filtering
    
    Parameters:
    -----------
    database_name : str
        Database to load ('HALLMARK', 'KEGG', 'REACTOME', 'GO_BP', 'GO_MF', 'GO_CC', 'ALL')
    category_filter : str
        Predefined category filter ('metabolism', 'immune', 'cancer', 'signaling', 'development', 'stress')
    custom_keywords : list
        Custom keywords to filter pathways
    exclude_keywords : list
        Keywords to exclude from pathways
    max_pathways : int
        Maximum number of pathways to return
    """
    print(f"\n--- Loading pathway databases: {database_name} ---")
    
    pathways = {}
    
    # Get keywords based on category filter
    if category_filter:
        filter_keywords = get_predefined_keywords(category_filter)
        if custom_keywords:
            filter_keywords.extend(custom_keywords)
    else:
        filter_keywords = custom_keywords
    
    if database_name == "HALLMARK":
        # MSigDB Hallmark pathways
        hallmark_pathways = load_msigdb_pathways(
            category='h.all', 
            keywords=filter_keywords, 
            exclude_keywords=exclude_keywords
        )
        pathways.update({f"HALLMARK_{k}": v for k, v in hallmark_pathways.items()})
        
    elif database_name == "KEGG":
        # KEGG pathways
        kegg_pathways = load_gseapy_pathways(
            library_name='KEGG_2021_Human',
            keywords=filter_keywords,
            exclude_keywords=exclude_keywords,
            max_pathways=max_pathways
        )
        pathways.update({f"KEGG_{k}": v for k, v in kegg_pathways.items()})
        
    elif database_name == "REACTOME":
        # Reactome pathways
        reactome_pathways = load_gseapy_pathways(
            library_name='Reactome_2022',
            keywords=filter_keywords,
            exclude_keywords=exclude_keywords,
            max_pathways=max_pathways
        )
        pathways.update({f"REACTOME_{k}": v for k, v in reactome_pathways.items()})
        
    elif database_name == "GO_BP":
        # Gene Ontology Biological Process
        go_bp_pathways = load_gseapy_pathways(
            library_name='GO_Biological_Process_2023',
            keywords=filter_keywords,
            exclude_keywords=exclude_keywords,
            max_pathways=max_pathways or 200  # Limit GO terms
        )
        pathways.update({f"GO_BP_{k}": v for k, v in go_bp_pathways.items()})
        
    elif database_name == "GO_MF":
        # Gene Ontology Molecular Function
        go_mf_pathways = load_gseapy_pathways(
            library_name='GO_Molecular_Function_2023',
            keywords=filter_keywords,
            exclude_keywords=exclude_keywords,
            max_pathways=max_pathways or 150
        )
        pathways.update({f"GO_MF_{k}": v for k, v in go_mf_pathways.items()})
        
    elif database_name == "GO_CC":
        # Gene Ontology Cellular Component
        go_cc_pathways = load_gseapy_pathways(
            library_name='GO_Cellular_Component_2023',
            keywords=filter_keywords,
            exclude_keywords=exclude_keywords,
            max_pathways=max_pathways or 100
        )
        pathways.update({f"GO_CC_{k}": v for k, v in go_cc_pathways.items()})
        
    elif database_name == "ALL":
        # Load multiple databases
        print("  Loading multiple databases...")
        
        # Core databases
        databases_to_load = [
            ("HALLMARK", 'h.all', None),
            ("KEGG", 'KEGG_2021_Human', 200),
            ("REACTOME", 'Reactome_2022', 150)
        ]
        
        for db_name, db_source, max_limit in databases_to_load:
            if db_name == "HALLMARK":
                db_pathways = load_msigdb_pathways(
                    category=db_source,
                    keywords=filter_keywords,
                    exclude_keywords=exclude_keywords
                )
            else:
                db_pathways = load_gseapy_pathways(
                    library_name=db_source,
                    keywords=filter_keywords,
                    exclude_keywords=exclude_keywords,
                    max_pathways=max_limit
                )
            pathways.update({f"{db_name}_{k}": v for k, v in db_pathways.items()})
    
    else:
        print(f"  ⚠️ Unknown database: {database_name}")
        return {}
    
    print(f"  ✅ Loaded {len(pathways)} pathways from {database_name}")
    
    # Apply overall limit if specified
    if max_pathways and len(pathways) > max_pathways:
        pathways = dict(list(pathways.items())[:max_pathways])
        print(f"  📊 Limited to {len(pathways)} pathways")
    
    return pathways

def perform_pathway_enrichment(gene_list, pathway_database, all_genes, p_threshold=0.05):
    """
    Perform pathway enrichment analysis using hypergeometric test
    
    Parameters:
    -----------
    gene_list : list
        List of genes of interest (significant DEGs)
    pathway_database : dict
        Dictionary of pathway_name: [gene_list]
    all_genes : list
        Background gene list (all tested genes)
    p_threshold : float
        P-value threshold for significance
    
    Returns:
    --------
    enrichment_results : DataFrame
        Enrichment analysis results
    """
    
    enrichment_results = []
    
    # Total number of genes in background
    total_genes = len(all_genes)
    # Number of significant genes
    sig_genes = len(gene_list)
    
    for pathway_name, pathway_genes in pathway_database.items():
        # Find overlap between pathway genes and all tested genes
        pathway_in_background = [g for g in pathway_genes if g in all_genes]
        pathway_size = len(pathway_in_background)
        
        if pathway_size < 5:  # Skip pathways with too few genes
            continue
            
        # Find overlap between significant genes and pathway
        overlap_genes = [g for g in gene_list if g in pathway_in_background]
        overlap_size = len(overlap_genes)
        
        if overlap_size == 0:  # Skip if no overlap
            continue
        
        # Hypergeometric test
        # Parameters: total_genes, pathway_size, sig_genes, overlap_size
        p_value = hypergeom.sf(overlap_size - 1, total_genes, pathway_size, sig_genes)
        
        # Calculate enrichment ratio
        expected = (pathway_size * sig_genes) / total_genes
        enrichment_ratio = overlap_size / expected if expected > 0 else 0
        
        enrichment_results.append({
            'pathway': pathway_name,
            'pathway_size': pathway_size,
            'overlap_size': overlap_size,
            'expected': expected,
            'enrichment_ratio': enrichment_ratio,
            'p_value': p_value,
            'overlap_genes': overlap_genes
        })
    
    # Convert to DataFrame and adjust p-values
    if enrichment_results:
        enrichment_df = pd.DataFrame(enrichment_results)
        enrichment_df = enrichment_df.sort_values('p_value')
        
        # Multiple testing correction
        _, enrichment_df['p_adjusted'], _, _ = multipletests(
            enrichment_df['p_value'], method='fdr_bh'
        )
        
        # Filter for significant pathways
        significant_pathways = enrichment_df[enrichment_df['p_adjusted'] < p_threshold].copy()
        
        return significant_pathways
    else:
        return pd.DataFrame()

def plot_pathway_enrichment(enrichment_results, title="Pathway Enrichment", max_pathways=20):
    """
    Create a dot plot for pathway enrichment results
    """
    if len(enrichment_results) == 0:
        return None
        
    # Take top pathways
    plot_data = enrichment_results.head(max_pathways).copy()
    plot_data = plot_data.sort_values('enrichment_ratio', ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create dot plot
    y_pos = range(len(plot_data))
    
    # Color by p-value, size by gene count
    colors = -np.log10(plot_data['p_adjusted'])
    sizes = plot_data['overlap_size'] * 20  # Scale for visibility
    
    scatter = ax.scatter(plot_data['enrichment_ratio'], y_pos, 
                        c=colors, s=sizes, alpha=0.7, 
                        cmap='Reds', edgecolors='black', linewidth=0.5)
    
    # Customize plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels([p.replace('_', ' ').replace('HALLMARK ', '') for p in plot_data['pathway']], 
                       fontsize=10)
    ax.set_xlabel('Enrichment Ratio', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add vertical line at enrichment ratio = 1
    ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
    
    # Add colorbar for p-values
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-Log10(Adjusted P-value)', fontsize=10)
    
    # Add size legend
    sizes_legend = [5, 10, 20]
    size_labels = ['5', '10', '20']
    legend_elements = [plt.scatter([], [], s=s*20, c='gray', alpha=0.7, edgecolors='black') 
                      for s in sizes_legend]
    legend1 = ax.legend(legend_elements, size_labels, title='Gene Count', 
                       loc='lower right', frameon=True, fontsize=8)
    ax.add_artist(legend1)
    
    # Add statistics text
    stats_text = f"Pathways shown: {len(plot_data)}\nTotal significant: {len(enrichment_results)}"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            verticalalignment='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    return fig

#%% Distance analysis
def identify_cholangiocyte_zones_graph(sample_adata, max_distance_um=30, scale_factor=0.5):
    """
    Identify cholangiocyte zones using graph-based approach (IMC-aligned)
    
    Parameters:
    -----------
    sample_adata : AnnData
        Single sample spatial data
    max_distance_um : float
        Maximum distance (μm) to connect cholangiocytes (default: 30μm from IMC)
    scale_factor : float
        Conversion factor from pixels to micrometers
    
    Returns:
    --------
    cholangiocyte_zones : dict
        Dictionary of cholangiocyte zone information
    """
    # Get cholangiocyte spots
    chol_mask = sample_adata.obs['cholangiocyte_enriched'] == True
    cholangiocyte_spots = sample_adata[chol_mask]
    
    if len(cholangiocyte_spots) < 2:
        print(f"    ⚠️  Insufficient cholangiocytes ({len(cholangiocyte_spots)}) for zone detection")
        if len(cholangiocyte_spots) == 1:
            # Create single zone for the one cholangiocyte
            if 'spatial' in cholangiocyte_spots.obsm.keys():
                coords = cholangiocyte_spots.obsm['spatial'] * scale_factor
            else:
                coords = np.column_stack([
                    cholangiocyte_spots.obs['x'].values * scale_factor,
                    cholangiocyte_spots.obs['y'].values * scale_factor
                ])
            
            return {
                'zone_0': {
                    'center': coords[0],
                    'cholangiocyte_spots': [cholangiocyte_spots.obs.index[0]],
                    'cholangiocyte_coords': coords,
                    'n_cholangiocytes': 1
                }
            }
        else:
            return {}
    
    # Get spatial coordinates (convert to micrometers)
    if 'spatial' in cholangiocyte_spots.obsm.keys():
        coords = cholangiocyte_spots.obsm['spatial'] * scale_factor
    else:
        coords = np.column_stack([
            cholangiocyte_spots.obs['x'].values * scale_factor,
            cholangiocyte_spots.obs['y'].values * scale_factor
        ])
    
    # Calculate pairwise distances between cholangiocytes
    distance_matrix = cdist(coords, coords, metric='euclidean')
    
    # Create adjacency matrix: connect spots within max_distance_um
    adjacency_matrix = (distance_matrix <= max_distance_um) & (distance_matrix > 0)
    
    # Convert to sparse matrix for connected components analysis
    sparse_adj = csr_matrix(adjacency_matrix)
    
    # Find connected components (cholangiocyte zones)
    n_components, labels = connected_components(
        csgraph=sparse_adj, 
        directed=False, 
        return_labels=True
    )
    
    print(f"    ✅ Identified {n_components} cholangiocyte zones using {max_distance_um}μm threshold")
    
    # Create zone dictionary
    cholangiocyte_zones = {}
    for zone_id in range(n_components):
        zone_mask = labels == zone_id
        zone_coords = coords[zone_mask]
        zone_spots = cholangiocyte_spots.obs.index[zone_mask].tolist()
        
        # Calculate zone center (centroid of cholangiocytes)
        center = np.mean(zone_coords, axis=0)
        
        cholangiocyte_zones[f'zone_{zone_id}'] = {
            'center': center,
            'cholangiocyte_spots': zone_spots,
            'cholangiocyte_coords': zone_coords,
            'n_cholangiocytes': len(zone_spots)
        }
        
        print(f"      Zone {zone_id}: {len(zone_spots)} cholangiocytes")
    
    return cholangiocyte_zones

def analyze_zone_microenvironment_distance_based(sample_adata, sample_abundance, zone_info, 
                                                distance_bins=[(0, 50), (50, 100), (100, 150)],
                                                scale_factor=0.5, zone_id=None, 
                                                external_tracker=None):
    """
    Analyze microenvironment using distance-based approach around cholangiocyte zone
    Modified to use external tracking and include cholangiocyte signature analysis
    
    Parameters:
    -----------
    sample_adata : AnnData
        Single sample spatial data
    sample_abundance : DataFrame
        Cell abundance data for this sample
    zone_info : dict
        Cholangiocyte zone center and information
    distance_bins : list
        Distance bins for analysis, e.g., [(0, 50), (50, 100), (100, 150)]
    scale_factor : float
        Conversion factor from pixels to micrometers
    zone_id : str
        Zone identifier for tracking
    external_tracker : dict
        External dictionary to track assignments across zones
    
    Returns:
    --------
    zone_microenv_data : dict
        Distance-based microenvironment composition including cholangiocyte signatures
    """
    center = zone_info['center']
    
    # Get all spot coordinates (convert to micrometers)
    if 'spatial' in sample_adata.obsm.keys():
        all_coords = sample_adata.obsm['spatial'] * scale_factor
    else:
        all_coords = np.column_stack([
            sample_adata.obs['x'].values * scale_factor,
            sample_adata.obs['y'].values * scale_factor
        ])
    
    # Calculate distances from all spots to zone center
    distances = np.sqrt(np.sum((all_coords - center)**2, axis=1))
    
    # Track assignments in external tracker (don't modify adata yet)
    if external_tracker is not None and zone_id is not None:
        # Initialize zone tracking if first zone
        if 'zone_assignments' not in external_tracker:
            external_tracker['zone_assignments'] = {}
            external_tracker['zone_distances'] = {}
            external_tracker['spot_indices'] = sample_adata.obs.index.tolist()
        
        # Store this zone's assignments and distances
        zone_assignments = np.full(len(sample_adata), 'unassigned', dtype='object')
        
        # Assign cholangiocytes first
        chol_mask = sample_adata.obs['cholangiocyte_enriched'] == True
        zone_assignments[chol_mask] = 'cholangiocyte'
        
        # Exclude cholangiocytes from microenvironment analysis
        non_chol_mask = sample_adata.obs['cholangiocyte_enriched'] == False
        
        max_radius = max([end for _, end in distance_bins]) if distance_bins else 0
        
        # Assign distance bins
        for bin_idx, (min_dist, max_dist) in enumerate(distance_bins):
            bin_mask = (distances >= min_dist) & (distances < max_dist) & non_chol_mask
            bin_label = f'bin_{bin_idx}_{min_dist}-{max_dist}um'
            zone_assignments[bin_mask] = bin_label
        
        # Assign beyond range
        if max_radius > 0:
            beyond_mask = (distances > max_radius) & non_chol_mask
            zone_assignments[beyond_mask] = 'beyond_range'
        
        # Store in external tracker
        external_tracker['zone_assignments'][zone_id] = zone_assignments
        external_tracker['zone_distances'][zone_id] = distances
    
    # Continue with original microenvironment analysis logic
    distance_assignment = np.full(len(sample_adata), 'unassigned', dtype='object')
    
    # Assign cholangiocytes first
    chol_mask = sample_adata.obs['cholangiocyte_enriched'] == True
    distance_assignment[chol_mask] = 'cholangiocyte'
    
    # Exclude cholangiocytes from microenvironment analysis (for distance bins)
    non_chol_mask = sample_adata.obs['cholangiocyte_enriched'] == False
    
    # Analyze each distance bin
    distance_analysis = {}
    total_microenv_spots = 0
    max_radius = max([end for _, end in distance_bins]) if distance_bins else 0
    
    # Define immune and stromal cell signatures
    cell_types = define_celltypes()
    immune_cell_types = cell_types["immune_cell_types"]
    stromal_cell_types = cell_types["stromal_cell_types"]
    cell_type_cols = immune_cell_types + stromal_cell_types

    # ADDED: Analyze cholangiocytes themselves (distance 0)
    chol_spots = sample_adata.obs.index[chol_mask]
    if len(chol_spots) > 0:
        # Get abundance data for cholangiocytes
        chol_abundance = sample_abundance[sample_abundance.index.isin(zone_info["cholangiocyte_spots"])]
        
        if len(chol_abundance) > 0:
            # Calculate signatures for cholangiocytes
            chol_immune_signature = chol_abundance[immune_cell_types].sum(axis=1).mean()
            chol_stromal_signature = chol_abundance[stromal_cell_types].sum(axis=1).mean()
            chol_immune_stromal_ratio = chol_immune_signature / (chol_stromal_signature + 1e-6)
            
            # Calculate individual cell type abundances for cholangiocytes
            chol_cell_abundances = chol_abundance[cell_type_cols].mean().to_dict()
            
            distance_analysis['cholangiocyte_distance_0'] = {
                'bin_index': -1,  # Special index for cholangiocytes
                'min_distance': 0,
                'max_distance': 0,
                'n_spots': len(zone_info["cholangiocyte_spots"]),
                'mean_distance': 0,
                'immune_signature': chol_immune_signature,
                'stromal_signature': chol_stromal_signature,
                'immune_stromal_ratio': chol_immune_stromal_ratio,
                'cell_abundances': chol_cell_abundances
            }
            
            print(f"        Zone {zone_id} cholangiocytes (distance 0): {len(chol_spots)} spots, I/S ratio: {chol_immune_stromal_ratio:.3f}")

    # Analyze regular distance bins (excluding cholangiocytes)
    for bin_idx, (min_dist, max_dist) in enumerate(distance_bins):
        # Define spots in this distance bin (excluding cholangiocytes)
        bin_mask = (distances >= min_dist) & (distances < max_dist) & non_chol_mask
        bin_spots = sample_adata.obs.index[bin_mask]
        
        # Assign distance bin label
        bin_label = f'bin_{bin_idx}_{min_dist}-{max_dist}um'
        distance_assignment[bin_mask] = bin_label
        
        if len(bin_spots) == 0:
            print(f"        Zone {zone_id} bin {min_dist}-{max_dist}μm: 0 spots")
            continue
        
        total_microenv_spots += len(bin_spots)
        
        # Get abundance data for spots in this bin
        bin_abundance = sample_abundance[sample_abundance.index.isin(bin_spots)]
        
        if len(bin_abundance) == 0:
            continue
        
        # Calculate signatures for this distance bin
        immune_signature = bin_abundance[immune_cell_types].sum(axis=1).mean()
        stromal_signature = bin_abundance[stromal_cell_types].sum(axis=1).mean()
        immune_stromal_ratio = immune_signature / (stromal_signature + 1e-6)
        
        # Calculate individual cell type abundances
        cell_abundances = bin_abundance[cell_type_cols].mean().to_dict()
        
        distance_analysis[bin_label] = {
            'bin_index': bin_idx,
            'min_distance': min_dist,
            'max_distance': max_dist,
            'n_spots': len(bin_spots),
            'mean_distance': distances[bin_mask].mean(),
            'immune_signature': immune_signature,
            'stromal_signature': stromal_signature,
            'immune_stromal_ratio': immune_stromal_ratio,
            'cell_abundances': cell_abundances
        }
        
        print(f"        Zone {zone_id} bin {min_dist}-{max_dist}μm: {len(bin_spots)} spots, I/S ratio: {immune_stromal_ratio:.3f}")
    
    # Assign spots beyond maximum range
    if max_radius > 0:
        beyond_mask = (distances > max_radius) & non_chol_mask
        distance_assignment[beyond_mask] = 'beyond_range'
        n_beyond = np.sum(beyond_mask)
        print(f"        Zone {zone_id} spots beyond max range (>{max_radius}μm): {n_beyond} spots")
    
    # Calculate overall zone microenvironment (all distance bins combined, still excluding cholangiocytes)
    if total_microenv_spots > 0:
        overall_mask = (distances <= max_radius) & non_chol_mask
        overall_spots = sample_adata.obs.index[overall_mask]
        
        if len(overall_spots) > 0:
            overall_abundance = sample_abundance[sample_abundance.index.isin(overall_spots)]
            
            overall_immune = overall_abundance[immune_cell_types].sum(axis=1).mean()
            overall_stromal = overall_abundance[stromal_cell_types].sum(axis=1).mean()
            overall_ratio = overall_immune / (overall_stromal + 1e-6)
            overall_cell_abundances = overall_abundance[cell_type_cols].mean().to_dict()
        else:
            overall_immune = overall_stromal = overall_ratio = 0
            overall_cell_abundances = {ct: 0 for ct in cell_type_cols}
    else:
        overall_immune = overall_stromal = overall_ratio = 0
        overall_cell_abundances = {ct: 0 for ct in cell_type_cols}
    
    # Create assignment summary (for this zone only)
    assignment_summary = pd.Series(distance_assignment).value_counts().to_dict()
    print(f"        Zone {zone_id} assignment summary: {assignment_summary}")
    
    zone_microenv_data = {
        'center': center,
        'distance_bins': distance_bins,
        'max_radius': max_radius,
        'total_microenv_spots': total_microenv_spots,
        'distance_analysis': distance_analysis,
        'assignment_summary': assignment_summary,
        'overall_immune_signature': overall_immune,
        'overall_stromal_signature': overall_stromal,
        'overall_immune_stromal_ratio': overall_ratio,
        'overall_cell_abundances': overall_cell_abundances
    }
    
    return zone_microenv_data

def finalize_distance_assignments(sample_adata, external_tracker, strategy='nearest_zone'):
    """
    Finalize distance assignments from external tracker and add to adata
    
    Parameters:
    -----------
    sample_adata : AnnData
        Sample data to modify
    external_tracker : dict
        External tracking dictionary with zone assignments
    strategy : str
        Strategy for resolving conflicts ('nearest_zone', 'first_zone', 'most_restrictive')
    
    Returns:
    --------
    sample_adata : AnnData
        Modified with final distance assignments
    """
    if 'zone_assignments' not in external_tracker:
        print("No zone assignments found in external tracker")
        return sample_adata
    
    n_spots = len(sample_adata)
    final_assignment = np.full(n_spots, 'unassigned', dtype='object')
    final_distance = np.full(n_spots, np.inf)
    nearest_zone_id = np.full(n_spots, 'none', dtype='object')
    
    zone_assignments = external_tracker['zone_assignments']
    zone_distances = external_tracker['zone_distances']
    
    print(f"  Finalizing assignments from {len(zone_assignments)} zones using '{strategy}' strategy...")
    
    if strategy == 'nearest_zone':
        # Assign each spot to the nearest zone's assignment
        for i in range(n_spots):
            min_dist = np.inf
            best_assignment = 'unassigned'
            best_zone = 'none'
            
            for zone_id in zone_assignments.keys():
                zone_dist = zone_distances[zone_id][i]
                if zone_dist < min_dist:
                    min_dist = zone_dist
                    best_assignment = zone_assignments[zone_id][i]
                    best_zone = zone_id
            
            final_assignment[i] = best_assignment
            final_distance[i] = min_dist
            nearest_zone_id[i] = best_zone
    
    elif strategy == 'first_zone':
        # Use first zone's assignment (order matters)
        first_zone = list(zone_assignments.keys())[0]
        final_assignment = zone_assignments[first_zone]
        final_distance = zone_distances[first_zone]
        nearest_zone_id = np.full(n_spots, first_zone, dtype='object')
    
    elif strategy == 'most_restrictive':
        # Use the most restrictive (smallest) distance bin assignment
        bin_order = ['cholangiocyte', 'bin_0_0-50um', 'bin_1_50-100um', 'bin_2_100-150um', 
                    'bin_3_150-200um', 'beyond_range', 'unassigned']
        
        for i in range(n_spots):
            best_assignment = 'unassigned'
            best_distance = np.inf
            best_zone = 'none'
            
            for zone_id in zone_assignments.keys():
                assignment = zone_assignments[zone_id][i]
                distance = zone_distances[zone_id][i]
                
                # Find the most restrictive assignment
                if assignment in bin_order:
                    current_priority = bin_order.index(assignment)
                    if assignment != 'unassigned':
                        if best_assignment == 'unassigned' or bin_order.index(best_assignment) > current_priority:
                            best_assignment = assignment
                            best_distance = distance
                            best_zone = zone_id
            
            final_assignment[i] = best_assignment
            final_distance[i] = best_distance
            nearest_zone_id[i] = best_zone
    
    # Add final assignments to adata
    sample_adata.obs['distance_assignment'] = final_assignment
    sample_adata.obs['distance_to_nearest_zone_center'] = final_distance
    sample_adata.obs['nearest_zone_id'] = nearest_zone_id
    
    # Create summary
    assignment_summary = pd.Series(final_assignment).value_counts().to_dict()
    print(f"  Final assignment summary: {assignment_summary}")
    
    return sample_adata

def analyze_sample_cholangiocyte_zones(sample_id, adata, abundance_df, 
                                     max_distance_um=30, 
                                     distance_bins=[(0, 50), (50, 100), (100, 150)],
                                     scale_factor=0.5,
                                     assignment_strategy='nearest_zone'):
    """
    Analyze cholangiocyte zones for a single sample using graph-based approach
    Modified to use external tracking for distance assignments
    
    Parameters:
    -----------
    sample_id : str
        Sample identifier
    adata : AnnData
        Complete spatial data
    abundance_df : DataFrame
        Cell abundance data
    max_distance_um : float
        Maximum distance to connect cholangiocytes (μm)
    distance_bins : list
        Distance bins for microenvironment analysis
    scale_factor : float
        Pixel to micrometer conversion factor
    assignment_strategy : str
        Strategy for finalizing assignments ('nearest_zone', 'first_zone', 'most_restrictive')
    
    Returns:
    --------
    sample_results : dict
        Zone analysis results for this sample
    """
    print(f"\n🔬 Analyzing cholangiocyte zones in sample: {sample_id}")
    
    # Get sample data
    sample_mask = adata.obs['sample_id'] == sample_id
    sample_adata = adata[sample_mask].copy()
    sample_abundance = abundance_df[abundance_df['sample_id'] == sample_id].copy()
    
    if len(sample_adata) == 0:
        print(f"  ❌ No data found for sample {sample_id}")
        return None
    
    # Reset abundance index to match adata indices
    sample_abundance.index = sample_adata.obs.index
    
    print(f"  Total spots: {len(sample_adata):,}")
    print(f"  Total cholangiocytes: {sample_adata.obs['cholangiocyte_enriched'].sum()}")
    
    # Identify cholangiocyte zones using graph-based approach
    cholangiocyte_zones = identify_cholangiocyte_zones_graph(
        sample_adata, max_distance_um=max_distance_um, scale_factor=scale_factor
    )
    
    if len(cholangiocyte_zones) == 0:
        print(f"  ❌ No cholangiocyte zones identified")
        return None
    
    # Initialize external tracker for distance assignments
    external_tracker = {}
    
    # Analyze microenvironment for each zone using distance-based approach
    zone_results = {}
    for zone_id, zone_info in cholangiocyte_zones.items():
        print(f"    Analyzing {zone_id}...")
        
        # Use modified function with external tracking
        microenv_data = analyze_zone_microenvironment_distance_based(
            sample_adata, sample_abundance, zone_info, 
            distance_bins=distance_bins, scale_factor=scale_factor,
            zone_id=zone_id, external_tracker=external_tracker
        )
        
        if microenv_data is not None:
            # Combine zone and microenvironment data
            zone_results[zone_id] = {
                **zone_info,
                **microenv_data,
                'sample_id': sample_id,
                'rfs_status': sample_adata.obs['RFS_status'].iloc[0]
            }
            print(f"      ✅ Zone microenv: {microenv_data['total_microenv_spots']} spots, overall I/S ratio: {microenv_data['overall_immune_stromal_ratio']:.3f}")
    
    # After all zones processed, finalize distance assignments
    print(f"  Finalizing distance assignments for {len(cholangiocyte_zones)} zones...")
    sample_adata = finalize_distance_assignments(
        sample_adata, external_tracker, strategy=assignment_strategy
    )
    
    return {
        'sample_id': sample_id,
        'rfs_status': sample_adata.obs['RFS_status'].iloc[0],
        'total_spots': len(sample_adata),
        'total_cholangiocytes': sample_adata.obs['cholangiocyte_enriched'].sum(),
        'n_zones': len(zone_results),
        'max_distance_um': max_distance_um,
        'distance_bins': distance_bins,
        'assignment_strategy': assignment_strategy,
        'zones': zone_results,
        'sample_adata': sample_adata,  # Keep for visualization
        'external_tracker': external_tracker  # Keep for debugging if needed
    }

def define_celltypes():
#     immune_cell_types = [

#     'q05cell_abundance_w_sf_CD4T_CXCL13',
#     'q05cell_abundance_w_sf_CD4T_IL7R',
#     'q05cell_abundance_w_sf_Treg_FOXP3',

#     'q05cell_abundance_w_sf_NK_FCGR3A',
#     'q05cell_abundance_w_sf_NK_GZMK',

# ]
    immune_cell_types = [
    'q05cell_abundance_w_sf_CD8T_GZMB',
    'q05cell_abundance_w_sf_CD8T_GZMK', 
    'q05cell_abundance_w_sf_CD8T_MKI67',
    # 'q05cell_abundance_w_sf_CD8T_MT',
    'q05cell_abundance_w_sf_CD8T_STAT1',
    'q05cell_abundance_w_sf_Macro_APOE',
    'q05cell_abundance_w_sf_Macro_FOLR2',
    # 'q05cell_abundance_w_sf_Macro_HSPH1',
    'q05cell_abundance_w_sf_Macro_S100A12',
    # 'q05cell_abundance_w_sf_Macro_SPP1',
    'q05cell_abundance_w_sf_Macro_VCAN',
    ]

    # macro_cell_types = [
    # 'q05cell_abundance_w_sf_Macro_APOE',
    # 'q05cell_abundance_w_sf_Macro_FOLR2',
    # 'q05cell_abundance_w_sf_Macro_HSPH1',
    # 'q05cell_abundance_w_sf_Macro_S100A12',
    # 'q05cell_abundance_w_sf_Macro_SPP1',
    # 'q05cell_abundance_w_sf_Macro_VCAN',
    # ]
    
    # dc_cell_types = [
    # 'q05cell_abundance_w_sf_cDC_CD1C',
    # 'q05cell_abundance_w_sf_cDC_CLEC9A',
    # 'q05cell_abundance_w_sf_pDC_LILRA4'
    # ]

    # Define stromal cell types (fibroblasts, CAFs, etc.)
    stromal_cell_types = [
        # 'q05cell_abundance_w_sf_CAF_CXCL12',
        # 'q05cell_abundance_w_sf_CAF_PDPN',
        'q05cell_abundance_w_sf_CAF_POSTN',
        'q05cell_abundance_w_sf_Fibro_COL3A1',
        # 'q05cell_abundance_w_sf_Fibro_MKI67',
        'q05cell_abundance_w_sf_Fibro_RGS5'
    ]

    return {
        "immune_cell_types":immune_cell_types,
        # "cd8t_cell_types":cd8t_cell_types,
        # "macro_cell_types":macro_cell_types,
        # "dc_cell_types":dc_cell_types,
        "stromal_cell_types":stromal_cell_types,
    }

def classify_zones_pcme(zone_df, method='percentile', distence_prefix=None, use_overall=True):
    """
    Classify cholangiocyte zones as PCME-I, PCME-S, or Intermediate
    
    Parameters:
    -----------
    zone_df : DataFrame
        Zone data with immune/stromal signatures
    method : str
        Classification method ('median', 'percentile', 'kmeans')
    distence_prefix : str
        Distance prefix for specific distance analysis
    use_overall : bool
        Use overall signatures vs distance-specific
    
    Returns:
    --------
    classifications : Series
        PCME classification for each zone (PCME-I, PCME-S, Intermediate)
    thresholds : dict
        Classification thresholds used
    """
    if use_overall:
        ratios = zone_df['immune_stromal_ratio'].values
        immune_vals = zone_df['immune_signature'].values
        stromal_vals = zone_df['stromal_signature'].values
        valid_indices = zone_df.index
    else:
        # Use specific distance bin for classification
        ratio_col = f'{distence_prefix}_immune_stromal_ratio'
        immune_col = f'{distence_prefix}_immune_signature'
        stromal_col = f'{distence_prefix}_stromal_signature'
        
        # Get valid data (non-null)
        valid_mask = (zone_df[ratio_col].notna() & 
                     zone_df[immune_col].notna() & 
                     zone_df[stromal_col].notna())
        
        ratios = zone_df.loc[valid_mask, ratio_col].values
        immune_vals = zone_df.loc[valid_mask, immune_col].values
        stromal_vals = zone_df.loc[valid_mask, stromal_col].values
        valid_indices = zone_df.index[valid_mask]
    
    if len(ratios) == 0:
        return pd.Series(['Intermediate'] * len(zone_df), index=zone_df.index), {}
    
    classifications = np.full(len(zone_df), 'Intermediate', dtype='object')
       
    if method == 'percentile':
        # Use 25th and 75th percentiles for clearer separation
        threshold_low = np.percentile(ratios, 25)
        threshold_high = np.percentile(ratios, 75)
        
        for i, idx in enumerate(valid_indices):
            ratio = ratios[i]
            if ratio >= threshold_high:
                classifications[zone_df.index.get_loc(idx)] = 'PCME-I'
            elif ratio <= threshold_low:
                classifications[zone_df.index.get_loc(idx)] = 'PCME-S'
            else:
                classifications[zone_df.index.get_loc(idx)] = 'Intermediate'
        
        thresholds = {'low': threshold_low, 'high': threshold_high, 'method': 'percentile'}
        
    elif method == 'kmeans':
        # Prepare features for clustering
        features = np.column_stack([immune_vals, stromal_vals, ratios])
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(features)
        
        # K-means clustering with 3 clusters
        kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(features_scaled)
        
        # Determine which cluster is PCME-I, PCME-S, and Intermediate
        cluster_ratios = []
        for i in range(3):
            cluster_mask = cluster_labels == i
            if np.sum(cluster_mask) > 0:
                cluster_ratios.append(ratios[cluster_mask].mean())
            else:
                cluster_ratios.append(0)
        
        # Sort clusters by mean I/S ratio
        sorted_clusters = np.argsort(cluster_ratios)
        pcme_s_cluster = sorted_clusters[0]      # Lowest I/S ratio
        intermediate_cluster = sorted_clusters[1] # Middle I/S ratio  
        pcme_i_cluster = sorted_clusters[2]      # Highest I/S ratio
        
        # Assign classifications
        for i, idx in enumerate(valid_indices):
            label = cluster_labels[i]
            if label == pcme_i_cluster:
                classifications[zone_df.index.get_loc(idx)] = 'PCME-I'
            elif label == pcme_s_cluster:
                classifications[zone_df.index.get_loc(idx)] = 'PCME-S'
            else:  # intermediate_cluster
                classifications[zone_df.index.get_loc(idx)] = 'Intermediate'
        
        thresholds = {
            'pcme_s_mean': cluster_ratios[pcme_s_cluster],
            'intermediate_mean': cluster_ratios[intermediate_cluster],
            'pcme_i_mean': cluster_ratios[pcme_i_cluster],
            'method': 'kmeans_3_clusters'
        }
    
    elif method == 'adaptive':
        # Adaptive method based on data distribution
        # Use IQR-based thresholds for better handling of outliers
        q25 = np.percentile(ratios, 25)
        q75 = np.percentile(ratios, 75)
        iqr = q75 - q25
        
        # Define thresholds based on IQR
        threshold_low = q25 - 0.5 * iqr
        threshold_high = q75 + 0.5 * iqr
        
        # Ensure we don't have too extreme thresholds
        threshold_low = max(threshold_low, np.percentile(ratios, 10))
        threshold_high = min(threshold_high, np.percentile(ratios, 90))
        
        for i, idx in enumerate(valid_indices):
            ratio = ratios[i]
            if ratio >= threshold_high:
                classifications[zone_df.index.get_loc(idx)] = 'PCME-I'
            elif ratio <= threshold_low:
                classifications[zone_df.index.get_loc(idx)] = 'PCME-S'
            else:
                classifications[zone_df.index.get_loc(idx)] = 'Intermediate'
        
        thresholds = {'low': threshold_low, 'high': threshold_high, 'method': 'adaptive_iqr'}
    
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'median', 'percentile', 'kmeans', 'adaptive'")
    
    # Print classification summary
    class_counts = pd.Series(classifications).value_counts()
    print(f"\nPCME Classification Summary ({method} method):")
    for class_name in ['PCME-I', 'Intermediate', 'PCME-S']:
        count = class_counts.get(class_name, 0)
        percentage = (count / len(classifications)) * 100
        print(f"  {class_name}: {count} zones ({percentage:.1f}%)")
    
    if method in ['median', 'percentile', 'adaptive']:
        print(f"  Thresholds: Low={thresholds['low']:.3f}, High={thresholds['high']:.3f}")
    elif method == 'kmeans':
        print(f"  Cluster means: PCME-S={thresholds['pcme_s_mean']:.3f}, "
              f"Intermediate={thresholds['intermediate_mean']:.3f}, "
              f"PCME-I={thresholds['pcme_i_mean']:.3f}")
    
    return pd.Series(classifications, index=zone_df.index), thresholds

def calculate_sample_pcme_metrics_zones(zone_df):
    """
    Calculate sample-level PCME metrics from zone classifications
    """
    sample_metrics = []
    
    for sample_id in zone_df['sample_id'].unique():
        sample_zones = zone_df[zone_df['sample_id'] == sample_id]
        total_zones = len(sample_zones)

        idx_ = sample_zones["pcme_classification"] == "Intermediate"
        sample_zones = sample_zones.loc[~idx_]
        
        if len(sample_zones) == 0:
            continue
            
        rfs_status = sample_zones['rfs_status'].iloc[0]
        
        # Count PCME types
        pcme_counts = sample_zones['pcme_classification'].value_counts()
        pcme_i_count = pcme_counts.get('PCME-I', 0)
        pcme_s_count = pcme_counts.get('PCME-S', 0)
        
        # Calculate percentages
        pcme_i_pct = (pcme_i_count / total_zones) * 100
        pcme_s_pct = (pcme_s_count / total_zones) * 100
        
        # Determine predominant type
        if pcme_i_count > pcme_s_count:
            predominant_type = 'PCME-I'
            dominance_score = pcme_i_pct
        elif pcme_s_count > pcme_i_count:
            predominant_type = 'PCME-S'
            dominance_score = pcme_s_pct
        else:
            predominant_type = 'Balanced'
            dominance_score = 50.0
        
        # Calculate zone-weighted averages
        weighted_immune = sample_zones['immune_signature'].mean()
        weighted_stromal = sample_zones['stromal_signature'].mean()
        weighted_ratio = sample_zones['immune_stromal_ratio'].mean()
        
        sample_metrics.append({
            'sample_id': sample_id,
            'rfs_status': rfs_status,
            'total_zones': total_zones,
            'pcme_i_count': pcme_i_count,
            'pcme_s_count': pcme_s_count,
            'pcme_i_percentage': pcme_i_pct,
            'pcme_s_percentage': pcme_s_pct,
            'predominant_type': predominant_type,
            'dominance_score': dominance_score,
            'weighted_immune_signature': weighted_immune,
            'weighted_stromal_signature': weighted_stromal,
            'weighted_is_ratio': weighted_ratio
        })
    
    return pd.DataFrame(sample_metrics)

#%% PCME expression analysis

def extract_pcme_region_expression(adata, zone_df, sample_results, max_radius=150):
    """
    Extract gene expression data from PCME-I and PCME-S regions
    
    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data
    zone_df : DataFrame
        Zone classifications
    sample_results : dict
        Sample analysis results
    max_radius : float
        Maximum radius (μm) to define PCME regions
    
    Returns:
    --------
    pcme_expression_data : dict
        Gene expression data for PCME regions
    """
    print(f"\n🧬 Extracting gene expression from PCME regions...")
    
    pcme_expression_data = {
        'PCME-I': {'spots': [], 'expression': []},
        'PCME-S': {'spots': [], 'expression': []}
    }
    
    # Scale factor for coordinate conversion
    scale_factor = 0.5
    
    for _, zone_row in zone_df.iterrows():
        pcme_type = zone_row['pcme_classification']
        if pcme_type not in ['PCME-I', 'PCME-S']:
            continue
            
        sample_id = zone_row['sample_id']
        zone_id = 'zone_' + str(zone_row['zone_id'].split('_')[-1]) 
        
        # Get sample data
        sample_mask = adata.obs['sample_id'] == sample_id
        sample_adata = adata[sample_mask]
        
        if len(sample_adata) == 0:
            continue
            
        # Get zone center from sample results
        if sample_id in sample_results and zone_id in sample_results[sample_id]['zones']:
            zone_info = sample_results[sample_id]['zones'][zone_id]
            zone_center = zone_info['center']
            
            # Get all spot coordinates (convert to micrometers)
            if 'spatial' in sample_adata.obsm.keys():
                all_coords = sample_adata.obsm['spatial'] * scale_factor
            else:
                all_coords = np.column_stack([
                    sample_adata.obs['x'].values * scale_factor,
                    sample_adata.obs['y'].values * scale_factor
                ])
            
            # Calculate distances to zone center
            distances = np.sqrt(np.sum((all_coords - zone_center)**2, axis=1))
            
            # Define PCME region (within max_radius, excluding cholangiocytes)
            # region_mask = (distances <= max_radius) & (sample_adata.obs['cholangiocyte_enriched'] == False)
            region_mask = (distances <= max_radius)
            region_spots = sample_adata.obs.index[region_mask]
            
            if len(region_spots) > 0:
                # Get expression data for this region
                region_expression = adata[region_spots].X.toarray()
                
                pcme_expression_data[pcme_type]['spots'].extend(region_spots.tolist())
                pcme_expression_data[pcme_type]['expression'].append(region_expression)
    
    # Combine expression data
    for pcme_type in pcme_expression_data:
        if pcme_expression_data[pcme_type]['expression']:
            combined_expression = np.vstack(pcme_expression_data[pcme_type]['expression'])
            pcme_expression_data[pcme_type]['combined_expression'] = combined_expression
            print(f"  {pcme_type}: {combined_expression.shape[0]:,} spots × {combined_expression.shape[1]:,} genes")
        else:
            print(f"  {pcme_type}: No expression data")
    
    return pcme_expression_data


def perform_pcme_differential_analysis(adata, pcme_expression_data, min_spots=10):
    """
    Perform differential gene expression analysis between PCME-I and PCME-S
    
    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data (for gene names)
    pcme_expression_data : dict
        Expression data for PCME regions
    min_spots : int
        Minimum number of spots required for analysis
    
    Returns:
    --------
    deg_results : DataFrame
        Differential expression results
    """
    print(f"\n📊 Performing PCME differential expression analysis...")
    
    # Check if we have enough data
    pcme_i_expr = pcme_expression_data['PCME-I'].get('combined_expression')
    pcme_s_expr = pcme_expression_data['PCME-S'].get('combined_expression')
    
    if pcme_i_expr is None or pcme_s_expr is None:
        print("❌ Insufficient expression data for analysis")
        return None
        
    if len(pcme_i_expr) < min_spots or len(pcme_s_expr) < min_spots:
        print(f"❌ Insufficient spots: PCME-I={len(pcme_i_expr)}, PCME-S={len(pcme_s_expr)}")
        return None
    
    print(f"Analyzing PCME-I ({len(pcme_i_expr)} spots) vs PCME-S ({len(pcme_s_expr)} spots)")
    
    # Perform differential expression for each gene
    deg_results = []
    gene_names = adata.var_names
    
    for i, gene in enumerate(gene_names):
        pcme_i_values = pcme_i_expr[:, i]
        pcme_s_values = pcme_s_expr[:, i]
        
        # Basic statistics
        mean_pcme_i = np.mean(pcme_i_values)
        mean_pcme_s = np.mean(pcme_s_values)
        
        # Fold change
        fold_change = mean_pcme_i / (mean_pcme_s + 1e-10)
        log2_fc = np.log2(fold_change) if fold_change > 0 else 0
        
        # Statistical test
        try:
            if np.std(pcme_i_values) == 0 and np.std(pcme_s_values) == 0:
                if mean_pcme_i == mean_pcme_s:
                    p_value = 1.0
                else:
                    p_value = 0.0
            else:
                _, p_value = stats.ttest_ind(pcme_i_values, pcme_s_values)
        except:
            p_value = 1.0
        
        # Expression frequency
        pct_pcme_i = (pcme_i_values > 0).mean() * 100
        pct_pcme_s = (pcme_s_values > 0).mean() * 100
        
        deg_results.append({
            'gene': gene,
            'mean_pcme_i': mean_pcme_i,
            'mean_pcme_s': mean_pcme_s,
            'fold_change': fold_change,
            'log2_fc': log2_fc,
            'p_value': p_value,
            'pct_pcme_i': pct_pcme_i,
            'pct_pcme_s': pct_pcme_s
        })
    
    deg_df = pd.DataFrame(deg_results)
    
    # Multiple testing correction
    _, p_adj, _, _ = multipletests(deg_df['p_value'], method='fdr_bh')
    deg_df['p_adj'] = p_adj
    
    # Sort by significance
    deg_df = deg_df.sort_values('p_adj')
    
    print(f"Differential expression analysis complete:")
    print(f"  Total genes tested: {len(deg_df)}")
    print(f"  Significant genes (padj < 0.05): {(deg_df['p_adj'] < 0.05).sum()}")
    print(f"  Upregulated in PCME-I (FC > 1.5): {((deg_df['fold_change'] > 1.5) & (deg_df['p_adj'] < 0.05)).sum()}")
    print(f"  Downregulated in PCME-I (FC < 0.67): {((deg_df['fold_change'] < 0.67) & (deg_df['p_adj'] < 0.05)).sum()}")
    
    return deg_df

def perform_pathway_enrichment_analysis(deg_df, significance_threshold=0.05):
    """
    Perform pathway enrichment analysis on PCME differential genes
    
    Parameters:
    -----------
    deg_df : DataFrame
        Differential expression results
    significance_threshold : float
        P-value threshold for significance
    
    Returns:
    --------
    enrichment_results : dict
        Pathway enrichment results
    """
    print(f"\n🛤️  Performing pathway enrichment analysis...")
    
    # Load pathway databases
    pathway_databases = load_pathway_databases()
    
    # Get significant upregulated and downregulated genes
    significant_genes = deg_df[deg_df['p_adj'] < significance_threshold]
    upregulated_genes = significant_genes[significant_genes['log2_fc'] > 0]['gene'].tolist()
    downregulated_genes = significant_genes[significant_genes['log2_fc'] < 0]['gene'].tolist()
    background_genes = deg_df['gene'].tolist()
    
    print(f"Genes for enrichment analysis:")
    print(f"  Upregulated (PCME-I enriched): {len(upregulated_genes)}")
    print(f"  Downregulated (PCME-S enriched): {len(downregulated_genes)}")
    print(f"  Background: {len(background_genes)}")
    
    enrichment_results = {}
    
    # Enrichment for upregulated genes (PCME-I enriched)
    if len(upregulated_genes) > 0:
        enrichment_up = perform_pathway_enrichment(
            upregulated_genes, pathway_databases, background_genes, p_threshold=0.05
        )
        enrichment_results['PCME-I_enriched'] = enrichment_up
        print(f"  PCME-I enriched pathways: {len(enrichment_up)}")
    
    # Enrichment for downregulated genes (PCME-S enriched)
    if len(downregulated_genes) > 0:
        enrichment_down = perform_pathway_enrichment(
            downregulated_genes, pathway_databases, background_genes, p_threshold=0.05
        )
        enrichment_results['PCME-S_enriched'] = enrichment_down
        print(f"  PCME-S enriched pathways: {len(enrichment_down)}")
    
    return enrichment_results

def create_volcano_plot(deg_results, fig_size=(10, 8)):
    """
    Create an enhanced volcano plot for PCME differential expression
    
    Parameters:
    -----------
    deg_results : DataFrame
        Differential expression results
    fig_size : tuple
        Figure size (width, height)
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis
    """
    # Prepare data
    x_vals = deg_results['log2_fc'].fillna(0)
    y_vals = -np.log10(deg_results['p_adj'].fillna(1))
    
    # Define significance thresholds
    fc_threshold = 0.5  # log2 fold change threshold
    pval_threshold = 0.05  # p-value threshold
    y_threshold = -np.log10(pval_threshold)
    
    # Color coding
    colors = []
    labels = []
    sizes = []
    
    deg_results = pd.DataFrame(deg_results.copy()) 

    for _, row in deg_results.iterrows():
        fc = row['log2_fc']
        pval = row['p_adj']
        
        if pval < pval_threshold and abs(fc) > fc_threshold:
            if fc > fc_threshold:
                colors.append('#E31A1C')  # Red for PCME-I enriched
                labels.append('PCME-I enriched')
                sizes.append(30)
            else:
                colors.append('#1F78B4')  # Blue for PCME-S enriched  
                labels.append('PCME-S enriched')
                sizes.append(30)
        elif pval < pval_threshold:
            colors.append('#FF7F00')  # Orange for significant but small FC
            labels.append('Significant')
            sizes.append(20)
        else:
            colors.append('#CCCCCC')  # Gray for non-significant
            labels.append('Non-significant')
            sizes.append(10)
    
    # Create figure
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Create scatter plot
    scatter = ax.scatter(x_vals, y_vals, c=colors, s=sizes, alpha=0.7, 
                        edgecolors='black', linewidth=0.3)
    
    # Add threshold lines
    ax.axhline(y=y_threshold, color='black', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=fc_threshold, color='black', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=-fc_threshold, color='black', linestyle='--', alpha=0.5, linewidth=1)
    
    # Label top genes
    significant_genes = deg_results[
        (deg_results['p_adj'] < pval_threshold) & 
        (abs(deg_results['log2_fc']) > fc_threshold)
    ].head(15)  # Top 15 most significant
    
    for _, gene in significant_genes.iterrows():
        if abs(gene['log2_fc']) > 1.0:  # Only label genes with substantial fold change
            ax.annotate(gene['gene'], 
                       (gene['log2_fc'], -np.log10(gene['p_adj'])),
                       xytext=(5, 5), textcoords='offset points', 
                       fontsize=8, alpha=0.8, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
    
    # Customize plot
    ax.set_xlabel('Log₂ Fold Change (PCME-I vs PCME-S)', fontsize=14, fontweight='bold')
    ax.set_ylabel('-Log₁₀ Adjusted P-value', fontsize=14, fontweight='bold')
    ax.set_title('PCME Differential Gene Expression', fontsize=16, fontweight='bold', pad=20)
    
    # Add statistics text box
    n_total = len(deg_results)
    n_sig = (deg_results['p_adj'] < pval_threshold).sum()
    n_up = ((deg_results['p_adj'] < pval_threshold) & (deg_results['log2_fc'] > fc_threshold)).sum()
    n_down = ((deg_results['p_adj'] < pval_threshold) & (deg_results['log2_fc'] < -fc_threshold)).sum()
    
    stats_text = f"Total genes: {n_total:,}\n" \
                f"Significant: {n_sig:,} ({n_sig/n_total*100:.1f}%)\n" \
                f"PCME-I enriched: {n_up:,}\n" \
                f"PCME-S enriched: {n_down:,}"
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
           verticalalignment='top', fontsize=10, fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.9))
    
    # Create custom legend
    legend_elements = [
        plt.scatter([], [], c='#E31A1C', s=50, label='PCME-I enriched (significant)'),
        plt.scatter([], [], c='#1F78B4', s=50, label='PCME-S enriched (significant)'),
        plt.scatter([], [], c='#FF7F00', s=30, label='Significant (small FC)'),
        plt.scatter([], [], c='#CCCCCC', s=20, label='Non-significant')
    ]
    ax.legend(handles=legend_elements, loc='upper right', framealpha=0.9, 
             edgecolor='black', fontsize=10)
    
    # Grid and styling
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Set reasonable axis limits
    max_y = np.percentile(y_vals,95) * 1.1
    max_x = np.percentile(x_vals,95) * 2
    ax.set_xlim(-max_x, max_x)
    ax.set_ylim(0, max_y)
    
    plt.tight_layout()
    return fig, ax

def create_pathway_enrichment_plots(enrichment_results, fig_size=(16, 10)):
    """
    Create side-by-side pathway enrichment bar plots
    
    Parameters:
    -----------
    enrichment_results : dict
        Dictionary containing enrichment results for up and down genes
    fig_size : tuple
        Figure size (width, height)
    
    Returns:
    --------
    fig : matplotlib figure
    """
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=fig_size)
    
    # Colors for different pathway categories
    pathway_colors = {
        'IMMUNE': '#E31A1C',      # Red
        'METAB': '#1F78B4',       # Blue  
        'ECM': '#33A02C',         # Green
        'SIGNAL': '#FF7F00',      # Orange
        'CELL_CYCLE': '#6A3D9A',  # Purple
        'STRESS': '#FB9A99',      # Light pink
        'OTHER': '#CCCCCC'        # Gray
    }
    
    def categorize_pathway(pathway_name):
        """Categorize pathway based on name keywords"""
        pathway_upper = pathway_name.upper()
        
        if any(keyword in pathway_upper for keyword in ['IMMUNE', 'T_CELL', 'B_CELL', 'INTERFERON', 'CYTOKINE']):
            return 'IMMUNE'
        elif any(keyword in pathway_upper for keyword in ['METABOLISM', 'GLYCOL', 'OXPHOS', 'FATTY_ACID']):
            return 'METAB'
        elif any(keyword in pathway_upper for keyword in ['ECM', 'COLLAGEN', 'MATRIX', 'ADHESION']):
            return 'ECM'
        elif any(keyword in pathway_upper for keyword in ['SIGNALING', 'PATHWAY', 'CASCADE']):
            return 'SIGNAL'
        elif any(keyword in pathway_upper for keyword in ['CELL_CYCLE', 'PROLIFERATION', 'DIVISION']):
            return 'CELL_CYCLE'
        elif any(keyword in pathway_upper for keyword in ['STRESS', 'HYPOXIA', 'RESPONSE']):
            return 'STRESS'
        else:
            return 'OTHER'
    
    # Plot PCME-I enriched pathways (upregulated)
    if 'PCME-I_enriched' in enrichment_results:
        pcme_i_data = enrichment_results['PCME-I_enriched'].head(15)  # Top 15
        
        if len(pcme_i_data) > 0:
            # Sort by enrichment ratio
            pcme_i_data = pcme_i_data.sort_values('enrichment_ratio', ascending=True)
            
            # Assign colors based on pathway category
            colors_i = [pathway_colors[categorize_pathway(pathway)] for pathway in pcme_i_data['pathway']]
            
            # Create horizontal bar plot
            y_pos = range(len(pcme_i_data))
            bars1 = ax1.barh(y_pos, pcme_i_data['enrichment_ratio'], 
                           color=colors_i, alpha=0.8, edgecolor='black', linewidth=0.5)
            
            # Clean pathway names for display
            clean_names = []
            for pathway in pcme_i_data['pathway']:
                # Remove common prefixes and make readable
                clean_name = pathway.replace('HALLMARK_', '').replace('KEGG_', '').replace('REACTOME_', '')
                clean_name = clean_name.replace('_', ' ').title()
                
                # Truncate long names
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
                
                clean_names.append(clean_name)
            
            ax1.set_yticks(y_pos)
            ax1.set_yticklabels(clean_names, fontsize=10, fontweight='bold')
            ax1.set_xlabel('Enrichment Ratio', fontsize=12, fontweight='bold')
            ax1.set_title('PCME-I Enriched Pathways\n(Immune-Infiltrated)', 
                         fontsize=14, fontweight='bold', color='#E31A1C')
            
            # Add p-value annotations
            for i, (_, row) in enumerate(pcme_i_data.iterrows()):
                p_val = row['p_adjusted']
                if p_val < 0.001:
                    p_text = '***'
                elif p_val < 0.01:
                    p_text = '**'
                elif p_val < 0.05:
                    p_text = '*'
                else:
                    p_text = f'{p_val:.2e}'
                
                ax1.text(row['enrichment_ratio'] + 0.1, i, p_text, 
                        va='center', fontsize=8, fontweight='bold')
            
            ax1.grid(axis='x', alpha=0.3, linestyle='--')
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
    
    # Plot PCME-S enriched pathways (downregulated)
    if 'PCME-S_enriched' in enrichment_results:
        pcme_s_data = enrichment_results['PCME-S_enriched'].head(15)  # Top 15
        
        if len(pcme_s_data) > 0:
            # Sort by enrichment ratio
            pcme_s_data = pcme_s_data.sort_values('enrichment_ratio', ascending=True)
            
            # Assign colors based on pathway category
            colors_s = [pathway_colors[categorize_pathway(pathway)] for pathway in pcme_s_data['pathway']]
            
            # Create horizontal bar plot
            y_pos = range(len(pcme_s_data))
            bars2 = ax2.barh(y_pos, pcme_s_data['enrichment_ratio'], 
                           color=colors_s, alpha=0.8, edgecolor='black', linewidth=0.5)
            
            # Clean pathway names for display
            clean_names = []
            for pathway in pcme_s_data['pathway']:
                clean_name = pathway.replace('HALLMARK_', '').replace('KEGG_', '').replace('REACTOME_', '')
                clean_name = clean_name.replace('_', ' ').title()
                
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
                
                clean_names.append(clean_name)
            
            ax2.set_yticks(y_pos)
            ax2.set_yticklabels(clean_names, fontsize=10, fontweight='bold')
            ax2.set_xlabel('Enrichment Ratio', fontsize=12, fontweight='bold')
            ax2.set_title('PCME-S Enriched Pathways\n(Stromal-Fibrotic)', 
                         fontsize=14, fontweight='bold', color='#1F78B4')
            
            # Add p-value annotations
            for i, (_, row) in enumerate(pcme_s_data.iterrows()):
                p_val = row['p_adjusted']
                if p_val < 0.001:
                    p_text = '***'
                elif p_val < 0.01:
                    p_text = '**'
                elif p_val < 0.05:
                    p_text = '*'
                else:
                    p_text = f'{p_val:.2e}'
                
                ax2.text(row['enrichment_ratio'] + 0.1, i, p_text, 
                        va='center', fontsize=8, fontweight='bold')
            
            ax2.grid(axis='x', alpha=0.3, linestyle='--')
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
    
    # Add pathway category legend
    legend_elements = [mpatches.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.8, 
                                         edgecolor='black', label=category)
                      for category, color in pathway_colors.items() if category != 'OTHER']
    
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.02), 
              ncol=len(legend_elements), framealpha=0.9, edgecolor='black', fontsize=10)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # Make room for legend
    
    return fig
