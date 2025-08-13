import pandas as pd
import numpy as np
from projects.scRNA_pipeline.packages.MetaTiME.metatime.mecmapper import scale
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

import warnings
warnings.filterwarnings('ignore')

from scipy import stats
from scipy.stats import hypergeom
from scipy.spatial.distance import cdist
from statsmodels.stats.multitest import multipletests
from gseapy import Msigdb
import gseapy as gp

#%% Visualize Data
def plot_sample_overview(metadata_df, figsize=(15, 10)):
    """
    Create sample overview plots
    
    Parameters:
    -----------
    metadata_df : DataFrame
        Sample metadata
    figsize : tuple
        Figure size
    """
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    fig.suptitle('Sample Overview - Step 1.1 Quality Control', fontsize=16, fontweight='bold')
    
    # 1. Sample distribution by group
    ax1 = axes[0, 0]
    group_counts = metadata_df['RFS_status'].value_counts()
    colors = sns.color_palette("Set2", len(group_counts))
    ax1.pie(group_counts.values, labels=group_counts.index, autopct='%1.0f%%', 
            colors=colors, startangle=90)
    ax1.set_title('Sample Distribution by Group')
    
    # 2. Number of spots per sample
    ax2 = axes[0, 1]
    sns.barplot(data=metadata_df, x='sample_id', y='n_spots', 
                hue='RFS_status', ax=ax2)
    ax2.set_title('Number of Spots per Sample')
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 3. Number of genes per sample
    ax3 = axes[0, 2]
    sns.barplot(data=metadata_df, x='sample_id', y='n_genes', 
                hue='RFS_status', ax=ax3)
    ax3.set_title('Number of Genes per Sample')
    ax3.tick_params(axis='x', rotation=45)
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 4. Spots distribution by group
    ax4 = axes[1, 0]
    sns.boxplot(data=metadata_df, x='RFS_status', y='n_spots', ax=ax4)
    ax4.set_title('Spots Distribution by Group')
    ax4.tick_params(axis='x', rotation=45)
    
    # 5. Sample composition table
    ax5 = axes[1, 1]
    ax5.axis('tight')
    ax5.axis('off')
    
    # Create summary table
    summary_table = metadata_df.groupby('RFS_status').agg({
        'sample_id': 'count',
        'n_spots': ['min', 'max', 'mean'],
        'n_genes': ['min', 'max', 'mean']
    }).round(0)
    
    # Flatten column names
    summary_table.columns = ['_'.join(col).strip() for col in summary_table.columns.values]
    summary_table = summary_table.reset_index()
    
    table = ax5.table(cellText=summary_table.values,
                     colLabels=summary_table.columns,
                     cellLoc='center',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    ax5.set_title('Summary Statistics by Group')
    
    # 6. Group comparison
    ax6 = axes[1, 2]
    group_comparison = metadata_df.groupby('RFS_status')['n_spots'].mean()
    bars = ax6.bar(group_comparison.index, group_comparison.values, 
                   color=colors[:len(group_comparison)])
    ax6.set_title('Average Spots per Group')
    ax6.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax6.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom')
    
    plt.tight_layout()
    return fig

def clean_cell_type_names(cell_type_cols):
    """Clean cell type names by removing prefix"""
    clean_names = []
    for col in cell_type_cols:
        clean_name = col.replace('q05cell_abundance_w_sf_', '')
        clean_names.append(clean_name)
    return clean_names

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
def calculate_spatial_distances(sample_adata, scale_factor):
    """
    Calculate distances between all spots in a sample
    
    Parameters:
    -----------
    sample_adata : AnnData
        Single sample spatial data
    
    Returns:
    --------
    distance_matrix : ndarray
        Distance matrix between all spots
    coordinates : ndarray
        X,Y coordinates of spots
    """
    # Get spatial coordinates
    if 'spatial' in sample_adata.obsm.keys():
        coordinates = sample_adata.obsm['spatial']
    else:
        # Use x,y coordinates from obs
        coordinates = np.column_stack([
            sample_adata.obs['x'].values, 
            sample_adata.obs['y'].values
        ])
    
    # Calculate pairwise distances
    distance_matrix = cdist(coordinates, coordinates, metric='euclidean') * scale_factor
    
    return distance_matrix, coordinates

def assign_distance_zones(sample_adata, distance_matrix, distance_bins):
    """
    Assign each spot to distance zones relative to cholangiocytes
    
    Parameters:
    -----------
    sample_adata : AnnData
        Single sample spatial data
    distance_matrix : ndarray
        Distance matrix between spots (n_spots x n_spots)
    distance_bins : list
        List of (start, end) distance tuples, e.g., [(0, 50), (50, 100), (100, 200)]
    
    Returns:
    --------
    zone_assignments : ndarray
        Zone assignment for each spot (-2: cholangiocyte, -1: unassigned, 0+: zone index)
    min_distances_to_chol : ndarray
        Minimum distance from each spot to nearest cholangiocyte
    """
    # Input validation
    n_spots = len(sample_adata)
    
    # Validate distance matrix shape
    if distance_matrix.shape != (n_spots, n_spots):
        raise ValueError(f"Distance matrix shape {distance_matrix.shape} doesn't match n_spots {n_spots}")
    
    # Validate distance_bins format
    if not isinstance(distance_bins, (list, tuple)) or len(distance_bins) == 0:
        raise ValueError("distance_bins must be a non-empty list of (min_dist, max_dist) tuples")
    
    for i, bin_range in enumerate(distance_bins):
        if not isinstance(bin_range, (list, tuple)) or len(bin_range) != 2:
            raise ValueError(f"distance_bins[{i}] must be a (min_dist, max_dist) tuple")
        if bin_range[0] >= bin_range[1]:
            raise ValueError(f"distance_bins[{i}]: min_dist must be < max_dist")
    
    # Initialize zone assignments
    zone_assignments = np.full(n_spots, -1, dtype=int)  # -1 means not assigned
    
    # Find cholangiocyte spots (handle both boolean and binary formats)
    chol_column = sample_adata.obs['cholangiocyte_enriched']
    if chol_column.dtype == bool:
        cholangiocyte_mask = chol_column.values
    else:
        # Assume binary (0/1) or numeric, convert to boolean
        cholangiocyte_mask = (chol_column > 0).values
    
    cholangiocyte_indices = np.where(cholangiocyte_mask)[0]
    
    # Handle case with no cholangiocytes
    if len(cholangiocyte_indices) == 0:
        print(f"  ⚠️  No cholangiocytes found in sample")
        min_distances_to_chol = np.full(n_spots, np.inf)  # Infinite distance
        return zone_assignments, min_distances_to_chol
    
    print(f"  ✅ Found {len(cholangiocyte_indices)} cholangiocyte spots")
    
    # For each spot, find minimum distance to any cholangiocyte
    min_distances_to_chol = np.min(distance_matrix[:, cholangiocyte_indices], axis=1)
    
    # Assign zones based on distance bins
    for zone_idx, (min_dist, max_dist) in enumerate(distance_bins):
        # Use >= for min and < for max to avoid overlaps
        mask = (min_distances_to_chol >= min_dist) & (min_distances_to_chol < max_dist)
        zone_assignments[mask] = zone_idx
        # n_assigned = np.sum(mask)
        # print(f"  Zone {zone_idx} ({min_dist}-{max_dist}): {n_assigned} spots")
    
    # Cholangiocytes themselves get special zone -2
    zone_assignments[cholangiocyte_indices] = -2
    
    # Report assignment statistics
    unique_zones, counts = np.unique(zone_assignments, return_counts=True)
    zone_stats = dict(zip(unique_zones, counts))
    
    print(f"  Zone assignment summary:")
    for zone, count in zone_stats.items():
        if zone == -2:
            print(f"    Cholangiocytes: {count} spots")
        elif zone == -1:
            print(f"    Unassigned: {count} spots")
        else:
            print(f"    Zone {zone}: {count} spots")
    
    return zone_assignments, min_distances_to_chol

def analyze_sample_spatial_zones(sample_id, adata, abundance_df, distance_bins, scale_factor = 0.22):
    """
    Analyze spatial zones for a single sample
    
    Parameters:
    -----------
    sample_id : str
        Sample identifier
    adata : AnnData
        Complete spatial data
    abundance_df : DataFrame
        Cell abundance data
    distance_bins : list
        Distance bin definitions
    
    Returns:
    --------
    sample_results : dict
        Results for this sample
    """
    print(f"\n📊 Analyzing sample: {sample_id}")
    
    # Get sample data
    sample_mask = adata.obs['sample_id'] == sample_id
    sample_adata = adata[sample_mask].copy()
    sample_abundance = abundance_df[abundance_df['sample_id'] == sample_id].copy()
    
    if len(sample_adata) == 0:
        print(f"  ❌ No data found for sample {sample_id}")
        return None
    
    print(f"  Total spots: {len(sample_adata):,}")
    print(f"  Cholangiocytes: {sample_adata.obs['cholangiocyte_enriched'].sum()}")
    
    # Calculate distances
    distance_matrix, coordinates = calculate_spatial_distances(sample_adata, scale_factor)
    
    # Assign distance zones
    zone_assignments, min_distances = assign_distance_zones(
        sample_adata, distance_matrix, distance_bins
    )
    
    # Add zone assignments to sample data
    sample_adata.obs['distance_zone'] = zone_assignments
    sample_adata.obs['min_distance_to_cholangiocyte'] = min_distances
    
    # Analyze zone composition
    zone_composition = {}
    
    for zone_idx, (min_dist, max_dist) in enumerate(distance_bins):
        zone_mask = zone_assignments == zone_idx
        zone_count = zone_mask.sum()
        
        if zone_count > 0:
            # Get abundance data for spots in this zone
            zone_spots = sample_abundance.iloc[zone_mask]
            
            # Calculate mean abundance for each cell type in this zone
            zone_cell_abundances = zone_spots[cell_type_cols].mean()
            
            zone_composition[f'Zone_{zone_idx}_{min_dist}-{max_dist}um'] = {
                'spot_count': zone_count,
                'percentage': (zone_count / len(sample_adata)) * 100,
                'cell_abundances': zone_cell_abundances.to_dict(),
                'mean_distance': min_distances[zone_mask].mean(),
                'std_distance': min_distances[zone_mask].std()
            }
            
            print(f"  Zone {zone_idx} ({min_dist}-{max_dist}μm): {zone_count} spots ({zone_count/len(sample_adata)*100:.1f}%)")
    
    # Cholangiocyte zone
    chol_mask = zone_assignments == -2
    chol_count = chol_mask.sum()
    if chol_count > 0:
        chol_spots = sample_abundance.iloc[chol_mask]
        zone_composition['Cholangiocytes'] = {
            'spot_count': chol_count,
            'percentage': (chol_count / len(sample_adata)) * 100,
            'cell_abundances': chol_spots[cell_type_cols].mean().to_dict(),
            'mean_distance': 0.0,
            'std_distance': 0.0
        }
        print(f"  Cholangiocytes: {chol_count} spots ({chol_count/len(sample_adata)*100:.1f}%)")
    
    sample_results = {
        'sample_id': sample_id,
        'rfs_status': sample_adata.obs['RFS_status'].iloc[0],
        'total_spots': len(sample_adata),
        'cholangiocyte_count': chol_count,
        'zone_composition': zone_composition,
        'coordinates': coordinates,
        'zone_assignments': zone_assignments,
        'min_distances': min_distances
    }
    
    return sample_results
