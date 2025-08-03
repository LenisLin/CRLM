import gc
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc

from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches


def explore_clinical_data(adata, variables=None, show_counts=True, show_percentages=True):
    """
    Explore clinical variables in AnnData object
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    variables : list, optional
        List of variables to explore. If None, explores all categorical variables
    show_counts : bool
        Whether to show raw counts
    show_percentages : bool
        Whether to show percentages
    
    Returns:
    --------
    dict : Dictionary containing summary statistics for each variable
    """
    
    if variables is None:
        # Auto-detect categorical and object columns
        variables = []
        for col in adata.obs.columns:
            if adata.obs[col].dtype in ['object', 'category'] or adata.obs[col].dtype.name == 'category':
                variables.append(col)
    
    results = {}
    
    print("="*80)
    print(f"CLINICAL DATA EXPLORATION - Total cells: {adata.n_obs:,}")
    print("="*80)
    
    for var in variables:
        if var not in adata.obs.columns:
            print(f"⚠️  Variable '{var}' not found in adata.obs")
            continue
            
        print(f"\n📊 {var.upper()}")
        print("-" * 50)
        
        # Get value counts
        value_counts = adata.obs[var].value_counts(dropna=False)
        percentages = adata.obs[var].value_counts(normalize=True, dropna=False) * 100
        
        # Store results
        results[var] = {
            'categories': value_counts.index.tolist(),
            'counts': value_counts.values.tolist(),
            'percentages': percentages.values.tolist(),
            'n_categories': len(value_counts),
            'missing_values': adata.obs[var].isna().sum()
        }
        
        # Display results
        summary_df = pd.DataFrame({
            'Category': value_counts.index,
            'Count': value_counts.values,
            'Percentage': [f"{p:.1f}%" for p in percentages.values]
        })
        
        print(summary_df.to_string(index=False))
        
        # Additional statistics
        print(f"\n📈 Summary:")
        print(f"   • Total categories: {len(value_counts)}")
        print(f"   • Missing values: {adata.obs[var].isna().sum()}")
        print(f"   • Most common: {value_counts.index[0]} ({value_counts.iloc[0]:,} cells, {percentages.iloc[0]:.1f}%)")
        
        if len(value_counts) > 1:
            print(f"   • Least common: {value_counts.index[-1]} ({value_counts.iloc[-1]:,} cells, {percentages.iloc[-1]:.1f}%)")
    
    return results

# Function to explore patient-level statistics
def explore_patient_level(adata, patient_id_col='Patient_ID', sample_id_col='Sample_ID'):
    """Explore patient and sample level statistics"""
    print("\n" + "="*80)
    print("PATIENT & SAMPLE LEVEL STATISTICS")
    print("="*80)
    
    # Patient level
    n_patients = adata.obs[patient_id_col].nunique()
    print(f"👥 Total patients: {n_patients}")
    
    # Sample level  
    n_samples = adata.obs[sample_id_col].nunique()
    print(f"🧪 Total samples: {n_samples}")
    
    # Cells per patient
    cells_per_patient = adata.obs.groupby(patient_id_col).size()
    print(f"\n📊 Cells per patient:")
    print(f"   • Mean: {cells_per_patient.mean():.0f}")
    print(f"   • Median: {cells_per_patient.median():.0f}")
    print(f"   • Range: {cells_per_patient.min():.0f} - {cells_per_patient.max():.0f}")
    
    # Samples per patient
    samples_per_patient = adata.obs.groupby(patient_id_col)[sample_id_col].nunique()
    print(f"\n🧪 Samples per patient:")
    print(f"   • Mean: {samples_per_patient.mean():.1f}")
    print(f"   • Median: {samples_per_patient.median():.1f}")
    print(f"   • Range: {samples_per_patient.min():.0f} - {samples_per_patient.max():.0f}")
    
    return {
        'n_patients': n_patients,
        'n_samples': n_samples,
        'cells_per_patient': cells_per_patient,
        'samples_per_patient': samples_per_patient
    }

# Set color
def color_setting():
    # Cancer Cell style color palettes
    cancer_cell_colors = {
        'major_types': {
            'T': '#E31A1C',
            'Epithelial': '#1F78B4', 
            'Plasma': '#33A02C',
            'B': '#FF7F00',
            'Myeloid': '#6A3D9A',
            'Stromal': '#B15928',
            'Encothelial': '#FB9A99',
            'Mast': '#A6CEE3'
        },
        'tissue_types': {
            'tumor': '#D62728',
            'blood': '#2CA02C',
            'normal': '#1F77B4',
            'border': '#FF7F0E',
            'Unknown': '#CCCCCC'
        },
        'treatment': {
            'Anti-PD1': '#1F77B4',
            'untreated': '#FF7F0E',
            'Anti-PD1 plus CapeOx': '#2CA02C',
            'Anti-PD1 plus Celecoxib': '#D62728',
            'Unknown': '#CCCCCC'
        },
        'microsatellite': {
            'MSI': '#E31A1C',
            'MSS': '#1F78B4',
            'Unknown': '#CCCCCC'
        },
        'response': {
            'pCR': '#2CA02C',
            'non_pCR': '#D62728',
            'Unknown': '#CCCCCC'
        },
        'treatment_stage': {
            'Pre': '#1F77B4',
            'On': '#FF7F0E', 
            'Post': '#2CA02C',
            'Unknown': '#CCCCCC'
        },
        'gender': {
            'M': '#4CAF50',
            'F': '#FF9800',
            'Unknown': '#CCCCCC'
        }
    }

    return cancer_cell_colors

# Figure 1: Overall cell type composition
def plot_cell_type_composition(adata, figsize=(8, 6), save_path=None):
    """Plot cell type composition pie chart and bar plot"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Data preparation
    cell_counts = adata.obs['Major_type'].value_counts()
    colors = [adata.uns['Major_type_colors'][ct] for ct in cell_counts.index]
    
    # Pie chart
    wedges, texts, autotexts = ax1.pie(cell_counts.values, 
                                      labels=cell_counts.index,
                                      colors=colors,
                                      autopct='%1.1f%%',
                                      startangle=90,
                                      textprops={'fontsize': 8})
    
    ax1.set_title('Cell Type Distribution\n(n=1,427,414 cells)', fontsize=8, fontweight='bold', pad=10)
    
    # Bar plot
    bars = ax2.bar(range(len(cell_counts)), cell_counts.values, color=colors)
    ax2.set_xlabel('Cell Types', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Number of Cells', fontsize=10, fontweight='bold')
    ax2.set_title('Cell Type Counts', fontsize=10, fontweight='bold')
    ax2.set_xticks(range(len(cell_counts)))
    ax2.set_xticklabels(cell_counts.index, rotation=45, ha='right', fontsize=8)
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K'))
    
    # Add count labels on bars
    for bar, count in zip(bars, cell_counts.values):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                f'{count/1000:.0f}K', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    return colors

# Helper function to safely handle categorical data
def safe_fillna_categorical(series, fill_value='Unknown'):
    """Safely fill NA values in categorical series"""
    if pd.api.types.is_categorical_dtype(series):
        # Add the category if it doesn't exist
        if fill_value not in series.cat.categories:
            series = series.cat.add_categories([fill_value])
        return series.fillna(fill_value)
    else:
        return series.fillna(fill_value)

# Figure 3: Cell type composition by tissue type
def plot_celltype_by_tissue(adata, figsize=(8, 5), save_path=None):
    """Plot cell type composition across tissue types"""
    # Prepare data
    tissue_series = safe_fillna_categorical(adata.obs['tissue'], 'Unknown')
    tissue_celltype = pd.crosstab(tissue_series, 
                                 adata.obs['Major_type'], normalize='index') * 100
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Stacked bar plot (percentage)
    colors = [adata.uns['Major_type_colors'][ct] for ct in tissue_celltype.columns]
    tissue_celltype.plot(kind='bar', stacked=True, ax=ax1, color=colors, 
                        legend=False, width=0.8)
    ax1.set_title('Cell Type Composition by Tissue', fontweight='bold', fontsize=8)
    ax1.set_xlabel('Tissue Type', fontsize=7, fontweight='bold')
    ax1.set_ylabel('Percentage (%)', fontsize=7, fontweight='bold')
    ax1.set_xticklabels(tissue_celltype.index, rotation=45, ha='right', fontsize=6)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=6)
    
    # Absolute counts heatmap
    tissue_celltype_counts = pd.crosstab(tissue_series, 
                                        adata.obs['Major_type'])
    
    im = ax2.imshow(tissue_celltype_counts.values, cmap='Reds', aspect='auto')
    ax2.set_xticks(range(len(tissue_celltype_counts.columns)))
    ax2.set_yticks(range(len(tissue_celltype_counts.index)))
    ax2.set_xticklabels(tissue_celltype_counts.columns, rotation=45, ha='right', fontsize=6)
    ax2.set_yticklabels(tissue_celltype_counts.index, fontsize=6)
    ax2.set_title('Cell Counts Heatmap', fontweight='bold', fontsize=8)
    
    # Add text annotations
    for i in range(len(tissue_celltype_counts.index)):
        for j in range(len(tissue_celltype_counts.columns)):
            count = tissue_celltype_counts.iloc[i, j]
            ax2.text(j, i, f'{count/1000:.0f}K' if count > 1000 else str(count),
                    ha='center', va='center', fontsize=5,
                    color='white' if count > tissue_celltype_counts.values.max()/2 else 'black')
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    return None

# Figure 5: Patient and sample overview
def plot_patient_sample_overview(adata, figsize=(8, 6), save_path=None):
    """Plot patient and sample level statistics"""
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    # Cells per patient
    cells_per_patient = adata.obs.groupby('patient').size()
    axes[0, 0].hist(cells_per_patient, bins=20, color='#1F77B4', alpha=0.7, edgecolor='black')
    axes[0, 0].set_title('Cells per Patient Distribution', fontweight='bold', fontsize=8)
    axes[0, 0].set_xlabel('Number of Cells', fontsize=7)
    axes[0, 0].set_ylabel('Number of Patients', fontsize=7)
    axes[0, 0].axvline(cells_per_patient.median(), color='red', linestyle='--', 
                      label=f'Median: {cells_per_patient.median():.0f}')
    axes[0, 0].legend(fontsize=6)
    
    # Study composition
    study_counts = adata.obs['study'].value_counts()
    axes[1, 0].pie(study_counts.values, labels=study_counts.index, autopct='%1.1f%%',
                  colors=plt.cm.Set3(np.linspace(0, 1, len(study_counts))))
    axes[1, 0].set_title('Study Composition', fontweight='bold', fontsize=8)
    
    # Summary statistics text
    axes[1, 1].axis('off')
    summary_text = f"""
    Dataset Summary:
    
    Total Cells: {adata.n_obs:,}
    Total Patients: {adata.obs['patient'].nunique()}
    Total Samples: {adata.obs['sample_id'].nunique()}
    
    Cell Types: {adata.obs['Major_type'].nunique()}
    Studies: {adata.obs['study'].nunique()}
    
    Median cells/patient: {cells_per_patient.median():.0f}
    Mean cells/patient: {cells_per_patient.mean():.0f}
    """
    axes[1, 1].text(0.1, 0.9, summary_text, transform=axes[1, 1].transAxes,
                   fontsize=7, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

# ============================================================================
# ANNOTATION FUNCTION FOR EACH CELL TYPE
# ============================================================================

def preprocessing_for_subtype_anno(adata, cell_type, marker_genes, figurePath, resolution=0.5):
    """
    Annotate subtypes for a specific major cell type
    
    Parameters:
    - adata: AnnData object with all cells
    - cell_type: string, major cell type to analyze ('T', 'Myeloid', etc.)
    - marker_genes: dict, marker genes for this cell type
    - figurePath: string, path to save figures
    - resolution: float, clustering resolution
    """
    
    print(f"\n🎯 Annotating {cell_type} cells")
    print("="*60)
    
    # Step 1: Extract cells of this type
    if cell_type == 'B_Plasma':
        cell_mask = (adata.obs['Major_type'] == 'B') | (adata.obs['Major_type'] == 'Plasma')
    else:
        cell_mask = adata.obs['Major_type'] == cell_type
    print(f"Found {cell_mask.sum():,} {cell_type} cells")
    
    if cell_mask.sum() < 1000:
        print(f"⚠️ Too few {cell_type} cells for detailed analysis")
        return adata
    
    # Create subset
    cell_subset = sc.AnnData(
        X=adata[cell_mask].raw.X,
        obs=adata[cell_mask].obs.copy(),
        var=adata.raw.var.copy(),
        obsm=adata[cell_mask].obsm.copy()
    )
    
    # Step 2: Basic preprocessing
    print("Preprocessing...")
    sc.pp.normalize_total(cell_subset, target_sum=1e4)
    sc.pp.log1p(cell_subset)
    sc.pp.highly_variable_genes(cell_subset, n_top_genes=2000)
    sc.pp.scale(cell_subset)
    sc.tl.pca(cell_subset)
    sc.pp.neighbors(cell_subset, n_neighbors=10, n_pcs=50,  use_rep='X_pca_harmony')
    # sc.pp.neighbors(cell_subset, n_neighbors=15)
    sc.tl.umap(cell_subset)
    
    # Step 3: Clustering
    print("Clustering cells...")
    sc.tl.leiden(cell_subset, resolution=resolution, key_added='leiden')
    
    print(f"Found {len(cell_subset.obs['leiden'].unique())} clusters")
    
    # Step 4: Calculate marker gene scores
    print("Calculating marker gene scores...")
    all_marker_genes = []
    
    for subtype, genes in marker_genes.items():
        available_genes = [g for g in genes if g in cell_subset.var_names]
        all_marker_genes.extend(available_genes)

        # available_genes = [g for g in genes if g in cell_subset.var_names]
        # if available_genes:
        #     print(f"  {subtype}: {available_genes}")
        #     sc.tl.score_genes(cell_subset, available_genes, 
        #                       score_name=f'{subtype}_score', use_raw=False)
        #     all_marker_genes.extend(available_genes)
        # else:
        #     print(f"  {subtype}: No genes found in dataset")
    
    # Remove duplicates
    all_marker_genes = list(set(all_marker_genes))
    
    # Step 5: Differential expression analysis
    sc.tl.rank_genes_groups(cell_subset, groupby="leiden", method="wilcoxon")
    
    # Step 6: Visualizations
    cell_figurePath = os.path.join(figurePath, f"{cell_type}_analysis")
    os.makedirs(cell_figurePath, exist_ok=True)
    
    # Cluster overview
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(cell_subset, color='leiden', ax=axes[0], show=False, 
               legend_loc='on data', title=f'{cell_type} Clusters')
    sc.pl.umap(cell_subset, color='Major_type', ax=axes[1], show=False, 
               title=f'{cell_type} Cells')
    plt.tight_layout()
    plt.savefig(os.path.join(cell_figurePath, f"{cell_type}_clusters.pdf"), 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    # Marker genes dotplot
    if len(all_marker_genes) > 0:
        sc.pl.rank_genes_groups_dotplot(
            cell_subset, groupby="leiden", standard_scale="var", n_genes=8, show=False
        )
        plt.tight_layout()
        plt.savefig(os.path.join(cell_figurePath, f"{cell_type}_markers_dotplot.pdf"), 
                    dpi=300, bbox_inches='tight')
        plt.show()
        
        # # Marker gene expression plots
        # if len(all_marker_genes) <= 20:  # Only plot if reasonable number
        #     fig = plt.figure(figsize=(20, 15))
        #     n_genes = len(all_marker_genes)
        #     n_cols = 5
        #     n_rows = (n_genes + n_cols - 1) // n_cols
            
        #     for i, gene in enumerate(all_marker_genes):
        #         ax = plt.subplot(n_rows, n_cols, i+1)
        #         sc.pl.umap(cell_subset, color=gene, ax=ax, show=False, 
        #                   frameon=False, size=8, color_map='Reds')
        #         ax.set_title(gene, fontsize=10, fontweight='bold')
        #         ax.set_xlabel('')
        #         ax.set_ylabel('')
        #         ax.set_xticks([])
        #         ax.set_yticks([])
            
        #     plt.tight_layout()
        #     plt.savefig(os.path.join(cell_figurePath, f'{cell_type}_marker_expression.pdf'), 
        #                 dpi=300, bbox_inches='tight')
        #     plt.show()
    
    ## Save processed data
    cell_subset.write_h5ad(os.path.join(cell_figurePath, f"{cell_type}_processed.h5ad"))

    del cell_subset
    gc.collect()

    return None