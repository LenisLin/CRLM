## Functions for deconvolution and analysis of spatial transcriptomics data

import scanpy as sc
import pandas as pd
import numpy as np

def adaptive_celltype_sampling(adata, celltype_col='cell_type', max_cells=10000, 
                             random_state=42, verbose=True):
    """
    Perform adaptive sampling on each cell type based on abundance.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing single-cell data
    celltype_col : str
        Column name in adata.obs containing cell type annotations
    max_cells : int
        Maximum number of cells to keep per cell type (threshold for downsampling)
    random_state : int
        Random seed for reproducible sampling
    verbose : bool
        Whether to print sampling information
        
    Returns:
    --------
    AnnData
        New AnnData object with adaptively sampled cells
    """
    
    # Check if celltype column exists
    if celltype_col not in adata.obs.columns:
        raise ValueError(f"Column '{celltype_col}' not found in adata.obs")
    
    # Get cell type counts
    celltype_counts = adata.obs[celltype_col].value_counts()
    
    if verbose:
        print(f"Original dataset: {adata.n_obs:,} cells, {len(celltype_counts)} cell types")
        print(f"Sampling threshold: {max_cells:,} cells per type")
        print("\nCell type distribution (before sampling):")
        for ct, count in celltype_counts.items():
            status = "DOWNSAMPLE" if count > max_cells else "KEEP ALL"
            print(f"  {ct}: {count:,} cells -> {status}")
    
    # Sample each cell type
    sampled_indices = []
    sampling_summary = []
    
    np.random.seed(random_state)
    
    for celltype in celltype_counts.index:
        # Get indices for this cell type
        celltype_mask = adata.obs[celltype_col] == celltype
        celltype_indices = np.where(celltype_mask)[0]
        
        original_count = len(celltype_indices)
        
        if original_count <= max_cells:
            # Keep all cells for this type
            selected_indices = celltype_indices
            final_count = original_count
            action = "kept_all"
        else:
            # Downsample this cell type
            selected_indices = np.random.choice(
                celltype_indices, 
                size=max_cells, 
                replace=False
            )
            final_count = max_cells
            action = "downsampled"
        
        sampled_indices.extend(selected_indices)
        sampling_summary.append({
            'cell_type': celltype,
            'original_count': original_count,
            'final_count': final_count,
            'action': action,
            'reduction_ratio': final_count / original_count
        })
    
    # Create new AnnData object with sampled cells
    adata_sampled = adata[sampled_indices].copy()
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(sampling_summary)
    
    if verbose:
        print(f"\nSampling completed!")
        print(f"Final dataset: {adata_sampled.n_obs:,} cells")
        print(f"Total reduction: {adata_sampled.n_obs / adata.n_obs:.2%} of original")
        
        print(f"\nDetailed sampling results:")
        downsampled_types = summary_df[summary_df['action'] == 'downsampled']
        kept_all_types = summary_df[summary_df['action'] == 'kept_all']
        
        print(f"  Cell types downsampled: {len(downsampled_types)}")
        print(f"  Cell types kept all: {len(kept_all_types)}")
        
        if len(downsampled_types) > 0:
            print(f"\nDownsampled cell types:")
            for _, row in downsampled_types.iterrows():
                print(f"  {row['cell_type']}: {row['original_count']:,} -> {row['final_count']:,} "
                      f"({row['reduction_ratio']:.1%})")
        
        if len(kept_all_types) > 0:
            print(f"\nCell types with all cells kept:")
            for _, row in kept_all_types.iterrows():
                print(f"  {row['cell_type']}: {row['final_count']:,} cells")
    
    # Add sampling information to adata
    adata_sampled.uns['sampling_info'] = {
        'method': 'adaptive_celltype_sampling',
        'max_cells_per_type': max_cells,
        'random_state': random_state,
        'summary': summary_df.to_dict('records'),
        'original_n_obs': adata.n_obs,
        'final_n_obs': adata_sampled.n_obs
    }
    
    return adata_sampled, summary_df


def plot_sampling_summary(summary_df, figsize=(12, 6)):
    """
    Plot sampling summary to visualize the effects of adaptive sampling.
    
    Parameters:
    -----------
    summary_df : pd.DataFrame
        Summary DataFrame returned by adaptive_celltype_sampling
    figsize : tuple
        Figure size for the plots
    """
    import matplotlib.pyplot as plt
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Plot 1: Before vs After cell counts
    cell_types = summary_df['cell_type'].values
    original_counts = summary_df['original_count'].values
    final_counts = summary_df['final_count'].values
    
    x_pos = np.arange(len(cell_types))
    width = 0.35
    
    ax1.bar(x_pos - width/2, original_counts, width, label='Original', alpha=0.7)
    ax1.bar(x_pos + width/2, final_counts, width, label='After Sampling', alpha=0.7)
    
    ax1.set_xlabel('Cell Types')
    ax1.set_ylabel('Number of Cells')
    ax1.set_title('Cell Counts: Before vs After Sampling')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(cell_types, rotation=45, ha='right')
    ax1.legend()
    ax1.set_yscale('log')  # Log scale for better visualization
    
    # Plot 2: Reduction ratios
    colors = ['red' if action == 'downsampled' else 'green' 
              for action in summary_df['action']]
    
    ax2.bar(x_pos, summary_df['reduction_ratio'], color=colors, alpha=0.7)
    ax2.set_xlabel('Cell Types')
    ax2.set_ylabel('Retention Ratio')
    ax2.set_title('Cell Retention Ratio by Type')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(cell_types, rotation=45, ha='right')
    ax2.axhline(y=1.0, color='black', linestyle='--', alpha=0.5)
    ax2.set_ylim(0, 1.1)
    
    # Add legend for colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='red', alpha=0.7, label='Downsampled'),
                      Patch(facecolor='green', alpha=0.7, label='Kept All')]
    ax2.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.show()

## Filter mitochondrial genes
def filter_mitochondrial_genes(adata):
    """
    Filter mitochondrial genes from the AnnData object.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing single-cell data
    
    Returns:
    --------
    AnnData
        New AnnData object with mitochondrial genes removed
    """
    adata.var['MT_gene'] = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obsm['MT'] = adata[:, adata.var['MT_gene'].values].X.toarray()
    return adata[:, ~adata.var['MT_gene'].values]  # Filter out mitochondrial genes