## Functions for pre-processing
import os
import gzip
import h5py
import warnings
from pathlib import Path

import scanpy as sc
import scanpy.external as sce
import anndata as ad
import scrublet as scr

import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
from sklearn.metrics import silhouette_score

import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

#%% Data Loading functions
# Function to find file with or without .gz extension
def find_file(data_dir, base_name):
    gz_file = data_dir / f"{base_name}.gz"
    regular_file = data_dir / base_name
    
    if gz_file.exists():
        return gz_file
    elif regular_file.exists():
        return regular_file
    else:
        raise FileNotFoundError(f"Neither {gz_file} nor {regular_file} found")

def load_seurat_data(data_dir, prefix="seurat_data", transpose_matrix=True):
    """
    Load Seurat exported data into scanpy AnnData object
    
    Parameters:
    -----------
    data_dir : str
        Directory containing the exported files
    prefix : str
        Prefix used when exporting from Seurat
    transpose_matrix : bool
        Whether to transpose the matrix (Seurat exports genes x cells, 
        scanpy expects cells x genes)
    
    Returns:
    --------
    adata : AnnData
        Scanpy AnnData object ready for analysis
    """
    
    # File paths - check for both compressed and uncompressed versions
    data_dir = Path(data_dir)
    
    # Find files (with or without compression)
    mtx_file = find_file(data_dir, f"{prefix}_matrix.mtx")
    features_file = find_file(data_dir, f"{prefix}_features.tsv")
    barcodes_file = find_file(data_dir, f"{prefix}_barcodes.tsv")
    metadata_file = data_dir / f"{prefix}_metadata.csv"
    
    # Check if metadata file exists
    if not metadata_file.exists():
        raise FileNotFoundError(f"Required file not found: {metadata_file}")
    
    print(f"Loading data from {data_dir}")
    print(f"Files found:")
    print(f"  - Matrix: {mtx_file.name}")
    print(f"  - Features: {features_file.name}")
    print(f"  - Barcodes: {barcodes_file.name}")
    print(f"  - Metadata: {metadata_file.name}")
    
    # 1. Load count matrix (scanpy automatically handles .gz files)
    adata = sc.read_mtx(mtx_file)
    
    # 2. Load gene names (features) - handle compression
    if features_file.suffix == '.gz':
        features = pd.read_csv(features_file, sep='\t', header=None, 
                              names=['gene_id', 'gene_symbol', 'gene_type'],
                              compression='gzip')
    else:
        features = pd.read_csv(features_file, sep='\t', header=None, 
                              names=['gene_id', 'gene_symbol', 'gene_type'])
    
    # 3. Load cell barcodes - handle compression
    if barcodes_file.suffix == '.gz':
        barcodes = pd.read_csv(barcodes_file, sep='\t', header=None, 
                              names=['barcode'], compression='gzip')
    else:
        barcodes = pd.read_csv(barcodes_file, sep='\t', header=None, 
                              names=['barcode'])
    
    # 4. Load metadata
    metadata = pd.read_csv(metadata_file, index_col='cell_barcode')
    
    # 5. Transpose matrix if needed (Seurat exports genes x cells)
    if transpose_matrix:
        adata = adata.T
    
    # 6. Set gene names and cell barcodes
    adata.var_names = features['gene_symbol'].values
    adata.var['gene_id'] = features['gene_id'].values
    adata.var['gene_type'] = features['gene_type'].values
    
    adata.obs_names = barcodes['barcode'].values
    
    # 7. Add metadata to observations
    # Align metadata with cell barcodes
    aligned_metadata = metadata.reindex(adata.obs_names, fill_value=np.nan)
    
    # Add all metadata columns to adata.obs
    for col in aligned_metadata.columns:
        adata.obs[col] = aligned_metadata[col].values
    
    # 8. Make variable names unique (in case of duplicates)
    adata.var_names_make_unique()
    
    # 9. Set basic info
    adata.raw = adata  # Save raw data
    
    print(f"Loaded AnnData object:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Metadata columns: {list(adata.obs.columns)}")
    
    return adata

#%% Data preprocessing functions
def preprocess_single_dataset(adata, study_name):
    """
    Complete preprocessing pipeline for a single dataset
    """
    print(f"\n🔄 Preprocessing {study_name}...")
    
    # Store raw data
    adata.raw = adata

    # Basic preprocessing for downstream analysis
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    print(f"  ✅ Preprocessing completed for {study_name}")
    return adata

#%% Preparation for integration
def merge_datasets_for_integration(preprocessed_datasets):
    """
    Merge all preprocessed datasets into a single AnnData object for integration
    """
    print(f"\n{'='*60}")
    print("🔗 Merging datasets for integration")
    print(f"{'='*60}")
    
    # Collect all datasets
    adata_list = []
    for study_name, adata in preprocessed_datasets.items():
        print(f"  Adding {study_name}: {adata.shape}")
        adata_list.append(adata)
    
    # Merge datasets
    print("  🔄 Concatenating datasets...")
    merged_adata = ad.concat(adata_list, axis=0, join='outer', merge='unique')
    
    # Make observation names unique
    merged_adata.obs_names_make_unique()
    
    # Add batch information for integration
    merged_adata.obs['batch'] = merged_adata.obs["orig.ident"].astype('category')
    
    print(f"  ✅ Merged dataset shape: {merged_adata.shape}")
    print(f"  📊 Studies included: {merged_adata.obs['batch'].value_counts()}")
    
    return merged_adata


#%% Batch correction methods
def downsample_stratified(adata, n_sample=50000, batch_key='batch', 
                         quality_weight=True, random_state=42):
    """
    Perform stratified downsampling with optional quality weighting
    
    Parameters:
    -----------
    adata : AnnData
        Input data
    n_sample : int
        Target number of cells to sample
    batch_key : str
        Column name for batch stratification
    quality_weight : bool
        Whether to weight sampling by cell quality
    """
    print(f"🎯 Downsampling from {adata.n_obs:,} to {n_sample:,} cells...")
    
    np.random.seed(random_state)
    
    # Get batch information
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)
    
    # Calculate cells per batch (proportional to original)
    batch_counts = adata.obs[batch_key].value_counts()
    batch_proportions = batch_counts / batch_counts.sum()
    
    # Ensure minimum representation per batch
    min_per_batch = max(500, n_sample // (n_batches * 4))  # At least 500 or 1/4 of average
    
    sampled_indices = []
    
    for batch in batches:
        batch_mask = adata.obs[batch_key] == batch
        batch_cells = np.where(batch_mask)[0]
        
        # Calculate target for this batch
        target_batch = max(min_per_batch, int(n_sample * batch_proportions[batch]))
        target_batch = min(target_batch, len(batch_cells))  # Don't exceed available
        
        if quality_weight and 'n_genes_by_counts' in adata.obs.columns:
            # Weight by cell quality (higher gene count = higher probability)
            quality_scores = adata.obs.loc[batch_mask, 'n_genes_by_counts'].values
            # Normalize to probabilities
            weights = quality_scores / quality_scores.sum()
            
            # Sample with replacement=False
            batch_sample = np.random.choice(
                batch_cells, 
                size=target_batch, 
                replace=False, 
                p=weights
            )
        else:
            # Random sampling
            batch_sample = np.random.choice(
                batch_cells, 
                size=target_batch, 
                replace=False
            )
        
        sampled_indices.extend(batch_sample)
        print(f"  • {batch}: {len(batch_sample):,} cells (original: {len(batch_cells):,})")
    
    # Create downsampled data
    adata_downsampled = adata[sampled_indices].copy()
    
    print(f"  ✅ Downsampled to {adata_downsampled.n_obs:,} cells")
    print(f"  📊 Batch distribution: {adata_downsampled.obs[batch_key].value_counts().to_dict()}")
    
    return adata_downsampled

def prepare_for_integration(adata, n_top_genes=2000, downsample_n=None, 
                          batch_key='batch', quality_weight=True):
    """
    Prepare data for batch correction with optional downsampling
    """
    print("🔬 Preparing data for integration...")
    print(f"  Initial shape: {adata.shape}")
    
    # Downsample if requested
    if downsample_n is not None and adata.n_obs > downsample_n:
        adata = downsample_stratified(
            adata, 
            n_sample=downsample_n, 
            batch_key=batch_key,
            quality_weight=quality_weight
        )
    else:
        print(f"  No downsampling needed (target: {downsample_n}, actual: {adata.n_obs})")
    
    # Store raw data if not already done
    if adata.raw is None:
        adata.raw = adata
    
    # Find highly variable genes across all datasets
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=n_top_genes, 
        batch_key=batch_key,
        subset=True
    )
    
    print(f"  • Selected {adata.n_vars} highly variable genes")
    
    # Scale data for integration
    sc.pp.scale(adata, max_value=10)
    
    # Compute PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    
    print(f"  ✅ Data prepared: {adata.shape}")
    return adata

def apply_bbknn_integration(adata, batch_key='batch'):
    """
    Apply BBKNN batch correction
    """
    print(f"\n🔗 Applying BBKNN integration on {adata.n_obs:,} cells...")
    
    adata_bbknn = adata.copy()
    
    try:
        # Apply BBKNN
        sc.external.pp.bbknn(
            adata_bbknn, 
            batch_key=batch_key,
            neighbors_within_batch=3,
            n_pcs=50
        )
        
        # Compute UMAP
        sc.tl.umap(adata_bbknn)
        
        print("  ✅ BBKNN integration completed")
        return adata_bbknn
        
    except Exception as e:
        print(f"  ❌ BBKNN integration failed: {str(e)}")
        return None

def apply_harmony_integration(adata, batch_key='batch'):
    """
    Apply Harmony batch correction
    """
    print(f"\n🎵 Applying Harmony integration on {adata.n_obs:,} cells...")
    
    adata_harmony = adata.copy()
    
    try:
        # Apply Harmony
        sc.external.pp.harmony_integrate(
            adata_harmony, 
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_pca_harmony'
        )
        
        # Compute UMAP on Harmony-corrected data
        sc.pp.neighbors(adata_harmony, use_rep='X_pca_harmony', n_neighbors=15)
        sc.tl.umap(adata_harmony)
        
        print("  ✅ Harmony integration completed")
        return adata_harmony
        
    except Exception as e:
        print(f"  ❌ Harmony integration failed: {str(e)}")
        return None

def apply_mnn_integration(adata, batch_key='batch'):
    """
    Apply MNN (Mutual Nearest Neighbors) batch correction
    """
    print(f"\n🔗 Applying MNN integration on {adata.n_obs:,} cells...")
    
    adata_mnn = adata.copy()
    
    try:
        # Apply MNN correction
        sc.external.pp.mnn_correct(
            adata_mnn,
            batch_key=batch_key,
            save_raw=True,  # Save raw data
            n_jobs=1,  # Use single core to avoid memory issues
            batch_categories=None  # Auto-detect batch categories
        )
        
        # The corrected data is stored in adata.X
        # Compute PCA on corrected data
        sc.tl.pca(adata_mnn, svd_solver='arpack', n_comps=50)
        
        # Compute neighbors and UMAP
        sc.pp.neighbors(adata_mnn, n_neighbors=15, n_pcs=50)
        sc.tl.umap(adata_mnn)
        
        print("  ✅ MNN integration completed")
        return adata_mnn
        
    except Exception as e:
        print(f"  ❌ MNN integration failed: {str(e)}")
        return None

def apply_scanorama_integration(adata, batch_key='batch'):
    """
    Apply Scanorama batch correction
    """
    print(f"\n🔄 Applying Scanorama integration on {adata.n_obs:,} cells...")
    
    adata_scanorama = adata.copy()
    
    try:        
        # Apply Scanorama
        sc.external.pp.scanorama_integrate(
            adata_scanorama, 
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_scanorama'
        )  
        
        # Compute neighbors and UMAP on corrected data
        sc.pp.neighbors(adata_scanorama, use_rep='X_scanorama', n_neighbors=15)
        sc.tl.umap(adata_scanorama)
        
        print("  ✅ Scanorama integration completed")
        return adata_scanorama
        
    except Exception as e:
        print(f"  ❌ Scanorama integration failed: {str(e)}")
        return None

def compute_integration_metrics(adata, batch_key='batch', n_sample=5000):
    """
    Compute metrics to evaluate integration quality
    """
    # Sample data if too large for metrics computation
    adata_metrics = adata.copy()
    if adata_metrics.n_obs > n_sample:
        print(f"    📊 Subsampling to {n_sample} cells for metrics computation")
        sc.pp.subsample(adata_metrics, n_obs=n_sample)
    
    # Get batch labels
    batch_labels = adata_metrics.obs[batch_key].astype('category').cat.codes
    
    # Compute silhouette score (lower is better for batch correction)
    try:
        sil_score = silhouette_score(adata_metrics.obsm['X_umap'], batch_labels)
        return {'silhouette_batch': sil_score}
    except:
        return {'silhouette_batch': np.nan}

def visualize_integration_results(adata_dict, batch_key='batch', figsize=(20, 12), 
                                n_cells_plot=10000):
    """
    Create comprehensive visualization of integration results
    """
    print(f"\n📊 Creating integration comparison plots (max {n_cells_plot:,} cells per plot)...")
    
    n_methods = len(adata_dict)
    fig, axes = plt.subplots(3, n_methods, figsize=figsize)
    
    if n_methods == 1:
        axes = axes.reshape(3, 1)
    
    # Plot batch effects
    for i, (method, adata) in enumerate(adata_dict.items()):
        if adata is not None:
            # Sample cells for visualization if too many
            if adata.n_obs > n_cells_plot:
                print(f"    📊 Subsampling {method} to {n_cells_plot:,} cells for visualization")
                adata_vis = adata[np.random.choice(adata.n_obs, n_cells_plot, replace=False)]
            else:
                adata_vis = adata
            
            # Batch plot
            sc.pl.umap(
                adata_vis, 
                color=batch_key, 
                ax=axes[0, i], 
                show=False, 
                frameon=False,
                title=f'{method} - Batch ({adata_vis.n_obs:,} cells)'
            )
                
            # study plot  
            sc.pl.umap(
                adata_vis, 
                color='study', 
                ax=axes[1, i], 
                show=False, 
                frameon=False,
                title=f'{method} - Study ({adata_vis.n_obs:,} cells)'
            )

            # Tissue plot  
            sc.pl.umap(
                adata_vis, 
                color='tissue', 
                ax=axes[2, i], 
                show=False, 
                frameon=False,
                title=f'{method} - Tissue ({adata_vis.n_obs:,} cells)'
            )
        else:
            axes[0, i].text(0.5, 0.5, f'{method}\nFailed', 
                          ha='center', va='center', transform=axes[0, i].transAxes)
            axes[1, i].text(0.5, 0.5, f'{method}\nFailed', 
                          ha='center', va='center', transform=axes[1, i].transAxes)
            axes[2, i].text(0.5, 0.5, f'{method}\nFailed', 
                          ha='center', va='center', transform=axes[2, i].transAxes)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def compare_all_integration_methods(adata_prepared, batch_key='batch'):
    """
    Apply all integration methods and compare results
    """
    print(f"\n{'='*80}")
    print(f"🚀 Comprehensive Batch Correction Comparison ({adata_prepared.n_obs:,} cells)")
    print(f"{'='*80}")
    
    # Store original (uncorrected) data with UMAP
    print(f"\n📍 Computing uncorrected UMAP on {adata_prepared.n_obs:,} cells...")
    adata_uncorrected = adata_prepared.copy()
    sc.pp.neighbors(adata_uncorrected, n_neighbors=15)
    sc.tl.umap(adata_uncorrected)
    
    # Dictionary to store results
    results = {'Uncorrected': adata_uncorrected}
    
    # Apply each method
    methods = [
        ('Harmony', apply_harmony_integration),
        ('BBKNN', apply_bbknn_integration),
        ('MNN', apply_mnn_integration),
        ('Scanorama', apply_scanorama_integration)
    ]
    
    for method_name, method_func in methods:
        try:
            result = method_func(adata_prepared.copy(), batch_key=batch_key)
            results[method_name] = result
        except Exception as e:
            print(f"❌ {method_name} failed: {str(e)}")
            results[method_name] = None
    
    # Compute metrics
    print(f"\n📊 Computing integration metrics...")
    metrics_df = pd.DataFrame()
    
    for method_name, adata_result in results.items():
        if adata_result is not None:
            metrics = compute_integration_metrics(adata_result.copy(), batch_key=batch_key)
            metrics['method'] = method_name
            metrics['n_cells'] = adata_result.n_obs
            metrics['n_genes'] = adata_result.n_vars
            metrics_df = pd.concat([metrics_df, pd.DataFrame([metrics])], ignore_index=True)
    
    # Print metrics
    print("\n📈 Integration Quality Metrics:")
    print("=" * 60)
    print(metrics_df[['method', 'silhouette_batch', 'n_cells', 'n_genes']].round(3))
    print("\nNote: Lower silhouette_batch score indicates better batch mixing")
    
    # Create visualizations
    fig = visualize_integration_results(results, batch_key=batch_key)
    
    # Print summary
    print(f"\n{'='*80}")
    print("🎯 Integration Summary:")
    successful_methods = [k for k, v in results.items() if v is not None]
    failed_methods = [k for k, v in results.items() if v is None]
    
    print(f"✅ Successful methods: {successful_methods}")
    if failed_methods:
        print(f"❌ Failed methods: {failed_methods}")
    
    print(f"{'='*80}")
    
    return results, metrics_df, fig

def save_integration_results(results, save_path="./integration_results/"):
    """
    Save integration results with downsampling info
    """
    import os
    os.makedirs(save_path, exist_ok=True)
    
    for method_name, adata in results.items():
        if adata is not None:
            filename = f"{save_path}/{method_name.lower()}_integrated_n{adata.n_obs}.h5ad"
            adata.write(filename)
            print(f"💾 Saved {method_name} results ({adata.n_obs:,} cells) to {filename}")

# Two-phase workflow functions
def quick_method_comparison(merged_adata, batch_key='batch', n_sample=20000, n_top_genes=2000):
    """
    Phase 1: Quick method comparison with small sample
    """
    print(f"\n🚀 PHASE 1: Quick Method Comparison (n={n_sample:,})")
    print("="*60)
    
    # Prepare data with downsampling
    adata_small = prepare_for_integration(
        merged_adata, 
        n_top_genes=n_top_genes,
        downsample_n=n_sample,
        batch_key=batch_key
    )
    
    # Compare methods
    results, metrics_df, fig = compare_all_integration_methods(adata_small, batch_key=batch_key)
    
    return results, metrics_df, fig

def assign_data_classes_for_obs(adata, numerical_cols, categorical_cols):
    """
    Automatically assigns appropriate data classes (dtypes) to columns in
    the .obs DataFrame of an AnnData object.

    This function identifies common numerical and categorical columns based on
    typical single-cell data conventions and converts their dtypes accordingly.
    - Numerical columns are converted to float.
    - Categorical columns are converted to pandas.Categorical dtype.
    - Columns not explicitly listed or inferred will retain their original dtype.

    Args:
        adata: The AnnData object whose .obs DataFrame needs type assignment.

    Returns:
        The modified AnnData object with updated column dtypes in .obs.
    """

    print("--- Starting data class assignment for adata.obs ---")

    # Process numerical columns
    for col in numerical_cols:
        if col in adata.obs.columns:
            original_dtype = adata.obs[col].dtype
            try:
                # Convert to numeric, coercing errors to NaN.
                # This handles cases where numbers might be stored as strings or mixed types.
                adata.obs[col] = pd.to_numeric(adata.obs[col], errors='coerce')
                # Convert to float to ensure consistency for numerical data
                adata.obs[col] = adata.obs[col].astype(float)
                if original_dtype != adata.obs[col].dtype:
                    print(f"  Converted numerical column '{col}' from {original_dtype} to {adata.obs[col].dtype}")
            except Exception as e:
                print(f"  Warning: Could not convert '{col}' to numeric. Error: {e}. Keeping original dtype.")
        else:
            print(f"  Numerical column '{col}' not found in adata.obs. Skipping.")

    # Process categorical columns
    for col in categorical_cols:
        if col in adata.obs.columns:
            original_dtype = adata.obs[col].dtype
            try:
                # For categorical columns, first ensure they are strings (especially for 'object' dtypes
                # that might contain mixed types like booleans or NaNs, as seen with 'Preoperative_Chemotherapy').
                # Then convert to 'category'. This handles 'nan' and 'False' becoming string categories 'nan' and 'False'.
                adata.obs[col] = adata.obs[col].astype(str).astype('category')
                if original_dtype != adata.obs[col].dtype:
                    print(f"  Converted categorical column '{col}' from {original_dtype} to {adata.obs[col].dtype}")
            except Exception as e:
                print(f"  Warning: Could not convert '{col}' to category. Error: {e}. Keeping original dtype.")
        else:
            print(f"  Categorical column '{col}' not found in adata.obs. Skipping.")

    for col in adata.obs.columns:
        if col not in numerical_cols and col not in categorical_cols:
            del adata.obs[col] 
            # If the column is not explicitly listed, we keep its original dtype
            print(f"  Removing original column '{col}'")

    print("--- Data class assignment complete ---")
    return adata