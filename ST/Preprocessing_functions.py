## Functions for preprocessing the data

def combine_stereo_process_info(adata, data):
    ## Assign the dimentinal reduction results
    # Access through the tool results structure
    try:
        print("Checking tl.result structure:")
        if hasattr(data, 'tl') and hasattr(data.tl, 'result'):
            
            # Copy PCA if available
            if 'pca' in data.tl.result:
                pca_result = data.tl.result['pca']
                # pca_result = pca_result.reindex(adata.obs_names)
                # adata.obsm['X_pca'] = pca_result.copy()
                adata.obsm['X_pca'] = np.array(pca_result)
                print("✓ Copied PCA from tl.result")
            
            # Copy UMAP if available
            if 'umap' in data.tl.result:
                umap_result = data.tl.result['umap']
                # umap_result = umap_result.reindex(adata.obs_names)
                # adata.obsm['X_umap'] = umap_result.copy()
                adata.obsm['X_umap'] = np.array(umap_result)
                print("✓ Copied UMAP from tl.result")
            
            # Copy neighbors if available
            if 'neighbors' in data.tl.result:
                neighbors_result = data.tl.result['neighbors']
                if 'connectivities' in neighbors_result:
                    adata.obsp['connectivities'] = neighbors_result['connectivities'].copy()
                    print("✓ Copied connectivities from tl.result")
                if 'nn_dist' in neighbors_result:
                    adata.obsp['distances'] = neighbors_result['nn_dist'].copy()
                    print("✓ Copied distances from tl.result")
            
            # Copy spatial coordinates
            if hasattr(data, 'position'):
                adata.obsm['spatial'] = data.position.copy()
                print("✓ Copied spatial coordinates")
        
        return adata

    except Exception as e:
        print(f"tl.result access failed: {e}")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import scanpy as sc

def identify_tma_cores(adata, eps=None, min_samples=50, standardize=True, plot=True):
    """
    Identify TMA cores using DBSCAN clustering on spatial coordinates
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with spatial coordinates in .obs['x'] and .obs['y']
    eps : float, optional
        The maximum distance between two samples for them to be considered 
        as in the same neighborhood. If None, will be estimated automatically.
    min_samples : int, default=50
        The number of samples in a neighborhood for a point to be considered 
        as a core point.
    standardize : bool, default=True
        Whether to standardize coordinates before clustering
    plot : bool, default=True
        Whether to create visualization plots
    
    Returns:
    --------
    None (modifies adata.obs inplace by adding 'tma_core' column)
    """
    
    print("=== TMA Core Identification using DBSCAN ===")
    
    # Extract spatial coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
        x_coords = coords[:, 0]
        y_coords = coords[:, 1]
    elif 'x' in adata.obs.columns and 'y' in adata.obs.columns:
        x_coords = adata.obs['x'].values
        y_coords = adata.obs['y'].values
    else:
        raise ValueError("No spatial coordinates found. Need either 'spatial' in obsm or 'x','y' in obs")
    
    # Prepare coordinate matrix
    coordinates = np.column_stack([x_coords, y_coords])
    print(f"Coordinate range: X=[{x_coords.min():.1f}, {x_coords.max():.1f}], Y=[{y_coords.min():.1f}, {y_coords.max():.1f}]")
    
    # Standardize coordinates if requested
    if standardize:
        scaler = StandardScaler()
        coordinates_scaled = scaler.fit_transform(coordinates)
        print("Coordinates standardized for clustering")
    else:
        coordinates_scaled = coordinates
    
    # Automatic eps estimation if not provided
    if eps is None:
        # Calculate distances to k-th nearest neighbor
        from sklearn.neighbors import NearestNeighbors
        k = min_samples
        neighbors = NearestNeighbors(n_neighbors=k)
        neighbors.fit(coordinates_scaled)
        distances, _ = neighbors.kneighbors(coordinates_scaled)
        distances = np.sort(distances[:, k-1], axis=0)
        
        # Use elbow method to estimate eps
        # Take the point where the distance increases most rapidly
        diff = np.diff(distances)
        eps = distances[np.argmax(diff)]
        print(f"Estimated eps: {eps:.4f}")
        
        if plot:
            plt.figure(figsize=(8, 4))
            plt.subplot(1, 2, 1)
            plt.plot(distances)
            plt.axhline(y=eps, color='r', linestyle='--', label=f'Estimated eps={eps:.4f}')
            plt.xlabel('Points sorted by distance')
            plt.ylabel(f'{k}-th nearest neighbor distance')
            plt.title('K-distance Graph for eps estimation')
            plt.legend()
    
    # Perform DBSCAN clustering
    print(f"Running DBSCAN with eps={eps:.4f}, min_samples={min_samples}")
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    cluster_labels = dbscan.fit_predict(coordinates_scaled)
    
    # Analyze results
    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise = list(cluster_labels).count(-1)
    
    print(f"Number of clusters found: {n_clusters}")
    print(f"Number of noise points: {n_noise}")
    print(f"Percentage of data clustered: {((len(cluster_labels) - n_noise) / len(cluster_labels)) * 100:.1f}%")
    
    # Handle noise points and relabel clusters
    if n_noise > 0:
        print(f"Warning: {n_noise} points classified as noise (outliers)")
        # Assign noise points to nearest cluster
        from sklearn.neighbors import NearestNeighbors
        
        # Get points that are not noise
        valid_mask = cluster_labels != -1
        valid_coords = coordinates_scaled[valid_mask]
        valid_labels = cluster_labels[valid_mask]
        
        if len(valid_coords) > 0:
            # Find nearest valid point for each noise point
            noise_mask = cluster_labels == -1
            noise_coords = coordinates_scaled[noise_mask]
            
            if len(noise_coords) > 0:
                nn = NearestNeighbors(n_neighbors=1)
                nn.fit(valid_coords)
                _, indices = nn.kneighbors(noise_coords)
                
                # Assign noise points to nearest cluster
                noise_labels = valid_labels[indices.flatten()]
                cluster_labels[noise_mask] = noise_labels
                print(f"Assigned {len(noise_coords)} noise points to nearest clusters")
    
    # Relabel clusters to start from 1 and be sequential
    unique_labels = np.unique(cluster_labels)
    label_mapping = {old_label: new_label + 1 for new_label, old_label in enumerate(unique_labels)}
    final_labels = np.array([label_mapping[label] for label in cluster_labels])
    
    # Add to adata
    adata.obs['tma_core'] = final_labels.astype(str)
    adata.obs['tma_core'] = adata.obs['tma_core'].astype('category')

    print(f"TMA cores identified and added to adata.obs['tma_core']")
    print(f"Core distribution:")
    core_counts = adata.obs['tma_core'].value_counts().sort_index()
    for core, count in core_counts.items():
        print(f"  Core {core}: {count:,} cells")

    # Visualization
    if plot:
        # Create consistent color mapping for all plots
        unique_cores = sorted(core_counts.index)
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_cores)))  # Use Set1 colormap for distinct colors
        core_color_map = dict(zip(unique_cores, colors))
        
        # Create color array for scatter plots
        color_array = [core_color_map[core] for core in adata.obs['tma_core']]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Plot 1: Original coordinates with clusters using consistent colors
        ax1 = axes[0, 0]
        for core in unique_cores:
            core_mask = adata.obs['tma_core'] == core
            ax1.scatter(x_coords[core_mask], y_coords[core_mask], 
                    c=[core_color_map[core]], label=f'Core {core}', 
                    s=0.8, alpha=0.7)
        ax1.set_xlabel('X coordinate')
        ax1.set_ylabel('Y coordinate')
        ax1.set_title('TMA Cores Identified by DBSCAN')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)
        
        # Plot 2: Core distribution histogram with consistent colors
        ax2 = axes[0, 1]
        bar_colors = [core_color_map[core] for core in core_counts.index]
        bars = ax2.bar(range(len(core_counts)), core_counts.values, color=bar_colors)
        ax2.set_xlabel('TMA Core')
        ax2.set_ylabel('Number of Cells')
        ax2.set_title('Cell Count per TMA Core')
        ax2.set_xticks(range(len(core_counts)))
        ax2.set_xticklabels([f'Core {core}' for core in core_counts.index], rotation=45)
        
        # Add value labels on bars
        for i, (bar, count) in enumerate(zip(bars, core_counts.values)):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(core_counts.values)*0.01,
                    f'{count:,}', ha='center', va='bottom', fontsize=9)
        
        # Plot 3: Spatial distribution with numbered core boundaries
        ax3 = axes[1, 0]
        for core in unique_cores:
            core_mask = adata.obs['tma_core'] == core
            ax3.scatter(x_coords[core_mask], y_coords[core_mask], 
                    c=[core_color_map[core]], label=f'Core {core}', 
                    s=0.8, alpha=0.7)
        ax3.set_xlabel('X coordinate')
        ax3.set_ylabel('Y coordinate')
        ax3.set_title('TMA Cores with Consistent Color Labels')
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)
        
        # Plot 4: Core centroid analysis with color coding
        ax4 = axes[1, 1]
        centroids = []
        for core in unique_cores:
            core_mask = adata.obs['tma_core'] == core
            centroid_x = x_coords[core_mask].mean()
            centroid_y = y_coords[core_mask].mean()
            centroids.append((centroid_x, centroid_y))
            
            # Plot core points with consistent color
            ax4.scatter(x_coords[core_mask], y_coords[core_mask], 
                    c=[core_color_map[core]], s=0.5, alpha=0.4, 
                    label=f'Core {core}')
            
            # Plot centroid with matching color but different marker
            ax4.scatter(centroid_x, centroid_y, s=200, marker='x', 
                    c=[core_color_map[core]], linewidths=4, 
                    edgecolors='black', linewidth=1)
            
            # Add core label with background for better visibility
            ax4.annotate(f'C{core}', (centroid_x, centroid_y), 
                        xytext=(8, 8), textcoords='offset points',
                        fontsize=12, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor=core_color_map[core], 
                                alpha=0.7, edgecolor='black'))
        
        ax4.set_xlabel('X coordinate')
        ax4.set_ylabel('Y coordinate')
        ax4.set_title('TMA Cores with Centroids (Color-Coded)')
        
        plt.tight_layout()
        plt.show()
        
        # Create a separate figure for core correspondence verification
        fig2, ax_verify = plt.subplots(1, 1, figsize=(12, 8))
        
        # Plot all cores with large, clearly labeled points
        for core in unique_cores:
            core_mask = adata.obs['tma_core'] == core
            sample_indices = np.where(core_mask)[0]
            
            # Sample points for clearer visualization if too many points
            if len(sample_indices) > 1000:
                sample_indices = np.random.choice(sample_indices, 1000, replace=False)
            
            ax_verify.scatter(x_coords[sample_indices], y_coords[sample_indices], 
                            c=[core_color_map[core]], s=3, alpha=0.8, 
                            label=f'Core {core} ({core_counts[core]:,} cells)')
        
        ax_verify.set_xlabel('X coordinate')
        ax_verify.set_ylabel('Y coordinate')
        ax_verify.set_title('TMA Core Correspondence Verification\n(Consistent Colors Across All Visualizations)')
        ax_verify.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax_verify.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # Print color mapping for reference
        print(f"\n=== Color Mapping Reference ===")
        for core in unique_cores:
            rgb = core_color_map[core][:3]
            print(f"Core {core}: RGB({rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f})")
        
        # Additional quality metrics
        print(f"\n=== Quality Metrics ===")
        
        # Calculate core compactness (average distance from centroid)
        compactness_scores = []
        for i, core in enumerate(unique_cores):
            core_mask = adata.obs['tma_core'] == core
            core_coords = np.column_stack([x_coords[core_mask], y_coords[core_mask]])
            centroid = centroids[i]
            distances = np.sqrt(np.sum((core_coords - centroid)**2, axis=1))
            avg_distance = np.mean(distances)
            std_distance = np.std(distances)
            compactness_scores.append(avg_distance)
            print(f"Core {core} compactness - Mean: {avg_distance:.1f}, Std: {std_distance:.1f}")
        
        # Calculate inter-core distances
        print(f"\n=== Inter-core Distances ===")
        for i in range(len(centroids)):
            for j in range(i+1, len(centroids)):
                dist = np.sqrt(np.sum((np.array(centroids[i]) - np.array(centroids[j]))**2))
                core_i = unique_cores[i]
                core_j = unique_cores[j]
                print(f"Distance between Core {core_i} and Core {core_j}: {dist:.1f}")
        
        # Calculate silhouette-like score for core separation
        print(f"\n=== Core Separation Quality ===")
        avg_compactness = np.mean(compactness_scores)
        min_inter_distance = float('inf')
        for i in range(len(centroids)):
            for j in range(i+1, len(centroids)):
                dist = np.sqrt(np.sum((np.array(centroids[i]) - np.array(centroids[j]))**2))
                min_inter_distance = min(min_inter_distance, dist)
        
        separation_ratio = min_inter_distance / avg_compactness if avg_compactness > 0 else 0
        print(f"Separation ratio (min inter-core distance / avg compactness): {separation_ratio:.2f}")
        print(f"Higher values indicate better core separation")

def optimize_dbscan_parameters(adata, eps_range=None, min_samples_range=None, target_clusters=9):
    """
    Optimize DBSCAN parameters to get approximately the target number of clusters
    """
    
    # Extract coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    else:
        coords = np.column_stack([adata.obs['x'].values, adata.obs['y'].values])
    
    # Standardize coordinates
    scaler = StandardScaler()
    coords_scaled = scaler.fit_transform(coords)
    
    # Define parameter ranges
    if eps_range is None:
        eps_range = np.linspace(0.1, 2.0, 20)
    if min_samples_range is None:
        min_samples_range = [20, 50, 100, 200]
    
    results = []
    
    print("Optimizing DBSCAN parameters...")
    for eps in eps_range:
        for min_samples in min_samples_range:
            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            labels = dbscan.fit_predict(coords_scaled)
            
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise = list(labels).count(-1)
            noise_ratio = n_noise / len(labels)
            
            results.append({
                'eps': eps,
                'min_samples': min_samples,
                'n_clusters': n_clusters,
                'n_noise': n_noise,
                'noise_ratio': noise_ratio,
                'cluster_diff': abs(n_clusters - target_clusters)
            })
    
    # Convert to DataFrame and find best parameters
    results_df = pd.DataFrame(results)
    
    # Filter out results with too much noise (>10%) or too few clusters
    good_results = results_df[(results_df['noise_ratio'] < 0.1) & (results_df['n_clusters'] >= 3)]
    
    if len(good_results) > 0:
        # Find parameters that give closest to target number of clusters
        best_params = good_results.loc[good_results['cluster_diff'].idxmin()]
        print(f"\nOptimal parameters found:")
        print(f"  eps: {best_params['eps']:.3f}")
        print(f"  min_samples: {best_params['min_samples']}")
        print(f"  n_clusters: {best_params['n_clusters']}")
        print(f"  noise_ratio: {best_params['noise_ratio']:.3f}")
        
        return best_params['eps'], int(best_params['min_samples'])
    else:
        print("No good parameters found with current ranges. Try expanding the search ranges.")
        return 0.5, 50

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def filter_non_coding_genes(adata, gene_name_col='real_gene_name', copy=False, verbose=True):
    """
    Filter out mitochondrial genes and other non-classic protein-coding genes
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    gene_name_col : str, default='real_gene_name'
        Column name in adata.var that contains gene symbols
    copy : bool, default=False
        Whether to return a copy or modify in place
    verbose : bool, default=True
        Whether to print filtering statistics
    
    Returns:
    --------
    AnnData or None
        Filtered AnnData object if copy=True, otherwise modifies in place
    """
    
    if copy:
        adata = adata.copy()
    
    if verbose:
        print("=== Gene Filtering Analysis ===")
        print(f"Starting with {adata.n_vars:,} genes and {adata.n_obs:,} cells")
    
    # Get gene names
    if gene_name_col in adata.var.columns:
        gene_names = adata.var[gene_name_col].astype(str)
    else:
        print(f"Warning: {gene_name_col} not found in adata.var.columns")
        print(f"Available columns: {adata.var.columns.tolist()}")
        gene_names = adata.var_names.astype(str)
    
    # Initialize filter mask (True = keep gene)
    keep_genes = pd.Series(True, index=adata.var_names)
    
    # 1. Mitochondrial genes
    mt_genes = gene_names.str.startswith('MT-') | gene_names.str.startswith('mt-')
    n_mt = mt_genes.sum()
    if verbose:
        print(f"\nMitochondrial genes (MT-*): {n_mt}")
        if n_mt > 0:
            mt_examples = gene_names[mt_genes].head(10).tolist()
            print(f"  Examples: {mt_examples}")
    keep_genes &= ~mt_genes
    
    # 2. Ribosomal genes
    ribo_genes = (gene_names.str.startswith('RPS') | 
                  gene_names.str.startswith('RPL') |
                  gene_names.str.startswith('Rps') | 
                  gene_names.str.startswith('Rpl'))
    n_ribo = ribo_genes.sum()
    if verbose:
        print(f"\nRibosomal genes (RPS*, RPL*): {n_ribo}")
        if n_ribo > 0:
            ribo_examples = gene_names[ribo_genes].head(10).tolist()
            print(f"  Examples: {ribo_examples}")
    keep_genes &= ~ribo_genes
    
    # 3. Long non-coding RNAs (lncRNAs)
    lnc_patterns = [
        r'LINC\d+',  # LINC genes
        r'.*-AS\d*$',  # Antisense genes
        r'.*-IT\d*$',  # Intronic transcripts
        r'.*-OT\d*$',  # Overlapping transcripts
        r'LOC\d+',  # LOC genes
        r'FLJ\d+',  # FLJ genes
        r'KIAA\d+',  # KIAA genes
        r'C\d+orf\d+',  # Chromosome open reading frames
    ]
    
    lnc_genes = pd.Series(False, index=gene_names.index)
    for pattern in lnc_patterns:
        pattern_match = gene_names.str.contains(pattern, regex=True, case=False)
        lnc_genes |= pattern_match
    
    n_lnc = lnc_genes.sum()
    if verbose:
        print(f"\nLong non-coding RNAs: {n_lnc}")
        if n_lnc > 0:
            lnc_examples = gene_names[lnc_genes].head(10).tolist()
            print(f"  Examples: {lnc_examples}")
    keep_genes &= ~lnc_genes
    
    # 4. microRNAs and small RNAs
    small_rna_genes = (gene_names.str.startswith('MIR') |
                       gene_names.str.startswith('SNOR') |
                       gene_names.str.startswith('SCARN') |
                       gene_names.str.startswith('U\d') |
                       gene_names.str.contains(r'^MIR\d+', regex=True) |
                       gene_names.str.contains(r'^SNORD\d+', regex=True))
    
    n_small_rna = small_rna_genes.sum()
    if verbose:
        print(f"\nSmall RNAs (miRNA, snoRNA, etc.): {n_small_rna}")
        if n_small_rna > 0:
            small_rna_examples = gene_names[small_rna_genes].head(10).tolist()
            print(f"  Examples: {small_rna_examples}")
    keep_genes &= ~small_rna_genes
    
    # 5. Pseudogenes
    pseudo_genes = gene_names.str.contains('P\d+$', regex=True) | gene_names.str.endswith('P')
    n_pseudo = pseudo_genes.sum()
    if verbose:
        print(f"\nPseudogenes: {n_pseudo}")
        if n_pseudo > 0:
            pseudo_examples = gene_names[pseudo_genes].head(10).tolist()
            print(f"  Examples: {pseudo_examples}")
    keep_genes &= ~pseudo_genes
    
    # Apply filtering
    n_removed = (~keep_genes).sum()
    n_remaining = keep_genes.sum()

    if verbose:
        print(f"\n=== Filtering Summary ===")
        print(f"Original genes: {adata.n_vars:,}")
        print(f"Genes to remove: {n_removed:,}")
        print(f"Remaining genes: {n_remaining:,}")
        print(f"Percentage removed: {(n_removed/adata.n_vars)*100:.1f}%")
        
        # Breakdown by category
        categories = {
            'Mitochondrial': mt_genes.sum(),
            'Ribosomal': ribo_genes.sum(),
            'Long non-coding': lnc_genes.sum(),
            'Small RNAs': small_rna_genes.sum(),
            'Pseudogenes': pseudo_genes.sum()
        }
        
        print(f"\nGenes removed by category:")
        for category, count in categories.items():
            if count > 0:
                print(f"  {category}: {count:,}")
    
    # Store information about removed genes
    adata.var['gene_category'] = 'protein_coding'
    adata.var.loc[mt_genes[mt_genes].index, 'gene_category'] = 'mitochondrial'
    adata.var.loc[ribo_genes[ribo_genes].index, 'gene_category'] = 'ribosomal'
    adata.var.loc[lnc_genes[lnc_genes].index, 'gene_category'] = 'lncRNA'
    adata.var.loc[small_rna_genes[small_rna_genes].index, 'gene_category'] = 'small_RNA'
    adata.var.loc[pseudo_genes[pseudo_genes].index, 'gene_category'] = 'pseudogene'
    
    # Filter the data
    adata_filtered = adata[:, keep_genes].copy()
    
    if verbose:
        print(f"\nFiltered data shape: {adata_filtered.shape}")
    
    if copy:
        return adata_filtered
    else:
        # Update original object
        adata._inplace_subset_var(keep_genes)
        print(f"✓ Filtering applied in place")