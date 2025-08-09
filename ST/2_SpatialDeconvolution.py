## Spatial Deconvolution of Spatial Transcriptomics Data

#%% Step 1: Import Libraries and define path
import os
import sys
import gc
from xml.sax.handler import property_interning_dict
from projects.Kidney_HE.MedSAM.train_one_gpu import show_box
import torch
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

from Deconvolution_functions import *

# create paths and names to results folders for reference regression and cell2location models
results_folder = '/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/analysis'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

if not os.path.exists(ref_run_name):
    os.makedirs(ref_run_name)
if not os.path.exists(run_name):
    os.makedirs(run_name)

#%% Step 2: Load and Preprocessing scRNA-seq data
## Load data
all_result_path_ = "/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA"
adata_ref = sc.read_h5ad(f"{all_result_path_}/final_integrated.h5ad")
adata_ref = sc.AnnData(
    X=adata_ref.raw.X,
    obs=adata_ref.obs.copy(),  # Keep current obs (cell metadata)
    var=adata_ref.raw.var.copy()  # Use raw var (gene metadata)
)

major_anno = sc.read_h5ad(f"/mnt/NAS_21T/ProjectData/IMC_CRLM/scRNA/major_anno_all.h5ad")

nonmalignant_result_path_ = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/scRNA/nonMalignant_Analysis/Nonmalignant_annotated.h5ad"
nonmalignant_cells = sc.read_h5ad(nonmalignant_result_path_)
nonmalignant_cells = nonmalignant_cells[nonmalignant_cells.obs['Sub_type'] != "Unknown"] ## Remove 'Unknown' cell type

malignant_result_path_ = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","Malignant_Analysis")
malignant_cells = sc.read_h5ad(os.path.join(malignant_result_path_, "malignant_epithelial_cells_annotated.h5ad"))
malignant_cells.obs['Malignant_type'] = malignant_cells.obs['Malignant_type'].str.replace('EC_', 'TC_')
malignant_cells = malignant_cells[malignant_cells.obs['Malignant_type'] != "Unknown"] ## Remove 'Unknown' cell type

## Combine major cell types
adata_ref.obs["Major_type"] = "Unknown"
adata_ref.obs.loc[major_anno.obs.index, 'Major_type'] = major_anno.obs['Major_type']
adata_ref.obs["Major_type"].value_counts()

## Combine sub cell types
adata_ref.obs["Sub_type"] = adata_ref.obs["Major_type"]
adata_ref.obs.loc[malignant_cells.obs.index, 'Sub_type'] = malignant_cells.obs['Malignant_type']
adata_ref.obs.loc[nonmalignant_cells.obs.index, 'Sub_type'] = nonmalignant_cells.obs['Sub_type']
adata_ref.obs["Sub_type"].value_counts()
adata_ref = adata_ref[~adata_ref.obs['Sub_type'].isin(["Plasma","T","B","NK","Endothelial","Unknown"])] ## Remove 'Unknown' cell type

adata_ref_back = adata_ref.copy()

## Downsample cells
adata_ref, summary = adaptive_celltype_sampling(
    adata_ref, celltype_col='Sub_type', max_cells=1200
)
plot_sampling_summary(summary)

del adata_ref.uns['sampling_info']  # Remove neighbors to avoid memory issues
adata_ref.write_h5ad(f"{results_folder}/downsampled_reference.h5ad") ## Save downsampled data

## Preprocess data
from cell2location.utils.filtering import filter_genes

## Load subset reference data
if "adata_ref" not in globals():
    adata_ref = sc.read_h5ad(f"{results_folder}/downsampled_reference.h5ad")

## Filter mitochondrial genes
adata_ref = filter_mitochondrial_genes(adata=adata_ref)

## Select genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.1, nonz_mean_cutoff=1.5)
adata_ref = adata_ref[:, selected].copy() # filter the object
adata_ref.X = adata_ref.X.astype(int)

#%% Step 3: Estimation of reference cell type signatures (NB regression)
# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='batch',
                        # cell type, covariate used for constructing signatures
                        labels_key='Sub_type'
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        # categorical_covariate_keys=['Method']
                       )

# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

# Train the regression model
mod.train(max_epochs=300,batch_size=2500)
mod.plot_history(20)

# Export the estimated cell abundance (summary of the posterior distribution).
# adata_ref = mod.export_posterior(
#     adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
# )
adata_ref = mod.export_posterior( # quantiles of the posterior distribution
    adata_ref, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05","q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)
# mod.plot_QC(summary_name="q05")

# Save anndata object with results
adata_ref.write(f"{ref_run_name}/sc.h5ad")

# export estimated expression in each cluster
adata_ref = sc.read_h5ad(f"{ref_run_name}/sc.h5ad")

name_of_per_cluster_mu_fg = "q05_per_cluster_mu_fg"
if name_of_per_cluster_mu_fg in adata_ref.varm.keys():
    inf_aver = adata_ref.varm[name_of_per_cluster_mu_fg][[f'{name_of_per_cluster_mu_fg}_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
    
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

#%% Step 4: Cell2location: spatial mapping

## Load 10x Visium data
rawST_datapath = "/mnt/public/lyx/IMC_HE_Merge/CRLM/figures/ST/raw"

slides = os.listdir(rawST_datapath)

for slide_ in slides:
    slide_name = slide_.split("_raw")[0] ## Get slide name
    adata_vis = sc.read_h5ad(os.path.join(rawST_datapath,slide_))
    adata_vis.obs["sample"] = slide_name

    ## Filter genes
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    if slide_name == "FDZS_A04932G3_bin50":

        tmas = adata_vis.obs['tma_core'].unique()

        for tma_ in tmas:
            adata_vis_ = adata_vis[adata_vis.obs['tma_core'] == tma_].copy()
            adata_vis_.X = adata_vis_.X.astype(int)
            cell2location.models.Cell2location.setup_anndata(adata_vis_, batch_key="sample")
            
            # create and train the model
            mod = cell2location.models.Cell2location(
                adata_vis_, cell_state_df=inf_aver,
                # the expected average cell abundance: tissue-dependent
                # hyper-prior which can be estimated from paired histology:
                N_cells_per_location=10, # 5 for Visium
                # hyperparameter controlling normalisation of
                # within-experiment variation in RNA detection:
                detection_alpha=20
                # refer: https://github.com/BayraktarLab/cell2location/blob/master/docs/images/Note_on_selecting_hyperparameters.pdf
            )
            mod.view_anndata_setup()

            # Training cell2location
            mod.train(max_epochs=30000,
                # train using full data (batch_size=None)
                batch_size=None,
                # use all data points in training because
                # we need to estimate cell abundance at all locations
                train_size=1
                )
            
            # export estimated cell abundance (absolute number of cells per location)
            adata_vis_ = mod.export_posterior(
                adata_vis_, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
            )

            # Save model and anndata
            mod.save(f"{run_name}/{slide_name}_{tma_}/", overwrite=True)
            adata_vis_.write(f"{run_name}/{slide_name}_{tma_}/sp.h5ad")
            
            # Load model and anndata for visualization
            # adata_vis = sc.read_h5ad(f"{run_name}/{slide_name}/sp.h5ad")
            # mod = cell2location.models.Cell2location.load(f"{run_name}/{slide_name}", adata_vis)

            # ===== MEMORY CLEANUP =====
            del mod, adata_vis_    
            gc.collect()

            torch.cuda.empty_cache()
            torch.cuda.synchronize()        
    else:
        cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

        # create and train the model
        mod = cell2location.models.Cell2location(
            adata_vis, cell_state_df=inf_aver,
            # the expected average cell abundance: tissue-dependent
            # hyper-prior which can be estimated from paired histology:
            N_cells_per_location=5, # 5 for Visium
            # hyperparameter controlling normalisation of
            # within-experiment variation in RNA detection:
            detection_alpha=20
            # refer: https://github.com/BayraktarLab/cell2location/blob/master/docs/images/Note_on_selecting_hyperparameters.pdf
        )
        mod.view_anndata_setup()

        # Training cell2location
        mod.train(max_epochs=30000,
            # train using full data (batch_size=None)
            batch_size=None,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1
            )

        # plot ELBO loss history during training, removing first 100 epochs from the plot
        # mod.plot_history(1000)
        # plt.legend(labels=['full data training'])

        # export estimated cell abundance (absolute number of cells per location)
        adata_vis = mod.export_posterior(
            adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
        )

        # Save model and anndata
        mod.save(f"{run_name}/{slide_name}/", overwrite=True)
        adata_vis.write(f"{run_name}/{slide_name}/sp.h5ad")
        
        # Load model and anndata for visualization
        # adata_vis = sc.read_h5ad(f"{run_name}/{slide_name}/sp.h5ad")
        # mod = cell2location.models.Cell2location.load(f"{run_name}/{slide_name}", adata_vis)

        # ===== MEMORY CLEANUP =====
        del mod, adata_vis    
        gc.collect()

        torch.cuda.empty_cache()
        torch.cuda.synchronize()

#%% Step 5: Visualising cell abundance in spatial coordinates

slides = os.listdir(run_name)

for slide_ in slides:

    print(f"Processing slide: {slide_}")

    ## Load slides
    slide_path = os.path.join(run_name, slide_,"sp.h5ad")
    adata = sc.read_h5ad(slide_path)

    # to adata.obs with nice names for plotting
    adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
    print(adata.obs.columns)

    # plot in spatial coordinates
    with mpl.rc_context({'axes.facecolor':  'black',
                        'figure.figsize': [4.5, 5]}):

        sc.pl.spatial(adata, cmap='magma',
                    # show first 8 cell types
                    color=['Hepatocyte', 'TC_Glycolysis', 'TC_Quiescent', 'TC_EMT',
                            'CD8T_GZMK', 'CAF_POSTN', 'EC_LYVE1', 'Fibro_COL3A1'],
                    ncols=4, size=1.3,
                    img_key='hires',
                    # limit color scale at 99.2% quantile of cell abundance
                    vmin=0, vmax='p99.2',show = False
                    )
        