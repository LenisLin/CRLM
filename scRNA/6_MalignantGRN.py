## GRN analysis of malignant epithelial cells using pySCENIC

import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose
from collections import Counter
from IPython.display import HTML, display

import scipy.io
import math
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss
from pyscenic.cli.utils import load_signatures

sc.settings.njobs = 12 # Set maximum number of jobs for Scanpy.

figurePath = os.path.join("/mnt/public/lyx/IMC_HE_Merge/CRLM/figures", "scRNA","Malignant_Analysis")
result_path = figurePath
save_path = os.path.join(figurePath,"GRN_analysis")
os.makedirs(save_path, exist_ok=True)

# Define Folder structure.
RESOURCES_FOLDERNAME = "/home/lenislin/Experiment/projects/zsz_CRC/pySCENIC/resources"
RESULTS_FOLDERNAME = save_path
FIGURES_FOLDERNAME = save_path
if not os.path.exists(RESULTS_FOLDERNAME):
    os.makedirs(RESULTS_FOLDERNAME,exist_ok=True)

# Auxilliary functions.
BASE_URL = save_path
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"

# Download Auxilliary data sets.
# Downloaded fromm pySCENIC github repo: https://github.com/aertslab/pySCENIC/tree/master/resources
HUMAN_TFS_FNAME = os.path.join(RESOURCES_FOLDERNAME, 'allTFs_hg38.txt')
# Ranking databases. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
RANKING_DBS_FNAMES = list(map(lambda fn: os.path.join(RESOURCES_FOLDERNAME, fn),
                       ['hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
                        'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather']))
# Motif annotations. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDERNAME, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

################### pyscenic ###################

if 'adata' not in locals():
    adata = sc.read_h5ad(os.path.join(result_path, "malignant_epithelial_cells_annotated.h5ad"))

## Results save
DATASET_ID="CRLM_Malignant"

METADATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.metadata.csv'.format(DATASET_ID))
EXP_MTX_QC_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.qc.tpm.csv'.format(DATASET_ID))
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(DATASET_ID))
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.motifs.csv'.format(DATASET_ID))
REGULONS_DAT_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.regulons.dat'.format(DATASET_ID))
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.auc.csv'.format(DATASET_ID))
BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.bin.csv'.format(DATASET_ID))
THR_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.thresholds.csv'.format(DATASET_ID))
ANNDATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.h5ad'.format(DATASET_ID))
LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.loom'.format(DATASET_ID))

adata.write_h5ad(ANNDATA_FNAME) # Categorical dtypes are created.
adata.to_df().to_csv(EXP_MTX_QC_FNAME)

## Load functions
def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    
    MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
    pd.set_option('display.max_colwidth', -1)
    display(HTML(df.head().to_html(escape=False)))
    pd.set_option('display.max_colwidth', MAX_COL_WIDTH)

def mapping_arr_to_matrix(exp_,colnames,meta_):
    # Create a mapping from labels to indices
    label_to_index = {label: idx for idx, label in enumerate(colnames)}

    # Find indices of the labels in characteristic_labels
    characteristic_labels = meta_.index
    col_indices = [label_to_index[label] for label in characteristic_labels if label in label_to_index]

    # Select the corresponding columns from the sparse matrix
    exp_ = exp_.tocsr()  # For efficient row slicing
    selected_columns = exp_[:, col_indices]

    # (Optional) Convert the result to a dense array for display
    # selected_columns_dense = selected_columns.toarray()
    
    return selected_columns

# REGULON CREATION
def derive_regulons(motifs, db_names=(#'hg19-tss-centered-10kb-10species', 'hg19-500bp-upstream-10species','hg19-tss-centered-5kb-10species'
        "hg38-limited-upstream500-tss-downstream100-full-transcript",
        "hg38-limited-upstream10000-tss-downstream10000-full-transcript"
                                      )):
    # motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]

    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]
    
    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0) 
                                                                      & ((motifs['Annotation'] == 'gene is directly annotated')
                                                                        | (motifs['Annotation'].str.startswith('gene is orthologous to')
                                                                           & motifs['Annotation'].str.endswith('which is directly annotated for motif')))
                                                                     ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))


def saveimg(fname: str, fig, folder: str=FIGURES_FOLDERNAME) -> None:
    """
    Save figure as vector-based PDF image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname))


def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

def clusterplot(meta, bin_mtx_, auc_mtx_, cellid_col = "cell_id", celltype_col = "sub_celltype", saveName = None):
        
    COLORS = [color['color'] for color in mpl.rcParams["axes.prop_cycle"]]

    cell_type_color_lut = dict(zip(meta[celltype_col].dtype.categories, COLORS))
    cell_id2cell_type_lut = meta.set_index(cellid_col)[celltype_col].to_dict()
    bw_palette = sns.xkcd_palette(["white", "black"])

    sns.set_style("whitegrid")
    fig = palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])
    saveimg(f'{saveName}-legend on_off.png', fig)

    sns.set(font_scale=0.8)
    fig = palplot(sns.color_palette(COLORS), meta[celltype_col].dtype.categories, size=1.0)
    saveimg(f'{saveName}-legend cell_type_colors.png', fig)

    sns.set(font_scale=1.0)
    sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})
    g = sns.clustermap(bin_mtx_.T, 
                col_colors=auc_mtx_.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
                cmap=bw_palette, figsize=(20,20))
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.set_xticks([])
    g.ax_heatmap.set_xlabel('Cells')
    g.ax_heatmap.set_ylabel('Regulons')
    g.ax_col_colors.set_yticks([0.5])
    g.ax_col_colors.set_yticklabels(['Cell Type'])
    g.cax.set_visible(False)
    g.savefig(os.path.join(FIGURES_FOLDERNAME, f'{saveName}-GRN clustermap.pdf'), format='pdf')

    return None

## STEP 1: Network inference based on GRNBoost2 from CLI
os.system(f"pyscenic grn {EXP_MTX_QC_FNAME} {HUMAN_TFS_FNAME} -o {ADJACENCIES_FNAME} --num_workers 8")

## STEP 2-3: Regulon prediction aka cisTarget from CLI
DBS_PARAM = ' '.join(RANKING_DBS_FNAMES)
os.system(f"pyscenic ctx {ADJACENCIES_FNAME} {DBS_PARAM} --annotations_fname {MOTIF_ANNOTATIONS_FNAME} --expression_mtx_fname {EXP_MTX_QC_FNAME} --output {MOTIFS_FNAME} --num_workers 12")

df_motifs = load_motifs(MOTIFS_FNAME)
df_motifs.head()

## STEP 4: Cellular enrichment aka AUCell (skip this step: https://github.com/aertslab/pySCENIC/issues/199)
# regulons = derive_regulons(df_motifs)
regulons = load_signatures(MOTIFS_FNAME)
with open(REGULONS_DAT_FNAME, 'wb') as f:
    pickle.dump(regulons, f)

with open(REGULONS_DAT_FNAME, 'rb') as f:
    regulons = pickle.load(f)

# AUCELL
df_tpm = pd.read_csv(EXP_MTX_QC_FNAME, index_col=0)
auc_mtx = aucell(df_tpm, regulons, num_workers=12) 
auc_mtx.to_csv(AUCELL_MTX_FNAME)

## OPTIONAL STEP 5 - Regulon activity binarization¶
auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0)

bin_mtx, thresholds = binarize(auc_mtx) 
bin_mtx.to_csv(BIN_MTX_FNAME) 
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME)

bin_mtx = pd.read_csv(BIN_MTX_FNAME, index_col=0)
thresholds = pd.read_csv(THR_FNAME, index_col=0).threshold

sort_auc_mtx = pd.DataFrame(data=np.sum(auc_mtx,axis=0))
sort_auc_mtx = sort_auc_mtx.sort_values(by=0, ascending=False)

sort_auc_mtx.to_csv(os.path.join(FIGURES_FOLDERNAME,'Sorted AUC of GRN in Malignant cells.csv'))
top_grn = [x for x in sort_auc_mtx.head(8).index]

fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(8, 4), dpi=100)
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

for i, grn in enumerate(top_grn):
    # Select the corresponding axis (ax1 to ax8)
    ax = axes[i]
    plot_binarization(auc_mtx, grn, thresholds[grn], ax=ax)

saveimg(os.path.join(FIGURES_FOLDERNAME,'GRN of Malignant cells.pdf'), fig)

## Create heatmap with binarized regulon activity
# adata = sc.read_h5ad(ANNDATA_FNAME)
# df_metadata = adata.obs
# df_metadata.head()
# df_metadata["cell_id"] = df_metadata.index

# meta = df_metadata[df_metadata['sub_celltype'].str.startswith(target_type, na=False)]
# meta['sub_celltype'] = meta['sub_celltype'].cat.remove_unused_categories()

# bin_mtx_ = bin_mtx.loc[meta.index,]
# auc_mtx_ = auc_mtx.loc[meta.index,]
# 
# clusterplot(meta=meta,bin_mtx_=bin_mtx_,auc_mtx_=auc_mtx_)

## STEP 6: Non-linear projection and clustering
adata = sc.read_h5ad(ANNDATA_FNAME)
auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0)
with open(REGULONS_DAT_FNAME, 'rb') as f:
    regulons = pickle.load(f)
f.close()

add_scenic_metadata(adata, auc_mtx, regulons)

rss = regulon_specificity_scores(auc_mtx, adata.obs.Malignant_type)
rss.head()
rss.to_csv(os.path.join(FIGURES_FOLDERNAME, "rss.csv"))

## CELL TYPE SPECIFIC REGULATORS - RSS
# List of features to plot
# features = [
#     'Mph_APOE', 'Mph_CCL20', 'Mph_S100A8', 'Mph_SPP1',
#     'Fibro_CCL11', 'Fibro_RGS5', 'Fibro_CXCL8', 'Fibro_MYH11'
# ]
features = [x for x in set(adata.obs["Malignant_type"])]
len_feat = math.ceil(len(features) / 2)

# Flatten the axes for easier indexing
sns.set(style='whitegrid', font_scale=1)
fig, axes = plt.subplots(2, len_feat, figsize=((len_feat + 0.5)*2, 8), dpi=300) # Create a 2x4 grid for 8 plots

axes = axes.flatten()

# Plot each feature
for ax, feature in zip(axes, features):
    plot_rss(rss, feature, ax=ax)
    ax.set_xlabel('')
    ax.set_ylabel('')

plt.tight_layout()
saveimg(os.path.join(FIGURES_FOLDERNAME,'GRN-RSS of Malignant cells.pdf'), fig)

## CELL TYPE SPECIFIC REGULATORS - Z-SCORE
df_obs = adata.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['Malignant_type']]
df_results = ((df_scores.groupby(by='Malignant_type').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
df_results[(df_results.Z >= 2.0)].sort_values('Z', ascending=False)

# Create pivot table with only top regulons
top_5_regulons = df_results.groupby('Malignant_type').apply(lambda x: x.nlargest(5, 'Z')).reset_index(drop=True)
df_heatmap_top5 = pd.pivot_table(data=top_5_regulons,index='Malignant_type', columns='regulon', values='Z')

# Plot heatmap
fig, ax1 = plt.subplots(1, 1, figsize=(12, 8))
sns.heatmap(df_heatmap_top5, ax=ax1, annot=True, fmt=".1f", linewidths=.7, 
            cbar=False, square=True, linecolor='gray', 
            cmap="YlGnBu", annot_kws={"size": 6})
ax1.set_ylabel('')
saveimg('Top5_Malignant_cell_regulons_heatmap.pdf', fig)
