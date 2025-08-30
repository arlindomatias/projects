# Import packages
import os
import scanpy as sc # pp: preprocessing; tl: tools; pl: plotting
import seaborn as sns
import numpy as np
import gseapy as gp
import pandas as pd
from anndata import AnnData
import scvi
import matplotlib.pyplot as plt

from gc import collect as clear
from scipy import stats
from matplotlib.pyplot import rc_context
from celltypist import models, annotate

# Global adjusts
sc.settings.set_figure_params(
    dpi_save = 300,
    facecolor = "white",
    format = "png") # Resolução
os.chdir("/home/arlindo/R/wd/SC") # Set working directory

# Import raw data
raw = sc.read_h5ad(
    "GSE216146_chemo_brain.h5ad",
    backed='r',
    as_sparse=()) # Importa arquivo bruto em modo de leitura
clear()

##################################
# Integration
## Identify interest samples
samples = raw.obs['biosample'].unique() # See all samples
print(raw[raw.obs['pathology'] == 'PN'].obs['biosample']) # Control samples
print(raw[raw.obs['pathology'] == 'CN'].obs['biosample']) # Chemo samples

## Select
samples = ['D20-6407','D21-2746','D21-2747','D20-6409','D21-2750','D21-2751'] # Definir amostra

## Create importing function
def pp(biosample):
    # Identificar células e importar dados
    data_id = raw.obs['biosample'].isin([biosample])
    adata = sc.AnnData(X=raw[data_id, :].raw.X.copy(),
                         obs=raw[data_id, :].obs.copy(),
                         var=raw[data_id, :].raw.var.copy())
    # Quality control
    adata.var['mt'] = adata.var.index.str.startswith('Mt-')
    adata.var['ribo'] = adata.var.index.str.startswith(('Rps', 'Rpl'))
    adata.var['hb'] = adata.var.index.str.contains("^Hb[ab]-")
    sc.pp.calculate_qc_metrics(adata,
                               qc_vars=['mt', 'ribo', 'hb'],
                               percent_top=None,
                               log1p=False,
                               inplace=True)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs.sort_values('total_counts')
    sc.pp.filter_cells(adata, min_genes=100)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 20]
    # Doublet detection
    doublets = adata.obs[adata.obs['predicted_doublet'] == True]
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]
    # Normalization
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Save
    adata.write_h5ad(f"{biosample}_data.h5ad")
    return adata

## Merge individual data
out = []
for sample in samples:
    out.append(pp(sample))
adata = sc.concat(out)

## Data evaluation
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_ribo")# Save changes to a new file

## Save progress
adata.obs.groupby('biosample')
adata.write_h5ad('combined.h5ad')
adata = sc.read_h5ad('combined.h5ad')
del raw
adata = adata.to_memory()

# Feature Selection
sc.pp.highly_variable_genes(
    adata,
    n_top_genes = 2000,
    batch_key="batch") # Filtrar apenas os n genes mais representativos

# Clustering and dimensionality reduction
## PCA
sc.pp.pca(adata, use_highly_variable=True)
sc.tl.pca(adata)

sc.pp.neighbors(adata) ## Neighbors

## Leiden (3 resolutions)
for res in [0.15, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    ) # "leiden_res_0.15", "leiden_res_0.50", "leiden_res_2.00"

sc.tl.umap(adata) ## UMAP

sc.tl.tsne(adata, n_pcs=20)## t-SNE
clear()

# Cell-type annotation
## Manual Annotation
marker_genes = {
    "Neuron": ["Snap25", "Rbfox3", "Syt1"],
    "Astrocyte": ["Gfap", "Aldh1l1", "Slc1a3"],
    "Microglia": ["Cx3cr1", "P2ry12", "Tmem119"],
    "Oligodendrocyte": ["Mog", "Plp1", "Mag"],
    "Endothelial": ["Pecam1", "Cldn5", "Flt1"]}

## Celltypist Annotation
### Prepare model
models.download_models()
models.models_description()
model = 'Mouse_Postnatal_DentateGyrus.pkl'

# Annotate
def cell_type(set = model, data = adata):
    model = models.Model.load(set)
    predictions = annotate(data, model=model, majority_voting=True)
    labels = predictions.predicted_labels
    data.obs['cell_type'] = labels['majority_voting'].values
    return data

adata = cell_type()
adata.obs['cell_type']
print(adata.obs['cell_type'].value_counts())

# Data visualization
sc.pl.highly_variable_genes(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca(
    adata,
    color=['pathology', 'batch'],
    dimensions=[(0, 1), (0, 1)],
    ncols=2,
    size=2,
)

sc.pl.umap(
    adata,
    color=['leiden_res_0.15', 'batch', 'pathology','cell_type'],
    wspace=0.5,
    # Setting a smaller point size to get prevent overlap
    size=2,
    ncols=2,
    #legend_loc="on data"
    )

sc.pl.tsne(adata,
           color=['leiden_res_0.15', 'batch', 'pathology', 'cell_type'],
           wspace=0.5,
           # Setting a smaller point size to get prevent overlap
           size=2,
           ncols=2,
           # legend_loc="on data"
           )

sc.pl.dotplot(adata,
              marker_genes,
              groupby=['pathology'],
              standard_scale="var",
              dendrogram=True)

with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata, color="Cdkn2a")

# Analysis
## Rank genes in clusters
sc.tl.rank_genes_groups(adata, 'leiden_res_0.15')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
markers = sc.get.rank_genes_groups_df(adata, None)
markers = pd.DataFrame(
    markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)])

## Counting cells
def map_condition(condition):
    if 'PN' in condition:
        return 'Control'
    else:
        return 'Chemotherapy'

adata.obs['pathology'] = adata.obs.pathology.map(map_condition)
num_tot_cells = adata.obs.groupby(['biosample']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
num_tot_cells

cell_type_counts = adata.obs.groupby(
    ['biosample', 'pathology', 'cell_type']).count()
cell_type_counts = pd.DataFrame(cell_type_counts)
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
cell_type_counts = cell_type_counts.rename(
    columns = {'initial_size_spliced': 'cell_type_counts'})

cell_type_counts['total_cells'] = cell_type_counts.biosample.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.cell_type_counts / cell_type_counts.total_cells

plt.ioff()
fig, ax = plt.subplots()
ax.plot = sns.boxplot(
    data = cell_type_counts,
    x = 'cell_type',
    y = 'frequency',
    hue = 'pathology')
plt.xticks(rotation = 45,
           rotation_mode = 'anchor',
           ha = 'right')
plt.show()

