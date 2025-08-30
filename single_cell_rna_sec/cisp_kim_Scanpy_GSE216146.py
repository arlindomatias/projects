# Import packages
import os
import scanpy as sc # pp: preprocessing; tl: tools; pl: plotting
import seaborn as sns
import numpy as np
import gseapy as gp
import pandas as pd
from anndata import AnnData

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
sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
sc.tl.pca(adata)

sc.pp.neighbors(adata) ## Neighbors

## Leiden (3 resolutions)
for res in [0.15, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    ) # "leiden_res_0.15", "leiden_res_0.50", "leiden_res_2.00"

sc.tl.umap(adata) ## UMAP

sc.tl.tsne(adata, n_pcs=20)## t-SNE

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
## Rank genes
sc.tl.rank_genes_groups(adata, 'leiden')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers_df = pd.DataFrame(markers) # Clusters' differential expression
markers_scvi = model.differential_expression(groupby = 'leiden')
markers_scvi



markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > .5)]
markers_scvi

sc.pl.umap(adata, color = ['leiden'], frameon = False, legend_loc = "on data")



sc.pl.umap(adata, color = ['EPCAM', 'MUC1'], frameon = False, layer = 'scvi_normalized', vmax = 5)
#, layer = 'scvi_normalized'



cell_type = {"0":"Macrophage",
"1":"Fibroblast",
"2":"CD4+ T-cell",
"3":"AT2",
"4":"AT1",
"5":"CD8+ T-cell",
"6":"Endothelial cell",
"7":"Plasma cell",
"8":"Macrophage",
"9":"AT2",
"10":"Fibroblast",
"11":"Fibroblast",
"12":"Macrophage",
"13":"Macrophage",
"14":"Airway epithelial",
"15":"Airway epithelial",
"16":"Monocyte",
"17":"Airway epithelial",
"18":"B-cell",
"19":"Aerocyte",
"20":"Airway epithelial",
"21":"Smooth muscle cell",
"22":"Cycling T/NK",
"23":"Neuronal cell",
"24":"Denditic cell",
"25":"Pericyte",
"26":"Fibroblast",
"27":"Erythroid-like",
"28":"Macrophage"
}



adata.obs['cell type'] = adata.obs.leiden.map(cell_type)

sc.pl.umap(adata, color = ['cell type'], frameon = False)

adata.uns['scvi_markers'] = markers_scvi
adata.uns['markers'] = markers

adata.write_h5ad('integrated.h5ad')

model.save('model.model')


# Analysis
adata = sc.read_h5ad('integrated.h5ad')

adata.obs.Sample.unique().tolist()



def map_condition(x):
    if 'cov' in x:
        return 'COVID19'
    else:
        return 'control'

adata.obs['condition'] = adata.obs.Sample.map(map_condition)
adata.obs



num_tot_cells = adata.obs.groupby(['Sample']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
num_tot_cells

cell_type_counts = adata.obs.groupby(['Sample', 'condition', 'cell type']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
cell_type_counts

cell_type_counts['total_cells'] = cell_type_counts.Sample.map(num_tot_cells).astype(int)

cell_type_counts['frequency'] = cell_type_counts.doublet / cell_type_counts.total_cells

cell_type_counts



import matplotlib.pyplot as plt

plt.figure(figsize = (10,4))

ax = sns.boxplot(data = cell_type_counts, x = 'cell type', y = 'frequency', hue = 'condition')

plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')

plt.show()

# Differential Expression
subset = adata[adata.obs['cell type'].isin(['AT1', 'AT2'])].copy()

#two options: SCVI or diffxpy

import diffxpy.api as de

subset.X = subset.X.toarray()

len(subset.var)


subset

sc.pp.filter_genes(subset, min_cells=100)

len(subset.var)

subset.obs = subset.obs.rename(columns = {'cell type':'cell_type'})

#subset = subset.raw.to_adata() #need to run this if you scaled/regress your data and have negative numbers

subset.obs

#if want to test between covid/non covid
# res = de.test.wald(data=subset,
#              formula_loc= '~ 1 + condition',
#              factor_loc_totest='condition'
#                   )


res = de.test.wald(data=subset,
             formula_loc= '~ 1 + cell_type',
             factor_loc_totest='cell_type'
                  )

dedf = res.summary().sort_values('log2fc', ascending = False).reset_index(drop = True)
dedf



subset.obs.cell_type.unique()

dedf['log2fc'] = dedf['log2fc']*-1
dedf = dedf.sort_values('log2fc', ascending = False).reset_index(drop = True)
dedf

dedf = dedf[(dedf.qval < 0.05) & (abs(dedf.log2fc) > .5)]
dedf

dedf = dedf[dedf['mean'] > 0.15]
dedf

genes_to_show = dedf[-25:].gene.tolist() + dedf[:25].gene.tolist() #top 25 and bottom 25 from sorted df

sc.pl.heatmap(subset, genes_to_show, groupby='cell_type', swap_axes=True)

#DE with scvi

model  = scvi.model.SCVI.load('model.model', adata)



model

scvi_de = model.differential_expression(
    idx1 = [adata.obs['cell type'] == 'AT1'],
    idx2 = [adata.obs['cell type'] == 'AT2']
    )

#any set of cells vs any set of cells
# scvi_de = model.differential_expression(
#     idx1 = [(adata.obs['cell type'].isin(['AT1', 'AT2'])) & (adata.obs.condition == 'COVID19')],
#     idx2 = [(adata.obs['cell type'].isin(['AT1', 'AT2'])) & (adata.obs.condition == 'control')]
#     )

scvi_de

scvi_de = scvi_de[(scvi_de['is_de_fdr_0.05']) & (abs(scvi_de.lfc_mean) > .5)]
scvi_de = scvi_de.sort_values('lfc_mean')
scvi_de



scvi_de = scvi_de[(scvi_de.raw_normalized_mean1 > .5) | (scvi_de.raw_normalized_mean2 > .5)]
scvi_de



genes_to_show = scvi_de[-25:].index.tolist() + scvi_de[:25].index.tolist() #top 25 and bottom 25 from sorted df

sc.pl.heatmap(subset, genes_to_show, groupby='cell_type', swap_axes=True, layer = 'scvi_normalized',
              log = True)

# GO Enrichment
gp.get_library_name()
# 'GO_Biological_Process_2021',
#'KEGG_2021_Human',

subset

enr = gp.enrichr(gene_list= dedf[dedf.log2fc > 0].gene.tolist(),
                 gene_sets=['KEGG_2021_Human','GO_Biological_Process_2021'],
                 organism='human', # don't forget to set organism to the one you desired!
                 outdir=None, # don't write to disk,
                 background = subset.var_names.tolist()
                )



enr.results

# Comparisons
sc.pl.violin(subset[subset.obs.cell_type == 'AT2'], 'ETV5', groupby='condition')

temp = subset[subset.obs.cell_type == 'AT2']

i = np.where(temp.var_names == 'ETV5')[0][0]

a = temp[temp.obs.condition == 'COVID19'].X[:,i]
b = temp[temp.obs.condition == 'control'].X[:,i]

stats.mannwhitneyu(a, b)

# Score gene signature



#gene signature, ie, input list of genes from user
with open('datp_sig.txt') as f:
    datp_sig = [x.strip() for x in list(f)]

sc.tl.score_genes(subset, datp_sig, score_name = 'datp')
subset.obs
sc.pl.violin(subset, 'datp', groupby='condition')



a = subset[subset.obs.condition == 'COVID19'].obs.datp.values
b = subset[subset.obs.condition == 'control'].obs.datp.values
stats.mannwhitneyu(a, b)

sc.pl.umap(subset, color = 'datp', vmax = 1)

#for thumbnail

adata
with rc_context({'figure.figsize': (8,8)}):
    sc.pl.umap(adata, color = ['cell type'], frameon = False, s = 5, legend_loc = 'on data',
              legend_fontsize=12, legend_fontoutline=2)



with rc_context({'figure.figsize': (8,8)}):
    sc.pl.umap(adata, color = ['MUC1'], frameon = False, layer = 'scvi_normalized', vmax = 5, s = 5)

