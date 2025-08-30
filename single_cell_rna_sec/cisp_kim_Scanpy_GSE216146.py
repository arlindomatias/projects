# Importar pacotes e funções
import os
import scanpy as sc
import pandas as pd
import seaborn as sns
import scvi
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp

from dtale import show as dt
from scipy.sparse import csr_matrix, issparse
from pandas import DataFrame as pdf
from gc import collect as clear
from scipy import stats
from matplotlib.pyplot import rc_context
from celltypist import models, annotate

# Ajustes
sc.settings.set_figure_params(dpi_save = 300, facecolor = "white", format = "png") # Resolução
os.chdir("/home/arlindo/R/wd/SC") # Diretório de trabalho
clear()

# Importar dados e isolar amostra
raw = sc.read_h5ad("GSE216146_chemo_brain.h5ad", backed='r', as_sparse=()) # Importa arquivo bruto em modo de leitura
print(raw[raw.obs['pathology'] == 'PN'].obs['biosample']) # Verificar IDs controle
biosample = 'D21-2746' # Definir a amostra
sample_id = raw.obs['biosample'].isin([biosample]) # ID selecionado
sample_data =  sc.AnnData(X=raw[sample_id, :].raw.X.copy(),
                         obs=raw[sample_id, :].obs.copy(),
                         var=raw[sample_id, :].raw.var.copy())

sample_data.write_h5ad(f"{biosample}_data.h5ad") # Salvar arquivo

df_var = pdf(sample_data.var)  # info dos genes
df_obs = pdf(sample_data.obs)  # info das células
counts_dense = sample_data.raw.X.toarray() if hasattr(sample_data.raw.X, "toarray") else sample_data.raw.X
counts_df = pd.DataFrame(counts_dense, index=sample_data.obs_names, columns=sample_data.var_names)

clear()

##################################
# Integration
## Identify interest samples
samples = raw.obs['biosample'].unique()
print(raw[raw.obs['pathology'] == 'PN'].obs['biosample']) # Grupo controle
print(raw[raw.obs['pathology'] == 'CN'].obs['biosample']) # Grupo chem

## Select
samples = ['D20-6407','D21-2746','D21-2747','D20-6409','D21-2750','D21-2751'] # Definir amostra

## Create importing function
def pp(biosample):
    # Identificar células e importar dados
    data_id = raw.obs['biosample'].isin([biosample])
    adata = sc.AnnData(X=raw[data_id, :].raw.X.copy(),
                         obs=raw[data_id, :].obs.copy(),
                         var=raw[data_id, :].raw.var.copy())

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

# Preprocessing
## Quality Control
### Mitochondrial genes
adata.var['mt'] = adata.var.index.str.startswith('Mt-')
### Ribosomal genes
adata.var['ribo'] = adata.var.index.str.startswith(('Rps', 'Rpl'))
### Hemoglobin genes
adata.var['hb'] = adata.var.index.str.contains("^Hb[ab]-")

### Calculate metrics and filter data
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=3)
adata.obs.sort_values('total_counts')
sc.pp.filter_cells(adata, min_genes=100)
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
upper_lim = 4500
adata = adata[adata.obs.n_genes_by_counts < upper_lim]
adata = adata[adata.obs.pct_counts_mt < 20]
adata = adata[adata.obs.pct_counts_ribo < 20]

### Visualize cleaned data
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_ribo")# Save changes to a new file
#adata.X = csr_matrix(adata.X)
adata.obs.groupby('biosample')
adata.write_h5ad('combined.h5ad')
adata = sc.read_h5ad('combined.h5ad')

# Feature Selection
sc.pp.highly_variable_genes(adata, n_top_genes = 2000, batch_key="batch") # Filtrar apenas os n genes mais representativos
sc.pl.highly_variable_genes(adata)

## Dimensionality Reduction
### PCA
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca(
    adata,
    color=['n_genes_by_counts', 'n_genes_by_counts', 'total_counts','total_counts'],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2)

### Neighbors
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color=["total_counts", "pathology"],
    # Setting a smaller point size to get prevent overlap
    size=2)

## Clustering
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])

## Re-assess QC
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3)

sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    size=2)

## Cell-type annotation
### Colors
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(
        sample_data, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data", size=2)
marker_genes = {
    "Neuron": ["Snap25", "Rbfox3", "Syt1"],
    "Astrocyte": ["Gfap", "Aldh1l1", "Slc1a3"],
    "Microglia": ["Cx3cr1", "P2ry12", "Tmem119"],
    "Oligodendrocyte": ["Mog", "Plp1", "Mag"],
    "Endothelial": ["Pecam1", "Cldn5", "Flt1"]
}
### Manual Annotation
sc.pl.dotplot(sample_data, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

### Celltypist Annotation
models.download_models(force_update = True)
models.models_description()
model = models.Model.load('Mouse_Isocortex_Hippocampus.pkl')  # Mouse
model

predictions = annotate(sample_data, model=model, majority_voting=True)
pred_labels_df = predictions.predicted_labels.copy()
sample_data.obs['cell_type_celltypist'] = pred_labels_df.iloc[:, 2].values
sc.pl.umap(sample_data, color='cell_type_celltypist', size=20, palette='tab20')
print(sample_data.obs['cell_type_celltypist'].value_counts()) # Conferir contagem por tipo



# Label Cells Types
sc.tl.leiden(adata, resolution = 1)
sc.tl.rank_genes_groups(adata, 'leiden')



#sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers

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

