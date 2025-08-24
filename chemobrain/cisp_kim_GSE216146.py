# Pacotes
import os
import scanpy as sc
import pandas as pd
import seaborn as sns
import scvi
import numpy as np
import gc
import dtale
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import gseapy as gp #this method requires internet connection
from scipy import stats
from matplotlib.pyplot import rc_context




# Verificar diretório de trabalho
os.chdir("/home/arlindo/R/wd")
print(os.getcwd())

# Converter df em panda
pdf = pd.DataFrame
# dtale.show(df) abrir o df em uma janela
dt = dtale.show

# Atalho para limpar a memória
clear = gc.collect
clear()

# Importar dados
raw = sc.read_h5ad("GSE216146_chemo_brain.h5ad", backed='r') # Importa arquivo bruto em modo de leitura
raw.obs['group_batch'] = raw.obs['pathology'].astype(str) + '_' + raw.obs['batch'].astype(str) # Subset de batch e tratamento
raw.obs['biosample'].unique() # Ver as amostras biológicas

# Identificar subsets isoladamente
used_samples_ids = raw.obs['pathology'].isin(['PN', 'CN']) # Identifica apenas as amostras de interesse
veh_ids = raw.obs['pathology'].isin(['PN']) # IDs veh
chem_ids = raw.obs['pathology'].isin(['CN']) # IDs chem
b1_ids = raw.obs['batch'].isin(['D20']) # IDs Bateria 1
b2_ids = raw.obs['batch'].isin(['D21']) # Ids Bateria 2

# Dados usados
adata = raw[used_samples_ids, :] # Apenas os dados de interesse
adata.obs['group_batch'].unique() # verifica se a tag de batch e tratamento está inclusa
adata.obs['group_batch'].value_counts() # mostra a quantidade de células em cada condição

clear()

# Doublet remove
veh_data = raw[veh_ids, :].to_memory() # até aqui não travou
sc.pp.filter_genes(veh_data, min_cells = 10) # Filtrar apenas genes presentes em n células ou mais
sc.pp.highly_variable_genes(veh_data, n_top_genes = 1000, subset = True, flavor = 'seurat_v3') # Filtrar apenas os n genes mais representativos
scvi.model.SCVI.setup_anndata(veh_data)
vae = scvi.model.SCVI(veh_data, n_hidden=64, n_latent=10, n_layers=1)
vae.train(max_epochs=200, batch_size=256)
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train(max_epochs=200, batch_size=256)
df = solo.predict()
df['prediction'] = solo.predict(soft = False)
#df.index = df.index.map(lambda x: x[:-9])
df.groupby('prediction').count()
df['dif'] = df.doublet - df.singlet
dt(df)
sns.displot(df[df.prediction == 'doublet'], x = 'dif')
doublet = df[(df.prediction == 'doublet') & (df.dif > 0.4)]
dt(doublet)

clear()

veh_data.obs['doublet'] = veh_data.obs.index.isin(doublet.index)
dt(veh_data.obs)
veh_data = veh_data[~veh_data.obs.doublet]
veh_data

# Preprocessing
veh_data.var
veh_data.var['mt'] = veh_data.var.index.str.startswith('mt-')
dt(veh_data.var)

ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
ribo_genes

veh_data.var['ribo'] = veh_data.var_names.isin(ribo_genes[0].values)
sc.pp.calculate_qc_metrics(veh_data, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
veh_data.var.sort_values('n_cells_by_counts')
sc.pp.filter_genes(veh_data, min_cells=3)
veh_data.obs.sort_values('n_genes_by_counts')
veh_data.obs.sort_values('total_counts')
sc.pp.filter_cells(veh_data, min_genes=200)
sc.pl.violin(veh_data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)

veh_df = pdf(veh_data.obs)

upper_lim = np.quantile(veh_data.obs.n_genes_by_counts.values, .98)
upper_lim = 350

veh_data = veh_data[veh_data.obs.n_genes_by_counts < upper_lim]
veh_data.obs

veh_data = veh_data[veh_data.obs.pct_counts_mt < 20]
veh_data = veh_data[veh_data.obs.pct_counts_ribo < 2]
veh_data

# Normalization
veh_data.X.sum(axis = 1)
sc.pp.normalize_total(veh_data, target_sum=1e4) #normalize every cell to 10,000 UMI
veh_data.X.sum(axis = 1)

sc.pp.log1p(veh_data) #change to log counts
dt(veh_data.obs)
veh_data.X.sum(axis = 1)
veh_data.raw = veh_data

# Clustering
sc.pp.highly_variable_genes(veh_data, n_top_genes = 1000)
veh_data.var
sc.pl.highly_variable_genes(veh_data)
veh_data = veh_data[:, veh_data.var.highly_variable]

sc.pp.regress_out(veh_data, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])
sc.pp.scale(veh_data, max_value=10)
sc.tl.pca(veh_data, svd_solver='arpack')
sc.pl.pca_variance_ratio(veh_data, log=True, n_pcs = 50)

sc.pp.neighbors(veh_data, n_pcs = 30)
sc.tl.umap(veh_data)
sc.pl.umap(veh_data)

sc.tl.leiden(veh_data, resolution = 0.5)
veh_data.obs

sc.pl.umap(veh_data, color=['leiden'])

################# Integration #################
def pp(csv_path):
    adata = sc.read_csv(csv_path).T
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat_v3')
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    df['prediction'] = solo.predict(soft=False)
    df.index = df.index.map(lambda x: x[:-2])
    df['dif'] = df.doublet - df.singlet
    doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]

    adata = sc.read_csv(csv_path).T
    adata.obs['Sample'] = csv_path.split('_')[2]  # 'raw_counts/GSM5226574_C51ctr_raw_counts.csv'

    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]

    sc.pp.filter_cells(adata, min_genes=200)  # get rid of cells with fewer than 200 genes
    # sc.pp.filter_genes(adata, min_cells=3) #get rid of genes that are found in fewer than 3 cells
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 2]

    return adata

out = []
for file in os.listdir('raw_counts/'):
    out.append(pp('raw_counts/' + file))

adata = sc.concat(out)

sc.pp.filter_genes(adata, min_cells = 10)

adata.X = csr_matrix(adata.X)

adata.write_h5ad('combined.h5ad')
adata = sc.read_h5ad('combined.h5ad')

adata.obs.groupby('Sample').count()



# sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset = True, layer = 'counts',
#                            flavor = "seurat_v3", batch_key="Sample") #no batch_key if one sample


scvi.model.SCVI.setup_anndata(adata, layer = "counts",
                             categorical_covariate_keys=["Sample"],
                             continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])



model = scvi.model.SCVI(adata)

model.train() #may take a while without GPU
adata.obsm['X_scVI'] = model.get_latent_representation()

adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)
sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.5)
sc.pl.umap(adata, color = ['leiden', 'Sample'], frameon = False)

adata.write_h5ad('integrated.h5ad')

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

