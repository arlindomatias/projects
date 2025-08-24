# Pacotes
import os
import scanpy as sc
import pandas as pd
import seaborn as sns
import scvi
import numpy as np
import gc
import dtale

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
