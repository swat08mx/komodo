#!pip install scanpy
#!pip install scvi-tools
#!pip install scikit-misc
#!pip install leidenalg
#!pip install gseapy
import gseapy as gp
import locale
from tqdm import tqdm
locale.getpreferredencoding = lambda: "UTF-8"
import scanpy as sc
import scvi
import seaborn as sns
import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
def pipeline(csv_path):
    adata = sc.read_csv(csv_path).T
    print("Read")
    sc.pp.filter_genes(adata, min_cells = 10)
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    df['prediction'] = solo.predict(soft = False)
    df['dif'] = df.doublet - df.singlet
    doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]

    adata = sc.read_csv(csv_path).T
    adata.obs['Sample'] = csv_path.split('_')[1]
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]


    sc.pp.filter_cells(adata, min_genes=200) #get rid of cells with fewer than 200 genes
    #sc.pp.filter_genes(adata, min_cells=3) #get rid of genes that are found in fewer than 3 cells
    ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 2]
    print("Function done")
    return adata

out=[]
lists=[]
path = "/home/ubuntu/uploads/"
pattern = r".*\.csv$"
#for file in tqdm(os.listdir(path)):
#  if re.search(pattern, file):
#    print("Matched")
#    lists.append(file)
#for names in lists:
#  out.append(pipeline('/home/ubuntu/' + names))
#data = sc.concat(out)

data = pipeline("/home/ubuntu/uploads/GSM5226577_C54ctr_raw_counts.csv")
print("File read")
sc.pp.filter_genes(data, min_cells=10)
data.X = csr_matrix(data.X)
data.obs.groupby('Sample').count()
data.layers['counts'] = data.X.copy()
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
data.raw = data
sc.pp.highly_variable_genes(data, n_top_genes=2000, subset = True, layer = 'counts', flavor= 'seurat_v3', batch_key="Sample")
scvi.model.SCVI.setup_anndata(data, layer='counts', categorical_covariate_keys=['Sample'], continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])
model = scvi.model.SCVI(data)
model.train()
data.obsm['X_scVI'] = model.get_latent_representation()
data.layers['scvi_normalized'] = model.get_normalized_expression(library_size=1e4)
sc.pp.neighbors(data, use_rep = 'X_scVI')
sc.tl.umap(data)
sc.tl.leiden(data, resolution=0.8)
sc.pl.umap(data, color = ['leiden', 'Sample'], frameon=False)
data.write_h5ad('integrated.h5ad')
sc.tl.rank_genes_groups(data, 'leiden')
sc.pl.rank_genes_groups(data, n_genes=20, sharey=False)
markers = sc.get.rank_genes_groups_df(data, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers_scvi = model.differential_expression(groupby='leiden')
markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > .5)]
sc.pl.umap(data, color = ['leiden'], frameon = False, legend_loc = "on data")
sc.pl.umap(data, color=['ABCA3', 'RTKN2', 'PECAM1'], frameon=False, layer = 'scvi_normalized')
dicto = {
    "0":"CD45",
    "1":"Macrophages",
    "2":"Pulmonary alveolar type I cells",
    "3":"Ionocytes",
    "4":"Pulmonary alveolar type II cells",
    "5":"Macrophages",
    "6":"NIL",
    "7":"NIL",
    "8":"NIL",
    "9":"NIL",
    "10":"Airway goblet cells",
    "11":"NIL",
    "12":"NIL",
    "13":"Ciliated cells",
    "14":"NIL"
}
data.obs['cell_type'] = data.obs.leiden.map(dicto)
sc.pl.umap(data, color = ['cell_type'], frameon=False)
data.uns['scvi_markers'] = markers_scvi
data.uns['markers'] = markers
model.save('model.model')
def map_contd(x):
  if 'cov' in x:
    return "COVID19"
  else:
    return 'control'
data.obs['condition'] = data.obs.Sample.map(map_contd)
num_tot_cells = data.obs.groupby(['Sample']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
cell_type_counts = data.obs.groupby(['Sample', 'condition', 'cell_type']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis=1) > 0].reset_index()
cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
cell_type_counts['total_cells'] = cell_type_counts.Sample.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.doublet / cell_type_counts.total_cells
plt.figure(figsize=(10, 4))
ax = sns.boxplot(data = cell_type_counts, x = 'cell_type', y = 'frequency', hue = 'condition')
plt.xticks(rotation = 35, rotation_mode= 'anchor', ha='right')
plt.show()
scvi_de = model.differential_expression(
    idx1 = [data.obs['cell_type'] == 'Macrophages'],
    idx2 = [data.obs['cell_type'] == 'CD45']
)
scvi_de = scvi_de[(scvi_de['is_de_fdr_0.05']) & (abs(scvi_de.lfc_mean) > .5)]
scvi_de = scvi_de.sort_values('lfc_mean')
scvi_de = scvi_de[(scvi_de.raw_normalized_mean1 > .5) | (scvi_de.raw_normalized_mean2) > .5]
genes_to_show = scvi_de[-15:].index.tolist() + scvi_de[:15].index.tolist()
subset = data[data.obs['cell_type'].isin(['CD45', 'Macrophages'])].copy()
sc.pl.heatmap(subset, genes_to_show, groupby='cell_type', swap_axes=True, layer='scvi_normalized', log=True)
enr = gp.enrichr(gene_list = scvi_de.index.tolist(),
                 gene_sets = ['KEGG_2021_Human'],
                 organism = 'human',
                 outdir=None,
                 background = subset.var_names.to_list())
print("DONE")
