
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad
import matplotlib
from scipy.io import mmwrite
from scipy.io import mmread

scv.set_figure_params()
adata = scv.read('~/mouse_479QC.loom', cache=True)

# generate metadta for single-cells 

pmeta = pd.read_csv('~/pmeta.txt' , delimiter = "\t")
use_cells=pmeta.iloc[:,0].tolist()

adata = adata[use_cells]

adata.obs.index.tolist() == use_cells

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
mode='dynamical'
scv.tl.recover_dynamics(adata)
scv.tl.velocity_graph(adata)

#adata.obsm['X_umap'] = pmeta[['umap1', 'umap2']].values
adata.obs['Cell_type'] = pmeta['celltype'].values



scv.tl.umap(adata)

## new umap
scv.pl.umap(adata, save= "umap.pdf", color='Cell_type',legend_fontsize='xx-small', color_map = 'Paired')
scv.pl.velocity_embedding(adata, basis='umap', arrow_length=2, arrow_size=1.5, dpi=300,color='Cell_type',legend_fontsize='xx-small',color_map = 'Paired', save="embd.pdf")
scv.pl.velocity_embedding_grid(adata, basis='umap', alpha=0.3, arrow_length=1, arrow_size=3, arrow_color="black", dpi=300,color='Cell_type',color_map = 'Paired',legend_fontsize='xx-small', save="embd_grid.pdf", figsize = (6,6))
scv.pl.velocity_embedding_stream(adata, basis='umap',color='Cell_type',legend_fontsize='xx-small',color_map = 'Paired',save="embd_stream.png")

#scv.pl.velocity(adata_oidx, basis='umap', var_names=['Runx1'], save=test_text+"markerg_velocity_newumap.pdf")
#scv.pl.scatter(adata_oidx, basis=list_of_genes, save=test_text + "markerg_scatter.pdf")
scv.pl.velocity_graph(adata, basis='umap',color='Cell_type',legend_fontsize='xx-small',color_map = 'Paired', save="velocity_graph.pdf")



## precomputed umap
scv.pl.umap(adata, save= "umap_precomput_atac.pdf", color='Cell_type',legend_fontsize='xx-small', color_map = 'Paired')
scv.pl.velocity_embedding(adata, basis='umap', arrow_length=2, arrow_size=1.5, dpi=300,color='Cell_type',legend_fontsize='xx-small',color_map = 'Paired', save="embd_precomput_atac.pdf")
scv.pl.velocity_embedding_grid(adata, basis='umap', alpha=0.3, arrow_length=1, arrow_size=3, arrow_color="black", dpi=300,color='Cell_type',color_map = 'Paired',legend_fontsize='xx-small', save="embd_grid_precomput_atac.pdf", figsize = (6,6))
scv.pl.velocity_embedding_stream(adata, basis='umap',color='Cell_type',legend_fontsize='xx-small',color_map = 'Paired',save="embd_stream_precomput_atac.png")

#scv.pl.velocity(adata_oidx, basis='umap', var_names=['Runx1'], save=test_text+"markerg_velocity_newumap.pdf")
#scv.pl.scatter(adata_oidx, basis=list_of_genes, save=test_text + "markerg_scatter.pdf")
scv.pl.velocity_graph(adata, basis='umap',color='Cell_type',legend_fontsize='xx-small',color_map = 'Paired', save="velocity_graph_precomput_atac.pdf")
