import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import warnings
import pandas as pd

scv.settings.verbosity = 3
cr.settings.verbosity = 2


adata = scv.read('~/mouse_479QC.loom', cache=True)
pmeta = pd.read_csv('~/pmeta.txt' , delimiter = "\t")
#genes_used=pd.read_csv('~/genes_used.csv')


adata.var_names_make_unique()
use_cells=pmeta.iloc[:,0].tolist()
use_genes=genes_used.iloc[:,1].tolist()


adata = adata[use_cells]

adata.obs.index.tolist() == use_cells
scv.pl.proportions(adata)


adata=adata[:, use_genes]

scv.pp.filter_and_normalize(adata, min_shared_counts=0, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

scv.tl.recover_dynamics(adata, n_jobs=8)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.tl.umap(adata)

#adata.obs['Cell_type'] = pmeta['savercat'].values
#adata.obsm['X_umap'] = pmeta[['umap1', 'umap2']].values
adata.obs['Cell_type'] = pmeta['celltype'].values

## celltype modification of mouse brain data
Cell_type_sub=[]
for ii in pmeta['savercat'].values:
  if ii in ['Neuron', 'Neuroblast', 'Glioblast','Radial glia']:
  Cell_type_sub.append(ii)
else:
  Cell_type_sub.append('NA')

adata.obs['Cell_type']=Cell_type_sub


scv.pl.velocity_embedding_stream(adata,save= "umap_emb_stream.png", basis='umap',legend_loc='best',size=35, legend_fontsize=12, title='', smooth=.8, min_mass=4,  color='Cell_type')



# louvain clustering 
## find terminal states
#cr.tl.terminal_states(adata, cluster_key='clusters', weight_connectivities=0.2)
cr.tl.terminal_states(adata, cluster_key='Cell_type', weight_connectivities=0.2)
cr.pl.terminal_states(adata, save='umap_terminal.pdf')

## find initial states
#cr.tl.initial_states(adata, cluster_key='clusters')
cr.tl.initial_states(adata, cluster_key='Cell_type')

cr.pl.initial_states(adata, discrete=True,save='umap_initial.pdf')

# compute fate map
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False, save='umap_lineages_seperate.pdf')

cr.pl.lineages(adata, same_plot=True, save='umap_lineages_cbn.pdf')

#### precompute umap
adata.obsm['X_umap'] = pmeta[['umap1', 'umap2']].values
scv.pl.velocity_embedding_stream(adata,save= "umap_emb_stream_precompute.png", basis='umap',legend_loc='best',size=35, legend_fontsize=12, title='', smooth=.8, min_mass=4,  color='Cell_type')


#cellrank.tl.compute_terminal_states(adata)
# louvain clustering at the bottom
## find terminal states
#cr.tl.terminal_states(adata, cluster_key='clusters', weight_connectivities=0.2)
cr.tl.terminal_states(adata, cluster_key='Cell_type', weight_connectivities=0.2)
cr.pl.terminal_states(adata, save='umap_terminal_precompute.pdf')

## find initial states
#cr.tl.initial_states(adata, cluster_key='clusters')
cr.tl.initial_states(adata, cluster_key='Cell_type')

cr.pl.initial_states(adata, discrete=True,save='umap_initial_precompute.pdf')

# compute fate map
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False, save='umap_lineages_seperate_precompute.pdf')

cr.pl.lineages(adata, same_plot=True, save='umap_lineages_cbn_precompute.pdf')
