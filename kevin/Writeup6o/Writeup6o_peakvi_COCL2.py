import os
import tempfile

import leidenalg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch

# see https://docs.scvi-tools.org/en/stable/tutorials/notebooks/atac/PeakVI.html

sc.set_figure_params(figsize=(4, 4), frameon=False)
torch.set_float32_matmul_precision("high")

adata = scvi.data.read_h5ad("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup6o/Writeup6o_all-data-atac_COCL2.h5ad")
adata

print("# regions before filtering:", adata.shape[-1])

# compute the threshold: 5% of the cells
min_cells = int(adata.shape[0] * 0.05)
# in-place filtering of regions
sc.pp.filter_genes(adata, min_cells=min_cells)

print("# regions after filtering:", adata.shape[-1])

scvi.model.PEAKVI.setup_anndata(adata)

print("Training PeakVI")
model = scvi.model.PEAKVI(adata)
model.train()

model.save("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup6o/Writeup6o_all-data-atac_COCL2_peakVI", overwrite=True)
# model = scvi.model.PEAKVI.load("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup6o/Writeup6o_all-data-atac_COCL2_peakVI", adata=adata)

PEAKVI_LATENT_KEY = "X_peakvi"

latent = model.get_latent_representation()
adata.obsm[PEAKVI_LATENT_KEY] = latent

print("Saving peakVI as CSV")
cell_names = adata.obs_names
df = pd.DataFrame(latent, index = cell_names)
df.to_csv('/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup6o/Writeup6o_all-data-atac_COCL2_peakVI.csv')

print("Plotting")
PEAKVI_CLUSTERS_KEY = "dataset"

# compute the k-nearest-neighbor graph that is used in both clustering and umap algorithms
sc.pp.neighbors(adata, use_rep=PEAKVI_LATENT_KEY)
# compute the umap
sc.tl.umap(adata, min_dist=0.2)

sc.pl.umap(adata, 
           color=PEAKVI_CLUSTERS_KEY,
           legend_loc='on data',
           legend_fontsize=12, 
           legend_fontoutline=2)
plt.savefig("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup6o/Writeup6o_peakVI_umap-COCL2_plot.png", dpi=300)

print("Done! :)")