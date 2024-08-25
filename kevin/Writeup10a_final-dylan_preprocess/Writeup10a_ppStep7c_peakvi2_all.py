import datetime
import io
import random
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import session_info
import torch

def capture_output(func):
    """
    Capture the output of a function that prints to consule

    Parameters:
    - func: Function that is run without parameters

    Returns:
    - The output of func()
    """
    original_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    func()
    sys.stdout = original_stdout
    captured_output = new_stdout.getvalue()
    new_stdout.close()
    return captured_output


date_of_run = datetime.datetime.now()
sessionInfo = capture_output(session_info.show)
torch.set_float32_matmul_precision("high")
torch.manual_seed(0)
random.seed(0)
np.random.seed(0)
scvi.settings.seed = 0

# see https://docs.scvi-tools.org/en/stable/tutorials/notebooks/atac/PeakVI.html
adata = scvi.data.read_h5ad("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/Writeup10a_ppStep7c_peakvi-prep_all.h5ad")

print("# regions before filtering:", adata.shape[-1])

# compute the threshold: 5% of the cells
min_cells = int(adata.shape[0] * 0.05)
sc.pp.filter_genes(adata, min_cells=min_cells)

print("# regions after filtering:", adata.shape[-1])

scvi.model.PEAKVI.setup_anndata(adata)

print("Training PeakVI")
model = scvi.model.PEAKVI(adata)
model.train()

model.save("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/Writeup10a_ppStep7c_peakvi_all", overwrite=True)
# model = scvi.model.PEAKVI.load("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup6o/Writeup6o_all-data-atac_CIS_peakVI", adata=adata)

PEAKVI_LATENT_KEY = "X_peakvi"

latent = model.get_latent_representation()
adata.obsm[PEAKVI_LATENT_KEY] = latent

print("Saving peakVI as CSV")
cell_names = adata.obs_names
df = pd.DataFrame(latent, index = cell_names)
df.to_csv('/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/Writeup10a_ppStep7c_peakvi_all_embedding.csv')

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
plt.savefig("/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/Writeup10a_ppStep7c_peakvi_all_umap.png", dpi=300)

print("Done! :)")