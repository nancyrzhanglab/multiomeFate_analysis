import anndata as ad
import pandas as pd
import numpy as np

adata = ad.read("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry.h5ad")

# Check if the layers exist before trying to delete them
layers_to_remove = ['ambiguous', 'spliced', 'unspliced']
for layer in layers_to_remove:
    if layer in adata.layers:
        del adata.layers[layer]


adata.X = adata.layers['matrix'].copy()
adata.X = adata.X.astype(np.float32)

del adata.layers['matrix']

adata.write("/home/users/kzlin/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry_kzlin-cleaned.h5ad")
a