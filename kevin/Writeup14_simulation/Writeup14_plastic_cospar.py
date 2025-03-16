# cospar virtual environment
import cospar as cs
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np

cs_input = sc.read_h5ad("/home/stat/nzh/team/kevinl1/nzhanglab/project/Multiome_fate/out/kevin/Writeup14/Writeup14_plastic-setting_seurat_CoSpar_prepared.h5ad")
pca_dat = pd.read_csv("/home/stat/nzh/team/kevinl1/nzhanglab/project/Multiome_fate/out/kevin/Writeup14/plastic-setting_embeddings.csv",index_col=0)

cs.logging.print_version()
cs.settings.verbosity = 2
cs.settings.data_path = "simulation_data_plastic"  # A relative path to save data. If not existed before, create a new one.
cs.settings.figure_path = "simulation_data_plastic_figure"  # A relative path to save figures. If not existed before, create a new one.

cs_input
cs_input.obs["ident"]

pca_dat

scale_dat = pd.read_csv("/home/stat/nzh/team/kevinl1/nzhanglab/project/Multiome_fate/out/kevin/Writeup14/plastic-setting_scale_data.csv",index_col=0)
scale_dat.shape
scale_dat
scale_dat.index
scale_dat.T.index

adata = sc.AnnData(scale_dat.T,obs = cs_input.obs)
adata.obs["cellID"] = adata.obs.index
adata.obsm['X_pca'] = pca_dat
adata.obsm['X_emb'] = cs_input.obsm['FT.COCL2.UMAP']
adata
adata.obsm['X_emb']

clone_data = pd.get_dummies(data =adata.obs[["cellID","assigned_lineage"]], dtype=float, columns=["assigned_lineage"])
clone_data =clone_data.set_index('cellID')
clone_data

meta = pd.read_csv("/home/stat/nzh/team/kevinl1/nzhanglab/project/Multiome_fate/out/kevin/Writeup14/plastic-setting_meta.csv",index_col=0)
meta

adata.obs = meta
adata
adata.obsm['X_emb'] = adata.obsm['X_emb'].to_numpy()

# run cospar
adata_original = cs.pp.initialize_adata_object(
    adata=adata,
    cell_names=scale_dat.T.index,
    time_info=adata.obs["dataset"],
    X_clone=clone_data,
    state_info=adata.obs["factor_lineage"],
    data_des="plastic_simulation",
    X_emb = adata.obsm['X_emb'],
    X_pca = adata.obsm['X_pca']
)
adata_original
cs.hf.check_available_choices(adata_original)

adata = cs.tmap.infer_Tmap_from_multitime_clones(
    adata_original,
    clonal_time_points=['day10_COCL2'],
    later_time_point='week5_COCL2',
    smooth_array=[20,15,10,5],
    sparsity_threshold=0.1,
    intraclone_threshold=0.2,
    max_iter_N=10,
    epsilon_converge=0.01,
    compute_new=True
)

cs.tl.fate_coupling(adata_original,source='X_clone')
cs.tl.fate_hierarchy(adata_original,source='X_clone')
cs.tl.clonal_fate_bias(adata_original, ['High', 'Low'])

cs.tl.progenitor(
    adata,
    selected_fates=["High", "Low"],
    source="intraclone_transition_map",
    map_backward=True,
    bias_threshold_A=0.5,
    bias_threshold_B=0.5,
    sum_fate_prob_thresh=0.2,
    avoid_target_states=True,
)

cs.tl.fate_bias(
    adata,
    selected_fates=["High", "Low"],
    source="intraclone_transition_map",
    pseudo_count=0,
)

adata
adata.obs.to_csv('/home/stat/nzh/team/kevinl1/nzhanglab/project/Multiome_fate/out/kevin/Writeup14/simulation_plastic_cospar_obs.csv')