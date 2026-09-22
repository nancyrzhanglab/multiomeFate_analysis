"""Run CoSPAR on a flat export written by `cospar_flat_io.R::export_for_cospar()`.

Usage (inside the `cospar` conda environment):

    python run_cospar.py EXPORT_DIR --t1 day10 --t2 week5 --fate-a High --fate-b Low

Reads EXPORT_DIR/{counts.mtx, cells.txt, genes.txt, obs.csv, X_pca.csv?} and
writes EXPORT_DIR/cospar_out/{obs.csv, dge_A.csv, dge_B.csv, run.json}.

Pipeline: initialize_adata_object -> (HVG + PCA when no X_pca.csv was exported)
-> infer_Tmap_from_multitime_clones(t1 -> t2) -> fate_bias(A vs B) ->
progenitor(A vs B) -> differential_genes(progenitor A vs progenitor B).
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as ssp
import scanpy as sc
import cospar as cs


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("export_dir")
    p.add_argument("--t1", required=True, help="earlier time_info label")
    p.add_argument("--t2", required=True, help="later time_info label")
    p.add_argument("--fate-a", required=True, help="state_info label of fate A at t2")
    p.add_argument("--fate-b", required=True, help="state_info label of fate B at t2")
    p.add_argument("--data-des", default=None,
                   help="cache key under cs.settings.data_path; defaults to the "
                        "export directory's basename. Must differ between datasets.")
    p.add_argument("--smooth", type=int, nargs="+", default=[20, 15, 10, 5])
    p.add_argument("--knn", type=int, default=20)
    p.add_argument("--sparsity", type=float, default=0.1)
    p.add_argument("--intraclone", type=float, default=0.2)
    p.add_argument("--max-iter", type=int, default=10)
    p.add_argument("--epsilon", type=float, default=0.01)
    p.add_argument("--bias-a", type=float, default=0.6,
                   help="fate_bias above this = progenitor of A")
    p.add_argument("--bias-b", type=float, default=0.4,
                   help="fate_bias below this = progenitor of B")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--n-pca", type=int, default=30,
                   help="PCs to compute when X_pca.csv is absent")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main():
    args = parse_args()
    d = args.export_dir
    np.random.seed(args.seed)

    cells = [l.rstrip("\n") for l in open(os.path.join(d, "cells.txt"))]
    genes = [l.rstrip("\n") for l in open(os.path.join(d, "genes.txt"))]
    X = ssp.csr_matrix(scipy.io.mmread(os.path.join(d, "counts.mtx")))
    obs = pd.read_csv(os.path.join(d, "obs.csv"), dtype=str, keep_default_na=False)
    obs.index = obs["cell_id"].values
    assert X.shape == (len(cells), len(genes)), "counts.mtx does not match cells/genes"
    assert list(obs.index) == cells, "obs.csv rows are not in cells.txt order"
    assert obs["time_info"].isin([args.t1, args.t2]).all(), "unexpected time_info label"
    assert obs.loc[obs["time_info"] == args.t2, "state_info"].isin(
        [args.fate_a, args.fate_b]).any(), "no t2 cell carries fate A or B"

    adata = sc.AnnData(X=X, obs=obs[["time_info", "clone_id", "state_info"]].copy())
    adata.var_names = genes

    pca_path = os.path.join(d, "X_pca.csv")
    X_pca = None
    if os.path.exists(pca_path):
        pca_df = pd.read_csv(pca_path, index_col=0)
        assert list(pca_df.index) == cells, "X_pca.csv rows are not in cells.txt order"
        X_pca = pca_df.to_numpy()

    # Long-format clone table -> X_clone via cospar's own helper, which aligns
    # on cell ID and refuses an all-zero result (the usual symptom of an ID mismatch).
    has_clone = obs["clone_id"].ne("NA") & obs["clone_id"].ne("")
    data_des = args.data_des or os.path.basename(os.path.normpath(d))
    cs.settings.data_path = os.path.join(d, "cospar_cache")
    cs.settings.figure_path = os.path.join(d, "cospar_cache", "figure")
    cs.settings.verbosity = 2

    adata = cs.pp.initialize_adata_object(
        adata=adata,
        cell_names=cells,
        time_info=obs["time_info"].values,
        state_info=obs["state_info"].values,
        X_pca=X_pca,
        X_emb=None if X_pca is None else X_pca[:, :2],
        data_des=data_des,
    )
    cs.pp.get_X_clone(adata,
                      clone_data_cell_id=list(obs.index[has_clone]),
                      clone_data_barcode_id=list(obs.loc[has_clone, "clone_id"]))

    if X_pca is None:
        cs.pp.get_highly_variable_genes(adata)
        cs.pp.get_X_pca(adata, n_pca_comp=args.n_pca)
        cs.pp.get_X_emb(adata)

    # time_ordering is lexicographic by default; force the intended order.
    cs.hf.update_time_ordering(adata, updated_ordering=[args.t1, args.t2])
    cs.hf.check_available_choices(adata)

    adata_t = cs.tmap.infer_Tmap_from_multitime_clones(
        adata,
        clonal_time_points=[args.t1],
        later_time_point=args.t2,
        smooth_array=args.smooth,
        CoSpar_KNN=args.knn,
        sparsity_threshold=args.sparsity,
        intraclone_threshold=args.intraclone,
        max_iter_N=args.max_iter,
        epsilon_converge=args.epsilon,
        compute_new=True,
    )

    fates = [args.fate_a, args.fate_b]
    src = "intraclone_transition_map"
    cs.tl.fate_bias(adata_t, selected_fates=fates, source=src, pseudo_count=0)
    cs.tl.progenitor(adata_t, selected_fates=fates, source=src,
                     bias_threshold_A=args.bias_a, bias_threshold_B=args.bias_b,
                     sum_fate_prob_thresh=0, avoid_target_states=True)

    bias_col = f"fate_bias_{src}_{args.fate_a}*{args.fate_b}"
    prog_a = f"progenitor_{src}_{args.fate_a}"
    prog_b = f"progenitor_{src}_{args.fate_b}"
    group_a = adata_t.obs[prog_a].astype(bool).values
    group_b = adata_t.obs[prog_b].astype(bool).values
    n_a, n_b = int(group_a.sum()), int(group_b.sum())
    print(f"progenitors: A={n_a}, B={n_b}")

    out = os.path.join(d, "cospar_out")
    os.makedirs(out, exist_ok=True)
    if n_a > 0 and n_b > 0:
        dge_a, dge_b = cs.tl.differential_genes(
            adata_t, cell_group_A=group_a, cell_group_B=group_b,
            FDR_cutoff=args.fdr, sort_by="ratio")
    else:
        dge_a = pd.DataFrame(columns=["gene", "ratio", "Qvalue"])
        dge_b = pd.DataFrame(columns=["gene", "ratio", "Qvalue"])
        print("WARNING: an empty progenitor group; no DE run", file=sys.stderr)
    dge_a.to_csv(os.path.join(out, "dge_A.csv"), index=False)
    dge_b.to_csv(os.path.join(out, "dge_B.csv"), index=False)
    adata_t.obs.to_csv(os.path.join(out, "obs.csv"))

    run = dict(vars(args))
    run.update({
        "cospar_version": cs.__version__,
        "scanpy_version": sc.__version__,
        "python_version": sys.version.split()[0],
        "n_cells_in_tmap_t1": int(len(adata_t.uns["Tmap_cell_id_t1"])),
        "n_cells_in_tmap_t2": int(len(adata_t.uns["Tmap_cell_id_t2"])),
        "n_progenitor_a": n_a, "n_progenitor_b": n_b,
        "n_dge_a": int(len(dge_a)), "n_dge_b": int(len(dge_b)),
        "fate_bias_column": bias_col,
        "progenitor_a_column": prog_a,
        "progenitor_b_column": prog_b,
        "pca_source": "exported X_pca.csv" if X_pca is not None else f"cospar HVG+PCA({args.n_pca})",
    })
    with open(os.path.join(out, "run.json"), "w") as f:
        json.dump(run, f, indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
