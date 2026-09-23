"""Run CoSPAR on a flat export written by `cospar_flat_io.R::export_for_cospar()`.

Usage (inside the `cospar` conda environment):

    python run_cospar.py EXPORT_DIR --t1 day10 --t2 week5 --fate-a High --fate-b Low

Reads EXPORT_DIR/{counts.mtx, cells.txt, genes.txt, obs.csv, X_pca.csv?} and
writes EXPORT_DIR/cospar_out/{obs.csv, dge_A.csv, dge_B.csv, run.json}.

Pipeline: initialize_adata_object -> (HVG + PCA when no X_pca.csv was exported)
-> infer_Tmap_from_multitime_clones(t1 -> t2) -> fate_bias(A vs B) ->
progenitor(A vs B) -> differential_genes(progenitor A vs progenitor B).

What CoSPAR computes, in one paragraph, so the calls below can be read against
it. The transition map T between the t1 and t2 cells is built from two
ingredients: the *barcode matrix* M (a t1-cell x t2-cell 0/1 matrix that ties
every t1 cell of a clone to every t2 cell of the same clone, with equal
weight - the barcode says which clone a cell came from, not which parent) and
two *within-time-point similarity smoothers* S_t1 and S_t2 (kNN graphs on the
expression embedding, one per time point). CoSPAR alternates between smoothing
T = S_t1 . M . S_t2 and re-sparsifying, so that cells that look alike at t1
get similar rows even when their barcodes differ; the t1-to-t2 expression
similarity is never used. The "intraclone" variant then restricts the smoothed
map back to barcode-consistent entries. The per-cell *fate bias* of a t1 cell
toward fate A is the share of its transition-map mass that lands on t2 cells
labeled A, normalized against the overall A/B mass ratio, so 0.5 means "no
preference". Cells outside the map (t1 cells whose clone has no t2 cells, when
`extend_Tmap_space` is left False) are *filled in* with 0.5, not NA - callers
that care must exclude such cells themselves before the export.
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
    # Defaults reproduce the CoSPAR paper's own simulation settings
    # (smooth_array 20 15 10 5, sparsity 0.1, intraclone 0.2), which are also
    # what the multiomeFate paper's earlier CoSPAR runs used.
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
    # smooth: the rounds of similarity smoothing. Each entry is the depth of
    # the kNN diffusion used in that round (larger = smoother); CoSPAR runs
    # one iteration per entry, so [20, 15, 10, 5] means four smooth-then-
    # sparsify rounds with progressively more local smoothing.
    p.add_argument("--smooth", type=int, nargs="+", default=[20, 15, 10, 5])
    # knn: neighbors per cell in the within-time-point similarity graphs
    # (S_t1, S_t2). This is the knob that decides how far a clone's fate
    # information spreads to transcriptomic neighbors.
    p.add_argument("--knn", type=int, default=20)
    # sparsity: after each smoothing round, transition-map entries below this
    # quantile-based threshold are zeroed, keeping the map sparse.
    p.add_argument("--sparsity", type=float, default=0.1)
    # intraclone: the threshold used when restricting the smoothed map back
    # to barcode-consistent (same-clone) entries; the fate bias below is read
    # from this intraclone map, not from the fully smoothed one.
    p.add_argument("--intraclone", type=float, default=0.2)
    # max-iter / epsilon: stop the smooth-sparsify alternation after this
    # many iterations or when the map changes by less than epsilon.
    p.add_argument("--max-iter", type=int, default=10)
    p.add_argument("--epsilon", type=float, default=0.01)
    # bias-a / bias-b: fate_bias is in [0, 1] with 0.5 = no preference;
    # cells above 0.6 are called progenitors of A, below 0.4 progenitors of
    # B, and the band between is left uncalled. Only the progenitor/DE step
    # uses these; the fate bias itself is exported unthresholded.
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

    # --- Read the flat export -------------------------------------------------
    # cells.txt / genes.txt fix the row and column order of counts.mtx;
    # obs.csv carries the three per-cell columns (time_info, clone_id,
    # state_info) in the same row order. Everything is cross-checked before
    # any CoSPAR call so an ID mismatch fails here, loudly, rather than as a
    # silently wrong clone matrix.
    cells = [l.rstrip("\n") for l in open(os.path.join(d, "cells.txt"))]
    genes = [l.rstrip("\n") for l in open(os.path.join(d, "genes.txt"))]
    # counts.mtx holds raw (not log-transformed) counts: CoSPAR's own HVG and
    # PCA steps, used only when no X_pca.csv was exported, assume raw counts.
    X = ssp.csr_matrix(scipy.io.mmread(os.path.join(d, "counts.mtx")))
    # keep_default_na=False: a clone literally named "NA" must stay a string;
    # unbarcoded cells are the literal string "NA" written by the R exporter.
    obs = pd.read_csv(os.path.join(d, "obs.csv"), dtype=str, keep_default_na=False)
    obs.index = obs["cell_id"].values
    assert X.shape == (len(cells), len(genes)), "counts.mtx does not match cells/genes"
    assert list(obs.index) == cells, "obs.csv rows are not in cells.txt order"
    assert obs["time_info"].isin([args.t1, args.t2]).all(), "unexpected time_info label"
    assert obs.loc[obs["time_info"] == args.t2, "state_info"].isin(
        [args.fate_a, args.fate_b]).any(), "no t2 cell carries fate A or B"

    adata = sc.AnnData(X=X, obs=obs[["time_info", "clone_id", "state_info"]].copy())
    adata.var_names = genes

    # --- Optional shared embedding -------------------------------------------
    # When the R side exported X_pca.csv (the embedding CYFER and the AED
    # also use), CoSPAR builds its similarity graphs on it and its own
    # HVG/PCA below is skipped entirely: the comparison holds the
    # representation fixed, so differences between methods are differences
    # in the model, not in preprocessing.
    pca_path = os.path.join(d, "X_pca.csv")
    X_pca = None
    if os.path.exists(pca_path):
        pca_df = pd.read_csv(pca_path, index_col=0)
        assert list(pca_df.index) == cells, "X_pca.csv rows are not in cells.txt order"
        X_pca = pca_df.to_numpy()

    # Long-format clone table -> X_clone via cospar's own helper, which aligns
    # on cell ID and refuses an all-zero result (the usual symptom of an ID mismatch).
    has_clone = obs["clone_id"].ne("NA") & obs["clone_id"].ne("")
    # data_des keys the on-disk similarity cache below cospar_cache/. Two
    # datasets sharing a key would silently reuse each other's similarity
    # matrices, so the caller must pass a unique --data-des per dataset.
    data_des = args.data_des or os.path.basename(os.path.normpath(d))
    cs.settings.data_path = os.path.join(d, "cospar_cache")
    cs.settings.figure_path = os.path.join(d, "cospar_cache", "figure")
    cs.settings.verbosity = 2

    # initialize_adata_object: stamps the AnnData with the fields every later
    # CoSPAR call reads (obs["time_info"], obs["state_info"], obsm["X_pca"],
    # obsm["X_emb"] for plots, uns["data_des"]). X_emb is only a 2-D plotting
    # convenience; the first two PCs are fine because nothing numerical reads it.
    adata = cs.pp.initialize_adata_object(
        adata=adata,
        cell_names=cells,
        time_info=obs["time_info"].values,
        state_info=obs["state_info"].values,
        X_pca=X_pca,
        X_emb=None if X_pca is None else X_pca[:, :2],
        data_des=data_des,
    )
    # get_X_clone: builds obsm["X_clone"], the cell x clone 0/1 membership
    # matrix, from the (cell_id, clone_id) pairs of the barcoded cells. This
    # matrix is the only t1-to-t2 link CoSPAR ever uses; a clone must appear
    # at both time points to contribute a constraint.
    cs.pp.get_X_clone(adata,
                      clone_data_cell_id=list(obs.index[has_clone]),
                      clone_data_barcode_id=list(obs.loc[has_clone, "clone_id"]))

    # Fallback preprocessing, used only when the caller exported no embedding:
    # CoSPAR's own highly-variable-gene selection, then PCA on those genes.
    if X_pca is None:
        cs.pp.get_highly_variable_genes(adata)
        cs.pp.get_X_pca(adata, n_pca_comp=args.n_pca)
        cs.pp.get_X_emb(adata)

    # time_ordering is lexicographic by default; force the intended order.
    cs.hf.update_time_ordering(adata, updated_ordering=[args.t1, args.t2])
    cs.hf.check_available_choices(adata)

    # --- Transition map -------------------------------------------------------
    # infer_Tmap_from_multitime_clones is the barcoded ("multi-time clones")
    # mode: the map is the alternation T = S_t1 . M . S_t2 described in the
    # module docstring, where M comes from X_clone and S_t1/S_t2 are kNN
    # similarity graphs on the (shared) PCA, one per time point. Only t1
    # cells whose clone has t2 cells enter the map's row space: with the
    # default extend_Tmap_space=False (not overridden here), t1 cells of
    # extinct/single-time clones are outside the map and later get the
    # filled-in fate bias of 0.5. A caller that considers those cells
    # informative (e.g. their clone size 0 IS the signal) must exclude or
    # handle them on its own side; this script does not decide that.
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

    # --- Fate bias ------------------------------------------------------------
    # fate_bias(A, B) on the intraclone map: for each t1 cell, the map mass
    # landing on t2 cells labeled A over the mass landing on A or B,
    # normalized so 0.5 = the population-average split. pseudo_count=0 means
    # no shrinkage toward 0.5; a cell whose mass is entirely on A gets 1.
    # The source is the intraclone map (the smoothed map restricted back to
    # same-clone entries), the choice CoSPAR's tutorials use for fate bias.
    fates = [args.fate_a, args.fate_b]
    src = "intraclone_transition_map"
    cs.tl.fate_bias(adata_t, selected_fates=fates, source=src, pseudo_count=0)
    # progenitor: thresholds the bias into two progenitor groups (above
    # bias-a -> A, below bias-b -> B); sum_fate_prob_thresh=0 keeps cells
    # with little mass on either fate, and avoid_target_states=True keeps t2
    # cells of the two fates themselves out of the progenitor groups.
    cs.tl.progenitor(adata_t, selected_fates=fates, source=src,
                     bias_threshold_A=args.bias_a, bias_threshold_B=args.bias_b,
                     sum_fate_prob_thresh=0, avoid_target_states=True)

    # CoSPAR writes its results as obs columns whose names embed the source
    # map and the fate pair (with a literal `*`); run.json records the exact
    # spellings so the R importer never has to reconstruct them.
    bias_col = f"fate_bias_{src}_{args.fate_a}*{args.fate_b}"
    prog_a = f"progenitor_{src}_{args.fate_a}"
    prog_b = f"progenitor_{src}_{args.fate_b}"
    group_a = adata_t.obs[prog_a].astype(bool).values
    group_b = adata_t.obs[prog_b].astype(bool).values
    n_a, n_b = int(group_a.sum()), int(group_b.sum())
    print(f"progenitors: A={n_a}, B={n_b}")

    # --- Differential expression between the progenitor groups ---------------
    # CoSPAR's own gene recipe (Wilcoxon on the raw counts, BH, ranked by
    # fold change). dge_A holds genes up in the A-progenitors. A caller
    # using the "matched" route (correlating genes with the fate bias in R)
    # will ignore these files; they are written for the recipe-faithful route.
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
    # The full obs table is the main output: it carries the fate-bias and
    # progenitor columns for every cell (t1 cells outside the map at 0.5).
    adata_t.obs.to_csv(os.path.join(out, "obs.csv"))

    # run.json: every argument, the package versions, the map dimensions and
    # the exact obs column names, so a run is auditable and the importer is
    # decoupled from CoSPAR's naming scheme.
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
