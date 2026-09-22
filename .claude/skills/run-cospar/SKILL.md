---
name: run-cospar
description: Use when benchmarking against CoSPAR (Wang et al. 2022) on lineage-barcoded single-cell data generated or held in R — exporting the dataset to flat files, running the transition-map inference and fate-bias / progenitor / differential-gene steps in the `cospar` conda environment, and reading the per-cell fate bias and implicated gene sets back into R. Also use when a cospar call aborts with "no multi-time clones", returns a constant 0.5 fate bias, or reuses a stale similarity-matrix cache.
---

# Run CoSPAR

CoSPAR (Coherent Sparse optimization) infers a cell-by-cell transition map between
two time points from a cell-state similarity graph plus clonal barcodes, then reads
off a per-cell **fate bias** toward one of two later-time-point fates. Genes are
implicated by differential expression between the biased progenitor groups, which is
the paper's own recipe (Wilcoxon, BH, FDR 0.05, ranked by fold change).

The hand-off between R and Python is a **directory of flat files**, the same
principle as `seurat-anndata-bridge`. This skill ships its own two templates because
the input here is a plain count matrix plus three per-cell columns, not a Seurat
object; reach for `seurat-anndata-bridge` only when a real Seurat object with
reductions and metadata has to cross over.

## Environment

The `cospar` conda environment on Kevin's laptop is the `COSPAR_ENV` location
(resolve the path through `CLAUDE_[name].md`); it holds cospar 0.4.1 on Python 3.9.
The source checkout is `COSPAR_SRC` (0.5.0). The API this skill uses is identical in
both. Invoke Python by absolute path rather than activating the environment:

```bash
"$COSPAR_ENV/bin/python" -c "import cospar as cs; print(cs.__version__)"
```

A fresh environment, when needed, follows the package's own recipe:
`conda create -n cospar python=3.9 && pip install cospar`. Done when the version
prints.

## Steps

1. **Copy both templates beside the analysis script**:
   `templates/cospar_flat_io.R` and `templates/run_cospar.py`. Copy, do not
   `source()` them out of the skill folder, so the result runs without the skill.

2. **Shape the three per-cell vectors in R.** Every cell at both time points gets:
   - `time_info`: exactly two labels, for example `"day10"` and `"week5"`.
   - `clone_id`: the barcode, `NA` when unbarcoded. **A clone constrains the map
     only if it has cells at both time points**; the exporter counts these and stops
     at zero.
   - `state_info`: the fate vocabulary. The convention from the paper's simulations
     is `"t1"` for every earlier cell and, for later cells, `"High"` when the clone
     is among the top-expanding ones and `"Low"` otherwise. Any two later-time-point
     labels work; they are what `--fate-a` / `--fate-b` name.

   Counts are non-negative and **not log-transformed**; cospar's HVG and PCA steps
   assume raw counts. Passing `pca_mat` makes cospar build its graph on your
   embedding and skip those steps, which is the right call when CYFER used the same
   embedding and the comparison should hold the representation fixed.

3. **Export** with `export_for_cospar(count_mat, time_info, clone_id, state_info,
   export_dir, pca_mat = NULL)`. Put `export_dir` under an `out/` tree, never in the
   code tree: the export plus cospar's similarity cache is tens of MB. Done when
   `manifest.json` reports the expected cell count and a positive
   `num_clones_multitime`.

4. **Run** from the shell:

   ```bash
   "$COSPAR_ENV/bin/python" run_cospar.py EXPORT_DIR --t1 day10 --t2 week5 \
       --fate-a High --fate-b Low --data-des NAME
   ```

   `--data-des` keys the on-disk similarity cache inside `EXPORT_DIR/cospar_cache/`;
   give every dataset its own name or a replicate silently reuses the previous graph.
   The defaults reproduce the paper's earlier runs (`smooth 20 15 10 5`,
   `sparsity 0.1`, `intraclone 0.2`, `max-iter 10`, `epsilon 0.01`), plus progenitor
   thresholds `0.6 / 0.4`. Done when the log ends with `progenitors: A=…, B=…` with
   both counts positive and `wrote …/cospar_out`. A zero count means the bias
   collapsed to one side; loosen the thresholds toward 0.5 or check that both fates
   are present at `t2`.

5. **Import** with `import_cospar_results(export_dir)`. It returns `fate_bias` (named
   over the earlier cells, in `[0, 1]`, **above 0.5 = toward fate A**), the two
   progenitor masks, `dge_a_df` / `dge_b_df`, and `run` with the versions and
   parameters. The implicated gene set for "associated with expansion" is
   `dge_a_df$gene` when fate A is the high-expansion label. Done when
   `length(fate_bias)` equals the number of earlier-time-point cells.

## Reading the outputs

- Cells outside the transition map are **0.5, not `NA`**, in the fate-bias column.
  The importer already restricts to `t1`; keep that restriction in any correlation.
- The Python column is `fate_bias_intraclone_transition_map_A*B` with a literal
  `*`. The importer reads with `check.names = FALSE` so the name survives; a plain
  `read.csv()` turns the `*` into `.`.
- `tl.progenitor` with `avoid_target_states = TRUE` is what prevents later-time
  cells from entering their own progenitor group.
- A smoke test lives in this skill's history: 1200 cells, 30 clones, 200 genes with
  10 causal genes recovered at Jaccard 1.0 in about 15 seconds. Runtime grows with
  the kNN graph, so a 10k-cell export takes minutes, not seconds.

## Prior usage

`kevin/Writeup14_simulation/*cospar*` and `kevin/Writeup17b_simulation-plots/
Writeup17b_barplot-*.R` are the earlier, notebook-style runs on the semi-synthetic
data. They fed cospar scaled data as `adata.X`, built the clone matrix with
`pd.get_dummies`, did the DE in R against fastTopics-denoised expression, and used
two contradictory winner/loser conventions (`< 0.5` in Writeup14, `> 0.6 / < 0.4`
in Writeup17b). This skill replaces all four choices; read those scripts for the
figure code, not the pipeline.
