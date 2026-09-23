# CLAUDE_kevin.md — Kevin's Context

> **Current state only.** Every narrative section below is updated *in place* each session — overwrite, don't append. The External Locations table is keyed by (location, machine): replace the matching row when a path changes, add a row only for a new pair, and delete rows that no longer exist. The append-only dated log lives in `HISTORY_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session writes this file and `HISTORY_kevin.md`. Other collaborators may read it for context but must not edit it.

## About Kevin
- Role in project: co-first author; owns the CYFER R package (`multiomeFate`, maintainer) and the `kevin/` analysis tree; runs the simulations
- Background: statistics faculty (UW Biostatistics); R package development, single-cell multiome methods
- Email: kzlin@uw.edu

## External Locations (per-machine paths)
Resolves the location names declared in the master `CLAUDE.md` → *External Locations* to real paths **on Kevin's machines**. One row per (location, machine) pair.

| Location name | Machine | Path | Notes |
|---|---|---|---|
| `CYFER_PKG` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/multiomeFate/` | branch `master`; source 1.0.3.000 and installed copy 1.0.3.0 match as of 2026-09-22; `devtools::install_local()` hangs here, use `R CMD INSTALL` from a scratch copy |
| `PAPER_REPO` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/multiome_fate_paper/` | also reachable as `/Users/kevinlin/Dropbox/...` (same folder); `paper_nbt.tex` is the newer source; no CLAUDE.md there |
| `PAPERS_DIR` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/papers/` | `multiome_fate_paper.pdf`, `simulation-planning.txt` |
| `OUT_ROOT` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/out/` | `Writeup_Simulations/sim1..7_*.rds` present; `Writeup14/` RData + CoSPAR obs CSVs present |
| `OUT_ROOT` | UW Biostat cluster (SLURM, account `biostat`) | `/home/users/kzlin/kzlinlab/projects/multiomeFate/out/kevin/` | what `sim3_heritability.slurm` and the sim scripts' `~/kzlinlab/...` paths resolve to |
| `OUT_ROOT` | Wharton HPC (legacy) | `/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/` | Writeup14 and earlier were run here; treat as archival |
| `COSPAR_SRC` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/cospar/` | version 0.5.0 checkout; tutorials in `docs/source/*.ipynb` |
| `COSPAR_ENV` | personal laptop (macOS) | `/opt/miniconda3/envs/cospar/` | cospar 0.4.1, Python 3.9.20, scanpy 1.10.3; smoke-tested 2026-09-22 |
| `COSPAR_ENV` | UW Biostat cluster | *(not present)* | no cospar environment built there yet |

## Project Status (as of 2026-09-23)
- Writeup21's trailblazing is complete through the `_r2` pass, which implements every decision from Kevin's round-1 review (`additional_context/responses-to-Writeup21-trailblazing_round1.txt`): 20-PC shared embedding (now the drivers' default), re-spaced Gini ladder {0.3–0.9}, extinct-clone t1 cells excluded before CoSPAR, split-half reliability removed, AED re-centring and seed scheme kept. Code in `kevin/Writeup21_new-simulations/` now also has `check_permutation_claude.R` (memo §7 label-permutation check) and `make_umaps_claude.R` (UMAPs of 6 example datasets); `run_cospar.py` (template + copy, identical) is heavily commented for Kevin's manual review.
- Output suffix scheme: `paste0("_r2", if(d_pca != 20) "_pca<d>", if(dry) "_dry")`; the round-1 files ("" = 10-PC, `_pca20`) stay on disk as the record. `_r2` RDS in `OUT_ROOT/Writeup21_new-simulations/`, CSVs in `csv/kevin/Writeup21/` (plus `check_permutation_r2.csv`, `umap_coords_r2.csv`), figures in `fig/kevin/Writeup21/`.
- `Writeup21_trailblazing-report_claude.{Rmd,html}` re-knitted 2026-09-23 for Kevin's review: `_r2` headline curves, UMAP section, explicit "lineage-DE has no per-cell fate score" note, fixed row-name leak in the realized tables, sanity-check section with the permutation table, problems-and-resolutions section.
- **Headline (`_r2`):** CYFER 0.90–0.93 across the Gini axis (0.52 at Gini 0.9, fate score still 0.82) and 0.79–0.97 on the AED axis. The CoSPAR exclusion helps it substantially on the Gini axis: near 0.70 through Gini 0.6 (was 0.28) and ~0.31 at 0.9 (was −0.13), erratic at 0.7–0.9 (SDs to 0.31); on the AED axis it moves ≤ +0.02 and CYFER/lineage-DE reproduce `_pca20` to 4 decimals, so CoSPAR's AED fall-off (0.76 → −0.03) is genuine.
- The "last easy level" now reads Gini 0.6 (one notch above the held 0.5; within two-replicate noise, so 0.5 stands unless the final round says otherwise) and AED 0.4 (confirmed).
- The memo `simulation_design_claude.md` is revised self-contained to the settled design and the `_r2` values; §9 is split into Settled (CoSPAR ≪ CYFER on 4B; uniform-link framing; noncausal) and Open (top-Gini regime, fixed Gini 0.5 vs 0.6, τ=0 top AED level, 300-weak-gene supplement). Memo revision rule: self-contained, between-round changes in `HISTORY_kevin.md`.
- Permutation check: metrics scatter around 0 with no consistent sign; CYFER's null spread is wide (−0.36 to 0.61 over 10 permutations) because under permuted labels its fit is an overfit random 20-PC direction (chance alignment SD ≈ 1/√20); a formal p-value would need hundreds of permutations.
- Installed `multiomeFate` is 1.0.3.0, matching the source, via `R CMD INSTALL` from a scratch copy (`devtools::install_local()` hangs on this laptop; Kevin approved this as the standing workaround).
- sim1–sim7 (`Writeup_Simulations/`) are complete with RDS outputs in `OUT_ROOT/Writeup_Simulations/` and CSVs in `csv/kevin/Writeup_Simulations/`.
- Uncommitted on branch `kevin`: everything from Writeup21 (scripts, CSVs, figures, memo and report revisions, the commented `run-cospar` template) and the two per-person files.

## Key Methodological Details
- CYFER's effective sample size is the number of clones; `cyfer()` refuses a fit when features ≥ training clones, so genes always go through an embedding (fastTopics on real data, PCA in synthetic sims). With 100 clones and 3 folds the cap is about 60 features, so 20 PCs is safe.
- `cell_imputed_score` is log10(expected progeny); `coefficient_vec` is natural-log; `lineage_imputed_count` is counts. `exp()` on the score is the documented mistake.
- The package accepts zero-count clones; every sim1–sim7 script filters them out. Writeup21 keeps them (`method_cyfer()`), and CYFER converged on every dataset with 0–33 extinct clones.
- Writeup21 generator: 100 clones × 10 t1 cells; 3,000 expected t2 cells via `β_0`, zeros kept, no cap; latent `d = 10`, `τ = s·sqrt(h²)`, `σ_w = s·sqrt(1 − h²)`, `g_l ~ Gamma(1.5, 1.5)` on the 9 non-causal coordinates (`spread_variation = "noncausal"`); `β = 1.5`; children `ρ = 0.8` with shared shift `||δ|| = 3s` and clone shift `τ_δ = 0.6`, both on non-causal coordinates; NB counts `θ = 10`, 2,000 genes, 100 causal. Seeds: replicate `r` of level `j` is `seed_base + 100j + r` (1000 for Gini, 2000 for AED); the generator's three stages use `10·seed + {0, 1, 2}`; calibration draws at `+50 + k`; permutation checks at `+80 + k`.
- Calibration: `h² = 1 − AED target`, then `s` bisected (geometric, 10 latent draws, tol 0.01) for the Gini target, then one full dataset's realized mean AED checked and `h²` shifted by the miss (up to 3 times, tol 0.03) with `s` re-bisected. `_r2` values: Gini axis (levels 0.3–0.9) `h² = 0.83, 0.675, 0.632, 0.6, 0.6, 0.6, 0.6` with `s` 0.37–2.06; AED axis `h² = 1, 0.832, 0.652, 0.45, 0.30, 0.106, 0` with `s` 0.62–1.25. Realized mean AED 0.40–0.43 across the Gini axis; realized Gini within ~0.06 of target in the mean of two replicates (per-replicate scatter ±0.05); floors: AED 0.1 realizes 0.12, AED 1.0 realizes 0.94.
- CoSPAR convention: t1 cells of extinct clones are excluded before the export (`method_cospar()`); their bias is `NA`, the gene correlation and fate-score diagnostic run over included cells, and `num_t1_excluded_cospar` is a replicate column. With `extend_Tmap_space = False` those cells would otherwise be filled in at bias 0.5.
- The causal latent axis is the weakest direction in gene space (100 genes × loading² ≈ 0.58 versus 2,000 × 0.0625 per non-causal coordinate), so it lands on PC 10–11; a 10-PC embedding drops it at the corners of both axes and a 15-PC one keeps it everywhere (R² 0.96–0.97). This is a property of the gene model, not of any method.
- Metric: per-gene vectors (CYFER and CoSPAR Spearman with the fate score over t1 cells; lineage-DE rank-biserial from ranks, U identical to `wilcox.test(x_H, x_L)$statistic`), Spearman with the truth `Spearman(x_g, Z_true)` over all 2,000 genes; Pearson and Jaccard (BH `q < 0.05`, p from the t-approximation for Spearman and the tie-corrected normal approximation for Wilcoxon) secondary. (The split-half reliability of the truth was removed by decision: its Spearman form is 0.2–0.7 by construction with 1,900 null genes.)
- CoSPAR: `state_info` `"t1"` / `"High"` / `"Low"` by the mean split; shared PCA passed as `X_pca` for all cells (t2 projected with the t1 centre and rotation); fate bias restricted to t1 cells; bulky export files deleted after import. CoSPAR's fate-score correlation with `Z_true` tops out near 0.7 even at the easy end and does not move with `τ_δ`.
- AED is the ratio of mean squared pairwise distances (twice the sum of per-coordinate variances), per clone, in the shared embedding; the paper's Methods formula is to be corrected to square.
- Real-data Gini anchors (0.64 at t1, 0.72–0.81 day 10, ~1 week 5) are single-time-point values on the same scale as the simulated axis.
- `Rscript` parses a script lazily: never edit a driver while it is running.
- The `_claude` suffix marks Claude-drafted files.

## Open Questions / Next Steps
1. Kevin reviews the re-knitted `Writeup21_trailblazing-report_claude.html` (UMAPs, `_r2` curves, permutation table) and the commented `run_cospar.py`.
2. Decide the fixed Gini for the final round: 0.5 (held; stands by default) versus 0.6 (the last easy level on the `_r2` curves under the CoSPAR exclusion); a 20-replicate question since CoSPAR is erratic at Gini 0.7–0.9 with 2 replicates.
3. Final round: 20 replicates on Hyak with the settled settings, which needs a `COSPAR_ENV` built there and a SLURM wrapper in the `sim3_heritability.slurm` pattern; the `τ_δ` supplementary row runs then. Optional if wanted: a many-permutation null for a formal permutation p-value.
4. Fix the paper's AED formula in `PAPER_REPO` to the squared version.
5. Commit the Writeup21 scripts, CSVs and figures on branch `kevin`.
