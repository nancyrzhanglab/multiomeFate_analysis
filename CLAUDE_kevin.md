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

## Project Status (as of 2026-09-22)
- Writeup21 is implemented and the trailblazing round has run on the laptop. Code in `kevin/Writeup21_new-simulations/`: `func_generate_claude.R`, `func_methods_claude.R`, `func_sweep_claude.R` (shared engine), `sim_gini_claude.R`, `sim_aed_claude.R`, `make_csvs_claude.R`, `plot_barplots_claude.R`, plus the two `run-cospar` templates copied beside them. Drivers take `WRITEUP21_DRY=1` and `WRITEUP21_PCA=<d>`; outputs carry `_dry` / `_pca<d>` suffixes. Per dataset about 1 min (CoSPAR 15–45 s of it); a 7 × 2 sweep is 11–13 min per axis; the two axes run concurrently.
- `kevin/Writeup21_new-simulations/Writeup21_trailblazing-report_claude.{Rmd,html}` is the knitted summary of the trailblazing round (reads the CSVs; knit with `RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`).
- Outputs: RDS and progress files in `OUT_ROOT/Writeup21_new-simulations/` (`sim_{gini,aed}_claude{,_pca20,_dry}.rds`), CSVs in `csv/kevin/Writeup21/` (summary, replicate_details, calibration per axis and pass, plus `check_tau_delta_easy_end_pca20.csv`), figures in `fig/kevin/Writeup21/` (Spearman, Jaccard, fate-score barplots and metric-vs-realized per axis and pass).
- **Headline result (20-PC pass, after the review fixes):** CYFER holds at 0.75–0.97 across both axes except the top Gini level (0.52, two replicates) while CoSPAR falls 0.67 → −0.13 (Gini) and 0.69 → −0.04 (AED) and lineage-DE 0.73 → 0.56 and 0.75 → 0.23; CYFER's fate-score correlation with the truth is 0.74–0.98 at every level. With the memo's 10-PC embedding CYFER instead collapses at the hard end of both axes because the causal latent axis falls to PC 11 (it loads on 100 genes, each non-causal coordinate on all 2,000); that pass is kept on disk (no suffix) as the record, the `_pca20` pass is the figure. The embedding dimension is the one design change to settle.
- Calibration now re-centres `h²` on the realized AED (memo §5); the Gini-0.2 level is unreachable at fixed AED 0.4 (mean AED 0.55 even at `h² = 1`, count noise dominates at `s = 0.21`), so the final ladder should start at 0.3. `"noncausal"` panned out (spread-size coupling within ±0.11); `τ_δ = 0.6` stands (CoSPAR insensitive to it at the easy end). Last easy level on the corrected curves: Gini 0.5 and AED 0.4, the provisional values.
- A fresh-context code review found no arithmetic bugs; its three fixes (AED re-centring, seed sub-streams `10·seed + {0,1,2}` with calibration draws at `+50`, constant CoSPAR bias reported as missing) are in. One fairness item is open for Kevin: t1 cells of extinct clones sit outside CoSPAR's map (`extend_Tmap_space = False`) and get a bias of exactly 0.5, contrary to memo §3.2.
- The memo `simulation_design_claude.md` still says 10 PCs everywhere and has not been revised for these results. Memo revision rule: each revision self-contained, between-round changes in `HISTORY_kevin.md`.
- Installed `multiomeFate` is 1.0.3.0, matching the source, via `R CMD INSTALL` from a scratch copy (`devtools::install_local()` hangs on this laptop).
- sim1–sim7 (`Writeup_Simulations/`) are complete with RDS outputs in `OUT_ROOT/Writeup_Simulations/` and CSVs in `csv/kevin/Writeup_Simulations/`.
- Uncommitted on branch `kevin`: everything from Writeup21 this session (the nine scripts, the CSVs, `fig/kevin/Writeup21/`) and the two per-person files; the retrofit, skills, memo, demo and responses files were committed earlier.

## Key Methodological Details
- CYFER's effective sample size is the number of clones; `cyfer()` refuses a fit when features ≥ training clones, so genes always go through an embedding (fastTopics on real data, PCA in synthetic sims). With 100 clones and 3 folds the cap is about 60 features, so 20 PCs is safe.
- `cell_imputed_score` is log10(expected progeny); `coefficient_vec` is natural-log; `lineage_imputed_count` is counts. `exp()` on the score is the documented mistake.
- The package accepts zero-count clones; every sim1–sim7 script filters them out. Writeup21 keeps them (`method_cyfer()`), and CYFER converged on every dataset with 0–33 extinct clones.
- Writeup21 generator: 100 clones × 10 t1 cells; 3,000 expected t2 cells via `β_0`, zeros kept, no cap; latent `d = 10`, `τ = s·sqrt(h²)`, `σ_w = s·sqrt(1 − h²)`, `g_l ~ Gamma(1.5, 1.5)` on the 9 non-causal coordinates (`spread_variation = "noncausal"`); `β = 1.5`; children `ρ = 0.8` with shared shift `||δ|| = 3s` and clone shift `τ_δ = 0.6`, both on non-causal coordinates; NB counts `θ = 10`, 2,000 genes, 100 causal. Seeds: replicate `r` of level `j` is `seed_base + 100j + r` (1000 for Gini, 2000 for AED); the generator uses `seed`, `seed + 1`, `seed + 2` for the three stochastic stages.
- Calibration: `h² = 1 − AED target`, then `s` bisected (geometric, 10 latent draws at seeds `seed_base + 100j + 50 + k`, tol 0.01) for the Gini target, then one full dataset's realized mean AED checked and `h²` shifted by the miss (up to 3 times, tol 0.03) with `s` re-bisected. At 20 PCs the Gini axis needs `h² = 1, 0.785, 0.66, 0.6, 0.6, 0.6, 0.6` (levels 0.2–0.9) and the AED axis `h² = 1, 0.83, 0.65, 0.45, 0.30, 0.11, 0`; calibrated `s` runs 0.21–2.06 on the Gini axis and 0.62–1.25 on the AED axis. Realized Gini and mean AED are within 0.03 of target except the unreachable Gini-0.2 AED.
- The causal latent axis is the weakest direction in gene space (100 genes × loading² ≈ 0.58 versus 2,000 × 0.0625 per non-causal coordinate), so it lands on PC 10–11; a 10-PC embedding drops it at the corners of both axes and a 15-PC one keeps it everywhere (R² 0.96–0.97). This is a property of the gene model, not of any method.
- Metric: per-gene vectors (CYFER and CoSPAR Spearman with the fate score over t1 cells; lineage-DE rank-biserial from ranks, U identical to `wilcox.test(x_H, x_L)$statistic`), Spearman with the truth `Spearman(x_g, Z_true)` over all 2,000 genes; Pearson and Jaccard (BH `q < 0.05`, p from the t-approximation for Spearman and the tie-corrected normal approximation for Wilcoxon) secondary. The Spearman split-half reliability of the truth is low by construction (most genes null); the Pearson version (0.58–0.93) is the usable reliability.
- CoSPAR: `state_info` `"t1"` / `"High"` / `"Low"` by the mean split; shared PCA passed as `X_pca` for all cells (t2 projected with the t1 centre and rotation); fate bias restricted to t1 cells; bulky export files deleted after import. CoSPAR's fate-score correlation with `Z_true` tops out near 0.7 even at the easy end and does not move with `τ_δ`.
- AED is the ratio of mean squared pairwise distances (twice the sum of per-coordinate variances), per clone, in the shared embedding; the paper's Methods formula is to be corrected to square.
- Real-data Gini anchors (0.64 at t1, 0.72–0.81 day 10, ~1 week 5) are single-time-point values on the same scale as the simulated axis.
- `Rscript` parses a script lazily: never edit a driver while it is running.
- The `_claude` suffix marks Claude-drafted files.

## Open Questions / Next Steps
1. Decide the shared-embedding dimension for the final round: 20 PCs (recommended; the re-centring keeps the AED axis on target there) versus changing the gene model so 10 PCs recover the latent space. Revise `simulation_design_claude.md` (self-contained) for this, the re-centring step as implemented, the seed scheme, and the CoSPAR extinct-clone sentence in §3.2 (or set `extend_Tmap_space = True` in `run_cospar.py`); record the changes in `HISTORY_kevin.md`.
2. Confirm the final-round fixed values (the corrected curves point to the provisional Gini 0.5 and AED 0.4), drop Gini 0.2 from the ladder (unreachable at fixed AED 0.4) and consider re-spacing the top of the Gini axis (the comparators die between 0.5 and 0.9; CYFER itself drops at 0.9).
3. Run the memo's §7 label-permutation check at the middle level of each axis (the other §7 quantities are already columns in the replicate CSVs).
4. Final round: 20 replicates on Hyak with the chosen settings, which needs a `COSPAR_ENV` built there and a SLURM wrapper in the `sim3_heritability.slurm` pattern; the `τ_δ` supplementary row runs then.
5. Fix the paper's AED formula in `PAPER_REPO` to the squared version.
6. Commit the Writeup21 scripts, CSVs and figures on branch `kevin`.
