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
| `CYFER_PKG` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/multiomeFate/` | branch `master`; source is 1.0.3.000, installed copy was 1.0.1.0 on 2026-09-22 (reinstall before sims) |
| `PAPER_REPO` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/multiome_fate_paper/` | also reachable as `/Users/kevinlin/Dropbox/...` (same folder); `paper_nbt.tex` is the newer source; no CLAUDE.md there |
| `PAPERS_DIR` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/papers/` | `multiome_fate_paper.pdf`, `simulation-planning.txt` |
| `OUT_ROOT` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/out/` | `Writeup_Simulations/sim1..7_*.rds` present; `Writeup14/` RData + CoSPAR obs CSVs present |
| `OUT_ROOT` | UW Biostat cluster (SLURM, account `biostat`) | `/home/users/kzlin/kzlinlab/projects/multiomeFate/out/kevin/` | what `sim3_heritability.slurm` and the sim scripts' `~/kzlinlab/...` paths resolve to |
| `OUT_ROOT` | Wharton HPC (legacy) | `/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/` | Writeup14 and earlier were run here; treat as archival |
| `COSPAR_SRC` | personal laptop (macOS) | `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/Nancy/multiomeFate/git/cospar/` | version 0.5.0 checkout; tutorials in `docs/source/*.ipynb` |
| `COSPAR_ENV` | personal laptop (macOS) | `/opt/miniconda3/envs/cospar/` | cospar 0.4.1, Python 3.9.20, scanpy 1.10.3; smoke-tested 2026-09-22 |
| `COSPAR_ENV` | UW Biostat cluster | *(not present)* | no cospar environment built there yet |

## Project Status (as of 2026-09-22)
- Repo retrofitted to the collaboration template: master `CLAUDE.md` rewritten with named external locations (no machine paths), `.githooks/` large-file hook installed and `core.hooksPath` set, `.gitignore` merged with the lab template plus `cospar_cache/` and `cospar_out/`, `additional_context/summary.md` created indexing the CoSPAR paper and the figure mockup.
- Two project-local skills written under `.claude/skills/`: `run-cyfer` (input contract, scales, identifiability floor, downstream recipes) and `run-cospar` (flat-file R export → Python runner in `COSPAR_ENV` → R import, with templates). The CoSPAR pipeline was smoke-tested end to end on a toy dataset (1200 cells, 30 clones, 200 genes): 10/10 causal genes recovered via CoSPAR's own progenitor-DE route in ~15 s.
- `kevin/Writeup21_new-simulations/simulation_design_claude.md` is the design memo for the two Figure 4 sweeps, revised three times against Kevin's written responses (`additional_context/responses-to-simulation-plan_round{1,2,3}.txt`; the between-round changes are logged in `HISTORY_kevin.md`, not in the memo). Core decisions are settled (see Key Methodological Details); §9 holds 2 open questions (keep the non-causal spread placement? confirm `τ_δ = 0.6`?), neither of which blocks implementation. `demo_aed_range_claude.R` (generation + PCA only, both spread placements, ~7.5 min) has been run; its CSV is in `csv/kevin/Writeup21/`. **No CYFER or CoSPAR fit has been run for Writeup21.**
- **Memo revision rule (Kevin, 2026-09-22):** each revision of `simulation_design_claude.md` must be self-contained, with no mention of previous rounds (no change log, no "Kevin's round-N answer" citations, no round-numbered question lists); the changes between rounds go into `HISTORY_kevin.md`.
- sim1–sim7 (`Writeup_Simulations/`) are complete with RDS outputs in `OUT_ROOT/Writeup_Simulations/` and CSVs in `csv/kevin/Writeup_Simulations/`.
- `fig/kevin/Writeup21/` exists and is empty.
- Uncommitted: the CLAUDE.md edits, the skills, the memo, the demo script and CSV, and all three responses files.

## Key Methodological Details
- CYFER's effective sample size is the number of clones; `cyfer()` refuses a fit when features ≥ training clones, so genes always go through an embedding (fastTopics on real data, PCA in synthetic sims).
- `cell_imputed_score` is log10(expected progeny); `coefficient_vec` is natural-log; `lineage_imputed_count` is counts. `exp()` on the score is the documented mistake.
- The package accepts zero-count clones; every sim1–sim7 script filters them out. Writeup21 keeps them by decision (they carry the signal at high Gini).
- Writeup21 generator, as settled: 100 clones × 10 t1 cells (equal); 3,000 expected t2 cells, zeros kept, no cap on clone size; latent `d = 10` with clone centres `τ = s·sqrt(h²)`, typical within-clone spread `σ_w = s·sqrt(1 − h²)`, and clone-varying spread `σ_l² = σ_w² g_l`, `g_l ~ Gamma(1.5, 1.5)` fixed across levels; a `spread_variation` switch puts `g_l` on all coordinates (`"isotropic"`) or on the 9 non-causal ones only (`"noncausal"`, the working default pending Kevin); `β = 1.5` on the first coordinate; progeny `ρ = 0.8` plus a t2 shift `δ + δ_l` (shared, large; clone-specific with `τ_δ = 0.6` fixed across levels in the trailblazing round, supplementary row `{0, 0.3, 0.6, 1.2}` in the final round); shared 10-PC PCA on log-normalized counts for CYFER, CoSPAR and AED.
- Writeup21 axes: t2-only Gini `Gini(Y_l)` (floor ~0.1 from Poisson noise, ceiling 0.99); mean *squared* AED over clones, which equals `1 − h²` to within a few hundredths and is bounded by 1, while the per-clone AED spans ~0.05–2 at the top level with `κ = 1.5`. Calibration: `h² = 1 − AED target` directly; `s` bisected for the Gini target (shrinking `s` lowers the Gini without moving the AED). One knob per axis, no iteration.
- Spread-size coupling (Jensen's inequality on the log-normal mean, `E[Y_l] ∝ exp(β² σ_l²/2)`): with isotropic clone-varying spread, diffuse clones expand more (Spearman(AED_l, Y_l) up to ~0.6 at the top AED level). Under `"noncausal"` it is ~0 (0.07 latent-only, within noise in the PCA demo) and the per-clone AED range is unchanged, but the Gini then *falls* along the AED axis at fixed `s` (0.57 → 0.39) instead of rising, so the calibrated `s` climbs to ~1.25 at the top AED level.
- Provisional middle point (`h² = 0.6`, t2 Gini 0.5): latent-only bisection gives `s ≈ 0.81`, `τ ≈ 0.63`, `σ_w ≈ 0.51` (`"noncausal"`); `τ_δ = 0.6` is that `τ`, rounded, held as an absolute value.
- Writeup21 metric: per-gene association vector per method (CYFER: Spearman with `Z_hat`; CoSPAR: Spearman with `fate_bias`, matched route only; lineage-DE: rank-biserial `2U/(n_H n_L) − 1` from the Wilcoxon on high-vs-low t1 cells, clones split at the *mean* t2 size), correlated by **Spearman** with the truth vector `Spearman(x_g, Z_true)` over all 2,000 genes. Pearson and Jaccard are secondary columns. The rank-biserial is preferred over `U` for cross-level comparability, sign, and effect-size interpretation; the Spearman headline is identical either way (memo §6.3).
- Two rounds: trailblazing (laptop, 7 levels × 2 replicates, provisional fixed values t2 Gini 0.5 and mean AED 0.4, `τ_δ = 0.6`; levels Gini `{0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9}`, mean AED `{0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0}`) to find the last level where all three methods do well, then the final round with those as the fixed values, re-spaced levels, 20 replicates on Hyak, and the `τ_δ` supplementary row. If CoSPAR fails even at the easy end, lower `τ_δ` first.
- Verified in `COSPAR_SRC`: the multi-time-clone map is `S_t1 · M · S_t2` (barcode link `M`, within-time similarity blocks); the t1-by-t2 similarity is never read, so a rigid t2 shift `δ` is inert and "t2 far from t1" is not barcoded CoSPAR's assumption. What breaks it: one-clone "High" fates at high Gini, uniform clone-level links under within-clone heterogeneity, ignoring extinct clones, and loss of cross-clone coherence at t2, dialled by `δ_l`.
- Genes are implicated by Spearman correlation with the fate potential, never from `β̂` (paper Methods); `β̂` lives on the PCs.
- CoSPAR fate bias is 0.5 (not NA) for cells outside the transition map; restrict to t1 cells before correlating. Clones with no t2 cells are single-time clones to CoSPAR (no barcode link, bias from smoothing only); not fatal.
- AED is the ratio of mean *squared* pairwise distances in the PCA embedding, per clone (twice the sum of per-coordinate variances, no pairs formed); the paper's Methods formula is to be corrected to square.
- Real-data Gini anchors (0.64 at t1, 0.72–0.81 day 10, ~1 week 5) are single-time-point values, now on the same scale as the simulated axis.
- The `_claude` suffix marks Claude-drafted files; `Writeup21` scripts follow it.

## Open Questions / Next Steps
1. Kevin answers the 2 questions in `simulation_design_claude.md` §9 (keep `spread_variation = "noncausal"`, suggested yes; confirm `τ_δ = 0.6`). Neither blocks implementation; defaults are stated in the memo.
2. Implement `func_generate_claude.R` (lift the generator from `demo_aed_range_claude.R` including the `spread_variation` switch, add t2 children and the `(h², s)` calibration by latent-only bisection), `func_methods_claude.R`, `sim_gini_claude.R`, `sim_aed_claude.R`, `make_csvs_claude.R`, `plot_barplots_claude.R` in `Writeup21_new-simulations/`, copying the two `run-cospar` templates beside them; drivers append a time-stamped progress `.txt` under `OUT_ROOT/Writeup21_new-simulations/`.
3. Trailblazing round on the laptop: 7 levels × 2 replicates per axis (~1–2 h sequential); CSVs to `csv/kevin/Writeup21/`, barplots to `fig/kevin/Writeup21/`; read off the last easy level per axis.
4. Reinstall `multiomeFate` from `CYFER_PKG` (source 1.0.3.000 vs installed 1.0.1.0) before any fit.
5. Final round: 20 replicates on Hyak with the chosen fixed values, which needs a `COSPAR_ENV` built there; the `τ_δ` supplementary row runs then.
6. Fix the paper's AED formula in `PAPER_REPO` to the squared version.
7. Commit the retrofit, the two skills, the memo, the demo and the three responses files on branch `kevin`.
