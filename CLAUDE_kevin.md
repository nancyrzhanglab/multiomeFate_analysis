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
- `kevin/Writeup21_new-simulations/simulation_design_claude.md` is the design memo for the two Figure 4 sweeps, revised once against Kevin's round-1 responses (`additional_context/responses-to-simulation-plan_round1.txt`). The core decisions are settled (see Key Methodological Details); §9 holds 8 new questions, the first two of which (pooled-Gini ceiling, AED range) change the level lists. **Nothing has been implemented or run.**
- sim1–sim7 (`Writeup_Simulations/`) are complete with RDS outputs in `OUT_ROOT/Writeup_Simulations/` and CSVs in `csv/kevin/Writeup_Simulations/`.
- `csv/kevin/Writeup21/` and `fig/kevin/Writeup21/` exist and are empty.
- Uncommitted: the `.DS_Store` deletions, the CLAUDE.md rewrite, the skills, the memo, and the responses file.

## Key Methodological Details
- CYFER's effective sample size is the number of clones; `cyfer()` refuses a fit when features ≥ training clones, so genes always go through an embedding (fastTopics on real data, PCA in synthetic sims).
- `cell_imputed_score` is log10(expected progeny); `coefficient_vec` is natural-log; `lineage_imputed_count` is counts. `exp()` on the score is the documented mistake.
- The package accepts zero-count clones; every sim1–sim7 script filters them out. Writeup21 keeps them by decision (they carry the signal at high Gini).
- Writeup21 design, as settled: 100 clones × 10 t1 cells (equal); ~3,000 t2 cells total with zeros; sweep `τ` for Gini and `σ_w/τ` for AED, bisecting so the other statistic stays fixed (pooled Gini 0.4, mean AED 0.5); progeny `ρ = 0.8` plus a large rigid t2 shift; shared 10-PC PCA on log-normalized counts for CYFER, CoSPAR and AED.
- Writeup21 metric: per-gene association vector per method (CYFER: Spearman with `Z_hat`; CoSPAR: Spearman with `fate_bias`, matched route only; lineage-DE: signed rank-biserial from high-vs-low t1 cells, clones split at the *mean* t2 size), correlated with the truth vector `Spearman(x_g, Z_true)` over all 2,000 genes. Jaccard is a secondary column only.
- Pooled Gini with equal t1 sizes `c` equals `Gini(Y) · mean(Y)/(c + mean(Y))`, a hard ceiling (~0.74 at 10 t1 cells and 3,000 t2 cells). Mean AED over clones is bounded by ~1 under a shared `σ_w` and has a count-noise floor.
- A rigid t2 shift `δ` reaches none of the three gene calls (all test t1 cells; CoSPAR smooths within time point), so it is not what breaks CoSPAR; what does is one-clone "High" fates at high Gini, clone-level bias at high AED, and ignoring extinct clones.
- Genes are implicated by Spearman correlation with the fate potential, never from `β̂` (paper Methods); `β̂` lives on the PCs.
- CoSPAR fate bias is 0.5 (not NA) for cells outside the transition map; restrict to t1 cells before correlating. Clones with no t2 cells are single-time clones to CoSPAR (no barcode link, bias from smoothing only); not fatal.
- The paper's AED is the ratio of mean *un-squared* pairwise distances in the PCA embedding, per clone, despite the name; Writeup21 uses that version and the paper will be fixed to match.
- Real-data Gini anchors (0.64 at t1, 0.72–0.81 day 10, ~1 week 5) are single-time-point values, not pooled; both are recorded per dataset so the figure can be relabelled.
- The `_claude` suffix marks Claude-drafted files; `Writeup21` scripts will follow it.

## Open Questions / Next Steps
1. Kevin answers the 8 new questions in `simulation_design_claude.md` §9: pooled-Gini ceiling / t2 total / whether 500 is a hard cap on clone size; AED range (mean ≤ 1 vs per-clone 0–2); CYFER per-gene statistic; lineage-DE signed statistic; Pearson vs Spearman outer correlation; rigid vs clone-specific `δ`; the fixed middle values; AED embedding dimension.
2. Then implement `func_generate_claude.R`, `func_methods_claude.R`, `sim_gini_claude.R`, `sim_aed_claude.R`, `make_csvs_claude.R`, `plot_barplots_claude.R` in `Writeup21_new-simulations/`, copying the two `run-cospar` templates beside them; drivers append a time-stamped progress `.txt` under `OUT_ROOT/Writeup21_new-simulations/`.
3. Pilot on the laptop: 7 levels × 2 replicates per axis (~1–2 h sequential), a pilot first to find the AED floor; CSVs to `csv/kevin/Writeup21/`, barplots to `fig/kevin/Writeup21/`.
4. Reinstall `multiomeFate` from `CYFER_PKG` (source 1.0.3.000 vs installed 1.0.1.0) before any fit.
5. Later: 20 replicates on Hyak, which needs a `COSPAR_ENV` built there.
6. Commit the retrofit, the two skills, the memo and the responses file on branch `kevin`.
