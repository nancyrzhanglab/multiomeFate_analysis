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
- `kevin/Writeup21_new-simulations/simulation_design_claude.md` holds the design memo for the two Figure 4 sweeps (Gini of clone size; within-clone heterogeneity / AED), with a proposed generator, calibration scheme, matched method pipelines, sanity checks, code plan, and 14 numbered questions for Kevin. **Nothing has been implemented or run.**
- sim1–sim7 (`Writeup_Simulations/`) are complete with RDS outputs in `OUT_ROOT/Writeup_Simulations/` and CSVs in `csv/kevin/Writeup_Simulations/`.
- Uncommitted: the `.DS_Store` deletions, the CLAUDE.md rewrite, and everything above.

## Key Methodological Details
- CYFER's effective sample size is the number of clones; `cyfer()` refuses a fit when features ≥ training clones, so genes always go through an embedding (fastTopics on real data, PCA in synthetic sims).
- `cell_imputed_score` is log10(expected progeny); `coefficient_vec` is natural-log; `lineage_imputed_count` is counts. `exp()` on the score is the documented mistake.
- The package accepts zero-count clones; every sim1–sim7 script filters them out. The Writeup21 memo recommends keeping them (they carry the signal at high Gini) pending Kevin's answer.
- Genes are implicated by Spearman correlation with the fate potential + BH, never from `β̂` (paper Methods). CoSPAR's own recipe is progenitor split (bias > 0.6 / < 0.4) then Wilcoxon + BH.
- CoSPAR fate bias is 0.5 (not NA) for cells outside the transition map; restrict to t1 cells before correlating. The obs column name carries a literal `*` (`fate_bias_..._High*Low`).
- The paper's AED formula is a ratio of mean *un-squared* pairwise distances despite the name; reviewer 2 flagged it. Under the memo's generator the squared version equals `1 − h²`.
- The earlier priming/plastic semi-synthetic sims reach Gini 0.26/0.28 only; real data run 0.64 (pre-treatment) to ~1 (week 5). The new sweeps must bracket that.
- The `_claude` suffix marks Claude-drafted files; `Writeup21` scripts will follow it.

## Open Questions / Next Steps
1. Kevin answers the 14 questions in `simulation_design_claude.md` §9 (Gini definition, mechanism for skew, lineage-DE at t1 vs t2, progeny inheritance `ρ`, truth-set definition, zero-count clones, scale/compute, embedding, gene-calling operating point, real-data AED/Gini anchors, decoupling, CoSPAR fate labels, control row, output conventions).
2. Then implement `func_generate_claude.R`, `func_methods_claude.R`, `sim_gini_claude.R`, `sim_aed_claude.R` in `Writeup21_new-simulations/`, copying the two `run-cospar` templates beside them.
3. Reinstall `multiomeFate` from `CYFER_PKG` (source 1.0.3.000 vs installed 1.0.1.0) before any fit.
4. Decide laptop (`mclapply`) versus Biostat cluster for the ~10 CPU-hour sweeps; the cluster has no `COSPAR_ENV` yet.
5. Commit the retrofit and the two skills on branch `kevin`.
