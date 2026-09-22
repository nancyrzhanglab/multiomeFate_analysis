# HISTORY_kevin.md — Kevin's Session Log

> **Append-only, ascending chronological order** (oldest at top, newest at the bottom). Add each session's dated entry to the END of this file. Never read at session startup — consulted only on demand for deep history. Current project state lives in `CLAUDE_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session appends to this file. Other collaborators may read it but must not edit or rewrite entries.

---

### 2026-05-11 (before this log existed — reconstructed from git)
- sim1–sim7 revision simulations written, run, and exported to CSV (commits `a9b9449` through `3d5144e`); `sim3_heritability.slurm` added for the Biostat cluster.
- The original `CLAUDE.md` carried absolute laptop paths and predated the external-locations convention.

### 2026-09-22 (Session 1 — retrofit, skills, Writeup21 design memo)
- Retrofitted the repo with `project-setup`: master `CLAUDE.md` rewritten around named external locations (`CYFER_PKG`, `PAPER_REPO`, `PAPERS_DIR`, `OUT_ROOT`, `COSPAR_SRC`, `COSPAR_ENV`); per-machine paths moved into `CLAUDE_kevin.md`; `.githooks/` + `core.hooksPath` installed; `.gitignore` merged with the lab template; `additional_context/summary.md` created.
- Decided the two method skills are project-local (`.claude/skills/run-cyfer`, `.claude/skills/run-cospar`) rather than global, because their contract is specific to this project's data shapes.
- Decided `run-cospar` ships its own flat-file templates instead of reusing `seurat-anndata-bridge`: the input is a count matrix plus three per-cell columns, and the earlier Writeup14 pipeline's choices (scaled data as `adata.X`, `pd.get_dummies` clone matrix, DE in R, two contradictory winner/loser thresholds) were all replaced by CoSPAR's own `get_X_clone` → `fate_bias` → `progenitor` → `differential_genes` route.
- Smoke-tested the CoSPAR round trip in `/opt/miniconda3/envs/cospar` (cospar 0.4.1): 1200 cells, 30 clones, 200 genes, 10 causal; Jaccard 1.0 against truth; similarity cache is ~38 MB per dataset, hence the new `.gitignore` entries.
- Found the installed `multiomeFate` (1.0.1.0) is behind the source (1.0.3.000); recorded as a pre-flight step in the `run-cyfer` skill.
- Wrote `kevin/Writeup21_new-simulations/simulation_design_claude.md`: one hierarchical latent generator with two knobs (`τ` between-clone, `σ_w` within-clone), calibrated by bisection to hit target Gini / AED levels while holding the other axis fixed; CYFER, CoSPAR (native and matched gene-calling routes), lineage-DE (t1 and t2 variants), and an oracle; Jaccard at BH q<0.05 plus top-k / AUROC.
- Open: the 14 design questions in the memo, chief among them the Gini definition (pooled vs t2-only), the mechanism generating skew, whether lineage-DE runs on t1 or t2 cells, how much progeny inherit their parent's state, and whether zero-count clones enter the CYFER fit.
- Open: the paper's "average-of-squared Euclidean distance" is written un-squared in the Methods; which version the resubmission keeps.
- Noted for later: `sim7` filters to `lineage_future_count > 0` and the package accepts zeros; at high Gini the filter discards most clones.
