# CLAUDE.md — multiomeFate Analysis Project

## Workflow Instructions
1. **Always enter plan mode** before starting any non-trivial task.
2. **Use superpower skills** where relevant: `/code-review` for code changes, `/r-style-guide` for any `.R` file.
3. **After every prompt**, run `/project-state`: refresh the current-state sections of the individual contributor's `CLAUDE_[name].md` in place, and append a dated entry to their `HISTORY_[name].md`. Write only non-obvious things; skip anything already in the code or git history.

## Project Context (High Level)
**Paper/Project**: Analysis repository for *"Resolution of Selection Versus Adaptation in Cellular Evolution"* (Chen, Lin, Schaff, et al.), which introduces **CYFER** (Clonal Fate Estimation by Exponential Regression): per-cell fate potential from single-cell multiome data paired with lineage barcoding.
**Authors**: Xinyi E. Chen, Kevin Z. Lin, Dylan Schaff (co-first); Robert Vander Velde, Christopher Cote, Sijia Huang; Andy J. Minn, Sydney M. Shaffer and Nancy R. Zhang (co-corresponding). Kevin owns this repo's `kevin/` tree; Chi-Yun's earlier scripts sit in `ChiYun/`.
**Goal**: Distinguish selection (pre-existing high-potential cells expand) from adaptation (surviving cells change state), and link each to molecular features.

**Submission status**: rejected at Nature Biotechnology (the version in `PAPER_REPO`), reviewed and rejected at Nature Genetics (January 2026; report and response letter in `PAPER_REPO/additional_context/`), now **being prepared for Nature Methods**. `kevin/Writeup_Simulations/` (sim1–sim7) answered the Nature Genetics critiques; the Nature Methods figure plan is the first seven slides of `additional_context/Mockup of figures of Nancy-Sydney paper.pptx`, whose **Figure 4** is the two new simulation sweeps designed in `kevin/Writeup21_new-simulations/`.

This is one of three sibling repositories (see External Locations): this analysis repo (GitHub `nancyrzhanglab/multiomeFate_analysis`, branch `kevin`), the R package `CYFER_PKG` (GitHub `nancyrzhanglab/multiomeFate`), and the manuscript `PAPER_REPO`. Read `CYFER_PKG/CLAUDE.md` for how CYFER is *implemented*; this file covers how it is *applied*.

## Repository Layout

```
multiomeFate_analysis/
  CLAUDE.md, CLAUDE_[name].md, HISTORY_[name].md
  .claude/skills/run-cyfer/       ← local skill: fitting CYFER and reading its outputs
  .claude/skills/run-cospar/      ← local skill: R → CoSPAR (Python) → R, with templates
  additional_context/             ← reference PDFs and the figure mockup; index in summary.md
  kevin/
    Writeup21_new-simulations/    ← CURRENT FOCUS: Gini and heterogeneity sweeps (design memo only so far)
    Writeup_Simulations/          ← sim1–sim7 revision simulations (Nature Genetics critiques)
    Writeup14_simulation/         ← semi-synthetic priming/plastic simulations + earlier CoSPAR runs
    Writeup17b_simulation-plots/  ← figure code for Writeup14 (Jaccard barplots, UMAPs, Gini)
    Writeup20_vignette-creating/  ← script that produced the package's bundled datasets
    Writeup*/                     ← numbered analysis writeups (real data; mostly run on the cluster)
    analysis_pipeline/
  csv/kevin/Writeup_Simulations/  ← flat CSVs exported by make_csvs.R (and Writeup14, Writeup18)
  fig/kevin/                      ← figures by writeup
  ChiYun/, demo/                  ← collaborator scripts and package demos
```

## External Locations
Folders this project depends on that live **outside** the project root. A path is external unless it can be written relative to this file's directory without `..`.

**No per-machine filesystem path appears in this file.** Each location is named here with its purpose; the path, and which machine it is on, is recorded per person in `CLAUDE_[name].md` under *External Locations (per-machine paths)*. Hostnames and URIs that are the same for everyone are fine here.

| Location name | Purpose | Copy semantics |
|---|---|---|
| `CYFER_PKG` | The `multiomeFate` R package source (CYFER). Install from it with `devtools::install_local()` before running any sim; the installed copy drifts behind the source | per-person copy (git clone; GitHub `nancyrzhanglab/multiomeFate`) |
| `PAPER_REPO` | Manuscript LaTeX (`paper.tex`, `paper_nbt.tex`), `fig/`, and `additional_context/` with the Nature Genetics review and response | per-person copy (Dropbox / local) |
| `PAPERS_DIR` | Paper PDF, `simulation-planning.txt` (reviewer critiques with planned responses), reviewer notes | per-person copy |
| `OUT_ROOT` | Large outputs: `Writeup_Simulations/*.rds` from sim1–sim7, `Writeup14/` RData and CoSPAR obs CSVs, future `Writeup21_new-simulations/` outputs and CoSPAR exports. Never tracked in git | per-person copy; a laptop copy and a cluster copy exist and drift |
| `COSPAR_SRC` | Checkout of the CoSPAR Python package (Wang et al. 2022), for reading the API and tutorials | per-person copy |
| `COSPAR_ENV` | A conda environment with `cospar` importable; the `run-cospar` skill invokes its `bin/python` by absolute path | per-person copy |

Refer to these by name in prose, code comments, and session notes. To resolve a name, read the current user's `CLAUDE_[name].md`; with no row, ask.

## Who Is Using This Session?
**Detect the current user** by running: `echo $USER`. This table maps each login to that person's **first-name** context file; it is the source of truth that `/brainstorm` and `/project-state` use to resolve `$USER` to the right filename (so the login `kevinlin` maps to `CLAUDE_kevin.md`, never `CLAUDE_kevinlin.md`).

| Username (login) | Current-state file (first name) | History archive |
|---|---|---|
| `kevinlin` | `CLAUDE_kevin.md` | `HISTORY_kevin.md` |

Chi-Yun and Sijia have contributed scripts (`ChiYun/`, `csv/kevin/Writeup14/sijia_*`) but have no row yet; they add one on their first session.

**File ownership.** Each row above names one person's files, and **only that person's session writes them.** Once `$USER` resolves to a first name, that is the only suffix you may create, edit, append to, rename, or delete — every other collaborator's `CLAUDE_[name].md`, `HISTORY_[name].md`, and `brainstorming_[name].md` is read-only. Read them for context when useful; never modify them, not even to fix a typo or add a cross-reference. This repository is shared, so an edit lands in the owner's working copy immediately and can overwrite state they wrote from their own machine. If a collaborator's file looks wrong or stale, say so instead of editing it. This master `CLAUDE.md` is the exception: it is shared and any collaborator may update it. The single override is the user, in the current turn, *directing you to write* that exact file ("add this to `CLAUDE_sarah.md`") — confirm once, then write it. Merely mentioning a collaborator, referring to their file, or being away does not qualify, and attribution does not launder the edit — a change stamped with your name and confined to two lines is still a write to a file you do not own. Record the decision here in the master `CLAUDE.md` and tell the user to contact the owner directly.

**Session startup — run this before any other work.** All four cases below are normal; none is an error to report back to the user.

1. Run `echo $USER` and look for a matching row.
2. **Row exists and the file exists** → read that person's `CLAUDE_[name].md` immediately. It holds the current project state and restores context in ~30 seconds. Do **not** read `HISTORY_[name].md` at startup; it is the append-only session log, consulted only on demand when deep history is needed.
3. **Row exists but `CLAUDE_[name].md` does not** → this is a collaborator's first session. Initialize `CLAUDE_[name].md` and `HISTORY_[name].md` from their templates via `/project-state`, then continue. Do not fall back to reading someone else's file, and do not treat the missing file as a reason to stop.
4. **No matching row** → ask the user their first name, add a row for them to the table above (this file is shared and any collaborator may update it), then initialize their files as in case 3. Never guess a first name from the login, and never write to a suffix you have not confirmed.

Adding a row is a project-level change, so record it the same way as any other (`/project-state`).

## Local skills
Two project-level skills live in `.claude/skills/` and are the reference for running the two methods; read them before writing any script that fits CYFER or calls CoSPAR.

| Skill | Fires when | Ships |
|---|---|---|
| `run-cyfer` | fitting `cyfer()` / `cyfer_finalize()`, scoring cells, ranking genes, computing the selection or adaptation index, or debugging a CYFER error | input contract, scales (`cell_imputed_score` is log10; coefficients are natural log), the identifiability floor (features < training clones), downstream recipes |
| `run-cospar` | benchmarking against CoSPAR on data held in R | `templates/cospar_flat_io.R` (export/import) and `templates/run_cospar.py` (transition map → fate bias → progenitor → differential genes), smoke-tested end to end |

The global `seurat-anndata-bridge` skill is the general R↔Python bridge; `run-cospar` carries its own lighter flat-file contract because its input is a count matrix plus three per-cell columns, not a Seurat object.

## CYFER: Core Method

**Model**: Per-cell fate potential `Z_i = β_0 + X_i^T β`
Clone count at t2: `Y_ℓ ~ Poisson(Σ_{i∈ℓ} exp(Z_i))`

**Loss (Poisson-log)**:
```
L(β) = (1/|clones|)[Σ_i exp(x_i^T β) − Σ_ℓ Y_ℓ log(Σ_{i∈ℓ} exp(x_i^T β))] + λ||β||²
```

**Optimization**: BFGS with analytical gradient, multiple random initializations.
**Regularization**: K-fold CV over a decreasing λ-sequence, folds hold out whole clones.

**Key outputs**:
- `fate_potential`: `cell_imputed_score`, log10(expected future progeny) per cell
- `selection_index`: Gini coefficient of estimated fate potentials (per clone, min 10 cells on real data; or population-wide)
- `adaptation_index`: fate-potential-weighted average expression distance t1→t2, normalized by a "no adaptation" distance

Neither index is in the package; `sim6_adaptation_index.R` and the `gini_coef()` below are the reference implementations. Genes are never read off `β̂`; they are implicated by Spearman correlation with the fate potential, BH-adjusted (paper Methods).

## Package API

```r
library(multiomeFate)

# Step 1: cross-validated fit
fit_res <- cyfer(
  cell_features        = X,          # matrix: cells × features, rownames + colnames, scaled, no constant column
  cell_lineage         = clone_vec,  # character vector of clone IDs
  lineage_future_count = lfc,        # named numeric vector: clone → t2 count (zeros allowed)
  lambda_initial       = 1,
  lambda_sequence_length = 10,
  num_folds            = 3,
  verbose              = 0
)

# Step 2: refit on full data at best lambda
final_fit <- cyfer_finalize(
  cell_features        = X,
  cell_lineage         = clone_vec,
  fit_res              = fit_res,
  lineage_future_count = lfc
)

# Step 3: per-cell fate potential (log10 scale) is returned directly
Z_hat <- final_fit$cell_imputed_score
# or by hand for new cells: log10(exp(X %*% coefficient_vec[-1] + coefficient_vec[1]))
```

**Preprocessing conventions** used by the sim scripts (the package does not enforce them):
- Filter to clones with `>= 2` cells for stable CV.
- The sim1–sim7 scripts also filter to `lineage_future_count > 0`; the package accepts zeros and Writeup21 plans to keep them (see the design memo).
- `num_folds = min(3, n_clones - 1)`, clamped to at least 2; the paper used 20 folds on real data.
- `scale()` the feature matrix before fitting.
- **Features must number fewer than the training clones** (`n_clones - ceiling(n_clones/K) >= p + 1`), so genes go through an embedding (fastTopics in the paper; PCA in synthetic sims) first.

**Defensive fitting pattern** (used across sim scripts, sometimes inlined):
```r
fit_cyfer_safe <- function(X, clone_labels, lfc) {
  valid <- names(lfc[lfc > 0])
  if (length(valid) < 3) return(NULL)
  keep <- which(clone_labels %in% valid)
  X_s <- X[keep, , drop=FALSE]; cl_s <- clone_labels[keep]; lfc_s <- lfc[valid]
  cs <- table(cl_s); valid2 <- names(cs[cs >= 2])
  if (length(valid2) < 3) return(NULL)
  k2 <- which(cl_s %in% valid2)
  X_s <- X_s[k2, , drop=FALSE]; cl_s <- cl_s[k2]; lfc_s <- lfc_s[valid2]
  nf <- min(3, length(unique(cl_s)) - 1); if (nf < 2) return(NULL)
  tryCatch({
    fr <- cyfer(X_s, cl_s, lfc_s, lambda_initial=1,
                lambda_sequence_length=10, num_folds=nf, verbose=0)
    cyfer_finalize(X_s, cl_s, fr, lfc_s)
  }, error = function(e) NULL)
}
```
Note: this pattern appears as `fit_cyfer` (sim5), `fit_cyfer_simple` (sim6), or inlined (sim1–4, sim7).

## Simulation Scripts (Writeup_Simulations/)

Seven scripts address specific reviewer critiques from the Nature Genetics review.
All use synthetic data + `multiomeFate` + `MASS` packages.

| Script | Addresses | Key Variable | n_replicates |
|--------|-----------|--------------|--------------|
| `sim1_growth_modes.R` | Rev 1a | growth model (exp/linear/logistic/power) | 20 |
| `sim2_rare_resistance.R` | Rev 1b, 9 | resistance fraction f ∈ {1%,2%,5%,10%,25%,50%} | 10 |
| `sim3_heritability.R` | Rev 1c, 2, 3 | h2_feature × h2_fate grid (5×4) | 10 |
| `sim4_power_analysis.R` | Rev 1d, 4, 5, 9 | #clones / cells-per-clone / capture rate | 10 |
| `sim5_barcode_dropout.R` | Rev 4, 5, 6 | p_overlap, p_cell_capture, min clone size | 10 |
| `sim6_adaptation_index.R` | Rev 8 | f_dying ∈ {10%,25%,50%,75%,90%,99%} | 10 |
| `sim7_sensitivity_specificity.R` | Rev 7, 1 | priming vs. plasticity scenario | 10 |

**Run**: `Rscript sim1_growth_modes.R` (single-core; see `README.md` there for runtimes ~5–40 min each). `sim3_heritability.slurm` is the cluster pattern (account `biostat`, partition `all-12c128g`). The scripts hardcode their output path; it resolves to `OUT_ROOT/Writeup_Simulations/` on whichever machine they run.

**Parallelization**: wrap replicate loops with `parallel::mclapply()` to speed up.

### RDS output structure (per script)

Each script saves a different structure — there is no uniform schema:

| Script | RDS keys |
|--------|----------|
| `sim1` | `detailed` (single-run metrics), `replicate_cors` (matrix), `summary`, `params` |
| `sim2` | `all_results` (list of per-f data frames), `summary`, `params` |
| `sim3` | `summary` (20-row grid), `grid`, `params` — **no replicate-level data stored** |
| `sim4` | `axis1_clones`, `axis2_cellsize`, `axis3_capture`, `params` |
| `sim5` | `partA` (clone dropout), `partB` (cell subsampling), `partC` (Gini vs clone size), `params` |
| `sim6` | `all_results` (nested list per f_dying), `summary`, `params` |
| `sim7` | `summary`, `permutation`, `all_results`, `null_fpr`, `params` |

## CSV Export (make_csvs.R)

`kevin/Writeup_Simulations/make_csvs.R` reads each RDS and writes flat CSVs to
`csv/kevin/Writeup_Simulations/`. Run: `Rscript make_csvs.R`.

| CSV file | Contents |
|----------|----------|
| `sim1_summary.csv` | Mean/SD/median correlation per growth model |
| `sim1_detailed_single_run.csv` | Correlation, Jaccard, Gini per model (single run) |
| `sim1_replicate_cors.csv` | Per-replicate correlation, one row per replicate |
| `sim2_summary.csv` | Mean AUROC, sensitivity, specificity per resistance fraction |
| `sim2_replicate_details.csv` | Per-(fraction, replicate) metrics |
| `sim3_summary.csv` | Per-(h2_feature, h2_fate) grid cell: CYFER vs naive correlation/Jaccard |
| `sim4_axis1_n_clones.csv` | Performance vs number of clones |
| `sim4_axis2_cells_per_clone.csv` | Performance vs cells per clone |
| `sim4_axis3_capture_rate.csv` | Performance vs barcode capture rate |
| `sim5_partA_clone_dropout.csv` | Correlation and Gini error vs p_overlap (random + size-biased) |
| `sim5_partB_cell_subsampling.csv` | Correlation and Gini error vs p_cell_capture |
| `sim5_partC_gini_vs_clone_size.csv` | Gini vs minimum clone size threshold |
| `sim6_summary.csv` | Adaptation index bias/correlation per f_dying |
| `sim6_replicate_details.csv` | Per-replicate adaptation index metrics |
| `sim7_summary.csv` | AUROC, AUPRC, sensitivity, specificity per scenario |
| `sim7_permutation_test.csv` | Permutation p-value for CYFER vs naive AUROC difference |
| `sim7_null_calibration_fpr.csv` | Observed false positive rate under the null |
| `sim7_replicate_details.csv` | Per-(scenario, replicate) full metrics for all three methods |

## Common Simulation Patterns

**Simulation parameters by script**:

| Script | n_cells | n_clones | n_features | n_causal |
|--------|---------|----------|------------|----------|
| sim1 | 1500 | 60 | 30 | 5 |
| sim2 | 1200 | 100 | 30 | 5 |
| sim3 | 1200 | 50 | 25 | 5 |
| sim4 | varies (K × n_per_clone) | varies | 20 | 4 |
| sim5 | 2500 (100×25) | 100 | 20 | 4 |
| sim6 | 1200 (60×20) | 20 | 20 | — |
| sim7 | 1200 | 60 | 50 | 10 |

Note: sim2 header comment says "n=2000 cells" but the code sets `n_cells <- 1200`.

**Hierarchical data generation** (used in sim1, sim3, sim4, sim5, sim7):
```r
library(MASS)
cc <- mvrnorm(n_clones, rep(0, n_features), sigma_between^2 * diag(n_features))
X  <- t(sapply(seq_len(n_cells), function(i) {
  cc[clone_ids[i], ] + mvrnorm(1, rep(0, n_features), sigma_within^2 * diag(n_features))
}))
X <- scale(X)
true_Z <- as.numeric(X %*% true_beta)
cell_future <- rpois(n_cells, exp(true_Z))
lfc <- tapply(cell_future, clone_labels, sum)
```

**sim3 heritability parameterization**:
```r
# h2 = sigma_between^2 / (sigma_between^2 + sigma_within^2)
# Fix sigma_between = 1, solve for sigma_within:
sigma_within <- sqrt((1 - h2) / h2)
# h2_feature: {0.05, 0.20, 0.50, 0.80, 0.95}
# h2_fate:    {0.10, 0.30, 0.60, 0.90}
```

**Gini coefficient** (the paper uses `dineq::gini.wtd`; identical on non-negative vectors):
```r
gini_coef <- function(x) {
  x <- pmax(x, 0); x <- sort(x); n <- length(x)
  if (n == 0 || sum(x) == 0) return(0)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}
```

**Jaccard index**:
```r
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
```

**AUROC** (Wilcoxon statistic, used in sim2 and sim7):
```r
auroc <- function(scores, labels) {
  n_pos <- sum(labels == 1); n_neg <- sum(labels == 0)
  if (n_pos == 0 || n_neg == 0) return(0.5)
  wilcox.test(scores[labels == 1], scores[labels == 0],
              alternative = "greater")$statistic / (n_pos * n_neg)
}
```

**sim7 p-values** are from `cor.test()` (CYFER) or `wilcox.test()` (naive), then scored
as `-log10(p)` for AUROC/AUPRC ranking. Includes three methods: CYFER, naive (Wilcoxon
by clone fate quartile), and clone-mean correlation.

**sim6 adaptation index**:
```r
adaptation_index <- function(d_vec, weights) {
  weights <- pmax(weights, 0)
  if (sum(weights) < 1e-10) return(mean(d_vec))
  sum(weights * d_vec) / sum(weights)
}
```
Three variants compared: CYFER-weighted, naive (uniform), oracle (true fate potential).

**Within-clone heterogeneity (AED)**, paper Methods: `AED_ℓ = Dist_ℓ / Dist_random`, the mean pairwise Euclidean distance among a clone's cells in the embedding at one time point, over the mean pairwise distance among all cells. The name says "squared" but the written formula does not square; reviewer 2 flagged it, and the Writeup21 memo asks which version to keep.

## Key Reviewer Critiques (Nature Genetics review)

1. **Rev 1a**: Test robustness to wrong growth model (linear, logistic, power-law)
2. **Rev 1b**: Test rare vs. common resistance fractions
3. **Rev 1c/2/3**: Simulate realistic phenotype heritability across lineages
4. **Rev 1d/4/5/9**: Power analysis — minimum clones, cells/clone, capture rate
5. **Rev 4/5/6**: Barcode overlap across timepoints; sampling bias effect on Gini
6. **Rev 8**: Dying cells may inflate adaptation index — show fate-weighting fixes this
7. **Rev 7**: Provide formal p-values and sensitivity/specificity for feature identification

Full reviewer quotes are in `kevin/Writeup_Simulations/README.md` and in `PAPERS_DIR/simulation-planning.txt`.

## Conventions
- New R or code files drafted by Claude carry a `_claude` suffix (`analysis_claude.R`) so a human reviews them before integrating.
- `additional_context/summary.md` is the index of reference material; update it whenever a PDF is added (see `/context-synthesis`).
- Do not run the simulations unless the user explicitly asks — they take 5–40 min each (sim1–sim7) and the Writeup21 sweeps will take hours.
- The `multiomeFate` package must be installed from `CYFER_PKG` before running any sim script, and the installed version checked against the source `DESCRIPTION` (the `run-cyfer` skill has the one-liner).
- Large outputs (RDS, RData, CoSPAR exports and caches) go under `OUT_ROOT`, never in the repo; `.gitignore` excludes `cospar_cache/` and `cospar_out/` as a backstop.
- Git: the large-file hook in `.githooks/` blocks commits ≥ 50 MB; each clone runs `git config core.hooksPath .githooks` once.
- Use "time point" (two words) in prose.

## Post-Prompt Update Instructions
After completing each user prompt, run `/project-state`. It will:
- **Refresh in place** the current-state sections of `CLAUDE_[name].md` (Project Status, Key Methodological Details, Open Questions / Next Steps).
- **Append a dated entry at the bottom** of `HISTORY_[name].md` recording new decisions, resolved/open questions, non-obvious code or LaTeX rationale, and empirical findings.

Do NOT record: things already in the LaTeX/code, git history, or reproducible from code.
