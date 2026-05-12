# CLAUDE.md — multiomeFate Analysis Project

## Project Overview

This is the analysis repository for the paper:
**"Resolution of Selection Versus Adaptation in Cellular Evolution"**
(Chen, Lin, Schaff, et al.) — submitted to Nature Genetics, currently under revision.

The paper introduces **CYFER** (Cell Fate via Exponential Regression), a method to estimate
per-cell fate potential from single-cell multiome data paired with lineage barcoding.

## Repository Layout

```
multiomeFate_analysis/
  kevin/
    Writeup_Simulations/    ← current focus: revision simulation scripts
    Writeup*/               ← prior analysis writeup directories (numbered)
    analysis_pipeline/
  csv/kevin/Writeup_Simulations/   ← flat CSVs exported by make_csvs.R
```

The **multiomeFate R package** lives at a sibling path:
```
/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate/
```

Reviewer critique notes and the paper PDF are at:
```
/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/papers/
```

## CYFER: Core Method

**Model**: Per-cell fate potential `Z_i = β_0 + X_i^T β`
Clone count at t2: `Y_ℓ ~ Poisson(Σ_{i∈ℓ} exp(Z_i))`

**Loss (Poisson-log)**:
```
L(β) = (1/|clones|)[Σ_i exp(x_i^T β) − Σ_ℓ Y_ℓ log(Σ_{i∈ℓ} exp(x_i^T β))] + λ||β||²
```

**Optimization**: BFGS with analytical gradient, multiple random initializations.
**Regularization**: K-fold CV over a decreasing λ-sequence.

**Key outputs**:
- `fate_potential`: log10(expected future progeny) per cell
- `selection_index`: Gini coefficient of estimated fate potentials
- `adaptation_index`: weighted-average expression distance t1→t2, weighted by fate potential

## Package API

```r
library(multiomeFate)

# Step 1: cross-validated fit
fit_res <- cyfer(
  cell_features        = X,          # matrix: cells × features
  cell_lineage         = clone_vec,  # character vector of clone IDs
  lineage_future_count = lfc,        # named integer vector: clone → t2 count
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

# Step 3: compute per-cell fate potential
Z_hat <- as.numeric(X %*% final_fit$coefficient_vec[-1]) +
         final_fit$coefficient_vec[1]
# coefficient_vec[1] is the intercept; [-1] are the feature coefficients
```

**Preprocessing notes**:
- Filter to clones with `lineage_future_count > 0` before fitting
- Filter to clones with `>= 2` cells for stable CV
- `num_folds = min(3, length(unique(clone_vec)) - 1)`
- Scale `X` with `scale()` before fitting

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

Seven scripts address specific reviewer critiques for the Nature Genetics revision.
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

**Run**: `Rscript sim1_growth_modes.R` (single-core; see README.txt for runtimes ~5-40 min each)

**Parallelization**: wrap replicate loops with `parallel::mclapply()` to speed up.

### RDS output structure (per script)

RDS files are saved to `/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup_Simulations/`.

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

**Gini coefficient**:
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

## Key Reviewer Critiques (Nature Genetics Revision)

1. **Rev 1a**: Test robustness to wrong growth model (linear, logistic, power-law)
2. **Rev 1b**: Test rare vs. common resistance fractions
3. **Rev 1c/2/3**: Simulate realistic phenotype heritability across lineages
4. **Rev 1d/4/5/9**: Power analysis — minimum clones, cells/clone, capture rate
5. **Rev 4/5/6**: Barcode overlap across timepoints; sampling bias effect on Gini
6. **Rev 8**: Dying cells may inflate adaptation index — show fate-weighting fixes this
7. **Rev 7**: Provide formal p-values and sensitivity/specificity for feature identification

Full reviewer quotes are in `README.txt` in the same directory.

## Working Notes

- The paper PDF and reviewer notes are in the Dropbox `papers/` folder (not Downloads — the sandbox cannot access `/Users/kevinlin/Downloads/`).
- The git branch for this work is `kevin`.
- Do not run the simulations unless the user explicitly asks — they take 5–40 min each.
- The `multiomeFate` package must be installed from source before running any sim script:
  ```r
  devtools::install_local("/path/to/git/multiomeFate")
  ```
