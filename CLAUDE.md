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

**Defensive wrapper** (used in all sim scripts):
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

## Simulation Scripts (Writeup_Simulations/)

Seven scripts address specific reviewer critiques for the Nature Genetics revision.
All use synthetic data + `multiomeFate` + `MASS` packages.

| Script | Addresses | Key Variable | n_replicates |
|--------|-----------|--------------|--------------|
| `sim1_growth_modes.R` | Rev 1a | growth model (exp/linear/logistic/power) | 20 |
| `sim2_rare_resistance.R` | Rev 1b, 9 | resistance fraction f ∈ {1%,2%,5%,10%,25%,50%} | 25 |
| `sim3_heritability.R` | Rev 1c, 2, 3 | h2_feature × h2_fate grid (5×4) | 15 |
| `sim4_power_analysis.R` | Rev 1d, 4, 5, 9 | #clones / cells-per-clone / capture rate | 25 |
| `sim5_barcode_dropout.R` | Rev 4, 5, 6 | p_overlap, p_cell_capture, min clone size | 20 |
| `sim6_adaptation_index.R` | Rev 8 | f_dying ∈ {10%,25%,50%,75%,90%,99%} | 20 |
| `sim7_sensitivity_specificity.R` | Rev 7, 1 | priming vs. plasticity scenario | 20 |

Each script saves an `.rds` file with `all_results` (replicate-level) and `summary` (aggregated).

**Run**: `Rscript sim1_growth_modes.R` (single-core; see README.txt for runtimes ~5-40 min each)

**Parallelization**: wrap replicate loops with `parallel::mclapply()` to speed up.

## Common Simulation Patterns

**Hierarchical data generation** (used in sim3, sim4, sim5):
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

**AUROC** (requires `pROC` or manual implementation):
```r
# Simple manual AUROC from ranked scores
auroc <- function(scores, labels) {
  n1 <- sum(labels == 1); n0 <- sum(labels == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  sum(rank(scores)[labels == 1]) / (n1 * n0) - (n1 + 1) / (2 * n0)
}
```

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
