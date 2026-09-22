---
name: run-cyfer
description: Use when fitting CYFER (the `multiomeFate` R package) to lineage-barcoded single-cell data — real or simulated — to estimate per-cell fate potential, when scoring cells or ranking genes from a CYFER fit, when computing the selection index (Gini) or adaptation index from one, or when a `cyfer()` / `cyfer_finalize()` call errors on shapes, scales, folds, or identifiability.
---

# Run CYFER

CYFER fits `Y_l ~ Poisson(sum_{i in l} exp(beta_0 + x_i' beta))` per clone `l`, using
features `x_i` measured at the *earlier* time point and clone sizes `Y_l` at the
*later* one. The output is a per-cell **fate potential**, `log10` of expected progeny.
Everything below is what the package's own documentation does not make obvious. The
package is the source of truth for signatures: `?cyfer`, `?cyfer_finalize`, and the
vignette `vignettes/simulation.rmd` in the `CYFER_PKG` location (resolve the path
through your `CLAUDE_[name].md`).

## 1. Load the current package

The installed copy drifts behind the source tree. Compare first:

```r
packageVersion("multiomeFate")                 # installed
read.dcf(file.path(CYFER_PKG, "DESCRIPTION"), fields = "Version")   # source
```

When they differ, reinstall from source before any fit
(`devtools::install_local(CYFER_PKG, force = TRUE)`), or `devtools::load_all(CYFER_PKG)`
for an interactive session. Done when the two versions print the same string.

## 2. Build the three inputs

| Input | Contract | Silent failure it prevents |
|---|---|---|
| `cell_features` | numeric **matrix**, cells x features, **rownames and colnames set**, `scale()`d, **no intercept or constant column** | unscaled features overflow `exp()` and the objective errors; a constant column errors in `cyfer_finalize()` |
| `cell_lineage` | character (factor is coerced), length `nrow(cell_features)`; if named, names must equal `rownames(cell_features)` **in the same order** | a permuted named vector is refused with a "DIFFERENT ORDER" error; an unnamed permuted vector is *not* caught |
| `lineage_future_count` | **named** numeric, names = clone IDs, values = later-time-point cell count | unnamed vector errors; clones absent from `cell_lineage` are dropped silently |

**Effective sample size is the number of clones, not cells.** The unpenalized end of
the lambda path needs `n_training_clones >= n_features + 1`, and `cyfer()` stops
otherwise. With `K` folds, `n_training_clones = n_clones - ceiling(n_clones / K)`. So
with 50 clones and 5 folds, at most 39 features. **Genes are never features**: reduce
to an embedding first (fastTopics or PCA, 10 to 30 dimensions), fit on the embedding,
and implicate genes afterwards (step 5).

Standard preprocessing every caller repeats:

```r
keep_clones <- names(lineage_future_count)[lineage_future_count > 0]
tab <- table(cell_lineage)
keep_clones <- intersect(keep_clones, names(tab)[tab >= 2])
idx <- which(cell_lineage %in% keep_clones)
X   <- scale(cell_features[idx, , drop = FALSE])
cl  <- cell_lineage[idx]
lfc <- lineage_future_count[keep_clones]
num_folds <- max(2, min(3, length(keep_clones) - 1))   # small sims; the paper used 20 on real data
```

Done when `nrow(X) == length(cl)`, `setequal(unique(cl), names(lfc))`, and
`length(lfc) - ceiling(length(lfc) / num_folds) >= ncol(X) + 1`.

## 3. Fit

```r
set.seed(10)
fit_res   <- multiomeFate::cyfer(
  cell_features = X, cell_lineage = cl, lineage_future_count = lfc,
  lambda_initial = 1, lambda_sequence_length = 10, num_folds = num_folds, verbose = 0)
final_fit <- multiomeFate::cyfer_finalize(
  cell_features = X, cell_lineage = cl, fit_res = fit_res, lineage_future_count = lfc)
```

- `lambda_initial` sets the top of a decreasing path that ends at exactly 0; `NA`
  derives it from the data, and that heuristic usually saturates at 101. Pass a
  number for reproducible comparisons across simulation settings.
- `seed_number` (default 10) governs fold assignment and random restarts. A wrapper
  that varies the replicate seed must also vary `seed_number`, or every replicate
  gets the same folds.
- The chosen lambda is the median held-out objective across folds; `lambda = 0` is a
  legitimate outcome on noiseless data, not a bug.
- Wrap the pair in `tryCatch(..., error = function(e) NULL)` inside replicate loops
  and count `NULL`s as non-converged; the CLAUDE.md pattern `fit_cyfer_safe()` is
  the reference form.

Done when `final_fit` has the four elements `cell_imputed_score`, `coefficient_vec`,
`lambda`, `lineage_imputed_count`, and `plot_trainTest(fit_res)` shows the test curve
minimum away from the path's top end.

## 4. Read the scales correctly

| Element | Scale | Use |
|---|---|---|
| `cell_imputed_score` | **log10**(expected progeny), one per input row | the fate potential the paper reports; the selection index and adaptation index take this |
| `coefficient_vec` | natural log, `Intercept` first | scoring new cells; never interpreted gene by gene |
| `lineage_imputed_count` | counts | calibration against `lineage_future_count` |

Identity that catches scale mistakes: `sum(10^cell_imputed_score[cl == l])` equals
`lineage_imputed_count[l]`. Linear predictor is `log(10) * cell_imputed_score`;
expected progeny is `10^cell_imputed_score`. Applying `exp()` to the score is the
documented mistake.

Scoring cells that were not in the fit (there is no `predict()` method):

```r
Z_new <- as.numeric(X_new %*% final_fit$coefficient_vec[-1]) + final_fit$coefficient_vec[1]
score_new <- log10(exp(Z_new)); names(score_new) <- rownames(X_new)
```

`X_new` must be scaled with the *training* centre and scale
(`attr(X, "scaled:center")`, `attr(X, "scaled:scale")`), and its columns must be in
the training order.

## 5. Downstream quantities

**Gene association.** With `Z_hat <- final_fit$cell_imputed_score`, test each gene's
expression (the original genes, not the embedding) against `Z_hat` with
`stats::cor.test()`; BH-adjust; call the set at `padj < 0.05`. Alternative used in the
paper's simulation plots: split cells at `Z_hat >= 0` into winners and losers and run
`wilcox.test()` per gene. `sim7_sensitivity_specificity.R` is the reference for the
correlation route and for how p-values are scored as `-log10(p)` for AUROC/Jaccard.

**Selection index.** Gini coefficient of `10^cell_imputed_score` over the cells of a
clone, or over all cells; the `gini_coef()` helper in the analysis repo's CLAUDE.md is
the definition. Clones need a minimum size (the paper used 10 cells).

**Adaptation index.** Fate-potential-weighted mean of each earlier-time-point cell's
expression distance to its later-time-point position, weights `10^cell_imputed_score`;
`sim6_adaptation_index.R` is the reference implementation. Neither index is in the
package.

**ANOVA decomposition.** `plot_anova()` reports between-clone / total variance of the
score as the "lineage effect" percentage. It rises mechanically as clones shrink, so
compare across settings only at matched clone counts.

## 6. Reference points

- `kevin/Writeup_Simulations/sim1`–`sim7` in this repo are the working examples of
  synthetic data, the safe-fit wrapper, and every metric above.
- Bundled data: `data("priming_simulation")` and `data("plastic_simulation")` in the
  package (7940 cells, 30 fastTopics features, 50 clones) run end to end in under a
  minute and are the fastest smoke test of an install.
- `CYFER_PKG/CLAUDE.md` holds implementation notes: fold construction holds out whole
  clones, `Intercept` is prepended internally, and the source files are named
  `lineage_cv*.R` rather than `cyfer*.R`.
