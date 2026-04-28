# Simulation Studies: CYFER Robustness and Validation

Revision simulations for *"Resolution of Selection Versus Adaptation in Cellular Evolution"*
(Chen, Lin, Schaff, et al.) — *Nature Genetics* revision.

Seven simulation scripts, each addressing one or more specific reviewer critiques. Every script is
self-contained and can be run independently after installing the `multiomeFate` R package (plus
`MASS` as a dependency for multivariate normal generation).

All simulations use synthetic data to allow ground-truth comparison, but are designed to reflect
the statistical structure of the real data (hierarchical expression, Poisson clone sizes,
exponential growth).

---

## Script Overview

| Script | Addresses | Key variable | Replicates |
|--------|-----------|--------------|-----------|
| `sim1_growth_modes.R` | Rev 1a | growth model (exp/linear/logistic/power) | 20 |
| `sim2_rare_resistance.R` | Rev 1b, 9 | resistance fraction *f* ∈ {1%, 2%, 5%, 10%, 25%, 50%} | 25 |
| `sim3_heritability.R` | Rev 1c, 2, 3 | *h²\_feature* × *h²\_fate* grid (5×4) | 15 |
| `sim4_power_analysis.R` | Rev 1d, 4, 5, 9 | #clones / cells-per-clone / capture rate | 25 |
| `sim5_barcode_dropout.R` | Rev 4, 5, 6 | *p\_overlap*, *p\_cell\_capture*, min clone size | 20 |
| `sim6_adaptation_index.R` | Rev 8 | *f\_dying* ∈ {10%, 25%, 50%, 75%, 90%, 99%} | 20 |
| `sim7_sensitivity_specificity.R` | Rev 7, 1 | priming vs. plasticity scenario | 20 |

---

## Script Details

### `sim1_growth_modes.R`

**Addresses:** Reviewer Comment 1a

> *"model (a) different growth modes, linear, exponential, others, and also fit the 'wrong'
> growth model to each scenario to show the fits work when the model is inaccurate."*

**Goal:** Assess whether CYFER (which assumes exponential/Poisson-log growth) is robust when the
true biological growth follows a different model.

**Scenarios tested:**

| Label | Data-generating process |
|-------|------------------------|
| (A) Exponential | `count_i ~ Poisson(exp(Z_i))` — CYFER's exact model |
| (B) Linear | `count_i ~ Poisson(Z_i)` |
| (C) Logistic | `count_i ~ Poisson(C / (1 + exp(-Z_i)))` — saturating growth |
| (D) Power-law | `count_i ~ Poisson(Z_i²)` — super-linear scaling |

**Design:** 1500 cells, 60 clones, 30 features (5 causal). CYFER always fits the exponential
model regardless of the true data-generating process. 20 Monte Carlo replicates.

**Metrics:** Pearson correlation (estimated vs. true fate potential), Jaccard index (causal
features), Gini coefficient of imputed fate potential.

**Expected result:** CYFER should recover monotone-equivalent rankings under all three alternative
models (B–D). The logistic model is most likely to attenuate performance at extreme *Z* values;
the power-law model degrades gracefully.

**Output:** `sim1_growth_modes_results.rds`

---

### `sim2_rare_resistance.R`

**Addresses:** Reviewer Comments 1b, 9

> *"different distributions of phenotype proportions (e.g. when resistance is very rare vs.
> moderately common)."*

**Goal:** Quantify CYFER's ability to identify rare "resistant" progenitor cells with
disproportionately high fate potential, across a range of resistance prevalence levels.

**Resistance fractions (*f*):** 1%, 2%, 5%, 10%, 25%, 50%

**Design:** 2000 cells, 100 clones, 30 features. Resistant cells: high expression of 5 causal
features (μ=2) and Z=+3. Non-resistant cells: low expression (μ=0) and Z=−1.5. 25 replicates
per prevalence level.

**Metrics:** AUROC (cell-level and clone-level), sensitivity/specificity at optimal Youden
threshold, Jaccard index (causal features).

**Comparison:** CYFER vs. naive (assign cells the log-count of their clone).

**Expected result:** CYFER maintains high AUROC (~0.9+) even at *f*=1%, because exponential
amplification of rare high-potential cells creates a strong clone-size signal.

**Output:** `sim2_rare_resistance_results.rds`

---

### `sim3_heritability.R`

**Addresses:** Reviewer Comments 1c, 2, 3

> *"realistically modelling phenotype heritability."*
> *"generating synthetic data to preserve phenotype heritability across lineages would be better
> validation of CYFER."*

**Goal:** Test CYFER performance across a 2D grid of heritability parameters.

**Heritability axes:**

| Parameter | Description | Values |
|-----------|-------------|--------|
| *h²\_feature* | Fraction of expression variance explained by shared clone-level programs | 0.05, 0.20, 0.50, 0.80, 0.95 |
| *h²\_fate* | Fraction of fate potential variance explained by expression features | 0.10, 0.30, 0.60, 0.90 |

**Design:** 1200 cells, 50 clones, 25 features (5 causal). Hierarchical generative model. 15
replicates per grid cell.

**Metrics:** Pearson correlation (estimated vs. true fate potential), Jaccard index (CYFER vs. naive).

**Expected result:** Performance increases with both *h²\_feature* and *h²\_fate*. CYFER advantage
over naive DE is greatest at low *h²\_feature*, where individual clone identities are poor proxies
for fate.

**Output:** `sim3_heritability_results.rds`

---

### `sim4_power_analysis.R`

**Addresses:** Reviewer Comments 1d, 4, 5, 9

> *"different cell numbers and barcoding efficacies (to give some estimates of power)."*
> *"further simulation to determine the minimum number of clones for stable performance."*

**Goal:** Determine the minimum experimental requirements (clones, cells per clone, barcode
capture rate) for reliable CYFER estimation.

**Three experimental design axes:**

| Axis | Fixed | Varied | Values |
|------|-------|--------|--------|
| Number of clones *K* | `n_per_clone = 20` | *K* | 10, 25, 50, 100, 200, 500 |
| Cells per clone | `K = 75` | `n_per_clone` | 3, 5, 10, 20, 50 |
| Barcode capture rate | `K = 75`, `n_per_clone = 20` | `p_capture` | 10%, 25%, 50%, 75%, 100% |

**Design:** 20 features (4 causal); 25 replicates per combination.

**Metrics:** Pearson correlation (fate potential), Jaccard index (features), Pearson correlation
(beta coefficients).

**Expected result:** Performance plateaus around K~50, n\_per\_clone~10, p\_capture~50%. This
provides empirical power recommendations for experimental design.

**Output:** `sim4_power_analysis_results.rds`

---

### `sim5_barcode_dropout.R`

**Addresses:** Reviewer Comments 4, 5, 6

> *"how stable barcode overlap is across timepoints."*
> *"How does sampling bias (especially in smaller clones) affect the estimation of the Fate
> Potential or Gini Coefficient?"*

**Goal:** Quantify the effect of incomplete sampling on fate potential estimation and Gini
coefficient measurement.

**Three parts:**

| Part | Description | Values tested |
|------|-------------|--------------|
| A — Clone-level dropout | Fraction of clones undetected at t2 (random vs. size-biased) | 10%, 25%, 50%, 75%, 100% overlap |
| B — Cell-level subsampling at t1 | Fraction of cells sequenced at early time point | 10%, 25%, 50%, 75%, 100% |
| C — Gini vs. min clone size | Minimum clone size cutoff when computing Gini of fate potentials | 1, 5, 10, 15, 20, 30 cells |

**Design:** 100 clones, 25 cells/clone, 20 features (4 causal); 20 replicates.

**Metrics:** Pearson correlation (fate potential), true vs. estimated Gini, Gini estimation error.

**Expected result:** The Gini coefficient is sensitive to small-clone bias (Part C), confirming
that the ≥10 cells cutoff used in the paper is important for unbiased estimation. Size-biased
dropout has larger impact on Gini than random dropout.

**Output:** `sim5_barcode_dropout_results.rds`

---

### `sim6_adaptation_index.R`

**Addresses:** Reviewer Comment 8

> *"one would still expect the inverse correlation between selection and adaptation; consider this
> counter-example: Suppose that a population contains 100 cells, with 1% 'primed' for resistance
> and the other 99% dying off... Can the authors offer any evidence that this temporal bias and
> actively dying cells are not driving the adaptation signal in the dataset?"*

**Goal:** Demonstrate that CYFER's fate-potential-weighted adaptation index correctly downweights
dying cells, preventing them from artifactually inflating the adaptation score.

**Scenario:** Each clone contains surviving cells (high Z, small expression change) and dying
cells (low/negative Z, large stress-induced expression change).

**Adaptation indices compared:**

| Method | Weights |
|--------|---------|
| CYFER-weighted | `w_i = exp(estimated fate potential)` |
| Naive unweighted | `w_i = 1` |
| Oracle-weighted | `w_i = exp(true fate potential)` |

***f\_dying* values tested:** 10%, 25%, 50%, 75%, 90%, 99% — 20 replicates per setting.

**Metrics:** Bias (estimated − true surviving-cell adaptation), per-clone adaptation correlation,
naive/CYFER bias ratio.

**Expected result:** Naive adaptation index is strongly biased upward when dying cells are
abundant (*f\_dying* ≥ 50%). CYFER-weighted index has near-zero bias across all *f\_dying*
values.

**Output:** `sim6_adaptation_index_results.rds`

---

### `sim7_sensitivity_specificity.R`

**Addresses:** Reviewer Comments 7, 1 (general)

> *"The differential expression analysis in Figure 3E–G is not clear: why is there no p-value
> associated with the CYFER analysis? Additionally, what is the specificity and sensitivity of
> these approaches in this DE analysis?"*

**Goal:** Provide comprehensive sensitivity/specificity analysis for feature identification in
priming and plasticity scenarios, with formal statistical tests.

**Scenarios:** Priming (low intra-clone variance) and Plasticity (high intra-clone variance).
1500 cells, 60 clones, 50 features (10 causal, 40 null). 20 replicates per scenario.

**Methods compared:**

| Method | Approach |
|--------|----------|
| CYFER | Pearson correlation of each feature with CYFER-estimated fate potential; *p*-values via correlation test |
| Naive | Wilcoxon test: cells from top-25% vs. bottom-25% clones (by observed t2 count) |
| Clone-mean | Correlate each feature's clone mean with log(clone count) |

**Metrics at α=0.05:** TP/FP/TN/FN, sensitivity, specificity, precision, F1, Jaccard, AUROC, AUPRC.

**Statistical tests:** Permutation test (100 permutations) comparing CYFER vs. naive AUROC;
null calibration (all β=0) to verify FPR ~α.

**Expected result:** In the plasticity scenario, CYFER substantially outperforms the naive
approach (higher AUROC, Jaccard, sensitivity) because it correctly attributes high fate potential
to rare outlier cells within each clone. Permutation test should confirm *p* < 0.05 for the
CYFER advantage in plasticity.

**Output:** `sim7_sensitivity_specificity_results.rds`

---

## Running the Simulations

### Requirements

- R ≥ 4.1
- `multiomeFate` package (install from source)
- `MASS` package

```r
# Install multiomeFate from source
devtools::install_local("/path/to/git/multiomeFate")

# Install MASS
install.packages("MASS")
```

### Execution

```bash
# From the command line
Rscript sim1_growth_modes.R
```

Or source interactively in R/RStudio.

### Expected runtimes (single core, approximate)

| Script | Runtime |
|--------|---------|
| `sim1` | 5–10 min (20 replicates × 4 models) |
| `sim2` | 10–20 min (25 replicates × 6 prevalences) |
| `sim3` | 15–25 min (15 replicates × 20 grid cells) |
| `sim4` | 20–40 min (25 replicates × 14 parameter values) |
| `sim5` | 10–20 min (20 replicates × 11 parameter values) |
| `sim6` | 5–10 min (20 replicates × 6 fractions) |
| `sim7` | 10–20 min (20 replicates × 2 scenarios + permutation tests) |

To speed up: reduce `n_replicates` (e.g., to 5) or `lambda_sequence_length` (e.g., to 5) at the
cost of noisier estimates. For parallelization, wrap replicate loops with
`parallel::mclapply()`.

---

## Output File Structure

Each script saves an `.rds` file containing:

- Full replicate-level results (for downstream analysis/plotting)
- Summary table (means and SDs across replicates)
- Parameter settings

```r
res <- readRDS("sim1_growth_modes_results.rds")
str(res)   # inspect structure
```

---

## Figures to Generate

| Script | Figure description |
|--------|-------------------|
| `sim1` | Line/bar plot of correlation (y) vs. growth model (x), with error bars |
| `sim2` | Line plot of AUROC (y) vs. resistance fraction *f* (x), CYFER vs. naive |
| `sim3` | Heatmap of correlation over *h²\_feature* × *h²\_fate* grid — CYFER, naive, and difference panels |
| `sim4` | Three panels (one per axis): correlation (y) vs. parameter value (x), with error bars; dashed line at *r*=0.6 |
| `sim5` | (A) Correlation and Gini error vs. *p\_overlap*, one line per dropout mechanism; (B) same for cell subsampling; (C) Gini error vs. min clone size |
| `sim6` | Line plot of adaptation index bias (y) vs. fraction dying (x) — CYFER, naive, oracle; naive bias increases rapidly while CYFER stays near zero |
| `sim7` | ROC curves (2×2 panel: priming/plasticity × CYFER/naive); Jaccard bar chart; permutation *p*-value table |
