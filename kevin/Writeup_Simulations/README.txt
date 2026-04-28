================================================================================
SIMULATION STUDIES: CYFER ROBUSTNESS AND VALIDATION
Revision simulations for "Resolution of Selection Versus Adaptation in
Cellular Evolution" (Chen, Lin, Schaff, et al.)
Created: 2026-04-27
================================================================================

This directory contains seven new simulation scripts, each addressing one or
more specific reviewer critiques from the Nature Genetics revision. Each script
is self-contained and can be run independently after installing the multiomeFate
R package (plus MASS as a dependency for multivariate normal generation).

All simulations use synthetic data to allow ground-truth comparison, but are
designed to reflect the statistical structure of the real data (hierarchical
expression, Poisson clone sizes, exponential growth).

================================================================================
SCRIPT OVERVIEW AND REVIEWER CORRESPONDENCE
================================================================================

------------------------------------------------------------------------
sim1_growth_modes.R
------------------------------------------------------------------------
ADDRESSES: Reviewer Comment 1a
REVIEWER QUOTE: "model (a) different growth modes, linear, exponential,
others, and also fit the 'wrong' growth model to each scenario to show
the fits work when the model is inaccurate."

GOAL: Assess whether CYFER (which assumes exponential/Poisson-log growth)
is robust when the true biological growth follows a different model.

SCENARIOS TESTED:
  (A) Exponential: count_i ~ Poisson(exp(Z_i)) — CYFER's exact model
  (B) Linear:      count_i ~ Poisson(Z_i) — linear proportionality
  (C) Logistic:    count_i ~ Poisson(C/(1+exp(-Z_i))) — saturating growth
  (D) Power-law:   count_i ~ Poisson(Z_i^2) — super-linear scaling

DESIGN: 1500 cells, 60 clones, 30 features (5 causal). CYFER always fits
the exponential model regardless of the true data-generating process.
Repeated over 20 Monte Carlo replicates.

METRICS:
  - Pearson correlation between estimated and true cell fate potential
  - Jaccard index of identified causal features
  - Gini coefficient of imputed fate potential

EXPECTED RESULT: CYFER should recover monotone-equivalent rankings under
all three alternative models (B-D), because any smooth monotone
transformation of the linear predictor Z preserves feature ranks. The
logistic model (saturation) is most likely to attenuate performance at
extreme Z values; the power-law model degrades gracefully.

OUTPUT: sim1_growth_modes_results.rds

------------------------------------------------------------------------
sim2_rare_resistance.R
------------------------------------------------------------------------
ADDRESSES: Reviewer Comment 1b; Reviewer Comment 9
REVIEWER QUOTE: "different distributions of phenotype proportions (e.g.
when resistance is very rare vs. moderately common)."

GOAL: Quantify CYFER's ability to identify rare "resistant" progenitor
cells that have disproportionately high fate potential, across a range
of resistance prevalence levels.

RESISTANCE FRACTION (f): 1%, 2%, 5%, 10%, 25%, 50%

DESIGN: 2000 cells, 100 clones, 30 features. Resistant cells have high
expression of 5 causal features (mu=2) and Z=+3. Non-resistant cells
have low expression (mu=0) and Z=-1.5. Cells randomly assigned to clones.
Repeated over 25 Monte Carlo replicates per prevalence level.

METRICS:
  - AUROC for cell-level identification (resistant vs. not)
  - AUROC for clone-level identification (clone contains resistant cell vs. not)
  - Sensitivity and specificity at optimal Youden threshold
  - Jaccard index of causal feature identification

COMPARISON: CYFER vs. naive (assign cells the log-count of their clone)

EXPECTED RESULT: CYFER maintains high AUROC (~0.9+) even at f=1%,
because exponential amplification of rare high-potential cells creates
a strong clone-size signal. The naive method (using clone-level counts
directly) has systematically lower cell-level AUROC, especially at rare
fractions where most cells in large clones are NOT the resistant ones.

OUTPUT: sim2_rare_resistance_results.rds

------------------------------------------------------------------------
sim3_heritability.R
------------------------------------------------------------------------
ADDRESSES: Reviewer Comments 1c, 2, 3
REVIEWER QUOTES:
  "realistically modelling phenotype heritability."
  "We were unsure how realistic the priming versus plasticity simulations
  were... this assignment entirely decouples any heritability that might
  exist in gene expression/accessibility."
  "generating synthetic data to preserve phenotype heritability across
  lineages would be better validation of CYFER."

GOAL: Test CYFER performance across a 2D grid of heritability parameters,
where heritability is defined both at the expression level (do sister
cells share gene expression?) and at the fate level (does expression
predict fate?).

TWO HERITABILITY AXES:
  h2_feature: Fraction of expression variance explained by shared
              clone-level programs (controls inter-clone vs intra-clone
              expression variance). Values: 0.05, 0.20, 0.50, 0.80, 0.95
  h2_fate:    Fraction of fate potential variance explained by expression
              features (controls how predictable fate is from expression).
              Values: 0.10, 0.30, 0.60, 0.90

DESIGN: 1200 cells, 50 clones, 25 features (5 causal). Hierarchical
generative model where cells within a clone share a clone-level expression
center (between-clone variance controlled by h2_feature), and fate
potential = beta^T X + noise (noise controlled by h2_fate). 15 replicates
per combination.

METRICS:
  - Pearson correlation of estimated vs. true fate potential (CYFER)
  - Jaccard index (CYFER vs. naive)

EXPECTED RESULT: Performance increases with both h2_feature and h2_fate.
Notably, CYFER should outperform naive DE even at low h2_feature, because
it explicitly estimates cell-level fate rather than relying on clone
identity. The advantage of CYFER over naive is greatest in low h2_feature
settings, where individual clone identities are poor proxies for fate.

OUTPUT: sim3_heritability_results.rds

------------------------------------------------------------------------
sim4_power_analysis.R
------------------------------------------------------------------------
ADDRESSES: Reviewer Comments 1d, 4, 5, 9
REVIEWER QUOTES:
  "different cell numbers and barcoding efficacies (to give some estimates
  of power)."
  "how stable barcode overlap is across timepoints and how clone size
  affects CYFER's estimation."
  "further simulation to determine the minimum number of clones for
  stable performance."

GOAL: Power analysis to determine the minimum experimental requirements
(clones, cells per clone, barcode capture rate) for reliable CYFER
estimation.

THREE AXES OF EXPERIMENTAL DESIGN:
  Axis 1 — Number of clones (K):
    Fixed n_per_clone=20; K varied: 10, 25, 50, 100, 200, 500
    Question: How many independent lineages does CYFER need?

  Axis 2 — Cells per clone (n_per_clone):
    Fixed K=75; n_per_clone varied: 3, 5, 10, 20, 50
    Question: How many cells per clone are needed for stable estimates?

  Axis 3 — Barcode capture efficacy (p_capture):
    Fixed K=75, n_per_clone=20; p_capture varied: 10%, 25%, 50%, 75%, 100%
    Simulates the fraction of clones successfully detected at t2.
    Question: How much barcode dropout can CYFER tolerate?

DESIGN: 20 features (4 causal); 25 replicates per combination.

METRICS:
  - Pearson correlation of estimated vs. true fate potential
  - Jaccard index of identified features
  - Pearson correlation between true and estimated beta coefficients

EXPECTED RESULT: Performance plateaus around K~50, n_per_clone~10,
p_capture~50%. Below these thresholds, coefficient estimation becomes
unstable. This provides empirical power recommendations for experimental
design.

OUTPUT: sim4_power_analysis_results.rds

------------------------------------------------------------------------
sim5_barcode_dropout.R
------------------------------------------------------------------------
ADDRESSES: Reviewer Comments 4, 5, 6
REVIEWER QUOTES:
  "how stable barcode overlap is across timepoints."
  "what is the fraction of barcodes that overlap between baseline, day 10,
  and week 5?"
  "How does sampling bias (especially in smaller clones) affect the
  estimation of the Fate Potential or Gini Coefficient?"

GOAL: Quantify the effect of two distinct types of incomplete sampling
on CYFER fate potential estimation and Gini coefficient measurement.

THREE PARTS:
  Part A — Clone-level dropout (barcode overlap at t2):
    Starting from complete data, a fraction of clones are undetected at t2.
    Two dropout mechanisms:
      (i)  Random: any clone equally likely to be missed
      (ii) Size-biased: small t2-clones more likely to be missed
    Overlap fractions tested: 10%, 25%, 50%, 75%, 100%

  Part B — Cell-level subsampling at t1:
    Only a fraction of cells are sequenced at the early time point.
    CYFER is fit using only the subsampled cells.
    Cell capture fractions: 10%, 25%, 50%, 75%, 100%

  Part C — Gini coefficient vs. minimum clone size threshold:
    Tests how applying a minimum-clone-size cutoff when computing the
    Gini coefficient of fate potentials affects estimation accuracy.
    Minimum sizes: 1, 5, 10, 15, 20, 30 cells.

DESIGN: 100 clones, 25 cells/clone, 20 features (4 causal); 20 replicates.

METRICS:
  - Pearson correlation of estimated vs. true fate potential
  - True Gini coefficient vs. estimated Gini coefficient
  - Gini estimation error |Gini_hat - Gini_true|

EXPECTED RESULT: The Gini coefficient is sensitive to small-clone bias
(Part C), confirming that the minimum clone size cutoff used in the
paper (>=10 cells) is important for unbiased Gini estimation. Clone-level
dropout (Part A) has moderate impact on fate potential correlation but
larger impact on Gini estimation, especially with size-biased dropout.

OUTPUT: sim5_barcode_dropout_results.rds

------------------------------------------------------------------------
sim6_adaptation_index.R
------------------------------------------------------------------------
ADDRESSES: Reviewer Comment 8
REVIEWER QUOTE: "one would still expect the inverse correlation between
selection and adaptation; consider this counter-example: Suppose that a
population contains 100 cells, with 1% 'primed' for resistance and the
other 99% dying off... at 5 weeks, the population would present a low
'adaptation' signature... Can the authors offer any evidence that this
temporal bias and actively dying cells are not driving the adaptation
signal in the dataset?"

GOAL: Demonstrate that CYFER's fate-potential-weighted adaptation index
correctly downweights dying cells (which may show large expression
changes due to stress), preventing them from artifactually inflating
the adaptation score.

SCENARIO: A clone contains two cell types:
  - Surviving cells (fraction f_survive): high Z, small expression change
    from t1 to t2 (stable molecular state)
  - Dying cells (fraction 1-f_survive): low/negative Z, large expression
    change (stress-induced state shift)

ADAPTATION INDICES COMPARED:
  (A) CYFER-weighted: w_i = exp(estimated fate potential)
  (B) Naive unweighted: w_i = 1 (uniform contribution)
  (C) Oracle-weighted: w_i = exp(true fate potential)

f_dying values tested: 10%, 25%, 50%, 75%, 90%, 99%
20 replicates per setting.

METRICS:
  - Bias = mean(estimated adaptation - true surviving-cell adaptation)
  - Correlation of per-clone adaptation index with true surviving-cell
    adaptation
  - Ratio of naive bias to CYFER bias (inflation factor)

EXPECTED RESULT: Naive adaptation index is strongly biased upward when
dying cells are abundant (f_dying >= 50%), because dying cells dominate
the unweighted average. CYFER-weighted index has near-zero bias across
all f_dying values, because fate-potential weighting suppresses the
contribution of cells with negative Z (expected to die).

OUTPUT: sim6_adaptation_index_results.rds

------------------------------------------------------------------------
sim7_sensitivity_specificity.R
------------------------------------------------------------------------
ADDRESSES: Reviewer Comment 7; Reviewer Comment 1 (general)
REVIEWER QUOTE: "The differential expression analysis in Figure 3E-G is
not clear: why is there no p-value associated with the CYFER analysis?
Additionally, what is the specificity and sensitivity of these approaches
in this DE analysis?"

GOAL: Provide comprehensive sensitivity/specificity analysis for feature
identification in both priming and plasticity scenarios, including formal
statistical tests comparing CYFER vs. the naive approach.

SCENARIOS: Priming (low intra-clone variance) and Plasticity (high
intra-clone variance). 1500 cells, 60 clones, 50 features (10 causal,
40 null). 20 replicates per scenario.

METHODS COMPARED:
  (A) CYFER: correlate each feature with CYFER-estimated fate potential;
      compute p-values using Pearson correlation test.
  (B) Naive: Wilcoxon test comparing cells from top-25% vs. bottom-25%
      clones (by observed t2 count).
  (C) Clone-mean: correlate each feature's clone mean with log(clone count).

METRICS AT alpha=0.05:
  - True Positives, False Positives, True Negatives, False Negatives
  - Sensitivity (recall), Specificity, Precision, F1, Jaccard
  - AUROC (area under ROC curve) for ranking features
  - AUPRC (area under precision-recall curve) for imbalanced setting

STATISTICAL TESTS:
  - Permutation test (100 permutations) comparing CYFER vs. naive AUROC
    to determine if the AUROC difference is statistically significant.
  - Null calibration (all beta=0) to verify false positive rate is ~alpha.

EXPECTED RESULT: In the priming scenario, both methods achieve similar
performance (high sensitivity and specificity). In the plasticity scenario,
CYFER substantially outperforms the naive approach (higher AUROC, Jaccard,
sensitivity) because it correctly attributes high fate potential to the
rare outlier cells within each clone, rather than treating all cells in
the clone as equivalent. The permutation test should confirm p < 0.05 for
the CYFER advantage in the plasticity scenario.

OUTPUT: sim7_sensitivity_specificity_results.rds

================================================================================
RUNNING THE SIMULATIONS
================================================================================

Requirements:
  - R >= 4.1
  - multiomeFate package (install from the git/multiomeFate directory)
  - MASS package (install.packages("MASS"))

To install multiomeFate from source:
  devtools::install_local("/path/to/git/multiomeFate")

Each script can be run from the command line:
  Rscript sim1_growth_modes.R

Or interactively in R/RStudio by sourcing the file.

Expected runtimes (approximate, single core):
  sim1: ~5-10 minutes  (20 replicates x 4 models)
  sim2: ~10-20 minutes (25 replicates x 6 prevalences)
  sim3: ~15-25 minutes (15 replicates x 20 grid cells)
  sim4: ~20-40 minutes (25 replicates x 14 parameter values)
  sim5: ~10-20 minutes (20 replicates x 11 parameter values)
  sim6: ~5-10 minutes  (20 replicates x 6 fractions)
  sim7: ~10-20 minutes (20 replicates x 2 scenarios + permutation tests)

To speed up: reduce n_replicates (e.g., to 5) or lambda_sequence_length
(e.g., to 5) at the cost of noisier estimates.

For parallelization, wrap the replicate loops with parallel::mclapply().

================================================================================
OUTPUT FILE STRUCTURE
================================================================================

Each script saves an .rds file with:
  - Full replicate-level results (for downstream analysis/plotting)
  - Summary table (means and SDs across replicates)
  - Parameter settings

Load results in R:
  res <- readRDS("sim1_growth_modes_results.rds")
  str(res)   # inspect structure

================================================================================
FIGURES TO GENERATE FROM THESE RESULTS
================================================================================

sim1: Line/bar plot of correlation (y) vs growth model (x), with error bars.
      Shows graceful degradation from exponential to alternative models.

sim2: Line plot of AUROC (y) vs resistance fraction f (x), with separate
      lines for CYFER vs. naive. Shows CYFER advantage at low f.

sim3: Heatmap of correlation (or Jaccard) over the h2_feature x h2_fate
      grid, with one panel for CYFER and one for naive, plus a difference panel.

sim4: Three panels (one per axis), each showing correlation (y) vs.
      parameter value (x), with error bars. Dashed line at r=0.6 threshold.

sim5: (A) Line plots of correlation and Gini error vs. p_overlap, one line
          per dropout mechanism.
      (B) Same for cell-level subsampling.
      (C) Bar chart of Gini error vs. minimum clone size threshold.

sim6: Line plot of adaptation index bias (y) vs. fraction of dying cells
      (x), with separate lines for CYFER, naive, and oracle. Should show
      naive bias increasing rapidly while CYFER stays near zero.

sim7: ROC curves (sensitivity vs. 1-specificity) for CYFER and naive, in
      priming and plasticity scenarios (2x2 panel). Bar chart of Jaccard
      indices. Table of permutation p-values.

================================================================================
