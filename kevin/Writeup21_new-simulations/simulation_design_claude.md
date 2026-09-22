# Writeup21: design for the Gini and heterogeneity simulation sweeps

Design memo for the two Figure 4 panels. Nothing here has been run. Decisions are
stated as settled; the points that still need Kevin's input are collected in
Section 9 and cross-referenced as **[Q1]**, **[Q2]**, ... where they arise.

## 1. What the two figures have to show

The Nature Methods mockup (`additional_context/Mockup of figures of Nancy-Sydney
paper.pptx`, slide 4) plans **Figure 4** as two panels, each a ladder of seven
simulated settings:

- **4A: varying clone-size skewness (Gini index).** How well each method recovers the
  genes associated with clonal expansion, as the Gini index of clone sizes rises.
- **4B: varying intraclonal heterogeneity (AED).** The same, as the within-clone
  average Euclidean distance rises.

Methods: CYFER, CoSPAR (Wang et al. 2022), and lineage-level differential expression
(the three columns of the paper's Table 1). Data: fully synthetic, two time points,
RNA only.

**The story both panels must tell is a gradient.** At the easy end of each axis
(low Gini: many clones expanded and many did not; low AED: every cell in a clone
looks alike) all three methods should do well. At the hard end (high Gini: one or two
lineages dominate the later time point; high AED: cells within a clone differ) only
CYFER should still work, because it models the expansion potential of every cell
rather than of clones or of smooth state transitions. The figure's content is where
along each axis the other two methods fall away.

Why the two comparators are expected to fail, in the terms the paper will use:

- **Lineage-DE** has to cut the clones into a "high" and a "low" group and pool all
  high clones together. When one or two clones dominate, the high group is one
  clone, and DE returns that clone's identity genes rather than the expansion
  program. When cells within a clone differ in fate (high AED), clone membership no
  longer predicts a cell's fate, so a clone-level split mislabels cells.
- **CoSPAR** links the two time points through barcodes only, and assumes that
  *within each time point* transcriptomic neighbours share fate. Its transition map
  is `S_t1 · M · S_t2`: `M` is the barcode link, which ties every t1 cell of a clone
  to every t2 cell of that clone equally, and `S_t1`, `S_t2` are within-time
  similarity smoothers (`COSPAR_SRC/cospar/tmap/_tmap_core.py`; the t1-to-t2
  similarity is never read). Three things break it:

  - When cells within a clone differ in fate, the barcode link cannot tell the
    clone-mates apart; the only way CoSPAR separates them is coherence across
    clones at t2, and that requires the t2 neighbourhood structure to mirror the
    t1 fate structure, which need not hold when each clone's descendants drift in
    their own direction.
  - When one or two clones dominate, the "High" fate is one clone's descendants
    and the bias is that clone's identity spread over its neighbours, lineage-DE's
    problem in another form.
  - A clone that went extinct is a single-time clone to CoSPAR, absent from the
    barcode link, whereas to CYFER its zero is data; at high Gini those are most
    clones.

  The assumption that a clone's t2 expression sits near its t1 expression belongs
  to expression-aligning methods (optimal transport, RNA velocity, CoSPAR's
  one-time-clone mode), not to barcoded CoSPAR, and the paper should not attribute
  it to CoSPAR.

The paper's current Methods say why this replaces the priming/plastic pair: those
two semi-synthetic datasets reach Gini 0.26 and 0.28 only, whereas the real data run
from Gini 0.64 before treatment to near 1 at week 5. The Nature Genetics reviewers
asked for exactly this sweep (reviewer 2 comment 3b, reviewer 3 comment 1b) and the
response letter committed to it.

## 2. The generative model

One generator serves both sweeps; each sweep moves one knob and recalibrates the
other so that the second statistic stays fixed (Section 5). Described bottom-up.

### 2.1 Latent cell state

Each earlier-time-point (t1) cell `i` in clone `l` has a latent state
`s_i ∈ R^d` (`d = 10`):

```
m_l   ~ N(0, τ² I_d)                          clone-level centre
s_i   = m_l + e_i,   e_i ~ N(0, σ_w² I_d)     within-clone deviation
```

`τ` and `σ_w` are the two knobs of the whole design. Between-clone spread `τ`
drives clone-size inequality; within-clone spread `σ_w` drives AED. Their ratio is
sim3's heritability (`h² = τ² / (τ² + σ_w²)`), so this is sim3's hierarchical model
with a gene layer on top and the two parameters exposed as the axes the figure
wants. `σ_w` is isotropic: heterogeneity moves the causal and non-causal
coordinates together, which is what real heterogeneity looks like and what AED
measures on all genes.

### 2.2 Fate potential

The first latent coordinate is the **expansion axis**:

```
Z_i = β_0 + β · s_i1
N_i ~ Poisson(exp(Z_i))           progeny of cell i at t2
Y_l = Σ_{i ∈ l} N_i               clone size at t2
```

`β` is fixed across the sweep (`β = 1.5`, as in sim7); `β_0` is solved at every
level so that the expected total number of t2 cells is the same across levels
(Section 2.6). One causal coordinate keeps the truth crisp and the axes
interpretable: `τ` moves the spread of clone means *along the causal axis*, `σ_w`
the within-clone spread along it. Between-clone spread `τ` is the primary knob for
Gini because its Gini is smoothly controllable by bisection; rare "jackpot" cells
within ordinary clones are not folded in, so 4A is not a re-run of sim2's
rare-resistance axis.

### 2.3 Genes

`G = 2000` genes. A gene program matrix `W ∈ R^{G × d}` with a sparse first column:

- 100 **expansion genes** load on `s_1` (50 positive, 50 negative loadings, magnitude
  drawn from `Unif(0.5, 1.0)`); zero loading on `s_1` for all other genes.
- Every gene may load on the other `d − 1` coordinates (dense, small loadings), so
  that PCA of the counts recovers a `d`-dimensional embedding and the non-causal
  coordinates are real structure, not noise.

Counts:

```
log μ_ig = a_g + W_g · s_i
Y_ig ~ NegBin(mean = L_i · μ_ig / Σ_g μ_ig, size = θ)     θ = 10, L_i ~ LogNormal
```

Negative-binomial rather than Poisson so that the "genes are noisy" part of the
problem is honest; library sizes `L_i` around 5,000; baselines `a_g` log-normal so
expression levels span the usual range. Nothing is tuned to make any method look
good, which is the point of synthetic rather than semi-synthetic data; clone
identity is inherited through `m_l`, so the reviewer's "you decoupled heritability"
objection does not apply.

The 100 loading genes are a generator device only. The truth the methods are scored
against is operational (Section 6.1), so no gene has to be declared "true" or
"false".

### 2.4 The later time point

Only CoSPAR needs t2 cells with expression (CYFER needs the clone counts `Y_l`;
lineage-DE splits clones on `Y_l` and tests t1 cells). The construction, step by
step:

1. Draw `s_i` for every t1 cell (Section 2.1) and `N_i` for every t1 cell
   (Section 2.2). `Y_l = Σ_i N_i` is the clone's t2 size, zeros included.
2. The t2 population is the multiset of children: parent `i` contributes `N_i`
   cells. No capture subsampling; the t2 cells CoSPAR sees are exactly the cells
   CYFER counts.
3. Each child of parent `i` gets its own latent state

```
s_child = ρ · s_parent + (1 − ρ) · m_l + sqrt(1 − ρ²) · σ_w · e_child + δ + δ_l
δ_l ~ N(0, τ_δ² I)   on the non-causal coordinates, one draw per clone
```

   with `ρ = 0.8`: a child takes after its specific parent more than after the
   clone's centre, but is not a copy. The t2 shift has two parts, both on the
   non-causal coordinates so that neither creates new expansion genes:
   - `δ`, a **shared** shift (a treatment-response program every surviving cell
     mounts), with its norm set so that t1 and t2 cells barely overlap in the top
     PCs, about three standard deviations of the t1 cloud;
   - `δ_l`, a **clone-specific** shift (each clone's descendants drift in their own
     direction), with spread `τ_δ`.
4. Counts for t2 cells come from the same gene model (Section 2.3).

`ρ`, `δ` and `δ_l` reach the gene calls only through CoSPAR, since CYFER and
lineage-DE never look at t2 expression. Within CoSPAR the two parts of the shift
do different things. The shared `δ` is inert: CoSPAR's t1→t2 link is the barcode
and its smoothing acts within each time point, so moving every t2 cell the same
way leaves the within-t2 neighbourhoods intact. It is kept for realism (the t2
population should not sit on top of the t1 population) and nothing is claimed for
it. The clone-specific `δ_l` is the knob that matters: at `τ_δ = 0` a clone's
descendants sit next to the descendants of transcriptomically similar cells from
other clones, so CoSPAR's cross-clone coherence at t2 can undo the uniform barcode
link and separate clone-mates that differ in fate; as `τ_δ` grows each clone's
descendants become their own t2 cluster, that coherence disappears, and CoSPAR is
left with the clone-level link alone. This is the mechanism by which "the later
time point need not be a smooth continuation of the earlier one" hurts CoSPAR, and
it is separate from the AED axis, which moves `σ_w` at t1.

The main sweeps fix `τ_δ` at a moderate value, comparable to the between-clone
spread `τ` at the middle level, so that clone-specific drift is present but not
extreme. A **supplementary row** moves `τ_δ` over `{0, τ/2, τ, 2τ}` at the middle
level of each axis, holding everything else fixed, to show CoSPAR degrading as t2
structure stops mirroring t1 while CYFER and lineage-DE do not move. **[Q6]**.

### 2.5 Sizes

| Quantity | Value | Why |
|---|---|---|
| clones `L` | 100 | sim2/sim5 scale; ~66 training clones at 3 folds, so up to ~60 CYFER features |
| t1 cells | 1,000: 10 per clone, equal | realistic for a pre-treatment barcoded population; 45 within-clone pairs per clone for AED |
| t2 cells | 3,000 expected total, fixed across levels via `β_0`; realized sizes from 0 up to a few hundred per clone | zeros kept; see **[Q1]** for the ceiling this places on pooled Gini and for the "up to 500" scale |
| genes | 2,000; 100 on the causal axis | enough for a per-gene correlation vector to mean something; small enough that CoSPAR runs in about a minute |
| latent `d` | 10 | PCA dims for CYFER, CoSPAR and AED |
| replicates | 2 per level for the laptop run; 20 later on Hyak | Section 8 |

## 3. Axis 1: Gini of clone sizes

### 3.1 Definition

The **pooled** Gini: `Gini(n_l^{t1} + n_l^{t2})` over the `L` clones, computed with
`gini_coef()` (identical to the paper's `dineq::gini.wtd` on non-negative vectors).
With every t1 clone at 10 cells, the pooled Gini is driven entirely by t2:

```
Gini(10 + Y) = Gini(Y) · mean(Y) / (10 + mean(Y))
```

because adding a constant leaves the mean absolute difference unchanged and raises
the mean. At 3,000 t2 cells over 100 clones, `mean(Y) = 30`, so the pooled Gini is
`0.75 × Gini(Y)` and **cannot exceed about 0.74** even when one clone holds
everything. The t2-only Gini is reported alongside for comparison with the paper's
real-data numbers, which are single-time-point values (0.64 at t1, 0.72–0.81 at day
10, near 1 at week 5). **[Q1]**.

### 3.2 The extreme end is a different regime

At the top level one or a few clones hold most t2 cells. Consequences:

- Lineage-DE's high group is one to three clones, 10–30 t1 cells against ~970. Its
  DE finds those clones' identity genes (their `m_l` on every coordinate), not the
  expansion genes. This is the intended failure.
- CoSPAR's "High" fate is one clone's descendants, and the coherence smoothing
  spreads that clone's bias over its transcriptomic neighbours. Most t1 cells belong
  to clones with no t2 cells at all; in CoSPAR's terms these are single-time
  clones, they contribute no barcode link to the transition map, and their fate bias
  comes entirely from smoothing over neighbours. This is not fatal to CoSPAR as long
  as some clones are observed at both time points (it is, at every level), and it is
  a fair picture of how CoSPAR behaves on the real data. The point to make in the
  paper is that CoSPAR treats an extinct clone as missing, whereas CYFER treats its
  zero as data.
- CYFER sees most clones with `Y_l = 0`. The package accepts zero-count clones and
  the fit keeps them: at high Gini the zeros carry the signal, and the filter to
  `Y_l > 0` used by sim1–sim7 would throw it away. This is a deliberate choice to
  state in the Methods.

Honest expectation for 4A: everyone is fine at the bottom, the comparators degrade
toward the top, and CYFER should degrade slowest if the zero-clone handling is
right. Nothing more is promised before it runs.

### 3.3 Levels

Seven targets on the pooled Gini, `{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}`, the top one
sitting just under the ceiling of Section 3.1. AED is held at its middle level
(Section 4.3) throughout. Whether the ceiling should instead be raised by
generating more t2 cells is **[Q1]**.

## 4. Axis 2: within-clone heterogeneity (AED)

### 4.1 Definition

The paper's clonal variability score, un-squared as written in the Methods
(`paper_nbt.tex`, "Clonal variability score calculation"), with the name left for
the paper to fix:

```
AED_l = Dist_l / Dist_random
Dist_l      = mean pairwise Euclidean distance among clone-l t1 cells in the PCA embedding
Dist_random = mean pairwise Euclidean distance among all t1 cells in the same embedding
```

The embedding is the same 10-PC PCA of log-normalized counts on all 2,000 genes
that the methods use (Section 6), computed on t1 cells only. Heterogeneity is
measured on all genes, causal and non-causal together; there is no separate
"non-causal coordinates only" experiment.

The sweep's x-axis is the **mean of `AED_l` over clones**. Under the generator the
expected squared within-clone distance is `2 d σ_w²` and the expected squared
overall distance is `2 d (τ² + σ_w²)`, so the un-squared AED is approximately
`sqrt(σ_w² / (τ² + σ_w²)) = sqrt(1 − h²)`, which is sim3's feature-heritability axis
in disguise. Two consequences for the range:

- The mean AED is **bounded above by about 1**: it reaches 1 when clones have no
  centre at all (`τ = 0`), because between-clone pairs are never closer, on average,
  than within-clone pairs. Individual clones exceed 1 only when their own spread
  exceeds the typical clone's, which needs `σ_w` to vary across clones. **[Q2]**
  asks how to reconcile this with a 0-to-2 axis.
- The mean AED has a **floor above 0**: at `σ_w = 0` the within-clone distance is
  pure count noise in the PCs, and so is part of the overall distance. The floor is
  found in a pilot run and the lowest level is set just above it.

### 4.2 Holding Gini fixed while AED moves

Raising `σ_w` widens the within-clone spread of `Z`, so `Σ_i exp(Z_i)` gets heavier
tails and the Gini rises on its own. To keep the panels unconfounded, the pooled
Gini is held at 0.4 (the middle of 4A's ladder) throughout 4B: at each AED level,
`τ` is bisected so the realized pooled Gini stays at 0.4 ± 0.03 (Section 5).
Realized Gini and realized AED are reported for every dataset and the metric is
also plotted against realized AED, as a check that the calibration worked.

### 4.3 Levels

Seven targets on the mean un-squared AED, spaced evenly between the pilot floor and
about 0.95, tentatively `{0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}` (in squared terms
`h²` from about 0.9 down to about 0.2, sim3's range). The middle level, AED 0.5
(`h² ≈ 0.75`), is the value held fixed during the Gini sweep. **[Q2]**, **[Q7]**.

## 5. Calibration

Both sweeps are stated in terms of a *realized* statistic (pooled Gini, mean AED),
not a generator parameter, so each level is calibrated by bisection:

1. **Gini sweep.** For a target pooled Gini, hold `σ_w` at the AED-0.5 value,
   simulate at a candidate `τ` (with `β_0` re-solved for 3,000 expected t2 cells),
   average the realized pooled Gini over 10 quick draws, bisect on `τ`. About 10
   generator calls per level; the generator alone runs in under a second, so this is
   negligible.
2. **AED sweep.** For a target mean AED, the ratio `σ_w / τ` sets the AED (Section
   4.1) while `τ` alone sets the Gini at fixed ratio, so bisect on `σ_w / τ` for the
   AED, then on `τ` for the Gini at that ratio, and iterate once. AED needs a PCA of
   the counts per draw, about a second each.

The calibrated `(τ, σ_w, β_0)` per level go into the RDS `calibration` table so the
levels are reproducible, and every replicate reports its realized Gini, AED, total
t2 cells and largest clone. The realized pooled Gini at a fixed `τ` varies across
replicates by about ±0.05 at `L = 100`; with two replicates per level that
variation is visible, which is fine for a pilot.

## 6. The three methods, matched fairly

### 6.1 The metric: a correlation of correlations

No method is asked to threshold genes into "called" and "not called". Instead each
method produces a **per-gene association with clonal expansion**, and the score for
the method is how well that vector agrees with the truth:

```
truth_g  = Spearman( x_g over t1 cells , Z_true )         g = 1..2000
m_g      = the method's per-gene statistic (below)
metric   = Pearson( m , truth ) over the 2,000 genes
```

`x_g` is the log-normalized expression of gene `g` in the 1,000 t1 cells. The truth
is operational: it is the association that a method with perfect knowledge of every
cell's fate potential would see on these very cells, so the ceiling is 1 by
construction and there is no oracle line. The cheap substitute for a ceiling is the
split-half reliability of `truth` (compute it on two random halves of the t1 cells
and correlate the halves), which says how much of the truth is recoverable from
1,000 cells at all. **[Q5]** asks whether the outer correlation should be Pearson
(dominated by the ~100 genes with large `|truth_g|`, which is what recovery of the
expansion program means) or Spearman (which also weighs the ordering among the
~1,900 near-zero genes). Both go in the CSV; one goes in the figure.

Jaccard at BH `q < 0.05` against the operational truth set (`truth` at `q < 0.05`)
is kept as a secondary column for continuity with the mockup's wording, not as the
headline.

### 6.2 Shared embedding

PCA (`d = 10`) on log-normalized t1 counts, computed once per dataset. CYFER uses
it as `cell_features`; CoSPAR gets it as `X_pca` so its similarity graph is built on
the same representation; AED is computed in it. Handing every method the same
embedding is the only way to attribute differences to the model rather than to
preprocessing. (The paper used fastTopics on real data; PCA is the generic choice
for a synthetic count matrix.)

### 6.3 Per-method statistics

**CYFER** (`run-cyfer` skill): fit on t1 PCs with clone counts `Y_l` *including
zeros*; `Z_hat = cell_imputed_score`; `m_g = Spearman(x_g, Z_hat)` over t1 cells.
This is the paper's Methods and sim7's `method_cyfer`. The fitted `β̂` lives on the
10 PCs, not on genes, so it is not the per-gene statistic; **[Q3]** covers the
alternative of projecting `β̂` back through the PCA loadings.

**Lineage-DE**: clones split into "high" and "low" by the **mean** of the t2 clone
sizes `Y_l` (high if `Y_l > mean(Y)`); when a few clones are unusually big the high
group is small and the low group is most clones, which is the intended behaviour
and is why the mean is used rather than the median. Then, using **t1 cells only**,
each gene is compared between all high-clone cells and all low-clone cells. The
per-gene statistic is a signed effect size: `m_g = 2·AUC_g − 1`, the rank-biserial
correlation from the Wilcoxon test (positive when high-clone cells express more).
**[Q4]** covers alternatives.

**CoSPAR** (`run-cospar` skill): `state_info` = `"t1"` for t1 cells; `"High"` for t2
cells of the high clones as defined by lineage-DE's mean split, `"Low"` for the
rest, so the two comparators see identical clone information;
`infer_Tmap_from_multitime_clones(t1 → t2)` on the shared PCA; `fate_bias(High, Low)`
on the intraclone map; `m_g = Spearman(x_g, fate_bias)` over t1 cells. This is the
"matched" route: it isolates the quality of CoSPAR's fate score from any
threshold. CoSPAR's own progenitor→DE recipe is not run. Cells outside the map
carry a bias of 0.5, not `NA`; the correlation is restricted to t1 cells.

**Score-level diagnostics** per dataset: `cor(Z_hat, Z_true)` and
`cor(fate_bias, Z_true)` over t1 cells, so a bad gene-level result can be traced to
the fate score rather than to the gene step.

## 7. Sanity checks before trusting a curve

- Split-half reliability of `truth` at every level (Section 6.1); a level where it
  drops well below 0.9 is one where no method can do well and is flagged as such.
- Label-permutation check at the middle level of each axis: shuffle clone labels
  before fitting and confirm every method's metric is near 0.
- Realized pooled Gini and mean AED within tolerance of their targets in every
  replicate; realized total t2 cells and largest clone size recorded.
- CYFER convergence (the `NULL` returns from the safe wrapper) per level.
- CoSPAR High and Low group sizes per level; if either is empty the bias has
  collapsed and that CoSPAR point is reported as missing rather than as 0.
- Heritability of the PCA embedding (`.anova_percentage`-style) per level, so 4B's
  x-axis can be cross-referenced to sim3.

## 8. Code plan, runtime, progress reporting, outputs

```
kevin/Writeup21_new-simulations/
  func_generate_claude.R      the generator (Sections 2, 5): one function, two knobs
  func_methods_claude.R       PCA / CYFER / CoSPAR export+import / lineage-DE / truth / metrics
  sim_gini_claude.R           axis 1 driver: calibrate, replicate, save RDS
  sim_aed_claude.R            axis 2 driver
  cospar_flat_io.R            copied from .claude/skills/run-cospar/templates
  run_cospar.py               copied from .claude/skills/run-cospar/templates
  make_csvs_claude.R          RDS -> csv/kevin/Writeup21/
  plot_barplots_claude.R      csv -> fig/kevin/Writeup21/
```

Each driver runs the R side in one process and shells out to the `cospar` conda
environment (`COSPAR_ENV`) per dataset via `system2()`. Per dataset: generator under
a second, PCA a second, CYFER about 20 s at 100 clones and 10 features, CoSPAR about
a minute at 4,000 cells, gene statistics a few seconds. Seven levels × 2 replicates
× two axes = 28 datasets, plus the optional `τ_δ` supplementary row (4 values × 2
replicates × two axes = 16 more), plus calibration, is about one to three hours on
the laptop, run sequentially so the progress file is readable.

**Progress reporting.** Each driver appends a time-stamped line to
`OUT_ROOT/Writeup21_new-simulations/progress_gini_claude.txt` (or `_aed_`) at every
stage: calibration of each level (target, calibrated parameters, realized value),
then for each level × replicate the start and end of generation, CYFER, CoSPAR and
the gene step, with elapsed time and a running estimate of time remaining. The same
lines go to the console. A glance at the file says how far along the run is.

**Outputs.** RDS to `OUT_ROOT/Writeup21_new-simulations/` (CoSPAR exports and
caches also stay there, never in the repo). Flat CSVs to `csv/kevin/Writeup21/`,
figures to `fig/kevin/Writeup21/`. RDS schema for both axes: `summary` (one row
per level × method with mean and SD of each metric), `replicate_details` (one row
per level × replicate × method), `calibration` (target, calibrated parameters,
realized mean), `params`.

**Figures.** Grouped barplots in the style of `Writeup17b_barplot-*.R`: one bar per
method at each of the seven levels, height the mean metric of Section 6.1 across
replicates with the individual replicates overplotted as points (an SD over two
replicates is not worth drawing), x-axis labelled with the target statistic and the
realized mean beneath it.

Later, once the pilot looks right: 20 replicates per level, `parallel::mclapply()`
over replicates or the SLURM pattern from `sim3_heritability.slurm` on Hyak, which
will first need a `COSPAR_ENV` built there.

## 9. Questions for Kevin

1. **Pooled Gini ceiling and t2 scale.** With 10 t1 cells per clone and 3,000 t2
   cells in total, the pooled Gini tops out near 0.74 (Section 3.1). The levels are
   set to `0.1–0.7` accordingly. If a higher top is wanted, the total t2 population
   must grow: 10,000 t2 cells give a ceiling near 0.9 at roughly three times the
   CoSPAR cost. Related: "lineages of size 0 through 500" was read as the scale of
   the t2 sizes, not a hard cap; at a fixed total of 3,000 the largest clone at the
   top level will hold well over 500 cells. If 500 is a cap, the total must shrink
   to about 1,000–1,500 and the pooled ceiling falls to about 0.5–0.6. Which is
   preferred: (a) 3,000 cells and levels up to 0.7, (b) a larger total, or (c) a
   hard cap at 500?
2. **AED range.** The mean AED over clones cannot exceed about 1 under any
   generator where clones share one `σ_w` (Section 4.1); the values above 1 in the
   real data are per-clone values for unusually spread clones. Should the x-axis be
   the mean AED on the reachable range (about 0.3 to 0.95, the default), or should
   `σ_w` vary across clones so that the per-clone distribution spans 0 to 2 at the
   top level, with the axis labelled by the mean?
3. **CYFER's per-gene statistic.** The plan uses `Spearman(x_g, Z_hat)` over t1
   cells, the paper's Methods. The alternative is a gene-level coefficient obtained
   by projecting the 10-PC `β̂` back through the PCA loadings (`W_pca · β̂`), which
   is closer to "the β's from CYFER" but is not what the paper does. Keep the
   Spearman version?
4. **Lineage-DE's per-gene statistic.** Now that no threshold is applied, lineage-DE
   needs a signed continuous statistic per gene: rank-biserial correlation from the
   Wilcoxon (the default; scale-free and comparable to a Spearman correlation), log
   fold change of mean log-normalized expression, or a Welch t-statistic?
5. **Outer correlation.** Pearson (default) or Spearman between the method's
   per-gene vector and the truth vector; over all 2,000 genes (default) or over the
   genes with the largest `|truth_g|`?
6. **The clone-specific t2 drift `τ_δ`.** The main sweeps fix `τ_δ` at about the
   between-clone spread `τ` of the middle level, and a supplementary row moves it
   over `{0, τ/2, τ, 2τ}` (Section 2.4). Is a moderate fixed value right for the
   main panels, or should the main panels use `τ_δ = 0` (CoSPAR at its best, so the
   Gini and AED axes alone do the work) with the drift shown only in the
   supplement? And is the supplementary row wanted at all in the pilot?
7. **The fixed values.** Pooled Gini 0.4 held during the AED sweep and mean AED 0.5
   held during the Gini sweep: are these the "middling" values wanted, or should
   they be chosen after the pilot as the last level where all three methods still
   do well?
8. **AED embedding dimension.** AED is computed in the same 10-PC embedding the
   methods use. Is 10 right, or should it match the number of PCs used for the
   real-data AED in Figure 2E?

## 10. What I am uncertain about, stated plainly

- Whether CoSPAR, given the same embedding, does nearly as well as CYFER on 4B. Its
  coherence prior assumes transcriptomic neighbours share fate, which is *true* under
  this generator (fate is a linear function of the embedding). If so, CYFER's edge is
  the extreme-Gini regime and the zero-clone handling rather than heterogeneity as
  such, and the framing should say that.
- How strongly the clone-specific drift `τ_δ` hurts CoSPAR in practice. The
  argument in Section 2.4 is from reading the map construction, not from running
  it; if CoSPAR turns out insensitive to `τ_δ`, then its weakness on these data is
  the uniform barcode link and the extinct clones, and the paper's framing should
  say that rather than "smooth continuum".
- Whether the pooled Gini is the right axis for the paper given its ceiling; the
  t2-only Gini is what the real-data numbers report, and both are recorded so the
  figure can be relabelled without re-running.
- Whether a structural set of 100 loading genes is too easy. Real expansion programs
  have weak effects on many genes; a version with 300 weak-loading genes would be a
  harder and more realistic supplement, and the operational truth handles it without
  any change to the metric.
