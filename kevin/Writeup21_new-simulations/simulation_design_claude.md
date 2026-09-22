# Writeup21: design notes for the Gini and heterogeneity simulation sweeps

Claude, 2026-09-22. A thinking-through, not an implementation. Nothing here has
been run. Questions for Kevin are collected in Section 9 and cross-referenced as
**[Q1]**, **[Q2]**, ... where they arise.

## 1. What the two figures have to show

The Nature Methods mockup (`additional_context/Mockup of figures of Nancy-Sydney
paper.pptx`, slide 4) plans **Figure 4** as two panels, each a ladder of about
seven simulated settings:

- **4A: varying clone-size skewness (Gini index).** Overlap of the genes each method
  calls with the true expansion genes, as the Gini index of clone sizes rises.
- **4B: varying intraclonal heterogeneity (AED).** Same metric, as the
  average-of-squared Euclidean distance within clones rises.

Methods: CYFER, CoSPAR (Wang et al. 2022), and lineage-level differential expression
(the three columns of the paper's Table 1). Metric: Jaccard index between the called
gene set and the true set. Data: fully synthetic, two time points, RNA only.

The paper's current Methods admit why this replaces the priming/plastic pair: those
two semi-synthetic datasets reach Gini 0.26 and 0.28 only, because each varies clones
along a single axis, whereas the real data run from Gini 0.64 before treatment to
near 1 at week 5. The Nature Genetics reviewers asked for exactly this kind of
sweep (reviewer 2 comment 3b, reviewer 3 comment 1b) and the response letter
committed to it. So the sweep is not decoration: it is the paper's answer to "when
does CYFER's advantage appear, and when does everything fail".

Two Table 1 claims are what the panels should make visible: CYFER "accounts for
exponential growth and selection" (4A) and "quantifies fate potential driven by
plasticity" (4B).

## 2. The generative model I would build

One generator serves both sweeps; each sweep moves one knob and recalibrates the
other to stay fixed. I describe it bottom-up.

### 2.1 Latent cell state

Each earlier-time-point (t1) cell `i` in clone `l` has a latent state
`s_i ∈ R^d` (`d = 10`):

```
m_l   ~ N(0, τ² I_d)               clone-level centre
s_i   = m_l + e_i,   e_i ~ N(0, σ_w² I_d)    within-clone deviation
```

`τ` and `σ_w` are the two knobs of the whole design. Between-clone spread `τ`
drives clone-size inequality; within-clone spread `σ_w` drives AED. Their ratio is
the heritability of sim3 (`h² = τ² / (τ² + σ_w²)`), so this is sim3's hierarchical
model with a gene layer on top and with the two parameters exposed as the axes the
figure wants rather than as a heritability grid.

### 2.2 Fate potential

The first latent coordinate is the **expansion axis**:

```
Z_i = β_0 + β · s_i1
N_i ~ Poisson(exp(Z_i))           progeny of cell i at t2
Y_l = Σ_{i ∈ l} N_i               clone size at t2 before sampling
```

`β` is fixed across the sweep (`β = 1.5`, as in sim7), `β_0` is solved so that the
expected total progeny matches a target population size. Making only one latent
coordinate causal keeps the truth set crisp (Section 2.4) and keeps the two axes
interpretable: `τ` moves the spread of clone means *along the causal axis*,
`σ_w` moves the within-clone spread along it. **[Q2]** asks whether this is the
mechanism Kevin wants for high Gini; alternatives are in Section 3.2.

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
problem is honest; library sizes `L_i` around 5,000. Baselines `a_g` from a
log-normal so expression levels span the usual range. Nothing here is tuned to make
any method look good, which is the point of synthetic rather than semi-synthetic
data; it is also why the reviewer's "you decoupled heritability" objection does not
apply, since clone identity is inherited through `m_l`.

### 2.4 Truth set

Structural definition: the 100 genes with non-zero loading on `s_1`. I would also
compute the **operational** truth (Spearman correlation of each gene with the true
`Z_i` over t1 cells, BH `q < 0.05`) and check it recovers the structural set at
Jaccard above 0.9 in every setting; where it does not, the effect size is too small
for *any* method and that level should be reported as such rather than as a method
failure. **[Q5]** asks which definition goes in the figure.

### 2.5 The later time point

The three methods need t2 cells for different reasons: CYFER needs only clone
counts; CoSPAR needs t2 cells with expression and a fate label; lineage-DE as Kevin
specified it tests t2 cells. So t2 cells must have expression. I would generate a
fixed number `n_2` of t2 cells by sampling progeny in proportion to `N_i` (this is
the capture step; `n_2 = 3000`), and give each progeny cell a latent state

```
s_child = ρ · s_parent + (1 − ρ) · m_l + sqrt(1 − ρ²) · σ_w · e_child + δ
```

with `ρ` the within-clone inheritance of the parent's deviation and `δ` a
t2-specific shift shared by all t2 cells (a treatment-response program on
non-causal coordinates, so it does not create new expansion genes). `ρ = 1`, `δ = 0`
means progeny are copies of their parents; `ρ = 0` means progeny regress to the
clone centre. **This knob decides how well lineage-DE at t2 works** (Section 4.3),
so the default matters. I would fix `ρ = 0.5` and `δ` a modest shift, and report
`ρ ∈ {0, 0.5, 1}` as a supplementary row. **[Q4]**.

### 2.6 Sizes

| Quantity | Value | Why |
|---|---|---|
| clones `L` | 100 | sim2/sim5 scale; leaves ~66 training clones at 3 folds, so up to ~60 CYFER features |
| t1 cells | 3,000 (30 per clone, equal) | AED needs pairs within clones; 30 gives 435 pairs per clone |
| t2 cells captured | 3,000 | fixed capture keeps CoSPAR's graph size constant across levels |
| genes | 2,000; 100 causal | enough for a BH-controlled gene call to mean something; small enough that CoSPAR runs in ~1 min |
| latent `d` | 10 | PCA dims for both CYFER and CoSPAR |
| replicates | 20 per level | sim1 used 20; SDs on Jaccard at 20 are about 0.03 |

Whether t1 clone sizes should also be unequal is **[Q1]**; the paper's pre-treatment
Gini of 0.64 says real t1 sizes are, and unequal t1 sizes make the "pooled" Gini
Kevin described meaningful.

## 3. Axis 1: Gini of clone sizes

### 3.1 Which Gini

Kevin's description: clone sizes across all time points at once, where a high value
means one clone has expanded far beyond the others at t2. Three candidate
definitions, all computable from the same data:

1. `Gini(n_l^{t1} + n_l^{t2})` over clones, the pooled size. This is the literal
   reading.
2. `Gini(n_l^{t2})`, the later-time-point inequality, which is what the paper reports
   rising to ~1 and what drives every method's difficulty.
3. `Gini` of the concatenated vector `(n_1^{t1}, ..., n_L^{t1}, n_1^{t2}, ..., n_L^{t2})`.

With equal t1 sizes, (1) is a damped version of (2); with unequal t1 sizes they
diverge. I would compute all three per dataset, sweep on (1) as specified, and
label the x-axis with realized values. **[Q1]**.

The paper's Gini is `dineq::gini.wtd`; the sim scripts use a hand-rolled
`gini_coef()`. They agree on non-negative vectors, so reuse `gini_coef()`.

### 3.2 What generates the inequality

The Gini of `Y_l = Σ_i exp(Z_i)` rises with three different things, and they are
not equivalent for the methods:

| Mechanism | Knob | Who it helps |
|---|---|---|
| (a) Clone means spread along the causal axis | `τ` | lineage-DE: high clones are uniformly high |
| (b) Steeper fate response | `β` | everyone, but the gene truth set gets larger effects |
| (c) Rare "jackpot" cells in otherwise ordinary clones | mixture on `e_i1` | CYFER only, in principle |

The realistic story in the paper is a mix of (a) and (c): a few clones are
uniformly primed, and some ordinary clones win through rare cells. A sweep on `τ`
alone (a) makes 4A partly a re-run of sim3's heritability axis; a sweep on (c)
alone makes 4A a re-run of sim2's rare-resistance axis. I would sweep `τ` as the
primary knob because it is the one whose Gini is smoothly controllable, keep `σ_w`
at the mid AED level, and **calibrate `τ` by bisection** to hit each target Gini
level in expectation (Section 5). **[Q2]** is whether Kevin wants (c) folded in.

### 3.3 The extreme end is a different regime

At Gini near 1, one or two clones hold nearly all t2 cells. Consequences:

- Lineage-DE compares one "high" clone with everything else, so it finds that
  clone's *identity* genes (its `m_l` on all coordinates), not the expansion genes.
  Its precision collapses even if recall is fine. This is the intended failure.
- CoSPAR's "High" fate is one clone's descendants, and its coherence smoothing then
  spreads that clone's bias over its transcriptomic neighbours. Also likely to fail,
  for the same reason wearing a different hat.
- CYFER sees most clones with `Y_l = 0`. Those zeros are informative (they say the
  cells are low), but every existing sim script filters to `Y_l > 0` before
  fitting, while the package itself accepts zeros. **At high Gini the filter throws
  away the signal.** I would keep zero-count clones in the fit and treat this as a
  deliberate choice to state in the Methods. **[Q6]**.

So the honest expectation for 4A is: everyone is fine in the middle, everyone
degrades at the top, and the figure's content is *how fast*. CYFER should degrade
slowest if the zero-clone handling is right. I would not promise more than that
before running it.

### 3.4 Levels

Seven targets: Gini ∈ {0.2, 0.35, 0.5, 0.6, 0.7, 0.8, 0.9} on definition (2), which
brackets the real data (0.64 at t1, 0.72–0.81 at day 10, ~1 at week 5). The exact
list depends on what the calibration can reach at `L = 100`; a Gini of 0.95 with
100 clones needs one clone to hold about 90% of cells, which is reachable but leaves
`n_2` almost entirely from one clone. **[Q7]** covers the level list and count.

## 4. Axis 2: within-clone heterogeneity (AED)

### 4.1 Definition, and a naming trap

The paper's Methods define the clonal variability score as

```
AED_l = Dist_l / Dist_random
Dist_l      = mean pairwise Euclidean distance among clone-l cells in the embedding
Dist_random = mean pairwise Euclidean distance among all cells
```

despite the name, the distances are **not squared** in the written formula, and
reviewer 2 (comment 20) already flagged the inconsistency. For the simulation I
would compute both the ratio of mean distances (as written) and the ratio of mean
squared distances, report which one the figure uses, and let the paper's Methods be
fixed to match. **[Q10]**. The embedding is the same PCA the methods use (Section 6),
at t1 only.

Under the generator, `E[squared distance within clone] = 2 d σ_w²` and
`E[squared distance overall] = 2 d (τ² + σ_w²)`, so the squared-distance AED is
`σ_w² / (τ² + σ_w²) = 1 − h²`, which is why this axis is sim3's feature-heritability
axis in disguise. Levels of AED map to `σ_w/τ` directly; the un-squared version is
a monotone transform of the same ratio.

### 4.2 Holding Gini fixed while AED moves

Raising `σ_w` widens the within-clone spread of `Z`, so `Σ_i exp(Z_i)` gets heavier
tails and the Gini rises on its own. If 4B does not hold Gini fixed, the two panels
are confounded. I would fix the target Gini at the middle level of 4A (about 0.6)
and, for each AED level, bisect on `τ` (or on `β_0` and `τ` jointly) so the realized
Gini stays at 0.6 ± 0.03. Report realized Gini and realized AED for every dataset
and plot Jaccard against the realized AED as a check that the calibration worked.
**[Q11]**.

### 4.3 Where the heterogeneity lives

Two very different things raise AED:

- **Heterogeneity on the causal axis** (`σ_w` on `s_1`). This is plasticity in the
  paper's sense: cells in one clone differ in fate potential. Lineage-level methods
  lose here because clone membership stops predicting cell fate.
- **Heterogeneity on non-causal axes** (`σ_w` on `s_2..s_d` only). AED rises, but
  every cell in a clone still has the same fate. Lineage-DE is unaffected; only the
  embedding gets noisier.

The primary sweep should move `σ_w` isotropically (both at once), because that is
what real heterogeneity looks like and it is what AED measures. A one-row
supplement moving only the non-causal coordinates would show that AED *per se* is
not what breaks lineage methods, which pre-empts a reviewer asking "isn't this just
noise". **[Q13]**.

### 4.4 The lineage-DE at t2 problem

Kevin specified lineage-DE as: split clones by expansion between the two time
points, then test t2 cells of high clones against t2 cells of low clones. Under
inheritance (`ρ = 1`), t2 cells are copies of the parents that expanded, so the
high clones' t2 cells are *already selected* for high `s_1`. That DE recovers the
expansion genes well even at high AED, because selection did the work that the
method cannot do at t1. In sim7 and in the paper's Figure 3, the naive DE used **t1
cells** of high versus low clones, and that is the version plasticity breaks.

I would run both: lineage-DE on t1 cells (the published comparator) and on t2 cells
(the one Kevin described). Under `ρ < 1` and a t2 shift `δ`, the t2 version drifts
toward calling the response program rather than the priming program, which is a
realistic failure worth showing. **[Q3]** and **[Q4]** together decide the default.

### 4.5 Levels

Seven targets on squared-distance AED: {0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95},
equivalently `h²` from 0.9 down to 0.05, which is sim3's range. Real-data anchors
would help: the values in Figure 2E and S3-1 should sit inside this range and I do
not have them. **[Q10]**.

## 5. Calibration

Both sweeps are stated in terms of a *realized* statistic (Gini, AED), not a
generator parameter. Two ways to handle that:

1. **Bisection per level.** For a target Gini, simulate at candidate `τ`, compute
   expected Gini over 10 quick draws, bisect. About 10 generator calls per level,
   negligible cost. Then draw the 20 replicates at the calibrated `τ`.
2. **Report realized values only.** Sweep the parameter, plot against the realized
   statistic. Simpler, but the x-axis becomes uneven and the levels are not
   nameable in a caption.

I would do (1) and *also* plot against realized values in a supplement, since the
realized Gini at a fixed `τ` varies across replicates by ±0.05 or so at `L = 100`.
The calibrated parameters go into the RDS `params` so the levels are reproducible.

## 6. The three methods, matched fairly

The Jaccard index depends as much on the gene-calling rule as on the score. To
separate "better fate score" from "different threshold", every method gets the same
final step where possible.

**Shared embedding.** PCA (`d = 10`, or `d` chosen by elbow) on log-normalized t1
counts, computed once per dataset. CYFER uses it as `cell_features`; CoSPAR gets it
as `X_pca` so its similarity graph is built on the same representation. Handing
both methods the same embedding is the only way to attribute differences to the
model rather than to the preprocessing. The paper used fastTopics; PCA is the
generic choice for a synthetic count matrix. **[Q8]**.

**CYFER** (`run-cyfer` skill): fit on t1 PCs, clone counts `Y_l` at t2 *including
zeros* (Section 3.3); per-gene Spearman `cor.test` against `cell_imputed_score` over
t1 cells; BH; call `q < 0.05`. This is exactly the paper's Methods and sim7's
`method_cyfer`.

**CoSPAR** (`run-cospar` skill): `state_info` = `"t1"` for t1 cells; `"High"` for t2
cells of the top-`k` clones by expansion ratio, `"Low"` for the rest;
`infer_Tmap_from_multitime_clones(t1 → t2)`; `fate_bias(High, Low)` on the
intraclone map; then **two gene-calling routes**:

- CoSPAR-native: `progenitor` (bias > 0.6 / < 0.4) → `differential_genes` (Wilcoxon,
  BH), the recipe in the CoSPAR paper. This is what a CoSPAR user would do.
- Matched: per-gene Spearman correlation of expression with `fate_bias` over t1
  cells, BH. This isolates the score's quality from the calling rule.

Report the native route in the main panel and the matched route in a supplement, or
the reverse; either way say which. **[Q9]** covers `k` and the thresholds. The
smoke test in the `run-cospar` skill recovered 10 of 10 causal genes on a toy, so
the pipeline itself is not the risk; the risk is CoSPAR doing *well* on a shared
embedding, which would be an honest result and a reason to lean on the extreme
Gini end where its High fate collapses.

**Lineage-DE**: clones ranked by `n_l^{t2} / n_l^{t1}` (or by `n_l^{t2}` when t1
sizes are equal); "high" = top quartile, "low" = bottom quartile (sim7's rule) or
top-`k` versus rest (Writeup14's rule); per-gene Wilcoxon between the two cell
groups; BH; `q < 0.05`. Run on t2 cells (Kevin's spec) and t1 cells (sim7's).
**[Q3]**, **[Q9]**.

**Oracle**: per-gene Spearman with the true `Z_i`. Its Jaccard is the ceiling and
should be reported on every panel as a grey line; where the ceiling itself drops,
the level is uninformative.

**Metrics per dataset**: Jaccard at `q < 0.05` (the headline), Jaccard at top-`|truth|`
genes (threshold-free), AUROC and AUPRC of the gene ranking, sensitivity and
specificity at `q < 0.05`, plus `cor(Z_hat, Z_true)` and `cor(fate_bias, Z_true)`
over t1 cells as score-level diagnostics. Sim7's helper functions already implement
all of these.

## 7. Sanity checks before trusting a curve

- Oracle Jaccard above 0.9 at every level; otherwise the truth is not recoverable
  and the level is dropped or flagged.
- Null run (`β = 0`) at the middle level: every method's false-positive rate at
  `q < 0.05` near 0.05 (sim7's calibration check).
- Realized Gini and AED within tolerance of targets for at least 18 of 20
  replicates.
- CYFER convergence count per level (the `NULL` returns from the safe wrapper).
- CoSPAR progenitor group sizes per level; a zero on either side means the bias
  collapsed, and that level's CoSPAR point is reported as missing rather than as
  Jaccard 0.
- Heritability of the PCA embedding (`.anova_percentage`-style) per level, so 4B's
  x-axis can be cross-referenced to sim3.

## 8. Code plan and runtime

```
kevin/Writeup21_new-simulations/
  func_generate_claude.R      the generator (Sections 2, 5): one function, two knobs
  func_methods_claude.R       CYFER / CoSPAR export+import / lineage-DE / oracle / metrics
  sim_gini_claude.R           axis 1 driver: calibrate, replicate, save RDS
  sim_aed_claude.R            axis 2 driver
  cospar_flat_io.R            copied from .claude/skills/run-cospar/templates
  run_cospar.py               copied from .claude/skills/run-cospar/templates
  make_csvs_claude.R          RDS -> csv/kevin/Writeup21_new-simulations/
```

Each driver runs the R side in one process and shells out to the `cospar` conda
environment per dataset via `system2()`. Per dataset: generator under a second,
CYFER 20–60 s at 100 clones and 10 features, CoSPAR about 1–2 min at 6,000 cells,
gene tests a few seconds. Seven levels × 20 replicates × two axes ≈ 280 datasets ≈
8–12 hours single-core, so `parallel::mclapply()` over replicates or the SLURM
pattern from `sim3_heritability.slurm`. **[Q7]**. RDS outputs go to the
`SIM_OUT` location; CoSPAR exports and caches are large and regenerable and stay
under `SIM_OUT`, never in the repo.

RDS schema (both axes): `summary` (one row per level per method with mean/SD of each
metric), `replicate_details` (one row per level × replicate × method), `calibration`
(target, calibrated parameter, realized mean), `params`.

## 9. Questions for Kevin

1. **Gini definition.** Pooled `n^{t1} + n^{t2}` per clone (literal reading),
   t2-only, or the concatenated vector? And should t1 clone sizes be unequal too
   (real pre-treatment Gini is 0.64)?
2. **Mechanism for high Gini.** Sweep between-clone spread `τ` on the causal axis
   (my default), or fold in rare jackpot cells within clones? The two make 4A a
   heritability story or a rare-resistance story respectively.
3. **Lineage-DE cells.** You specified t2 cells. Run the t1 version too (it is the
   comparator in the current Figure 3 and in sim7)? Which one goes in the main panel?
4. **What t2 cells look like.** How much of the parent's state does a progeny cell
   keep (`ρ`), and is there a t2-wide response shift `δ`? This decides whether t2
   lineage-DE is nearly an oracle (`ρ = 1`, `δ = 0`) or realistically confounded.
5. **Truth set.** Structural (100 genes with non-zero loading) or operational
   (correlation with true `Z`)? I would use structural and verify the two agree.
6. **Zero-count clones in CYFER.** Keep them (my recommendation; the package accepts
   them and they carry the signal at high Gini) or filter to `Y_l > 0` as every
   existing sim script does?
7. **Scale.** Seven levels per axis as in the mockup? 20 replicates? Laptop with
   `mclapply` or Hyak SLURM? Total budget is roughly 10 CPU-hours as sized.
8. **Embedding.** PCA on log-normalized counts, shared by CYFER and CoSPAR (my
   default), or fastTopics to match the paper's real-data analysis?
9. **Gene-calling operating point.** BH `q < 0.05` for every method as the headline,
   with top-`k` and AUROC as supplements? For CoSPAR, native progenitor-DE route or
   the correlation route matched to CYFER in the main panel? Top-`k` or quartile
   split for defining high/low clones, and what `k`?
10. **Real-data anchors.** What AED values does the real data span (Figure 2E, S3-1),
    and is the squared or un-squared version the one the paper will keep? The
    sweep should bracket the real range on both axes.
11. **Decoupling.** Hold Gini fixed at ~0.6 while sweeping AED, and hold AED at its
    middle level while sweeping Gini, by recalibrating? Or let the second statistic
    float and just report it?
12. **CoSPAR fate labels.** Define `High` by the same clone split lineage-DE uses
    (so the two comparators see identical information), and pass CoSPAR the same PCA?
13. **Control row for 4B.** Include the "heterogeneity on non-causal axes only"
    variant as a supplement?
14. **Outputs.** CSVs under `csv/kevin/Writeup21_new-simulations/` like
    `Writeup_Simulations`, and figure style from `Writeup17b_barplot-*.R`
    (grouped bars by method) or line-plus-ribbon over the seven levels?

## 10. What I am uncertain about, stated plainly

- Whether CoSPAR, given the same embedding, does nearly as well as CYFER on 4B. Its
  coherence prior assumes transcriptomic neighbours share fate, which is *true* under
  this generator (fate is a linear function of the embedding). If so, CYFER's edge
  on 4B is the extreme-Gini regime and the zero-clone handling, not heterogeneity as
  such, and the framing should say that.
- Whether a Gini near 0.9 is reachable with 100 clones without `n_2` becoming one
  clone; the calibration will tell, and more clones (200) may be needed at the top.
- Whether a structural truth set of 100 genes is too easy. Real expansion programs
  have weak effects on many genes; a version with 300 weak-loading genes would be a
  harder and more realistic supplement.
