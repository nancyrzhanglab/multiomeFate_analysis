# Writeup21: design for the Gini and heterogeneity simulation sweeps

Design memo for the two Figure 4 panels. Nothing here has been fitted; the only
code run is the generation-only demo of Section 4.1. Decisions are stated as
settled. The points that still need Kevin's input are collected in Section 9
and cross-referenced as **[Q1]**, **[Q2]** where they arise; neither blocks
implementation, and the memo states the working default for each.

## 1. What the two figures have to show

The Nature Methods mockup (`additional_context/Mockup of figures of Nancy-Sydney
paper.pptx`, slide 4) plans **Figure 4** as two panels, each a ladder of seven
simulated settings:

- **4A: varying clone-size skewness (Gini index).** How well each method recovers the
  genes associated with clonal expansion, as the Gini index of clone sizes at the
  later time point rises.
- **4B: varying intraclonal heterogeneity (AED).** The same, as the within-clone
  average squared Euclidean distance rises.

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

One generator serves both sweeps; each sweep moves one knob and holds the other so
that the second statistic stays fixed (Section 5). Described bottom-up.

### 2.1 Latent cell state

Each earlier-time-point (t1) cell `i` in clone `l` has a latent state
`s_i ∈ R^d` (`d = 10`):

```
m_l   ~ N(0, τ² I_d)                                   clone-level centre
s_i   = m_l + e_i,   e_i ~ N(0, diag(σ_l1², ..., σ_ld²))   within-clone deviation
σ_lk² = σ_w² · g_l   for the clone-varying coordinates k,  σ_lk² = σ_w² otherwise
g_l   ~ Gamma(shape = κ, rate = κ)                      one draw per clone, mean 1
```

The within-clone spread **varies across clones**: `σ_w` is the typical spread and
`g_l` (mean 1, shape `κ = 1.5`, fixed across all levels of both sweeps) makes
some clones tight and some diffuse, which is what puts the per-clone AED on a
0-to-2 range (Section 4.1). `κ = ∞` gives every clone the same spread and is
kept in the generator as a reference setting only.

**Where the clone-varying factor `g_l` acts** is a generator switch,
`spread_variation`:

- `"isotropic"`: `g_l` multiplies every coordinate, causal and non-causal alike.
  Heterogeneity then moves the expansion axis together with everything else,
  and a diffuse clone has a wider spread of fate potential and so a larger
  *expected* clone size. This is Jensen's inequality for the log-normal mean:
  `Z_i` is normal within a clone with variance `β² σ_l1²`, and

  ```
  E[Y_l | m_l, σ_l1] = 10 · exp(β_0 + β m_l1 + β² σ_l1² / 2)
  ```

  so two clones with the same centre but different spreads differ in expected
  size by the factor `exp(β² σ_l1² / 2)`: the jackpot cells on the far side of
  the expansion axis outweigh the cells on the near side. Called the
  **spread-size coupling** below; the demo measures it (Section 4.1).
- `"noncausal"`: `g_l` multiplies the `d − 1` non-causal coordinates only, and
  the causal coordinate keeps the shared spread `σ_w` in every clone. A clone's
  spread on the expansion axis, and with it its expected t2 size, is then the
  same for tight and diffuse clones, so **AED and clone size are independent by
  construction**. The per-clone AED still spans 0 to 2 (Section 4.1), because
  nine of the ten coordinates carry the clone-to-clone variation.

The working default is `"noncausal"`; whether to keep it is **[Q1]**, decided on
the demo results in Section 4.1. Either way the *mean* within-clone spread on
the causal coordinate is `σ_w` and grows along the AED axis, so the axis moves
within-clone fate heterogeneity as intended; the switch only decides whether the
clone-to-clone *variation* in spread reaches the fate potential.

The two knobs are reparameterized as a **scale** and a **share**:

```
τ   = s · sqrt(h²)          s   overall latent scale
σ_w = s · sqrt(1 − h²)      h²  between-clone share of variance = τ² / (τ² + σ_w²)
```

`h²` is sim3's heritability. The AED depends on `h²` alone; the Gini depends on
`s` as well, because `s` sets the spread of the fate potential across cells
(Section 2.2). Section 5 shows that this makes each axis a one-knob calibration.

### 2.2 Fate potential

The first latent coordinate is the **expansion axis**:

```
Z_i = β_0 + β · s_i1
N_i ~ Poisson(exp(Z_i))           progeny of cell i at t2
Y_l = Σ_{i ∈ l} N_i               clone size at t2
```

`β` is fixed across the sweep (`β = 1.5`, as in sim7); `β_0` is solved at every
level so that the expected total number of t2 cells is 3,000 given the drawn
latent states (`β_0 = log(3000 / Σ_i exp(β s_i1))`). One causal coordinate keeps
the truth crisp and the axes interpretable. Because `Z` is linear in `s_i1`, the
spread of `exp(Z)` across cells, and with it the inequality of the `Y_l`, grows
with the overall latent scale `s`: that is what the Gini sweep moves. Rare
"jackpot" cells within ordinary clones are not folded in as a separate device, so
4A is not a re-run of sim2's rare-resistance axis.

### 2.3 Genes

`G = 2000` genes. A gene program matrix `W ∈ R^{G × d}` with a sparse first column:

- 100 **expansion genes** load on `s_1` (50 positive, 50 negative loadings, magnitude
  drawn from `Unif(0.5, 1.0)`); zero loading on `s_1` for all other genes.
- Every gene may load on the other `d − 1` coordinates (dense loadings,
  `N(0, 0.25²)`), so that PCA of the counts recovers a `d`-dimensional embedding
  and the non-causal coordinates are real structure, not noise.

Counts:

```
log μ_ig = a_g + W_g · s_i
Y_ig ~ NegBin(mean = L_i · μ_ig / Σ_g μ_ig, size = θ)     θ = 10, L_i ~ LogNormal
```

Negative-binomial rather than Poisson so that the "genes are noisy" part of the
problem is honest; library sizes `L_i` around 5,000 (`sdlog = 0.3`); baselines
`a_g ~ N(0, 1)` on the log scale so expression levels span the usual range.
Nothing is tuned to make any method look good, which is the point of synthetic
rather than semi-synthetic data; clone identity is inherited through `m_l`, so the
reviewer's "you decoupled heritability" objection does not apply.

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
s_child = ρ · s_parent + (1 − ρ) · m_l + sqrt(1 − ρ²) · σ_l ∘ e_child + δ + δ_l
δ_l ~ N(0, τ_δ² I)   on the non-causal coordinates, one draw per clone
```

   with `ρ = 0.8`: a child takes after its specific parent more than after the
   clone's centre, but is not a copy (`σ_l ∘ e_child` is the coordinate-wise
   product with the clone's spread vector of Section 2.1). The t2 shift has two
   parts, both on the non-causal coordinates so that neither creates new
   expansion genes:
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
it is separate from the AED axis, which moves `σ_l` at t1.

**The main sweeps fix `τ_δ = 0.6`** in the trailblazing round, the same absolute
value at every level of both sweeps. **[Q2]**. The value is the between-clone
spread `τ` at the provisional middle point of the two sweeps (t2 Gini 0.5 at mean
AED 0.4, Section 8): a latent-only calibration at `h² = 0.6`, `κ = 1.5` gives
`s ≈ 0.81`, hence `τ = s · sqrt(0.6) ≈ 0.63` and `σ_w ≈ 0.51`, under the
`"noncausal"` default (`s ≈ 0.76`, `τ ≈ 0.59` under `"isotropic"`). So `δ_l`
scatters the clones' t2 centres by about as much as their t1 centres already
differ, and by a little more than one within-clone standard deviation per
coordinate: clone-specific drift is present but does not swamp the t1 structure.
A fixed absolute value rather than "the level's own `τ`" keeps the t2 drift
identical across levels, so a change along either axis is attributable to that
axis alone. If CoSPAR does far worse than expected even at the easy end of either
axis, `τ_δ` is the first thing to lower (Section 10). A **supplementary row**
moves `τ_δ` over `{0, 0.3, 0.6, 1.2}` (`0`, `τ_δ/2`, `τ_δ`, `2τ_δ`) at the
chosen middle level of each axis, holding everything else fixed, to show CoSPAR
degrading as t2 structure stops mirroring t1 while CYFER and lineage-DE do not
move; it runs in the final round (Section 8), once the middle levels are known,
and its grid is re-centred on the final `τ_δ` if that value changes.

### 2.5 Sizes

| Quantity | Value | Why |
|---|---|---|
| clones `L` | 100 | sim2/sim5 scale; ~66 training clones at 3 folds, so up to ~60 CYFER features |
| t1 cells | 1,000: 10 per clone, equal | realistic for a pre-treatment barcoded population; 45 within-clone pairs per clone for AED |
| t2 cells | 3,000 expected total, fixed across levels via `β_0`; realized clone sizes from 0 upward with **no cap** (over 800 in one clone at the top Gini level) | zeros kept: at high Gini they carry the signal |
| genes | 2,000; 100 on the causal axis | enough for a per-gene correlation vector to mean something; small enough that CoSPAR runs in about a minute |
| latent `d` | 10 | PCA dims for CYFER, CoSPAR and AED, throughout |
| replicates | 2 per level in the trailblazing round; 20 in the final round | Section 8 |

## 3. Axis 1: Gini of clone sizes

### 3.1 Definition

The **t2-only** Gini: `Gini(Y_l)` over the `L = 100` clones, zeros included,
computed with `gini_coef()` (identical to the paper's `dineq::gini.wtd` on
non-negative vectors). This is the statistic the paper's real-data numbers report
(single-time-point values: 0.64 at t1, 0.72–0.81 at day 10, near 1 at week 5), so
the simulated axis and the real anchors are on the same scale.

Its range is essentially the whole unit interval. The floor is Poisson counting
noise: with `Z` nearly constant across cells every clone is `Poisson(30)` and the
Gini is about 0.10 (the demo of Section 4.1 measured 0.12–0.16 at its smallest
latent scale). The ceiling is one clone holding everything, Gini `1 − 1/L = 0.99`.
The pooled Gini `Gini(10 + Y_l)` over t1 and t2 cells together, which at these
sizes is `0.75 × Gini(Y)` and tops out near 0.74, is recorded as a secondary
column so the figure could be relabelled without re-running, but it is not the
axis.

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

Seven targets on the t2-only Gini for the trailblazing round,
`{0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9}`: the bottom just above the Poisson floor,
the top in the one-clone-dominates regime, denser in the range the real data
occupy. The final round may re-space them around where the comparators fall away.
AED is held at its fixed value (Section 8) throughout.

## 4. Axis 2: within-clone heterogeneity (AED)

### 4.1 Definition

The paper's clonal variability score, **squared** as its name says (the Methods
formula in `paper_nbt.tex`, "Clonal variability score calculation", currently
reads un-squared and will be corrected to match):

```
AED_l = SqDist_l / SqDist_random
SqDist_l      = mean pairwise squared Euclidean distance among clone-l t1 cells in the PCA embedding
SqDist_random = mean pairwise squared Euclidean distance among all t1 cells in the same embedding
```

The embedding is the same 10-PC PCA of log-normalized counts on all 2,000 genes
that the methods use (Section 6), computed on t1 cells only. Heterogeneity is
measured on all genes, causal and non-causal together; there is no separate
"non-causal coordinates only" experiment. The squared form is computed without
forming pairs: the mean pairwise squared distance among `n` points is twice the
sum of the per-coordinate sample variances.

The sweep's x-axis is the **mean of `AED_l` over clones**. Under the generator the
expected within-clone squared distance for clone `l` is `2 Σ_k σ_lk²` and the
overall one is about `2 d (τ² + σ_w²)`, so

```
AED_l  ≈ (1 − h²) · g_l                       (isotropic)
AED_l  ≈ (1 − h²) · (1 + (d − 1) g_l) / d     (noncausal)
mean_l AED_l ≈ 1 − h²                          (either)
```

The mean is sim3's feature-heritability axis exactly, and is **bounded by 1**
(reached when clones have no centre at all, `τ = 0`), while individual clones
exceed 1 in proportion to their own `g_l`. That is why the per-clone spread
needs clone-varying `σ_l`, and why the axis is labelled by the mean while the
per-clone distribution is what spans 0 to 2.

**Demo** (`demo_aed_range_claude.R`, generation and PCA only, no fitting; results
in `csv/kevin/Writeup21/demo_aed_range_claude.csv`). Per-clone squared AED at
latent scale `s = 1`, `κ = 1.5`, one seed, under both placements of `g_l`. The
second seed agrees to within the sampling noise of the maximum, which is set by
the largest `g_l` draw. The last column is the spread-size coupling, the
Spearman correlation between a clone's AED and its t2 size over the 100 clones
(its sampling standard error is about 0.1).

`spread_variation = "isotropic"`:

| `h²` | mean AED | min | median | 95th pct | max | t2 Gini | Spearman(AED_l, Y_l) |
|---|---|---|---|---|---|---|---|
| 0.90 | 0.12 | 0.01 | 0.10 | 0.27 | 0.38 | 0.56 | −0.07 |
| 0.75 | 0.28 | 0.02 | 0.23 | 0.62 | 0.83 | 0.55 | −0.02 |
| 0.50 | 0.52 | 0.03 | 0.43 | 1.14 | 1.33 | 0.55 | 0.20 |
| 0.25 | 0.74 | 0.04 | 0.65 | 1.60 | 1.67 | 0.61 | 0.35 |
| 0.10 | 0.87 | 0.04 | 0.77 | 1.74 | 1.86 | 0.65 | 0.50 |
| 0.00 | 0.95 | 0.05 | 0.86 | 1.84 | 1.99 | 0.71 | 0.64 |

`spread_variation = "noncausal"`:

| `h²` | mean AED | min | median | 95th pct | max | t2 Gini | Spearman(AED_l, Y_l) |
|---|---|---|---|---|---|---|---|
| 0.90 | 0.12 | 0.02 | 0.09 | 0.27 | 0.36 | 0.57 | −0.14 |
| 0.75 | 0.28 | 0.03 | 0.23 | 0.62 | 0.80 | 0.55 | −0.11 |
| 0.50 | 0.52 | 0.06 | 0.44 | 1.14 | 1.31 | 0.50 | −0.10 |
| 0.25 | 0.74 | 0.04 | 0.65 | 1.58 | 1.68 | 0.46 | −0.02 |
| 0.10 | 0.87 | 0.05 | 0.78 | 1.74 | 1.86 | 0.43 | −0.07 |
| 0.00 | 0.96 | 0.06 | 0.85 | 1.88 | 1.99 | 0.39 | −0.02 |

What the demo established:

- **The mean tracks `1 − h²` to within a few hundredths at every level, under
  either placement**, so the count-noise floor is negligible with 2,000 genes in
  10 PCs (mean 0.12 at `h² = 0.9`, against 0.10 predicted). The level list can be
  set from `h²` directly.
- **With a shared spread (`κ = ∞`) the per-clone AED stays within 0.6–1.4 even at
  the top**; with `κ = 1.5` it spans about 0.05 to 2.0 at the top level (up to
  2.5 on the other seed) and about 0.01 to 0.4 at the bottom. `κ = 1` (an
  exponential) reaches 2.3–2.7 at the top but crowds more clones near zero (5th
  percentile 0.07–0.10 against 0.15–0.20 at `κ = 1.5`); `κ = 3` tops out near
  2.1–2.2. So `κ = 1.5`, fixed across levels, gives the 0-to-2 range wanted.
- **The per-clone AED range is the same under both placements.** At the top
  level the `"noncausal"` placement spans 0.06–1.99 (seed 10) and 0.11–2.50
  (seed 20) against 0.05–1.99 and 0.12–2.47 for `"isotropic"`; the medians and
  5th–95th percentile ranges agree to within 0.05 at every level. Moving `g_l`
  off the causal coordinate costs nothing on the axis that matters.
- **The spread-size coupling is real under `"isotropic"` and gone under
  `"noncausal"`.** Isotropic: diffuse clones are bigger at t2, with the Spearman
  between `AED_l` and `Y_l` about 0 at the bottom and 0.5–0.65 at the top.
  Non-causal: −0.02 to −0.21 across levels and seeds, within the sampling noise
  of a Spearman over 100 clones. A latent-only check over 50 draws (no counts, no
  PCA) puts the coupling at 0.60 (SD 0.08) isotropic against 0.07 (SD 0.12)
  non-causal at the top level, and 0.26 against 0.01 at `h² = 0.5`; the
  residual 0.07 is the sample-level version of the same Jensen effect (a clone
  whose 10 cells happen to scatter more on the causal axis expands a little
  more) and is present under any placement.
- **The two placements move the Gini in opposite directions along the AED axis
  at fixed `s`.** Isotropic: 0.55 at the bottom rising to 0.71 at the top,
  because the coupling fattens the tail. Non-causal: 0.57 falling to 0.39, the
  same as with a shared spread (0.58 to 0.39), because a clone's size sums 10
  cells whose causal coordinates become independent draws as `h² → 0`, which
  averages the inequality out. Either way the Gini is a side effect and the
  calibration of Section 5 removes it; under `"noncausal"` the calibrated `s`
  rises along the AED axis (from about 0.68 to 1.25 for a Gini of 0.5, against
  0.67 to 0.82 under `"isotropic"`, from the latent-only calibration).
- **The latent scale `s` is the Gini's knob under either placement.** At the
  top AED level with `κ = 1.5` and `"noncausal"`, the t2 Gini goes from 0.36–0.39
  at `s = 1` to 0.22–0.25 at `s = 0.7`, 0.17–0.19 at `s = 0.5` and 0.12–0.14 at
  `s = 0.3`, while the mean AED stays at 1.0 and the per-clone range at 0.2 to 3.
  So any Gini above the Poisson floor is reachable at any AED level by moving
  `s`, and the AED does not move with it.

For reference, the un-squared mean AED at the same settings runs 0.33 → 0.95 over
the same `h²` range: the squared version has the wider dynamic range at the easy
end, which is one more reason to prefer it.

### 4.2 Holding Gini fixed while AED moves

The t2 Gini moves on its own along the AED axis (Section 4.1: down from 0.57 to
0.39 under the `"noncausal"` default, up from 0.55 to 0.71 under `"isotropic"`).
To keep the panels unconfounded, the t2 Gini is held at its fixed value
throughout 4B: at each AED level, `s` is bisected so the realized Gini stays
within ±0.03 of the target (Section 5). Realized Gini and realized mean AED are
reported for every dataset and the metric is also plotted against realized AED,
as a check that the calibration worked.

### 4.3 Levels

Seven targets on the mean squared AED for the trailblazing round, evenly spaced
from just above the floor to the top, `{0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0}`,
which is `h² ∈ {0.9, 0.75, 0.6, 0.45, 0.3, 0.15, 0}`. At the top level the
per-clone AED spans about 0.05 to 2; the barplot's x-axis carries the mean, and
the per-clone 5th–95th percentile range is printed beneath it so a reader sees
the 0-to-2 spread. The final round may re-space them.

## 5. Calibration

Both sweeps are stated in terms of a *realized* statistic (t2 Gini, mean squared
AED), not a generator parameter. The demo showed the two statistics separate
cleanly onto the two knobs, so no iteration is needed:

1. **AED is set directly.** Mean squared AED ≈ `1 − h²` regardless of `s`, so each
   AED level sets `h² = 1 − target`; the realized mean is verified per dataset and
   the level re-centred only if it is off by more than 0.03.
2. **Gini is bisected on `s`.** At the level's `h²`, simulate at a candidate `s`
   (with `β_0` re-solved for 3,000 expected t2 cells), average the realized t2
   Gini over 10 quick draws of the latent states and clone sizes (no counts
   needed), and bisect on `s` until the mean is within 0.01 of the target. About
   10 generator calls per level; the latent generator alone runs in milliseconds.

For the **Gini sweep**, `h²` is fixed at the fixed AED value and `s` is bisected
for each Gini target. For the **AED sweep**, `h²` moves and `s` is bisected at
each level so the Gini stays at its fixed value. At the provisional middle point
(Gini 0.5, `h² = 0.6`) the latent-only calibration gives `s ≈ 0.81` under the
`"noncausal"` default, which is where the `τ_δ` of Section 2.4 comes from.

The calibrated `(h², s, τ, σ_w, β_0)` per level go into the RDS `calibration`
table so the levels are reproducible, and every replicate reports its realized
Gini, mean AED, per-clone AED quantiles, total t2 cells, largest clone and number
of extinct clones. The realized t2 Gini at a fixed `s` varies across replicates by
about ±0.05 at `L = 100` (visible in the demo's two seeds, more so at the top);
with two replicates per level that variation is visible, which is fine for the
trailblazing round.

## 6. The three methods, matched fairly

### 6.1 The metric: a correlation of correlations

No method is asked to threshold genes into "called" and "not called". Instead each
method produces a **per-gene association with clonal expansion**, and the score for
the method is how well that vector agrees with the truth:

```
truth_g  = Spearman( x_g over t1 cells , Z_true )         g = 1..2000
m_g      = the method's per-gene statistic (below)
metric   = Spearman( m , truth ) over all 2,000 genes
```

`x_g` is the log-normalized expression of gene `g` in the 1,000 t1 cells. The truth
is operational: it is the association that a method with perfect knowledge of every
cell's fate potential would see on these very cells, so the ceiling is 1 by
construction and there is no oracle line. The cheap substitute for a ceiling is the
split-half reliability of `truth` (compute it on two random halves of the t1 cells
and correlate the halves), which says how much of the truth is recoverable from
1,000 cells at all. The outer correlation is Spearman over all genes; the Pearson
version goes in the CSV as a secondary column, since it weights the ~100
large-`|truth_g|` genes more heavily and may separate the methods differently.

Jaccard at BH `q < 0.05` against the operational truth set (`truth` at `q < 0.05`
from `cor.test`) is kept as a secondary column for continuity with the mockup's
wording, not as the headline; each method's `q` comes from the test that produces
its statistic (`cor.test` for CYFER and CoSPAR, the Wilcoxon for lineage-DE).

### 6.2 Shared embedding

PCA (`d = 10`) on log-normalized t1 counts, computed once per dataset. CYFER uses
it as `cell_features`; CoSPAR gets it as `X_pca` so its similarity graph is built on
the same representation; AED is computed in it. Ten PCs throughout, including the
AED. Handing every method the same embedding is the only way to attribute
differences to the model rather than to preprocessing. (The paper used fastTopics
on real data; PCA is the generic choice for a synthetic count matrix.)

### 6.3 Per-method statistics

**CYFER** (`run-cyfer` skill): fit on t1 PCs with clone counts `Y_l` *including
zeros*; `Z_hat = cell_imputed_score`; `m_g = Spearman(x_g, Z_hat)` over t1 cells.
This is the paper's Methods and sim7's `method_cyfer`. The fitted `β̂` lives on the
10 PCs, not on genes, and is not projected back to genes: correlation among genes
would let the projected coefficient land on the wrong genes, whereas the Spearman
version depends only on the fate potential being well estimated.

**Lineage-DE**: clones split into "high" and "low" by the **mean** of the t2 clone
sizes `Y_l` (high if `Y_l > mean(Y)`); when a few clones are unusually big the high
group is small and the low group is most clones, which is the intended behaviour
and is why the mean is used rather than the median. Then, using **t1 cells only**,
each gene is compared between all high-clone cells (`n_H` of them) and all
low-clone cells (`n_L`). The per-gene statistic is the **rank-biserial correlation
from the Wilcoxon rank-sum test**, written out explicitly:

```
x_H = x_g over the n_H high-clone t1 cells;  x_L = x_g over the n_L low-clone t1 cells
U_g = number of (i in H, j in L) pairs with x_H[i] > x_L[j], counting ties as 1/2
    = stats::wilcox.test(x_H, x_L)$statistic        # R's "W" is exactly this U for the first argument
m_g = 2 · U_g / (n_H · n_L) − 1                    # in [−1, 1]; equals 2·AUC − 1
```

`m_g` is positive when high-clone cells express gene `g` more, zero when the two
groups are exchangeable, and scale-free, so it is comparable to the Spearman
statistics of the other two methods. In code it is computed for all 2,000 genes
at once from ranks (`U_g = Σ_{i ∈ H} rank(x_g)[i] − n_H (n_H + 1) / 2`, with
average ranks for ties) and cross-checked against `stats::wilcox.test` on a
handful of genes; the Wilcoxon `p`-value from the same call feeds the secondary
Jaccard column.

*Why `m_g` rather than `U_g`, given that the test's null distribution and
`p`-value are built on `U_g`.* Within one dataset the two are the same ranking:
`m_g` is an increasing affine function of `U_g` with `n_H` and `n_L` fixed
across genes, so `Spearman(m, truth) = Spearman(U, truth)` exactly, and the
headline metric does not depend on the choice. The `p`-value is a function of
`U_g` and the two group sizes, so nothing is lost either: the same
`wilcox.test` call supplies both, `m_g` as the effect size and `p` for the
Jaccard column. The reasons to record `m_g` rather than `U_g` are the ones a
reader of the CSV will care about:

- **Comparability across levels.** `U_g` lives on `[0, n_H · n_L]`, and `n_H`
  changes with the level (at the top of the Gini axis the high group is one to
  three clones, 10–30 cells; at the bottom it is a third of the clones), so
  `U_g` from two levels are on different scales. `m_g` is on `[−1, 1]` at every
  level, which is what the Pearson secondary column and any pooled plot need.
- **A sign and a zero.** `m_g` is centred at 0 under exchangeability and its
  sign says which group is higher, the same convention as the Spearman
  statistics of CYFER and CoSPAR and of the truth vector; `U_g` is centred at
  `n_H n_L / 2`, and a two-sided `p`-value (or `−log10 p`, as sim7 used) folds
  the sign away, so it could only be correlated with `|truth_g|`, a different
  and weaker target.
- **Effect size, not evidence.** `m_g = 2·AUC − 1` (Cliff's delta) answers "how
  separated are the two groups", which is the counterpart of the Spearman
  statistics; `p` answers "how sure are we", which conflates the separation with
  `n_H · n_L` and so would reward the levels with balanced groups for reasons
  unrelated to gene recovery. For a fixed `(n_H, n_L)` the two-sided `p` is a
  decreasing function of `|m_g|` up to the tie correction in the variance, so a
  signed `−log10 p` would rank the genes the same way within a dataset as
  `m_g`; `m_g` is simply the version that is also interpretable across
  datasets.

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
- Realized t2 Gini and mean squared AED within tolerance of their targets in every
  replicate; realized total t2 cells, largest clone, number of extinct clones and
  the per-clone AED quantiles recorded.
- Spearman between `AED_l` and `Y_l` per dataset, so the spread-size coupling of
  Section 2.1 is visible in the outputs rather than assumed (near 0 under the
  `"noncausal"` default; 0.5–0.65 at the top AED level if `"isotropic"` is
  chosen instead).
- CYFER convergence (the `NULL` returns from the safe wrapper) per level.
- CoSPAR High and Low group sizes per level; if either is empty the bias has
  collapsed and that CoSPAR point is reported as missing rather than as 0.
- Heritability of the PCA embedding (`.anova_percentage`-style) per level, so 4B's
  x-axis can be cross-referenced to sim3.

## 8. Code plan, two rounds, progress reporting, outputs

```
kevin/Writeup21_new-simulations/
  demo_aed_range_claude.R     DONE: generation-only demo behind Section 4.1
  func_generate_claude.R      the generator (Sections 2, 5): one function, knobs (h², s, κ, spread_variation, τ_δ)
  func_methods_claude.R       PCA / CYFER / CoSPAR export+import / lineage-DE / truth / metrics
  sim_gini_claude.R           axis 1 driver: calibrate, replicate, save RDS
  sim_aed_claude.R            axis 2 driver
  cospar_flat_io.R            copied from .claude/skills/run-cospar/templates
  run_cospar.py               copied from .claude/skills/run-cospar/templates
  make_csvs_claude.R          RDS -> csv/kevin/Writeup21/
  plot_barplots_claude.R      csv -> fig/kevin/Writeup21/
```

**Two rounds.** The fixed value of each axis is defined as *the last level at
which all three methods still do reasonably well*, which is only known after
seeing the curves, so:

1. **Trailblazing round**, on the laptop: both sweeps at 7 levels × 2 replicates,
   with **provisional fixed values: t2 Gini 0.5 held during the AED sweep, mean
   squared AED 0.4 (`h² = 0.6`) held during the Gini sweep**, and `τ_δ = 0.6`
   (Section 2.4). Its outputs are the curves at provisional settings plus, per
   axis, the last easy level.
2. **Final round**: the fixed values are set to the levels found in round 1, the
   level lists are re-spaced if the comparators fall away in a narrow window,
   and both sweeps run at 20 replicates on Hyak (which first needs a `COSPAR_ENV`
   built there). If the new fixed values change the shape of the curves, a
   second laptop pass at 2 replicates comes before the Hyak run. The `τ_δ`
   supplementary row (Section 2.4) runs at the chosen middle levels in this
   round.

Each driver runs the R side in one process and shells out to the `cospar` conda
environment (`COSPAR_ENV`) per dataset via `system2()`. Per dataset: generator under
a second, PCA a second, CYFER about 20 s at 100 clones and 10 features, CoSPAR about
a minute at 4,000 cells, gene statistics a few seconds. Seven levels × 2 replicates
× two axes = 28 datasets, plus calibration, is about one to two hours on the laptop
for the trailblazing round, run sequentially so the progress file is readable; the
supplementary row adds 16 datasets in the final round.

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
per level × replicate × method, with the realized statistics and diagnostics of
Section 7), `calibration` (target, calibrated parameters, realized mean),
`params`.

**Figures.** Grouped barplots in the style of `Writeup17b_barplot-*.R`: one bar per
method at each of the seven levels, height the mean metric of Section 6.1 across
replicates with the individual replicates overplotted as points (an SD over two
replicates is not worth drawing), x-axis labelled with the target statistic and the
realized mean beneath it (for 4B also the per-clone 5th–95th percentile range).

## 9. Questions for Kevin

1. **Keep the clone-varying spread on the non-causal coordinates only?** The
   generator's `spread_variation` switch (Section 2.1) defaults to `"noncausal"`,
   and the demo (Section 4.1) says the change is free on the axis that matters
   and removes the confound: the per-clone AED range is unchanged (0.05–2 at the
   top level, medians and percentiles within 0.05 of the isotropic version at
   every level), the mean still equals `1 − h²`, and the spread-size coupling
   drops from a Spearman of 0.5–0.65 to about 0 (0.07 in the latent-only check,
   the irreducible sample-level Jensen effect). The costs: (a) the Gini now falls
   rather than rises along the AED axis at fixed `s`, so the calibrated `s`
   climbs to about 1.25 at the top level (a larger within-clone spread of `Z`,
   hence bigger jackpot cells, at a Gini that is nonetheless held at 0.5); (b)
   the clone-to-clone *variation* in AED is then carried by fate-irrelevant
   coordinates, so a clone's own AED no longer says anything about how
   heterogeneous its fate potential is, only the level's mean does. **Suggestion:
   keep `"noncausal"`.** The figure's claim is about heterogeneity per se, and
   under `"isotropic"` a reviewer could say that the AED axis is partly a
   clone-size axis (diffuse clones are the big clones, so lineage-DE's high group
   is also the diffuse group) and that the per-clone AED distribution shown under
   the barplot is partly a size distribution; under `"noncausal"` the Methods can
   state that AED and clone size are independent by construction. Cost (b) is
   cosmetic because the analysis is at the dataset level, and cost (a) is handled
   by the calibration. If Kevin prefers the isotropic version, the switch is one
   argument and the demo results for it are in the same CSV.
2. **`τ_δ = 0.6` for the trailblazing round.** Section 2.4 fixes the
   clone-specific t2 drift at the between-clone spread `τ` of the provisional
   middle point, as an absolute value shared by every level of both sweeps
   (about one within-clone standard deviation per coordinate). Confirm, or name a
   different value; the supplementary row's grid `{0, 0.3, 0.6, 1.2}` scales with
   it.

## 10. What I am uncertain about, stated plainly

- Whether CoSPAR, given the same embedding, does nearly as well as CYFER on 4B. Its
  coherence prior assumes transcriptomic neighbours share fate, which is *true* under
  this generator (fate is a linear function of the embedding). If so, CYFER's edge is
  the extreme-Gini regime and the zero-clone handling rather than heterogeneity as
  such, and the framing should say that.
- How strongly the clone-specific drift `τ_δ` hurts CoSPAR in practice. The
  argument in Section 2.4 is from reading the map construction, not from running
  it. **If CoSPAR performs far too poorly even at the easy settings of either
  sweep, `τ_δ` will be decreased**: the main panels are meant to show the Gini
  and AED axes doing the work, not the drift. Conversely, if CoSPAR turns out
  insensitive to `τ_δ`, its weakness on these data is the uniform barcode link
  and the extinct clones, and the paper's framing should say that rather than
  "smooth continuum".
- Whether the spread-size coupling of Section 2.1 would read as a feature or as a
  confound to a reviewer. The `"noncausal"` default removes it from the data;
  **[Q1]** decides whether that stays. Under `"noncausal"` the per-clone AED
  variation is fate-irrelevant by construction, which is the honest price of the
  independence and should be said in the Methods in one sentence.
- At the top AED level `τ = 0`: clones have no centre and clone identity at t1 is
  spread alone. Lineage-DE's two groups then differ in variance, not in mean, so
  its statistic should sit near 0 for every gene; that is the intended failure, but
  it means the top level is qualitatively different from the rest of the ladder,
  not just harder, and the text should say so. Under `"noncausal"` the calibrated
  `s` is largest at this level (about 1.25), so the within-clone spread of `Z` is
  also largest there; the level is extreme in two ways at once.
- Whether a structural set of 100 loading genes is too easy. Real expansion programs
  have weak effects on many genes; a version with 300 weak-loading genes would be a
  harder and more realistic supplement, and the operational truth handles it without
  any change to the metric.
