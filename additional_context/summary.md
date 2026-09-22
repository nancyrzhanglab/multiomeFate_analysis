# additional_context/ — Reference Material

This folder holds reference material for the multiomeFate analysis project (CYFER).
Read this file before opening any file here — it exists so that future sessions do
not re-read the same sources.

The material falls into 2 categories:
- **Foil methods**: methods CYFER is benchmarked against in the simulations.
- **Figure planning**: the working plan for the Nature Methods resubmission.

The Nature Genetics reviewer report and response letter are **not** here; they live in
`PAPER_REPO/additional_context/`, and their simulation-relevant points are distilled
in `PAPERS_DIR/simulation-planning.txt` and `kevin/Writeup_Simulations/README.md`.

Entries are keyed by `generate-paper-id` citation key. Superseded material is marked
`[SUPERSEDED by X]` rather than deleted.

Last updated: 2026-09-22

---

## Foil methods

### wang2022cospar — Wang et al. (2022)

**Full citation:** Wang, S.-W., Herriges, M. J., Hurley, K., Kotton, D. N., & Klein,
A. M. (2022). CoSpar identifies early cell fate biases from single-cell transcriptomic
and lineage information. *Nature Biotechnology*, 40, 1066–1074.
doi:10.1038/s41587-022-01209-1. File: `wang_2022_cospar.pdf`.

**What it does:** Infers a cell-by-cell transition map `T(t1, t2)` between two time
points from a cell-state similarity graph plus lineage barcodes, by iterating
{enforce observed intra-clone transitions → smooth over the kNN graph → sparsify →
renormalize}. From the map it derives per-cell fate probabilities, a two-fate **fate
bias** in `[0, 1]` (0.5 = neutral, also the value for cells outside the map), fate
potency (entropy), fate coupling and hierarchy.

**Key mathematical ideas:**
- Objective `min_T ‖T‖₁ + α‖LT‖₂` subject to the clonal constraint; `L` is the graph
  Laplacian. Sparsity: each cell reaches few states. Coherence: transcriptomic
  neighbours share fates. The user-facing knobs are the smoothing depths
  (`smooth_array`) and the sparsity threshold, not `α`.
- Clones seen at both time points fully constrain the problem
  (`infer_Tmap_from_multitime_clones`). Clones seen only at t2 need a joint estimate
  of the initial clone matrix, initialized by optimal transport or HVG pseudo-clones
  (`infer_Tmap_from_one_time_clones`).
- **Fate-associated genes**: threshold the fate bias into two progenitor groups, then
  Wilcoxon rank-sum with BH at FDR 0.05, ranked by fold change
  (`tl.progenitor` → `tl.differential_genes`). No CoSPAR-specific gene statistic.

**Temporal & methodological context:** 2022; benchmarked on Weinreb 2020 LARRY
hematopoiesis (days 2/4/6), Biddy reprogramming, and iPSC lung differentiation.
Robust to as few as 30 barcodes, to barcode homoplasy, and to clonal dispersion.
Stated limits: learns only the *average* fate bias of observed states, cannot
separate division-rate from differentiation-rate bias, and depends on the similarity
metric and smoothing depth.

**Project relevance:**
- **Foil / alternative method**: yes. One of the three columns of the paper's Table 1
  and the "existing method" in the planned Figure 4 sweeps. Its coherence prior
  assumes neighbours share fate, which is exactly what high within-clone
  heterogeneity on the fate axis violates.
- **Prior use**: `kevin/Writeup14_simulation/*cospar*` ran it on the semi-synthetic
  priming/plastic data (supplementary Fig S3-3); the `run-cospar` local skill
  replaces that pipeline.
- **What it does not model**: exponential growth, selection, or lineage extinction
  (Table 1's CYFER-only rows).

## Figure planning

### mockup2026figures — Figure mockup for the Nature Methods resubmission

**File:** `Mockup of figures of Nancy-Sydney paper.pptx` (67 MB, so git-ignored and Dropbox-only; 67 slides; the first
7 are the planned main figures, no speaker notes).

**Slide-by-slide (text verbatim where short):**
- **Fig 1** — "CYFER addresses skewed clonal growth and intraclonal heterogeneity
  that hinder treatment-response prediction." C–D clone-size evenness histograms and
  quantification; E–F intraclonal heterogeneity histograms; G CYFER overview. The
  problem statement moves ahead of the experiment.
- **Fig 2** — Clonal growth and diversity across treatments and time points. D
  histogram of Gini indices, "comparable or even more extreme than other people's
  data"; E violin of intraclonal heterogeneity.
- **Fig 3** — "CYFER estimates future clone size across conditions": A–D one example
  per setting (simulated data, hematopoiesis, cancer treatment from Schaff et al.,
  the authors' own data); E barplot of predicted-vs-observed clone size accuracy
  across datasets. Schaff et al. is a *new* external benchmark.
- **Fig 4** — "CYFER recovers fate-associated genes across conditions. A) Simulation
  of varying clone size skewness and comparing the overlap of DEGs found by our and
  existing methods versus ground truth. B) Simulation of varying intraclonal
  heterogeneity and comparing overlap of DEGs." Two rows of ~7 thumbnails each.
  **This is the Writeup21 deliverable.**
- **Fig 5** — "CYFER identifies priming and plastic clonal growth" (real data).
- **Fig 6** — "Example clones on selection driven or adaptation driven in dabtram
  and cisplatin."
- **Fig 7** — untitled full-page figure (molecular programs / clinical, carried over).

**Project relevance:** Fig 4 fixes the design target for `kevin/Writeup21_new-simulations/`:
two seven-level sweeps (Gini of clone sizes; within-clone heterogeneity), metric =
DEG overlap with ground truth, methods = CYFER, CoSPAR, lineage-level DE.
