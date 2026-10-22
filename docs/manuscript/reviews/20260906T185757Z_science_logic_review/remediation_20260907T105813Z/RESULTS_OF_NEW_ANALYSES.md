# New analyses run in response to the review, 2026-09-07

Four analyses were run. Two settle a claim the manuscript was making, one
strengthens the strongest surviving result, and one explains a null the review
raised. Every number below was read from its producing file in the same session.

---

## 1. Are the two gene classes distinguishable? **No, at this n.**

`RNA-seq/results/class_difference_test/cdt-20260907T105813Z/` · job 21489259 ·
prespecification sealed and hash-verified at runtime · full reasoning in
`VERDICT.md`.

The manuscript's organizing sentence was that the genetic and disease-state maps
are "largely non-overlapping". They are — and that cannot fail: 92.9% of a
437-gene set is expected to miss a 1,066-gene set in a 14,920-gene background by
arithmetic, against 95.4% observed. This test replaced it with one that can.

| | hepatocyte specificity | protein stage response, residualised |
|---|---:|---:|
| matched pairs | 242 | 54 |
| SMD after matching | +0.005 | −0.001 |
| mean paired difference | −0.204 | −0.037 |
| sign split | 123 / 119 (p = 0.85) | 24 / 30 (p = 0.50) |
| Wilcoxon p | 0.088 | 0.288 |
| corrected paired null p | **0.064** | 0.154 |
| minimum detectable paired d | 0.180 | 0.381 |

**Verdict: not distinguishable.** The prespecified rule also required the
direction to match the prediction recorded in advance (colocalized genes *more*
hepatocyte-specific); the observed sign is the opposite.

Three things this turned up that matter more than the verdict:

- **The prespecified permutation null was anti-conservative** and its p = 0.0036
  must not be quoted. Relabelling genes and rematching destroys the covariate
  imbalance that defines the real comparison: the classes differ in expression by
  SMD 1.576, so only 242 of 415 colocalized genes have a partner in the caliper
  at all, whereas relabelled pseudo-classes match easily. Measured, that null has
  sd 0.0684 against the statistic's own standard error of 0.1091 — **1.59× too
  narrow.** The correct paired null (randomise which member of each realised pair
  carries the label) has sd 0.1097 and gives p = 0.064.
- **The obvious version of the protein comparison is circular and would have
  produced a confident wrong claim.** Raw, it gives d = −0.554, p = 0.0002.
  Disease-state genes were *defined* by their bulk RNA response and protein
  tracks RNA. Residualising on |bulk logFC| removes it entirely (p = 0.154).
- **The one large, reproducible difference between the classes is
  methodological**: colocalized genes sit 4.3 log2 units higher in liver
  expression, because colocalization requires a well-powered eGene.

Positive control passed (hepatocyte vs fibroblast markers separate by 2.67 log2,
predicted > 2.0 in writing beforehand), so the null is real and not a broken
measurement.

**Manuscript consequence: "separable" cannot be used of these two classes, and
"complementary"/"largely non-overlapping" cannot be presented as findings.**
Applied to `01_intro.md` and `fig5_discussion.md`.

---

## 2. Coding-mass asymmetry: the three owed checks, **all passed**

`GWAS/finemapping/results/fig2_robustness/20260907T105813Z/` · job 21489350 ·
62 refits, seed 20260817, 10,000 permutations each, reusing the adopted
producer's clustering, estimand and permutation unchanged. The reference arm
reproduces the adopted ratio 6.5494 exactly before anything varies.

| check | fits | ratio range | worst p |
|---|---:|---|---:|
| leave-one-study-out, direct | 13 | **4.99 – 9.58** | 8.0e-4 |
| leave-one-study-out, enzyme | 18 | 5.53 – 7.12 | 2.0e-4 |
| leave-one-cluster-out, direct | 24 | 5.57 – 6.84 | 3.0e-4 |
| window size 0.25–5 Mb | 5 | 5.23 – 6.55 | 2.0e-4 |
| ancestry-matched LD blocks | 1 | **6.50** (33 vs 588 clusters) | 1.0e-4 |

The check that mattered was leave-one-study-out on the direct class, because
MVP_NAFLD_AMR alone contributes 95 of the 208 direct credible-set instances.
**Dropping it raises the ratio to 9.58 rather than removing it.** No single
study, cluster, window or clustering rule carries the result.

This is now the best-supported quantitative result in the repository.

---

## 3. A third of the program registry is mislabelled, and the cause is measurable

`Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/`
· `FINDINGS.md`.

The review found the two featured programs are labelled hepatocyte and are not.
The question was how many. Scoring all 117 programs' fixed members against
donor-level per-lineage pseudobulk, z-scored across lineages:

**80 of 117 (68%) match their declared lineage. 37 (32%) do not.**

Validated first against the donor-paired gold standard — the index independently
returns Fibroblasts for ECM/IGFBP7 and Cholangiocytes for BICC1, and the script
exits non-zero if it does not.

The mismatches are systematic: **17 of the 37 are programs discovered outside
hepatocytes whose content is hepatocyte**, several named for it outright
("Hepatocyte xenobiotic (CYP2C)" labelled fibroblast; "Bile-acid transport
(ABCB11)" and "Plasma lipoprotein (APOE)" labelled cholangiocyte; "Plasma
protein (APOA1)" labelled macrophage *and* separately fibroblast).

The mechanism is ambient RNA, measured not asserted:

| program set | n | median decontX ambient burden |
|---|---:|---:|
| non-hepatocyte label, hepatocyte content | 17 | **0.432** |
| non-hepatocyte label, content matches | 60 | **0.061** |

Mann-Whitney one-sided p < 1e-5, rank-biserial **+0.998** — the sets barely
overlap.

**Consequences.** Figures 4C–4E label programs by content with the discovery
subset kept in the legend as provenance. The two featured programs are the
general case, not exceptions. And the vocabulary null (26 of 114
fibrosis-supported against 18.4 expected, p = 0.134) now has a cause: a third of
the vocabulary does not measure what its name says. That is an explanation, not
a rescue — the registry still does not beat matched random sets on fibrosis. The
registry itself is unchanged: nothing was rediscovered, reweighted or renamed,
and no disease outcome was read.

---

## 4. The 113/27 versus 114/26 program counts: reconciled, not superseded

Recomputed directly. The unweighted vocabulary arm additionally admits T-cell
module 9 (`hotspot_tcells_eaa9c262b58d9323`, weight coverage 0.783), giving
114 tested in a 228-test family against the weighted arm's 113 in 226. **The 26
supported are a strict subset of the 27.** The one program that differs is
cholangiocyte module 18, "CFTR / ion transport", at weighted q = 0.0369 and
unweighted q = 0.0541 — one borderline call either side of 0.05, moved by the
scoring and family size. The extra tested program is supported in neither arm.
`RESULTS.md` called this pair superseded and `HIGH_IMPACT_PATTERNS.md` quoted it
live; both now carry the reconciliation.

---

## 5. The F0→F1 contrast is almost entirely control-versus-disease

`figures/candidates/f0-reference-arms-20260907T105813Z/analysis/f0_arms/` ·
job 21489343 · `FINDINGS.md`.

F0 is 104 participants: **31 source controls and 73 F0 disease biopsies**,
nested in cohort (GSE130970's F0 is entirely control; GSE162694's is entirely
disease, because its controls carry `fibrosis_stage = NA` and are dropped
upstream). F1–F4 are effectively clean.

| contrast | A: all F0 (n=104) | B: disease-only F0 (n=73) | C: all F0 + control-status term (n=104) |
|---|---:|---:|---:|
| **F0→F1** | **1,471** | **7** | **258** |
| F1→F2 / F2→F3 / F3→F4 | 548 / 522 / 244 | 548 / 522 / 244 | 532 / 506 / n/a |
| F2 / F3 / F4 vs F0 | 3,904 / 6,320 / 6,424 | 3,122 / 4,949 / 4,336 | 3,238 / 5,728 / 5,038 |

**The collapse is the confound, not the smaller n.** Arm C keeps all 104
participants and only adds a control-status term: 1,471 becomes 258, so 82% is
removable at unchanged n. Arm B leaves 7. Arm C is a partial adjustment because
control status is nested in cohort; the gap between 258 and 7 is what the
nesting hides. All-gene effects still correlate (rho 0.775), so this changes
what is resolvable, not the direction.

**Two claims retired.** "Fibrosis remodeling expands early and then stabilizes"
(`HIGH_IMPACT` row 20) rested on 1,471 → 548 → 522 → 244; on disease-only F0
that is **7 → 548 → 522 → 244**, which is not an early expansion. And the
GNMT and MAT1A F0→F1 stage effects do not survive: GNMT q = 0.0085 → 0.216
(disease-only) and 0.141 (control-status term at unchanged n); MAT1A
q = 1.3e-4 → 0.105 and 0.027. Point estimates move 15–20%. **The pooled
844-participant disease-versus-control result for GNMT/MAT1A/CYP2C19 is
untouched** — different contrast, still valid.

**Defect fixed in passing.** `build_pi_stage_extensions.R:112` builds its t
critical value with an unnamed `qt()` vector and indexes it by gene name, so all
257,070 confidence intervals in the adopted stage-extension output are NA. The
rerun names the vector and fails closed if any interval is missing.

---

## 6. Which trait class drives the colocalized map (and the headline it supports)

`GWAS/finemapping/results/coloc_trait_scope/cts-20260907T105813Z/`.

`HIGH_IMPACT_PATTERNS.md` row 12 recorded 465/473 enzyme versus 22 direct and
marked it stale on the adopted release. The adopted gene-level table carries
`driving_trait`, so it is directly recomputable:

| set | genes | strongest evidence from a liver-enzyme trait | from direct disease / liver fat |
|---|---:|---:|---:|
| SuSiE PP.H4 >= 0.5 (primary) | 462 | **448 (97.0%)** | **14 (3.0%)** |
| ABF PP.H4 >= 0.5 | 823 | 758 (92.1%) | 65 (7.9%) |

By trait: ALT 150, AST 164, GGT 134; NAFLD 8, PDFF 5, NASH 1.

**Scope caveat, checked:** the file holds one row per gene (18,977 rows, no
duplicated identifier, one missing id), so `driving_trait` is each gene's
*strongest* colocalization. This measures which trait class supplies each gene's
best evidence, not that a gene lacks weaker evidence from the other class. An
earlier version of this analysis reported a "both = 0" column; that is true by
construction and was removed rather than reported.

### The argument this completes

Three verified numbers now form one chain:

1. Direct MASLD and liver-fat credible sets carry **22.5%** of their fine-mapped
   posterior probability in protein-altering or canonical-splice variants
   against **3.4%** for liver enzymes; robust across all 62 refits.
2. Protein-altering variation acts on protein function, not transcript
   abundance, so it need not produce an expression QTL and cannot be recovered
   by expression-based colocalization.
3. The colocalized map behaves exactly as that predicts: **448 of 462** genes
   get their strongest evidence from a liver-enzyme trait, **14** from a direct
   disease or liver-fat trait.

Consequence: expression-based colocalization applied to liver traits returns
predominantly a **liver-enzyme regulatory map**, and the most disease-proximal
part of inherited MASLD risk is the part it is least able to see. The sparse
overlap with the disease-state transcriptome is a consequence of this, not a
finding — which is what the original headline was trying, and failing, to be.

Step 2 is a standard molecular-genetics inference, not something tested here;
steps 1 and 3 are measured. Woven into `abstract.md`, `01_intro.md`,
`draft/fig2.md`, `fig5_discussion.md`, and added to `PAPER.md` as the central
quantitative claim, with the retired framings and the "PIP mass"/"coding mass"
jargon added to the forbidden-language table.
