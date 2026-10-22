# MASLD model registry

**Every model this project has built, its state, and what may be claimed for it.**
One row per model. Update this file when a model is built, changed, released or retired.

Last updated 2026-09-08 20:35 UTC (**masld-liver-chromatin-state-v1.2 SHIPPED**: v1.1 + per-region H3K27ac profile head, primary form RRR-offset + cis, in-cohort +0.1413, transfers to 39 multiome donors' ATAC at p 0.001 with the shuffled pipeline at null; `test_release.sh` ALL CHECKS PASSED. Phase 3b: no context term advances (well leakage). Zero-shrinkage RRR instability found and fixed.)

---

## Released — a stranger can run these

| model | what it does | input | released at | tests |
|---|---|---|---|---|
| **masld-severity-v1.2** ⚠contract refuses GSE193066 — RESOLVED 2026-09-08, card §4.2 amended, bound unchanged | orders bulk liver RNA-seq by fibrosis severity (`latent_severity`, bitwise = v1.1) **and by NASH-CRN activity (`latent_activity`)**; prints a label-free shift diagnostic | counts, 26,629-gene axis | `release/masld-severity-v1.2/` | `test_release.sh` |
| **masld-liver-chromatin-state-v1.1** | **development-grade**: ten RNA-predicted directions of one cohort's H3K27ac covariance (heads bitwise = v1) + a fixed **steatosis-associated chromatin axis** (OOF partial +0.636); exact-recipe accuracies k1 +0.596 / k2 +0.612; ⛔ no direction reproduces externally, k3/k5/k10 are sex | counts, 42,163-gene axis | `release/masld-liver-chromatin-state-v1.1/` | `test_release.sh` ALL CHECKS PASSED (job 21503182); manifest generated after the genome-wide and cis/trans lines; final test job queued |
| masld-liver-chromatin-state-v1 | superseded by v1.1 (`SUPERSEDED.md`); ⛔ do not quote its card | counts, 42,163-gene axis | `release/masld-liver-chromatin-state-v1/` | `test_release.sh` |
| **masld-liver-chromatin-state-v1.2** ⭐ | v1.1 (ten state heads + steatosis head, bitwise) **+ a per-region H3K27ac PROFILE head** (`--profile`): primary form RRR (shrinkage 1e-7, rank 20) fixed + a ridge on the ≤10 nearest genes within 100 kb on its residual; in-cohort OOF skill **+0.1413 [+0.1125, +0.1686]** (+0.022 over RRR every fold; the gene block is LOCAL: +0.036 over random genes); **transfers** to 39 GSE296875 donors' snATAC: reliable-set mean ρ +0.256 vs within-well null +0.071 (p 0.001), +0.251 well-residualised, shuffled-pairing pipeline at null; four alternative forms ship (`states_cis` ties in cohort and is more robust under raw transport). Nothing nonlinear beats it. | counts, 42,163-gene axis; marginal transport needs ≥10 samples | `release/masld-liver-chromatin-state-v1.2/` | `test_release.sh` ALL CHECKS PASSED (states bitwise v1.1; every transfer ρ reproduced to 0.0 from the release code; sources sha-verified) |
| masld-severity-v1.1 | superseded by v1.2 (frozen; `SUPERSEDED.md`) | counts, 26,629-gene axis | `release/masld-severity-v1.1/` | `test_release.sh` |
| **masld-panel-v1** | same ordering; ⛔**qPCR deployability NOT established** | counts **or Cq** | `release/masld-panel-v1/` | `test_release.sh`, 4/4 |
| **masld-liver-protein-severity-v1.3** | orders liver proteomes by fibrosis | protein quant | `release/masld-liver-protein-severity-v1.3/` | `test_release.sh` |

### masld-liver-chromatin-state-v1 — RNA to chromatin state (2026-09-07) — ⛔ CARD NUMBERS SUPERSEDED, v1.1 in preparation
Ten PLS components of GSE267145 H3K27ac (residualised on the per-donor concentration block), each
predicted from bulk RNA by a ridge on an 81-dim frozen PCA basis. 29 MB.
- ⛔⛔ **The v1 card's accuracies were copied from a different arm than the shipped heads** (handoff
  finding 1). Exact-recipe out-of-fold partials (`chromatin-state-exact-recipe-eval-20260907T230940Z`,
  aligned components, given grades+sex+14 RNA descriptors): k1 **+0.596** [0.40, 0.74] (card 0.666),
  k2 **+0.612** [0.39, 0.75] (card 0.770), k4 +0.604, k5 +0.472, k6 +0.448, k7 +0.438, k9 +0.399,
  k10 +0.347, k8 +0.297, k3 +0.250. All ten survive BH q<0.05; k3 fails the family-wise max-statistic
  null (p 0.113). Whole-matrix reconstruction of the SHIPPED low-rank form **+0.045** [0.030, 0.058]
  (card quoted +0.1194/+0.1030, which belong to unshipped models). Shipped `oof_residual_sd` wrong for
  every component (k2 0.527 shipped vs 0.866 exact).
- ⛔⛔ **"Orthogonal to histology" was an alignment artefact.** Aligned, k2 correlates +0.55 with
  steatosis, +0.50 lobular, +0.43 ballooning, +0.31 fibrosis: it is the activity axis. k5 (−0.40), k3
  (−0.27), k1 (−0.23) carry sex in training. Only k1–k3 have a stable identity across folds (|cos| ≥ 0.85);
  k4–k10 are re-assigned between folds.
- ⭐ Ten unsupervised in-fold RNA PCs match the supervised ridge for k1–k7 (within ±0.12); the model's gain
  is only in k8–k10, the least stable components. Ridge penalties sit at the grid ceiling.
- ⭐ **Fixed steatosis-associated head (Step 2 target, same PCA-81 ridge): +0.636 [0.47, 0.76]**, the one
  motivated target; ships in v1.1. Female-only (n=85): k1 +0.695, k2 +0.638.
- ⭐ RNA-predictable chromatin covariance is high-dimensional: 26/30 aligned K=30 components have OOF lower
  bounds > 0; stable subspace 8.3 dims at K=10, 22.2 at K=30. "Ten states" was an arbitrary cut.
- ⛔⛔ **EXTERNAL UTILITY NOT DEMONSTRATED** (`chromatin-state-external-utility-20260907T224715Z`, sealed,
  8 RNA cohorts, 847 samples, 0/10 axes reproducible): k1 test-retest ICC **+0.008** on 58 repeat-biopsy
  pairs; no axis beats a random direction in the same 81-PC space (null 0.33 ± 0.20); **k3, k5, k10 are
  sex axes** externally (|ρ| 0.43–0.53); k2, k7, k8 correlate with fibrosis externally (median |ρ|
  0.40–0.43); k7, k8 track age; ten-axis geometry does not transfer (−0.05 vs +0.26–0.49 for random axes).
- ✅ **Review 2026-09-08 (`review-quantification-sensitivity-20260908T102607Z`): the external negative SURVIVES.** Scoring the
  same 726 samples from STAR/featureCounts vs kallisto through v1.1 agrees at ρ 0.90–1.00 per direction (median k1 0.968,
  k2 0.992) although the two quantifications agree only at 0.56 per sample gene-wise: the 81-PC projection is
  quantification-insensitive, so the negative is not a pipeline artefact. ICC re-derived from ANOVA sums to 2e-16; pairing
  rebuilt independently (58 pairs); second-seed null identical. ⭐ Reference row the lane lacked: the release's own PC4 has
  test-retest ICC **+0.81** and five PCs exceed 0.60 vs the best chromatin direction 0.63 (a sex axis) — reproducible
  donor directions EXIST in this space; the chromatin-covariance directions are just not among them.
- ✅ **Tuning/implementation review 2026-09-08 (`review-parameter-sensitivity-20260908T103644Z`,
  `review-cross-implementation-20260908T104115Z`): the external negative is NOT a tuning or coding artefact.**
  Representation size (NPC 10/20/40/81/n−1) and number of directions (K 3/5/10) were swept, each judged against a
  random-direction null **regenerated in its own representation**: no configuration's best direction survives BH
  (smallest q 0.148). ⛔ My first pass compared a 10-dim configuration against an 81-dim null and printed a false
  "tuning artefact" flip on a 0.0013 margin — repaired in `prespec/ADDENDUM_01.md`. ⭐ In 6 of 8 configurations the
  most reproducible direction is k10, which correlates +0.46 to +0.59 with **sex**.
  ⭐ **The estimator provably cannot matter for this endpoint**: shipped ridge, a grid extended to 1e4, and PLS
  regression give partials identical to 6 dp — raw predictions differ by up to 2.15, but within each fold they are
  affine transforms (after the endpoint's within-fold standardisation: 6.8e-15, zero rank changes). Not scale-invariant
  for the profile endpoint.
  ⭐ Two unused cohorts (GSE213621 n=367, GSE167523 n=98) reproduce every association with more power: k2/k7 vs
  fibrosis +0.41/+0.34, k3/k5 vs sex −0.60/−0.35, k7/k8 vs age +0.49/+0.34, k1 vs NASH +0.51.
  ✅ Cross-implementation: k2 partial via sklearn/statsmodels |Δ| 0.0000; ten ICCs via MixedLM max |Δ| 0.0122
  (k6 boundary clamp), ordering preserved; lane C's cis arm via sklearn RidgeCV inside its bootstrap interval.
- ⛔ In-sample fit 0.92–0.99 means nothing. Within-batch ordering only; one cohort (n=99); development
  data after ≥26 executions.
- ⭐ **Genome-wide (94,026 regions, sealed, `corgi-genomewide-h3k27ac-20260907T195320Z`): RNA→chromatin (nested RRR)
  participant-differential skill +0.119, Corgi +0.0001 (real = shuffled context); release claim STANDS under the sealed
  rule and under fixed training-fold calibration (ADDENDUM_01: +0.119 vs −0.0004, [0.092, 0.144]).** ⛔ 'ours' there is the
  unshipped RRR form; the shipped rank-10 model recovers +0.045.
- ⭐⭐ **Cis vs trans (`chromatin-cis-trans-20260907T231147Z`, sealed, run 2 is the record): a ridge on each region's ≤20
  neighbouring genes (±100 kb) recovers +0.044 [0.036, 0.051] of participant-level H3K27ac — as much as the global
  81-PC form (+0.048) — and the two ADD: +0.082 together, increment +0.034 [0.025, 0.042].** Cis skill is concentrated
  (1,675 regions > 0.3 vs 170 for trans). ⛔ The shipped model has no gene-level weights and leaves the local component
  unrecovered; this is why raw-gene forms (+0.10–0.12) beat PCA forms (+0.045–0.048) on reconstruction all along.
  ⭐ Next build that is actually justified: a sparse cis head (~96,460 × ≤20 weights, ~15 MB), not a new model family.
  ⭐ 2026-09-08: BUILT nested (`chromatin-cis-head-20260908T102817Z`): +0.1426; beats RRR by +0.0232 paired — see the 2026-09-08 profile-head lane below. ⛔ The ~15 MB-vs-32 GB size argument is dead: RRR stores two factors (~6.6 MB).
  ⛔ Figure 7's 'signal is diffuse/tissue-level' is half right: an additive part is expression-local; Corgi recovers neither.
- **Read as: a development-grade predictor of one cohort's RNA–chromatin covariance, not ten liver
  chromatin states.** v1.1 (heads bitwise identical + steatosis head + corrected card) replaces this entry.

### masld-severity-v1.2 — the parent instrument, with an activity head (2026-09-07)
- ⚠⚠ **2026-09-08 (`review-cross-implementation-20260908T104115Z`): the shipped sparsity contract admits every training
  cohort (172/172, 112/112, 93/93, 78/78, 66/66) but refuses **127 of 164 GSE193066 samples** — the very cohort whose
  test-retest ICC 0.7658 is quoted on the card. The STAR quantification passes at a LOWER rate (23/160), so this is not a
  kallisto artefact: GSE193066 is shallower (median library 1.85e7 vs 2.3–4.4e7) and the bound is depth-sensitive. The
  reliability figure is not wrong, but a user re-running that validation through the released path today would be stopped.
  **Owner decision pending**; the measured bound was not changed.
`latent_severity` is bitwise identical to v1.1 (asserted by `test_release.sh` on GSE268273); every
v1.1 number below stands. Added:
- **`latent_activity`** — `NAS_head_464` from `executions/two-axis-gate-20260904T230306Z` exported as
  primal weights (w = Aᵀα; reproduces the gate's predictions to 5.7e-13,
  `executions/severity-v12-activity-head-20260907T162936Z`). LOCO on **464 SAMPLES** from four cohorts
  (GSE240729 carries no NAS). Sealed gate **PASS at the minimum 2 of 4** (margins +0.0134, +0.0068 over the
  paired MDE); vs NAS: GSE193066 first 0.4367 (n=106), GSE225740 0.6054 (93), GSE281797 0.4900 (91),
  GSE267145 0.7987 (99, frozen-weight scoring under the 2026-09-02 carve-out; this number stays out of
  `release/`). ⛔ **No two-axis claim**: Condition B (increment over a 4-column clinical floor on
  GSE281797) failed, 0.169 vs MDE 0.226. ⛔ Ordering only; not a subtype.
- **Shift diagnostic, NOT a gate** — T², Q, per-PC displacement, ribosomal/mitochondrial z on the
  release plane. ⛔ The prespecified T²/Q **refusal bound was measured and NOT shipped**: it admitted
  GSE281797 (0/94; displacement 6.55 SD, equal to in-contract GSE276114) and refused GSE225740 (73/93)
  for coverage 0.605. `executions/severity-v12-joint-contract-20260907T162640Z/RESULTS.md`; card §6d.
  ⭐ GSE281797's ribosomal shift is +3.5 SD against the five-cohort pool (+15 SD was against the
  99-sample classifier training set) — the classifier failure was a single-cohort envelope problem.

### masld-severity-v1.1 — the parent instrument (superseded by v1.2; numbers unchanged)
Kernel ridge on an EM errors-in-variables latent Kleiner target, 26,629 genes, per-sample quantile
map onto a frozen reference.
- Out-of-cohort Spearman **0.6631** on GSE268273 (109 participants), **69.6%** of a 0.9523 ceiling.
- External frozen AUROC F>=3: **0.9187 / 0.9002 / 0.8290**.
- Test-retest ICC **0.7658** vs the Kleiner label's own **0.1886** — the only comparison a tissue
  model can legitimately win.
- ⛔ It is a WITHIN-BATCH RANKER. The absolute severity scale is conclusively not recoverable.
- ⚠ The EM latent contributes **+0.0008** to ordering (median, 492 matched configs). Its value was
  claimed for reliability, not rank accuracy.
- ⚠ Raw score units do NOT transport between this release and a refit of the same recipe: sd ratio
  to **1.85x**, driven by axis coverage, 6-18% of within-batch pairs ordered oppositely.

### masld-panel-v1 — the deployable panel ⭐
20 gene pairs over **20 wells**; ships a pair list, coefficients and an intercept, nothing else.
- Retention vs the full 26,629-gene model: **0.973 / 1.007 / 0.974** (GSE268273 109, GSE213621 299,
  GSE276114 177 participants).
- Above the detection floor the score is **exactly a weighted sum of delta-Cq** — verified to
  **6.7e-16**, and invariant to a global cycle offset to 1.1e-15.
- Refuses below 0.70 pair coverage. EGFLAM is load-bearing (10 of 20 denominators).
- ⚠ Informative range is **k <= 30**; at k=100 a log-ratio panel is not separable from a random one.
- ⚠ The log-ratio identity is exact only ABOVE the floor; the T2 guard was left FAILING, not weakened.
- ⛔ An ORDERING, not a stage. `interval_calibration = "NONE"`.
- ⛔⛔ **NOT DEMONSTRATED RUNNABLE ON qPCR — 7 of 20 genes are at or below detection in liver.**
  Median library CPM over 9 cohorts: OR2I1 **0.06**, STMN2 **0.22** (a numerator in 4 of 20 pairs),
  GPR158 **0.28** (<1 CPM in 92-100% of samples in EVERY cohort), NWD2 0.34, FAP 0.57,
  ENSG00000271239 0.59 (**92.5%** zeros), MOXD1 0.91. OR2I1 and ENSG00000271239 are not measured
  at all in 2 of 9 cohorts. ⭐EGFLAM, the obvious suspect, is FINE (2.52 CPM, 3.4% max zeros).
  The 6.7e-16 delta-Cq identity is a statement about the ALGEBRA, not about whether the
  transcripts amplify; no thermocycler has run this panel.
  See `executions/panel-abundance-audit-20260905T140637Z/RESULTS.md`.
- ⚠ **OR2I1 swings 1,700-fold across cohorts** (median 0.00 in GSE193066 to **103.79** in
  GSE281797) while carrying the 2nd-largest coefficient. Not biology. ⭐A training-pool
  abundance filter does NOT catch it — only a cross-cohort one does.
- ✅ **The panel is NOT a depth artifact.** Score is on CPM-of-panel, so scaling cancels
  algebraically above the floor; partial Spearman removing depth/entropy/top-50-share moves the
  median cohort by only ~-0.03 and 2 of 9 go UP. ⚠Except **GSE281797**: 0.3709 -> 0.2484, landing
  on the repbench `library_complexity` control (0.2374). Treat its panel numbers as unresolved.
- ⚠ An abundance-constrained refit (released selection code, one line changed) is **within noise**:
  two of three cohorts move in OPPOSITE directions with CIs excluding zero (GSE213621 -0.0940,
  GSE276114 +0.0494, GSE268273 -0.0008), mean far inside the MDE 0.0667. **No superiority claim
  either way**; 4 thresholds x 2 filters tried, uncorrected, not prespecified.

### masld-liver-protein-severity-v1.3
LOO Spearman 0.7291, 80.1% of its 0.9107 ceiling, n=58, single cohort.
⛔ No external replicate exists or can. Characterised, not validated.

---

## Tested and REJECTED — each by its own prespecified kill line

| model | task | result | why it died |
|---|---|---|---|
| ordinal cumulative-logit on RNA | fibrosis | −0.0076, MDE 0.0265 | best head on protein; does not carry to RNA |
| anchor regression (cohort anchor) | fibrosis | like-for-like delta **0.0000** | curve flat, argmax AT gamma=1 which IS ridge |
| pairwise-difference, L2 loss | fibrosis | — | **algebraically IS** a within-cohort-centred ridge (identity to 2.6e-14) |
| pairwise logistic (RankNet) | fibrosis | −0.0131, CI [−0.0235, −0.0027] | measurably worse, not merely indistinguishable |
| absolute conformal (percentile, ordinal) | fibrosis | coverage failure | a LOCATION problem; one arm gave singleton rate 1.000 at coverage 0.000 |
| DEBIAS-M port | — | not built | for a linear head the identifiability constraint makes it a no-op at deployment |

⭐ **Every rejected model added machinery to the existing estimand. Both survivors CHANGED the
estimand** — absolute stage to within-batch ordering, dense standardised model to scale-free contrast.

---

## Survived but not released

| model | state |
|---|---|
| **pairwise conformal** | ✅ a genuine split-conformal object. Pooled coverage 0.8902 over 386 participants, call rate 0.706, accuracy 0.845. Not packaged as a release. |
| **NAS activity head** | ⚠ Condition A passed 2 of 4 at the sealed minimum; Condition B FAILED on power (+0.1689 at 2.09 SE, needs 2.80). **n ~ 163 would certify** [Estimate]; GSE281797 is n=91. |

---

## Under active examination — DO NOT QUOTE

| model | state |
|---|---|
| **Geneformer V2-104M** (Apache-2.0) | ⛔ **two of my own runs disagree.** L12 mean-pool: 0.3618 vs twin 0.3160 (+0.0458) in the readout sweep, but 0.3122 vs 0.3156 (−0.0034) in `score_neural.py`. Same seed, same folds. **Unresolved — no number from this arm is usable until it is.** |
| TranscriptFormer (MIT), UCE (CC-BY-4.0) | not run |

**What IS established about Geneformer:** the forward pass is correct — masked tokens are recovered
at rank **3-24 of 20,275** against a chance baseline of ~10,137, so the pretrained weights are real
and working. ⚠ **The encoder sees at most 2,046 of 26,629 axis genes = 7.7% of the axis.**
⛔ **Exposure: 0 of 15 bulk cohorts reach `clean_declared`** — all `unknown`. Any foundation arm is
DESCRIPTIVE-ONLY, and the pre-registered exposure-stratified check is NOT EXECUTABLE because that
stratum is empty. The random-init twin is therefore the load-bearing control.

---

## Development lane 2026-09-08 — GSE296875 donor x lineage ATAC, and what context adds

`executions/gse296875-lineage-peak-matrix-20260908T105743Z` (job 21517427, COMPLETE) and
`executions/gse296875-context-regulatory-increment-20260908T110722Z` (main pass COMPLETE, null running).

**The peak x unit matrix now exists.** 273 donor x lineage units x 306,706 peaks, counted from the
deposited fragment files. The counting convention was decided by measurement against the authors' own
peak assay over the same 68,398 nuclei, and it went **against** the sealed prediction:

| convention | median across-peak Spearman vs the authors' donor matrix | total counts |
|---|---:|---:|
| **fragment overlap** (adopted) | **+0.99987** (min +0.99961) | 340,943,240 |
| Tn5 cut sites | +0.98928 | 626,487,979 |

Reference total 341,091,266, so the adopted matrix lands 0.043 % below it. The authors' peak assay
counts fragments overlapping peaks, not transposition events, despite Tn5-adjusted coordinates.

**Context increment, 2,000 genes, 39 donors, out of fold (⚠ null and sensitivities still running —
DO NOT QUOTE until they land).** Skill is measured against the training-fold mean, so the baseline
scores exactly 0:

| model | skill | increment over the previous model |
|---|---:|---|
| M1 promoter accessibility only | +0.160 | — |
| M2 20 cis peaks, hepatocyte lineage | +0.286 | **+0.127** [+0.083, +0.162] |
| M3 + the same peaks in all lineages (**context**) | +0.296 | +0.009 [−0.004, **+0.027**] |
| M4 + donor transcriptional state | +0.378 | **+0.082** [+0.073, +0.091] |
| M5 + state x accessibility interaction | +0.410 | **+0.032** [+0.025, +0.044] |

Two sealed predictions failed: lineage context does **not** clear zero (P2), and the interaction
**is** resolvable at n=39 (P3 predicted it would not be). The interaction's null has not returned.

⛔ **A defect found and repaired before any outcome was read.** The prespec's endpoint (pooled
out-of-fold Spearman) contradicted its own guard G4 (the baseline must score 0). The training-fold
mean differs per fold, so pooling gave the **baseline itself** |rho| up to 0.608. The endpoint was
replaced by the skill definition lanes A and C already use. `prespec/ADDENDUM_01.md` has it.

**RESULT 2026-09-08 (null + sensitivities in; `RESULTS.md`).** Frozen-fold skill vs training-fold mean: M1 promoter +0.160,
M2 cis peaks +0.286, M3 + lineage context +0.296, M4 + donor state +0.378, M5 + state × cis interaction +0.410.
Lineage context does NOT add (+0.009 [−0.004, +0.027]; P2 failed). The interaction DID add on the frozen folds
(+0.032 [+0.025, +0.044]; state-permutation null max +0.010, 0/200) — **and collapses under leave-one-well-out:
every model +0.02–0.03, interaction +0.003.** ⛔⛔ The 10x well is a batch shared by RNA and ATAC (same nuclei;
median per-feature well R² 0.24–0.26 vs 0.17 chance) and every frozen-fold test donor has well-mates in training.
**No context term advances; nothing becomes a release component; decision 4 is not triggered.** ✅ Checked 2026-09-08: the
teacher–student lane B already ran leave-one-well-out in GSE296875 and its RNA→ATAC whole-matrix skills move by < 0.006 (RRR
+0.021 → +0.018; teacher–student +0.008 → +0.010), so the '−0.013 replicated in GSE296875' line is well-robust. The leakage is
specific to accessibility→expression models on the same nuclei; no other quoted GSE296875 number rests on the frozen donor folds.

## Development lane 2026-09-08 — are any REPRODUCIBLE RNA directions chromatin-predictive?

`executions/chromatin-reproducible-directions-20260908T110210Z`. Directions fitted on GSE193066
test-retest pairs only, frozen, applied to GSE267145. ⚠ final rerun in flight; the numbers below are
from the run whose guards all passed.

- Reproducible directions exist and are **not** outlier artefacts: the best release coordinate PC4
  reads ICC **+0.812** with a drop-one-participant range of [+0.804, +0.851] and Spearman +0.793;
  the top fitted direction reads +0.958 with drop-one [+0.956, +0.964]. Only 1 of 81 coordinates and
  0.23 % of random directions exceed 0.8, so a high ICC here is genuinely rare.
- **The most reproducible non-degenerate direction is a library-complexity axis** (Gini −0.479,
  top-5 % share −0.507, entropy +0.353), consistent with complexity being this project's recurring
  confound. The single highest-ICC direction is refused as degenerate; the next is sex and age
  (−0.513, −0.429).
- **The reproducible directions carry essentially no chromatin.** Whole-matrix out-of-fold skill:
  ten reproducible directions **+0.0004** [−0.0064, +0.0070] against the 81-PC trans arm **+0.048**
  (reproduced from lane C to 4.4e-16) and the in-fold PLS arm +0.089.
- Across the 81 fixed release coordinates, ICC rank vs chromatin-partial rank is **+0.217** — the two
  properties are largely unrelated rather than opposed.

⛔ **Three defects of my own, each caught by a guard and each repaired on the record**
(`prespec/ADDENDUM_01.md`, addenda 01-03):
1. The permuted-fit null returned held-out ICC **+0.916**. Re-evaluated on **shuffled held-out pairs**
   — where no participant identity exists — it still read **+0.944**. The fit was exploiting the
   35-dimensional null space of a rank-46 within-pair scatter in 81 dimensions. Fitting is now
   restricted to 40 dimensions and the shuffled-held-out control is a guard that can fail; the null's
   control now reads +0.022.
2. Pooling five folds' held-out scores imported fold-to-fold offsets as between-participant variance
   (the same artefact as the Phase 3b endpoint above, found the same day). Each fold is now
   standardised on its own training pairs, and the within-fold ICC is reported beside the pooled one.
3. The S3a bootstrap p had **no floor**: with 400 draws it returned q = 0.0, but the smallest
   attainable BH q at rank 1 in a 100-member family was 0.5, so the q<0.05 criterion was unearnable
   and its True verdict was empty. p is now floored at 1/(B+1) and B raised to 10,000.

---

## Resolved 2026-09-08 — severity-v1.2's sparsity contract and GSE193066

`executions/severity-contract-gse193066-20260908T122902Z` (job 21521650). **The bound was not moved.
The card was corrected.**

The shipped scorer admits **37/164** GSE193066 samples (kallisto) and **23/160** (STAR), and
GSE193066 is the cohort §4.2 cites for test-retest reliability. Only **3 of 58 pairs** have both
biopsies in contract, so the figure cannot be restated on in-contract input at all.

What the measurement settled:

- The card's **ICC 0.7658 reproduces exactly** (|diff| 0.0000) through the released path with
  `--acknowledge-out-of-contract`, and acknowledged scores are **bitwise identical** to v1.1's on
  in-contract samples (0.00e+00) — the flag changes admission, not arithmetic.
- **New:** the same computation on an independent STAR quantification gives ICC **+0.7610**. The
  reliability result reproduces across quantifiers, which §4.2 did not previously show.
- On §6b's **deposited** `thin` curve, the median GSE193066 sample sits between the reference and
  p = 0.3, a measured **1.1 %** Spearman loss; the worst sits between p = 0.1 and p = 0.03. Every one
  of those levels was scored **acceptable** by the sealed rule for `thin`, whose own L* is p = 0.03
  (T1 0.38199). The shipped 0.2544 comes from the **`mask`** mechanism — §6b's own words,
  "adversarial rather than realistic" — plus the 1.10x floor clause.

So the cohort is refused by a bound calibrated on the adversarial mechanism, at a tie-structure level
the realistic mechanism measures as costing about 1 %. The bound stays: it came from a rule sealed
before its curve existed, and loosening it after seeing which cohort it excludes is the failure that
rule prevents. Revisiting it needs a **new sealed prespec** that fixes the mechanism before looking,
and that is an owner decision.

⛔ **Two things I had previously reported about this were wrong and are corrected:** the bound is
**not** a range-based "looks like training data" gate — it came from a sealed degradation curve; and
the ICC was **not** computed before the bound existed — **v1.1 ships the same bound**.

---

## Development lane 2026-09-08 — the H3K27ac profile head: what wins, and against which reference

- **GPU nonlinearity lane** (`chromatin-nonlinear-profile-head-20260908T113726Z`, L40S, sealed + labelled
  addendum; guards vs lane C 2.4e-15 and vs lane A 4.7e-15). No neural arm beats its matched linear twin:
  mlp_states10 − linear_states10 = **−0.0014 [−0.0054, +0.0028]**; on 81 PCs the MLP is worse by 0.047.
  ⭐ **Nested-rank RRR on raw genes is the best global form, +0.1194**, beating the 81-PC per-target ridge
  (+0.0480) by **+0.0714 [+0.0549, +0.0872]** — 96,460 independent penalty choices from ~78 participants
  cost more than a shared low-rank constraint buys. ⛔ Its ranks sat on the grid edge ([12,12,12,11,11] of
  1–12); the cis-on-RRR lane widens the grid.
- **Nested sparse cis head** (`chromatin-cis-head-20260908T102817Z`, job 21529239, `RESULTS.md`): whole-matrix
  OOF skill **+0.1426 [+0.1142, +0.1694]**; the nested procedure chose (100 kb, ≤10 genes, **10 predicted
  states** as the global term) in 4/5 folds — the 10 states beat 81 PCs inside a per-region fit by a third of
  the inner error. nested − trans +0.0946 [+0.0761, +0.1116]; reliable set 4,847 regions, **distal**
  (promoter 0.096 vs 0.137; P4 failed in the opposite direction); P2 failed on the upside (the "raw-gene
  ceiling" +0.103 was a per-target ridge, not a ceiling). Reproduction of lane C: 0.00e+00 after fixing a
  seed-literal bug that had made two runs fail at 5e-3 (my verification script had hard-coded the wrong
  lane's seed and "passed"). **Sealed rule: ship v1.2.**
- ⭐⭐ **Post-hoc, labelled, paired on identical folds/targets (`out/posthoc_paired_vs_gpu_lane.json`):
  cis head − RRR = +0.0232 [+0.0150, +0.0310], positive in all five folds** (+0.024 +0.011 +0.019 +0.027
  +0.035). The number a card may quote for "what the cis head adds" is **+0.023 over RRR**, not +0.095 over
  the per-target ridge. cis head − linear states10 = +0.054; RRR − states10 = +0.031.
- **Running / submitted**: `chromatin-cis-on-rrr-20260908T140356Z` (job 21531974; cis on top of RRR with a
  recalibration-matched reference and a shuffled-cis null — its sealed G2 could not pass by its own arm
  definition, Addendum 01; wide (shrinkage, rank) grid, Addendum 02; shared-penalty variants).
  `chromatin-profile-transfer-gse296875-*` (external transfer of the all-99 profile head to 39 multiome
  donors' ATAC, within-well permutation null, shuffled-fit control). `chromatin-discrepancy-map-*`
  (21531731) and `chromatin-measurement-decision-*` (21531732) on the cis-head OOF profile.
- **Measurement decision** (`chromatin-measurement-decision-20260908T113422Z`, job 21532145, `RESULTS.md`): with a
  budget to measure chromatin in 20 % of a fold's participants, the model's own **cis–trans disagreement** picks the
  participants holding **0.229** of the out-of-fold profile error vs **0.192** for a random pick (null sd 0.0095, max of
  1,000 draws 0.226; paired +0.032 [+0.005, +0.056]); histology 0.191, library quality 0.194, predicted extremity 0.211
  (marginal). A triage statement about prediction error, not biology. The within-grade-1 endpoint favours extremity rules
  but is uncorrected over five rules (smallest p 0.015) — suggestive only.
- **Discrepancy map** (`chromatin-discrepancy-map-20260908T113252Z`, job 21532144, `RESULTS.md`): per-participant
  unexplained share `D` ranges **0.59–1.40** of the training-mean baseline and replicates across disjoint region halves
  (+0.997, a low bar for a ratio over 48k regions). ⛔ The predicted concentration confound is **weak**: strongest
  descriptor |ρ| 0.136, whole descriptor block **13.6 %** of `D`'s variance (P2 failed); no grade or sex survives. 86 % of
  a replicable participant-level model error is unexplained; FRiP/ChIP enrichment is not in the fixture. The script's
  "descriptor-driven" decision line is corrected in its RESULTS.md. Not offered as biology.
- **Exact shipped recipe** (`chromatin-v12-exact-shipped-recipe-20260908T163444Z`, job 21532488): the fixed (100 kb, ≤10 genes,
  states) configuration — the one a release ships — scores **+0.1429 [+0.1140, +0.1695]** out of fold, +0.0003 above the nested
  procedure (fold 0's 50-kb choice was the worse one) and **+0.0235 [+0.0153, +0.0312] over RRR**, positive in all five folds;
  folds 1–4 reproduce the deposited nested predictions to 4e-7. **RRR + cis as it would ship: +0.1193 — equal to RRR alone
  (+0.0001 [−0.007, +0.007]) and below states + cis by −0.024 in every fold.** ⛔ Ten neighbouring genes add +0.054 on top of
  ten states but NOTHING on top of a rank-12 RRR, which cannot encode region-specific local effects — so the "cis" block may be
  supplying GLOBAL information a per-region fit cannot get from 81 PCs. `chromatin-cis-locality-control-*` (job 21534430,
  sealed: cis vs random genes vs other-region cis vs 1–5 Mb same-chromosome) decides whether the card may say "cis" at all.
- ⭐⭐ **Locality control** (`chromatin-cis-locality-control-20260908T*`, job 21534430, `RESULTS.md`): **the block IS local.**
  states10 + 10 nearest genes +0.1429; + 10 random genes +0.1071; + another region's cis genes +0.1072; + 10 genes 1–5 Mb
  away on the same chromosome +0.1091; states only +0.1180. **cis − random +0.0359 [+0.0316, +0.0397]**, every fold ≥ +0.032;
  random genes COST 0.011 vs no genes; per-region advantage decays with TSS distance (+0.027 < 1 kb → +0.015 at 20–50 kb).
  The card may say "cis"; the local component is +0.036 over random genes. The RRR + cis null result is a fitting artefact
  (one strong feature absorbs a per-region ridge), not evidence against locality.
- **cis-on-RRR** (`chromatin-cis-on-rrr-20260908T140356Z`, job 21532398, `RESULTS.md`; guards ≤ 1e-14): sealed Q2 fails —
  RRR-as-a-per-region-FEATURE + 20 cis genes = RRR (+0.0005 [−0.006, +0.007]). ⛔ But re-fitting the RRR prediction per region
  COSTS −0.017 [−0.023, −0.012] and the cis genes add back +0.018 [+0.013, +0.022] (20 shuffled cis genes −0.023): the null is a
  fitting artefact of the feature form. Wide (shrinkage, rank) grid +0.001 over lane A's (grid edge immaterial); shared vs
  per-region penalty +0.001 (no difference); 81 PCs on top of RRR −0.033. ⭐ Untested until now: RRR as a fixed OFFSET with a cis
  ridge on its residual — `chromatin-rrr-offset-cis-*` (chained after the transfer job) measures it in cohort and in transfer.
- ⭐⭐ **External transfer** (`chromatin-profile-transfer-gse296875-20260908T161209Z`, job 21532238, `RESULTS.md`; 7 guards, all pass;
  nothing fitted on GSE296875): the all-99 profile head applied to 39 multiome donors' snRNA and tested against the same donors'
  snATAC at 83,627 mapped regions. **states + cis: mean per-region ρ +0.205 on the 4,260 reliable regions vs a within-well null
  of +0.054 ± 0.024 (p 0.001); +0.106 on all regions (null +0.034); +0.196 after well residualisation; the shuffled-pairing
  pipeline sits at its null (p 0.47 all regions; +0.051 vs +0.015 ± 0.021 on the reliable set, p 0.05).** The cis genes
  transfer (+0.041 [+0.036, +0.045] over states alone). Transfer ρ tracks in-cohort skill (decile +0.056 → +0.206). ⛔ Ordering
  depends on transport: marginal (z-score within batch) RRR+cis +0.262 > RRR +0.221 > states+cis +0.205 > states +0.164; raw
  transport states+cis +0.175 > RRR+cis +0.140, RRR collapses to +0.100. T4 (states+cis ≥ RRR externally) failed (−0.016, n.s.).
  ⛔ The within-well null has a POSITIVE mean (well = shared batch) — the global null is not the right null here.
- ⭐ **RRR-offset + cis** (`chromatin-rrr-offset-cis-*`, jobs 21535502 / 21537530, `RESULTS.md`): in cohort **+0.1431 [+0.1141, +0.1695]**,
  +0.0237 [+0.0192, +0.0279] over RRR in every fold, +0.0278 over ten random genes on the residual, **a dead tie with states + cis
  (+0.0002 [−0.0066, +0.0076])**. Re-fitting the 10 states on RRR's residual HURTS (−0.013; union −0.015 vs offset + cis). The local
  component is additive to the global one when the global prediction is held fixed. Transfer half running (job 21537530); the sealed
  rule (ties within 0.005 → the better-transferring form is primary) will decide between states + cis and RRR-offset + cis.
- ⛔⛔ **Numerical defect caught by a reproduction guard on the transfer set**: the transfer lane's all-99 RRR chose shrinkage 0; the
  SVD solution then weights the training null direction by 1/σ ≈ 1e14 — ~0 in sample, ~1e-2 evaluation-order noise on NEW samples
  (two implementations differed by 7.7e-3 in per-region ρ; states + cis reproduced to 0.0). Every RRR-based number is being recomputed
  with a shrinkage floor of 1e-7 (`chromatin-stable-rrr-forms-*`, job 21537620), which also re-applies the sealed primary rule.
- ✅ **v1.2 SHIPPED** (`chromatin-state-v12-build-20260908T172230Z/build_v12.py`; `release/masld-liver-chromatin-state-v1.2/`,
  `test_release.sh` ALL CHECKS PASSED). Stable lane (`chromatin-stable-rrr-forms-*`, shrinkage floor 1e-7): in cohort states + cis
  +0.1429 vs RRR-offset + cis +0.1413 (tie, −0.0016 [−0.0081, +0.0050]); transfer (reliable set, marginal) +0.205 vs **+0.256**
  (+0.051 [+0.029, +0.071]) → **primary `rrr_offset_cis`** by the sealed rule. RRR-feature + cis transfers 0.006 better still but is
  0.022 worse in cohort. Raw transport reverses the ordering (states + cis +0.175 > offset +0.129 > RRR +0.099). Profiles ship as
  float64 (float32 flipped donor ranks in 2 of 83,627 regions). Measurement-decision and discrepancy numbers on the card are from
  reruns on the shipped form's OOF profile (0.227 vs 0.192 random; D 0.58–1.41, descriptors 14 %). Shipped-form reliable set 4,364
  (states + cis 4,847, overlap 3,146), both in `predictable_regions.bed`. v1.1 carries `SUPERSEDED.md`.

## Binding claim limits, for every model above

- ⛔ **No comparison to FIB-4.** A tissue model needs a biopsy and can never be a screening test.
  SomaSignal 0.90 for F>=3 and LSM-VCTE 0.85 need no biopsy and we do not beat them.
- ⛔ **No performance claim over the incumbent is certifiable**: MDE 0.0667 at K=14 arms on 585
  pooled participants; reaching 0.05 needs ~1,041 participants, a cohort that does not exist.
- ⛔ No stage, no interval, no diagnostic, prognostic or screening claim.
- The only genuinely external published liver-transcriptome staging number is Kamzolas 2026,
  AUC 0.80 (0.72-0.89), n=58. Govaere 2020's 0.93 is in-sample AND includes T2DM as a predictor.

## Where the detail lives

`executions/transportable-substrate-20260904T214259Z/RESULTS.md` — arm-by-arm record and six defects.
`executions/liver-repbench-synthesis/REPORT.md` — the representation x task matrix.
`executions/panel-abundance-audit-20260905T140637Z/RESULTS.md` — panel gene abundance, the floor,
the depth confound, and the abundance-constrained refit. Reproduces the released 0.6451/0.6429/
0.6567 and the 0.9523/0.9739/0.7918 ceilings before making any new claim.
`FINDINGS_2026-09-05.md` — reader-facing summary. ⚠ **Stale**: predates the exposure determination,
the MLM validation, and the Geneformer discrepancy above.

## Development lane 2026-09-07 — RNA carries chromatin variation beyond histology (GSE267145, development only)

`executions/chromatin-beyond-histology-20260907T161521Z`, sealed prespec + addendum 01, seed 20260907.
No weight leaves the directory (GSE267145 derivative-weight terms). Discovery PC5 is a discovery object
(all 99 donors); 26 executions have touched this cohort — development evidence, not validation.
- **Zero-fit (Step 1):** frozen `NAS_head_464` (never saw GSE267145 or chromatin) vs H3K27ac PC5:
  marginal 0.707; partial given steatosis + ballooning + lobular + fibrosis + sex **0.434 [0.181, 0.629]**,
  Freedman–Lane 0/5000, MDE(80%) 0.293; + 14 RNA descriptors 0.354 [0.090, 0.600]; + 7 chromatin
  descriptors + BayesPrism hepatocyte fraction **0.426 [0.163, 0.672]**; + GOCC_LIPID_DROPLET score 0.382.
  Not concentration, not composition, not "finer steatosis" by that proxy. The v1.1 fibrosis score
  carries far less (0.24 given histology). ⛔ PC index 5 is not a stable object: the fixture's own PC5
  correlates 0.33 with the deposited PC5 (neighbouring eigenvalues 4.6 / 3.7 / 2.7 %).
- **Nested CV (Step 2):** in each training fold, after residualising regions on the 4-descriptor
  concentration block, the component most partial-associated with steatosis is index **3** in all five
  folds (EV 0.044–0.056), loading stability mean |cos| **0.711** (pairs 0.585–0.819; fixed index 5 only
  0.353); pooled target vs deposited PC5 0.672, vs steatosis 0.626. RNA ridge (λ by inner CV) predicts it
  out of fold: partial given histology + sex **0.598 [0.437, 0.713]**; given + RNA descriptors
  **0.552 [0.351, 0.704]** (primary); given + chromatin descriptors + hepatocyte fraction 0.579; given +
  lipid score 0.548. In-fold RNA PCs 1–5 + histology reach 0.574 (ridge − PCs +0.02 [−0.16, +0.19]);
  frozen NAS head + histology 0.374 (ridge − head +0.18 [−0.08, +0.45]); descriptor-only model ≈ 0
  (ridge − descriptors +0.55 [0.19, 0.83]). Nested-rank RRR (rank 11–12) beats rank 4 on identical
  folds: skill 0.1194 vs 0.0988, **+0.0206 [0.0136, 0.0276]**. Within steatosis grade 1 (n=43), top-8 vs
  bottom-8 by RNA prediction differ by 1.36 SD of measured target (0/2000 random draws ≥). ⚠ ridge λ hit
  the grid floor (c = 1e-4) in every fold. **Regenerated 500-draw null** (RNA rows permuted against
  chromatin + histology, ridge with inner CV refitted per draw): mean +0.005, sd 0.122, **q97.5 0.229**,
  max 0.445; observed 0.552, **0/500**, MDE(80%) 0.343. **Sealed decision rule → proceed** (stability ✓,
  lower bound > 0 ✓, above null q97.5 ✓, not absorbed by lipid score ✓). Full report: `RESULTS.md` in
  the lane. ⛔ Development evidence: one cohort, discovery-motivated hypothesis, age/BMI unmeasured, no
  second H3K27ac cohort exists. **Exploratory H1**: the RNA correlates of the axis are the SREBP2 sterol-biosynthesis program (Hallmark CHOLESTEROL_HOMEOSTASIS the only set at BH q≤0.05, +0.349; HMGCS1/SQLE/IDI1/MVD/HMGCR/NSDHL/FDFT1 top genes; MYLIP/IDOL opposite); 0/16 BayesPrism fractions; descriptors |rho|≤0.26; no ribosomal/MT genes. Candidate only — statin exposure (unrecorded) induces exactly these genes.

## Development lane 2026-09-07 (evening) — what does NOT beat a linear model on raw genes

Four families were tested against the same GSE267145 target (held-out H3K27ac, five frozen folds,
concentration-residualised log2 CPM), each under its own sealed prespec. **All four lost.**

| family | execution | result |
|---|---|---|
| **Teacher–student (paired training)** | `paired-teacher-student-20260907T183136Z` | ⛔ **−0.0570 [−0.0722, −0.0442]** vs the best direct arm; replicated at −0.0128 [−0.0290, −0.0027] in GSE296875. The pairing control shows the chromatin half IS informative (destroying it costs a further 0.04): averaging it with an RNA-side score the student can already reproduce **dilutes the target**. |
| **MOFA+** | same | ⛔ **−1.7446** on its own 5,000-region panel — 1.74x worse than the training mean — vs +0.1006 for nested RRR on the identical panel. No amortised encoder; a new sample's factors are least-squares-solved against the RNA loadings and that does not transfer. |
| **BulkFormer-37M embeddings** | `rna-representation-benchmark-20260907T192029Z` | ⛔⛔ **The pretrained embedding loses to its own RANDOM-INIT TWIN on all 5 aggregations and all 3 seeds**: pretrained skill +0.009/+0.030/+0.012/+0.020/+0.013 (mean/max/median/all/interested) against twins +0.049 to +0.071. Component partials likewise (+0.19–0.26 vs +0.32–0.40). Pretraining does not merely fail to help on this task; it discards participant variation a random projection of the same architecture keeps. |
| **Corgi, zero-shot AND fine-tuned** | `corgi-gse267145-h3k27ac-comparator-20260907T185333Z`, `corgi-film-h3k27ac-finetune-20260907T191151Z` | ⛔ Zero-shot participant-differential skill **+0.00028** with the real context indistinguishable from a **shuffled** one (+0.017 [−0.018, +0.053]), while region-level Spearman is +0.476 — it knows which regions are active, not who the sample is. ⛔⛔ **Fine-tuning does not fix it**: the ladder's `film_plus_head` rung (5,869,229 trainable params, backbone frozen, 4 epochs, loss monotone on all 5 folds) moved skill by **+0.00016**. Ours on the identical regions: **+0.1186**. Same ordering against raw log2 CPM. |

⭐ **What wins:** ridge/nested-rank RRR on raw log2 CPM genes (whole-matrix skill **+0.1194**), and for
the *component* endpoint an in-fold 81-dim PCA is better still (mean component partial **+0.5164** vs
+0.4765 for 42,163 raw genes; k1 +0.724, k2 +0.792). Tuned kernel ridge on the same PCs is slightly
worse (+0.5039), consistent with every earlier nonlinearity result in this project.

⛔⛔ **A defect of mine, corrected and on the record.** Every "ours vs Corgi" number in the first
scoring run was computed against a **broken reconstruction of our own model** (lane A's per-fold
component predictions multiplied by the ALL-99 PLS loadings — different bases, arbitrary SVD signs).
It scored −0.095 where the real model scores +0.103. Three reported comparisons were retracted and a
sealed prediction recorded as FAILED is in fact satisfied by ~400x. Correction:
`corgi-gse267145-h3k27ac-comparator-20260907T185333Z/CORRECTION_01.md`; the true out-of-fold
reconstruction is deposited at `paired-teacher-student-.../out/lane_a_oof_reconstruction.npz` under a
guard that reproduces the published skills to 1e-15. ⭐ **Rule: a derived artifact must reproduce a
number its producer already published, or abort.**

`corgi-genomewide-h3k27ac-20260907T195320Z` extended the comparison to **all 94,026 covered autosomal regions**
(8 L40S shards, forward strand as declared): ours (nested RRR) **+0.1191**, raw-gene ridge +0.1039, Corgi actual
**+0.0001**, shuffled +0.0002, nearest-wrong +0.0001, cohort-mean exactly flat; ours − Corgi profile ρ +0.286
[0.252, 0.319]; region-level ρ of region means +0.468. All four sealed predictions held; **release claim STANDS**.
`ADDENDUM_01.md` (post hoc, handoff finding 4): training-fold-mean centring with training-fixed scale agrees
(+0.1187 vs −0.0004). ⛔ The shipped rank-10 model's own genome-wide number is +0.045.
⛔ Corgi jobs must request `--gres=gpu:l40s:1`: the project's `rotary_emb` overlay is compiled for
sm_89 only and fails on b6k with "no kernel image is available for execution on the device".

---

## Candidates 2026-09 (not adopted)

Added 2026-09-26; revised 2026-09-29. Candidate models from the 2026-09-23 program; none is
released or adopted, and none has been applied to a Resource number. The A1 identity check shows
that three sets of Resource participant counts are overstated or cannot yet be checked
(`docs/STATUS.md` open defects). Record: `docs/technical/NEW_MODELS_PROGRAM_2026-09-23.md`.
Execution paths are under `executions/`.

| model | question | state | numbers (source) |
|---|---|---|---|
| **A. In-biopsy allelic (Model A-final)** | Is the reference liver cis allelic effect observable in MASLD biopsy RNA, and does it differ in proportion across recorded fibrosis stage (F0-1 / F2 / F3-4, cross-sectional)? | **in progress**. Design chosen 2026-09-26 (κ/τ version); margin Δ = 0.125 per stage step. Gates 0, A, B and C not yet passed; no tag allelic ratio has been read against stage. Claims limited to a cross-sectional, severity-associated difference in tissue-level allelic ratio | A1 identity: 1,466 libraries, 1,349 genetic individuals, including 105 GSE193066 individuals that are provisional until the contamination-aware rerun (`in-biopsy-allelic-a1-full-20260924T025858Z/identity/identity_summary.json`); identity defects in the 844-participant bulk fit, GSE213621 and GSE193066 are `docs/STATUS.md` open defects |
| **B. Noise-aware aspects** | How many separable molecular axes do the four histology components carry once each score is a noisy ordinal readout? | **closed** 2026-09-26 at the SVI result (user decision). Candidate qualifier only: component scores are not read as separate molecular readouts; component attribution `untestable` under PRESPEC_B §7 (minimum detectable effect above the declared limit in GSE281797 for all four components); in paper text the resolved aspect is recorded as indeterminate (`docs/PAPER.md:530`). Held-out likelihood preferred K = 2 over K = 1 (gain positive in 5 of 7 cohorts; one of the five is GSE193066, gain 0.067, where 132 of 164 libraries meet the possibly-mixed rule; 4 of 6 without it); the axes were not characterized (NUTS did not converge, no loadings); no sealed test | K = 2 over K = 1: 0.0280 nats per participant record, p 0.005 (1/201, the floor of 200 draws), 1,007 participant records (`aspect-latent-b3-tests-corr-20260924T043421Z/summary.json`, `observed.json`), at most 975 people by genotype (GSE213621 25, GSE135251 4 and GSE130970 3 repeat libraries; at most 969 if the provisional GSE193066 calls hold; GSE267145 not in the identity panel; `aspect-latent-b2-inputs-20260924T003835Z/b_inputs_corr.tsv.gz` joined to `in-biopsy-allelic-a1-full-20260924T025858Z/identity/`); K = 3 p 1.0, K = 4 p 0.398, rejected (same file); original counts 0.0217, p 0.010 (`aspect-latent-b3-tests-orig-20260926T014659Z/summary.json`); test 3 did not detect component-specific signal, p 0.465 to 0.981, BH q 0.981, with the observed value below the K = 1 null mean for three labels (`aspect-latent-b5-test3-corr-20260926T140056Z/summary.json`); NUTS max R-hat 4.55 (`aspect-latent-b4-final-20260926T135812Z/nuts_K2_summary.json`) |
| **C. H&E to molecular programs** | Can program scores be read from GTEx postmortem liver H&E beyond fat content and recorded postmortem covariates? | **done** (sealed test run once, 2026-09-25). Claim limit: "readable from GTEx postmortem H&E beyond recorded covariates", not a biopsy readout. Image exposure: Phikon-v2 pretraining included GTEx slides (images only; no RNA target seen) | 16/79 programs pass all four criteria (development and sealed, primary and residualized; `he-molecular-c2-20260924T004615Z/c5_sealed_test/readability_four_criteria.tsv`), including an immediate-early (FOS) program, so residualization does not remove all postmortem signal; sealed 45 donors, 61/79 primary, 30/79 residualized (`c5_sealed_test/c5_summary.json`); embeddings predict ischemic time R² 0.62 and Hardy class R² 0.35 in 179 development donors (`c4_run3/c4_summary.json`) |
| **C external test (TCGA-LIHC non-tumor, Amendment 6)** | Do the 16 GTEx-passing programs carry to liver that was not collected postmortem? | **done** (run once, 2026-10-03). Primary endpoint not met. Claim limit: not demonstrated in this cohort; not `tested_negative`. Image exposure: Phikon-v2 pretraining included TCGA slides | 50 donors; median partial Spearman 0.204 (p 0.103), residualized 0.231 (p 0.001); 1/16 programs meets the rule in both versions (Myeloid / immune, DOCK2); hypoxia negative control 0.313 (`he-molecular-c6-tcga-20261003T133755Z/c6_external_test/c6_summary.json`) |
| **F mouse diet models vs human program stage trend (PRESPEC_F)** | Which fixed human programs does each public mouse diet model move in the human stage direction, and is that specific to diet class? | **done** (2026-10-03, candidate, not adopted). Primary endpoint not met: no diet-class specificity. Public mouse data only | 17 contrasts, 16 testable; same-class minus different-class similarity −0.012 (p 0.479); median Spearman with human trend 0.726 (ceiling 0.801); 2/16 at BH q<0.05; 27 of 62 stage-associated human programs untestable in mouse (`preclinical-fidelity-20261003T153626Z/f2_fidelity_r2/f2_summary.json`) |

## Candidates 2026-10 (not adopted)

Added 2026-10-17 by user decision ("yes to all"; `executions/model-program-20261015T223347Z/decisions_pending.md`, last section). Candidate models from model-program-20261015T223347Z, Track B. None is adopted, and none has been applied to a Resource number, figure or the Gene Catalog. P = `executions/model-program-20261015T223347Z`.

| model | question | state | numbers (source) |
|---|---|---|---|
| **B1. Pairwise relative conformal severity calls (installed as `release/masld-severity-v1.3/`, package masld-pairwise-conformal-v1)** | Within one bulk cohort, can a pool-calibrated threshold on a refit severity score give ordering calls for biopsy pairs, abstaining when the gap is small, with coverage near 0.90? | **comparable to masld-severity-v1.2** (2026-10-16, candidate, not adopted; installed as a release candidate 2026-10-17, does not supersede v1.2). All prespecified lines pass. The refit score orders pairs no better than v1.2; what B1 adds is the calibrated abstention threshold. No finite-sample guarantee: q is calibrated on LOCO out-of-fold margins in the 449-unit pool and applied to the full-pool refit. Coverage is empirical and marginal over one discordant within-cohort pair, not P(correct given a call). Coarse classes, within-cohort only, cross-sectional stage | GSE213621 (confirmatory; 257 genetic individuals; F0F1 < F2 < F3F4; 21,986 discordant pairs): coverage 0.9016 (participant-clustered 95% 0.8736-0.9278), call rate 0.7863, call accuracy 0.8749 (clustered lower bound 0.8360); same-class pair call rate 0.6780; 0.80-level coverage 0.8114, but q at 0.80 is about -0.003, so that rule calls nearly every pair. B1 minus v1.2 at matched call rate -0.0106 (95% CI -0.0197 to -0.0021), comparable under the 0.02 margin; re-estimate cohorts +0.0027 (-0.0113 to +0.0168). Re-estimate on seen cohorts: GSE225740 coverage 0.9027 (87 units), GSE281797 0.9186 with call rate 0.5995 (91 units) (`P/trackB/B1/run-20261015T232535Z/results.json`, `P/trackB/B1/RESULT.md`) |
| **B2. Geneformer V2-104M on bulk rank-value input vs random-init twins and linear ridge** | Do frozen Geneformer V2-104M L12 mean-pool embeddings of bulk rank-value input order liver biopsies by fibrosis stage better than identical random-weight encoders, and than a ridge on the 26,629-gene axis? | **no call vs twins (lambda shift); worse than ridge** (2026-10-16, candidate, not adopted). Exposure unknown for all four cohorts, so every number is descriptive. The object is bulk rank-value input to a single-cell-trained encoder; the encoder sees at most 2,046 of the 26,629 axis genes. Bounds only this task. Nothing is packaged | 454 genetic individuals, 4 cohorts, leave-one-cohort-out, 2,048 tokens: S(L) 0.6442; Delta = S(P) minus mean S(twins) -0.08, 97.5% CI [-0.24, +0.24], no call (a one-step lambda shift moves Delta by more than 0.01); P minus L -0.44, 97.5% CI [-0.59, -0.17], worse; MDE(Delta) 0.3077 (`P/trackB/B2/stage2-20261016T014110Z/eval/results.json`). GSE281797 frozen transfer (secondary, n = 94; 93 of 94 outside the v1.1 sparsity contract): Delta +0.18 [-0.04, +0.38], no call (lambda shift); P minus L -0.15 [-0.35, +0.06], inconclusive (`P/trackB/B2/transfer-gse281797-20261016T065452Z/collect/transfer_results.json`, `P/trackB/B2/RESULT.md`) |
| **B3. Donor-mean foundation-model cell embeddings vs linear pseudobulk baseline (snRNA)** | Do donor means of frozen single-cell foundation-model embeddings (TranscriptFormer TF-Sapiens, scGPT whole-human; Geneformer V2-316M descriptive) order donors by cross-sectional stage better than a linear pseudobulk baseline? | **inconclusive in T1; no headline** (2026-10-16, candidate, not adopted). Unit is the biological donor. T1 headline rule not met for either decision arm; "comparable" was unreachable for every delta, so no result reads as "no difference". UCE 33-layer and 4-layer dropped at the forward check (no MASLD result). T2 is a reanalysis of an opened set (GSE289173; second use of its labels), development-grade, every class descriptive. Bounds only donor-mean embeddings with the prespecified PCA + ridge readout | T1 (GSE202379 + GSE244832, 55 evaluation donors, stage Healthy / Steatosis / Steatohepatitis): S baseline 0.5367, TF-Sapiens 0.0859 (delta_base -0.4508, inconclusive), scGPT 0.0791 (delta_base -0.4576, inconclusive); delta_twin TF 0.0827, scGPT -0.1540, both inconclusive (`P/trackB/B3/readout-20261016T052514Z/main/scores.tsv`, `deltas.tsv`). T2 (GSE289173 FL series, 26 donors): S baseline 0.5986, TF-Sapiens 0.1707 (exposure unknown; delta_base -0.4279, inconclusive), scGPT 0.0236 (delta_base -0.5751, worse, descriptive), Geneformer 0.1280 (`P/trackB/B3/t2-readout-20261016T103328Z/main/scores.tsv`, `deltas.tsv`, `P/trackB/B3/RESULT.md`). Label note (2026-10-17, user decision): T1 labels rebuilt from source differ for 2 of 18 GSE244832 donors (Steatosis/Steatohepatitis, one class each) and for 1 GSE202379 donor, with 1 source donor absent from the roster. Under the source labels every T1 class stays inconclusive, but the scGPT delta_base refit decision-level upper bound rises from +0.0146 to +0.0860, so no exclusion of a +0.02 scGPT gain over pseudobulk may be stated (`P/deep-20261017T022651Z/raw/B3_perdonor/run-20261017T035409Z/main/relabel_deltas.tsv`) |
| **B4. Cross-assay H3K27ac transport (rank-4 RRR head vs PCA ridge)** | Does the frozen rank-4 RNA-to-H3K27ac head, trained on GSE267145 CUT&RUN, predict per-donor H3K27ac count composition in GSE128072 ChIP donors better than the frozen PCA ridge recipe? | **worse** (2026-10-15, candidate, not adopted). Development evidence: the 17 receiving donors were pre-exposed by jobs 22026740 and 22026751. Label holds under all four scorings (F12, F9, F6, MASK). Bounds only these two frozen recipes, this deviance scoring and these donors; not a statement about RRR vs ridge as methods. Training and receiving donors are different cohorts (bariatric biopsies vs transplant-donor livers); overlap implausible [Inference], no genotype check | Delta (rank-4 minus PCA skill, F12) −0.068, BCa 95% [−0.113, −0.014], sign-flip p 0.020; n = 17 donors; PCA skill −0.408, rank-4 −0.475, both below the pooled training reference in 17/17 donors; size depends on the floor (F6: −0.111). Autosomal-only Delta −0.054 [−0.101, +0.004] includes 0; receiving 11/17 male vs training 14/99 male. Development D1 (GSE267145, n = 99) −0.032 [−0.048, −0.011]. Reproduced by two independent recomputes (`P/trackB/B4/runs/main_20261015T230900Z_22261752/receiving_results.json`, `P/trackB/B4/RESULT.md`, `P/trackB/B4/review1_20261015T235223Z/recompute.json`, `P/trackB/B4/review2_20261015T235500Z/recompute.json`) |
| **B5. GSE268273 missing-component RNA program prediction (RRR vs ridge)** | Does a reduced-rank map predict 8 held-out program scores from 714 panel inputs better than full ridge? | **comparable, degenerate** (2026-10-15, candidate, not adopted). RRR chose rank 8 (full rank) in 5/5 outer folds and the secondary split, so it is the same model as ridge; the rank question was not tested. Acceptable-error line untestable (Stage 0 power 0.764 < 0.8). Recovery not shown to be program-specific. Spent development cohort | RRR = ridge skill 0.777, 95% [0.705, 0.836]; Delta 0.000 [0, 0], H2 p 1.0; n = 109 participants, 5-fold cross-fit. Technical-only ridge 0.564; panel genes add 0.216 [0.159, 0.317]. Matched gene sets: p_spec 0.178 (17/100 ≥ 0.777, median 0.760). Permutation p 0.001 (floor). Ridge penalty at grid floor 1e-7 in 3/5 folds (`P/trackB/B5/runs/main-20261015T231819Z/results.json`, `selections.tsv`, `P/trackB/B5/runs/stage0-20261015T231819Z/stage0.json` (power 0.764), `P/trackB/B5/RESULT.md`) |
