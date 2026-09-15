# MASLD model results dashboard

**Last refreshed: 2026-09-02 EDT (multi-modal day).** Long-form history: `CHANGELOG_DETAIL.md`.
Frozen evaluator output files are the metric authority; live SLURM is the scheduler authority.
Goal: a lineage-specific, uncertainty-aware map of cross-sectional liver remodeling across
fibrosis, NAS components and a descriptive continuum — not a clinical classifier.

## 2026-09-02 — THE MULTI-MODAL ANSWER

### ✅ What survives its own nulls
| claim | evidence |
|---|---|
| **RNA weights read PROTEOMES** | ρ +0.6128, perm p 1e-4, permuted-weight control −0.004; **survives complexity control at +0.5382 partial** |
| **RNA and protein share a SEVERITY direction** | concordance **0.4523**, **no shared participants**, 3,941 symbols; **label-null audited, 5.8 sd, max 0.3669** |
| **Chromatin carries ordering RNA does not** | partial **+0.3194** vs RNA-given-chromatin **+0.0130** — two lanes, two estimators |
| H3K27ac + RNA both clear complexity on NAS | +0.1925 and +0.2378, both excluding zero |
| Liver proteins add over complexity | **+0.2448 / +0.2139** on *live* gates (⛔ not the dead-gate +0.2413) |
| plasma → NAS over clinical covariates | +0.2791 [+0.118, +0.434] |
| ⛔ **The liver "confound" IS FIBROGENESIS** | 22 ECM / 5 collagens up, hepatocyte enzymes down — **tested, not inferred** |

### ⛔ What fails, with mechanisms
**Fusion gains nothing anywhere** — both modalities over-rank F0 and under-rank every higher stage *in
the same direction on the same participants*; that shared error IS the 0.945 residual correlation.
**A second tissue HURTS** (liver+scWAT −0.1091 [−0.181, −0.051]). **Cross-tissue concordance fails its
correct null** (primary p 0.254, not 0.0005). **MASH/MASL transfer fails the training-label null** in
both arms despite AUROC 0.82. **No peripheral tissue substitutes for biopsy.** scWAT vacuous, oWAT not
evaluable, ATAC bounded/confounded, methylation has no endpoint.

### ⭐ Modality inventory
**ONE externally validated modality** (bulk RNA, 3 cohorts, 585 participants) **+ four single-cohort
demonstrations.** GSE267145 is the **only** cross-modal pairing with an endpoint — and it **forbids
shipping weights**.

---

## 2026-09-01 — SHIPPED, and four things below this line are RETRACTED

⭐⭐ **`release/masld-severity-v1/` exists and a stranger can run it.** Independently verified by an
agent that did not build it (jobs 21328834, 21328840): OOD Spearman **0.663119** reproduced end to end
through `score.py`, filter enforcement bitwise identical, four refusals fire, manifest 10/10.
Shipped arm is **per-sample quantile normalization** — the only construction that is coherent
(same transform at train and inference) *and* runnable on one sample.

| claim | status after 2026-09-01 |
|---|---|
| "reproduces ~4x better than the histology label" | ⛔ **RETRACTED as a general claim.** The label's own ICC is **0.62** in PXD051911, not 0.19. Survives only as a within-sample contrast at a stated stage variance. |
| σ_e = 0.6433 as *the* Kleiner constant | ⛔ **Does not replicate.** PXD051911 gives 0.4655–0.5092, one arm formally `discrepant`. Ship the **range**; 0.6433 is the **anti-conservative** end. |
| "transfers ABOVE native" | ⛔ **"Not worse than native."** Paired delta +0.0509 CI [−0.028, +0.137] includes zero. |
| "predictive intervals are calibrated" | ⛔ Raw coverage **0.404** at nominal 50%. Ordering is the deliverable. |
| Lane D residualized control | ⛔ **Unidentified** — the paired delta reverses sign with the residualization scope (−0.119 pooled, +0.128 within-etiology, shift 0.253 > the effect). |
| Lane C, the non-invasive route | ⛔ **Negative.** Routing plasma→predicted-liver→severity is *worse* than direct; liver 0.699 vs plasma 0.050 on the same 48. |

⭐ **The recurring confound, found in three independent lanes *(⚠ later SEVEN — see the modality push)* today**: library-complexity descriptors
(Shannon entropy, Gini, top-N count shares). ⛔ An earlier attribution to `albumin_fraction` is
**retracted** — dropping albumin costs 0.048, not 0.34. **Put them in the baseline before claiming
biology.**

⭐⭐ **What survives and is the real contribution**: **the instrument transfers to THREE independent
external cohorts, 585 participants** — ρ/ceiling **0.68–0.79** across two platforms and three fibrosis
coding schemes, matching a ridge trained inside each cohort (every paired interval covers zero,
permutation p 0.0001). And **σ_e is far more transportable than ICC** — 1.38x versus 3.2x (pooled A4 arm, 0.6119/0.1886 = 3.244; the A1 arm gives 3.302, which is the 3.3x quoted in the detail section below — same numbers, different arm) across two
deposits at an essentially unchanged denominator.

---

## Status at a glance

**Updated 2026-08-31.** The pattern below is the campaign's main result, and it is consistent across
every lane scored so far. ⚠ **Read the 2026-09-01 block above first — it retracts four claims that
appear below.**

### ✅ What works: translating between two readouts of ONE observed state

| lane | result | n | replication |
|---|---|---:|---|
| **RNA → chromatin** | **+0.13345** (H3K27ac) / **+0.12985** (ATAC), p 0.000 both, 5/5 folds | 99 participants / 39 donors | ⭐ **recurs across two cohorts AND two assays**; ⚠ disease-dependence is **GSE267145 ONLY** — see correction below |
| Zonation from expression | accuracy **0.7895** vs chance 0.333; MAE 0.2281 vs 0.8889 | 19 participants | one cohort only |
| Graded fibrosis beyond composition | C1 **+0.17329** [0.0785, 0.2683] ⚠ *(quoted elsewhere as [0.0694, 0.2771]; `verdict.json` carries the former)* frozen, LOCO ⚠ | 521 samples | 5 graded cohorts |
| NAS transfer | Spearman **0.756**, survives leak hunt | 99 | 1 external cohort |

### ⛔ What fails: predicting the effect of a PERTURBATION from sequence

| lane | result | why it is a real negative |
|---|---|---|
| **Variant → gene eQTL effect** | Pearson **0.8669 → 0.0210** under honest LD-block splitting | leakage was carrying it; only 0.57% of held-out rows share a gene with training |
| Signed effect direction | DNA-LM heads **0.5016** / **0.5015**; allele-feature ridge 0.5033 | at chance vs a **measured** 0.4998, against ceilings to 1.0000 |
| **S-LDSC, 117 programs** | **0 of 117** pass BH on the coefficient, 5 traits | 18 "enrichment hits" include **5 depletions**, one *negative*; largest tau_z +2.788 |
| Enformer on MPRA | **−0.0492** vs an allele-identity ridge, 0/5 seeds | loses to knowing which allele it is |
| Architectural complexity | fusion, routing, multi-task all give **exact zeros** | routed system is bitwise the identity map |
| **Response to intervention** | diet-arm effect **reverses sign in 5/5 readouts**; within-arm p 0.374–0.888 | pooled "tracking" is the treatment arm; surgery contrast rests on **1 participant** |
| **DNA-LM embeddings on variant→gene** | **inert** — permuting the embeddings *raises* the score (+0.2126 → +0.2269) | the whole score is the distance block; embedding-only \|z\| is −0.0001 |

⭐ **A third contrast, from the encoder field test:** pretrained single-cell encoders beat a
50-dimensional PCA of raw counts by **+0.045–0.054 with a linear probe** and only **+0.004–0.017 once
the head has one hidden layer** — and on clean held-out studies that PCA is a **dead tie with a
liver-specific encoder trained here on 102 donors** (+0.0004 [−0.0057, +0.0070]). ⛔ The two encoders
that look best pooled, Geneformer and TranscriptFormer, are **exactly the two with no clean fold**.
**Much of what reads as pretraining value is head capacity plus unaudited exposure.**

⭐ **The thesis this supports:** these models interpolate well **within an observed molecular state**
and do not predict **what an intervention or a sequence change will DO to it**. Two independent
contrasts now carry it:
- RNA predicts a held-out donor's chromatin profile in **two cohorts and two assays**, while sequence
  cannot predict the sign of a cis-eQTL effect **above chance**.
- The same severity models that grade a biopsy **do not move correctly when the patient is treated** —
  the pooled appearance of tracking is the treatment arm, and the one stratum that can separate them
  reverses sign.

⭐ Stated plainly: **grading a state is not the same capability as predicting a change to it**, and
this campaign now has direct evidence of the gap in both the genetic and the therapeutic direction.

⚠ **The honest qualifier, and it is now narrow.** On MPRA — a *reporter* assay — hyenadna beats the
ridge (**+0.0820**, 5/5 seeds) and nucleotide_transformer does too (**+0.0537**, 5/5). But the ridge
on **allele identity alone** already reaches 0.1342 of an absolute 0.2162, so even there the task is
dominated by *which allele it is*, not by what the surrounding sequence does. And on the endogenous
variant→gene task the same two models are **inert**: permuting the embeddings does not hurt the
score. ⭐ **The one place sequence models add anything is a synthetic reporter readout, not the
endogenous effect the Resource cares about.**

### ⚠ Two qualifications on the strongest positive

**1. The disease-dependence is normal-vs-diseased, not graded severity.** The RNA→chromatin map does
not weaken across the severity boundary, it **stops working** — cross-stratum transfer falls to or
below the trivial baseline. But dropping the 24 **normal** livers attenuates the effect by **59%**
(Δ +0.1457 → +0.0595, p 0.24). Within diseased liver only, the residual is positive in all arms and
never resolvable (MDE 0.158–0.201 vs points 0.059–0.120). ⛔ **Graded severity is `untestable`, not
negative.** ⛔ Composition is not separable from the map — LOW is 24 normal + 25 NAFL, HIGH is 3 NAFL
+ 47 NASH — so write *the mapping differs between normal and diseased liver*, never within-cell-type
rewiring. ✅ Sex was chased and excluded: it is a complete confound in the pre-registered split, yet
female-only makes the effect **larger** (+0.2482).

**2. Part of the +0.13345 rides on non-biological structure.** Components 2 and 4 are recognisable
liver biology (high-variance genes, 39/150 and 38/150 in Hotspot programs; component 4 tracks
activity at ρ = −0.615). But components 1 and 3 — including the **largest variance share, 0.165** —
load on low-abundance, low-variance features, **3/150 and 2/150 in any program**, topped by
unannotated ENSG ids, lncRNAs and olfactory receptors.

### ✅ That confound is now resolved — and it narrowed the claim

The signed-eQTL negative was confounded with label noise. **It is now `tested_negative`, not
`untestable`**: per-row sign reliability `Phi(|Beta|/SE)` averages **0.7818**, giving ceilings of
0.7878 (flat) to 1.0000 (at |z|≥5), and **every stratum sits at measured chance with intervals
crossing zero headroom**. The embeddings stay **inert at high confidence** — permuting them is ≥
intact in 9 of 10 cells and *helps* in several with intervals excluding zero. A positive control
confirms the representation works: it recovers substitution class at NT AUC **0.9444**.

⛔ **But a third possibility, which I had not named, bounds the claim.** The eQTL labels are
**marginal associations with no fine-mapping**. Under one causal variant per (gene, LD block), **at
most 1.9–4.3% of held-out rows can be causal**, so the ceiling for *any* sequence model is
**0.5055–0.5211** — inside measurement resolution. Much of the task is asking a sequence model to
predict the direction of **LD tagging**, which is not sequence-determined.

⛔ **"Sequence models cannot predict eQTL direction" is NOT a licensed sentence.** The licensed one
is: *on marginal, un-fine-mapped cis-eQTL associations from this panel, DNA-LM embeddings carry no
directional information, and the design cannot separate a model failure from a target that is not
sequence-determined at the variants tested.* The experiment that separates them is the same head on
**eQTL-fine-mapped** variants — the GWAS credible sets supplied the variant list, but the eQTL side
was never fine-mapped.

### ⛔ Three numbers quoted in this file have NO producing file, or the wrong one

Found 2026-08-31 while building Figure 7, by an agent that refused to place any number it could not
trace. **Rule 1 of the plan is that every number in the manuscript needs a producing file on disk**;
these did not meet it, and the summary prose was the only place they existed.

1. **Harmonized C1 `+0.2333 [0.1187, 0.3450]` has no producing file.** The frozen ladder's
   `verdict.json` carries **C1 = +0.17329 [0.0785, 0.2683]**. The status table above has been
   corrected to the traceable value. ⛔ Anywhere `+0.233` appears in older prose in this file, it is
   **unsourced** and must not be carried into the manuscript without locating its producer.
2. **GSE83452 "diet-only MDE 0.138–0.215" is not in the results file.** The producing `results.json`
   has only `P2_mde_80pct` = **0.1001–0.1766**, with **no arm qualifier**. The intervention verdict
   does not depend on it — the sign reversal (surgery ρ +0.5029 vs diet −0.3663) is unambiguous in
   the file — but the quoted range was invented precision.
3. **Multi-head arm `−0.001160`** is absent from `routed_multihead_eval/results.json`.

⚠ Also corrected: the tie-structure ceilings are **0.7918** (GSE267145) and **0.9558** (GSE130970),
and in `two_axis_activity_fibrosis_external_activation.json` the `graded_primary_*` keys are
**GSE130970** while `omitted_level_arm_*` is **GSE193066** — a first extraction that keyed on the key
name rather than each entry's own `endpoint` string had these swapped.

⭐ **The lesson: a summary file is a report, not a source.** Numbers that live only in prose drift,
acquire false precision, and cannot be re-derived. Figure 7's build reads all 27 producing files
directly and records their sha256 in `tables/PROVENANCE.tsv`.

### Standing limits

| | |
|---|---|
| Best confirmed external transfer | NAS Spearman **0.756** |
| Project-reserved model | none; GSE289173 outcomes unopened |
| Unified architecture | **not justified** — specialist complementarity never established |
| Weight release | GSE267145 / GSE296875 runs need written permission; scVI bundle is **atlas-only** |

**Three earlier check failures, three different causes** — none was a method failure:

| Lane | n | Check | Cause |
|---|---:|---|---|
| GSE49541 fibrosis | 72 | FAIL 0/4 | Source had **4** advanced cases; every design linearly separable |
| GSE83452 NASH | 148 | FAIL 3/4 | Paired bootstrap bound −0.047; the one PASS turns on **1 participant** |
| NAS component sum | 99 | FAIL 1/3 applicable | Transfer **works**; margin between two good models is only +0.021 |

## Model-building day — 2026-08-30

Five lanes launched against the approved plan. First results below; two campaigns still running.

### ✅ Liver zone prediction PASSES — the pipeline check, not a discovery

`executions/gse105127-zone-prediction-design-20260830/`, job 21247207, 61 s.
19 participants × 3 laser-microdissected zones, participant-held-out, folds settled without outcomes.

| metric | value | 95% CI | chance |
|---|---:|---|---:|
| **mean absolute zone error** (primary) | **0.2281** | 0.1228–0.3509 | 0.9181 |
| accuracy (secondary) | 0.7895 (45/57) | 0.7018–0.8772 | 0.3333 |

Confusion matches the shape predicted *before* the run: 12 errors, **11 one zone away, exactly 1
inverting the axis**. Ends read cleanly (CV 18/19, PP 16/19); the middle hesitates (IZ 11/19).
Pre-registered band was 0.67–0.78 predicting MAE 0.2222.

⛔ **This licenses "the GSE105127 RNA processing and join are sound" and nothing biological.**
Zonation is among the strongest gradients in the liver and the source study selected sections where
it is visible. Close to a positive control, and it behaved like one. A failure would have been more
informative, pointing upstream at the quantification.

⛔ **A p-value was withdrawn.** A constrained mode exploiting the one-section-per-zone design returned
57/57, MAE 0.0000, p 1.6e-15. Controls killed it: **replacing every probability with 1/3 also scores
1.0000**, because rows are ordered CV/IZ/PP in all 19 participants, so a solver breaking ties by
identity reproduces the truth knowing nothing. Two controls survive and support a narrower claim with
no p-value: shuffling row order within participant still gives 1.0000 over 200 repeats, and the
worst-case assignment margin is 1.71 log units.

⛔ Only zonation cohort in the project. A negative here is **`indeterminate`, never `tested_negative`** —
recorded before any fit existed.

### ⭐ Variant → program substrate built; the LD split is the load-bearing number

`executions/variant_program_substrate-20260830T000000Z/`, 277 files.

- **Max cross-split r² 0.00992** (379-sample panel) and 0.00976 (489-sample), **zero** cross-fold pairs
  ≥0.01 across 8.8 M and 8.5 M examined. Calibrated against each panel's own noise floor: pairs >4.5 Mb
  apart give mean r² 0.00266, matching 1/(n−1)=0.00265. Residual linkage is **below the level either
  panel can distinguish from sampling noise**.
- ⭐ **Per-gene mean beta scores Pearson 0.902 with the block left in and 0.059 held out.** A
  variant-level split would have inherited that 0.902 and none of it real.
- Join rate **3,366/3,721 program genes (90.5%)**, 6,593/7,093 membership rows (93.0%), all by exact
  HGNC symbol; zero rescued by alias, GENCODE v19/v49, or ENSG routes.
- ⛔ **Sign is unpredicted by anything available.** Distance → |z| Spearman 0.174; distance → signed beta
  **0.010**. Per-gene mean beta held out 0.059. All apparent signal is magnitude.
- Nearest-TSS picks the largest-effect gene in a hub **23.1%** of the time vs a within-hub permuted null
  of 7.0% ± 0.16. That is the number a gene-assignment model must beat.
- ⛔ **Ancestry loss has a mechanism**: 10,751 of 41,428 credible-set variants absent from the panel —
  AFR **42.8%**, EAS 26.9%, AMR 16.6%, EUR 9.0%. **Zero surviving variants have European MAF <1%**, so
  the panel carries a European MAF floor and drops exactly what is common in African-ancestry populations.
- ⛔ Power is set by blocks, not variants: 263 genome-wide, median 14 per program, **74 of 117 below 20**.

### ⛔ GSE83452 paired intervention — a well-designed NULL on the only human intervention data

`executions/gse83452-paired-intervention-design-20260830/`, jobs 21247259 and 21247284.
60 participants biopsied twice a year apart, **32 diet / 28 bariatric surgery**, all 120 arrays present.
Unit is the participant; the measurement is the within-participant delta on the fixed-projection axis,
nothing refit on this cohort. 130 of 139 loading genes present.

| | value |
|---|---:|
| mean delta, surgery (n=28) | −0.5971 |
| mean delta, diet (n=32) | −0.1060 |
| **arm difference** | **−0.4911** |
| bootstrap 95% interval | −1.1491 to +0.1938 |
| permutation p (10,000) | 0.1603 |
| delta sd across all 60 | 1.3261 |
| pre-registered floor (0.5 sd) | 0.6631 — **not cleared** |

Direction as pre-registered; magnitude indistinguishable from zero. 21 of 28 surgery participants fell
against 14 of 32 on diet.

⛔ **Regression to the mean INFLATES this comparison — the opposite of what the prespecification argued.**
Baseline axis score is −1.852 for surgery against −2.292 for diet, so surgery starts more severe with
further to fall, and delta correlates with baseline at **r = −0.58**. Post-hoc stratified permutation
(labelled as such): p 0.275 within baseline strata, residual arm difference **−0.1152**, residual
permutation p 0.696. **Only 23.5% of the raw difference survives adjustment**; roughly three quarters
travels with starting position. Surgery is over-represented in the most severe stratum, 9 vs 6.

⛔ **A control fired and the primary was withdrawn by the pre-registered rule — but nothing is hidden.**
The timepoint-scramble check flips each participant's delta, which destroys the effect itself, so it is a
sign-flip significance test of the quantity under study rather than an independent artifact probe. It
fires whenever the effect is not significant, which the primary already reports. The withdrawal is
correct by the rule and is **not** evidence of an artifact; all three statistics agree.
✅ Check A is what a control should look like: observed −0.4911 against a random-gene-set null of
+0.0078, sd 0.1065, **0 of 1,000 as extreme**.

⚠ **Informed secondary, NOT pre-registered** — the transition table by arm was seen before the design was
fixed, and the agent disclosed it unprompted. n=54, binary transition only, since no NAS or fibrosis
stage exists for this series: resolvers −0.783 (n=20), persisters −0.001 (n=17).

⛔ **No NAS and no fibrosis stage exist for GSE83452.** The strong question — does the readout track
*graded* histological change — cannot be asked on this cohort at all.

**Transferable lesson, in the agent's words:** *check a confound on the axis you are measuring, not on a
correlated label.* The binary histology said diet was sicker; the continuous readout says surgery is.
The two imbalances run in opposite directions, so a confound argument built on the label pointed the
wrong way with full confidence.

### ⛔ Variant → program: the gene-level model does not beat distance, and the program layer is closed

`executions/variant_program_substrate-20260830T000000Z/`, README Part 2. Fit on folds 0–3, evaluated once
on fold 4; intervals bootstrap the 54 linkage blocks in fold 4.

| target | distance | model | model − distance |
|---|---:|---:|---|
| absolute z | 0.241 | 0.195 | **−0.044** [−0.085, +0.003] |
| absolute beta | 0.237 | 0.277 | +0.050 [−0.105, +0.203] |

Hub top-1 is a tie: model 0.2163 vs nearest-TSS 0.2138 against a permuted null of 0.0654. **Sign is not
predicted**: Spearman 0.0027 and AUC **0.493** — below chance — on 172,947 held-out rows. Tail results
(AUC 0.599 at |z|>4, 0.627 at p<5e-8) have intervals spanning 0.5 and rest on 32–35 blocks; the block
bootstrap killed them.

⛔ **Splitting by LD block is very nearly splitting by GENE.** Only 164 of 2,587 held-out genes appear in
training; **3.0% of held-out rows** belong to a gene with any training data; 142 of 7,220 training genes
span more than one block. Every gene-level summary feature is unavailable at prediction time.

⛔ **Program aggregation: NEGATIVE, and structurally out of scope.** Same-program gene pairs correlate no
more than different-program pairs at matched genomic separation: observed +0.0165, shuffled −0.0009 ±
0.0393, **z 0.44, p 0.657**, vacuity check z 1.28 (null usable). The specified variance-ratio test came
back degenerate at 1.0002 and its own vacuity check caught it — permuted-within-hub data still beat
shuffled at p 0.02, which is the sparsity pattern, not signal. Statistic reported, **not used**.
⭐ **The reason the lane closes is not the p-value.** 453,266 same-program gene pairs exist and **609
(0.13%) are observable in cis at all**. Whether a program responds coherently to a variant is not a
question cis-eQTL data can answer for 99.87% of the pairs. That needs trans effects or a perturbation
readout; no modelling or extra variants fix it.

⛔ **Power was computed before any null ran**: median 14 blocks per program, **24 of 117 named untestable**
(<10 blocks), 95 of 117 cannot detect below 0.5. Only the pooled test over 263 blocks has power, MDE 0.172.

**The bar for the sequence lane, recorded before it runs: 0.172.** Beating a gene-level sign model at
0.003 is not the test — a sequence result of 0.05 would look like a sixteen-fold gain and be
indistinguishable from nothing.

### ⭐ The 117 programs are effectively 112 independent tests

From eigenvalues of the annotation correlation matrix over all 9,997,231 reference SNPs: Li–Ji **112.00**,
Cheverud–Nyholt 116.78, mean pairwise correlation 0.0112, median −0.0031. Despite 1.91 programs per gene,
once 100 kb windows are laid over the genome the annotations are nearly uncorrelated. **No case for a
lighter multiple-testing correction.**

### ⛔ No MASLD endpoint can carry the heritability question

S-LDSC needs h² z >7 and attenuation <0.2. Only **ukbb_alt (16.4, 0.088)** and **ukbb_ast (12.3, 0.103)**
clear it. pdff 6.35 and ghouse_cirr 6.29 fall just short. finngen_nafld 4.36 / 0.275, ghodsian_nafld
2.83, finngen_cirr 2.33 all excluded; both MVP enzyme arms show attenuation **0.376/0.377**, ~38% of mean
signal attributed to confounding. **The analysis answers a question about serum ALT and AST, not MASLD**,
and that sentence travels in three places including a column on every data row.

### ⛔ Putting baselines on B6K strands the L40S candidate

`config/tasks/cell_state_mapping.toml` lists **scvi_baseline and scanvi_baseline as mandatory baselines**,
not candidates; TF-Sapiens is the candidate. Its claim needs gain ≥0.02 over the strongest baseline,
currently hvg_pca_knn 0.90748, which TF-Sapiens clears by 0.0314. **If scVI or scANVI lands above 0.91887
the margin closes** — a verdict decided by numbers formally incomparable to it (TF-Sapiens is L40S/cu124;
cu124 has no sm_120 SASS, so that arm cannot execute on B6K at all).

⚠ **Runtime change dominates architecture change ~19×**: version effect cu124→cu128 max_rel_difference
**0.566**, architecture effect L40S→B6K **0.0305**. So re-running to escape a pool difference perturbs the
numbers by more than the pool difference. Any cu128 re-run is a **new evaluation**, never a relocation.

✅ **Recorded decision**: the five HVG/PCA baselines stay on their Intel CPU numbers. Re-running them on
B6K would cross Intel→AMD (a boundary this cluster already documented via SCAN/AVX-512), hold a B6K GPU
idle 8.5 h to run sklearn, and move the strongest baseline for silicon reasons. They have zero CUDA
exposure, which is the entire mechanism the comparator rule guards.

⛔ `celltypist_baseline` is not merely unscored, it is **undesigned** — `train_locally` with normalization,
solver, regularization, class balancing, seed, majority-vote policy and calibration all unfrozen. The
eight-member baseline set is 5 done, 2 in flight, 1 not designed.

### ⛔ MPRA claim narrowed twice

The `variant_to_regulation` numbers are real — hyenadna delta_ridge 0.21622, **gain +0.08201** over the
allele-identity ridge 0.13422, 5/5 seeds; nucleotide_transformer 0.18792, +0.05371. These come from
`executions/model-cpu-train-605-21099008` and must **not** be pooled with
`config/dna_lm_family_reconciliation.json`, a different task (`gse281364_mpra_development_sequence_transfer`),
different heads (linear_ridge/two_layer_gelu), 1,000 bootstrap, baseline 0.17708, all gains negative.

⛔ But that run's own receipt records `mandatory_task_native_baselines_complete: false` with five missing,
including MPRALegNet. **Defensible claim: beats the allele-identity control in an incomplete baseline set.**
Not "beats the task-native baselines."
⛔ And **"5/5 seeds" means five candidate refits each beat one fixed reference** — the control carries
`seed_behavior = deterministic_repeated_for_schema_only`. Not five-versus-five.
⛔ hyenadna at CI low 0.0274 is the only robust arm; nucleotide_transformer's 0.0011 touches zero and it is
flagged `restricted`.

⛔ **MPRALegNet, deltaSVM and gkm-SVM are terminally blocked**, so the mandatory baseline set on that task
**cannot be completed by running anything**. Only sequence_CNN and sequence_transformer remain runnable and
their authorisation flags are all false, awaiting the PI.

### ⛔ Four guards found incapable of catching their own target

A recurring defect worth naming, alongside last week's six unusable nulls:
1. A `DatasetManifest` widened to accept arbitrary tables cannot detect an unvalidated field.
2. A 50% coverage floor cannot detect a reference that lost half a chromosome.
3. An artifact sweep reading only `{path, sha256, size_bytes}` tables reported **154 intact, 0 broken**
   while the registry was still failing on standalone `*_sha256` keys it never read.
4. A label-permutation null that a **uniform-probability model passes at p 1.6e-15**.
Two of the four were proposed by me; the fourth got as far as producing a number.

## Graded-fibrosis instrument — 2026-08-28/29

Substrate `executions/pooled_graded_fibrosis_substrate-20260828T112808Z`; primary
`executions/pooled_graded_fibrosis_ladder-20260828T114911Z` (held-back prespec, reproduced
bit-for-bit before any sensitivity arm). **521 bulk RNA SAMPLES**, 5 cohorts each carrying a full
graded Kleiner 0-4 scale, leave-one-cohort-out. `donor_key_exists=false`: the sequencing sample is
the unit and no participant count exists. Every number below was independently re-derived by a
verifier agent that did not produce it.

### The headline: expression carries graded fibrosis beyond composition

| rung | sees | frozen | harmonized |
|---|---|---:|---:|
| a | hepatocyte fraction alone (1 df) | 0.3869 | 0.3771 |
| b | eligible-subcomposition CLR (6 df) | 0.2786 | 0.2917 |
| c | composition + lineage-residualized expression | 0.5602 | 0.6103 |
| d | plain bulk expression, matched support | 0.6879 | 0.6951 |

**C1 = c − a = +0.1733 [0.0694, 0.2771] frozen, +0.2333 [0.1187, 0.3450] harmonized**
(paired change +0.0600 [+0.0147, +0.1076]). Partialling hepatocyte rank out of both rung d and
stage leaves 0.6219; library size leaves 0.6825; both together 0.6128. Tie ceiling ~0.965, so 0.688
is not a ceiling output file. Measured MDE for C1 is ~0.136.

C2 = c − d = −0.1278 frozen → −0.0848 harmonized, **below its 0.09304 floor**: the verdict drops
V4 → V3, so only *plain bulk is not worse* survives, not *plain bulk is better*. C3 dies entirely
(interval covers zero) and **changes sign** under reference-supported lineage families.

### Diffuseness — NARROWED 2026-08-29, the first statement of this was overstated

**Correct claim (supported):** *once a model is fitted*, feature identity buys almost nothing,
because the fibrosis signal is recoverable from very many feature sets.
- **20,000 random H3K27ac regions score 0.7869** against a variance-selected set's 0.7945;
  **20 of 50 draws at or above observed**. 100 random regions already reach 0.683.
- Supervised expansion 139 → 39,798 genes buys **+0.0201** (p 0.059).
- Sharper form: a random 139-gene set's PC1 correlates **0.981** with the global PC1 of all
  39,798 genes, and that global axis alone reaches macro |rho| **0.1587**. **About a third of the
  best published signature's association is available from any 139 genes.**

⛔ **NOT supported — do not write "no signature beats matched chance".** Under a pooled five-cohort
*unsupervised membership* test with a permutation-calibrated null
(`executions/published_masld_signature_panel-20260829T122548Z`), **14 of 19 clean published or
curated sets beat a size- and expression-matched null at BH q<=0.05**, Kamzolas among them
(**0.5535** vs null median 0.1588, q95 0.1786, p 0.0005, sign-consistent in 5/5 cohorts).

⚠ **Two different tests must be reconciled before either enters a figure.** The earlier
"Holm p 0.0819 in both decisive cohorts" is a *per-cohort* test on *two* cohorts using the held-back
fixed-projection loadings; the above is a *pooled five-cohort* membership test. They are not the
same claim and were briefly conflated here.

**What passing this test actually means:** it indexes **fibrogenesis, not MASLD specificity**.
Hallmark EMT (0.438), NABA core matrisome (0.401), Reactome ECM organisation (0.391) and collagen
formation (0.378) — none built for liver disease — beat the null about as convincingly as the
purpose-built signatures. ⛔ And the WikiPathways set literally named *non-alcoholic fatty liver
disease* (152 genes) is the **worst arm in the panel**: 0.1089 vs a null median of 0.1031,
**p 0.443**, indistinguishable from matched random genes. Curation naming a disease is not evidence
it indexes it.

⛔ **Expression×variance matching does not achieve exchangeability.** For **22 of 41 arms zero of
1000** matched random sets reached the signature's PC1 explained-variance fraction (kamzolas145
0.540 vs null median 0.237). Published signatures are co-expressed modules; matched random sets are
not, which is what broke the uncalibrated test.

### The axis is readable in chromatin (development-grade)

GSE267145, 99 samples, participant-held-out, target = frozen out-of-fold RNA prediction of NAS.
H3K27ac recovers the RNA ordering at **rho 0.7945 [0.690, 0.867]**; **0.5552 [0.392, 0.690]**
survives residualizing the binary disease call; the two modality positions still agree at
**0.4867 [0.307, 0.632] after conditioning on the true NAS score**. Floors (i) and (iv) NOT cleared.
`champion_eligible=false`, single cohort.

### It ranks a batch; it cannot place a patient

`axis_raw` is **not** a shared coordinate: stage-matched samples sit **1.386–1.874 Kleiner stages
apart** (p 5.0e-05), 35–47% of the axis's whole 0→4 range. Cause in code —
`02_build_unsupervised_axes.R:172-197` applies `limma::normalizeQuantiles()` **per cohort** while
the loading and centre are shared. Pure shift, not scale; case mix explains 3.4%.
**Label-free fix works**: subtract each batch's own q10 → out-of-cohort offset 0.862 → **0.151**
stages, and pooled Spearman *rises* 0.541 → 0.597. But anchor RMSE is **2.848 stages at n=1**,
0.352 at n=20, 0.217 at n=40; a severe-only batch mis-anchors by up to 2.141 stages. The released
pipeline deposits no frozen quantile reference, so it cannot place a new sample at all.

### Uncertainty: ranks well, resolves stage poorly

First conformal layer in the tree. Nominal-90% half-width **±1.478 Kleiner units** — 69% of the
0-4 scale, failing to exclude **2.94 of 5 stages**. At matched coverage only **7.8% narrower** than
the trivial marginal interval; rungs a and c are *wider* than it. Adaptivity absent (width CV 0.023).
Conditional coverage by predicted stage 0.861 / 0.773 / 0.684. Meets the claimable CV+ level
(1−2α) in 5/5 cohorts, nominal in 3/5.

### The continuum is redundant inside the model, and wins below n=160

Added as one column at matched support its coefficient is **exactly zero in all 5 folds**
(delta +0.000000). Not weakness — **redundancy**: it ranks 2nd–4th of 39,799 columns on the
screening statistic, enters at support 2–5, and is dropped by support 20–27. An elastic net
reconstructs it from the held-out cohort at **rho 0.990–0.996**.
**Cold-start crossover n = 160** labelled training samples (per-cohort 80/80/80/160); at n = 0 the
supervised model is undefined, not zero. Cross-cohort sd is *worse* for the continuum (0.109 vs
0.060) and is **structured** — best where the supervised model is worst.

### GSE240729 resolved, and the two signature tests reconciled — 2026-08-29

**GSE240729 is not a bad cohort; it has the cleanest fibrosis gradient in the pool.**
Five of six anomalies reduce to one operand: C1 = rho(c) − rho(a), and there rung a (composition)
is **+7.45 SD** above the other four while rung c is −0.95 and **rung d is −0.22, entirely ordinary**
(0.6757 in a pool range 0.6355–0.7919). C1 is negative because **composition over-performs**.
Carrier is the **Fibroblast share** (rho 0.711 with stage, pool max vs 0.215–0.474); hepatocyte
fraction is its closure mirror (rho −0.928; partial rho(hep, stage | fibroblast) −0.6415 → **+0.070**).
Not deconvolution: every deconvolution-FREE probe is also at pool maximum there (stellate 0.696,
Reactome collagen 0.679, NABA 0.638, Moylan 0.695, Kamzolas 0.705).
It **is** a damaged-RNA assay, perfectly separated on read-length-free metrics — mapped/input length
**0.9226 vs 0.9973–1.0070** (AUC 0.000), deletion rate **0.13% vs 0.01–0.02%** (AUC 1.000), mismatch
0.555% vs 0.16–0.35%, gene-assigned reads **0.203 vs 0.465–0.788** — but the composition signal
retains **80.2%** after conditioning on nine technical metrics.
⛔ **FFPE and cohort are perfectly confounded**; the only other FFPE liver set in the tree has no
deposited staging. **Action: retain in expression arms; report every composition-dependent contrast
(C1, rung a) as a declared stratum with a drop-cohort arm.**

**The two signature tests are reconciled: the flip is the NULL, not the score.**
Swapping only the null loading source (foreign 135-donor discovery PCA → PC1 on the target) moves
p **0.0430/0.0410 → 0.0010/0.0010**; swapping only the observed score does **not** flip it
(→0.0500/0.0719). Any 139 genes pushed through the discovery PCA already reach |rho| ~0.35 with
fibrosis before a sign is chosen, versus 0.159 for a random set's own PC1 on our substrate — so
test A benchmarks against a baseline already handed most of the answer by the discovery covariance.
⭐ **A de novo 139-gene set selected on our own data also fails test A** (p 0.1009): that null is hard
for every gene identity, not just the inherited one. Aggregation is a second independent factor —
keeping A's null but using one pooled macro statistic gives p 0.003.
✅ **Test B is NOT a coherence test**: within arms, PC1 explained variance correlates **negatively**
with the B statistic (median −0.153, 40/41 arms), and across arms the raw +0.435 becomes **−0.598**
after conditioning on member-gene fibrosis association. Kamzolas 145 ranks **1st of 19** clean arms.
⚠ **Calibration asymmetry**: B's own validity block disqualifies its raw p; A was never calibrated.

**q10 anchor prespecification held-back** (`eb4471b9…`), with clause J5 stating plainly that 0.10 was
selected on performance and specifying out-of-sample confirmation instead. **No label-free case-mix
guard exists** — four shift-invariant shape statistics caught 0/10, 2/10, 2/10, 2/10 severe
truncations — so the guard is label-based. **n ≥ 20 minimum** (RMSE ≤ 0.459 stages worst-case),
**hard refusal at n ≤ 9** (structural: at type-7 interpolation the 0.10 quantile IS the sample
minimum when n ≤ 10). A new 20-sample batch carries **0.588 Kleiner stages** combined uncertainty;
⛔ quoting the 0.151 LOCO figure for a new small batch is prohibited.

⛔ **UNVERIFIED, do not propagate**: the GSE240729 continuum figures repeated in earlier notes
(anchor rho 0.7100 clearing both matched nulls, 0.7100 → 0.4231 under partialling, "unsupervised
PC1 0.715 → supervised 0.587", offset +2.878, anchor failure +2.141) could **not be located in the
executions tree**; an independent plain within-cohort SVD PC1 there gives **0.105**. Different
construction; must not be conflated or quoted until a producing output file is found.

### The routed multi-modal system is a no-op, and the oracle forecloses the whole rule class — 2026-08-29

Held-back prespec `executions/routed_multihead_prespec-20260829T214500Z` (`c3402418…`);
evaluation `executions/routed_multihead_eval-20260829T230000Z`. 521 bulk RNA SAMPLES,
5 cohorts, leave-one-cohort-out.

**The routing rule selects the reference component in 5 of 5 folds, and the routed
prediction vector is bitwise identical to the supervised elastic net on all 521 samples**
(delta_R = 0.0000000000). This is the third exact zero in the campaign.

⭐ **The decisive number is the oracle.** A router with perfect post-hoc knowledge of which
component wins on each sample gains **+0.0071360**, which is **13.0× below the pre-declared
floor of 0.09303604**. An oracle is an upper bound over every attainable routing rule, so
this forecloses the class rather than reporting one failed attempt.

⛔ **The mechanism is structural, not a tuning failure.** The only component that beats the
reference anywhere is the zero-df continuum in GSE240729 (0.7287 vs 0.7073), and **all four
training cohorts prefer the supervised model by 0.11–0.12** — so any rule fitted on training
data must choose wrongly in the one cohort where choosing differently would help.

⛔ **The multi-prediction arm died on coverage, before any modelling.** Fibrosis 521/521
across 5 cohorts; NAS total 250/521 in 2 cohorts, never both multi-cohort-trained and
out-of-cohort evaluable; NAS components 78/521 in 1 cohort; **disease status 0 informative
samples**. Label arithmetic recurred: NAS − (steatosis + lobular + ballooning) = 0 on 78/78,
and `GSE135251.Stage` is exactly I(fibrosis ≥ 3).

The multi-head arm returned **−0.001160** with coupling μ* selected at 0 in 2 of 5 folds.

⛔ **Sixth unusable null.** Uniform-random routing places the all-reference solution at its
99th percentile, because averaging components of unequal quality makes the best one extreme.
⛔ **Disclosed deviation:** 1 guard of 82 failed on a ~9.8e-8 arithmetic slip in the prespec's
own oracle constant, and the run continued instead of raising.

### Corrections issued this cycle

- **"The continuum does not detect disease" — RETRACTED and reversed.** It rested entirely on
  PRJNA512027, retired 2026-05-15 for an L0/S0 library-prep × disease confound
  (all 34 controls L-batch, all 116 inflammation/fibrosis S-batch, rho 0.6153). The axis detects at
  **AUC 0.8638 / 0.8929 / 0.8239 / 0.8130** in four cohorts whose scores were already frozen.
  The composition "winner" there was reading batch (CLR vs batch 0.72–0.77, vs disease within
  batch 0.000).
- **The deconvolution split was a conda-env accident.** `12_run_bayesprism.R` branches on
  `requireNamespace("InstaPrism")`; two cohorts ran under `rnaseq` and got the reference-updated
  **final** posterior, three under `music_deconv` and got the **initial** posterior.
- **The frozen 7-lineage eligibility family is an artifact** of deriving the rule on 3 of 5 cohorts;
  applied honestly to all five it yields 3. pDCs were kept on a **134-cell** reference profile while
  Fibroblasts and Endothelial were excluded.
- **Two nulls specified for this campaign were unusable** and were caught by producing agents rather
  than reported as passes: one vacuous (the known-failing set cleared it 5/5), one a point mass
  (all 250 permutations returned exactly the observed value). Any new null must be shown to reject a
  known failure before a pass from it is trusted.

## Overnight results — 2026-08-28

**Four calibrated negatives, one recovered output file. Every direction the network model
needed is now tested.**

### Lane 1 — paired contrasts on the frozen GSE267145 resamples ✅ closed
⛔ **The modality split does not survive calibration.** Fibrosis rank 4/4
sign-consistent but **0/4 exclude zero**; the **activity half is absent** (−0.012,
mixed signs); one fibrosis endpoint reverses. **The only clean separation is a
zero-inflation shrinkage output file** — fibrosis is 71/99 zeros, and chromatin wins
ordinal MAE by sitting nearer the mode, the same mechanism that lets a constant
baseline beat every fitted model. **Do not write the modality split as a finding.**
Fusion is **F2 not F1** — no evidence either way; late fusion is indistinguishable
from its best parent on **all nine** metrics. The real result is narrower and was
declared in advance: **the entire penalty lives in one prediction head**
(ordinal MAE −0.468 [−0.590,−0.343] against regression MAE +0.003).
`model-check-1102-21183907` · counts unordered: M 6/36, MM 1/9, F 9/54, B 58/70.

### Lane 2 — program coherence along the continuum ✅ closed
⛔ **The 43-system graph cannot rewire** — every edge is IDF-weighted gene-membership
overlap; no sample enters its construction. Adjacent nodes share a median Jaccard of
0.108 (8 genes), so donor-correlation edges would be mechanically confounded; the 116
programs are **79.4% pairwise gene-disjoint**, which is what licenses the program-level
object instead.

**Non-transferability negative.** The two **pre-registered** focal programs show a
real, non-artifactual coherence increase in **GSE162694 only**. **GSE213621 (n=361)
and GSE135251 (n=214), the two largest, are silent on both.** Declared primary pools
to **−0.0012 ± 0.0196**.

✅ **Two mechanisms excluded by measurement, not argument.** The rank-one deflation
cannot produce a tercile contrast — `Cov(x|t)` is **constant in t**, so it lowers
coherence identically in every band; measured floor ≤0.005 across all 40 cells,
against a prediction recorded before the output was read and then falsified. And
**linear** differential expression does not manufacture a contrast either: the
rotation null preserves `Cov(x,t)` exactly, genes' means vary with t as much as in
real data, and the contrast is still zero. Only the **nonlinear** channel remains.

⛔ Two nulls, not interchangeable: split-permutation is inferential (stricter);
rotation bounds the deflation channel only (Gaussianises marginals). Where they
diverged, the wider was taken.

✅ **CLOSED — the dissenting cohort is not attributable to composition-independent signal.** GSE162694
residualized, declared primary: **ductular observed +0.2263 against a
composition-conditional null mean of +0.2206, z = +0.13, p = 0.90.** The one cell
that survived every other null is produced by composition. Stromal is ~half explained
with a nominal residual (1 of 40 tests). ⛔ Phrase it as *"conditioning on deconvolved
lineage proportions reproduces the effect"*, not *"composition causes it"* — the
proportions come from the same matrix `t` is projected from, so the test over-removes
and is conservative. `hac-coherence-{mvp,cal,n1,n2}-20260828*`.

### Also closed
- **CCC**: per-donor LIANA already run correctly — **0 of 13 LR pairs** pass donor FDR.
- **SCENIC**: five cell types, all n=18, **0 significant in every one**.

### ✅ Recovered
**GSE202379 steatosis and activity grades**, frozen additively —
`gse202379-saf-recovery-21183702`. ⛔ Steatosis has **four** levels at n=40.
One ordering bug in the fibrosis harvest with **eight** victims: 7 invented stages,
1 overwritten real grade (P98), 1 agreeing by luck (P30). Blast radius traced and
**contained** — never reached the canonical donor metadata or the 117 frozen programs.


## Network-model redirect (2026-08-28)

The target changed twice, both on user correction. It is **not** per-gene histology
labels. It is a **navigable network**: how variants, genes and gene-networks affect
each other in MASLD, within cell types, and how Hotspot modules form at varying
degrees. Plan: `~/.claude/plans/review-this-handoff-doc-glimmering-stonebraker.md`.

### ⛔ RETRACTED: the modality split does NOT survive calibration

I built this plan around "chromatin reads fibrosis, expression reads activity",
citing its agreement with Track A's 80%/73% composition split as two independent
lines converging. **Paired contrasts on the frozen resamples do not support it.**

| contrast | H3-RNA | interval | pairs |
|---|---:|---|---|
| fibrosis_cumulative spearman | +0.110 | [-0.056, +0.293] | 4/4 positive, **0/4 exclude zero** |
| fibrosis_regression spearman | +0.105 | [-0.056, +0.283] | 4/4 positive, **0/4 exclude zero** |
| **nash_crn_component_sum spearman** | **-0.012** | [-0.107, +0.094] | **signs MIXED** |
| fibrosis_group3 macro-F1 | -0.064 | [-0.190, +0.061] | **4/4 RNA-favouring** |
| fibrosis_cumulative ordinal MAE | **+0.213** | **[+0.106, +0.319]** | **4/4 exclude zero** |

⛔ **The only clean separation is a zero-inflation output file.** Fibrosis is **71 of 99
zeros**; on ordinal MAE the H3 models predict 0.79-0.81 against RNA's 0.99-1.07 and a
true mean of 0.455 — **H3 wins by sitting nearer the zero-inflated mode**, the same
mechanism that lets a constant-0 baseline beat every fitted model. **The metric that
separates the modalities is shrinkage-driven; the rank metrics immune to shrinkage
are inside noise.** The activity half is absent on rank entirely, and one fibrosis
endpoint reverses.

**Fusion is F2, not F1.** Neither arm's interval excludes zero on the primary, and
late fusion is indistinguishable from its best parent on **all nine** metrics.
⛔ The real result is narrower and was declared in advance: **the entire penalty
lives in one prediction head** — `block_pca` vs the same parent is **-0.468
[-0.590, -0.343]** on ordinal MAE and **+0.003** on regression MAE. Carry "the
ordinal expected-value head misbehaves under the joint block", **never** "naive
fusion hurts".

Output files: `model-check-1102-21183907-paired-contrasts` (`91f15c88…`), frozen
resamples `fb7ea4b2…` verified before outcomes were read, 99/99 point estimates
reproduced, 0 refits. Counts, unordered: M 6/36, MM 1/9, F 9/54, B 58/70.

### (superseded) The one strong signal, already paid for
A **completed, audited** GSE267145 campaign (11 models × 5 endpoints × 99 paired
RNA+H3K27ac participants) shows **chromatin reads fibrosis, expression reads
activity**:

| endpoint | best H3K27ac | best RNA |
|---|---:|---:|
| fibrosis (Spearman) | **0.513** | 0.405 |
| NAS component sum (Spearman) | 0.734 | **0.776** |

⚠ **CORRECTED: fusion loses to the BEST parent, not to both lanes.** On stage3
macro-F1 both fusion arms (0.6310, 0.6328) sit **above** three single-modality
models, and `block_pca` beats its own RNA parent. The defensible claim is
"fusion does not beat the best parent" — it is worst of eleven only on
fibrosis_cumulative MAE. Track A predicts the modality split independently
(fibrosis 80% cell-intrinsic, activity 73% compositional), so each modality wins
where its measurement matches the biology.
⚠ The crossing must be tested as a **paired interaction**, not read off marginal
point estimates. ⛔ And the baseline's fibrosis Spearman is **`not_estimable`**
(constant prediction, 0 of 10,000 valid replicates) — **not** 0.0000; the 0.0047
figure is its NAS Spearman, a different endpoint.

### Three directions tested and NULL — do not re-run hoping otherwise
- **CCC**: per-donor LIANA already run correctly (117 donor parquets, pseudoreplication
  fix applied). **0 of 13 stage-progressive LR pairs pass donor FDR**; the
  donor-collapsed stage model output is a **0-byte file**.
- **SCENIC / regulons**: five cell types, all n=18, **0 significant in every one**.
- **The 43 molecular systems**: **0/43** beat a connected-subgraph null, and the
  **median system sits below its own null's median**. "33/43" did no filtering —
  plain BH passes 43/43.

### ⛔ Traps found tonight
1. **All ragged N+1 headers have one cause** — an unnamed pandas index at column 0.
   `pd.read_csv` aligns it; `awk`/`cut`/`csv.reader`/R `read.table` silently return
   the neighbouring column. Also why the lineage roster is 16, not 15.
2. **Differential correlation is confounded with differential expression.** A
   permutation null changes no gene's mean while **14,651/23,231 transcripts move
   along the axis** — the mean shift is scored as rewiring. ✅ Residualize on the axis
   first; if the effect collapses, it was expression.
3. **Never put the permissive null on the decisive test.** `system_01` clears a
   size-matched uniform null at q=0.0086 and is **rejected** by a connected-subgraph
   null. Null choice, not effect size, decided it.
4. **The 43-system graph's edges are annotation overlap** — no sample enters its
   construction, so they cannot rewire. Edge-level work is also unpowered (needs
   Δr 0.92 per-cohort).
5. **`343_harvest_documented_fstage.py` has ONE ordering bug with EIGHT victims** —
   the disease-status label map runs before the SAF regex and returns on first hit.
   7 invented stages plus **P98, where `Healthy control` overwrites a present S1A1F1**.
   P30 hits it too and agrees by luck.
6. ⚠ **The continuum-vs-stage advantage is smaller than headlined** — 2.79 compares a
   marginal against an adjusted coefficient; like-for-like 1.45–1.57;
   **reliability-corrected 0.87–1.31, which crosses 1.0**. And 2.021 is one cohort
   (n=361) against 1.15 in the other (n=109) — a ratio scaling with n is partly
   mechanical.

### ✅ Free assets found
- **Two focal Hotspot programs are PRE-REGISTERED** in the frozen prespec —
  `ductular_injury` and `stromal_ecm`. Primary at m=2, no multiplicity argument.
- **BayesPrism cell-type-specific EXPRESSION exists** (posterior Z, 13 matrices) —
  ⛔ but only **4,341 genes**, 13 non-matching types, and STAR-era stale.
- ✅ **GSE202379 steatosis and activity grades RECOVERED and frozen** —
  `executions/gse202379-saf-recovery-21183702` (`682b992e…`, 21 tests, additive:
  `343...py` untouched). n=40: S `0:1 1:17 2:18 3:4` · A `0:1 1:8 2:4 3:18 4:9` ·
  F `0:3 1:9 2:12 3:12 4:4`. ⛔ **Steatosis has FOUR levels at n=40** (an S0 appears).
  ⛔ One cohort — axes, not power; joint coverage **39** for expression, **37** for CCC.


## Showcase model — design checks

A check here decides whether a model is *definable*, not whether it scores well.
Narrative in `CHANGELOG_DETAIL.md`.

| Stage | Question | Outcome |
|---|---|---|
| 0 | Are steatosis / ballooning / inflammation separable? n=99 | **STOP** — one axis |
| 0b | How many axes at n=99? | **`ONE_AXIS_ACTIVITY`** (fibrosis underpowered) |
| 0c | Fibrosis axis at n=180? | **`TWO_AXIS_ACTIVITY_FIBROSIS`** |
| 1 | Does the assignment transfer externally? | **`TRANSFERS_WITHOUT_SPECIFICITY`** |
| 3 | Composition-mediated or within-cell? | blocked on proportion-variant provenance |

**Stage 0.** C1 ballooning\|inflammation 0.818; **C2 participation ratio 1.42 of 3.0 vs a 2.0 floor
— decisive, the number to defend**; C3 0.825/0.872/0.872. `model_fitted: false`.

**Stage 0b/0c.** Partial-association BH counts against a matched permutation count null:

| direction | n=99 (0b) | n=180 (0c) |
|---|---:|---:|
| activity \| fibrosis | 805 | **4,388** |
| fibrosis \| activity | **0** | **1,305** |

✅ **0c proves 0b's fibrosis null was thinness, not biology** — same axis, better substrate,
split-half **0.3574 → 0.4786**, marginal BH count **8 → 3,321**. The pre-declared
`not_independent_or_underpowered` guard fired correctly and was vindicated independently.

**Stage 1.** Assignment cells: activity-only **4,187**, fibrosis-only **1,104**, both **201**,
neither **56,448**. S1 met 3/3; S2 met 1/2.

| set (GSE130970) | own-axis | other-axis | reading |
|---|---:|---:|---|
| activity_only | **0.0319** | 0.0071 | **4.5:1 preference**, not specific |
| fibrosis_only | ~0.062 | below null | **axis-specific**, `tested_negative`, powered |

⛔ The binary verdict hides this: activity genes prefer their own axis 4.5:1, and the failing margin
(0.0048) sits at the MDE (0.005). **Catalog: fibrosis may claim specificity; activity claims
association plus a preference ratio and must not claim specificity.**

**Substrate.** GSE135251 n=180 trains (both axes native); GSE130970 n=78 is the graded external arm
(fibrosis 25/28/9/14/2, all five levels, join >99.9%); GSE193066 n=106 participants is
activity-only (fibrosis omits F4 → **a different instrument**, never meta-analysed).

| cohort | n | fibrosis marginal | **tie-structure ceiling** |
|---|---:|---|---:|
| GSE130970 | 78 | 25/28/9/14/2 | **0.956** |
| GSE193066 | 106 | 6/37/37/26 | **0.949** |
| *GSE267145* | 99 | **71**/15/9/4 | **0.792** |

### ⛔ Traps recorded on this lane
1. **Compute the tie-structure ceiling before modelling** — the largest rho *any* untied predictor
   could reach, fixed by the marginal. GSE267145 caps at **0.792** because 71 of 99 sit at stage 0;
   one dominant tie group costs more than being 21% smaller. Ceiling > MDE means detectable, **not**
   that an effect exists.
2. **Fibrosis at n=99 was underpowered, not empty** — max \|partial r\| **0.3767** below a
   familywise floor of **0.4502**; its best gene was undetectable. Marginal count was **8**.
3. **A scale that omits a level is a different instrument** — separate arm, never meta-analysed.
   Approved standard; GSE267145 (0-3, no F4) is not a recode hazard but a different measurement.
4. **Spearman-Brown fails on two vectors from one cohort** — returned \|rho\| **1.50-1.68**. Errors
   are shared. Use a measured shared-sample floor (mean 0.124, p95 0.312).
5. **Cross-axis vector comparison is reliability-confounded** — fibrosis split-half 0.357 vs 0.479-0.553;
   4 of 6 contrasts sit on their attenuation prediction, so separation holds **vs steatosis only**.
6. **Bootstrap df must use DISTINCT units** — resamples hold ~62.8 of 99, so `df=n-3` gave a median
   BH count of **484.5** against **65.0** correct (observed 163). Symptom: observed below its own
   2.5th percentile. Plain-correlation bootstraps are unaffected. Fix: jackknife.
7. **Count nulls are heavy-tailed** (median 0, p95 0, max ~9,000) — use the exceedance p.
8. **N+1 ragged headers are systemic in this tree — four file families so far.** Read by name,
   GSE193066's `fibrosis stage` returns **sex** and `nafld activity score` returns **fibrosis**, and
   nothing downstream looks wrong. Resolve by measured position; validate contents against the
   ordinal claimed.
9. **The eligible lineage family is SIX**, derived on canonical/kallisto in the two kallisto arms:
   Hepatocytes, Macrophages, Mono+mono derived, cDC2s, T cells, pDCs. ⛔ An earlier "five"
   (with Plasma cells and cDC1s) came from the **wrong variant**; only three overlap.
   ⛔ **Only 3 of 6 are variant-stable** (Hepatocytes, Macrophages, pDCs) — the backup-variant
   family is those three alone. Flips: Mono+mono 1.000 vs 0.528, cDC2s 1.000 vs 0.489,
   T cells 0.991 vs 0.361. ⛔ A **pooled** three-cohort rule gives just **2**, losing four to
   STAR's detection floor. Rule is per-sample detectability at ≥90%, not a median floor
   (which gives 13/7/7). ⛔ Fibroblasts (**0.227**) and endothelium (**0.361**) collapse in the
   training cohort and are fine in both external arms.
10. **Two BayesPrism variants exist per cohort** (`_star_backup`) — **not one row is identical in
    any cohort**, max differences **0.5594** and **0.5086** on a quantity confined to [0,1], on a
    byte-identical roster. Resolved: **canonical = kallisto** (77,078 genes, key `TSPAN6`),
    **`_star_backup` = STAR** (37,6xx genes, key `A1BG`); the frozen axis is kallisto, so canonical
    is the only choice consistent with the expression. ⛔ Settled by **gene universe and sample-count
    agreement, not the filename** — and the 180-row backup matching our frozen n **pointed at the
    wrong file**.
11. ⛔ **The effective-floor difference is a PIPELINE property, not a cohort one — this supersedes
    W2.** Values in (0, 1e-8): kallisto cohorts **997**, human STAR cohorts **0**, mouse 178.
    GSE240729's 7.19e-08 floor is a property of STAR, not of GSE240729. **So a detection-floor rule
    evaluated per cohort is really evaluated per quantification** — the asymmetric-denominator shape.
    ⛔ The three axis arms are **not** on one quantification: GSE135251 and GSE130970 are kallisto,
    **GSE193066 is STAR**. An independent second ground for the separate-arm rule. ✅ The fibroblast
    contrast survives — 22.7% vs 100% is kallisto-vs-kallisto.
12. ⛔ **A manifest must record what it EXCLUDED, not only what it included.** A
    `*/*_bayesprism_proportions.tsv` glob silently drops `_star_backup.tsv` — a filter that reads
    like a listing, accurate about what it pins and silent about what it omitted.
13. **Proportion row sums are not bitwise 1.0** — only 44/216, 7/78, 19/164 exact. Use a tolerance.
14. **The GPU dispatcher enforces an exact queue-item key set** — no `disabled_reason` field.

**Evidence states.** NAS-component separability in GSE267145 = **`tested_negative`**; aspect
**attribution** = **`indeterminate`**; **SAF is settled by neither**. `validate_resource_scope.py` PASS.

⛔ **In-cohort at training. The fibrosis axis has exactly one external arm, n=78, no assertable
donor key.** External support is **set-level, never per-gene**.

**Two Cobolt queue items disabled** (052 v8, 057 v9-generic6; 50 tasks) — separation-doomed, and v9 is
on record as never evaluable. Not re-frozen; 67 items remain enabled.

## Endpoint results

### Fibrosis / NASH / NAS — the goal's core lane

| Model / lane | Cohort, endpoint | Result | Read |
|---|---|---:|---|
| **gene_median_ridge** | GSE135251→GSE267145, NAS sum, n=99 | **Spearman 0.756** | ✅ verified working transfer; null p95 0.165 |
| gene_median_pca_ridge | same | 0.735 | all 4 learned models clear their nulls |
| per_array_rank_(pca_)ridge | same | 0.662 / 0.638 | |
| training_mean | same | 0.000 | constant null |
| gene_median_elastic_net | GSE135251→GSE83452, NASH, n=148 | macro-F1 **0.6224** | beats prevalence 0.4127 decisively |
| gene_median_pca_elastic_net | same | 0.6017 | strongest comparator (post-hoc) |
| per_array_rank_elastic_net | GSE267145→GSE49541, fibrosis, n=72 | AUPRC 0.567, p=0.072 | candidate; **below its own null mean** |
| gene_median_elastic_net | same | AUPRC 0.692, **p≈7.5e-4** | real signal but **post-hoc**; predicts 0 positives, macro-F1 0.357 |
| H3-only linear SVM | GSE267145 in-cohort, NOR/NAFL/NASH | macro-F1 0.691 | best 3-state source baseline |
| RNA nearest centroid | GSE267145 in-cohort, CRN sum | Spearman 0.776 | within-cohort; **different estimand** from transfer |
| Prespecified RNA SVM | → GSE260666, n=16 | macro-F1 0.419 | weak external transfer |
| **Lineage composition → fibrosis** | 3 bulk cohorts, n=350 samples | Hepatocytes **rho −0.439** (q 1e-13) | ✅ operator **validated**; a positive control, not a discovery |
| ↳ Macrophages | same | +0.225 → **−0.165** under closure | ❌ closure output file; sign **reverses** |
| ↳ pDCs | same | +0.133 | eligible, below effect floor, null |
| Donor × lineage pseudobulk | GSE296875, n=37–39 | **0/12 survive BH** | power ≈ **0** at realistic effect sizes |

**Power, frozen.** 80% power under BH needs **~150 donors** for a moderate effect (AUROC 0.70 /
rho 0.30). GSE296875 has 37–38; pooled with GSE202379, 77. The world ceiling for public liver
single-cell with donor histology is 86 donors, embargoed. This is why the lineage lane moved to
a bulk deconvolution substrate.

### Variant → regulatory activity (GSE281364 MPRA, 1,033 elements, 5 seeds)

| Model | Score | Gain vs allele ridge | Read |
|---|---:|---|---|
| **HyenaDNA delta** | **0.2162** | **+0.0820 [+0.0274, +0.1389]**, 5/5 seeds | strongest completed regulatory result |
| Sei native 40-class | 0.1939 | +0.0596, CI crosses zero | promising, not separated |
| Nucleotide Transformer delta | 0.1879 | +0.0537 [+0.0011, +0.1108] | positive; restricted comparator |
| Caduceus delta | 0.1449 | +0.0107, crosses zero | |
| *allele-identity ridge* | *0.1342* | *reference* | strong simple baseline |
| gkm-SVM / deltaSVM | 0.126 / 0.127 | negative | did not adopt |
| DNABERT-2 delta | 0.0964 | negative | |
| Enformer SAD / SAR | 0.083 / 0.085 | **negative, 0/5 seeds** | fails on MPRA transfer |
| MPRALegNet | — | — | **terminally blocked**; scoring forbidden |

### Chromatin / RNA-ATAC (GSE296875)

| Model | Result | Read |
|---|---:|---|
| **Observed-ATAC GLM** | deviance skill **0.1627 [0.1520, 0.1733]** | strongest control; all 5 lineages positive |
| Shuffled-modality control | 0.0561 | retained to expose shortcut structure |
| Observed-ATAC-only | 0.0313 | much weaker |
| RNA-only / masked | 0.0176, CI crosses zero | **no established RNA-conditioned skill** |
| PeakVI transfer smoke | **−0.0559**, 0/5 lineages | negative |
| MultiVI smoke | 0.79% vs training mean | below the 5% check |
| scGLUE paired / StabMap | MRR 0.0444 / 0.0281 | both below linear CCA 0.0503 |
| scBasset preliminary | −0.4876 | negative screen |
| Corgi frozen head | 2.5e-7 | zero incremental skill |
| BPNet / ChromBPNet / CNN / transformer | rectangle in progress | **no partial ranking permitted** |
| Cobolt v12 | 25/25 GPU + 150/150 CPU controls verified | **unscored**; 1,000-nucleus smoke fixture |

### Cell-state mapping (broad identity, not disease state)

| Model | 50k donor-balanced | Study-balanced | Eligible to be the best model? |
|---|---:|---:|---|
| TF-Sapiens MLP | **0.9389 [0.9275, 0.9486]** | 0.9260 | ❌ exposure `unknown` |
| hvg_pca_knn | 0.9075 [0.8934, 0.9190] | 0.8178 | ✅ best eligible |
| linear SVM / logistic | 0.9054 / 0.9022 | — | ✅ |
| elastic_net | 0.8898 | **0.8835** | ✅ best on transfer |
| nearest centroid | 0.7016 | — | ✅ |

**⚠ elastic_net's study-balanced edge is carried by one 149-cell study.** Post-hoc concentration
diagnostic over the 7 studies: elastic_net's gain over linear_svm (**+0.0215**) and over logistic
(**+0.0238**) both **reverse sign** when GSE189600 — 149 cells, 2 donors, ~0.3% of the set — is
dropped (−0.0251 and −0.0227). It is the most influential single study in **10 of 20** pairwise
comparisons. An unweighted macro-average over studies is the right choice for a transfer question;
the margin must simply never be quoted without saying what carries it.

**The ranking flips between metrics** — kNN leads within-study and is worst of the four on
cross-study transfer. Paired on the same donors: kNN−elastic_net **+0.0169 [+0.0070, +0.0266]**;
TF-Sapiens−elastic_net **+0.0416 [+0.0323, +0.0513]**. At 1,000 cells Geneformer V2 316M reaches
0.9522 (+0.0181 [+0.0044, +0.0320]). These models predict **cell identity**, not disease remodeling.

## Cross-cutting methods findings

1. **The source-only selection rule is unreliable inside a seed-width.** GSE83452's top two source
   models differ by 0.00037 against a seed sd of 0.00074; GSE49541's best exceeded no null at all.
   The NAS lane — gap 0.0347 against a spread of 0.00077 — selected correctly. The rule works when
   models are genuinely separated and is arbitrary when they are not.
2. **Prevalence is never the AUPRC null.** Measured random-scorer nulls sit *above* prevalence in
   all three cohorts (0.474 vs 0.444; 0.712 vs 0.703). Scoring against prevalence previously
   produced six false positives. Use each model's own realised score-vector tie structure.
3. **Bootstrap intervals must be paired.** Overlapping marginal intervals are not evidence of
   absence; `bootstrap_distributions.npz` draws are index-aligned across scoring runs.
4. **Vacuous check conditions inflate headlines.** A deterministic fit satisfies a seed-direction
   condition for free. Checks now record `not_applicable` and report *met / applicable*.
5. **Report p-values at their resolution floor.** 3/10,000 is not 3.0e-4; state p ≤ 1e-4.
7. **Ask what would satisfy a check condition vacuously, before freezing it.** Four instances of
   one shape: a seed condition on a deterministic fit; an effect floor scored while also filtering
   upstream; a dominance guard on a closed composition (the one that got through); and a condition
   whose subject the evaluator read from a hardcoded constant rather than the check. A condition the
   analysis never computes is the same defect in a different costume — it reports absent rather than
   failing — so every evaluator must assert its computed condition set equals the check's registered set.
9. **A global permutation null is too easy when the evaluation set has structure.** Permuting
   *within* declared strata preserves the nuisance and destroys only the signal. On a fixture where
   the score carries no within-stratum information, the global null gives **p ≤ 1e-4** and the
   stratified null **p = 0.976** — same data, same observed value. Now a procedural requirement:
   where strata exist, report both. `evaluators.auprc_reference.stratified_permutation_reference`.
10. **A pooled gain must be reported with where it comes from.** Leave-one-stratum-out concentration
   is mandatory alongside any gain claim, and is never pass/fail — there is no defensible threshold
   for "too concentrated". `evaluators.stats.gain_concentration`, `not_applicable` below 4 strata.
11. **`all([])` is `True` — a check with zero applicable conditions can PASS having tested nothing.**
   Applying `not_applicable` honestly can empty the applicable set, and a boolean aggregation over
   an empty list then reads as success. Note the earlier fix *created* this trap: before honest
   `not_applicable`, the same check would have reported a misleading `met/total`. Every evaluator must
   handle 0/0 explicitly. Third verdict state: **`NO_APPLICABLE_CONDITIONS`** — not PASS (nothing was
   demonstrated), not FAIL (that would assert the outputs were tested and found wanting).
12. **A row count is not a unit count.** The Zhu substrate is 282 variants x 2 cell models x 2 stimuli;
   both 1,128 and 564 were row counts of a 4x expansion, and 17 "DAV rows" were 5 distinct variants.
   Every power number quoted before the distinct-unit check was wrong. Same shape as SRR-runs-as-donors
   and the 273-to-125 donor collapse. Report distinct units beside row counts and say which the
   arithmetic uses.
13. **Unit alignment with external evidence is a registration-time constraint.** It is invisible until
   an external corroboration lane is attempted, by which point outputs are frozen at the wrong unit.
8. **On a closed composition, "different component" is not "independent component."** A dominance
   guard must require the qualifying non-dominant component to survive subcomposition. Testing that
   two components are distinct objects does not test that they vary independently.
6. **Measure a column's contents, never trust its name.** Bulk metadata and proportions files carry
   an R `write.table(row.names=TRUE)` offset (header N, rows N+1), and a column labelled `SYMBOL`
   is 46% ENSG. Assert the N+1 invariant and fail closed.

## Live compute (volatile — refresh from SLURM)

| Lane | State |
|---|---|
| GPU rectangle | 4 running / 21 dependency-held; **queue backlog is 0** (all enabled items claimed) |
| Remaining GPU work | 25 jobs / 345 logical tasks |
| Dead dependencies | none |
| W2 composition lane | **complete.** Verdict + interpretive record + forward-only check v2 frozen; v1 byte-unchanged |
| W3 prioritization lane | **closed as a recorded negative.** Substrate cannot support it; no check registered by ruling |

## Metric authorities

- **NAS transfer** `executions/model-eval-934-21134193/evaluation` — verified by `verify-nas-independent-audit-20260826`
- **GSE83452** `executions/model-eval-890-21130318-nash_vs_nafl-eval` — verified by `verify-gse83452-independent-audit-20260826`
- **GSE49541** `executions/model-eval-866-21130184/evaluation` — verified by `verify-gse49541-independent-audit-20260826`
- **Power statement** `executions/gse296875-power-statement-21129430`
- Cell state `RESULTS_2026-08-23_cell_state_mapping.md`; TF-Sapiens `executions/model-cpu-score-504-21083585/`
- Classical 50k `executions/model-cpu-train-503-21097161`, scored by `model-cpu-score-502-21097344`
- DNA-language `executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv`
- Enformer/Sei `executions/model-cpu-train-604-21097189/evaluation/family_native_summary.tsv`
- Observed multiome `executions/model-cpu-score-328-21100332/receipt.json`

## Change log index

One line per entry; full text in `CHANGELOG_DETAIL.md`.

| Date | Entry |
|---|---|
| 08-29 | **GSE240729 resolved and the signature tests reconciled.** Five of six anomalies are one operand: rung a is +7.45 SD there while rung d is −0.22, so **composition over-performs**; carrier is the fibroblast share (0.711, pool max) and every deconvolution-free probe is also at pool max — it has the **cleanest fibrosis gradient in the pool**. Real damaged-RNA assay (mapped/input 0.923, deletion 0.13%, AUC 0.000/1.000) but composition retains 80.2% after conditioning on nine technical metrics. FFPE vs cohort **cannot be separated**. Signature tests: the A/B flip is entirely the **null geometry** — a de novo set selected on our own data also fails A (p 0.101). Test B is not a coherence test (within-arm EVR correlates −0.153 with the statistic). q10 prespec held-back with the performance-selection included; **no label-free case-mix guard exists**. **UNVERIFIED and withdrawn: the GSE240729 continuum figures (0.7100, PC1 0.715) have no locatable producing output file.** |
| 08-29 | **NARROWING: "no signature beats matched chance" is refuted; the diffuseness claim is rewritten.** Pooled five-cohort unsupervised membership test with a permutation-calibrated null: **14 of 19 clean published sets beat a size- and expression-matched null**, Kamzolas at 0.5535 vs null median 0.1588. The supported claim is that feature identity buys little **once a model is fitted**, plus the sharper fact that a random 139-gene set PC1 is the global liver axis (|r| 0.981) worth macro |rho| 0.1587 on its own. **Passing indexes fibrogenesis, not MASLD** — Hallmark EMT 0.438 and NABA matrisome 0.401 beat the WikiPathways NAFLD set (0.1089, p 0.443, the worst arm). Two tests (per-cohort held-back-projection vs pooled membership) were briefly conflated and must be reconciled before either is figured. |
| 08-29 | **Third distinct null failure mode recorded: degenerate in DIRECTION.** Under permuted grades the specified matched null fires at 24-38% for real signatures vs 4.3% for random controls, because a random set PC1 has \|r\| 0.981 with the global PC1 — the null occupies one direction with sd ~0.015 and any coherent module aimed elsewhere clears it. Repaired by permutation calibration with the preference order recorded before the run. Joins the vacuous matched null and the point-mass permutation null in the ledger. |
| 08-29 | **Graded-fibrosis instrument consolidated.** C1 +0.1733 → +0.2333 after deconvolution harmonization; C2 falls below floor (V4→V3); C3 dies. Severity is **diffuse** — curated 139 genes fail a matched null, 20k random H3K27ac regions match a selected set, 139→39,798 genes buys +0.020. Chromatin recovers the RNA ordering at 0.795 (0.487 after conditioning on true NAS). `axis_raw` offset **1.39–1.87 stages**, label-free q10 anchor fixes it but n=1 RMSE is 2.85 stages. Conformal 90% interval spans 69% of the scale. Continuum coefficient exactly zero inside the model; cold-start crossover n=160. |
| 08-29 | **RETRACTION: "the continuum does not detect disease" is reversed.** It rested on PRJNA512027, retired 2026-05-15 for an L0/S0 library-prep × disease confound. The axis detects at AUC 0.81–0.89 in four cohorts with already-frozen scores; the composition comparator there was reading batch. |
| 08-26 | **W2 scored: check PASS 4/4 applicable but NOT a composition finding.** Macrophages reverses sign under closure (+0.225 → −0.165); dominance guard satisfied by the output file it excluded |
| 08-26 | W2 collapse filter cut the testable family 16→3; Cholangiocytes (rho 0.426, below floor in 96.5% of GSE135251) would otherwise have been the headline |
| 08-26 | NAS transfer **verified** — 0.756 survives leak hunt: not library size, not a global axis, not the gene set alone, not the 24 tied zeros |
| 08-26 | Determinism threshold weakness recorded — 0.999 correlation is blind to spread; deliberately not retrofitted |
| 08-26 | W2 lineage lane scoped; CLR rejected (819 values <1e-8), rank primary adopted; shared-reference dependency named |
| 08-26 | All 3 unadjudicated bulk-pool cohorts fail the **NASH** requirements — scope caveat: decides nothing about fibrosis endpoints |
| 08-26 | NAS lane scored: transfer works, check FAIL 1/3 applicable; the seed fix fired on its first real use |
| 08-26 | GSE135251 included — two narrow admissions; `positive_required` met only via the documented-equivalent clause |
| 08-26 | GSE83452 scored + verified: FAIL 3/4; one condition turns on a single participant |
| 08-26 | GSE49541 scored + verified: FAIL 0/4; source had 4 advanced cases |
| 08-26 | Power statement frozen: ~150 donors needed for a moderate effect |
| 08-26 | GPU backlog is 0; the three "stale" 1,033-task bundles were already complete |
| 08-25 | Independent audit: every headline number reproduces; three *descriptions* did not |
| 08-25 | Model census closed 136/136; pin cascade 16→0 |
| 08-25 | Cobolt v12 25/25 on B6K — first bundle to survive the resource separation |
| 08-25 | GSE83452 + GSE49541 summarised through their inductive separations |
| 08-25 | First phenotype result: a comprehensive, well-characterised **negative** |
| 08-25 | 50k classical baseline complete; ranking flips between metrics |
| ← | earlier entries in `CHANGELOG_DETAIL.md` |

## Update requirements

- Update only the affected row when a fit, evaluator, transfer, missing piece or held-back status changes.
- Keep distinct endpoints in distinct tables. There is no universal leaderboard.
- Never report partial-rectangle metrics. Preserve negative results and terminal failures.
- No number enters this file without independent re-derivation by a different agent.
- Add a one-line change-log entry here; put the detail in `CHANGELOG_DETAIL.md`.
| 2026-08-27 | **Two evaluation conventions adopted from the oFM paper (arXiv:2608.24688), and the retrospective already found something.** Frozen at `executions/model-check-904-21178510-evaluation-conventions` (`3068cef1333a2ed012241cb729bfde3f71dbb9b0cbcc95f544c08e5cbdd77de0`), additive and **forward-only** for enforcement; every read input verified byte-identical before and after. **(1) Stratified permutation nulls, a PROCEDURAL requirement.** Their null reshuffles within propensity-score bins rather than globally; every null this campaign had computed was global. Demonstration regenerated inside the freeze job: on a fixture where the score carries **no** within-stratum information about the label, the global null gives **p <= 1e-4** (resolution floor) and the stratified null **p = 0.976** — same data, same observed AP 0.784. The control matters as much: with strata assigned at random the two nulls converge (0.5162 vs 0.5165), proving the gap comes from structure and not the new code path. Guards fail closed — non-contributing strata are counted not dropped, and all-degenerate raises rather than returning a meaningless p-value. **(2) Gain-concentration, a MANDATORY REPORTED diagnostic, never pass/fail** — no defensible threshold exists for "too concentrated". Leave-one-stratum-out rather than Cochran's Q, because it is what actually caught GSE189600; `not_applicable` with **no numbers** below 4 strata, the k=3 lesson applied at design time. **Retrospective on the 50k cell-state lane (post hoc, revises nothing):** 7 studies, 20 comparisons, all applicable. **elastic_net's gain over linear_svm (+0.0215) and over logistic (+0.0238) both REVERSE SIGN when GSE189600 is dropped** (-0.0251, -0.0227). That study is 149 cells from 2 donors, ~0.3% of the evaluation set, and is the most influential single stratum in **10 of 20** comparisons. This extends the 2026-08-25 audit rather than repeating it: there the margin collapsed 0.066 -> 0.0004; here it inverts. **Not adopted, recorded so nobody re-adopts them:** the oFM architecture (scale-dependent at 1.67M patients), its longitudinal framing (needs repeated observations we do not have), and **its use of prevalence as the average-precision chance line** — this project established prevalence is not the AP null, and their strata reach n=20 where that bias is material. **A methodology error of this session is also recorded:** a first suite run compared a whole-tree `unittest discover` against `PREEXISTING_FAILURES`, an empty frozenset calibrated for three modules, and reported 35 unrelated pre-existing failures as regressions. Re-run at the correct scope: **127 tests, `new_failures: []`**. The wider tree does carry pre-existing failures outside the baseline's scope, so "the debt register is empty" describes three modules, not the tree. |
| 2026-08-27 | **Stratified nulls applied to the completed microarray lanes: both hold. No conclusion changes.** `executions/model-check-905-21178664-stratified-nulls` (`e6b60c7e72ff6ab72869aa050474313be74740d7ffbd3e85a43e022ff9e6bffb`), post hoc, revises nothing, 20,000 permutations, every read input byte-identical before and after. **The stratum is SCAN DATE**, a genuine technical batch for a microarray assay - shared reagent lot, instrument state, operator - chosen because it is real, not because it would produce an effect. GSE49541 has **3 scan days** (29/25/18 arrays); GSE83452 has **13** at baseline, clustering into two campaigns 14 months apart (Oct-Nov 2013, Jan 2015). **A prediction was made before running and it held.** Batch is not confounded with the outcome in either lane - GSE49541 is 44%/40%/48% advanced across its three days, and GSE83452 is **70% NASH in both the 2013 and 2015 campaigns** - so the global null had nothing to wrongly credit and the stratified null should barely move. It barely moved: `models_whose_conclusion_changes_at_0.05` is **empty**, and every shift is <= 0.011 with most under 0.002. GSE49541 `gene_median_elastic_net` 0.00110 -> 0.00045; GSE83452 `gene_median_elastic_net` 0.00030 -> 0.00245, `gene_median_linear_svm` 0.00015 -> 0.00135, `per_array_rank_elastic_net` 0.03355 -> 0.02320. Observed AP values reproduce the frozen ones independently (0.6923, 0.8269, 0.8394). **So the p-values reported for these two lanes stand as published.** That is the useful outcome: the convention was applied where a real stratum exists and returned reassurance rather than a correction. **Scope limit, stated rather than glossed:** the stratified null as implemented covers average precision and Spearman. The 50k cell-state lane's endpoint is donor-class-balanced macro-F1, which needs its own permutation machinery, and its labels are evaluator-only by construction (`prediction_tables_contain_observed_labels: false`) - so that lane received the **concentration** diagnostic instead, which is where its real problem was found. The NAS lane's target is a single cohort whose evaluator-only labels carry no batch variable, so no stratum exists to apply. |

---

## RNA-conditioned H3K27ac: the first clean cross-modal positive — 2026-08-30

**GSE267145, n=99 participants, 5-fold participant-held-out, reduced-rank regression at rank 4.**
Job 21248730, `executions/gse267145-crossmodal-ceiling-20260830/results_rna2h3_v2_20260830T170000`.

Predict a held participant's H3K27ac profile from their bulk RNA. This is the direction
`paired_bulk_rna_h3k27ac_task.toml` permits; held H3K27ac is evaluator-only.

| quantity | value |
|---|---|
| participant-macro multinomial profile deviance skill | **+0.13345** |
| permutation p vs destroyed-pairing null (100 perms) | **0.000** |
| destroyed-pairing null | mean −0.04968, sd 0.01805, **max −0.00376** |
| folds positive | **5 of 5**, range +0.11742 to +0.14886 |
| observed-pair retrieval (separate lane) | MRR 0.3455, top-1 0.1616 vs chance 0.0505 |

Rank 4 was taken from the cross-covariance ceiling, not chosen by search: components 5 and 6
fail the pairing null and a top-k method cannot step over them. So there is no selection
inflation to discount.

The fold range is tight and every fold is positive, so the **concentration diagnostic is clean** —
this is not one cohort-fold carrying a pooled number.

⚠ **What this does and does not say.** A deviance skill of 0.133 means 13.3% of the excess
deviance over a training-mean profile is explained. Modest, real, and highly consistent. It is
**one cohort**; by this project's own standard that is not a verdict. The second substrate is
GSE296875 (39 donors, same-nucleus multiome), and the question there is whether the relationship
*recurs*, not what a pooled estimate would be. The zero-inflated-fibrosis caveat that retracted
the earlier "chromatin reads fibrosis" claim does **not** apply here: this lane predicts H3K27ac
profiles from RNA and never touches a fibrosis label.

⛔ **The retrieval lane is reported separately and cannot support an RNA-conditioned claim.**
Top-1 at 3.2x chance is a statement about matching observed pairs, not about prediction.

### The gate had to be rewritten first, and the rewrite is the story

The first run (21247860) stopped at its own gate and printed:

```
null skill: mean -0.04968  sd 0.01805  2.5% -0.07641  97.5% -0.01290  max -0.00376
GATE FAILED - the pipeline finds skill in unpaired data
```

**All 100 permutations scored below zero.** The message asserts the opposite of the line directly
above it. The condition was `abs(null_mean) < 0.02` — a two-sided bound on a quantity with a
known one-sided expectation, since a reduced-rank fit on scrambled pairings *must* land below a
training-mean baseline. It failed because the null was too clean.

The fix was not to loosen it. It is now one-sided on positive skill, **and** carries the separate
criterion the two-sided version was clumsily reaching for:

> When the null does not centre on zero, **"beats the null" and "is useful" are different claims.**
> With the null at −0.05, a result of −0.02 beats destroyed pairing while still losing to the mean
> profile. Both are reported, and the reading of all four combinations (A/B/C/D) is printed by the
> script *before* the primary is computed.

This run came out **A**: p < 0.05 and skill > 0, the only combination that licenses the sentence
"RNA predicts held-participant H3K27ac".

`scripts/91_assert_the_gate.py` tests the rewritten gate by constructing its failure on purpose —
the identity permutation, where pairing is preserved. Result: **fires on preserved pairing
(+0.13345), passes on destroyed (max −0.01470). A demonstrably live guard, not a check that
happened to return green.**

⚠ One honest wrinkle: the assertion script needs the paired skill to distinguish "guard is dead"
from "nothing was there to leak", so it prints the primary before step 2's gate formally opens it.
The A/B/C/D reading was committed in code before either step ran and the gate keys only on the
null, so no branch could be rationalised — but the ordering is not the clean one.

---

## Two dormant cohorts activated — 2026-08-30

`executions/cohort-activation-20260830T195335Z/`. Both ran in under 3 s; no SLURM job was warranted.

### GSE193066 — 58 within-participant pairs recovered

**106 participants, 164 sample rows, 58 pairs.** The only bulk cohort in the project with both a
real donor key and within-participant repeats.

The title convention: paired participants are `HUnafld035_1` / `HUnafld035_2`; **unpaired titles
carry no suffix at all.** The correct rule is `re.sub(r'_\d+$', '', Sample_title)` applied to
*every* title. The naive rule — strip `_2`, look up the bare id — matches **0 of 58**, because the
bare id is not itself a title in the file.

Five independent confirmations that the 58 are not spurious: exactly one 1st and one 2nd biopsy per
pair; sex constant across all 106; the published `HUn106.gct` header equals exactly the 106
first-biopsy titles and `HUn164.gct` equals all 164; the deposit's own prose says "the 58 high-risk
patients with paired biopsy tissues"; and a pre-existing trajectory file names the same 58 with zero
disagreement either way.

⛔ **A selection limit that travels with every paired contrast: all 58 paired participants are
high-risk at first biopsy**, while the 48 unpaired split 35 low / 13 high. The pairs are not a
random sample of the 106.

⛔ **No interval field exists.** The deposit carries seven keys and none is a date. Per-sample age
gives median 2 years, range 1-7 over 57 pairs. HUnafld080 records age 46 → 41, impossible; sex and
biopsy labels agree, so it is a deposit inconsistency, pinned as an enumerated exception so a *new*
negative delta still raises. Fibrosis delta across the 58 pairs is **−0.0345** (28 unchanged, 15
down, 15 up).

### GSE202379 — all three SAF aspects recovered, and steatosis really does have four levels

At n=40 graded (47 donors total):

| aspect | levels | counts |
|---|---|---|
| Steatosis S | **4** | S0=1, S1=17, S2=18, S3=4 |
| Activity A | 5 | A0=1, A1=8, A2=4, A3=18, A4=9 |
| Fibrosis F | 5 | F0=3, F1=9, F2=12, F3=12, F4=4 |

⛔ **The single S0 donor is P30.** A 1-3 steatosis scale silently drops exactly that donor — the
failure is invisible because it looks like a missing value, not an error.

n by use: **40 graded → 39 expression** (P70 absent from the substrate) **→ 37 cell-cell** (P62 and
P67 fail LIANA's two-lineages-at-30-cells gate, a library-depth property rather than missing
phenotype).

The source fix in `Analysis/SingleCell/scripts/343_harvest_documented_fstage.py` widens the regex
from `S\d+A\d+F(\d)` to `S(\d)A(\d)F(\d)` and emits `saf_S`/`saf_A`/`saf_F`/`f_stage_source`.
`F_stage_documented` is **value-identical pre- and post-fix on all 59 samples**, verified by running
both module versions side by side — so recovering two aspects changes no existing number.

⚠ A separate label-map precedence defect was deliberately **not** fixed: P98's real F1 is
overwritten to F0, P30 hits the same short-circuit and agrees by luck, and 7 donors carry an F with
no SAF string. Fixing it would silently move a frozen value; it is now *visible* through
`f_stage_source` instead. The frozen `donor_fstage_documented.tsv` was not rewritten.

⚠ These are three axes on one cohort at n=39/37 — **axes, not power.**

### Guards

GSE193066: 10 demonstrably live, 1 verified positively, 1 live only under a changed constant
(`unified_metadata.csv` measured at header NF=9 / row NF=9, so it is **not** ragged and the ragged
branch cannot be exercised on it). GSE202379: 5 demonstrably live, 2 verified positively, 1 live
only under a changed constant. 31/31 unit tests pass.

Violations constructed on purpose that fired: the naive `_2` rule (0 pairs → raised); the
participant key permuted under seed 20260830; GEO fibrosis shifted +1; the 106-column GCT header
compared against 164 titles; a lowercased SRR key set; the unstripped `GSE202379_<id>` join; reading
`group(1)` of the three-group pattern as F; and a 1-3 steatosis scale. **Every zero join raises and
names the rule; none degrades to `applicable: false`.**

---

## Liver zonation from expression: positive — 2026-08-30

**GSE105127, 19 participants x 3 zones, 57 libraries, participant-held-out folds.**
Job 21249083 (4 min), `executions/gse105127-zonation-lane-20260830T195218Z/`.

| quantity | value | pre-registered bar |
|---|---|---|
| accuracy | **0.7895** (45/57), CI95 [0.702, 0.877] | chance 0.3333, 0.05 bar 0.4737, MDE 0.546 |
| mean absolute zone error | **0.2281**, CI95 [0.123, 0.351] | chance 0.8889 |
| refit-permutation p (200 perms) | **0.00498** — the floor at 200 | null mean 0.329, q95 0.422 |
| exact within-participant p | 1.008e-11 | |

Clears every pre-registered bar. The bars were fixed, and the `indeterminate` label for a negative
was written down, **before the result was seen** — this is the project's only zonation cohort, so a
negative here could never have been separated from insufficient power.

⭐ **The ordinal structure is respected: 11 one-step errors and only 1 two-step.** That is why mean
absolute zone error is the right readout and accuracy alone would have thrown it away. Confusion:
CV 18/21, IZ 11/14, PP 16/22.

### The row-order trap was live, and it would have produced a perfect fake

⛔ **All 19 participants have their rows ordered CV/IZ/PP, and row position ALONE scores accuracy
1.0.** The trap named in advance was not hypothetical in this dataset — it was sitting there. Rows
were shuffled before fitting and predictions were shown invariant to shuffling and to
within-participant reversal.

⛔ **A second, independent vacuous null was caught and the claim withdrawn.** The *constrained*
assignment arm (one library per participant-zone cell, solved as an assignment problem) reaches
accuracy 1.0 — but a **uniform-probability model also reaches 1.0 under it, at p = 1.6e-15.** That
null cannot reject a known-empty model, so it proves nothing and **no p-value is quoted for the
constrained arm.** The free-argmax null passes the same control correctly: a uniform model scores
0.333 at p ≈ 1.0.

This is the seventh unusable null in this campaign, and the first one caught by running the
known-empty-model control *as a precondition* rather than discovering it afterwards.

### Guards

**15 demonstrably live, 2 verified positively, 0 dead.** Six violations were constructed on purpose
and all six fired, including the decisive one: a one-hot row-position model — literally the
tie-break-by-row-identity solver the trap describes — was fed to the invariance guard and it fired,
reporting accuracy 1.0.

⚠ **What this does not license.** One cohort, n=19 participants. The result says liver zone is
recoverable from bulk expression in these participants; it is not a validated zone classifier, and
there is no second cohort in this project to test transfer against.

---

## `sequence_native_regulatory` registration: closed-negative — 2026-08-30

The belief that Enformer, Sei, Borzoi and the DNA language models read as "blocked" *because* the
task definition is unregistered is **wrong**, and registering it would surface nothing.

1. **The existing scores belong to a different, already-registered task.** Enformer, Sei and four DNA
   LMs are scored under `task_id = variant_to_regulation`, which is active. They are already in this
   file's MPRA table. Nothing was hidden.
2. **Borzoi has no scores at all** — nine execution dirs, all admission/fixture/preflight. Its status
   is `blocked_terms`: "Official weight terms are undeclared." **That is a licence block; no registry
   edit can lift it.** Same for `scooby_onek1k`. `evo` has no score either.
3. **Zero executions carry `task_id = sequence_native_regulatory`** — 0 of 1,684 task_spec files.
   Registering it would surface nothing that exists.
4. `config/evaluation/sequence_native_regulatory_task.toml` already exists, written today and
   deliberately parked **outside** `config/tasks/` with a `registration_note` giving the reason.

⛔ **The cost of registering anyway, traced in a scratch copy:** `Registry.load` rejects the parked
file on an unknown `evidence` field — and `TaskSpec` has no `notes` field, so the usual
"unknown tables go under `notes`" workaround does not apply. After schema repair the loader fails
*deeper*: all 10 baseline models fail task-compatibility because none declares the task in
`supported_tasks`, and fixing that means mutating `config/models/*.toml`, which **22 overlay records
pin by SHA-256 as `base_registry_unchanged: true`.**

The other four unregistered tasks (`histone_regulatory`, `methylation_regulatory`,
`observed_multiome`, `protein_transport`) hit the same two-layer blocker. All four left alone.

**A genuine `sequence_native_regulatory` result needs a scoring run, not a registry edit.**

⚠ Do not mix blockings: the DNA-LM numbers in `config/dna_lm_family_reconciliation.json`
(ridge 0.1771) use 1,033 outer-locus blocks; the MPRA table above uses 239 long-range blocks.

---

## GSE268273: already quantified, and it CANNOT be pooled with the five-cohort graded set — 2026-08-30

`executions/gse268273-evaluator-phenotype-20260830T205701Z/`. Acquisition and quantification
completed **2026-08-25** (plan `model-data-075-21082574`, campaign hash `fb284966…f8d3a`); the whole
DAG ran through admit, independent audit and readiness freeze. ⛔ `AGENT_HANDOFF_2026-08-25.md` is
**stale** — it lists the quantify array as RUNNING. Nothing was re-downloaded.

On disk: 968 FASTQs, **512,881,099,727 bytes (477.66 GiB, byte-exact to the plan)**, all 109
`rsem.genes.results` at 60,719 genes, and a 109 x 14,078 counts matrix.

### Counts: 824 runs / 109 samples / 109 participants

Confirmed four independent ways — 109 unique GSM, 109 unique patient codes, 109 unique BioSamples,
109 unique SRX. Runs are technical partitions: **680 single-end from 73 participants, 144 paired-end
from 36**, with `--paired-end` genuinely used for those 36 (tallied from the per-participant
receipts, not from a flag). ✅ All four metadata tables measured for raggedness first — **none are
ragged**; the endemic pandas-index defect does not apply here.

### ⭐ A well-powered graded cohort: fibrosis ceiling 0.9523 at n=109

F0 36, F1 41, F2 13, F3 14, F4 5 — largest tie group 41 (37.6%). **Not** capped the way GSE267145 is
at 0.792. Per stratum: IMID n=69 ceiling **0.9401**, classic n=40 ceiling **0.9664**. This is the
best-powered graded fibrosis substrate in the project after the 521-sample pool.

### ⛔ The reference does not match the six locally reprocessed bulk cohorts

| | GSE268273 (and GSE105127) | the six bulk cohorts |
|---|---|---|
| genome | `refdata-gex-GRCh38-2024-A` **p14 + patches** | `fasta/genome.fa`, 10x **primary** assembly |
| quantifier | **RSEM 1.3.1** + STAR 2.7.10b | **featureCounts 2.1.1** + STAR 2.7.11b |
| strandedness | **unstranded** | **reverse-stranded `-s 2`** |

GENCODE v49 is shared; nothing else is. Three consequences before any pooling: RSEM
`expected_count` is EM-distributed and fragment-scaled while featureCounts emits integer read/pair
counts (this project's own read-vs-fragment note records a **~1.98x** gap); unstranded vs
reverse-stranded diverges on antisense-overlapping genes; and p14-with-patches vs primary assembly
changes the target set. **Treat GSE268273 as a separate arm, not a sixth graded cohort.**

⛔ **The delivered matrix is 14,078 genes, not the transcriptome** — that axis is inherited from the
source's own filtered supplementary table crosswalked to v49, a *source-chosen* gene set. The full
60,719-gene RSEM output is retained per participant, so widening it costs a re-consolidation, not
~90 CPU-hours of realignment.

### Two phenotype defects found and repaired

⛔ **Case-inconsistent yes/no.** `hypertension` and `t2d/ir` both ship as `No`/`no`/`yes`, so
grouping on the raw string yields **three strata instead of two**. After casefold: hypertension 58/51
(**42 rows changed**), t2d/ir 40/69 (**50 rows changed**). ✅ `obesity` was already clean at 78/31
with **0 rows changed** — which is what makes the other two trustworthy rather than a blanket
transform that would have "fixed" everything indiscriminately.

⛔ **Comma decimal**: one `hba1c` is `"5,9"` and becomes NaN silently under `pd.to_numeric`.
Repaired to 109 usable. GSM8289030 genuinely lacks HDL/LDL — 108 usable for those two only.

### Guards

**Checksum: demonstrably live.** Re-hashed all 968 files against ENA md5s (968/968 OK). Before
reporting that pass, flipped byte `bf`→`40` at offset 83,888,407 of a copy, asserted the sha256
changed and size was preserved, and required the guard to fire — clean `OK`, corrupted
`MD5_MISMATCH`.

**Reproduction, not digest: demonstrably live.** Rebuilt the 109 x 14,078 matrix from the 109
retained RSEM outputs — bitwise identical, 0 cells differing, all 11 PAR_Y duplicate-sum groups
recovered; perturbing one cell by 1.0 makes it fail.

⚠ The prior chain's audit receipts assert their own counts; everything above was recomputed from
FASTQs, RSEM outputs and the GEO deposit directly.

---

## The P98 delta is material — and it exposed that Fig 3D is run-level — 2026-08-30

`executions/labelmap_correction_delta_20260830T210030Z/`. Jobs 21250896 / 21250897.
**Nothing adopted**: `donor_fstage_documented.tsv` still `9548296f…`, `docs/` unmodified, no consumer
repointed.

| Fig 3D hepatocyte carrier share | frozen (ledger) | P98-corrected | Δ |
|---|---|---|---|
| F0→F1 | **14.5%** (607/4180) | **19.8%** (829/4191) | **+5.26 pp** |
| F1→F2 | **25.8%** (943/3662) | **22.2%** (814/3660) | **−3.51 pp** |
| F2→F3 | 20.2% | 20.2% | 0.00 |
| F3→F4 | 54.3% | 54.4% | +0.12 pp |

⛔ **The narrative shape changes.** `14.5 → 25.8 → 20.2 → 54.3` (rise, dip, jump) becomes
`19.8 → 22.2 → 20.2 → 54.4` (flat ~20-22% plateau, then jump). **The F1→F2 peak disappears.** The
F2→F3 lineage clause survives verbatim (Endo 28.1→27.9, Chol 27.6→27.5, Fib 24.1→24.4, Hep 20.2).

✅ Both arms proven to differ by sha256 at every stage, and the frozen arm reproduces the sealed
117-program artifact to max |Δβ| **1e-15** and the released hep shares exactly. Internal control:
512's label-independent `primary_effects` is bit-identical, max |Δβ| **0.000e+00**.

### ⛔ The real finding: one donor moves it 5 pp because the chain is RUN-LEVEL

`343e` sets `n_donors = len(sub)` over per-(sample, dataset) rows → **8/8/12/11/19 = 58 runs**.
`348b` then prints `[Hepatocytes] donors: 58` and treats every run as an independent donor.
**The true donor tally is 6/8/12/11/9 = 46** (verified against `donor_pairing.csv`: 46 donors,
67 runs). So P98 is 1/8 of the F0 arm instead of 1/6, and **part of the +5.26 pp swing is
pseudoreplication, not label correction.**

This project already knew — `project-fstage-58-is-runs-not-donors-2026-08-08` — and the chain still
runs at run level. `512_hotspot_v2_donor_refit.R` IS donor-collapsed; `348b` / `WS3_01` are not.
**User decision 2026-08-30: fix the run-vs-donor collapse first, then apply P98 on top.**

### ⚠ The correction perturbs contrasts P98 is not in

`348b` fits one joint `~ F_factor + dataset + sex` model, so eBayes moderation shifts everywhere.
Endothelial F4_vs_F3 significant genes **1456→1352**; T cells F4_vs_F2 **672→555**. Across all
604,572 shared rows **zero have identical logFC** (max |ΔlogFC| 2.54, median 0.017). Any "only the
F0/F1 arms move" assumption is wrong.

### Chain B, run-level, both labels

117 program betas (`score ~ F_stage_documented`, n=46 donor-collapsed): |Δβ| median **0.0064**, max
0.0778, **0/117 unchanged**, **0/117 change sign**; median |β| 0.109 so the median relative shift is
**5.6%**. BH q<0.05 **32→35** (3 gained, 0 lost, all borderline). HC3 q<0.05 **9→9**.

⛔ `crn_transition_scores.tsv` moves hardest: only the two F0/F1 transitions change (F2→F3 and F3→F4
bit-identical), but **5 gained / 4 lost** at q<0.05 and individual cells swing far —
macrophages_16 F0→F1 z **7.13 → 2.76** (q 0.00004→0.045), macrophages_13 z −4.79→−2.58
(q **0.00006 → 0.060**, loses significance), macrophages_6 q **0.00045 → 0.096**.

### ⛔ Two pre-existing breakages, unrelated to the label

1. **Current-HEAD `512_hotspot_v2_donor_refit.R` + `lib_donor_collapse.R` fail on unchanged canonical
   inputs** — `Pooled-cell donor(s) missing metadata: SRR17375011, …`, confirmed with an
   all-canonical control run. `lib_donor_collapse.R` went 8,037→9,822 bytes on **Aug 13, after the
   Aug 7 seal**. Chain B ran on sealed-era versions from git (`44bd7c9`, `3d99022`), which is why the
   frozen arm reproduces the sealed artifact exactly. **Needs its own fix.**
2. **WS3_01's coarse axis no longer reproduces its released values** — SH **0.483→0.376**, Steatosis
   **0.256→0.192** — because the canonical `pseudobulk_de/` tables drifted since Jul 2.

---

## Fig 3D at donor level: pseudoreplication was SUPPRESSING the effect — 2026-08-30

`executions/fig3d_donor_level_20260830T214514Z/`. Nothing adopted; no canonical file touched.
n = **46 biological donors** (F0 6 / F1 8 / F2 12 / F3 11 / F4 9 frozen; 5 / 9 corrected).

⚠ The model is `~ F_factor` **with no covariates** — `dataset` is constant (GSE202379 only) and
`sex` is 100% NA here, so 348b's `has_dataset` / `has_sex` switches are both FALSE. The eBayes point
stands regardless: one joint fit per cell type, so every contrast moves.

| unit | label | F0→F1 | F1→F2 | F2→F3 | F3→F4 |
|---|---|---|---|---|---|
| run (n=58) | frozen | 14.522 | 25.751 | 20.168 | 54.256 |
| run (n=58) | corrected | 19.780 | 22.240 | 20.168 | 54.374 |
| **donor (n=46)** | **frozen** | **15.224** | **26.085** | **20.649** | **54.414** |
| **donor (n=46)** | **corrected** | **28.371** | **22.295** | **20.808** | **54.197** |

⛔ **My hypothesis was backwards. Pseudoreplication did not inflate the label swing — it suppressed
it by 7.89 pp.** Run-level swing +5.259 pp; **donor-level swing +13.147 pp**, 2.5x larger. At donor
level P98 is 1/6 of the F0 arm rather than 1/8. The *level* effect on the frozen label alone is small
(+0.703 pp), so the ledger row barely moves — **what changes is its sensitivity to one donor.**

⚠ **And +13.15 pp is not "hepatocytes gained".** At donor level the macrophage F0 group is exactly
{P98, PHL1, PHL2} = **3 donors, sitting exactly on `MIN_DONORS_PER_STAGE = 3`**. Moving P98 leaves 2,
so 348b **drops the entire macrophage F1_vs_F0 contrast** and macrophages leave the denominator
(24.83% → 0.00%). Like-for-like on the lineage set common to both arms the swing is **+8.118 pp**
(20.252 → 28.371); ~5.03 pp is dropout. **This number hinges on one donor crossing an arbitrary
floor.**

F2→F3 lineage, donor level (frozen → corrected): Endo 29.313→28.809, Chol 26.527→26.713,
Fib 23.511→23.671, Hep 20.649→20.808. Macrophages and T cells contribute 0 at F2→F3 in all four arms.

### ⛔ `crn_transition_scores.tsv` was run-level, and most of its significance was pseudoreplication

**184 run units → 64 donor units.** q<0.05 rows **176 → 42** (frozen) and **177 → 39** (corrected).
F0→F1 group sizes 13/11 → **6/7**.

- macrophages_16 F0→F1: run z 7.128→2.760; **donor z 5.362→2.246, q 0.0000→0.5172 — significance
  lost outright.**
- macrophages_13: run q 0.0001→0.0600; **donor q 0.0213→0.5172.**

⛔ 3 of 64 donors (GSE244832 D01/D02/D06) had runs whose `F_stage_inferred` disagreed (D02 runs
assigned 3,3,1,2,2,3,3) — **the same donor was previously entered into two different transition
groups at once.** Resolved by majority vote at collapse.

✅ **512's `documented_fstage_sensitivity.tsv` was already donor-collapsed** (line 694 hard-asserts
58 records → 46 donors), so the 117 betas are unchanged: median |Δβ| 0.006384, max 0.077833,
**0/117 sign flips**, q<0.05 32→35.

### ✅ `lib_donor_collapse.R` HEAD is an improvement, and the blocker was a stale artifact

HEAD (`e85072df`) is a **strict superset** of the sealed lib (`57eedf02`): +48 Liver_Atlas runs → 19
donors, **0 of 225 shared run ids disagree**, GSE202379 sub-map byte-identical.

⛔ The real blocker: `hotspot_modules/donor_collapse/donor_scores_all_weighted.tsv` is dated
**2026-07-12**, the lib changed **2026-08-13**. 11 of its 85 "donor" ids are raw Liver_Atlas run ids,
so 512's `assert_joinable` correctly refused. `511_pooled_cell_donor_scores.py` on disk already
mirrors the HEAD map — **only its output was never regenerated.** Rerunning 511 unchanged fixes it.

⛔ **A second defect found:** today's `512_delta.R` had a hardcoded `PRIMARY_SCORE_FILE`, dropping the
canonical `Sys.getenv` override — so the first HEAD attempt silently re-read the stale artifact and
reproduced the original error. Restored to the canonical form.

✅ **The lib choice is immaterial to every documented-F-stage number:** sealed-lib+stale-scores vs
HEAD-lib+regenerated-scores give **max |Δbeta| = 0.000e+00, max |Δq| = 0.000e+00**, identical
n_donors on all 117 programs.

### Guards

Run-level arms re-run as a reproduction guard: `routing_summary_by_transition.csv` and 510's
`crn_transition_scores.tsv` **bit-identical** to today's for both labels; routing categoricals
identical with numeric drift ≤3.3e-13 and 348b logFC ≤3.6e-13 (BLAS noise; the label effect is
max 2.54). Four arms proven to read four distinct (level, input-sha) pairs.

⚠ **A ragged-file guard fired on a FALSE alarm:** R's `strsplit` drops trailing empty fields, so a
row ending in NA looks one field short. Fixed to count separators instead (`awk` confirms uniform
NF=72). The mirror of the usual failure — a guard that fires when nothing is wrong.

---

## ⭐ RNA→chromatin RECURS in GSE296875 — a second cohort, a different assay — 2026-08-30

`executions/gse296875-crossmodal-recurrence-20260830T193000Z/`. Job 21251685, 3:54:48, 1.32 GB.
**39 donors**, 68,398 same-nucleus RNA+ATAC nuclei pseudobulked to donor level, 1,000 permutations.

| | GSE267145 | GSE296875 |
|---|---|---|
| assay | bulk RNA → bulk **H3K27ac** | same-nucleus RNA → **ATAC** |
| unit | 99 participants | **39 donors** |
| skill | **+0.13345** | **+0.12985** |
| permutation p | 0.000 | **0.000** |
| folds positive | 5/5 | **5/5** (+0.0439 to +0.1711) |
| null max | −0.00376 | **−0.00430** |
| pre-committed outcome | **A** | **A** |

Donor jackknife (df 38, 39 units): **CI95 [0.0981, 0.1630]**, se 0.0166. The two point estimates
differ by **0.0036**, well inside that interval.

⭐ **This is recurrence across instruments, not replication.** Different assay (H3K27ac vs ATAC),
different resolution (bulk vs same-nucleus), different cohort, different n. That makes it *stronger*
evidence that the relationship is real than a like-for-like repeat would be — and it licenses **no
shared effect size**. Report the two separately.

✅ **Rank 4 was transferred from the GSE267145 lane, not fitted here** (`rank_source: "transferred
from the GSE267145 lane"`). A descriptive sweep shows skill still climbing at rank 8 (0.1474), and
the file marks it `descriptive_rank_sweep_is_not_a_selection: true`. **The reported number uses the
pre-committed rank, not the best one.**

### Two secondary arms, both surviving

- **Leave-one-well-out**: skill **+0.14206**, **7 of 8 wells positive**, jackknife CI95
  [0.1078, 0.1778]. So the effect is **not well batch structure** — the obvious confound in a
  multiome with 8 wells.
- **Hepatocyte-only**: skill **+0.10299**, 5/5 folds, CI95 [0.0669, 0.1410]. Survives restriction to
  one lineage, so it is not purely composition.

### The gate is demonstrably live, and its threshold is reported

`clean_gate_passed: true` on 200 destroyed-pairing permutations (max **−0.00226**, mean −0.0916),
and **all 200 permutations actually moved rows** (mean 38.08 of 39 donors relocated) — the
permutation was verified to permute. `leaky_gate_fired: true` on the identity permutation.
⭐ The file also records `smallest_p97_5_threshold_that_would_fire_on_the_leaky_arm: 0.1298`, i.e.
**how much margin the gate had**, rather than only that it passed.

### The barcode trap was real and was avoided

⛔ **2,784 bare barcodes are shared across wells** (up to 3 wells sharing one), and a bare-barcode
join would have **mis-paired 5,625 nuclei**. Joined on `(well_id, raw_barcode)`, verified unique:
68,398 unique pairs for 68,398 nuclei against only 65,557 unique bare barcodes.

✅ Join verified by head-to-head cosine, which **can fail and was shown to**: correct-partner median
cosine **0.9825** vs wrong-partner **0.0333**, with **0 nuclei where the wrong partner wins** and 0
below 0.5. The tamper relocated 58,068 nuclei and was asserted not to be a no-op.
✅ Reproduction guard: 8,308,427 values recomputed across 227 units, **0 disagreeing**, and a
deliberate 1-value tamper was detected.

⛔ **`controlled_data_policy_dataset_used: true`** — this run used GSE296875, so any weight release
from it needs written permission. The GSE267145 result is separately affected. An atlas-only model
remains the permissive case.

---

## Rank sensitivity: the RNA→H3K27ac result is not rank-specific — 2026-08-30

`executions/gse267145-crossmodal-ceiling-20260830/results_ranksweep_20260830T220000Z/`.
Job 21251831, 1:42:32. Eight ranks, each with its own one-sided 100-permutation destroyed-pairing
gate. ✅ **Reproduction guard: rank 4 returned +0.13345 to within 3.4e-06** of the reported value, so
this harness agrees with the run it extends.

| rank | skill | p | gate | folds+ | null max | retrieval top-1 |
|---:|---:|---:|:--|:--|---:|---:|
| 1 | +0.06326 | 0.000 | ok | 5/5 | +0.01738 | 0.0404 |
| 2 | +0.10570 | 0.000 | ok | 5/5 | +0.01334 | 0.0909 |
| 3 | +0.11749 | 0.000 | ok | 5/5 | −0.00710 | 0.1313 |
| **4** | **+0.13345** | **0.000** | **ok** | **5/5** | −0.01888 | 0.1616 |
| 6 | +0.14674 | 0.000 | ok | 5/5 | −0.02418 | 0.2121 |
| 8 | +0.15693 | 0.000 | ok | 5/5 | −0.00765 | 0.2323 |
| 12 | +0.16900 | 0.000 | ok | 5/5 | −0.03716 | 0.2929 |
| 16 | +0.17490 | 0.000 | ok | 5/5 | −0.04743 | 0.3131 |

**8 of 8 gates pass; every rank gives p = 0.000 and 5/5 folds positive.** The finding does not
depend on the choice of 4.

⛔ **Rank 16 scores highest (+0.17490) and reporting that is FORBIDDEN.** Rank 4 was fixed in advance
for a structural reason (components 5 and 6 fail the pairing null and a top-k method cannot step over
them); picking the sweep maximum afterwards is selection on the outcome. The prohibition is written
into the script's own output beside the winning rank so it cannot be made by accident.

⚠ **Skill is still climbing monotonically at rank 16, and so is retrieval top-1 (0.0404 → 0.3131).**
The same pattern appears independently in GSE296875, still rising at rank 8 (0.1474). So **the
reported +0.13345 is conservative** — there is more shared structure than rank 4 extracts. Whether
the later components are trustworthy is a separate question the pairing null already answered NO for
components 5 and 6 in the ceiling analysis, which is exactly why 4 was chosen. **Do not resolve that
tension by raising the rank; resolve it by re-deriving the ceiling if it matters.**

---

## Corrected DNA-LM embeddings complete, and the repair is confirmed downstream — 2026-08-30

Job 21250329, 1:33:20 (both models), `variant_program_substrate-20260830T000000Z/embeddings/`.

| model | variants | dim | tokens/seq | strands | context fed |
|---|---|---|---|---|---|
| nucleotide_transformer | 30,667 | 1,280 | 1,000 | forward | 5,994 bp |
| **hyenadna** | 30,667 | 256 | 6,001 | **both** | 6,000 bp |

**All 8 assertions green on both**, both bound to prespecification `8037d00b…`, 10 windows dropped
for non-ACGT (the two known assembly-gap regions). ⭐ **hyenadna completed for the first time** — it
had never finished in any prior attempt.

✅ **The indel repair verified independently, measured in the embeddings rather than the windows:**
comparing against the pre-repair run at identical row order, **exactly 33 rows changed their delta —
all 33 indels, 0 non-indels**, out of 30,667. That reproduces the repair agent's count from the
downstream side, which is a different measurement than the one it made.

The leverage finding is unchanged after repair: indel delta L2 median **1.2977** vs SNV **0.0874**
(**14.8x**), and **91 of the 100 largest deltas** are still indels (0.31% of rows). All finite.

⛔ **These embeddings are usable; the LABELS in this lane are not.** 1,215 rows at those same 33
variants in `tables/variant_gene_training_long.parquet` still carry the wrong `beta_alt` sign,
because `32_reorient_to_reference.py` used the same broken gate. **Fit heads from
`variant_program_substrate_rebuild-20260830T160000Z`, which has these right.**

---

## S-LDSC of the 117 programs: a clean negative on the coefficient — 2026-08-31

`GWAS/ldsc/program_heritability_20260830/results/`. Job 21247851, 11:05:47. 585 regressions
(117 programs x 5 traits), baselineLD v2.2 + a lineage-matched background annotation.

**0 of 117 programs pass BH on the coefficient (tau) for ANY trait.**

| trait | tau BH<0.05 | enrichment BH<0.05 | max tau_z |
|---|---:|---:|---:|
| ukbb_alt | **0** | 2 | +1.695 |
| ukbb_ast | **0** | 14 | +1.556 |
| pdff | **0** | 0 | +2.788 |
| ghouse_cirr | **0** | 0 | +1.299 |
| finngen_obesity (specificity control) | **0** | 2 | +2.396 |

The largest tau_z anywhere is **+2.788**, on pdff, which is itself a `secondary` trait (h2 z 6.35,
below the recommended bar of 7). Nothing survives multiplicity.

### ⭐ The run is its own demonstration that enrichment is the wrong statistic

The 18 "enrichment hits" refute themselves on inspection:

- **5 of the 18 have enrichment < 1** — that is **depletion** counted as a pass, and one
  (`finngen_obesity macrophages::15`) has a *negative* point estimate of **−0.1998**. `hepatocytes::4`
  on ukbb_ast has enrichment **0.266 at BH p = 9.7e-14**, the most significant number in the table,
  and it means the annotation carries *less* heritability than its SNP share.
- **tau_z for the 18 ranges from +1.56 down to −3.06.** Two point significantly the *opposite* way
  from their enrichment. None reaches BH significance on tau (min BH p 0.965).
- The **background** annotations show the same: `background::hepatocytes` on **finngen_obesity — a
  deliberate non-liver specificity control** — has enrichment p **3.4e-4** with tau_z 0.83. A large
  annotation is enriched almost automatically.

✅ This is why the pre-registered read was `Coefficient_z-score`, not `Enrichment`. Had the headline
been taken from enrichment, this run would have produced 14 "significant liver programs" on AST and
2 on an obesity control.

### ⛔ What this negative does and does not mean

The run's own `HEADLINE.txt` states it: **no MASLD disease endpoint clears the heritability bar**, so
these results answer a question about **serum ALT and AST**, not about MASLD. `ukbb_alt` (h2 z 16.4)
and `ukbb_ast` (12.3) are the only properly powered traits; `pdff` (6.35) and `ghouse_cirr` (6.29)
are labelled `secondary` and below the bar of 7.

So: **the 117 programs do not carry detectable liver-enzyme heritability enrichment beyond a
lineage-matched background.** A null on the disease endpoints would be uninformative — those GWAS
are too small to have shown anything either way. This is `tested_negative` for the enzyme traits and
`untestable` for the disease endpoints, not a single negative verdict.

⚠ The programs were derived from single-cell expression, not from genetics, so this is not a
contradiction — it says common-variant heritability for liver enzymes is not concentrated in them.

---

## ⛔ Intervention response: NEGATIVE — the apparent tracking is the treatment arm — 2026-08-31

`executions/gse83452-intervention-response-20260831/`. **60 paired biopsies one year apart, 32 diet /
28 surgery**, 171 participants, 231 samples. The only human intervention substrate in this project.

**Every model in this campaign is scored on association with a histology grade. This is the test that
asks whether the score MOVES when the same person is treated. It does not.**

⛔ **First, a scope limit found by reading the raw deposit rather than assuming.** Parsing all 231
`characteristics_json`, the only histology field at either visit is `liver status`
(NASH / no NASH / undefined). **There is no NAS, fibrosis, steatosis or ballooning at follow-up.**
So sign is testable and **magnitude is not** — no magnitude claim is made anywhere in this lane.

⛔ **No fitted model had ever produced a follow-up score.** The existing prediction bundle is stamped
`timepoint: baseline`, `followup_arrays_read: false`, `n_predictions: 152`. The scoring everyone
relied on had never touched the second biopsy.

| readout | pooled ρ | **surgery** | **diet** | within-arm p | diet-only effect |
|---|---:|---:|---:|---:|---:|
| gene_median_elastic_net | +0.123 | **+0.503** | **−0.366** | 0.770 | +0.059 |
| gene_median_linear_svm | +0.225 | **+0.492** | **−0.355** | 0.888 | +0.066 |
| gene_median_pca_en | +0.087 | **+0.458** | **−0.351** | 0.571 | +0.046 |
| per_array_rank_en | +0.152 | **+0.380** | **−0.204** | 0.555 | +0.061 |
| continuum | +0.258 | **+0.391** | **−0.125** | 0.374 | +0.416 |

⭐ **The structural fact that decides it.** Within the 37 baseline-NASH participants, surgery is
**14 resolved / 1 persisted** and diet is **6 / 16**. *Resolution is nearly a restatement of arm.* So
the pooled contrast cannot separate "tracks histology" from "responds to surgery". **Diet is the only
stratum that can separate them — and there the effect REVERSES SIGN in 5 of 5 readouts.** Holding arm
composition fixed in the null flattens everything to p 0.374–0.888.

⛔ **The surgery arm is structurally untestable: its contrast rests on one participant** (1 persisted
of 15). The pre-specified concentration diagnostic is unusable for the same reason — it credits
51–82% of the effect to surgery, whose contrast is that single person. A within-arm permutation
replaces it.

**Verdict: negative and underpowered.** ⛔ **CORRECTED 2026-09-02: the range "0.138–0.215" appears
NOWHERE in the producing file and is retracted.** `results.json` carries `P2_mde_80pct` per arm:
**0.1001** (gene_median_elastic_net), **0.1206** (linear_svm), **0.1024** (pca_elastic_net),
**0.1766** (per_array_rank_elastic_net) — so the correct range for the comparable arms is
**0.1001–0.1766**. ⚠ The logit variants run 0.4168–0.7592 and `continuum_fixed_projection` is
**1.1673**, so any quoted range must name its arms.
✅ **The verdict is unchanged and in fact TIGHTER**: observed effects of 0.046–0.066 sit below even the
smallest MDE (0.1001), so this cohort cannot exclude an effect of that size — **and the observed sign
is wrong.** The claim that these models track intra-participant histological change is unsupported, and
the pooled appearance of tracking is the treatment arm.

### ⛔ Corrections to earlier work in this campaign

- **The regression-to-the-mean argument does not transfer between readouts.** The −1.852 / −2.292
  baseline imbalance is a **continuum-axis property**; the four probability models are *balanced* at
  baseline. What does replicate is r(baseline, delta) = **−0.518 to −0.584 on all five**. Adjustment
  leaves 42–77% of each effect, all intervals spanning zero. ⭐ **Check a confound on the specific
  readout you are using, not on a sibling that shares the cohort.**
- **A tie-structure ceiling was computed with `argsort` instead of its inverse** — correct on sorted
  input, **0.05 instead of 0.849** on real data. Caught before the primary.
- **A tamper test fired for only 1 of 4 models** because it perturbed a gene outside three models'
  1,000 selected features. Replaced with one gene each model actually uses; then all four fired.
- ⚠ **The specificity null is weak and is reported as such.** Gene-scramble p = 0.000–0.003 looks like
  a pass but is not scale-matched: scrambled models produce near-zero deltas, so it mostly tests
  whether the model moves at all. Same class of error as the timepoint-scramble it replaced. **It
  carries no weight in the verdict.**

✅ Reproduction: all four models reproduce the frozen `predictions.tsv` to **1.11e-16**; continuum
deltas match the 2026-08-30 run at **max abs diff exactly 0**. Matrix confirmed as the frozen oligo
path, not SCAN. ⚠ The five readouts are **not independent** — mean off-diagonal ρ 0.489, and the
three `gene_median` models sit at 0.71–0.94, so "5 of 5 reverse" is closer to 2–3 independent looks.

---

## ⛔ DNA-LM embeddings contribute NOTHING to variant→gene prediction — 2026-08-31

`executions/dnalm_heads-20260831T183647Z/`. Prespecification `cba90896f5ef…` written before any fit.
Join **29,641 of 29,651** label variants (the 10 missing are exactly the receipt-dropped assembly-gap
windows). Labels read only from the rebuild. Paired bootstrap, 2,000 resamples of the **55 held-out
LD blocks**, identical draws for every method.

| task | baseline | best embedding head | paired delta |
|---|---|---|---|
| \|z\| | distance ridge **0.2344** | hyenadna delta+dist 0.2126 | **−0.0218 [−0.0329, −0.0100]** *worse* |
| \|z\| | distance ridge 0.2344 | NT delta+dist 0.2049 | **−0.0294 [−0.0421, −0.0162]** *worse* |
| \|z\| | distance ridge 0.2344 | 128-PC + dist interaction | +0.0001 / +0.0005, **both cross zero** |
| \|z\| | raw distance 0.2309 | embedding **alone** | **−0.2345 / −0.2409** |
| signed beta (AUC) | 0.5030 | 0.5062 | −0.0015 / −0.0014, **cross zero** |
| hub top-1 | nearest-TSS 0.2488 | 0.2400–0.2438 | −0.0064 [−0.0169, +0.0031], **crosses zero** |

**No configuration beats the baseline anywhere.** Full-dimension heads are measurably *worse*; the
only configuration that does not lose (heavy PCA plus regularisation) reproduces the distance
baseline to within ±0.002 — i.e. it recovers the baseline by discarding the embedding.

### ⭐ G7 is the result, not a failed guard

| guard | intact | after permutation |
|---|---|---|
| **G6** permute the **labels** | +0.2126 / +0.2049 | **−0.0800 / −0.0549** — harness detects signal loss |
| **G7** permute the **embeddings** | +0.2126 / +0.2049 | **+0.2269 / +0.2203** — score goes **UP** |

⭐ **Destroying the correspondence between each variant and its embedding did not hurt the model.**
Read against G6, which shows the harness *can* detect signal loss, that is decisive: **the entire
score comes from the distance block.** The embedding is not weak here, it is inert.

Embedding-only heads confirm it: |z| at **−0.0001 / −0.0065**, signed effect AUC **0.4951 / 0.4996**,
and hub top-1 **0.0908 against a within-hub null of 0.0911 ± 0.0035** — structurally, since a
variant-level feature is constant across every gene in a hub.

⚠ hyenadna does beat NT on |z| (**+0.0077 [+0.0048, +0.0107]**, excludes zero), so the MPRA ordering
reappears — **but both sit below the baseline, so it ranks two losers.**

### ⛔ A join defect that would have silently inverted 18,730 variants

The npz `variant_id` is `chrom:pos:allele1:allele2` in **source-panel order**, while the windows were
built from reference-gated `ref_allele`/`alt_allele`. **These disagree on 19,756 of 30,667 rows.**
Separately, **1,026 positions carry two rows under mirrored ids whose arrays are bitwise identical**,
so `chrom:pos` alone is ambiguous.

⛔ Joining on the id string would have looked **~37% successful** and then inverted `alt_minus_ref` on
**18,730 variants**. A naive row-order join would have mis-assigned **21,788 of 29,641**. The correct
join (chrom:pos + reference-gated alleles) gives orientation **+1 on 29,641/29,641, zero swapped**.

### Guards, each shown able to fail

G1 identifier join 0 → 214 mismatches under a 1% key permutation. G2 positional join **21,788**
mis-assignments (informative, not vacuous). G3 orientation moves by exactly 500 when 500 alleles are
swapped. **G4 LD split: 0 violations across 83,477 cross-fold pairs, max cross-fold r² 0.0099948 —
the constraint is binding, not slack — against 8,687,237 under a by-variant split.** G5 reproduces
all six recorded baselines to **0.00000**. Gram reduction matches a literal X'WX to 9.65e-15;
a tampered centring blows it to 3.46e3.

### Caveats stated by the producer

⚠ Ridge is a **linear** head, so this bounds linear readout only. ⚠ Embeddings are mean-pooled over
1,536 bp and per-token states were never written, so **no positional head was testable**. ⚠ Paired
intervals resolve to about ±0.013, so a real gain above ~+0.015 would have been visible.
⚠ `e_refalt` looked like +0.0353 on fold 4 but flips sign across folds (+0.0145, −0.0262, −0.0179,
+0.0132, +0.0353) — **noise, killed by the 5-fold view before the interval spoke.**

---

## The signed-eQTL negative resolved: `tested_negative`, but the licensed claim is NARROW — 2026-08-31

`executions/eqtl_sign_reliability-20260831T195702Z/`. Prespecification `5b7e5a80fd9e78fb…` hashed
before computing. Held-out fold 4: 111,426 rows / 5,930 variants / **55 LD blocks**.

### World B is refuted — the recorded signs are not coin flips

`r = Phi(|Beta|/SE)` per row: **mean 0.7818**, median 0.7891, below 0.6 on only 17.5% of rows.

| ceiling | all | \|z\|≥1 | ≥2 | ≥3 | ≥5 |
|---|---|---|---|---|---|
| flat prior | 0.7878 | 0.9429 | 0.9949 | 0.9998 | 1.0000 |
| empirical Bayes (π₀ = 0.1492) | 0.6375 | 0.7534 | 0.9073 | 0.9882 | 1.0000 |
| full oracle (MC) | 0.8774 | 0.9878 | 0.9999 | 1.0000 | 1.0000 |

✅ Chance was **measured**, not assumed — within-block label permutation, 2,000 draws, 0.4998–0.5019
across every band, so sign is not clustered in LD blocks.

### Every stratum sits at chance, against ceilings up to 1.0000

| stratum | held-out rows / blocks | best embedding AUC | measured chance |
|---|---|---|---|
| all | 111,426 / 55 | 0.5016 [0.4901, 0.5133] | 0.4998 |
| \|z\|≥1 | 47,777 / 54 | 0.5019 [0.4830, 0.5199] | 0.5003 |
| \|z\|≥2 | 17,113 / 49 | 0.4985 [0.4822, 0.5143] | 0.4999 |
| \|z\|≥3 | 8,862 / 39 | 0.4940 [0.4581, 0.5344] | 0.5003 |
| \|z\|≥5 | 3,847 / 31 | 0.4900 [0.4688, 0.5121] | 0.5019 |

**Every interval crosses zero headroom.** Every stratum cleared its power gate from counts *before*
any AUC was read.

⭐ **The embeddings stay inert at high confidence** — the number this test existed to produce.
Permuting the variant→embedding map is **≥ intact in 9 of 10 model × stratum cells**, and in several
the paired delta *excludes zero in the wrong direction*: on |z|, hyenadna **−0.0152 [−0.0235,
−0.0068]**, NT **−0.0105 [−0.0189, −0.0020]**. Permuting does not merely fail to hurt — it helps.
**So the original result is not a label-noise artefact.**

✅ **Positive control: the representation is not broken.** The pooled delta recovers *substitution
class* on held-out folds at NT median AUC **0.9444**, hyenadna 0.8304, both above all 40 shuffle
draws. hyenadna sits at chance on exactly the four self-complementary substitutions (A>T, T>A, C>G,
G>C) — the expected signature for a both-strands model.

### ⛔ Two corrections to what this file previously said

1. **The numbers were misattributed.** `0.5033 [0.4855, 0.5269]` is the **allele-feature ridge**
   (substitution dummies, EAF, indel length, CpG flags) — a *tabular* model, not a DNA-LM head. The
   DNA-LM heads' own values are **0.5016** (NT) and **0.5015** (hyenadna). And **0.5016 is the
   majority-class fraction, not an AUC floor** — chance AUC is 0.5 whatever the class balance.
   Earlier text in this file treating 0.5016 as a floor was wrong.
2. ⛔ **A third world exists that neither the brief nor this file had named, and it bounds the
   claim.** The eQTL labels are **marginal association statistics with no fine-mapping column**.
   Under one causal variant per (gene, LD block), **at most 1.9–4.3% of held-out rows can be causal**,
   so the ceiling for *any* sequence model is **0.5055–0.5211** — inside the measurement resolution
   unless there are ≥2 (at |z|≥2) to ≥5 (all) independent signals per locus.

### ⭐ The sentence this licenses, and the one it does not

✅ **Licensed:** *on marginal, un-fine-mapped cis-eQTL associations from this panel, DNA-LM embeddings
carry no directional information, and the design cannot separate a model failure from a target that
is not sequence-determined at the variants tested.*

⛔ **NOT licensed:** "sequence models cannot predict eQTL direction." Most variants in a marginal
association are tagging rather than causal, so the task as posed partly asks a sequence model to
predict the direction of **LD tagging**, which is not sequence-determined at all.

**The experiment that separates them: the same head on eQTL-fine-mapped credible-set variants.**
The GWAS credible sets supplied the variant list, but the eQTL side was never fine-mapped.

⚠ One instability worth noting: `b1_dist`/`m1_dist_emb` fall **below chance in all five folds at
|z|≥3**, driven by fold 0 (0.3170 at ≥3, 0.1192 at ≥5) — the upstream/downstream indicator is
fold-unstable. `e_emb` never shows it and straddles chance everywhere, which is why the
embedding-only negative is the clean one.

⚠ **G6 was vacuous on the signed target and is reported as such**, not as a pass: intact is already
at chance, so there is nothing for a label permutation to collapse. It is live on |z|
(+0.3445 → −0.2063).

---

## RNA→chromatin is DISEASE-DEPENDENT — but it is normal-vs-diseased, not graded severity — 2026-08-31

`executions/crossmodal-severity-stratification-20260831T190048Z/`. `PRECOMMIT.md` byte-unchanged at
`fbb72586…`. Δ = difference-in-differences on the paired per-donor home advantage
`h = s_own − s_other`, against a label-permutation null that refits at the observed group sizes.

| arm | pre-reg | n low/high | Δ | z | perm p | MDE | verdict |
|---|---|---|---:|---:|---:|---:|---|
| GSE267145 activity≥3 | **yes** | 49/50 | **+0.2087** | +4.81 | **0.0000** | 0.1239 | **DISEASE_DEPENDENT** |
| GSE267145 activity≥4 | no | 64/35 | +0.1457 | +3.48 | 0.0050 | 0.1154 | DISEASE_DEPENDENT |
| GSE267145 female-only ≥3 | no | 49/36 | **+0.2482** | +5.52 | **0.0000** | 0.1287 | DISEASE_DEPENDENT |
| GSE267145 **NOR dropped**, ≥4 | no | 40/35 | **+0.0595** | +1.06 | 0.24 | 0.1580 | **UNDERPOWERED** |
| GSE267145 female-only, NOR dropped | no | 34/27 | +0.1200 | +1.77 | 0.09 | 0.2014 | UNDERPOWERED |
| GSE296875 fibrosis_any | yes | 22/15 | +0.0548 | +1.30 | 0.18 | 0.1290 | UNDERPOWERED |
| GSE296875 steatosis≥5% | yes | 18/20 | +0.0161 | +0.38 | 0.70 | 0.1258 | UNDERPOWERED |

⭐ **The 2×2 is the headline: the map does not weaken across the severity boundary, it stops
working.** LOW donors score **+0.0487** under a LOW-trained map and **−0.0364** under a HIGH-trained
one; HIGH donors **−0.0072** vs **+0.1165**. Cross-stratum transfer lands at or below the trivial
mean-profile baseline.

✅ **Sex was chased and excluded.** Sex is a *complete* confound in the pre-registered split
(49F/0M low, 36F/14M high) — but female-only makes Δ **larger** (+0.2482), so sex is not the driver.
GEO submission order shows no stratum association (Mann-Whitney p = 1.00 RNA, 0.77 H3K27ac).

### ⛔ The load-bearing caveat: 24 normal livers carry the effect

Two arms share the same cut (≥4), the same 35 HIGH donors and the **identical** matched training size
m = [28,28,26,28,30]. The only difference is whether the 24 **NOR** (normal) livers sit in the LOW
pool. Δ falls **+0.1457 → +0.0595, a 59% attenuation**, and the LOW-trained model's skill on HIGH
donors *improves* from +0.032 to +0.071.

**The 24 normal livers are what break the transfer.** Within diseased liver only the residual is
positive in all three arms but never resolvable — MDE 0.158–0.201 against point estimates
0.059–0.120. ⛔ **Graded severity is `untestable` here, not negative.**

⛔ **Composition is not separable.** LOW is 24 normal + 25 NAFL; HIGH is 3 NAFL + 47 NASH. Bulk
RNA→bulk H3K27ac cannot distinguish "the map inside a cell type changed" from "the cell mixture
changed." ✅ **Write it as *the RNA→chromatin mapping differs between normal and diseased liver*.
Do NOT write within-cell-type rewiring.**

### ⚠ Part of the parent +0.13345 rides on non-biological structure

The rank-4 map splits cleanly in two:
- **Components 2 and 4** load on high-abundance, high-variance genes (SD percentile 88.1 / 91.2 vs 50
  for all genes), 39/150 and 38/150 in Hotspot programs, and are recognisable liver biology — CXCL2,
  GADD45G/B, NFKBIA, RHOB, HES1; JUNB, FASN, TRIB1, ETS2, TFF3, IER2. Top programs *Immediate-early
  AP-1 (FOS)*, *Activated fibroblast (MMP19)*, *Plasma lipoprotein (APOE)*, *Lipogenic hep-like
  (MLXIPL)*. Component 4 tracks activity at **ρ = −0.615**.
- ⛔ **Components 1 and 3 — including the LARGEST variance share (0.165) — load on low-abundance,
  low-variance features** (percentiles ~39 / ~37), **3/150 and 2/150 in any program**, topped by
  unannotated ENSG ids, lncRNAs and olfactory receptors.

In GSE296875 the biology-loaded component has ρ = +0.047 with steatosis; **no component there is
severity-linked.**

### Guards, and two method corrections they forced

Seven guards, each fired on a violation built for it (job 21288303): exact reproduction of both
parent numbers to **1.1e-16 / 8.9e-16**; ragged-row and unnamed-index TSV guards; endpoint-join loss
on a tampered copy (sealed sources verified untouched); the one-sided within-stratum gate; a
null-reference check on a content-free 80/19 split; and ⭐ **detection of a map difference as small as
48 of 96,460 peaks rewired (Δ +0.256)** — so the UNDERPOWERED verdicts above are measured negatives,
not blindness.

⭐ **Two corrections were forced by guards BEFORE any real severity split was read**
(`PRECOMMIT_AMENDMENT_01/02.md`, with failed jobs 21286015 and 21287171 as evidence):
1. The pooled donor mean of `h` returned **−0.046 on a content-free split**, because unequal group
   sizes reweight a constant model offset. Replaced by the difference-in-differences.
2. The bootstrap-around-zero was replaced by a label-permutation null, because **Δ has a non-zero
   null mean** — leave-fold-out shrinks each stratum's pool by its own test count.

⚠ **The producer pre-committed `constitutive` as the likeliest outcome and was wrong**, and says so.
⚠ One gate criterion (`null_max < 0.05`) failed on the HIGH-stratum arm at 500 permutations
(+0.0552, +0.0500) while `p97.5` passed comfortably (−0.0023, −0.0014) — **that criterion tightens as
N_PERM grows** and the parent lane used 100 draws. That arm was declared non-decisive regardless.

---

## The flagship variant→gene→program chain is dead at hop 1 — 2026-08-31

Recording this explicitly rather than leaving it implied across three separate results.

The plan's flagship lane was **variant → signed effect on every gene in its cis window → aggregate to
the 117 programs**. It carried a pre-registered trap: *"a program score is a weighted sum of its
member genes' betas, so if gene-level prediction is good, program-level follows by arithmetic"* — the
same defect as NAS = steatosis + ballooning + lobular. The plan therefore required the program-level
claim to beat a membership-shuffled null holding gene-level accuracy fixed.

**That test is now moot, because hop 1 does not clear.** Three independent measurements:

1. Honest LD-block splitting collapses gene-level prediction from Pearson **0.8669 to 0.0210**.
2. DNA-LM embeddings are **inert**: permuting the variant→embedding map *raises* the score, in 9 of
   10 model x stratum cells, while permuting labels collapses it.
3. Signed direction sits at **measured chance in every |z| stratum**, against ceilings up to 1.0000.

⭐ **Aggregating inert per-gene predictions cannot produce an informative program-level prediction**,
so the membership-shuffled null has nothing to adjudicate. The lane is reported as **gene-level
only, and negative** — exactly the fallback the plan pre-specified.

⚠ **Two things this does NOT say.** It does not say the programs are not genetically relevant — that
is a different question, tested separately by S-LDSC (0/117 on the coefficient) and now by the
reverse TWAS direction. And per the sign-reliability analysis, the licensed claim about hop 1 remains
narrow: **on marginal, un-fine-mapped cis-eQTL associations from this panel**, embeddings carry no
directional information, and the design cannot separate model failure from a target that is not
sequence-determined at the variants tested.

---

## ⛔ PENDING EDIT to `docs/PAPER.md` — a claim the evidence now contradicts — 2026-08-31

**Not yet applied.** `docs/` edits abort running campaign bundles (the resource-firewall collision),
and three jobs were live when this was found. Apply when the queue is clear, then re-run
`scripts/manuscript/validate_resource_scope.py`.

**Current text, `docs/PAPER.md` around line 431** (written 2026-08-30, before the endogenous results):

> Whether a pretrained model improves on a task-native baseline depends on the task. It does so for
> variant-to-regulatory-activity, **where the signal is local and sequence-encoded**, and it does not
> for tissue-level severity transfer.

⛔ **The clause "where the signal is local and sequence-encoded" is an explanation the evidence now
refutes.** It generalises a single result on a *synthetic reporter* assay into a mechanism, and the
endogenous test says the opposite:

- **MPRA (reporter):** hyenadna +0.0820, NT +0.0537, both 5/5 seeds — but an allele-identity ridge
  alone reaches **0.1342 of an absolute 0.2162**, so even there the task is dominated by which allele
  it is.
- **Endogenous variant→gene:** the same two models are **inert**. Permuting the variant→embedding map
  *raises* the score (+0.2126 → +0.2269) in 9 of 10 cells, while permuting labels collapses it.
  Signed direction sits at measured chance in every |z| stratum against ceilings up to 1.0000.

**Proposed replacement:**

> Whether a pretrained model improves on a task-native baseline depends on the task. On a synthetic
> reporter assay it does, by a small margin over a baseline that uses allele identity alone. On the
> endogenous variant-to-gene effect it does not: the embeddings are inert, in that permuting the
> correspondence between a variant and its embedding does not degrade the score. It also does not
> improve tissue-level severity transfer. These are statements about measured transfer on these
> cohorts, not about model quality in general and not about clinical use.

⚠ **The endogenous claim must carry its own scope limit** wherever it appears: the eQTL labels are
**marginal, un-fine-mapped associations**, and at most 1.9–4.3% of held-out rows can be causal, so
the design cannot separate a model failure from a target that is not sequence-determined at the
variants tested.

### Two other claims to re-check in the same pass

1. *"The limits on tissue-level severity transfer are stated as measured bounds rather than as failed
   attempts"* — still true and now better supported (routing oracle +0.0071 against a pre-declared
   floor of 0.0930), but the **intervention** result adds a second, different limit: severity models
   do not respond correctly to treatment. Consider whether that belongs here.
2. *"Model scores ... reported beside a per-cohort record of what each model saw during its own
   training. Where that record is unresolved the score is reported as uninterpretable"* — ✅ upheld
   and now realised in panel 7C, which draws TF-Sapiens' top macro-F1 (0.9389) hollow in control grey
   because its exposure is `unknown`.

---

## ⭐ Panel 7F: the top-scoring encoder CANNOT be cleanly compared — 2026-08-31

`executions/encoder-field-score-21314196/evaluation/`. Jobs 21313780 (common head, 20:25) and
21314196 (scoring, 1:13). **102 donors, 50,000 identical rows, 5 cell classes, study-held-out outer
folds, 10,000-resample donor cluster bootstrap with shared multiplicities across all models.**

**Seven feature blocks through one shared head** — the comparison nobody had actually run.

### Raw ranking, and why it is misleading

| model | family | exposure | clean / confounded / unresolved folds | macro F1 |
|---|---|---|---|---|
| geneformer_v2_316m::mlp | pretrained | target_label_unexposed | **0 / 3 / 2** | **0.9438** |
| tf_sapiens::mlp | pretrained | **unknown** | **0 / 0 / 5** | 0.9389 |
| scimilarity_v1_1::mlp | pretrained | **encoder_seen** | 3 / 2 / 0 | 0.9381 |
| uce_4l::mlp | pretrained | target_label_unexposed | 3 / 2 / 0 | 0.9377 |
| **scanvi_baseline** | **liver-specific, trained here** | trained on outer-training studies only | **5 / 0 / 0** | 0.9338 |
| scvi_liver_latent::mlp | liver-specific | as above | **5 / 0 / 0** | 0.9303 |
| hvg_pca_task_native::mlp | **task-native, no pretraining** | no corpus | **5 / 0 / 0** | 0.9265 |

⭐ **The best-scoring model has ZERO clean held-out studies.** Geneformer is `encoder_seen` on
GSE136103, GSE185477 and Liver_Atlas, and `unknown` on the other four — so **every fold is either
confounded or unresolved, and its 0.9438 cannot be compared to anything.** The verdict string in the
file is literally `no_clean_held_out_study_exists_for_this_encoder`. Same for TF-Sapiens, whose
exposure is `unknown` on all five folds.

### Clean-fold paired deltas against the liver-specific encoder

| model | clean studies | delta vs liver encoder | verdict |
|---|---:|---|---|
| geneformer::linear | **0** | — | **no clean held-out study exists** |
| tf_sapiens::linear | **0** | — | **no clean held-out study exists** |
| uce_4l::linear | 5 | **+0.0063 [+0.0004, +0.0095]** | above |
| scimilarity::linear | 5 | **+0.0127 [+0.0057, +0.0190]** | above |
| uce_33l::linear | 5 | +0.0038 [−0.0013, +0.0078] | crosses zero |
| hvg_pca_task_native::linear | 7 | **−0.0880 [−0.0933, −0.0807]** | below |
| uce_33l::mlp | 5 | **−0.0924 [−0.0961, −0.0892]** | **below** |
| uce_4l::mlp | 5 | **−0.0927 [−0.0964, −0.0899]** | **below** |
| scimilarity::mlp | 5 | +0.0087 [+0.0009, +0.0153] | above |
| hvg_pca_task_native::mlp | 7 | +0.00002 [−0.0050, +0.0060] | crosses zero |

⭐ **Where a clean comparison is possible, the margins are small and the head matters more than the
encoder.** UCE beats the liver encoder by **+0.006** with a linear head and **loses by −0.093** with
an MLP — a swing far larger than any gap between encoders. A liver-specific model trained here on
102 donors is within ~0.01 of billion-parameter foundation encoders.

⛔ **SCimilarity is the honest caveat**: it clears the liver encoder by +0.0127, and its checkpoint
state is **`encoder_seen`**. It is reported, not excluded, with that label attached.

### Why this panel is the paper's point

The contract commits to *"reporting every score beside an audit of what each model already saw during
training so that apparent skill is separated from prior exposure."* This run does exactly that, and
the audit **changes the ranking**: the two highest raw scores are the two that cannot be cleanly
compared at all. Exposure is resolved **per study**, not per model — Geneformer is `encoder_seen` on
three studies and `unknown` on four — so a single model-level label would have hidden it.

✅ Bootstrap: 10,000 resamples, donors as independent units (`cells_are_independent_units: false`),
paired on the same draws with shared multiplicities across every model, seed 20260824.
✅ `clinical_claim_allowed: false` recorded in the receipt. `calibration_fit: false`.

### 7F continued: against the task-native baseline, the field's advantage is mostly the HEAD

The comparison that matters most was missing from my first reading of the files. Same shared head,
same 50,000 rows, delta against a **50-dimensional PCA of raw counts (HVG+PCA)**:

| encoder | linear head | **two-layer MLP head** |
|---|---|---|
| geneformer_v2_316m | +0.0493 [+0.0404, +0.0599] | +0.0173 [+0.0119, +0.0230] |
| uce_4l | +0.0537 [+0.0453, +0.0638] | +0.0112 [+0.0047, +0.0165] |
| transcriptformer | +0.0458 [+0.0377, +0.0555] | +0.0124 [+0.0056, +0.0199] |
| scimilarity_v1_1 | +0.0482 [+0.0390, +0.0586] | +0.0117 [+0.0054, +0.0181] |
| uce_33l | +0.0448 [+0.0361, +0.0550] | **+0.0040 [−0.0029, +0.0095] crosses zero** |
| **scvi liver** | +0.0311 [+0.0251, +0.0388] | **+0.0039 [−0.0020, +0.0095] crosses zero** |

⭐ **The whole field's apparent advantage is largely a statement about the linear head.** Give the
head one hidden layer and a **50-dimensional PCA of raw counts matches the liver-specific encoder
outright**, and closes most of the gap to every pretrained encoder. Pretraining buys +0.045–0.054
over a linear probe and only +0.004–0.017 once the head can bend.

### Clean studies only, MLP head, against the liver encoder

⛔ **CORRECTED 2026-09-01.** A "GSE189600 dropped" variant of this table circulated with intervals
attached. **Those intervals exist in no producing file.** `clean_fold_paired_deltas_vs_liver_encoder.tsv`
computes its clean-studies delta as an unweighted mean over all clean studies **with GSE189600
included**. The dropped-study *point estimates* do reproduce (+0.0156, +0.0004, −0.0025, −0.0021),
but their intervals would require re-aggregating `bootstrap_distributions.npz` — reconstruction, not
reading. The panel plots the file's own values and names GSE189600 as the driver.

**The file's actual values, which are what panel 7F block c shows:**

| model | delta vs liver encoder | reading |
|---|---|---|
| **scimilarity_v1_1** | **+0.0087 [+0.0009, +0.0153]** | the **only** pretrained encoder clearing zero on a clean comparison |
| HVG+PCA task-native | **+0.0000 [−0.0050, +0.0060]** | **a dead tie with the liver encoder** |
| uce_33l | **−0.0924 [−0.0961, −0.0892]** | below |
| uce_4l | **−0.0927 [−0.0964, −0.0899]** | below |
| geneformer, transcriptformer | — | **no clean held-out study; no number reportable** |

⚠ Block c's metric is **study-balanced within the clean set**, because that is the only form its
producing file provides; blocks a and b are donor-class-balanced. A donor-class-balanced version with
GSE189600 dropped needs a re-aggregation run, not a reconstruction.

⭐ **The two encoders that look best pooled are exactly the two with no clean fold** — Geneformer
(+0.0135 vs liver) and TranscriptFormer (+0.0086). ⛔ And SCimilarity, the one pretrained encoder
that does clear zero cleanly, is itself `encoder_seen` in aggregate and marked
`ineligible_aggregate_development_encoder_seen` for sealed use.

✅ **Reproduction controls, not digest checks.** The loop reproduced **60 already-frozen shards
(Geneformer + TranscriptFormer) at max absolute difference 0.0**, temperatures and inner-fold losses
identical; the liver-latent row placement reproduced the frozen `scvi_baseline` at **2.2e-16**, and
reversing the reference rows (same multiset) broke it by **1.0**. The scorer matches two prior
independent scorers on 14 metrics to <1e-12. 4 guards demonstrably live, 3 verified positively.

### Caveats for the caption

⚠ **Report the donor-weighted metric, not study-balanced.** Study-balanced macro-F1 is dominated by
GSE189600 (149 cells, 2 donors), where UCE's MLP collapses to **0.53** and drags its study-balanced
score to 0.864 while its donor-weighted score is 0.938.
⚠ The liver encoder's 3-seed ensemble varies **encoder and head** seed while the frozen blocks vary
head only — a mild advantage to the liver encoder.
⚠ **Score does not track embedding width**: 128-d SCimilarity (0.9381) beats 1280-d UCE-33L (0.9305),
so width is a weak confound.

⛔ Three claims in my brief were wrong and are corrected here: 21247811 **did** compute metrics (the
`metrics_calculated: false` flags are in the *fit* job's receipts); TranscriptFormer was already
trained **and scored**, its matrix living at `model-training-500-21080428/embeddings/` where a
row-order-hash search misses it; and a **fifth** encoder exists on the identical rows, SCimilarity
v1.1 at 128-d.

---

## Reverse direction: 0 of 117 programs carry a coherent signed genetic signal — 2026-08-31

`executions/program_to_locus_twas-20260831T000000Z/`, `FINDINGS.md`, `SHA256SUMS` over 62 files.
Prespecification sealed; an addendum records two statistics added after sealing, with the original
hash preserved.

**Join rate 0.8979** — 3,341 of 3,721 program genes carry an OTTERS Z, matching the expected ~0.90
(losses: 290 ENSG absent from the TWAS panel, 90 symbols absent from GENCODE hg19). Programs asserted
at 117 / 3,721 / 7,093 from the five **named** lineage dirs; the 15-file glob was run on purpose and
rejected (20,605 rows / 4,788 genes / 406 programs).

✅ **The null holds AND rejects** — the hard part of this lane. On a known-empty Z field
P(p<0.05) = **0.0519 / 0.0525 / 0.0514**, KS uniformity p 0.21–0.91. With an effect planted in one
program it rejects: power **0 → 0.38 → 0.98** across δ = 0, 1.0, 2.0.

**Programs passing the pre-specified test: 0 of 117.** Signed statistic, BH q<0.05 within 117
programs, under both nulls, on either primary trait.

⭐ **Two independent genetic tests now agree.** S-LDSC gives **0 of 585** program×trait cells on tau;
TWAS aggregation is a **different estimand on a different substrate** and also returns nothing.
⛔ Neither implies the other, and the write-up must say so.

### ⛔ Three methodological findings that change how these results read

1. **The genome-rotation null (N1) is ANTI-CONSERVATIVE for liver gene sets.** On synthetic pool sets
   Q/N1 rejects at **0.091** against a nominal 0.05, while Q/N3 gives 0.0485 — because program genes
   carry more TWAS signal than the genome average (pooled mean χ² **4.225 vs 3.668**, p = 8.9e-3).
   **Anything read off N1 alone is inflated.**
2. **The max-|Z| statistic cannot reach significance in this design at all.** A permutation null
   reuses the same Z values, so its p-value floor is ~n_genes/n_universe, which exceeds the BH
   threshold for **116 of 117** programs. ⛔ **Its zeros are `not_run`, not negative.**
3. ⛔ **Correction to a number I put in the brief.** "Median 14 LD blocks, 24 of 117 untestable" comes
   from the **forward** variant→gene lane and counts components of the credible-set *variant* LD
   graph. This lane's **gene-side** blocks are median **32** on ukbb_alt. On pdff, 75 of 116 programs
   span fewer than 10 blocks, which is why pdff cannot answer the question whatever it returns.

### One post-hoc hit that survives everything, reported as what it is

**cholangiocytes::5 / ukbb_alt** on the count statistic: q = 0.0023 (N1) and 0.0070 (N3); at
B = 5×10⁶, **K = 9 against a null mean of 1.73, p = 3.1e-5**; its hot genes span **8 distinct 5 Mb
blocks** so it is not one GWAS locus (block-collapsed p = 1.0e-4); and it replicates in the
independent P0.001 model (p = 2.2e-5).

⛔ **But its signed statistic is only nominal (p = 0.035) and the nine genes split 6+/3−.** The honest
reading is *an unusual concentration of individually strong ALT associations across eight loci*
(STARD10, RBP4, GPX1, NUPR1, DDT/CHCHD10, SDC1, ZNHIT1, SCAND1) — **not a coordinated program-level
effect.** The other two candidate cells died: `tcells::3` is MHC-driven (block-collapsed p = 0.20),
and `cholangiocytes::16` sits on ukbb_ggt, which has no heritability diagnostic at all.

### Two self-inflicted traps the producer caught

⛔ A first power simulation **planted the effect in every program at once**, which makes the N3 null
invariant by construction and made N3 look powerless; redone one program at a time it tracks N1.
⛔ At B = 2,000 the smallest attainable p (5.0e-4) sits **above** the BH threshold (4.27e-4), so the
simulation *could not have rejected a perfect signal*. There is now a guard that refuses to run
unless the floor is 4x below the threshold, and it fires on that exact configuration.

✅ **Scope limit honoured**: `program_credible_set_summary.tsv` gives a credible set per program —
median **1,470** fine-mapped variants, range 189–19,056 — **never a named causal variant.**
✅ **11 of 11 guards demonstrably live**, each violated on purpose with the tamper asserted to change
the bytes.

---

## Zonation deep null: p = 5.0e-05, and it is a measurement now — 2026-08-31

`executions/gse105127-zonation-deepnull-20260830T220000Z/`. Job 21278341, **7:52:46**, all
**20,000 / 20,000** permutations completed, `stopped_on_time_budget: false`.

| | 200 permutations | **20,000 permutations** |
|---|---|---|
| accuracy p (one-sided) | 0.004975 | **4.99975e-05** |
| MAE p (one-sided) | 0.004975 | **4.99975e-05** |
| null mean / q95 | 0.329 / 0.422 | **0.33266 / 0.45614** |

⭐ **The earlier 0.004975 was exactly 1/201 — the floor at 200 permutations, not a measurement.**
The new value is 1/20,001, so it is *still* at its floor: **not one of 20,000 within-participant
permutations reached the observed accuracy of 0.7895.** The difference is that the floor is now
100x lower, and the reported number is bounded by the budget rather than by anything about the data.
Report it as **p < 1e-4**, not as a point value.

The null is stable and matches the shallow run (mean 0.329 → 0.33266, q95 0.422 → 0.45614), so the
deeper run confirms rather than revises. Verdict unchanged: **positive**.

⚠ Still one cohort, n=19 participants. The `indeterminate`-on-a-negative label was pre-registered and
never applied, because the result was positive — see the Figure 7E field-selection bug, where that
conditional was mistakenly plotted as the outcome.

### The three engineering faults this job cost, all mine

Attempt 1 (21251830) ran the full 20 h walltime and wrote **zero bytes**:
1. ⛔ **I dropped `OMP_NUM_THREADS`/`MKL_NUM_THREADS`/`OPENBLAS_NUM_THREADS=8`** when writing a fresh
   sbatch from the validated one, so BLAS oversubscribed the allocation. With them restored the cost
   is **1.41 s/perm** — my original 1.2 s estimate was nearly right; the environment was the problem.
2. ⛔ stdout was block-buffered, so the log was empty for 20 h.
3. ⛔ Everything serialised at the end, so the timeout lost all of it.

✅ The rebuilt job carries a **wall-clock budget inside the loop** that stops cleanly and reports the
count *achieved*, plus partial nulls written every 25 permutations. It needed neither: the estimate
was right once the environment was. ⭐ **The fix was not a better estimate — it was making the job
degrade instead of vanish when the estimate is wrong.**


---

## ⛔ CORRECTION 2026-09-01: the RNA→chromatin disease-dependence is SINGLE-COHORT

I have repeatedly reported disease-dependence as a two-cohort finding. **It is not.** From
`executions/crossmodal-severity-stratification-20260831T190048Z/results/ALL_ARMS_SUMMARY.json`:

| arm | Δ | p | verdict |
|---|---|---|---|
| GSE267145 activity≥3 (pre-registered) | +0.2087 | 0.0000 | DISEASE_DEPENDENT |
| **GSE296875 fibrosis_any** | +0.0548 | **0.18** | **UNDERPOWERED_untestable** |
| **GSE296875 steatosis≥5%** | +0.0161 | **0.70** | **UNDERPOWERED_untestable** |
| GSE267145 NOR-dropped | +0.0595 | 0.24 | UNDERPOWERED_untestable |

**Both pre-committed GSE296875 arms are `UNDERPOWERED_untestable`.** The disease-dependence rests on
GSE267145 alone, and within it on 24 normal livers. ✅ What *does* recur across two cohorts and two
assays is the **base result** (+0.13345 / +0.12985) — that part stands. Do not carry the
disease-dependence into any transfer pitch as a replicated finding.

## ⛔ CORRECTION: "opaque keys" was a policy label, not a data fact

GSE267145 H3K27ac regions are literal `chr1:96466-96736` strings on **GRCh38.p14 primary contigs**,
widths 165–4,387 bp — the column is even named `opaque_source_feature_key`. What is unresolved is the
0-based/1-based convention, and it is **unresolvable from sources**: the crosswalk records
`one_based_inclusive_bounds_compatible: true` AND `zero_based_half_open_bounds_compatible: true`, and
no primary source states it. `interval_overlap` and `coordinate_conversion` are `blocked_actions`, so
the lane is **fail-closed by governance, not blocked by genomics**.

⭐ The ambiguity is **±1 bp against peaks 165–4,387 bp wide** — any overlap join is stable to <1%
either way. The honest amendment is to **compute under both conventions and report both**; if they
agree, the convention never mattered.
⚠ GSE296875's coordinates are **also unaudited** — `start_1based` is an assertion in a parse regex,
and its own TOML still lists the peak-coordinate crosswalk as an open blocker. Both sides are
unresolved; only one is fail-closed about it.

## ✅ [RESOLVED LATER] PROVENANCE: `0.7918` had no producing file — now it does

The GSE267145 fibrosis tie ceiling **0.7918** appears **only** in this file's prose. It recomputes
correctly to **0.791771** from the frozen marginal (71/15/9/4) in `endpoint_summary.json` — the same
routine reproduces GSE130970's stored 0.95584 to 6 dp — but it must be given a producing artifact
before it is used anywhere. The stored source says only "roughly 0.79" in prose.

Similarly, the retired-QWK claim (`0.638 → 0.225 under corrected units`) appears only in a sealed
prespec with **no producing artifact for either number**; every other `0.638` in the tree is a
Spearman. Verify before citing.

---

## ⚠ [SUPERSEDED — σ_e IS A RANGE] SHIPPED (A): the Kleiner measurement-error constant — 2026-09-01

⛔ **THIS SECTION'S CONSTANT IS SUPERSEDED. σ_e = 0.6433 is ONE deposit, not "the" constant.**
PXD051911 gives **0.5092 [0.304, 0.694]** (27 at-surgery pairs) and **0.4655 [0.365, 0.548]** (30
follow-up pairs, formally **discrepant**, ratio CI [0.537, 0.956] excluding 1). **Ship the range
0.4655–0.6433.** ⛔ **0.6433 is the ANTI-CONSERVATIVE end** for ceiling-proximity and disattenuation
claims — the two things this project most often asserts. Every number below is correct *for
GSE193066* and must not be quoted as general. See the full retraction further down.

`executions/kleiner-reliability-20260901T122243Z/`. Prespecification `d5b39ef2…`, **verified at
runtime**. 58 paired repeat biopsies, 106 participants, median ~2 y apart. Source series matrix
sha256 `f35dae7e…`.

### ⛔ The deliverable is σ_e, NOT an ICC — and that correction is the contribution

| quantity | value |
|---|---|
| **single-read error, σ_e** | **0.6433 stages**, CI95 **[0.517, 0.760]** |
| error variance | 0.4138 stages² |
| exact agreement | 28/58 = 0.483 |
| within ±1 | 52/58 |
| ICC(1,1) *in this sample* | 0.1886 |

⭐ **Why σ_e and not ρ.** ICC and kappa depend on the fibrosis-stage variance of the sample they are
measured in. In this one deposit **the same σ_e yields ρ = 0.19, 0.28 or 0.47** depending only on
which stage variance is used as the denominator (paired-set one-way model, paired-set first biopsy,
full-106 first biopsy; classical range correction gives 0.400). **Publishing a single ρ would be
publishing a property of our sample.** σ_e is on the stage scale and transports.

**How anyone uses it:**
```
rho      = 1 - 0.414 / Var(Kleiner stages in your sample)
r_latent = r_observed / sqrt(rho)          ceiling: r_observed <= sqrt(rho)
beta_latent = beta_observed / rho
```
ρ by stage SD: 0.7→0.159, 0.9→0.491, 1.0→0.588, 1.1→0.659, 1.2→0.714, 1.3→0.756, 1.5→0.817.
⭐ A worked example that bites: a reported r = 0.5 against Kleiner in a population like this one
implies r_latent = 1.15 — **impossible**, so either the population has more stage variance or the
result needs re-examining.

### ⛔ Three honesty constraints the producer surfaced, and they bind

1. **The ICC is not separated from zero.** Parametric CI95 **[−0.070, 0.424]**, p = 0.075;
   permutation p = 0.103, observed at the 89.7th percentile; bootstrap **9.7% of draws at or below
   zero**. ⛔ **0.189 must never be quoted as an established value** — I was about to do exactly that.
2. **The symmetry test supports noise but cannot exclude real progression.** Sign test **p = 1.0**
   (15 up / 15 down), Wilcoxon p = 0.742, mean Δ CI [−0.276, +0.207]. But the mean-delta arm has 80%
   power only for **0.34 stages over 2.49 y = 0.135 stages/year**, which **cannot exclude published
   NASH progression rates.** ⭐ So the argument rests on the **sign and shape** of the distribution,
   not on the precision of the mean, and must be written that way.
3. **Interval-dependence tests are an underpowered null**, labelled as such per the prespec, not read
   as support (Spearman MDE at 80% power is 0.36; observed |ρ| ≤ 0.14).

✅ **The pairing is validated by a positive control that could have failed**: age ICC **0.978**,
age r **0.993**, sex concordant **58/58**. Wrong pairings would not track age and sex.

### Scope limits recorded on the artifact

⛔ **No inter-reader substrate exists in this project** — no deposit records who scored a slide, how
many did, or whether anything was re-read. This is a **test–retest** statistic across two physically
distinct biopsies: it folds reader variability, **biopsy sampling variability**, and any real change
into one term and cannot separate them. We can never claim Kleiner "agrees with itself better than
two pathologists agree."
⛔ **The paired set is selected**: all 58 are high-risk at first biopsy (the 48 unpaired split 35 low
/ 13 high), and their stage variance is **0.74x** that of the full 106.
⛔ **Stage 4 appears in 1 of 164 deposited samples** — the top of the scale is truncated.
⛔ NAS is descriptive only: its delta is **asymmetric** (34 down / 7 up), real improvement plus
regression to the mean, so it cannot carry the reliability argument.

**Guards: 12 demonstrably live, 1 verified positively, 1 live only under a changed constant.**

⭐ **What this replaces.** `docs/RESULTS.md` (~line 745) currently assumes a Kleiner reliability band
of 0.6–0.9 and says so explicitly — an assumption that collapses a 2.79x claim to 0.87–1.31. The
measured range-corrected band is **[0.400, 0.470]**, *below* the assumed floor, with σ_e as the
transportable form. **That edit is pending and must be applied to `docs/RESULTS.md`.**

---

## ⭐ OOD TRANSFER WORKS: a model trained elsewhere recovers 85–99% of native performance — 2026-09-01

`executions/full-axis-rebuild-and-ood-transfer-20260901T122932Z/`. Job 21323999 (2:18).
Target **GSE268273, 109 participants, real donor key**, tie-structure ceiling **0.952322**.

✅ **Firewall honoured:** `target_predictions_written_and_hashed_before_target_outcomes_opened = true`,
predictions sha256 `cb263414…`. Predictions were frozen before any outcome was read.

| source head | transfer ρ | native in-cohort ρ | paired delta | excludes zero? | **% of in-cohort signal recovered** |
|---|---|---|---|---|---|
| linear SVR | 0.5219 | 0.5637 | **−0.0073** | **No** | **98.7%** |
| cumulative ordinal logistic | 0.5219 | 0.6021 | −0.0802 | **No** | **86.7%** |
| elastic net | 0.5142 | 0.6064 | −0.0922 | **No** | **84.8%** |

**Primary transfer Spearman 0.5743** against a ceiling of 0.9523.

⭐ **This is the result a released tool needs.** A model fitted on one cohort, applied to 109
participants it has never seen, **loses nothing that excludes zero** against a model trained natively
in the target — and it does so **across a full pipeline change**: source cohorts are
featureCounts 2.1.1 / reverse-stranded / 10x primary assembly; GSE268273 is RSEM 1.3.1 / unstranded /
GRCh38.p14-with-patches. Only GENCODE v49 is shared.

⚠ Layout matters and is reported: paired-end n=36 gives ρ 0.620 (ceiling 0.9636) against the
single-end arm — report the stratum, not just the pooled number.
⛔ **CORRECTED 2026-09-02 — the AUROC 0.901 / AUPRC 0.644 previously printed here was NOT this
transfer's number.** It is `results.native__cumulative_ordinal_logistic.advanced_f3_f4` in
`job2_results.json` — the **native in-cohort baseline**, and specifically the best of four native
heads (ridge 0.8836, elastic net 0.8626, linear SVR 0.8322). Printed directly under "Primary transfer
Spearman 0.5743", it read as the transfer AUROC. **The shipped instrument's own value on this cohort
and contrast is AUROC 0.9187** (DeLong CI95 [0.8571, 0.9803], AUPRC 0.7813), measured through
`release/masld-severity-v1/score.py` in
`executions/standardized-metric-panel-20260902T130923Z/results/panel.json`.
✅ The null statement was and remains correct: **prevalence (0.174) is NOT the AUPRC null; the
random-scorer mean is 0.206** (p95 0.293). The metric panel measures the same quantity by the same
construction and gets 0.2069.

### Four guards, all demonstrably live, each with its own failure text

| guard | on real input | on a constructed violation |
|---|---|---|
| tie-structure ceiling | 0.952322, matches the frozen receipt to 6 dp | 0.950725 if one F4 were relabelled F0 |
| paired bootstrap indices shared | sha256 `2e93c436…` | a different index set per arm — would make the interval marginal, not paired |
| not scoring its own training rows | **0** GSE267145 and **0** GSE135251 ids in target | 109 when the target id set is compared with itself |
| target row alignment | ρ **0.5743** | **0.1944** on row-shuffled outcomes |

### The input rebuild that made it possible

| | delivered axis | full quantification |
|---|---|---|
| GSE268273 coverage of the 42,163-gene model axis | 13,710 = **32.5%** | **42,163 = 100.0%** |
| genes gained by rebuilding | — | **+28,453** |

✅ **GSE268273 confirmed non-voom and bitwise reproducible** from 109 retained `rsem.genes.results`:
`bitwise_identical: true`, max abs difference **0.0**, **0 cells differing**. The earlier "32.5%
coverage" and "global limma-voom" concerns both described the *delivered* axis, not the data.

⛔ **A rebuild trap found and avoided.** The delivered 14,078-gene matrix does **not** equal a naive
stable-ID aggregation of the full RSEM output — **3 genes, 327 cells differ**. The cause is
pseudoautosomal genes (P2RY8, CSF2RA, ASMTL-AS1) whose X and PAR_Y rows are **elementwise identical**:
a naive stable-ID sum double-counts them (`delivered_over_stable_rule_ratio = 0.5`). The delivered
matrix is correct; a naive rebuild would have silently doubled three genes.

⚠ An amendment was written mid-run and discloses its own timing: job 21323950 was cancelled at 11m47s
**with no results file and no target outcome ever opened**, so the change (source-side PCA scaling and
a narrowed SVR grid, both purely computational) could not have been influenced by a result.

---

## ⚠ [HEADLINE SUPERSEDED] SHIPPED (B): a severity score, reproducibility claim NARROWED — 2026-09-01

⛔ **This section's original title was "reproduces ~4x better than the histology label". That is
retracted as a general claim**: the Kleiner label's own ICC is **0.62** in PXD051911, not 0.19, so the
comparison is a **within-sample contrast at a stated stage variance (0.5091 stages²)** and nothing
more.

⚠ **This marker does NOT certify the rest of the section.** A second claim below is also superseded:
**"Predictive intervals are calibrated: nominal 50% → empirical 0.478"** is **RETRACTED** — measured
coverage is **0.404** on the shipped quantile arm and **0.330** on the earlier hybrid arm. Any σ_e
quoted below as a single constant is likewise superseded by the range 0.4655–0.6433.
**Everything else below is as measured.**

`executions/latent-severity-20260901T124751Z/`. Trained on the **521-sample** graded pool; reproducibility
validated on the **58 GSE193066 participant pairs**, which the model never saw.

| arm | test–retest ICC | **ΔICC vs histology** | ΔICC CI95 | validity ρ | gate |
|---|---:|---:|---|---:|:--|
| **histology (the label itself)** | **0.1886** | — | — | — | — |
| b1_stararm | **0.7745** | **+0.586** | [0.307, 0.848] | 0.510 | ✅ |
| b1_pooled_std | 0.7293 | +0.541 | [0.258, 0.813] | 0.508 | ✅ |
| **b1_latent** | **0.7167** | **+0.528** | **[0.245, 0.809]** | 0.509 | ✅ |
| ridge_reference | 0.6909 | +0.502 | [0.211, 0.788] | 0.530 | ✅ |
| b1_kallisto | 0.6458 | +0.457 | [0.150, 0.756] | 0.473 | ✅ |

**The model score reproduces across repeat biopsies roughly four times better than the Kleiner stage
does, and the paired ΔICC excludes zero** (99.98% of draws > 0 for the reference arm). Predictive
intervals are **calibrated**: nominal 50% → empirical 0.478 [0.434, 0.522]; nominal 80% → 0.798
[0.761, 0.832].

### ⭐ The control that makes the claim non-trivial: `adv_sex_only`

| adversarial arm | ICC | ΔICC | validity ρ | gate |
|---|---:|---:|---:|:--|
| **adv_sex_only** | **1.0000** | **+0.811** | **−0.229** | ❌ **FAILS** |
| adv_constant_plus_noise | 0.3126 | +0.124 | −0.053 | ❌ |
| adv_random_signature | 0.1446 | −0.044 | 0.004 | ❌ |
| adv_library_size | 0.0729 | −0.116 | 0.086 | ❌ |

⭐ **Sex is perfectly reproducible across repeat biopsies — ICC 1.0, the best of anything tested —
because it does not change.** If reproducibility alone were the claim, a sex-only predictor would win
outright. It fails on validity (ρ = −0.229), and the gate is **directional by design**: the file notes
that *"a two-sided 'excludes zero' test passes a sex-only predictor whose association with stage is
significantly NEGATIVE, which is why it is not used."* **Reproducibility is necessary and nowhere near
sufficient; the pairing of a reproducibility gain with a directional validity gate is the whole design.**

### ⭐ It is not merely "continuous beats ordinal"

Discretising the model score to the label's own coarseness **keeps the advantage**:

| discretised arm | ICC | ΔICC | CI95 |
|---|---:|---:|---|
| b1_latent, equal-frequency bins | 0.6122 | **+0.424** | [0.120, 0.732] |
| ridge_reference, marginal-matched bins | 0.5619 | +0.373 | [0.057, 0.691] |
| b1_latent, marginal-matched bins | 0.5281 | +0.340 | [0.019, 0.654] |

Both still exclude zero. The gain is not an artifact of comparing a continuous score against a 5-level
ordinal.

### ⛔ NEGATIVE: composition-invariance (B2) does NOT help — my hypothesis was wrong

| B2 arm | ICC | ΔICC | CI95 | verdict |
|---|---:|---:|---|---|
| b2_clrsub | 0.5000 | +0.311 | [0.005, 0.592] | barely excludes zero |
| b2_stararm | 0.4396 | +0.251 | [−0.056, 0.532] | **crosses zero** |
| b2_clr16 | 0.3426 | +0.154 | [−0.176, 0.449] | **crosses zero** |
| b2_kallisto | 0.3240 | +0.135 | [−0.214, 0.457] | **crosses zero** |

**Every composition-invariant arm is worse than its plain counterpart, and three of four cross zero.**
I predicted that removing composition would improve test–retest reproducibility by stripping biopsy
sampling variability. The opposite happened: **composition carries reproducible severity signal, and
removing it costs more than the sampling noise it removes.** Report as a clean negative on the
invariance hypothesis, not as a tuning failure.

### Scope limits recorded on the artifact

⛔ **Cross-sectional only** — "the intervention arm reverses sign 5/5, which forbids any monitoring,
longitudinal or treatment-response claim. The 58-pair analysis is a REPRODUCIBILITY test, not a
progression measurement."
⛔ **Units**: 521 **SAMPLES** for the pool (`donor_key_exists=false`; participant/donor forbidden, 844
forbidden); **58 PARTICIPANTS** for validation.
⛔ **No cell-fraction claim** — composition enters only as a nuisance basis.
⛔ **No inter-reader claim** — σ_e is test–retest across two physically distinct biopsies and folds
reader variability, biopsy sampling variability and real change together.

Guards: **4 demonstrably live, 4 passed, 0 failed.**

---

## ⚠ [HEADLINE SUPERSEDED] The instrument transfers, and a small panel keeps most of it — 2026-09-01

⛔ **Original title said "transfers ABOVE native". Retracted**: the paired delta against native is
**+0.0509 CI [−0.028, +0.137], which includes zero** (MDE 0.118). The defensible claim is **not worse
than native**. The native baseline is itself a 10-fold in-cohort ridge estimate with its own sampling
error, and must be named whenever the comparison is quoted.

`executions/panel-and-ood-20260901T135507Z/`. Prespec `d671007f…` **verified at runtime**. Firewall:
predictions sealed at 14:04:33.427 UTC, target outcomes opened 14:04:33.431 — hashed before opened,
predictions sha256 `72b85ecc…`.

### The latent model out-of-cohort: ρ 0.649 against a native baseline of 0.609

| | value |
|---|---|
| transfer Spearman on GSE268273 (109 participants) | **0.6492** |
| native in-cohort 10-fold participant CV ridge | 0.6094 |
| **fraction of in-cohort signal recovered** | **106.5%** |
| paired delta | +0.0400, **does not exclude zero** |
| ρ / tie-structure ceiling (0.952322) | 0.682 |

⭐ **A model trained on 521 samples from five other cohorts, applied to 109 participants it has never
seen, performs at least as well as a model trained natively inside the target** — across a full
quantification change (featureCounts / reverse-stranded / primary → RSEM / unstranded /
GRCh38.p14-with-patches). That is the property a released tool needs, and it is stronger than the
generic-head result (85–99%) that preceded it.

✅ It **reused** the prior execution's axis, PAR-duplication fix, bootstrap indices and native
baseline rather than rebuilding them, and documented an int64-vs-int32 cast on the bootstrap receipt
as lossless.

### ⭐ A FIVE-gene panel retains ~83% of the reproducibility gain

**`smallest_k_whose_dICC_excludes_zero` = 5** — the smallest size tested — and
`largest_k_whose_dICC_does_not_exclude_zero` = **None**. Every panel size tested clears zero.

| k=5 panel | value |
|---|---|
| test–retest ICC on the 58 pairs | **0.6272** |
| **ΔICC vs histology (0.1886)** | **+0.4386** [+0.121, +0.719], 99.63% of draws > 0 |
| validity ρ | 0.5072 [0.350, 0.640], directional gate **passes** |
| full-transcriptome model for comparison | ΔICC +0.528 |

**Genes: FAP, NALCN, STMN2, SOX9, ITGBL1.** ⭐ This is interpretable liver-fibrosis biology, not an
opaque signature — **FAP** is the canonical activated-fibroblast marker and a fibrosis drug target,
**SOX9** marks ductular reaction, **ITGBL1** is ECM/fibroblast. At k=10 it adds **THBS2** and
**FBLN5**, both established fibrosis biomarkers, plus CXCL6 and MOXD1.

⭐ **Five genes run on qPCR or a targeted panel in any pathology lab.** That is the difference between
a result and something a trial can use.

### The controls that make the panel credible

- ⛔ **Selection was done inside the folds**, not on the full data — `fold_vs_fulldata_gene_overlap`
  shows 3–5 genes shared per cohort fold rather than an identical set, so the classic
  select-then-report-CV leak is excluded.
- ✅ **`n_sex_linked_genes_in_panel = 0`**, and score-vs-sex ρ = **−0.166**. A small panel is
  especially vulnerable to accidentally encoding a stable nuisance — this one does not.
- ✅ **Random-panel null**: 10 size-matched random panels give ΔICC mean −0.017, p95 0.182, and the
  real panel **beats all of them**.
- ✅ The directional validity gate is retained, with its note that a two-sided test would pass
  `adv_sex_only`.

Guards on both runs: all demonstrably live, none failed.

---

## ⭐⭐⭐ THE HEADLINE: the score predicts which biopsy is wrong — 2026-09-01

`executions/latent-vs-label-20260901T153025Z/`. 58 participants, 30 model–histology discordant.

**At biopsy 1 the model says X and histology says Y. Does biopsy 2 move toward X?** Yes, over and
above regression to the mean.

| | RTM-only | **RTM + model score** |
|---|---|---|
| r² for Δ fibrosis | 0.4836 | **0.6452** |
| | | **Δr² = +0.1616** |

| statistic | value |
|---|---|
| β (model score, per SD) | **+0.4174** |
| CI95 | **[+0.2159, +0.5992]**, excludes zero |
| permutation p (10,000) | **1.0e-04** — the floor |
| MDE at 80% power | 0.2743 (observed 0.417 clears it) |
| content-free random-gene score | **+0.0881** (ICC 0.145, validity +0.004) |
| jackknife over 58 participants | **all positive**, max shift 0.042 |

⭐ **When the model and the biopsy disagree, the NEXT biopsy tends to move toward the model.** This is
the strongest available evidence that the reproducibility gain is not mere stability — the score is
tracking the latent quantity, and the disagreement is the biopsy's error. It is verified
**prospectively by an independent measurement**, and it is clinically actionable: it says which
biopsies to distrust.

⚠ **The registered statistic is the Y1-adjusted one and the honesty is preserved**: direction-only
AUC **0.729** [0.515, 0.906] adjusted, but **0.449 unadjusted** and raw sign concordance **0.500** —
both reported and both labelled *not the registered statistic*. The claim rests on the regression,
not the AUC.
⛔ **This is not a monitoring claim.** It is about which read is wrong, not about progression, and
magnitude is uninterpretable because the deposit records no interval.

### ⛔ Test 1 is confounded — and the agent found it mid-analysis and disclosed it

**All 58 paired participants are high-risk by the PLS signature, which is derived from BIOPSY 1
EXPRESSION.** Selection is therefore on a function of the biopsy-1 transcriptome, which restricts the
biopsy-1 expression range and attenuates any correlation involving Y1 — and our model *is* an
expression model. Var(Y1) 0.577 vs Var(Y2) 0.449.

✅ **Only the Ybar-vs-Y2 comparison is clean**; Ybar-vs-Y1 and Ybar-vs-mean-single inherit the
confound. Recorded as **discovered during analysis, not pre-registered**.
✅ Averaging does improve the reference as predicted: ρ single read **0.1886** → two-read average
**0.3174**.
⚠ The paired set spans stages 1–3 only (Y1 marginal 0/12/23/23/0; Y2 0/8/34/15/1) — no stage 0, one
stage 4.

### Test 3 — null, as declared in advance

`"expected_result": "NULL - declared before computing"`. Uncertainty vs |Δ fibrosis| ρ = **−0.203**
[−0.450, +0.064], does not exclude zero, MDE 0.371. Change across the pair is largely noise, so this
was not achievable in principle. Reported, not hunted.

### Reliability weighting — a no-op, and correctly predicted analytically before running

`"prediction_made_before_running": "a no-op: per-cohort rho differs only through stage variance, so
optimal weights (1/sigma2_e,c) are equal"`. ⭐ Since σ_e is constant across cohorts, per-cohort ρ
varies **only** through the stage-variance denominator — so inverse-variance weights are uniform and
weighting cannot help. **Weighting by ρ would have been wrong.** Per-cohort ρ 0.679–0.751.

Guards: 4 demonstrably live, 1 verified positively, 0 failed.

---

## ⭐ Lane D — no etiology-specific loss detected, VERIFIED and RESCOPED — 2026-09-01

`executions/etiology-transfer-fibrosis-distinctness-20260901T135543Z/`. GSE276114, 177 participants,
one source-native 3-level fibrosis outcome, ONE quantification, three etiologies. Etiology is the
only thing that changes across cells, which is the whole reason this cohort can ask the question.
Prespecification sealed before any transfer cell was computed; stage-1 substrate pinned by sha256
(`c96cfd1a...`) and verified at runtime.

**Substrate:** MASLD 81 / CVH 82 / ARLD 14; outcome levels 39 / 24 / 114; 19,681 genes after symbol
collapse and all-zero drop.

**Decisive stage-matched arm** (marginals forced identical across etiology at 13/11/42, n=66 per
etiology, matched tie ceiling **0.8545**):

| cell | transfer ρ | native ρ | paired Δ | interval |
|---|---:|---:|---:|---|
| CVH→MASLD ridge | +0.690 | +0.690 | −0.0013 | [−0.0779, +0.0912] |
| CVH→MASLD elastic net | +0.679 | +0.676 | +0.0037 | [−0.0792, +0.0817] |
| CVH→MASLD linear SVR | +0.676 | +0.638 | +0.0368 | [−0.0735, +0.2143] |
| CVH→MASLD ordinal | +0.680 | +0.702 | −0.0206 | [−0.0958, +0.0618] |
| MASLD→CVH ridge | +0.686 | +0.653 | +0.0318 | [−0.0358, +0.1172] |
| MASLD→CVH elastic net | +0.647 | +0.648 | −0.0081 | [−0.1275, +0.1225] |
| MASLD→CVH linear SVR | +0.682 | +0.563 | +0.1154 | [−0.0168, +0.2448] |
| MASLD→CVH ordinal | +0.634 | +0.671 | −0.0394 | [−0.1373, +0.0490] |

**8 of 8 straddle zero, in both directions, across four heads.** Per the sealed prespec this reads as
*transfers*: the fibrosis signal is etiology-general and the instrument reaches further than claimed.

⛔ **ARLD enters no cell, and this is degeneracy rather than low power.** 14 of 14 sit at F4, so the
outcome has one level and Spearman is undefined at every effect size. No n would make the cell askable.

### ⚠ The load-bearing weakness, flagged before the verifier reports

The argument needs the pipeline to be *capable* of showing a loss, and the batch control is what
supplies that: an etiology-mixed cross-run-block split loses skill (early→late Δ −0.154 to −0.236,
CI excludes zero, both etiologies on both sides). **But that arm trains on n=34 and tests on n=129,**
while its reverse trains on 129, tests on 34, and has a native Spearman of only 0.113. A training-size
difference and not a block effect could produce the whole drop. If the batch control does not
establish detectability, the etiology null is just a null, and the honest label may be
`UNDERPOWERED_untestable` rather than *transfers*.

⚠ Second concern: the quoted intervals are over **200 subsampling draws**, and the results file itself
records `interval_is_over_subsampling_draws_not_a_participant_bootstrap: true`. That may not be the
right uncertainty for a claim about participants.

An independent verifier that did not produce this result is checking both, plus participant leakage
between source and target, the ARLD exclusion, and the four guards. **Do not propagate this to
`docs/` or any figure until that returns.**

---

## ⚠ Lane F — aspect decomposition at n=377 stays indeterminate — 2026-09-01

`executions/aspect-axis-four-arm-21324229/`. Prespec digest guard `G_PRESPEC_DIGEST` fired on a
deliberate tamper (rc=1), so it is demonstrably live.

| arm | outcome | axes | per-aspect status |
|---|---|---:|---|
| A | ONE_AXIS | 1 | steatosis supported; ballooning and fibrosis indeterminate; lobular not independent |
| B | INDETERMINATE | 0 | all four untestable |
| C | TWO_AXES | 2 | fibrosis and steatosis supported; ballooning and lobular indeterminate |
| D | ONE_AXIS | 1 | fibrosis supported; steatosis indeterminate; ballooning and lobular untestable |

**The four arms disagree, and the disagreement is the finding.** These are four different instruments
with four different label vocabularies, so no arm is the answer and a naive pool was always
forbidden. Best case is arm C: fibrosis and steatosis separate; ballooning and lobular inflammation
never do in any arm. The earlier n=99 `ONE_AXIS_ACTIVITY` and n=180 `TWO_AXIS_ACTIVITY_FIBROSIS`
readings were power artifacts, but n=377 does not resolve to three independent axes either.
`docs/ROADMAP.md:41` keeps its `indeterminate` state; this carries no claim.

---

## ⛔ GSE249997 is NOT a validation target — correction to my own recommendation — 2026-09-01

I proposed it as "the cleanest external validation target for the severity instrument" on the
strength of `exposure_status = "target_label_unexposed"`. **Both halves of that were wrong.**

1. ⛔ **It is not on disk.** A full-tree scan finds exactly one file:
   `config/datasets/gse249997.toml`. `automatic_download = false`, and nothing matching PRJNA1051582
   exists anywhere.
2. ⛔ **It has no severity endpoint.** Its own registry blocker says so verbatim: *"until then this
   cohort has no stage endpoint."* Early-MASLD only, no healthy controls, no stage range, and 77 GEO
   samples against 70 participants in the publication, so the samples are not known to be independent
   people.

⭐ **The lesson: `target_label_unexposed` is necessary but nowhere near sufficient.** It says a model
never saw the labels. A validation target also needs labels to exist. I conflated *unexposed* with
*usable* and would have funded a lane that could not have produced a number.

The locally processed bulk cohorts are exactly ten: GSE126848, GSE130970, GSE135251, GSE162694,
GSE167523, GSE174478, GSE193066, GSE213621, GSE240729, PRJNA512027 (retired, library-prep confound).

---

## ⭐ PXD051911's plasma arm is the non-invasive substrate — the crosswalk problem is absent, not solved — 2026-09-01

Measured directly from `data/PXD051911/meta_data.txt` (222 rows, 154 unique `patient_name`):

| arm | rows | participants |
|---|---:|---:|
| liver proteome | 58 | 58 |
| **plasma proteome (baseline)** | **143** | **143, already 1:1** |
| plasma follow-up | 41 | 41 |
| liver ∩ plasma, person-level | — | **48** |

- ⭐ **Complete histology on all 143**: Kleiner F0 39 / F1 69 / F2 25 / F3 7 / F4 3, NAS 143/143,
  SAF No_MASLD 39 / MASL 74 / MASH 30. That is a **wider stage range than the liver arm**, which has
  no F4 at all (F0 13 / F1 31 / F2 11 / F3 3).
- ⭐ **All 41 follow-up participants carry their own second Kleiner stage and their own second plasma
  proteome**, and all 41 are also in `initial_sample`. A second paired-visit cohort.
- ⭐ **No crosswalk is required.** The Zeybel route stalled at `AMBIGUOUS` because it had to join two
  separate deposits. Here plasma and liver are the same people in the same deposit with an
  authoritative join in `meta_data.txt`. The problem that could not be solved by inference is simply
  not present.

✅ **The registry blocker is satisfied, not overridden.** It excluded the plasma arm because "184 rows
cover 143 participants with 41 repeat-visit rows". The baseline column is already one row per person;
the follow-up is a **separate column** that must never be stacked into the same matrix.

⛔ **Two constraints cap this lane and both are structural.** (1) No paired RNA-protein relationship
may be inferred from this accession under any lane, so nothing here may feed the RNA severity
instrument. (2) No protein champion claim is permitted — `protein_transport` needs an external
protein cohort and zero are eligible; the only other proteomic source is CODEX imaging proteomics
behind a data-sharing agreement. Development-only.

⚠ Scope limits to carry: 112 F / 31 M in a bariatric/obesity cohort, not representative; F3+F4 is
10 of 143 so the tie ceiling will bind hard; the paired interval brackets a real surgical
intervention, so arm 3 canNOT be read as a noise-dominated test-retest the way GSE193066 was.

---

## 🔍 Independent verification of Lane D — my own suspicion REFUTED, wording rescoped — 2026-09-01

`executions/verify-etiology-transfer-refutation-20260901T000000Z/` (jobs 21327145, 21327150). A
verifier that did not produce the result was briefed to refute it and defaulted to "refuted" on
anything it could not confirm.

### ⛔ My load-bearing objection was WRONG, and the measurement says so

I claimed the batch control was a training-size artifact (train n=34, test n=129). The verifier built
the control I should have specified: a **fixed** 40-participant late-block test set (marginal 11/6/23,
ceiling 0.8867), every arm trained at **n=34**, 200 draws, pipeline verbatim.

| arm, all trained at n=34 | ridge | EN | SVR | ordinal |
|---|---:|---:|---:|---:|
| CROSS34 (the real 34 early-block participants) | 0.5470 | 0.5861 | 0.6324 | 0.6448 |
| WITHIN34m (34 from LATE at the same 2/5/27 marginal) | 0.7550 | 0.7511 | 0.7537 | 0.7537 |
| **block effect, size AND marginal fixed** | **−0.2116** | −0.1669 | −0.1234 | −0.1117 |

Decomposition of the published −0.232 / −0.236 / −0.051 / −0.154:
pure training size (n=34 vs n=89, same block) −0.031 to −0.036; training marginal vs random at n=34
−0.007 to −0.022; **residual genuine block effect −0.11 to −0.21**.
⭐ **Training size accounts for about 16% of the reported loss, not the loss.** The batch control does
establish what it was used to establish.

⚠ Two caveats the verifier stated rather than hid: WITHIN34m draws from `late` without an etiology
constraint (~55% MASLD vs the early block's 71% CVH), so the measured block effect carries a residual
etiology-composition term; and at n_test=40 the individual draw intervals straddle zero, with the
medians unanimous in sign across four heads.

### ⛔ The interval defect is real and structural

Of the 132 participants entering a matched draw, **66 have sampling fraction exactly 1.0** — all 42
MASLD F4, all 13 CVH F0-2, all 11 CVH F3 are in EVERY draw and contribute **zero** variance. Two draws
share ~85% of their participants. **Exactly half the substrate contributes nothing to the quoted
interval.** The producer's own `interval_is_over_subsampling_draws_not_a_participant_bootstrap: true`
flag is correct and the number must never be read as a confidence interval.

### ✅ But `UNDERPOWERED_untestable` is too strong — the right label is a BOUNDED NULL

MDE at 80% power from the primary arm's stage-stratified participant bootstrap (10,000 paired reps):
**0.101–0.177 Spearman** — ⛔ **SUPERSEDED by V3, see the correction below. The true MDE is
0.185–0.362.** This proxy figure came from the primary arm's bootstrap standing in for the matched
arm's and understated it by about 2x. The argument built on it — that the arm is powered to detect a
block-sized loss — does not survive.

⚠ **Unasked and it cuts toward the claim:** in every etiology cell the transfer trains on MORE data
than the native it is compared against (82 vs 73, 81 vs 74; decisive 66 vs 53), which is why some cells
report `fraction_of_native_recovered > 1`. ⛔ **Nobody may claim transfer beats native.** Removing the
advantage moves the delta by 0.000–0.030 (median ~0.006) and leaves it scattered around zero and
sign-inconsistent, so the bias does not manufacture the null.

### ✅ Everything else re-derives

Substrate sha256 identical; 177 participants; 81/82/14; 19,983 − 302 all-zero = 19,681 (**0 duplicate
symbols, so the "symbol collapse" was a no-op**); outcome 39/24/114; matched marginal 13/11/42 realized
identically in all 200 draws. Leakage **zero** everywhere: MASLD∩CVH = 0, early∩late = 0, max
source∩target across 200 draws = 0. ARLD honest: 14/14 at one level, zero rows in any cell — one
footnote, the all-zero gene filter runs over all 177 rows, so ARLD does influence *which* 302 genes are
dropped; label-free and immaterial, but "enters no cell" is not literally true at that step.

### ⛔ One deposited guard is DEAD

`stage_matched_marginals_are_identical` **cannot fail**. `rng.choice(..., replace=False)` fixes the
per-stratum count by construction, so `marg_ok` can never become False; 2,000 trials produced 0
reachable mismatches, and over-requesting a stratum raises ValueError before the check runs. Its
"constructed violation" is the unmatched CVH marginal, a *different* quantity. The asserted property is
true and was verified independently — but this is another instance of
[[feedback-guard-must-be-able-to-fail-2026-08-30]] and it shipped inside a lane I described as
carefully guarded. Two further guards (`stage1_audit_pin`, `arld_enters_no_cell`) deposit **hardcoded
literals** for their violation readings rather than executed branches, though both underlying aborts
are real.

### ⭐ A scope hole I never suspected, tested and closed

The verifier suspected the transfer was **transductive** — z-scoring the target on the target's own
statistics would silently remove the etiology mean shift and require the whole target cohort in hand.
Re-standardizing on SOURCE statistics costs **−0.004 to +0.001** across all 8 cells, every interval
straddling zero. The transfer is genuinely prospective.

### The wording that may be published

> In GSE276114, at a stage-matched marginal (13/11/42, n=66 per etiology), a MASLD-trained instrument
> reaches Spearman 0.63–0.69 against a tie ceiling of 0.8545 in chronic viral hepatitis and vice versa,
> with no loss detectable at this power. **The data are consistent with the transfer losing up to
> ~0.20 Spearman, about 30% of native skill** (MDE80 0.185–0.362 by direct participant bootstrap);
> smaller losses are not excluded. This is a weak null stated as a bound, not equivalence.

⛔ **Not** "transfers with no loss." Scope: one cohort, one quantification, one 3-level source-native
outcome, two etiologies, ARLD untestable by construction, native baselines trained on less data.

### ⛔⛔ PRIORITY QUESTION — the source study may already own this claim

The verifier flags that the source study for GSE276114 (**Zeybel 2025, PMID 39889710**) is itself
titled "...advanced fibrosis in chronic liver disease **across etiologies**". If the cross-etiology
point is theirs, then per `CLAUDE.md` — *"Do not present source-owned public findings as discoveries of
this paper"* — Lane D is a reproduction on their data and not our result. **This must be settled before
Lane D enters any figure, claim or `docs/` authority.** Checking now.

⚠ Still running: V2 at 100/200 (max |median difference| 0.003 against the published decisive arm,
converging), V3 not started — V3 would give the matched arm its own participant-level SE instead of
borrowing the primary arm's. The verifier notes the matched arm is smaller so its true MDE is likely
somewhat LARGER, strengthening the underpowered direction, and explicitly labels that as inference
rather than measurement.

### ✅ PRIORITY RESOLVED — Lane D's DESIGN is novel; its DIRECTION is not, and must be cited

Source study checked directly, 2026-09-01: **Zeybel et al. 2025, *Cell Reports Medicine*,
"Integrative proteo-transcriptomic characterization of advanced fibrosis in chronic liver disease
across etiologies"**, PMID 39889710, PMC11866494, GSE276114, 330 individuals (40 healthy + 290 with
histologically characterized fibrosis from viral, alcohol or metabolic causes).

⛔ **They never ran a cross-etiology transfer.** All three etiologies were **pooled into single
training sets**; models trained on 70% of samples with stratified 5-fold cross-validation in the
discovery cohort; reported AUROC 0.98 (≥F3) and 0.99 (cirrhosis). The words *transfer*, *generalize*
and *held-out etiology* appear nowhere in connection with etiology. Their only etiology-wise analysis
is differential expression.

⚠ **But they already state the direction**, verbatim: *"The transcriptome profiles did not exhibit
clear distinctions among these etiology-based groups."*

**Therefore:**
- ✅ **Novel to us:** the train-in-one / test-in-another design, the stage-matched marginal that removes
  the composition confound, and a bounded null with a stated MDE (**0.185–0.362 Spearman**, V3).
- ⛔ **Not novel:** the direction of the answer. Cite Zeybel 2025 as anticipating it. Per `CLAUDE.md`,
  *"Do not present source-owned public findings as discoveries of this paper."* Lane D contributes the
  quantitative, prospective, confound-controlled version of an observation the source already made
  descriptively.

⭐ **Free external check on our substrate build.** Their published composition reproduces ours exactly:
they report ARLD n=14, MASLD n=42, CVH n=58 among cirrhotics and F0-2 39 / F3 25 / F4 114. Our
independent derivation found ARLD n=14 all at a single outcome level, matched-draw pools of
MASLD stage2 = 42 and CVH stage2 = 58, and outcome levels 39/24/114.
⚠ **One discrepancy to pin before publication: they report 178 liver transcriptomes, we derive 177,
and their F3 count is 25 against our 24.** One sample. Resolve it.

⚠ **Relevant to the PXD051911 plasma lane:** Zeybel 2025 also built **non-invasive plasma predictors**
of advanced fibrosis (132 circulating proteomic signatures, AUROC 0.98/0.99). That space is occupied.
Our plasma lane is a different cohort (PXD051911, bariatric/obesity, n=143), a different estimand
(graded 5-level Kleiner rather than binary ≥F3) and a different question (does routing through
predicted liver state beat direct plasma→severity). It must be differentiated on exactly those axes and
must not be pitched as a first non-invasive predictor.

---

## ⛔⛔ TWO unannotated confounds in GSE276114, both verified from the GEO record — 2026-09-01

Surfaced while resolving a bookkeeping discrepancy (178 vs 177). Neither was in view when Lane D was
designed. Both are read directly from `GSE276114_series_matrix.txt.gz`, fetched 2026-09-01.

### ✅ First, the sample count is resolved and costs us nothing

**"Liver sample 95" was withheld at submission.** GEO issued 177 contiguous accessions
(GSM8491342–GSM8491518) for 177 submitted samples; the count matrix header carries 178 tab fields =
1 `Gene.name` + 177 samples, titled "Liver sample 1".."Liver sample 178" with exactly one gap and no
duplicates. Our design merges 177 both / 0 left-only / 0 right-only with zero field mismatches, so
**our build drops nothing** — the all-zero filter drops 302 genes and never a sample.
The paper's 178 is internally consistent (Fig 1B 39+25+114; Table 1 72+14+58+34).
⭐ **Lane D is provably insensitive**: the matched marginal is min(MASLD, CVH) per stage, CVH binds at
F0-2 (13) and F3 (11) and MASLD binds at F4 (42), so restoring sample 95 (MASLD F3 [Inference], from
strict etiology block-contiguity plus arithmetic) leaves **13/11/42, n=66 per etiology, unchanged**.

⚠ **Trap to record: cite Fig 1G, never Table 1, for etiology counts.** Table 1 breaks HCC out as its
own class (MASLD 72 / ARLD 14 / CVH 58 / HCC 34); GEO folds HCC into the underlying etiology
(81/82/14). They reconcile exactly — 82 = 72 + 10 and 82 = 58 + 24, with 10 + 24 = 34 — but anyone
checking our numbers against Table 1 concludes we are wrong by nine and twenty-four.

### ⛔ Confound 1 — 34 HCC patients, folded into their etiology, UNIDENTIFIABLE

`!Series_summary` says the cohort spans *"from early-stage fibrosis to cirrhosis and **associated
HCC**."* That is the **only** occurrence of HCC/carcinoma/tumour/cancer in the entire series matrix.
GEO deposits exactly three sample characteristics — `tissue`, `disease group`, `disease` — and
**no per-sample cancer annotation exists**. The 34 split 10 MASLD / 24 CVH, so MASLD is ~12% cancer
patients and CVH ~29%.

⛔ The cross-etiology contrast is therefore partly a contrast in **HCC proportion**. Stage-matching
removed the stage confound and does nothing about this one.
⛔ Worse for the instrument: HCC is strongly associated with cirrhosis and F4 is 114 of 177, so any
tumour or tumour-adjacent signal concentrates exactly where the model has the most data.

### ⛔⛔ Confound 2 — EXPLANT versus BIOPSY, with a different RNA kit for each, also unannotated

`!Series_overall_design`: *"178 snap-frozen **explant** liver samples **and/or biopsies**..."*
`!Sample_extract_protocol_ch1`: *"...RNeasy **mini** kit for explant liver samples (74104) and RNeasy
**micro** kit (74004) for **biopsy** samples."* — and that identical string is deposited for **every**
sample, so the field cannot distinguish them either.

⛔ **This is a tissue-sampling difference AND a library-prep difference, and both are almost certainly
confounded with stage.** An explant is a whole liver removed at transplant, so explants are
end-stage and overwhelmingly F4; biopsies carry the earlier stages. **A model that appears to predict
fibrosis in this cohort may be predicting explant-versus-biopsy**, which is nearly a deterministic
function of the outcome.
⚠ This is the same failure mode that retired PRJNA512027 (library-prep x disease confound) — see
[[reference-prjna512027-retired-libprep-confound-2026-08-28]].

⭐ **Note what this does and does not touch.** The stage-matched transfer arm forces identical stage
marginals across etiology, so if sampling type is a function of stage it is matched too and the
*transfer* comparison survives. What is exposed is the weaker premise underneath it — that the model
measures fibrosis at all in this cohort rather than sampling modality. ⚠ And the verifier's measured
run-block effect (−0.11 to −0.21) may itself BE a sampling-type or kit effect rather than a batch
effect, which would change what the batch control means.

**Both questions are dispatched to the verifier. Lane D stays out of `docs/` and out of every figure
until they return.**

---

## ⭐ The delta-shape check — CORRECTED: it is a FIBROSIS contrast, and the decisive version is INTERNAL — 2026-09-01

From `executions/pxd051911-plasma-severity-20260901T191009Z`, 41 participants with two staged visits.

| cohort | n pairs | context | down / same / up | mean Δ | reading |
|---|---:|---|---|---:|---|
| GSE193066 | 58 | observation, ~2 y | 15 / 28 / 15 | **−0.034** | symmetric → noise-dominated |
| **PXD051911 fibrosis** | 41 | brackets bariatric surgery | **5 / 24 / 12** | **+0.171** | asymmetric UP |
| **PXD051911 NAS** | 41 | same | **15 / 20 / 6** | **−0.707** | asymmetric DOWN |

⭐ **Activity improves while fibrosis does not.** Coherent for a bariatric cohort — steatosis,
ballooning and lobular inflammation respond to weight loss while fibrosis lags or progresses.

⭐ **This is a genuine external check on [[project-kleiner-measurement-error-constant-2026-09-01]].**
We read GSE193066's 15-up/15-down symmetry as the signature of a measurement-noise-dominated process.
If our procedure manufactured symmetry, it would have manufactured it here too — and it did not. Two
paired-biopsy cohorts, two different signatures, **each in the direction its own clinical context
predicts.**
⚠ Not interchangeable: PXD051911's interval brackets surgery and GSE193066's does not, so neither
validates the other's magnitude. **Only the contrast in SHAPE is informative.**

### ⛔ Correction — Arm 3 canNOT test intervention response, and I was wrong to ask for it

I asked the lane to flag loudly if the paired arm showed the model tracking real surgical improvement.
**It structurally cannot.** Testing model *response* needs the score at visit 1 against the score at
visit 3, and **visit-3 plasma sits entirely in one batch covering only 14 of the 143 baseline
samples.** A within-participant plasma change is `UNTESTABLE`, not negative. The lane barred it in its
prespec (section 3, rule 1) *before* computing anything.
Arm 3b asks the strictly weaker question — does the visit-1 score forecast the visit-2 stage beyond
RTM — and a positive there is forward validity of a cross-sectional score. ⛔ **It would not contradict
[[project-intervention-response-negative-2026-08-31]] and must never be written as if it did.**

### ⛔ Two further framings of mine, withdrawn

1. ⛔ **"Our bottom-heavy case mix is the harder, more realistic distribution and that is a feature."**
   Withdrawn — it reads as spinning a null in the same artifact that reports our own incremental result
   as untestable. Report the distribution, the tie ceiling (**0.9284**, so bottom-heaviness costs ~7
   points of attainable Spearman) and the MDE, and let those carry the point without the adjective.
2. ✅ **Better sentence than the convention I supplied**: in our 143 baseline participants **≥F3 is 10
   positives and F4 is 3**, so a binary AUROC comparable to Zeybel's is *arithmetically unformable* on
   our side, not merely bad practice. One participant moves the estimate materially.

**Arm 1 label, assigned mechanically from the deposited verdict rather than by judgment:**
`REPRODUCTION_ATTEMPT_increment_over_age_sex_BMI_UNTESTABLE_at_this_n`. Beating a permutation null only
reproduces the already-published fact that plasma carries fibrosis information. Reported as a
reproduction and as a control for Arm 2, never as a finding.

### Prior-art verification granularity

✅ **Verified by direct fetch of PMC11866494** (2026-09-01): pooled training across all etiologies, 70%
of samples, stratified 5-fold CV; AUROC 0.98 (≥F3) / 0.99 (cirrhosis); **no** cross-etiology transfer
and no occurrence of transfer/generalize/held-out etiology in that connection; liver transcriptome
n=178 with cirrhotic ARLD 14 / MASLD 42 / CVH 58; the "did not exhibit clear distinctions" quote;
title, Cell Reports Medicine, 2025, PMID 39889710.
⚠ **[Unverified]**, from a search-result summary only: the n=330 cohort total and the "132 circulating
proteomic signatures" figure. The constraint on our framing rests on the verified items.

### ✅ The explant-versus-biopsy worry is NOT supported — measured — 2026-09-01

Run directly on `data/GSE276114/GSE276114_raw.count.txt.gz` (19,983 x 177, joined 1:1 to the design on
`sample_number`, 177/177 asserted). Features **technical only**, no biological gene selection: library
size, detected genes, zero fraction, MT fraction (13 MT genes), and each sample's top-50 count share.

| predictor | Spearman vs fibrosis stage |
|---|---:|
| technical-only ridge, 10-fold CV x 5 seeds | **0.163** (per-seed 0.183/0.168/0.114/0.198/0.150) |
| MT fraction alone | 0.241 |
| 2-means technical cluster (sizes 145/32) | 0.080, against an **MDE of 0.211** at n=177 |
| **native biological ridge, for scale** | **0.72–0.79** |

Permutation null for the technical model: mean −0.017, 97.5% +0.164, **p = 0.028**. The technical
cluster is independent of `disease_group` (χ² p = 0.148) and of etiology (p = 0.186), and **both
clusters span sample_number 1–178**, so it does not align with the early/late blocks either.

⭐ **A small technical association with stage is real (p = 0.028) — expected, since sicker tissue
degrades — but it cannot manufacture a 0.72–0.79 signal.** The specific failure mode, *the model is
really predicting explant-versus-biopsy*, is not supported.

⚠ **Two limits, stated rather than glossed.** (1) Five crude proxies. This shows that whatever
separates explants from biopsies is not captured by them and is not strong enough to generate the
stage signal — **not** that the sampling types are balanced. (2) **MT fraction is the one feature that
genuinely tracks stage** (ρ +0.261, p = 4.4e-4, above the MDE; medians 0.041 / 0.042 / 0.061 across
F0-2 / F3 / F4) and it is ambiguous between biology (mitochondrial dysfunction is central to MASH) and
tissue quality. It also differs by etiology: **MASLD median 0.066 vs CVH 0.048.**

⭐ **New control commissioned because of that last number.** MT fraction differs by etiology AND tracks
stage, which is exactly the nuisance that could carry part of a cross-etiology transfer. The decisive
arm is to add a **technical-features-only baseline** on the same folds, matched draws and bootstrap
indices, and report the transfer as a paired delta over that as well as over the native baseline. The
claim then becomes "not distinguishable from native **and** not explained by library quality." If the
technical baseline transfers as well as the expression model, we want to find that ourselves.

---

## ⭐ HCC confound RESOLVED — real, unremovable, and much smaller than it first looked — 2026-09-01

### The finding that de-escalates it: the tissue is PERITUMORAL

Read from the paper's Methods directly, not inferred:
> **"For peritumor HCC samples, the samples were collected at least 1 cm distant from the tumor."**
> **"The fibrosis stage in the peritumoral liver tissue of HCC patients ... were histologically
> evaluated by expert liver pathologists based on Kleiner et al."**

⭐ **The sequenced tissue is not tumour, and the outcome is staged on the same tissue that was
sequenced.** Exclusion criteria removed "any malignancy other than HCC". Fig 1B: 178 = CLD 144 +
CLD-with-HCC 34.

✅ **The data agree with the protocol.** An attempt to identify the 34 from expression fails in an
informative way: **AFP is uniformly low (median logCPM 1.21, max 2.44, ZERO samples above 3)**, GPC3
flags 2, and no marker in a 17-gene HCC/proliferation panel yields anything like 34 split 24/10.
That is exactly what peritumoral tissue ≥1 cm from tumour should look like. **There is no measurable
tumour signal in this tissue to learn.** My "the model may be learning cancer" worry is substantially
lowered — though not to zero, since a field effect is not excluded.

### Identification: impossible. Every route exhausted

GEO series matrix (three characteristics only), **BioSample** (SAMN43432520/21/620 checked directly),
SRA (BioSample-derived), Supplementary Data S1–S7 (all analysis outputs, no per-sample table),
Document S1 (PMC returns **HTTP 520**, not retrievable), Mendeley 6brkvh3f97 (1,460 assays x 218
"Subject N", zero phenotype rows). ✅ **Our design file discarded nothing — it kept every informative
field GEO has.**

### Exposure, quantified from the paper's own numbers

> "70.6% have CVH-associated HCC, and 29.4% have MASLD-associated HCC" = **24 CVH-HCC, 10 MASLD-HCC,
> ARLD-HCC = 0.**

| group | HCC exposure in the deposited 177 |
|---|---|
| CVH | **24 of 82 = 29.3%** (certain) |
| MASLD | **9 or 10 of 81 = 11.1–12.3%** (depends on the withheld sample) |

⭐ **This independently CONFIRMS the sample-95 attribution.** Totals including HCC are MASLD 72+10 = 82,
CVH 58+24 = 82, ARLD 14, sum 178; GEO deposits MASLD 81 / CVH 82 / ARLD 14, so **the withheld sample
must be MASLD**. Upgrade from [Inference] to confirmed — now supported by the paper's etiology
arithmetic as well as by block contiguity.

### ✅ The matched arm does NOT make it worse

Under the worst-case assumption that every HCC is F4 (biologically realistic — HCC arises in cirrhosis
and staging was on peritumoral tissue): MASLD 10/66 = **15.2%** (raw 12.3%, so the matched arm
*enriches* MASLD HCC by dropping F0-2 from 26 to 13); CVH expected 17.4/66 = **26.4%** (raw 29.3%,
slightly depleted). **Gap narrows from 17.0 pp to 11.2 pp.** Arithmetic under a stated assumption, not
a measurement — per-sample HCC is unavailable.

⚠ Weak supporting evidence, labelled weak: a peritumoral field effect mimicking fibrosis should inflate
apparent skill in the HCC-richer group, so CVH→MASLD should transfer *worse* than native. Observed
matched deltas −0.024 to +0.031, all straddling zero. Consistent with, not proof of, a small effect.

**Required wording:** the contrast is between two clinically defined groups that differ in etiology
**and** in HCC burden (~12% vs ~29% with concurrent HCC, tissue sampled peritumorally ≥1 cm from
tumour and staged on that same tissue). It is not a clean etiology contrast and cannot be made into one
from public data.

### ✅ V2 landed — EXACT reproduction

**Max |difference| between recomputed and deposited medians is 0.000000** for both transfer and native
across all 8 cells at the full 200 draws. Size-matched deltas (source cut 66→53 to match the native's
own median training n): **no cell excludes zero**; largest shift is elastic_net MASLD→CVH
(−0.008 → −0.039). **The training-size advantage does not manufacture the null.**

---

## ⛔ Only 2 of 9 severity cohorts can be affirmatively cleared of folded-in cancer — 2026-09-01

Metadata scan across all `!Series_*`/`!Sample_*` lines for hcc/carcinoma/tumor/cancer/malignant/
neoplasm/adjacent/resection: **zero hits** in GSE126848, GSE130970, GSE135251, GSE162694, GSE167523,
GSE213621, GSE240729; 1 boilerplate hit in GSE174478; 2 in GSE193066.

⛔⛔ **The method's limit, and it is the whole point: GSE276114's own GEO metadata and abstract are
equally silent about its 34 HCC patients. A zero-hit metadata scan is exactly the evidence that FAILED
to catch the cohort we are worried about.** Clearing a cohort requires reading its source paper's
cohort table.

Against open-access full texts (6 of 9):
- ✅ **GSE193066 — affirmatively clean**: explicit **"106 HCC-naive NAFLD patients."** (⚠ deliberately
  enriched for HCC *risk*, 71 of 106 PLS-high-risk — a selection characteristic, not prevalent cancer.)
- ✅ **GSE167523 — near-affirmatively clean**: liver cancer appears only as an *incident follow-up*
  outcome, which requires patients to be cancer-free at biopsy.
- ⚠ **GSE240729 — needs a targeted follow-up.** One ambiguous sentence: *"Considering that HCC or
  cirrhosis is often an exclusion criterion in biomarker studies, these findings further support the
  biomarker panel's broad applicability."* That could imply their cohort did **not** exclude HCC. GEO
  carries only fibrosisscore F0–F4 (6 at F4). Unresolved.
- SILENT (cancer in discussion only): GSE130970, GSE162694, GSE174478, GSE213621.
- SILENT (not open access, abstract only): GSE126848, GSE135251.

**Verdict: no positive evidence of another folded-cancer cohort, but 7 of 9 are CANNOT-DETERMINE.**

⚠ **Species trap found in passing:** GSE213621's series title and `overall_design` describe the
paper's **mouse** arm while its 368 deposited samples are titled "human liver sample ...". Our build is
correct, but this is the same family as the Table-1 trap.

---

## ⛔ CORRECTION to the delta-shape check — my framing overstated it, twice — 2026-09-01

### ⛔ NAS does NOT contrast. Only fibrosis does.

Pulled from both deposits rather than retyped
(`kleiner-reliability-20260901T122243Z/results/kleiner_reliability.json`):

| cohort | fibrosis | NAS |
|---|---|---|
| GSE193066 (n=58) | 15 down / 28 same / 15 up, mean **−0.034**, sign p **1.000** | **34 down / 17 same / 7 up, mean −1.121** |
| PXD051911 (n=41) | 5 down / 24 same / 12 up, mean **+0.171** | 15 down / 20 same / 6 up, mean −0.707 |

⛔ **NAS falls in BOTH cohorts.** My sentence — "two cohorts, two different signatures, each in the
direction its clinical context predicts" — is **true for fibrosis and false for NAS**. Any claim in
that shape must be scoped to fibrosis explicitly or it overstates the evidence.

### ⭐ The decisive version was already INTERNAL to GSE193066, and it is stronger

In the **same 58 participants**, with the **same pairing**, the **same delta definition** and the
**same sign test**: fibrosis 15/15 (p = 1.000) against NAS 34 down / 7 up.

⭐ **The procedure reported a large asymmetry in one outcome and none in the other, on identical
participants and identical pairing.** A pairing defect, or any artifact that manufactures symmetry,
would have flattened both. It flattened neither by construction. **So the fibrosis symmetry is
informative rather than structural — and that conclusion needs no second cohort at all.**

### What PXD051911 actually adds — narrower, but only a second cohort can supply it

It replicates that capability under a **completely different pairing mechanism**. GSE193066 pairs by
regex on sample titles (`re.sub(r'_\d+$', '', title)`); PXD051911 pairs by an explicit metadata join on
`patient_name` between `initial_sample` and `follow_up_sample` rows, with no title parsing anywhere.
That retires the narrower worry that **the title regex specifically** was mispairing participants —
which the internal NAS check cannot address, because it shares the pairing.

**Correct order to write it: lead with the internal fibrosis-vs-NAS dissociation, and use PXD051911 as
the pairing-mechanism replication behind it.**

### What the deposit explicitly refuses to license (`results/DELTA_SHAPE_CROSS_COHORT.json`)

- ⛔ Does **not** show the GSE193066 pairing was correct — only that the procedure *can* report asymmetry.
- ⛔ Does **not** confirm the noise reading. Symmetry excluded a net upward drift, never balanced real change.
- ⛔ Does **not** license comparing magnitudes. `sd_delta`, σ_e and ICC are not commensurable across the
  two: different intervals, different selection (all 58 GSE193066 pairs high-risk at first biopsy), and
  PXD051911 brackets surgery.
- ⛔ Does not revise the constant **as of that date**. ⚠ **LATER SUPERSEDED**: σ_e 0.6433 [0.5170, 0.7600] is the GSE193066 value only; the shipped quantity is the range **0.4655–0.6433**.

### One prior-art constraint I had not drawn

Zeybel's AUROCs are **pooled-etiology discovery-cohort** numbers under 70/30 with 5-fold CV, and no
cross-etiology transfer is claimed. So the occupied ground is *"plasma carries fibrosis information and
can be fitted to a binary threshold WITHIN a cohort"* — **exactly what our Arm 1 reproduces**. It is
not a held-out-cohort transfer result and our lane must not describe it as one.

---

## ⛔ V3 RETRACTS the verifier's own power verdict — the etiology null is WEAKER than reported — 2026-09-01

Direct participant-level bootstrap of the matched arm (300 replicates, stratified by etiology x stage,
**CV folds assigned by unique participant** so a duplicated participant cannot sit in both a training
fold and its own test fold). Reproduction guard: max |diff| from deposited medians = **0.000000**.

| cell | boot Δ | 95% CI | SE | subSD | ratio | **MDE80** | % native |
|---|---:|---|---:|---:|---:|---:|---:|
| CVH→MASLD ridge | +0.0124 | [−0.124, 0.188] | 0.0748 | 0.0432 | 1.73 | **0.2096** | 30.4% |
| CVH→MASLD elastic_net | −0.0084 | [−0.150, 0.103] | 0.0668 | 0.0419 | 1.59 | 0.1870 | 27.7% |
| CVH→MASLD linear_svr | +0.0796 | [−0.121, 0.323] | 0.1158 | 0.0765 | 1.51 | 0.3244 | 50.8% |
| CVH→MASLD ordinal | −0.0112 | [−0.140, 0.135] | 0.0661 | 0.0426 | 1.55 | 0.1852 | 26.4% |
| MASLD→CVH ridge | +0.0282 | [−0.142, 0.256] | 0.0988 | 0.0433 | 2.28 | **0.2768** | 42.4% |
| MASLD→CVH elastic_net | +0.0099 | [−0.173, 0.181] | 0.0919 | 0.0704 | 1.31 | 0.2576 | 39.7% |
| MASLD→CVH linear_svr | +0.1171 | [−0.141, 0.378] | 0.1292 | 0.0666 | 1.94 | 0.3620 | 64.3% |
| MASLD→CVH ordinal | −0.0163 | [−0.199, 0.157] | 0.0921 | 0.0487 | 1.89 | 0.2581 | 38.5% |

⛔ **The published subsampling interval is 1.31x to 2.28x TOO NARROW.** The structural argument (66 of
132 participants at sampling fraction 1.0) is now quantified.
⛔ **MDE80 is 0.185–0.362 = 26–64% of native skill**, not the 0.101–0.177 (14–28%) taken from the
primary-arm proxy. **The proxy understated it by about 2x.**
⛔ **Worst-case loss still consistent with the data: −0.1995 Spearman.**

**What this retracts.** The verifier had argued against `UNDERPOWERED_untestable` on the grounds that
the arm is powered to detect a loss the size of a run-block shift (0.11–0.21). For the prespecified
primary head (ridge) MDE80 is **0.210 and 0.277 — at or above the block effect.** That argument does
not survive. It had been labelled inference rather than measurement, and the measurement moved it far
enough to change the conclusion, not merely the number.

⭐ **Where the label lands.** Still not `UNDERPOWERED_untestable` — the arm produces a real, non-null
transfer (0.63–0.69 against a 0.8545 ceiling, permutation p = 1e-4). But **"no loss" is wrong and even
"no loss at this power" needs the number attached: the data are consistent with the transfer losing up
to ~0.20 Spearman, about 30% of native skill.** A weak null, stated as a bound.

---

## ✅ V6 — the explant/biopsy confound is largely cleared, in the way that mattered

**Procurement, from Methods:** three routes — transplantation (explant), resection surgery, and
diagnostic percutaneous core biopsy. Two RNA kits, mini for explants and micro for biopsies. ⭐ **But
library prep is IDENTICAL for every sample** (Illumina Stranded Total RNA + Ribo-Zero Plus). So this is
an extraction-kit and tissue-type confound, **not** the library-chemistry confound that retired
PRJNA512027 — weaker in kind than that precedent.

- **Is sampling type recoverable?** A 2-group structure exists but is weak and unlabellable: KMeans k=2
  on 10 purely distributional features splits 117/60, with Cramér's V 0.208 vs stage, ρ 0.204,
  χ² p 0.021. **Nothing like the deterministic explant-equals-F4 relationship I feared.**
- ⭐ **Does the block effect coincide with it? NO, and this is clean.** Multivariate CV AUROC predicting
  late-vs-early block from all 14 technical features = **0.4941, chance.** Per-feature 0.50–0.61; only
  mt_fraction reaches p=0.045 uncorrected, 1 of 14 tests, not surviving correction. **The run-block
  boundary is technically invisible, so the batch control still means what we read it as meaning.**
- ✅ **My rescue argument's premise was TESTED, not assumed.** Within F4 only (n=114) the technical
  cluster is etiology-independent (MASLD 0.524 vs CVH 0.535, χ² p=0.195), so stage-matching does
  approximately technical-match. ⚠ Honest limit: this can only detect a route difference if the
  technical cluster tracks sampling route, and that cannot be verified because the labels do not exist.
  **State it as "no detectable difference", never "no difference".**
- ⚠ And granting the premise rescues only the transfer comparison's **internal validity**, not the
  interpretation. If the skill were technical, "the signal transfers across etiology" would collapse to
  "a kit signature transfers across etiology" — true and vacuous.

---

## ⭐⭐ The finding that matters more than the kit: summary statistics recover most of the model in CVH

Within-etiology CV Spearman against fibrosis stage (**ignore the pooled 0.507 — it is inflated by
pooling etiologies with different stage marginals plus all-F4 ARLD**):

| feature set | MASLD | CVH |
|---|---:|---:|
| all 14 technical summaries | 0.418 | **0.527** |
| distributional only (kit-sensitive) | 0.406 | 0.093 |
| **expression model, for scale** | **0.719** | **0.623** |

⭐ **In CVH, 14 summary statistics recover 0.527 of the expression model's 0.623 — about 85%.** The
kit-sensitive distributional features are weak there (0.093), so **this is not a kit artifact.** The
signal comes from two compositional features with straightforward biological readings:

| feature | ρ overall | MASLD | CVH | reading |
|---|---:|---:|---:|---|
| `albumin_fraction` | **−0.320** | −0.22 | −0.37 | falls with fibrosis — hepatocyte depletion |
| `mt_fraction` | **+0.261** | +0.22 | +0.48 | rises with fibrosis |

⛔ **This is a real statement about what the instrument measures: much of it is coarse compositional
shift, not a fibrogenesis program.** It is consistent with
[[project-expression-beyond-composition-graded-fibrosis-2026-08-28]] (C1 +0.17329 beyond composition —
composition explains much, not all) and with the project's standing position on de-emphasizing
hepatocyte-intrinsic signal. **Worth recording independently of Lane D.**

### ⛔ CORRECTION: albumin is NOT the driver — library COMPLEXITY is

Measured under one protocol (ridge, 10-fold, 5 seeds, pooled n=177), by the verifier:

| feature set | n feat | pooled ρ |
|---|---:|---:|
| the original 5 | 5 | 0.1497 |
| distributional | 10 | 0.4423 |
| all 14 | 14 | 0.4879 |
| **all 14 minus albumin** | 13 | **0.4395** |

⛔ **Dropping `albumin_fraction` costs 0.048, not 0.34.** The 5-feature set reproduces exactly
(0.1497 vs 0.163 — no data discrepancy). The gap is five *distributional* descriptors: top-10 and
top-500 shares, Shannon entropy, Gini, and the median/IQR of logCPM. **The reading "the signal comes
from two compositional features" above is wrong and is retracted.**

⚠ The **85% figure is within-etiology CV in CVH only** (0.527/0.623). Against the *transfer* Spearman
the technical baseline reaches **52–70%**. Always state which denominator.

---

## ✅ Cancer-exposure audit of the nine severity cohorts — 0 contaminated, 1 cleared, 8 undetermined — 2026-09-01

`executions/cohort-cancer-exposure-audit-20260901T200423Z/`. Method: **only a source publication's
cohort description, inclusion/exclusion criteria or participant table can clear a cohort.** GEO
family.soft.gz + PubMed esummary + PMC/Europe PMC full text + publisher PDFs. GEO web pages are
recaptcha-blocked to `wget` here; the FTP SOFT route works and was used.

| verdict | cohorts |
|---|---|
| CONTAMINATED | **none** |
| AFFIRMATIVELY_CLEAN | GSE193066 |
| CANNOT_DETERMINE | GSE126848, GSE130970, GSE135251, GSE162694, GSE167523, GSE174478, GSE213621, GSE240729 |

⭐ **The positive control proves the method's necessity.** Run against known-contaminated GSE276114
(34 HCC, 19%): **sample-level cancer-keyword hits = 0**, series-level = 0. The only word that flags it
is **`explant`** — 177 occurrences — which is a tissue-source term, not a cancer term. **A keyword scan
would clear a cohort that is one-fifth cancer patients.**

⭐ **`explant` is the tell, not cancer vocabulary.** Transplant tissue implies end-stage disease and
routine HCC surveillance. Grep for procurement route, not for oncology words.

### ✅ GSE240729 alarm RESOLVED — it does not implicate the deposited data

The flagged sentence is preceded by *"The presence of HCC did not influence the performance of the
panel either"*, so HCC patients **are** present in that study — but the paragraph is explicitly scoped
to **the n=128 serum TESTING cohort** (68 ANCHOR, 22 BARICO, 38 Leiden UMC), which has no RNA-seq and
**contributes zero samples to GSE240729**. GEO holds only the 67 FFPE liver-biopsy translation-cohort
samples (F0 9 / F1 17 / F2 25 / F3 10 / F4 6). [Inference], strong: panel performance was evaluated
only in the serum cohorts.

### ⛔ GSE167523's near-clearance is OVERTURNED

My earlier reading — that liver cancer appears only as an incident follow-up outcome, so patients must
have been cancer-free at biopsy — **does not transfer**. The incidence analysis is scoped to *"all 164
patients ... at Kaizuka City Hospital"*, while the **98 RNA-sequenced samples come from Sendai Kousei
Hospital, 2016–2018**. Different arms of the same 311-patient multicentre study. ⚠ Residual
reassurance only: the four enrolment criteria describe a diagnostic, non-oncologic biopsy indication,
and no cancer exclusion is stated for any of the 311.

### ⛔⛔ NEW TRAP — GSE193066's PARENT SuperSeries IS an HCC resource

GSE193066 is clean: *"106 **HCC-naive** NAFLD patients (164 biopsy tissues...)"* and *"106
non-cirrhotic HCC-naive patients ... during which 6 patients developed HCC"* — incident, not prevalent.

⛔ **But it is a SubSeries of SuperSeries GSE193084**, whose sibling SubSeries **GSE192959, GSE193080
and GSE200460 are HCC cohorts**: the PLS-NAFLD derivation set is *"48 patients who underwent curative
... radiofrequency ablation (RFA) for early-stage NAFLD-related HCC"*, and tissue validation set 2 is
*"59 **HCC-experienced** patients with NAFLD who underwent complete tumor resection."*

⛔ **Never ingest GSE193084 or any sibling SubSeries.** This matters more than the others because
GSE193066 is the cohort behind **both** the Kleiner σ_e constant and the Test 2 forward-validity
result. Widening to the SuperSeries would silently import HCC patients into the two results the
campaign leans on hardest.

### What the eight undetermined verdicts actually mean

Their source publications state exclusion criteria that **simply do not mention malignancy** — not that
cancer patients were included. **Per-sample cancer status is recoverable for 0 of 9 cohorts**; none
carries a cancer field in GEO. So the substrate has an unquantified, probably small, cancer exposure
that cannot be measured from public data and should be stated as a limitation rather than resolved.

---

## ⛔ The panel-size curve cannot select a panel size — 2026-09-01

Read directly off `executions/panel-and-ood-20260901T135507Z/results/ood.json`, not off summary prose.

| k | transfer ρ | paired Δ | excludes zero |
|---:|---:|---:|:--|
| 5 | 0.6535 | 0.0443 | no |
| 10 | 0.6713 | 0.0619 | no |
| **15** | **0.6864** | 0.0770 | **yes** |
| 20 | 0.6740 | 0.0646 | no |
| 25 | 0.6762 | 0.0668 | no |
| **30** | 0.6804 | 0.0709 | **yes** |
| 40 | 0.6725 | 0.0630 | no |
| 50 | 0.6259 | 0.0166 | no |
| 100 | 0.6512 | 0.0423 | no |
| 200 | 0.6411 | 0.0322 | no |

⛔ **The curve is non-monotone and the two arms that exclude zero are not adjacent.** k15 and k30 clear
while k20 and k25, which sit between them, do not. Every MDE in the column is ~0.105–0.120 and every
point estimate is below it except those two. **This is the shape of noise around a flat curve, not a
size–performance relationship.** Selecting k from this table would be selecting on the noise.

⛔ **Consequence for the release**: a small panel may be shipped as a *demonstration* that most of the
signal survives aggressive gene reduction. It may not be shipped as a recommended assay size, and no
text may imply that 5, 15 or 30 was chosen because the curve favoured it.

⚠ Also note every panel arm reports `fraction_of_in_cohort_signal_recovered` above 1.0 (1.03–1.13).
Transfer exceeding the native in-cohort baseline is a statement about the *baseline*, not only about
the model — the native number is a 10-fold in-cohort ridge at 0.609434, and it is itself an estimate
with sampling error. **"Transfers above native" must always be written with the baseline named.**

---

## 🔍 Scrutiny of THE HEADLINE — it survives, and two of its sentences do not — 2026-09-01

Re-derived from `executions/latent-vs-label-20260901T153025Z/results/test2_partial.json` by a reader
who did not produce it. **The result stands.** β +0.4174 [0.2159, 0.5992], perm p 1.0e-04, jackknife
all-positive, content-free control +0.0881. Three separate things push the estimate DOWN, not up:

1. `shared_error_caveat`, in the artifact's own words: *"M1 and Y1 come from the SAME biopsy, so the
   biopsy-sampling component of the error is shared and cancels from the residual."* M1 therefore
   carries biopsy 1's own sampling error, and biopsy 2 samples elsewhere, which drags β toward zero.
2. `selection_robustness` is correct that selection on a **regressor** does not bias the coefficient —
   so the PLS-selection confound that damages Test 1 does not damage Test 2.
3. The paired set spans stages 1–3 with one stage 4, compressing the outcome range.

⛔ **But that same caveat narrows the claim, and the narrowing was not carried into the prose.** If
the sampling component cancels, this test can only detect the **reading** component of the error. So
the demonstrated statement is *the model corrects the slide read*, **not** *the model corrects the
needle*. Whether it also beats biopsy sampling variability is **untested here** — and biopsy sampling
is the part clinicians most want fixed. σ_e folds reader, sampling and real change into one term; this
test reaches only one of the three.

⛔ **"Clinically actionable: it says which biopsies to distrust" is not earned.** The direction-only
statistic a clinician would actually use is the pre-registered Y1-adjusted **AUC 0.7289, CI95
[0.5150, 0.9060]** — a lower bound of 0.515 is a hair above coin-flip. Unadjusted AUC is **0.4489**
and raw sign concordance is **0.500**. The regression coefficient is the sound claim; the AUC is the
clinical claim, and it is barely separated from chance at n=30 discordant. **Write the regression,
drop the actionability sentence, and report the AUC with its interval whenever direction is mentioned.**

### ⛔ The campaign's single largest structural risk, stated plainly

**σ_e = 0.6433, the ICC comparison 0.7167 vs 0.1886, and this headline are all the same 58 pairs from
the same deposit, GSE193066** — selected high-risk, stages 1–3, one stage 4 in 164 deposited samples.
Three flagship results, one cohort, n=58. No amount of internal validation fixes that.
✅ A second repeat-biopsy substrate has now been located and is under test — see the PXD051911
replication lane. Until it reports, every one of the three should be written as single-cohort.

---

## ✅ V8 — the technical confound is REFUTED by measurement: expression wins 8/8 — 2026-09-01

`executions/verify-etiology-transfer-refutation-20260901T000000Z/v8_results.json`. Job 21327444
(52:28, exit 0). Same matched draws (seeds SEED+5000+d), same heads, same CV rule, same code path;
only the design matrix changes, so every delta is paired.

| technical feature set | median technical transfer | median Δ (expression − technical) | cells expression wins, interval excludes 0 |
|---|---:|---:|---:|
| the original 5 | 0.2566 | **+0.4205** | **8/8** |
| distributional 10 | 0.3981 | **+0.2668** | **8/8** |
| all 14 | 0.3971 | **+0.2728** | **8/8** |

Honest participant bootstrap (300 reps, ridge, grouped folds), all four excluding zero:

| set | direction | expression | technical | Δ | CI95 |
|---|---|---:|---:|---:|---|
| 5 | CVH→MASLD | 0.678 | 0.117 | +0.556 | [0.313, 0.847] |
| 5 | MASLD→CVH | 0.683 | 0.266 | +0.408 | [0.077, 0.853] |
| 14 | CVH→MASLD | 0.678 | 0.298 | +0.385 | [0.103, 0.693] |
| 14 | MASLD→CVH | 0.683 | 0.386 | +0.283 | [0.067, 0.716] |

⭐ **My library-quality objection is dead, and it died by measurement rather than argument.** The
cross-etiology transfer is not explained by technical variation.

⚠ **Two things to disclose anyway.** The technical baseline is **not negligible** — 52–70% of the
expression transfer's Spearman from 10–14 numbers. And it is **asymmetric**: technical transfers into
CVH (0.35–0.47) far better than into MASLD (0.14–0.36), tracking `mt_fraction`'s ρ +0.48 in CVH vs
+0.22 in MASLD.

## ⚠ V7 — residualization, and why it settles less than it appears to

| variant | median transfer | median Δ vs native | cells excluding zero |
|---|---:|---:|---:|
| published, no residualization | 0.6795 | +0.0012 | **0/8** |
| residualized on 10 distributional | 0.4512 | −0.1188 | **1/8** |
| residualized on all 14 | 0.2442 | −0.2113 | 4/8 |

⛔ **"Residualization renders the transfer lossy" is NOT supported.** In the clean 10-feature variant
**7 of 8 intervals cover zero**, and the sole exclusion is `CVH_to_MASLD__linear_svr` (−0.284) — the
learner already recorded as unstable (SVR on unscaled PCs is superlinear in C), which is also the most
extreme arm in all-14 (−0.382) and the lowest transfer anywhere (0.046). One fragile learner dissenting
in the direction its own instability predicts is not a finding.

⛔ **The all-14 variant cannot carry the claim** — entropy, Gini and top-N shares are computed *from*
the counts and are partly biological (a cirrhotic liver genuinely has a more skewed transcriptome).
**0.231 bounds the technical contribution from above; it does not estimate it.** Never write "34% of
the signal is technical."

⛔ **Residualization costs power, so it corroborates nothing either.** CI half-widths go from ~±0.09
published to ±0.20 residualized (⛔ I first wrote ±0.25, which was the MAX of eight arms, not the median). The residualized test is a *weaker* bounded null, not support.

⚠ **Unchecked**: the residualization was fitted pooled across all 177 participants. If the technical
summaries correlate with etiology — and the mt_fraction asymmetry says they do — pooled residualization
strips genuine between-group structure and biases the transfer comparison *toward* finding no
difference. Flagged to the verifier; recorded as unchecked, not as benign.

✅ **What survives, and it is the right reframing**: the transferable part of the signal is
disproportionately the coarse compositional axis, so Lane D is **generality of a composition-dominated
severity axis, not of a fine-grained fibrogenesis program.**

---

## 📌 DECISION (user, 2026-09-01): Lane D is BENCHMARK-ONLY — it never enters `docs/`

**Not `docs/RESULTS.md`, not `docs/PAPER.md`, not any figure, not as a validated update.** Scope
reasoning: `docs/` holds the Resource paper's authorities, and a modeling result about etiology
generality belongs to a different paper. The evidence is unchanged; only the destination is settled.

**Second decision: the flagship write-up waits on the PXD051911 replication.** σ_e = 0.6433, the ICC
0.7167-vs-0.1886 gain, and the "predicts which biopsy is wrong" headline all rest on the same 58 pairs
from GSE193066. Nothing is written up until the ~27 independent at-surgery pairs report.

### Lane D FINAL wording, for the benchmark record (corrected, verbatim)

**(1)** In GSE276114, much of what a bulk expression severity model appears to measure is coarse
compositional rather than fine-grained transcriptional variation: fourteen per-library summary
statistics recover Spearman 0.40–0.53 against fibrosis stage within etiology, against 0.62–0.72 for
the full expression model, and projecting ten purely distributional library summaries out of every
gene removes 16–34% of both native and cross-etiology performance, the split between the two arms
depending on whether the residualization is fitted across the cohort or within etiology.

**(2)** Downstream of that, at a stage-matched marginal (F0–2 13, F3 11, F4 42; n = 66 per etiology),
a MASLD-trained model predicts fibrosis stage in chronic viral hepatitis at Spearman 0.63–0.69 against
a matched tie ceiling of 0.854, and a CVH-trained model does the same in MASLD, with the paired delta
against a native model covering zero in all eight direction-by-learner arms.

**(3)** This is a bounded null rather than an equivalence — a participant-level bootstrap gives 80%
power only for a loss of 0.19–0.36 Spearman (26–64% of native skill), so a loss of up to approximately
0.20 Spearman is not excluded — and the technically residualized re-test is **unidentified rather than
merely imprecise**: its paired delta reverses sign with the scope of the residualization (−0.119 fitted
across the cohort, +0.128 fitted within etiology; median paired shift 0.253, positive in every arm),
because pooled fitting removes structure shared between the etiologies while within-etiology fitting
removes precisely the within-group variation a native model exploits, so neither specification is
privileged and the residualized comparison cannot adjudicate the etiology question in either direction.

**Accompanying note (methods/supplement).** Under distributional residualization the paired point
estimates are consistently negative — seven of eight arms negative, median paired delta −0.119 — but
seven of eight intervals cover zero, so a loss is not established and this comparison cannot
adjudicate it. The single exclusion is the linear-SVR arm (−0.284, 95% interval [−0.540, −0.029]), the
learner recorded elsewhere as unstable (SVR on unscaled principal components is superlinear in C),
which is also the most extreme arm under full residualization and has the lowest transfer of any cell
(0.046). Dropping both SVR arms leaves the median paired delta unchanged at −0.119, so the negative
direction is not SVR-driven, but no arm other than SVR excludes zero and the result must not be read
as one of four learners reaching significance. **Paired statistics are medians of per-draw differences
throughout; differences of marginal medians are not reported, as the two diverge here (−0.119 against
−0.004).**

⛔ **These figures are specific to residualization fitted across the full cohort.** Refitting within
etiology reverses them: **zero of eight arms negative, median paired delta +0.128**, three arms
excluding zero in the POSITIVE direction. Neither specification is privileged, and the pair is
reported together rather than either quoted alone.

**Scope limit (a), HCC burden.** Concurrent hepatocellular carcinoma differs between the compared
groups: 24 of 82 participants with chronic viral hepatitis (29.3%) and 9 or 10 of 81 with MASLD
(11.1% or 12.3%). The one-participant range follows from a discrepancy in the deposited record: the
source study reports 34 participants with HCC in a 178-participant cohort (24 CVH, 10 MASLD), but only
177 were deposited to GEO, and the omitted participant is MASLD of unknown HCC status. Per-sample HCC
status is not deposited and could not be recovered from the series matrix, BioSample, SRA, the
published supplementary data, or the associated Olink deposit, so it cannot be excluded or adjusted
for. Tissue from participants with HCC was peritumoral, collected at least 1 cm from the tumour and
staged histologically on that same tissue, and shows no AFP elevation (maximum log2 CPM 2.44, no
sample above 3). AFP is a single marker and does not exclude field effects in peritumoral tissue; the
peritumoral sampling and the ≥1 cm margin are the primary basis for this limit.

**Scope limit (b), sampling route and kit.** Sampling route (transplant explant, surgical resection,
percutaneous core biopsy) and RNA extraction kit are not annotated per sample; library preparation was
identical for all participants. A technical-summary model does not reproduce the fibrosis signal (the
expression model exceeds it in all eight arms, paired delta +0.18 to +0.55) and does not track the
run-block boundary used as a positive control (cross-validated AUROC 0.494), but sampling-type balance
between etiologies is undetected rather than demonstrated.

### ⛔ Three defects the verifier corrected — two of them in text I had ALREADY ACCEPTED

The paired-vs-marginal error recurred twice more after the first catch, in passages I had signed off:

| where | wrong | correct |
|---|---|---|
| Correction-2 reasoning | difference of medians **−0.0037** | median of paired differences **−0.1188** |
| cell count | "6 of 8 negative" | **7 negative, 1 positive** |
| sentence (1) | "a third of **both**" | **a quarter native (0.155), a third transfer (0.231)** |
| sentence (3) | "triples ... ±0.25" | **more than doubles (2.27x) ... ±0.20** |

⭐ **My own ±0.25 was the LARGEST arm, not the median** — I re-derived the eight half-widths myself:
0.1558, 0.1825, 0.1922, 0.1948, 0.2091, 0.2200, 0.2292, 0.2555 → **median 0.2020** against a published
median of **0.0888**, ratio **2.27x**. I would have written the max as if it were typical.

⭐ **Getting (1) right makes the paragraph internally consistent rather than accidentally
contradictory**: residualization harms transfer (0.231) MORE than native (0.155), which is the same
direction as the −0.119 paired delta. The marginal phrasing had hidden that agreement.

⚠ A third level of the same issue: the deposited `delta_median` is the median over 200 draws of the
per-draw difference and diverges from (transfer_median − native_median) **per cell** by up to 0.0165.
**Paired statistics never reduce to differences of marginals here.**

### ⛔ The pooled-residualization deposit note is RETRACTED — V9 refutes it

**Replacement text, verbatim:** *The technical residualization was fitted across all 177 participants
at once. Whether this removes genuine between-etiology signal was tested directly by refitting it
within each etiology (group means added back) and repeating the decisive arm on identical draws. The
two specifications disagree by more than the effect under study: the paired delta is −0.119 (seven of
eight arms negative) fitted pooled and +0.128 (zero of eight negative) fitted within etiology, a median
paired shift of 0.253 present in every cell. The sign of the conclusion therefore follows the pooling
choice rather than the data, and the residualized re-test is not a usable control; it is excluded from
adjudicating the etiology comparison, and the unresidualized result stands on its own. A marginal
check — no technical feature differing between etiologies after BH-FDR correction, etiology predicted
from the ten distributional features at cross-validated AUROC 0.558 — did not anticipate this, and a
marginal association should not be used to bound an analysis-specification effect.*

✅ **Control that makes the comparison clean**: the unresidualized variant reproduces V7 **exactly,
max |difference| 0.000000** across all eight cells, so the machinery is sound.
⛔ Residualization also reverses **which arm it hurts**: pooled costs transfer 33.9% and native 25.5%;
within-etiology costs transfer 16.1% and native 33.5%.

### ⛔ CITATION: the deposit is **Yang**, not Zeybel

Zeybel M is the **last** author. Every other cohort in this project is mapped by *first* author
(Govaere, Pantano, Suppli), so "Zeybel 2025" is inconsistent and wrong in a manuscript.

> Yang H, Atak D, Yuan M, et al. Integrative proteo-transcriptomic characterization of advanced
> fibrosis in chronic liver disease across etiologies. Cell Rep Med. 2025;6(2):101935.
> doi:10.1016/j.xcrm.2025.101935. PMID 39889710; PMCID PMC11866494.

They already report the direction: *"The transcriptome profiles did not exhibit clear distinctions
among these etiology-based groups."* **The transfer DESIGN is ours; the direction is theirs.**

---

## ✅ RESOLVED: prospective standardization is FREE on the proper axis — my alarm was an artifact

⛔ **This section's original headline was "THE REPRODUCIBILITY CLAIM DOES NOT SURVIVE PROSPECTIVE
STANDARDIZATION". That was wrong and is retracted.** The collapse was a numerical artifact of running
on the unfiltered 42,163-gene axis, not a property of the instrument. Kept below because the
find→diagnose→fix sequence is the useful part.

### The resolution, from job 21328574 on the detection-filtered release axis (26,629 genes)

| arm | GSE268273 ρ | GSE193066 ICC | ΔICC vs histology |
|---|---:|---:|---|
| transductive, latent | 0.6553 | 0.7391 | +0.5505 [+0.277, +0.816] |
| **prospective, latent** | **0.6543** | **0.7557** | **+0.5671 [+0.295, +0.832]** |
| hybrid, latent (**selected**) | 0.6601 | 0.7521 | +0.5635 [+0.289, +0.827] |

**Paired delta prospective − transductive: ICC +0.0158 [−0.029, +0.061], OOD −0.0009 [−0.028, +0.027]
— both cover zero.** Prospective standardization costs nothing. **The tool is shippable, retains
ICC ≈ 0.75, and a stranger with one sample can run it.**

### ⭐ The mechanism, and it is mundane

`results/mechanism.json`. Gene `ENSG00000235699` has **train_sd_log2cpm = 3.3e-07** — silent in the
training pool, nonzero in GSE193066. Dividing by that SD gives **z = 1,459,136**, and that one gene
contributes **+509.9** to a score whose GSE268273 range is [0.87, 2.83].

| | GSE193066 | GSE268273 | training pool |
|---|---:|---:|---:|
| max abs standardized value | **1.46e6** | 299 | 13.0 |
| score range | [0.024, **511.7**] | [0.87, 2.83] | — |

4,217 genes sit below the 10th-percentile training SD. They carry **8.3% of the weight mass** but
**24.5% of the contribution mass in GSE193066** against 3.1% in GSE268273. Every one of the top
blow-up genes is flagged `on_release_axis: false` — **the release axis already excluded them.**

✅ **The detection filter beats SD-flooring**, which was the obvious alternative: floors at p1/p5/p10
give ICC 0.505 / 0.509 / 0.569, all *worse* than filtering's 0.7557. Flooring retains noisy genes with
a rigged denominator; filtering removes genes that should never have been on the axis.

⛔ **My "range" hypothesis was also wrong.** I proposed that GSE268273 survived because it spans F0–F4
while GSE193066 is narrow. The real reason is that GSE193066 simply had more genes that are silent in
the training pool. Nothing to do with stage range.

✅ **Measured axis-coverage threshold 0.30** (conditions hold contiguously to 0.20; next point up for
margin). Selection was made from the training pool alone, before outcomes.

---

## [SUPERSEDED — see resolution above] The original alarm, 2026-09-01

`executions/release-instrument-20260901T211017Z/`, job 21328299 (FAILED at 00:02:53 in the dual→primal
conversion, but the measurement completed first and is in `logs/jobA-21328299.out`).

The published instrument standardizes features **within cohort** (`fit_latent_severity.py:65`,
`zscore_within(X, groups)`). That is **transductive** — it uses the target batch's own per-gene mean
and SD. A user with one biopsy cannot do it. Replacing it with a prospective standardization (per-gene
constants frozen from the 521-sample pool) gives:

| arm | test–retest ICC, 58 pairs | ΔICC vs histology |
|---|---:|---|
| **transductive (as published)** | **0.7208** / 0.6818 | +0.532 [+0.256, +0.799] |
| **prospective (frozen pool moments)** | **0.0016** / 0.0009 | −0.187 |
| **hybrid** | **0.0040** / 0.0027 | −0.185 |

Paired delta prospective − transductive: **−0.58, CI [−0.835, −0.191], excludes zero in every cell.**

⛔ **Meanwhile the out-of-cohort Spearman on GSE268273 barely moves**: −0.017 to −0.033, every CI
covering zero. **The two metrics disagree completely about whether the swap matters.**

⭐ **This is not a leakage artifact — and that makes it worse, not better.** GSE193066 contributes
**0 samples** to the 521-sample training pool (verified: GSE135251 172, GSE162694 112, GSE174478 93,
GSE130970 78, GSE240729 66). So scoring GSE193066 prospectively **is** the realistic simulation of a
stranger running the tool, and under that condition the test–retest ICC is approximately zero.

⛔ **"Reproduces about four times better than the histology label" is a property of the standardization,
not of the instrument**, and it cannot ship in its current form.

⛔ **The agent's release-selection rule chose on OOD Spearman** (prospective 0.632 vs hybrid 0.656) — a
metric nearly invariant to the thing under test — and therefore selected an arm with ICC 0.0040.
**A selection rule that cannot see the failure mode is the wrong rule.**

### Working hypothesis, under test — not adopted

Spearman within a cohort is invariant to any monotone transform applied to all samples alike, and a
*global affine* change leaves ICC unchanged too, so the collapse cannot be a global rescaling. It has
to be **per-gene**: transductive divides each gene by its GSE193066 SD, prospective by the pool SD, and
that reweights the score rather than transforming it. If so, the reason GSE268273 survives is **range**
— it spans F0–F4 so severity still dominates the ranking, whereas GSE193066 is selected high-risk,
stages 1–3, stage variance 0.74x the full 106, narrow enough for admitted technical variance to swamp
the between-participant spread.

Two diagnostics requested: (i) correlate transductive vs prospective scores within GSE193066 — near
-identical refutes the hypothesis; (ii) decompose the prospective score's variance into
between-participant and within-pair terms, since ICC≈0 can mean the between term collapsed **or** the
within term exploded, and those have different fixes.

### The untried construction that could rescue it

**Per-sample quantile normalization to a frozen reference** derived once from the 521 pool: prospective
by construction, needs nothing but the single sample being scored, and removes exactly the per-gene
scale sensitivity implicated above. Now added as a fourth arm.

**If nothing prospective recovers the ICC**: the tool ships as a severity *ranking* instrument for
cohorts with a full stage range, and the model card states that test–retest reliability under the
shipped standardization was measured at approximately zero. Not "reduced". Not "attenuated".

---

## ⛔⛔⛔ THE KLEINER CONSTANT DOES NOT CLEANLY REPLICATE — and the ICC framing is unsafe — 2026-09-01

`executions/kleiner-reliability-replication-pxd051911-20260901T211358Z/`. Independent deposit
PXD051911, 57 patients with ≥2 Kleiner grades, split into four arms. Guards **12 passed, 0 failed**
(10 demonstrably_live, 1 verified_positively, 1 live_only_under_a_changed_constant).
Pairing verified by positive control: **age ICC 0.9981, sex concordance 27/27**.

| arm | n | σ_e | CI95 | ICC(1,1) | stage var | verdict |
|---|---:|---:|---|---:|---:|---|
| **A1** initial→at-surgery | 27 | **0.5092** | [0.304, 0.694] | 0.6228 | 0.687 | **`consistent_but_underpowered`** |
| **A2** initial→follow-up, no surgery | 30 | **0.4655** | [0.365, 0.548] | 0.6117 | 0.558 | ⛔ **`discrepant`** |
| A3 post-surgical | 11 | 0.4264 | [0.213, 0.564] | 0.5918 | 0.446 | — |
| A4 pooled mixture | 41 | 0.4553 | [0.366, 0.530] | 0.6119 | 0.534 | — |
| **source GSE193066** | 58 | **0.6433** | [0.517, 0.760] | **0.1886** | ~0.51 | — |

⛔ **A1 is a failure to reject, not a replication** — the artifact says so itself. The source point sits
inside A1's CI, but A1's point does *not* sit inside the source CI, and the ratio CI [0.437, 1.156]
falls outside the pre-declared equivalence margin [0.75, 1.333]. **MDE at n=27: only a ratio ≥1.58x or
≤0.61x is detectable, so a genuine 25% difference is invisible.**
⛔ **A2 is formally `discrepant`**: ratio CI [0.5366, 0.9555] **excludes 1**.

### ⭐⭐⭐ The result that matters: ICC moved 3.3x, σ_e moved 1.38x

**The Kleiner label's own ICC is 0.19 in GSE193066 and ~0.62 in PXD051911.** Same label, same scale,
different stage variance — exactly the dependence σ_e was chosen to avoid. **This is the σ_e-over-ICC
argument demonstrated empirically instead of asserted, and it is now the strongest contribution in
the lane.**

⛔ **Consequence for the flagship, and it is severe.** Any sentence of the form *"Kleiner staging has
ICC 0.19"* or *"the instrument reproduces better than the histology label"* is **contradicted by an
independent deposit**. Safe only as a **within-sample contrast at a stated stage variance**.
Combined with the release finding (model ICC 0.7208 transductive → **0.0016 prospective**), the
"reproduces ~4x better than the label" claim needs **both** the transductive standardization **and**
GSE193066's particular stage variance. **Neither transports.**

### ✅ Two competing explanations for the σ_e gap were tested and BOTH fail

1. **Stage composition does not explain it.** Standardising A2 to the source's stage mix moves σ_e to
   **0.2232**, *away* from 0.6433, widening the gap rather than closing it.
2. **Real change does not explain it.** A2 shows significant net progression (11 up / 17 same / 2 down,
   sign p = 0.0225), and MSW estimates σ²_e + σ²_D/2, so real change can only **inflate** it — yet A2
   is already *below* source. Removing the net shift lowers it further, to 0.4143.

⚠ **My brief to the agent had the arm structure wrong.** I assumed follow-up meant post-bariatric. The
agent correctly separated A2 (no surgery sample; net progression) from A3 (post-surgical; 1 up / 3
down, improvement). The split is the agent's, not mine.

⚠ **A1 vs A2 as primary is genuinely unclear**: A1's CI is 76.5% of its point estimate, A2's is 39.2%,
and A2 is the only arm returning a verdict. A1 has the cleaner inter-read interval; A2 has precision.
Both reported; recommendation requested.

### What the shipped constant should be — my own analysis, pending the producer's recommendation

⛔ **First, an independence trap I walked into and had to back out of.** The four PXD arms are **not**
independent studies: A4 (41) = A2 ∪ A3, and A3 (11) is a **subset of A1's patients** (the 11 with all
three sample types). Only **A1 (27) and A2 (30) are disjoint**. Pooling all four inflates the weight of
the PXD deposit and manufactures precision. Any meta-analysis here uses A1 and A2 only.

On the three disjoint arms: fixed-effect σ_e **0.5472**, random-effects **0.5403** CI95 [0.428, 0.682],
Q = 5.27 on 2 df (p = 0.072), I² **62.0%**.

⛔ **But that random-effects number should not ship either.** The heterogeneity is *between deposits*
(A1/A2 ratio is only 1.094, while the deposit contrast is ~1.3x), and with **n_deposits = 2 a
between-deposit variance cannot be estimated**. A τ² fitted on two deposits dressed as three arms is a
number with no sampling justification behind it.

✅ **Recommended form: a RANGE with the mechanism named, not a pooled point.**

| quantity | GSE193066 | PXD051911 |
|---|---:|---:|
| σ_e | 0.6433 [0.517, 0.760] | 0.4655–0.5092 |
| error variance | 0.4138 | 0.2167–0.2593 |
| **ρ at Var(stages) = 1.0** | **0.586** | **0.741–0.783** |
| ceiling √ρ | 0.766 | 0.861–0.885 |
| r_latent multiplier | 1.306 | 1.130–1.161 |

⭐ **Even as a range this beats what it replaces.** `docs/RESULTS.md` previously assumed ρ **0.6–0.9**
and said outright that this was a guess. The measured two-deposit range at Var = 1.0 is **ρ 0.586–0.783**
— it excludes the top third of the assumed band and shifts the whole interval down. **A measured range
that narrows and relocates an admitted guess is still a contribution;** it is simply not the single
universal constant the plan hoped for.

✅ **What a user should be told to do**: measure σ_e in their own setting if they have repeat reads,
and otherwise use the range and report which end they assumed. **Never borrow an ICC.**

### ⭐ What survives the two negatives, stated plainly

**Survives intact:**
- **The instrument transfers.** OOD ρ 0.649 vs native 0.609 on GSE268273, and this is the one metric
  that is **robust to the standardization swap** (−0.017 to −0.033, every CI covering zero).
- **σ_e is far more transportable than ICC** — 1.38x versus 3.3x across deposits. This is now measured
  in two deposits rather than argued, and it is the methodological contribution.
- **The disattenuation machinery**, re-parameterised by a range.
- **Lane D**, as a bounded null, benchmark-only.

**Does not survive in general form:**
- "Reproduces ~4x better than the histology label" — needs both the transductive standardization and
  GSE193066's stage variance.
- σ_e = 0.6433 as *the* constant.
- The reproducibility gain as a property of the shipped tool — pending the per-sample quantile arm.

**Survives as a within-sample result, with conditions stated:** the 58-pair headline that the score
predicts which biopsy is wrong. It must carry the standardization and the stage variance.

### ⛔⛔ WHICH END OF THE σ_e RANGE IS CONSERVATIVE — I had this backwards, and so does the project

From `posthoc_pooling_options.json`. **A larger σ_e gives a SMALLER ρ, hence a LOWER ceiling √ρ and a
LARGER disattenuated r_latent = r_obs/√ρ. Both of those FLATTER a predictor.**

| claim shape | which σ_e end favours it |
|---|---|
| "our model approaches the measurement ceiling" | **large σ_e (0.6433)** — lowers the ceiling toward the model |
| "the true correlation is higher than it looks" | **large σ_e (0.6433)** — inflates r_latent |
| conservative choice for both | **small σ_e (0.4655)** |

⛔ **The project has been using 0.6433 — the anti-conservative end — for exactly those two claim
shapes.** Every ceiling-proximity and disattenuation statement made so far is flattered by the choice
of deposit, and none of them said so.

⛔ **The worked example I recorded is fragile and must be dropped or rewritten.** "A reported r = 0.5
implies r_latent = 1.15, which is impossible" holds only at the large end: at stage SD 0.7 the
disattenuated value is **1.2679 at σ_e 0.6433 but 0.6695 at σ_e 0.4655** — entirely possible. The
"impossible" punchline was an artifact of picking the anti-conservative deposit.

✅ **Every disattenuation statement must now name the σ_e it used and the direction of the bias.**

### ✅ The producer independently reached my pooling conclusion, with better reasons

It rejected meta-analytic pooling for two reasons stronger than mine:
1. **DerSimonian–Laird τ² is severely downward biased and unstable at k = 2–3**, so the random-effects
   CI is untrustworthy however it lands (its own runs: k=2 gives I² 75.8%, RE σ_e 0.554 [0.404, 0.760];
   k=3 gives I² 58.1%, RE 0.5435 [0.442, 0.668]).
2. **Random effects assumes the deposits are exchangeable draws from a population of settings** — denied
   here by a composition difference at **Fisher p = 7.5e-07** (F0/F1 vs F2+) and by the unverifiable
   reading-independence question.

It also confirmed my A4 ICC arithmetic independently (0.6119 quoted, 0.6119 computed).

---

## ⛔ The chromatin residual carries NO biology beyond technical variation — outcome B — 2026-09-01

`executions/chromatin-residual-severity-20260901T211335Z/`. Job 21328587, COMPLETED 00:02:11.
GSE267145, **99 participants**. Fibrosis tie ceiling **0.791771**, reproduced three independent ways
and matching the prose value to 6 dp — so `0.7918` finally has a producing file.

| statistic | value | CI95 |
|---|---:|---|
| **D1** = ρ(RNA+residual) − ρ(RNA) | **+0.1656** | [+0.0177, +0.3342] |
| **D2** = ρ(RNA+tech+residual) − ρ(RNA+tech) | **−0.0729** | [−0.1591, +0.0068] |
| MDE(D1) at 80% power | 0.2248 | — |

⛔ **D2 is the decisive comparison and it is null.** Once technical summaries are in the baseline, the
chromatin residual adds **nothing** — the point estimate is negative and the interval covers zero.

⛔ **What the residual actually predicts is library complexity**, not chromatin biology. Of 21
technical features, 4 pass BH q<0.05 against the residual-only model prediction, and they are all
H3K27ac distributional descriptors: `shannon_entropy` ρ **+0.668**, `gini` **−0.659**,
`frac_counts_top1000` **−0.633**, `frac_counts_top100` **−0.576**. Agreement between the technical-only
and residual-only predictions is ρ **+0.513**.

⚠ **D1 excludes zero but sits BELOW its own MDE (0.166 vs 0.225)** — a detection at roughly 50% power,
so the effect size is subject to winner's curse and should not be quoted as a magnitude.
⚠ **It is rank-fragile**: D1 = +0.180 at K=5, +0.166 at K=10, **+0.003 at K=20**.

### ✅ The producer disclosed a post-hoc gate amendment — and the disclosure is validated by the outcome

Its gate threshold **failed**, and it wrote the amendment **after seeing it fail**, recording that as a
defect in its own prespecification: the constants 0.02/0.10 were borrowed from a lane whose statistic
is a *bounded* profile-deviance skill that shrinks toward baseline, whereas D1 is a *difference of two
Spearman correlations* on ~78 training rows, whose null legitimately has variance. A threshold on the
null's upper percentile therefore demands a **degenerate** null and conflates "the null is centred
wrong" (a defect) with "the null has spread" (necessary). The sealed prespec fixed the gate's **form**
(one-sided on positive increment) but not its constants.

⭐ **The amendment is credible precisely because it opened the door onto a NEGATIVE.** Motivated
loosening produces a positive; this one let through outcome B. No reported number moves either way —
D1, D2, their intervals and the permutation p are computed against the same full null regardless; the
gate only decides whether the primary opens.

### ⭐⭐ THE CONVERGENT PATTERN — three independent lanes *(⚠ later SEVEN — see the modality push)*, same confound, one session

1. **GSE276114 etiology**: ten distributional library summaries recover much of the expression model,
   and residualizing them out removes a quarter of native and a third of transfer performance.
2. **The release standardization**: swapping per-gene scaling from cohort-specific to pool-frozen
   collapses test–retest ICC from 0.7208 to 0.0016 — admitting technical variance swamps a narrow
   between-participant spread.
3. **The chromatin residual** (here): the apparent severity increment is H3K27ac library complexity,
   and it disappears once technical features are in the baseline.

⛔ **Library-complexity descriptors — entropy, Gini, top-N count shares — are the recurring confound in
this project, not albumin or mitochondrial fraction.** Every lane that finds an apparent signal must
now include them in the baseline before claiming biology. Note they are partly biological (a cirrhotic
liver genuinely has a more skewed transcriptome), so they **bound** the technical contribution from
above and never estimate it.

---

## ⭐⭐⭐ SHIPPED: `release/masld-severity-v1/` — a stranger can run it — 2026-09-01

`executions/release-instrument-20260901T211017Z/`. Jobs 21328573 (6:52, 3.4 GB) and 21328574 (5:07,
2.1 GB), 12 CPU / 180 GB / partition `cpu`. Prespec `4ae1562b…`, amendments `819f00ed…` / `cf702191…`,
all verified at runtime. **Guards 9/9 passed, all demonstrably_live, 0 failed.**

```
python score.py --counts my_counts.tsv --out scores.tsv
bash test_release.sh          # asserts OOD Spearman 0.660132 end to end
```

⭐ **The smoke test re-derives 0.660132 on 109 PARTICIPANTS through `score.py` from raw files** — the
released code path, not the training path. G1 reproduces the published `latent_42k` 0.649204 to
**1.5e-7**, proving the harness was scoring the published instrument.

### Input contract
**Unversioned ENSG, 26,629 genes, fixed order.** Trailing `.7` stripped; duplicate stable ids
**summed** (the PAR lesson — P2RY8 / CSF2RA / ASMTL-AS1 would double-count and corrupt the library
denominator). **Counts only.** The scorer raises on CPM/TPM (identical column sums), on logged input
(max < 25, since log2CPM ≤ 19.9), on a zero join, on coverage below the measured **0.30** threshold,
and on four HCC accessions held in a live denylist. ✅ **Each refusal was verified to have ACCEPTED
before it was fixed** — the first CPM check was a library-size heuristic that let CPM straight through.

### Performance, with the honest phrasing
| endpoint | value |
|---|---|
| OOD Spearman, 109 participants | **0.660132** = 69.3% of the 0.952322 ceiling |
| vs native 0.609434 | paired delta **+0.0509 CI [−0.028, +0.137] — INCLUDES ZERO**, MDE 0.118 |
| test–retest ICC, 58 participants | 0.7521 |

⛔ **"Not worse than native", NOT "beats native".** This retires the earlier "transfers ABOVE native"
framing everywhere it appears.
⛔ **The prospective-vs-transductive claim is bounded, not null**: MDE 0.039–0.041 Spearman and
0.064–0.081 ICC, so the correct statement is *any cost is below 0.04 Spearman / 0.07 ICC*, which is all
this design resolves.

### ⛔ The intervals are NOT calibrated — disclosed unprompted
Raw coverage on GSE268273: **0.330** at nominal 50%, 0.752 at 80%, 0.927 at 90%, 0.982 at 95%.
Under-covers in the middle, over-covers in the tails, and the score carries an uncorrected per-cohort
offset a caller cannot compute. **This retreats from the earlier internal claim of 0.478 at nominal
50%.** ⭐ **The ordering is the deliverable; the absolute value is not a Kleiner stage.**

### ⛔ Cancer exposure — my brief was WRONG and the correction is worse
I wrote "1 of 9 training cohorts is affirmatively cancer-free." The nine-cohort list is the project's
**evaluation** substrate. The model trains on **five** (GSE130970, GSE135251, GSE162694, GSE174478,
GSE240729) and **all five are CANNOT_DETERMINE**. GSE193066 is the only affirmatively clean cohort and
it is the **validation** set, excluded by a leakage guard — independently confirmed, 0 rows in
`pool_samples.tsv`. **0 of 5 training cohorts are affirmatively cancer-free.**

### Amendment 02 is outcome-informed and says so on its face
The decision to change the axis followed the observed collapse. The remedy uses only the inherited
training-pool detection rule (≥1 count in ≥20% of samples in every training cohort, computed on the
521 training samples); **no validation outcome enters it.** Correct handling of an unavoidable
post-hoc decision.

⚠ Primal weights reproduce the dual to **1.33e-14** — this certifies **numerical safety and that the
deposited vector is the one in use**, NOT two independent estimators agreeing. The card says so.
⚠ The k5 panel is a demonstration arm and is **not in the release tree**. No panel size was selected;
the card records the k-curve's non-monotonicity.

### ✅ INDEPENDENTLY VERIFIED: the release runs, checked by an agent that did not build it

Job 21328683, `executions/independent-release-verification/`. `test_release.sh` run from a clean shell
on a compute node:

```
scored 109 samples
axis coverage 1.0000 (26629/26629); threshold applied 0.30 (measured 0.30)
stripped a version suffix from 42163 identifiers
== ALL CHECKS PASSED          SMOKE_TEST_EXIT=0
```

It ingests versioned gene ids as a stranger's file would, strips them, applies the coverage gate, and
reproduces the expected OOD Spearman through `score.py` alone. It also emits an unprompted scope
warning on stderr, so a caller who never opens the model card is still told the absolute level carries
an uncorrected per-cohort offset. ⭐ **The "a stranger can run it" claim is verified, not accepted.**

⚠ **A trap I walked into with a memory already on file.** The first attempt (21328621) reported
`COMPLETED` exit 0 and produced **nothing readable** — `--output` pointed at `/scratch`, which is
**node-local**, so the log was written to the compute node's own disk. A job that reports success and
leaves no evidence is indistinguishable from one that did nothing. Always give SLURM an absolute
`/gpfs` output path.

### ⛔ I OVERSTATED the `kleiner.py` defect — corrected after running it instead of reading it

I claimed a user in a low-variance population would be told a legitimate result is impossible. **Wrong.**
The demo explicitly scopes itself to "a population like the 58-participant paired set, whose Kleiner
stage variance is 0.5100", names the escapes (more stage variance, an optimistic r, averaged reads),
and closes on the correct lesson that r = 0.50 is uninterpretable until the stage variance is on the
page. **The variance handling is right and that criticism is withdrawn.**

**The real defect is narrower: the punchline flips on an undisclosed choice of deposit.**

| σ_e | ρ at var 0.51 | ceiling | r_latent from r = 0.50 |
|---|---:|---:|---:|
| **0.6433** GSE193066 | 0.1886 | 0.4343 | **1.1513 — "IMPOSSIBLE"** |
| **0.4655** PXD051911 | 0.5751 | 0.7583 | **0.6594 — unremarkable** |

Same population, same reported correlation, opposite verdicts. ✅ **Fix is to print BOTH rows, not to
delete the example** — the contrast is a sharper teaching point than the original and is the honest
state of the evidence.

### ⭐⭐ Job C resolves the standardization question — per-sample quantile wins on every metric

Job 21328736 (11:45). Added the fourth arm and the variance decomposition.

| arm, release axis 26,629 | in-pool LOCO | OOD | **ICC** | s²_between | s²_within |
|---|---:|---:|---:|---:|---:|
| hybrid, latent | 0.6637 | 0.6601 | 0.7521 | 0.1482 | 0.0488 |
| prospective, latent | 0.6840 | 0.6543 | 0.7557 | 0.2488 | 0.0804 |
| **quantile, latent (CHOICE)** | **0.6872** | **0.6631** | **0.7658** | 0.2458 | 0.0752 |

**Hybrid was dropped from eligibility.** Per-sample quantile normalization to a frozen reference is
best on in-pool LOCO, OOD *and* ICC, and it is the only construction that is fully prospective with no
reference to a training distribution beyond fixed quantiles — a caller with **one sample** can run it.

### ⭐ The variance decomposition answers the mechanism question exactly

| arm | s²_between | s²_within | ICC |
|---|---:|---:|---:|
| full 42,163, prospective | 22.23 | **14,219.72** | 0.0016 |
| full 42,163, hybrid | 8.62 | **2,235.47** | 0.0038 |
| full 42,163, **quantile** | 0.3554 | **0.3688** | **0.4908** |
| release 26,629, quantile | 0.2458 | 0.0752 | 0.7658 |

⭐ **The within-pair term EXPLODED by five orders of magnitude; the between-participant term did not
collapse.** That was the diagnostic question — ICC≈0 can mean either, and they have different fixes.
It was the within term, so the fix is to stop admitting the offending genes, which is what the axis
filter does.

⭐ **Quantile normalization fixes it independently of the axis filter** (ICC 0.4908 on the *unfiltered*
axis, against 0.0016 prospective and 0.0038 hybrid) because it never divides by a per-gene SD at all.
Belt and braces: the shipped arm has both.

### ⭐ Why OOD Spearman never noticed the failure

Score correlation, transductive vs prospective:

| axis | Pearson | Spearman |
|---|---:|---:|
| full 42,163 | **+0.0925** | **+0.8492** |
| release 26,629 | +0.9844 | +0.9787 |

**On the broken axis the scores are nearly uncorrelated in Pearson terms but still rank-preserved.**
A handful of blown-up samples destroy the linear relationship while leaving most of the ordering
intact — so a rank metric survives while a variance-ratio metric dies. ⛔ **This is why OOD Spearman
was the wrong selection criterion, and why an arm was nearly shipped with ICC 0.004.**

---

## ⛔⛔ V9: the residualized delta REVERSES SIGN with the residualization scope — 2026-09-01

Job 21328403 (01:05:51). `v9_results.json`.

| variant | median transfer | median native | **median delta** | half-width | cells excluding 0 |
|---|---:|---:|---:|---:|---:|
| none, published | 0.6795 | 0.6620 | **+0.0012** | 0.089 | 0 |
| pooled, as in V7 | 0.4512 | 0.4549 | **−0.1188** | 0.202 | 1 |
| **within-etiology** | 0.5725 | 0.4052 | **+0.1279** | 0.145 | **3** |

**A 0.253 swing and a sign change, produced by nothing but the scope of the residualization.**
⚠ **0.253 is the median of the per-cell PAIRED differences; 0.247 is the difference of the two median
deltas.** I first quoted 0.247 — violating, within the hour, the rule I had just recorded that paired
statistics never reduce to differences of marginals here.

### ⛔ The conservatism argument was BACKWARDS and I propagated it

The verifier argued — and I recorded — that pooled residualization biases toward finding NO loss while
V7 moved toward MORE loss, so the residualized result is conservative. **The measurement says the
opposite: pooled PRODUCES the apparent loss that within-etiology residualization removes.** Pooled
−0.119, within +0.128. Retracted.
⚠ Two agents agreed on an unmeasured direction, which made it *less* likely to be checked, not more.

### ⭐⭐ The methodological finding — worth more than the etiology answer

The deposit reads: CV AUROC predicting etiology from the ten distributional features is **0.5581**, and
*"near 0.5 means the pooled/within distinction is second order."* **Its own experiment refutes that
bound — a 0.558 classifier produced a 0.253 swing in the delta.**

⛔ **A weak classification AUROC does NOT bound the regression-relevant variance removed by
residualizing on those features.** Classification asks whether a direction separates two groups;
residualization deletes the entire projection onto it. A direction that barely discriminates can still
carry a large share of the variance a regression uses. **Never use a classification AUROC as a bound
on a residualization's impact.**

### ✅ Why this strengthens the lane rather than weakening it

Neither residualized arm is identified. Within-etiology residualization removes exactly the
within-etiology variation a native model exploits — which is why native falls furthest there (0.405 vs
pooled 0.455 vs unresidualized 0.662) and why transfer then appears to beat it. Pooled removes shared
structure and tilts the other way.

⭐ **The sign of the residualized delta is a function of an analysis choice with no principled
resolution, so the residualized re-test is UNIDENTIFIED and cannot adjudicate the etiology question in
either direction.** That is stronger than "too imprecise" and replaces that clause in sentence (3).
**The published unresidualized result — delta +0.0012, 0 of 8 cells excluding zero — stands untouched.**

---

## ⛔ LANE C — the non-invasive route does NOT work, and the gap is large — 2026-09-01

`executions/pxd051911-plasma-severity-20260901T191009Z/`. Job 21327263, COMPLETED 02:55:30.
**Guards 13/13 met (11/11 demonstrably_live, 3 verified_positively, 2 live_only_under_a_changed_constant).**
G13 and G14 both **fired their violations**, so the verdict rule is demonstrably able to return every
one of its categories.

### Arm 2 — routing through predicted liver state does not help (n = 48)

| model | Spearman |
|---|---:|
| plasma → severity, direct | **0.0497** |
| plasma → predicted liver → severity, **routed** | **−0.1306** |
| **liver expression → severity, direct** | **0.6989** |

Paired delta routed − direct **−0.1804** [−0.394, +0.024], MDE 0.2982 → `UNDERPOWERED_untestable`.

⛔ **Routing is worse than direct, not better.** The composed model
*plasma → predicted liver state → latent severity* was the plan's highest-upside idea and it does not
survive contact with the data.
⭐ **The number that matters is the third row.** On the same 48 participants, liver expression predicts
severity at **0.699** while plasma protein predicts it at **0.050**. The signal is in the tissue; the
plasma proxy is not carrying it.
⚠ Only **7 of 48** have plasma at visit 2 and liver at visit 1, so the routing arm is thin by
construction. And note plasma reached 0.2542 against SAF at n=143 in Arm 1 — **0.0497 at n=48 is a
small-sample figure and should not be quoted as plasma's ceiling.**

### Arm 3b — the plasma score does NOT predict which biopsy is wrong

β **−0.0431** [−0.236, +0.062], permutation p **0.667**, ΔR² +0.0047, MDE 0.2131 →
`UNDERPOWERED_untestable`. The headline result that holds for the liver instrument does **not** transfer
to plasma. Guard `G3b_visit1_only` verified positively.

### Arm 3a — a third σ_e figure, and it is NOT independent

Kleiner on the 41 paired visits: exact 0.585, within-1 **1.000**, mean delta +0.171,
**σ_e 0.4445, ICC 0.6119**, sign p 0.143.
⚠ **This is the same substrate as the replication lane's A4** (41 pairs, σ_e 0.4553, ICC 0.6119) — not
independent corroboration. ✅ **But two agents with separate implementations returned ICC 0.6119
identically**, which is a reproducibility check worth having.
⛔ NAS is far noisier: **σ_e 1.2067**, ICC 0.4880 — reinforcing that NAS stays descriptive.
⚠ C3 BMI coupling +0.237 (Kleiner), +0.365 (NAS): the retest signal is partly tracking weight change.

### What Lane C reduces to
**The plasma proteome carries severity information (Arm 1, p = 0.0091) but cannot be shown to add over
three free clinical covariates, cannot be usefully routed through predicted liver state, and does not
inherit the liver instrument's forward-validity result.** Every arm that could have been positive is
`UNDERPOWERED_untestable` at n = 48–143. **This is a bounded negative, not a demonstration of absence.**

---

## ⭐ Encoder benchmark release — build verified, two findings worth keeping — 2026-09-01

`executions/encoder-benchmark-release-20260901T221202Z/`. Jobs 21328838 (baselines, 3:31) and 21328843
(build+verify, 1:59), both COMPLETED exit 0.

### ⭐⭐ "Cells are not biological replicates" — now a measured factor, not a principle

Guard `G08_DONOR_IS_THE_UNIT`, probe `hvg_pca_task_native::two_layer_mlp`, 2,000 resamples:

| interval | width |
|---|---:|
| donor-clustered (102 donors) | **0.022396** |
| cell-weighted (50,000 cells) | **0.001061** |
| **ratio** | **21.11x** |

⛔ **Bootstrapping 50,000 cells instead of 102 donors narrows the interval by a factor of 21.**
The project has carried "cells are not biological replicates" as a rule; this is the first time it has
been *quantified* on this substrate. Any cell-level interval on this data overstates precision ~21-fold.

### ✅ The exposure problem is enforced structurally, not documented

Guard `G11_NO_HARMFUL_NUMBER_SHIPPED`: **4 model/head combinations have zero clean held-out studies** —
`geneformer_v2_316m` (both heads) and `transcriptformer_tf_sapiens` (both heads). Every `.tsv`, `.json`
and `.md` in the release was scanned: **0 offending lines, 0 offenders.** A model with no clean study
ships **no number at all**, which is what makes the leaderboard readable rather than misleading.

### Verdicts stable across 8 seeds, zero flips

| model | linear head | MLP head |
|---|---|---|
| `hvg_pca_task_native` (the 50-d PCA baseline) | below reference | **crosses zero** |
| `scimilarity_v1_1` | above reference | above reference |
| `uce_33l` | crosses zero | **below reference** |
| `uce_4l` | **above reference** | **below reference** |

⭐ **The head changes the verdict, and in opposite directions for the same encoder.** `uce_4l` is above
reference with a linear head and below it with an MLP; `uce_33l` crosses zero with linear and is below
with MLP. ⛔ **A leaderboard that fixes one head would report a different winner than one that fixes the
other** — which is the reusable content of this benchmark, and the reason the harness ships both.
⚠ `uce_4l` (4 layers) beats `uce_33l` (33 layers) on the linear head — depth is not buying anything here.

---

## ⭐ THE NESTED INCREMENT IS POSITIVE — plasma adds over clinical covariates — 2026-09-01

`executions/pxd051911-plasma-severity-20260901T191009Z/results/DELTA_NULL_DIAGNOSTIC.json`.
Job 21328169, COMPLETED 01:13:43. **n = 143 participants.** This is the corrected test: the plasma
block's rows are permuted while y, the clinical covariates and the folds are held fixed, and **both**
arms of a nested pair are refit every draw — `arm A: ridge on [clinical, plasma_perm]`,
`arm B: ridge on [clinical]`, with `lam_p = inf` in the grid so arm B is literally arm A restricted.

| outcome | nested increment | CI95 | null mean ± sd | perm p | registered non-nested delta |
|---|---:|---|---:|---:|---|
| Kleiner fibrosis | +0.1916 | [−0.020, +0.398] | −0.0605 ± 0.1005 | **0.0135** | +0.1971 [−0.033, +0.428] |
| **NAS** | **+0.2791** | **[+0.118, +0.434]** | −0.0705 ± 0.0699 | **0.0005** | +0.2470 [+0.042, +0.444] |
| SAF | +0.1900 | [−0.007, +0.388] | −0.0481 ± 0.0659 | **0.0025** | +0.1363 [−0.100, +0.373] |

⭐⭐ **NAS is a clean positive on both statistics**: the increment excludes zero **and** beats its
permutation null. **The plasma proteome carries NAS information beyond age, sex and BMI at n = 143.**

⭐ **The null is centred where theory says it must be** — **−0.048 to −0.071**, near zero with the
slight negative offset from arm A's extra noise columns. ⛔ The rejected construction (permuting plasma
against an *unpermuted* stored clinical arm) would have centred at **−0.118** and manufactured
significance for all three. Both the mis-specification and its magnitude were predicted before the run.

✅ **Kleiner and SAF behave exactly as PRE-REGISTERED**: small permutation p with a CI covering zero is
**not a contradiction**. The deposit committed to the reading in advance — *the increment exceeds a
noise control within this cohort, and the cohort-level increment is still not estimated precisely
enough to exclude zero.* A permutation conditions on the observed 143 participants; a bootstrap does
not. Having that sentence sealed beforehand is why it is reportable rather than arguable.

⛔ **Scope.** This is Arm 1, cross-sectional, n = 143. It does **not** rescue Arm 2 (routing plasma
through predicted liver state is *worse* than direct: −0.131 vs +0.050, with liver itself at 0.699) or
Arm 3b (plasma does not inherit the forward-validity result: β −0.043, p 0.667). **The defensible claim
is that plasma protein adds to free clinical covariates for NAS — not that it replaces a biopsy.**

---

## ⭐⭐⭐ EXTERNAL VALIDATION: three cohorts, 585 participants, all NOT WORSE THAN NATIVE — 2026-09-01

`executions/external-validation-expansion-20260901T190000Z/`. Job 21329195, COMPLETED exit 0, **1:01**,
MaxRSS 0.62 G. Prespec `b67a1c64…` verified at runtime. **Frozen released weights scored through
`release/masld-severity-v1/score.py` as a caller would — not the fitting code.**

| cohort / arm | n | ρ | ceiling | **ρ/ceiling** | native | Δ vs native | CI95 | perm p |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| GSE268273 (shipped ref) | 109 | 0.6631 | 0.9523 | 0.696 | 0.6094 | **+0.0537** | [−0.025, +0.137] | — |
| **GSE213621 PRIMARY** 3-level | **299** | 0.6382 | 0.9423 | 0.677 | 0.6784 | −0.0402 | [−0.099, +0.017] | 0.0001 |
| GSE213621 sensitivity 4-level | 367 | 0.7473 | 0.9659 | 0.774 | 0.7510 | −0.0037 | [−0.044, +0.039] | 0.0001 |
| **GSE276114 PRIMARY** 3-level | **177** | 0.6739 | 0.8483 | **0.794** | 0.6853 | −0.0114 | [−0.080, +0.054] | 0.0001 |

⭐ **The normalized figure is the headline: ρ/ceiling 0.68–0.79 across three cohorts, two platforms,
three fibrosis coding schemes and axis coverage from 0.63 to 1.00.** That stability across instrument
variation is a stronger generality statement than any single ρ.

⛔ **"Matches an in-cohort baseline across three cohorts" — NOT "beats" one.** All three new deltas are
**negative** (−0.040, −0.004, −0.011) against the shipped cohort's +0.054. **One favourable point
estimate out of four does not license `fraction_of_native > 1` as a general claim**, and the release's
own wording was corrected accordingly.
⛔ **The registered primary is the weakest arm**: Δ −0.040 with a lower bound of **−0.099**, so a loss
of up to ~0.10 is not excluded there. Stated plainly rather than softened.

### ✅ The Control ruling, turned from an argument into a number

Adding the 68 controls moves the **instrument** +0.109 (0.638 → 0.747) and **native** +0.073
(0.678 → 0.751). So the control-vs-disease contrast is easy for both and the delta survives either way
— but the **absolute** figure is inflated by ~0.11. ⛔ **Headlining the 4-level arm would have published
0.75 for what is a 0.64 result on diseased participants.** The harder arm was registered as primary and
the instrument cleared it.
⚠ Dropping Control cost almost nothing in power (MDE 0.1471 → 0.1482), so **there was never a power
argument for keeping it — only a flattering one.**

### Integrity, five guards all live
Prespec pinned and aborting on mismatch; **all 10 released-artifact hashes verified before a single
sample was scored** (closing the receipt-staleness problem); phenotype joins complete and one row per
participant; crosswalk verified 1:1 with **zero axis genes reachable from more than one symbol**;
**denylist live** — handed GSE193084, scorer exited rc=2.
⚠ GSE213621 **pre-merges F0F1 and F3F4**: three effective severity bins from 299 participants, not five.
⚠ MDEs (0.083, 0.059, 0.095) resample test participants against fixed predictions and do not resample
the native model's training, so they are a **lower bound on the true SE** — same estimator as the
release, comparable by construction.

### What this supports, every number measured through the released path
**"Out of cohort the instrument reaches Spearman 0.64–0.75 against tie ceilings of 0.85–0.97 — 68% to
79% of the maximum any score could reach — in three cohorts totalling 585 participants that the model
never trained on, spanning MASLD, chronic viral hepatitis and alcohol-related disease, and in each it
matches a ridge trained inside that cohort."**

---

## ⭐ What the instrument is WORTH, in the unit clinicians use — 2026-09-02

⚠ **POST-HOC and descriptive.** Computed by me on the already-opened external-validation scores in
`executions/external-validation-expansion-20260901T190000Z/`. **Not prespecified, not a registered
result.** Participant bootstrap, 2,000 draws, seed 20260902.

**Question**: rank participants by model score, biopsy from the top down, and stop at 80% sensitivity
for advanced fibrosis. How many biopsies does that take?

| cohort | n | prevalence | fraction biopsied for 80% sens | **biopsies avoided vs random** | CI95 |
|---|---:|---:|---:|---:|---|
| **GSE213621** | 299 | 31.8% (F3F4) | **45.2%** [38.1, 51.8] | **34.8%** | [28.2, 41.9] |
| GSE276114 | 177 | 78.0% (F3/F4) | 66.7% [61.0, 70.6] | 13.3% | [9.4, 19.0] |

Both exclude zero. At the top decile in GSE213621 the PPV is **76.7%** against a 31.8% base rate —
**2.41x enrichment**.

⭐ **Utility is a function of prevalence, and the two cohorts show it cleanly.** In a 32%-prevalence
population — the realistic screening setting — ranking by score avoids about a third of biopsies. In a
78%-prevalence population, already selected for advanced disease, it avoids 13%, because almost
everyone is a case and there is nothing left to enrich. **The instrument is worth most exactly where
biopsy burden is worst.**

### ⛔ The caveat that decides how much this is worth

**"Versus a random ordering" is the wrong comparator and no clinician uses it.** The real comparator is
FIB-4, NFS or FAST — cheap, already standard, and none of them are available for these cohorts. **This
calculation shows the score carries actionable ranking information; it does NOT show it beats existing
non-invasive scores, and that is the comparison a reader will want.** Until it is run against FIB-4 on
a cohort carrying the inputs, this is an upper bound on the incremental value, not an estimate of it.

⚠ GSE213621 pre-merges F0F1 and F3F4, so "advanced" there is a coarse bin.
⚠ 80% sensitivity is a choice; the whole curve is in the producing calculation, not just this point.

---

## ⚠ FLAGSHIP REPLICATION IN LIVER PROTEIN: INCONCLUSIVE, underpowered by construction — 2026-09-02

`executions/pxd051911-liver-latent-vs-label-20260902T002426Z/`. Job 21329333. Cross-modality attempt:
RNA in GSE193066 → **protein** in PXD051911.

✅ **The protein instrument works**: leave-one-participant-out Spearman(score, V1 Kleiner) **0.7291**
over 58 (Pearson 0.7196), against a prespecified floor of 0.20. Scores were **sealed and sha256-hashed
before any V2/V3 Kleiner vector was constructed** (`e3e93c54…`, re-verified at runtime).

| arm | n | **MDE** | β | CI95 | perm p | ΔR² | verdict |
|---|---:|---:|---:|---|---:|---:|---|
| **A at-surgery (~1 y), PRIMARY** | 26 | **0.4198** | +0.2661 | [−0.054, +0.520] | 0.0378 | +0.103 | BOUNDED NULL |
| B follow-up (2–5 y) | 21 | **0.6026** | +0.3527 | [−0.059, +0.784] | 0.0048 | +0.118 | BOUNDED NULL |

⛔⛔ **The MDE EXCEEDS the original effect (0.4174) in both arms.** This test could not have detected an
effect the size of the one it was trying to replicate, let alone half of it. **It is uninformative, not
a refutation.** The artifact says so in its own verdict string rather than leaving it to a reader.

⛔ **My instruction cost the power.** I required the arms be reported separately — correct on the
science, since a follow-up pair contains real change and is a different quantity — but it split n=36
into **26 + 21** and roughly doubled the MDE. A pooled analysis would have had n=36 and MDE ≈0.35, still
marginal but far better. **The right design was pooled with an arm interaction as primary** (the agent
ran it as an H1 secondary: arm-A β +0.3278, interaction −0.0864 over 36 participants).

⚠ **The content-free control matters here and cuts against arm A.** 200 random-protein signatures give
mean +0.028, 2.5–97.5 pct **[−0.310, +0.292]**, max |β| 0.386. **Arm A's +0.2661 sits INSIDE that
range** — indistinguishable from a random signature. Arm B's +0.3527 sits outside its own control's
97.5th percentile (+0.264). So what little signal there is lives in the arm that **contains real
change**, which is the arm least able to speak to reading error.

⚠ Same permutation-vs-bootstrap split as the plasma lane, and the same reading applies: a small
permutation p beside a CI covering zero means the statistic beats a noise control **within this
sample** while the population effect stays unestimated.

**Net: the flagship remains single-cohort, single-modality.** Directionally consistent point estimates
(+0.266, +0.353 against the original +0.417) are worth recording and are not evidence.

---

## 🗺 THE MULTI-MODAL EVIDENCE MAP — what exists, what is being built, what is missing

**Written 2026-09-02 to answer one question: does this project have a multi-modal result, or an RNA
result with modality-flavoured footnotes?** Status as of now.

### ✅ ESTABLISHED — cross-modal PREDICTION, already replicated

| claim | effect | n | replication |
|---|---|---:|---|
| **RNA → H3K27ac** profile, held participant | **+0.13345** skill, p 0.000, 5/5 folds | 99 participants | ⭐ **recurs in a second cohort AND a second assay** |
| **RNA → ATAC** | **+0.12985**, p 0.000, 5/5 folds | 39 donors | GSE296875 |

⭐ **This is the strongest genuinely multi-modal result in the project and it is already done.** RNA
predicts a held-out participant's chromatin state, across two cohorts and two assays.
⛔ Its disease-dependence is **single-cohort** (GSE267145 only; both GSE296875 arms are
`UNDERPOWERED_untestable`) — do not describe that part as replicated.
⛔ Inferred chromatin **is** a deterministic function of RNA — ship for annotation/interpretation
transfer only, never as new signal. The **residual** carries nothing beyond technical variation
(D2 −0.073, CI covers zero).

### ✅ ESTABLISHED — severity instrument, ONE modality

| modality | instrument | external validation |
|---|---|---|
| **bulk RNA** | `release/masld-severity-v1/`, shipped + independently verified | ⭐ **3 cohorts, 585 participants**, ρ/ceiling 0.68–0.79 |

### 🔨 UNDER CONSTRUCTION — the modalities that would make it multi-modal

| modality | substrate | status |
|---|---|---|
| **H3K27ac severity** | GSE267145, **99 participants PAIRED with RNA** | ⚠ job 21329564 **FAILED** at 00:00:28 on a guard; superseded by `chromatin-severity-instrument-20260902T011438Z` — **completed, see the NAS result below** |
| **liver protein severity** | PXD051911, 58, **already at LOO Spearman 0.7291** | being packaged; plus RNA-weights→protein cross-modal transfer test |
| **modality inventory** | everything on disk | eligibility table: endpoint, n participants, ceiling, verdict — published before modelling |

### ⚠ ESTABLISHED but NOT an instrument

- **plasma protein → NAS**, nested increment **+0.2791 [+0.118, +0.434]**, p 0.0005, n=143. Real, and
  currently a finding rather than a packaged tool.
- **Zonation from expression**: accuracy 0.7895 vs chance 0.333, 19 participants, one cohort.

### ⛔ WHAT IS MISSING, stated plainly

1. ⚠ **[RESOLVED LATER IN THIS FILE]** No cross-modal AGREEMENT test had been run *at the time of writing* — do an RNA-derived and a chromatin-derived
   severity score rank the same 99 participants the same way? GSE267145 is the only substrate where
   this is possible and it is now in flight.
2. **Every non-RNA instrument is single-cohort.** Protein n=58 one deposit; chromatin n=99 one deposit.
   ⛔ **Three single-cohort instruments called "multi-modal" is weaker than two honest ones**, and the
   paper must say which modality is *validated* and which are *demonstrations*.
3. **The flagship "predicts which biopsy is wrong" is single-cohort, single-modality**, and the protein
   replication was uninformative (MDE 0.42–0.60 against an effect of 0.417).
4. **No modality other than RNA has been tested out-of-cohort at all.**

### ⛔ METHYLATION (GSE105127) — REJECT for a severity instrument, verified from source

Checked against the GEO SOFT record and the source article, **not** a summary file:
- `!Sample_characteristics_ch1` over all **114** sample records: `tissue: liver` (114),
  `isolation: laser capture microdissection` (114), `hepatic zone` PP/IZ/CV (38 each).
  **114 = 19 participants x 3 zones x 2 assays.** No disease, stage, NAS, fibrosis or BMI field on any
  sample.
- Source article **PMC6175862**: the string **"fibrosis" appears 0 times in the entire paper.** It was
  never scored. NAS appears only as **group medians with ranges** (NC 0, HO 0, STEA 3 [1–3], EARLY
  3 [2–4]) — **no per-participant value is published anywhere.**

⚠ **One correction to my framing**: it is not endpoint-free. There is an ordered **4-level participant
group** (NC 4 / HO 5 / STEA 5 / EARLY NASH 5 = 19) in `!Series_overall_design`. But that is a different
endpoint from fibrosis stage, and the group label is all that exists.

⛔ **The 4-level arm is not worth running and was declined.** n=19, 20,000 permutations, seed 0:
tie ceiling **0.9688** (ties are NOT the constraint), null mean −0.0010 sd 0.2359,
**MDE 0.5879** = 0.607 of ceiling. **The binding constraint is n=19.** An arm would need ρ ≈ 0.59 to
clear — at or above what this project's best *in-cohort* fibrosis instruments reach (0.57–0.69).
⭐ **Combined with the toml's own pre-registered limit** (only zonation cohort, no escalation target, so
a negative is indeterminate and never `tested_negative`), **this arm can return a positive or an
indeterminate but never a usable negative.** An asymmetric test that cannot fail informatively is not
worth compute.

⚠ Housekeeping: `config/datasets/gse105127_zonated_rna_rrbs.toml` still carries
`role/status = "blocked"` and `admission_blocking = true`, which is **stale** — but its `notes` body is
current to 2026-08-30 and already records the completed CpG crosswalk and RSEM quantification. **The
status keys are stale, not the body.**

---

## 🗺 MODALITY ELIGIBILITY TABLE — ONE validated modality, four single-cohort demos — 2026-09-02

`executions/modality-eligibility-scouting-20260901T200000Z/modality_eligibility_table.json`. Criteria
sealed `3d5aaf23…` **before measurement**. Scouting only — no modelling, no compute beyond file reads.

### ✅ ELIGIBLE
| modality | cohort | n participants | endpoint | ceiling | MDE |
|---|---|---:|---|---:|---:|
| **bulk RNA** | 5 train + **3 external** | 521 + **585** | fibrosis | 0.85–0.97 | **validated** |
| plasma proteomics | PXD051911 | **143** | fibrosis / NAS | 0.9284 / 0.9792 | 0.253 / 0.249 |
| H3K27ac CUT&RUN | GSE267145 | **99** | fibrosis 4-lvl / **NAS comp-sum** | 0.7918 / **0.9871** | 0.246 / 0.248 |
| liver proteomics | PXD051911 | 58 | fibrosis / NAS | 0.9107 / 0.9840 | 0.347 / 0.344 |
| scWAT proteomics | PXD051911 | 58 | fibrosis | 0.8975 | 0.348 |

### ⚠ UNDERPOWERED
oWAT 27 (MDE **0.5038**) · plasma follow-up 41 (**not independent** — same participants as the 143) ·
snATAC GSE296875 37 (**binary** any/none, already `UNDERPOWERED_untestable`).

### ⛔ REJECTED
- **RRBS methylation GSE105127 — NO ENDPOINT.** 114 samples carry **exactly three distinct
  characteristic values**, all zone assignments. ⚠ **Its config lists `histopathology` as a modality —
  that is the ZONE call, not a severity score.** A skimmer will think an endpoint exists.
- **PXD052937 — UNRECOVERABLE, closed.** 72 runs in four groups (7/17/38/10); zero hits across all ten
  ConditionSetup columns for any disease keyword; run ids internal with no participant key. ⭐ **And the
  inference is not even self-consistent: the source is a 64-patient BINARY study, the deposit has 72
  runs in FOUR groups.** Neither count nor group number matches. 62 GB, closed.
- single-cell/snRNA (102/26/18/12 — none lists pathology scores), spatial (n=4; two blocked),
  permission-blocked (**`fnih_86_liver` n=86 would otherwise be the richest multi-modal source**),
  seven perturbation cohorts whose units are cells/pools/animals rather than participants.

### ⛔ CROSS-MODAL PAIRING IS ESSENTIALLY UNIQUE
**GSE267145 (RNA + H3K27ac, 99 participants) is the ONLY genuine cross-modal pairing WITH an endpoint.**
GSE296875 pairs snRNA+snATAC on 37 but the endpoint is binary and untestable; GSE105127 pairs RNA+RRBS
on 19 with no endpoint; GSE244832 pairs on 18 with no endpoint.
⛔ **AND TWO OF THE FOUR TISSUES WERE NOT ON DISK.** `data/PXD051911` holds only liver and plasma
quantifications. The `scWAT_proteomics_filename` and `oWAT_proteomics_filename` columns are **`.raw`
INSTRUMENT filenames, not measurements** — nothing in the repo quantified them. **I claimed four
tissues were available on the strength of metadata columns; two of them were one missing download away
from not existing.** The multitissue lane located `scwat_protein_quant.txt` (6988 x 58) and
`owat_protein_quant.txt` (6988 x 27) on **PRIDE PXD051911**, downloaded them into its execution's
`source/`, sha-pinned them, and confirmed both join bijectively (58/58 and 27/27).
✅ **Check that a column names a measurement, not an instrument file, before counting it as a modality.**

⚠ **Cross-TISSUE pairing within PXD051911 is PARTICIPANT-wise, NOT SAMPLE-wise — my earlier framing
was wrong and is corrected here.** Measured from the metadata by the alignment lane:

| tissue | n | visit composition |
|---|---:|---|
| liver | 58 | **100% V1** |
| scWAT | 58 | 32 V1 + **26 V2** |
| oWAT | 27 | **100% V2** |
| plasma | 143 | 136 V1 + 7 V2 |

⛔ **Of the 56 liver∩scWAT shared participants, only 31 are same-visit**; 25 are cross-visit and **7 of
the 56 carry a DIFFERENT Kleiner grade in the two rows.** ⛔ **All 25 liver∩oWAT shared participants are
cross-visit — ZERO same-visit pairs**, because oWAT is drawn entirely at bariatric surgery.

**Consequence**: the liver→scWAT weight transfer I called "the easiest real win" has a usable n of
**31**, not 56 — MDE **0.4850** against a tie ceiling of 0.8975, so the honest expectation is a bounded
null. **The liver→oWAT arm cannot be a claim at all.**

✅ **What survives without pairing**: per-protein severity *associations* can be computed in each tissue
independently (liver 58, scWAT 58) and the association vectors correlated over the shared feature
space. **That concordance quantity needs no paired participants and sidesteps the visit problem
entirely.**

### ⛔ SHIPPING CONSTRAINT on the chromatin instrument
`config/datasets/gse267145_znf469_human_liver.toml`: *"Do not release derivative weights materially
trained on this source until model-specific source terms and legal review pass."* The
coordinate-independent paired-bulk lane **is** activated (TaskSpec `8bf05fe6…`), so measurement is in
scope — **but a chromatin release is not.** Sequence, interval, motif, source-outcome-selection and
external-claim uses also remain blocked.

### ⚠ Three stale records found
1. **`data/REGISTRY.md` lists GSE276114 as "Proteomics (SomaScan), Govaere"** — it is Yang/Zeybel bulk
   RNA-seq. Wrong on modality, platform and author.
2. **PXD051911 is under-registered**: its config lists only liver/scWAT/oWAT. **The plasma arm (143,
   complete Kleiner and NAS) is absent** — the largest non-RNA asset we have is unregistered.
3. GSE105127's `histopathology` field is the zone call (see above).

---

## ⭐⭐⭐ CHROMATIN AND RNA BOTH CARRY SEVERITY BEYOND TECHNICAL — on NAS, n=99 — 2026-09-02

`executions/chromatin-severity-instrument-20260902T011438Z/results_addendum_e_*/`. GSE267145, **99
participants, paired RNA + H3K27ac**, PCA rank-reduced per block to stabilise the alpha selection.

### NAS component sum (ceiling 0.9871) — the well-powered endpoint

| block | ρ | ρ/ceiling |
|---|---:|---:|
| T library complexity | 0.5346 | 0.5415 |
| C H3K27ac | 0.7408 | 0.7504 |
| R RNA | 0.7667 | 0.7767 |
| **RC RNA + H3K27ac** | **0.7926** | **0.8030** |

**Nested increments over the technical baseline, BOTH excluding zero:**

| increment | point | MDE | excludes 0 |
|---|---:|---:|:--|
| **D3 chromatin over complexity** | **+0.1925** | 0.1394 | ✅ |
| **D5 RNA over complexity** | **+0.2378** | 0.1615 | ✅ |

⭐⭐ **Both modalities carry real severity signal beyond library complexity.** This is the first time in
this project that a non-RNA modality has cleared the complexity gate on a nested increment — the
confound that killed the chromatin residual, and that beat the primary modality in the fibrosis arm,
does **not** explain the NAS result.

### Fibrosis (ceiling 0.7918, 72% of participants at stage 0) — everything is a bounded null

| increment | point | MDE | excludes 0 |
|---|---:|---:|:--|
| D3 chromatin over complexity | +0.0096 | 0.1247 | ❌ |
| D5 RNA over complexity | +0.0314 | 0.1338 | ❌ |

⭐ **The endpoint choice decided the result, exactly as the ceiling analysis predicted.** Same
participants, same assay, same folds — a ceiling of 0.9871 versus 0.7918 is the difference between two
clear positives and two bounded nulls.

### ⛔ "Chromatin beats RNA" does NOT survive to the better endpoint

On fibrosis, chromatin led marginally (0.4157 vs 0.2690) and that ordering was stable across three alpha
rules. **On NAS the ordering reverses: RNA 0.7667 against chromatin 0.7408.** So the earlier striking
marginal result is **fibrosis-specific** and must not be reported as a general statement that chromatin
outperforms RNA.

### ✅ The CT anomaly is resolved — it was regularisation instability

CT − C is now **−0.0590** (fibrosis, MDE 0.161) and **−0.0137** (NAS, MDE 0.098), both comfortably
within noise. Rank reduction fixed it. ⚠ **But the producer's own predicted mechanism failed**: repeated
inner CV narrowed the alpha range in only **1 of 7 blocks** (S1 VOID), so *rank reduction* fixed it,
**not** better alpha selection. Its pre-committed adjudication returned **DIAGNOSIS_WRONG** on that
basis and it reported that against itself.
⚠ S5 also voided: the rank-20 RC arm scores below rank-10 on fibrosis but not NAS — higher rank is not
uniformly worse, so the overfitting story is incomplete.

⚠ **Fusion's gain is small and unqualified**: RC 0.7926 against R alone 0.7667 is **+0.026**, far below
the ~0.14 MDEs in play. **Do not report RNA+chromatin as beating RNA** without a nested RC−R increment
and its interval.
✅ `G25_NO_WEIGHTS_DEPOSITED` confirms no release tree and no weights, respecting GSE267145's
derivative-weight restriction.

## ⛔ ATAC (GSE296875, n=39): no severity signal detectable — and it does NOT refute the H3K27ac result

`executions/gse296875-atac-severity-instrument-20260901T012940Z/`. Job 21329615 (7:18).

| endpoint / CV | tech | ATAC | tech+ATAC | nested increment | p |
|---|---:|---:|---:|---:|---:|
| fibrosis 3-level, donor-grouped | −0.291 | −0.268 | −0.260 | +0.031 | 0.281 |
| **steatosis %, leave-one-well-out** | −0.041 | −0.339 | −0.001 | **+0.040** | 0.227 |
| steatosis %, donor-grouped | +0.233 | −0.299 | −0.000 | −0.233 | 0.975 |

⛔⛔ **DO NOT read the negative ρ as inverted biology — and my first framing of why was wrong.**
I wrote that ρ = −0.3 is noise against "a null sd near 0.16 centred at zero." **The measured null is not
centred at zero.** Leave-one-group-out on a group-structured outcome is a *deterministic
anti-correlation generator*: hold out well7 (5/5 positive), train on a negative-enriched complement,
predict LOW for true-HIGH donors. Under a **permuted** outcome with no signal present at all,
ATAC-alone still centres at **ρ = −0.3155 (sd 0.2005)**.

So the observed −0.6135 sits **1.5 sd low against its own artifact-aware null, p 0.9745** — inside the
null's lower tail, not evidence of anything. ⭐ **Raw ρ is uninterpretable under leave-one-group-out;
only the null-referenced quantity means anything.**
⚠ The increment nulls are *narrower* than their components (sd 0.108–0.159 vs 0.20), because a
difference of two positively correlated statistics has less variance than either.

### ⛔ This is a BOUNDED NULL, not a refutation of chromatin

Three reasons it cannot speak to the GSE267145 positive:
1. **Underpowered by construction** — MDE ≈ 0.42 against a ceiling of 0.85–0.99. The H3K27ac increment
   that cleared was +0.1925; **an effect that size was never detectable here.**
2. **Different endpoint.** GSE267145's positive is on the **NAS component sum** (ceiling 0.9871).
   GSE296875 has **no NAS and no fibrosis stage** — its `forbidden_claims` names both — only binary
   any-fibrosis and steatosis percent.
3. **Fibrosis in this cohort is confounded with well** (well-only ρ +0.4775, above its own MDE, with
   all 39 donors nested in 8 wells), so its fibrosis arms were pre-declared void regardless.

⭐ **So chromatin severity remains a single-cohort, single-assay result.** H3K27ac clears the complexity
gate on NAS at n=99; ATAC cannot test it at n=39 on a different endpoint. **Cross-assay support for
chromatin is unavailable, not absent** — and the honest sentence names which.

### Sealed verdicts, and a prediction failure that voided the primary

- **fibrosis (both arms): UNINTERPRETABLE** on two independent grounds — the well-only predictor
  (+0.4775) exceeds the model's own ρ, and prediction P3 failed. Fisher OR **31.5, p 0.00032**; 60% of
  positives vs 4.5% of negatives sit in wells 7–8.
- **steatosis: BOUNDED NULL**, bounded at **MDE 0.3156** under the leave-one-well-out scheme actually
  used (observed −0.0012, p 0.1309, ceiling 0.9863, MDE/ceiling 0.320). The one clean endpoint
  (well KW p 0.4211).
  ⚠ **The 0.4115 I first recorded came from a free permutation carrying no CV structure** — the LOWO
  null is shifted downward (mean −0.2182, sd 0.1911) so the absolute threshold falls with it. ⛔ **My
  prediction that leave-one-well-out would cost power was WRONG: it TIGHTENS the bound by ~0.10.**
- ⛔ **"Fibrosis is untestable by any scheme" is NOT supported** — I proposed that framing and the MDE
  does not approach the ceiling in any arm. The correct, narrower statement: **fibrosis here is
  uninterpretable by EITHER scheme for two DIFFERENT reasons** — donor-grouped lets the well confound
  help the model (P9 confirmed), and leave-one-well-out neutralises the confound as a predictor but
  injects the anti-correlation artifact into the score. **Neither is a power problem.**

### ⛔ AUPRC was quoted against prevalence — the exact trap, and it was self-corrected
First reported as "AUPRC tech 0.6177, ATAC 0.2728 (prevalence 15/37 = 0.405)". **Prevalence is not the
null.** Correct references: **random-scorer AUPRC null mean 0.4582, sd 0.0835, 95th 0.6115** —
reproducing the project's documented 0.459 and sitting well *above* prevalence; and a model-based
permutation null under LOWO of tech 0.3633 / ATAC 0.3448 / tech+ATAC 0.3605.
Re-read correctly: **ATAC's 0.2728 is below its own null mean, p 0.9705.** Conclusion unchanged, stated
reason wrong.

### ⚠ Technical descriptors are the ONLY thing predicting fibrosis here — and adding ATAC destroys it
Tech alone AUPRC **0.6177** against null 0.3633, nominal p 0.0250 (LOWO); 0.6451, p 0.0200
donor-grouped. **tech+ATAC falls to 0.3271, p 0.5272 — a nested increment of −0.2906.**
✅ **But nothing survives multiple testing.** The config mandates BH q 0.05 within the
`GSE296875_histopathology_secondary_family`; across **m = 12** the smallest p is 0.0200 against a
rank-1 threshold of 0.00417, and **every q ≥ 0.1219.** So this is a nominal p only — best read as the
well confound expressed through depth and complexity descriptors.

### ✅ All three admission blockers cleared
Barcode-to-donor join frozen and hashed (`5d3d9466…`, 68,398 unique well-namespaced barcodes,
`bare_barcode_is_a_unique_key: false`); endpoint rows and missingness masks frozen from the
authoritative supplement (`8356c200…` matching the config exactly, arms run separately at n=37 and
n=38, never intersected); donor-grouped folds and the LOWO sensitivity sealed in prespec `771a3e8a…`
and `sha256sum -c` verified at runtime in step 0 of every job.
⚠ **Known deviation, recorded not hidden**: the config also specifies `paired_donor_cluster_bootstrap`
at 10,000 replicates. Permutation p-values and null mean/sd were reported; **bootstrap CIs were not
computed.** Since nothing clears BH, CIs would not change the reading — but it is an unmet config
requirement and is logged as such.

⚠ **P3 predicted the fibrosis increment in [−0.15, +0.20]; observed −0.5620.** By the producer's own
sealed rule that **VOIDS** the primary rather than qualifying it — it had not anticipated the
leave-one-group-out geometry above. ✅ **P9 was correct in all 3 blocks** (donor-grouped beat LOWO
every time), so the well-confound reading is confirmed rather than assumed.

### ⚠ SIDE FINDING — an existing arm elsewhere splits on the confounded variable

`executions/crossmodal-severity-stratification-20260831T190048Z/` has a `fibrosis_any` arm on this
cohort that **splits donors on exactly this well-confounded variable** (same OR 31.5, p 0.00032), so
that stratification is close to a well split. It returned `UNDERPOWERED_untestable`, so **no live claim
rests on it** — but a positive there would have been uninterpretable too, and ⛔ **anyone proposing to
power it up by adding donors from these same wells would be chasing batch.** The steatosis≥5%
stratification is clean (chi² p 0.6223).

⭐ **What would settle chromatin-as-a-pillar**: not more work on GSE296875, but an ATAC cohort whose
batch structure is **crossed with severity rather than nested in it**, at n well above 37.

---

## ⭐⭐⭐ CROSS-MODAL TRANSFER: RNA-trained weights read PROTEOMES — 2026-09-02

`release/masld-liver-protein-severity-v1/`. **Independently verified** by an agent that did not build
it (job 21329671): 58 samples scored, max diff **4.954e-09**, **5/5 refusals fired**, exit 0.

| arm | Spearman vs Kleiner |
|---|---:|
| **T1 released RNA weights applied to protein z-scores** | **+0.6128** |
| T1 with leave-one-out standardisation | +0.6142 |
| T2 literal: RNA's own mu/sd applied to protein log2 | +0.6711 |
| **T3 permuted-weight control (200 draws)** | **−0.0043** (sd 0.1762) |
| T4 protein-native leave-one-out | +0.7291 |
| MDE at n=58 | 0.3608 |

Permutation p **0.00010** (10,000 draws), null mean +0.00035, sd 0.13357.

⭐⭐ **A model trained on transcriptomes predicts fibrosis from protein abundance with no refitting** —
**84% of the protein-native score from 19.6% of the RNA axis** (5,221 of 26,629 genes, carrying 14.35%
of sum|w|). **The permuted-weight control at −0.0043 is what makes this a result**: it isolates the
gene-to-weight pairing from the shape of the protein data.
✅ **Not leakage** — no Kleiner label enters T1, T2 or T3 at any point, and T4 is leave-one-out, so they
are directly comparable.

⭐ **This is the strongest genuinely multi-modal claim in the project.** It says RNA and protein *share
the severity axis*, which is a stronger statement than having an instrument in each.

### ⛔ Two limits that must stay prominent

⛔ **The shipped RNA tool CANNOT score a proteome at all — and the reason is stronger than coverage.**
`masld-severity-v1/score.py` refuses with **`ZERO GENE JOIN: none of the 26,629 axis genes matched your
7,096 identifiers`** — an *identifier* refusal (the axis is unversioned ENSG; the protein matrix is
symbols), not a coverage refusal. **Two separate barriers, and I had conflated them:**
1. **Identifier**: the released scorer rejects symbol input outright. The transfer required a symbol
   crosswalk applied **outside the released path**.
2. **Coverage**: even after crosswalking, only **5,221 of 26,629** axis genes are measured in the
   protein axis (**19.61%**, carrying **14.35%** of sum|w|) — below the release's own measured 30%
   threshold.

**So the claim is about the WEIGHTS, and is further from "the shipped tool does this" than I first
wrote.** ⚠ T3's permuted-weight control has `max_abs` **0.5041** and `pct97_5_abs` 0.3818; T1's 0.6128
exceeds both, so it clears even the control's most extreme draw.
⛔ **No external validation**: n=58, PXD051911 only, against the RNA instrument's three cohorts and 585
participants. The two releases are not at parity and the card must say so.

### ⛔ A release defect I found by running it, which the builder had not
`test_release.sh:5` reads `PY="${PYTHON:-python3}"` and **defaults to a `python3` with no pandas**, so
a stranger running `bash test_release.sh` gets `ModuleNotFoundError` and exit 1. I hit it directly.
The RNA release defaults to the absolute path of a known-good interpreter. **A smoke test that fails
out of the box on the machine it was built on is worse than none** — it teaches the first user that
the release is broken. Fix pending re-verification.

---

## ⛔ FOUR-TISSUE PROTEOMICS: only liver works, and a second tissue HURTS — 2026-09-02 (smoke)

`executions/pxd051911-four-tissue-severity-20260902T012712Z/`. Smoke job 21329624 (14:45, MaxRSS
**249 MB**) ran the full pipeline at reduced draws — **ρ point estimates are deterministic and final**,
intervals are noisy pending the 10,000-draw run.

| tissue | n | ρ | ceiling | ρ/ceil | perm p | **complexity-only ρ** | nested increment |
|---|---:|---:|---:|---:|---:|---:|---|
| **liver** | 58 | **0.7291** | 0.9107 | **0.801** | 0.010 | **+0.4835** | **+0.2413 [+0.069, +0.415]** |
| scWAT | 58 | 0.0917 | 0.8975 | 0.102 | 0.687 | +0.0248 | +0.0703 [−0.293, +0.396] |
| oWAT | 27 | 0.3282 | 0.8835 | 0.372 | 0.299 | +0.1275 | +0.2083 [−0.196, +0.481] |
| plasma | 143 | 0.2120 | 0.9284 | 0.228 | 0.080 | +0.1406 | +0.1465 [−0.103, +0.350] |

⛔ **scWAT is VACUOUS** (ρ 0.092 below the 0.20 instrument floor, p 0.687). **oWAT is the bounded null
promised in advance** (0.328 below its own MDE 0.517). **Plasma sits just under its MDE.**
✅ P1 confirmed: liver reproduces **0.7290918599842632 to the 16th digit**.

### ⛔ The complexity gate bites the ONE tissue that works

**Ten intensity descriptors carrying no protein identity reach ρ 0.4835 in liver — two-thirds of the
full 0.7291.** The nested increment is +0.2413 with an interval excluding zero, so proteins do add real
signal — **but the liver instrument is far less protein-specific than its headline suggests.**
⚠ The producer's own caveat, kept: the descriptors are deterministic functions of the matrix and the
protein block is ~5,300 of 5,310 joint dimensions, so this increment measures **added predictive
value, not independence from complexity.**

### ⛔⛔ SECOND CV-GEOMETRY ARTIFACT OF THE DAY — and it changes how to read every ρ here

**Liver's content-free null centres at −0.4785, not 0.** Under heavy shrinkage a leave-one-out
prediction collapses toward the LOO training mean `(n·ȳ − y_i)/(n−1)`, which is **exactly anti-monotone
in the held-out label**, driving Spearman toward −1.

⛔ **ρ must be compared to the empirical permutation null, NEVER to zero**, and **a negative ρ in a weak
tissue is this artifact, not anti-correlation with fibrosis.**
⛔ **The parametric MDE understates the real bar by roughly 2x at n=58.** Every "bounded null at MDE X"
in this project computed parametrically at small n is therefore **optimistic about its own power.**

⭐ This is the **second independent discovery today** that a CV scheme manufactures anti-correlation —
the ATAC lane found leave-one-**group**-out does it under group-structured outcomes (null centred
−0.3155), this lane finds leave-one-**out** does it under heavy shrinkage. Different mechanisms,
same lesson: **the null must be measured, not assumed.**

### The three cross-tissue questions, all answered NO

**A. Agreement — no, or unresolved.** liver~scWAT same-visit n=31: raw +0.313, partial given Kleiner
+0.359, p 0.085 against MDE 0.485. liver~plasma n=41: raw +0.162, partial +0.066, p 0.33.
scWAT~plasma −0.003. **No pair clears its MDE.** ⛔ **The "technical or batch structure" explanation for the +0.36 is RETRACTED** — it was the
producer's reading, I propagated it, and it does not survive checking. **Acquisition batch is not
associated with Kleiner in any tissue** (Kruskal p 0.16 liver, 0.72 scWAT, 0.41 oWAT), and **liver and
WAT were acquired in separate campaigns** (liver 2019/2020, WAT 2020/2021) with no shared batch
structure. Age/sex/BMI adjustment moves it essentially not at all (+0.359 → +0.353 at n=31;
+0.294 → +0.297 at n=56). ✅ **So it is neither batch nor demography. At n=31 with MDE 0.485 and
p 0.085 the honest answer is UNRESOLVED, not explained.**

**B. Fusion — a second tissue HURTS.** liver+scWAT **−0.1091 [−0.195, −0.056]**, interval excludes zero
and is signed negative; same-visit n=31 −0.1062; liver+oWAT **−0.1927 [−0.397, −0.007]**; liver+plasma
flat (+0.016 / −0.006).

**C. Substitution — clearly no.** Liver 0.77–0.82 of ceiling against scWAT −0.016/+0.077 and plasma
+0.016/+0.061; delta −0.65 to −0.79, every interval excluding zero. **Neither peripheral tissue
substitutes for a liver biopsy at these n.**

⚠ Two predictions FAILED and were flagged rather than qualified: P3 (scWAT 0.15–0.45, observed 0.092)
and P9 (the null-centring above). P3 failed **toward** the null so it cannot manufacture a positive.

### ⭐ MATRIX-FIDELITY POSITIVE CONTROLS — oWAT WITHDRAWN, scWAT narrowed

The producer added a control I had not asked for and it changes two of the four arms. **Sex
classification** is the decisive one: it is unaffected by range restriction or by whether a tissue is
relevant to liver disease, so it measures whether a matrix supports *any* participant-level instrument.

| tissue | n | sex AUC | perm null mean (sd) | p | → fat_pct | → age | → BMI |
|---|---:|---:|---:|---:|---:|---:|---:|
| liver | 58 | **0.9553** | 0.4982 (0.0794) | 0.0005 | **+0.574** | **+0.604** | +0.103 |
| scWAT | 58 | 0.7829 | 0.4978 (0.0818) | 0.0010 | +0.277 | +0.225 | −0.139 |
| oWAT | 27 | 0.7000 | 0.5046 (0.1293) | **0.120** | — | +0.215 | −0.160 |

⛔ **oWAT is WITHDRAWN, not a bounded null.** It **fails its own positive control** — it cannot classify
sex. So ρ 0.328 against Kleiner cannot distinguish "no fibrosis signal in omental adipose" from "at
n=27 this matrix supports no instrument at all." **The correct line is: not evaluable at this n.**
⭐ A bounded null asserts the design *could* have detected an effect; this design cannot be shown to
detect anything.

⚠ **scWAT's null is REAL but narrower than "adipose does not track fibrosis."** The matrix is **not
broken** — it detects sex at p 0.0010, decisively above its own null. So ρ 0.092 against Kleiner is a
genuine null, **at a demonstrated participant-level fidelity of AUC 0.78 against liver's 0.96.**
⚠ **Range restriction is excluded as the explanation**: age sd 13.6 vs 13.5 (ranges 18–68, 19–68) and
BMI sd 6.88 vs 6.94 on the same 58 participants. **The difference is the matrix.**

⚠ **DISCLOSED: a threshold was moved after seeing a number.** The prespec said liver ≥ 0.90 and
scWAT/oWAT ≥ 0.85 means "sound matrix". **scWAT came in at 0.783, below its own pre-declared line**,
and the producer is *not* calling it uninterpretable — arguing the permutation test against its own
null (p 0.0010) is the principled criterion and 0.85 was arbitrary. **That reasoning is sound and the
threshold move is still a threshold move.** ⛔ **Both readings must travel with any scWAT claim** — its
pre-registered fidelity bar was not met.

⚠ Liver's 10,000-draw permutation null: mean **−0.0875**, sd 0.2240, p 0.0104; complexity increment
holds at **+0.2413 [+0.0830, +0.4198]**. Note this differs from the **content-free** null at −0.4785 —
permuted labels and random weights are different constructions and are not interchangeable.

---

## ⛔ FUSION: no measurable gain — and the complementarity gate is empirically disproven — 2026-09-02

`executions/crossmodal-complementarity-20260902T013839Z/`. Development run (B=200); full 10,000-draw
run on 21329697 sharpens intervals only, point estimates are deterministic given the seed.

### ⭐⭐ A synthetic positive control the gate would have THROWN AWAY

**C1** is built so two blocks carry genuinely **independent** components of the same outcome. Its stack
beats the train-selected best single by **+0.2821, CI95 [+0.1360, +0.4084], MDE 0.2125** — a real,
interval-excluding, *detected* complementarity. **Its residual correlation is 0.8217, ABOVE the 0.8
threshold.** Enforcing the gate as a stop rule would have discarded a complementarity the very same
machinery then measured.

⛔ **And every pair that PASSES the gate passes because its second modality has an out-of-fold Spearman
at or below zero** — plasma −0.012, scWAT −0.193. **Both failure directions are now demonstrated
empirically, not argued from algebra.** The criterion in
`config/complementarity_readiness_audit.toml` is defective and must not be used.

| pair | n | ρ_A | ρ_B | r_resid | floor(c=0) | reachable? |
|---|---:|---:|---:|---:|---:|:--|
| GSE267145 RNA × H3K27ac | 99 | 0.288 | 0.420 | 0.945 | 0.852 | NO |
| PXD liver × plasma (48) | 48 | 0.694 | **−0.012** | 0.719 | 0.720 | passes — *because B is null* |
| PXD liver × scWAT (31) | 31 | 0.665 | **−0.193** | 0.759 | 0.710 | passes — *because B is null* |
| **C1 synthetic positive** | 99 | 0.363 | 0.378 | **0.822** | 0.841 | **NO — yet complementarity is real** |
| C2 duplicate (negative ctrl) | 99 | 0.363 | 0.364 | 0.9999 | 0.848 | NO |

### The stack itself: no gain on the only powered substrate

GSE267145, n=99: H3K27ac **0.4199**, RNA 0.2880, **stack 0.4142**, train-selected best single 0.4199.
**Paired delta −0.0058, MDE 0.0470, CI95 [−0.0439, +0.0259].** NNLS kept RNA at non-zero weight in only
**2 of 5 folds**. Consistent with the chromatin lane's own D1 +0.0992 at MDE 0.2595.

### ⛔⛔ A NEW TRAP: a small MDE on a paired delta is NOT "adequately powered"

The protein pairs return MDEs of **0.018–0.065**, and the producer's sealed prediction P-12 (that all
would exceed 0.25) was **wrong**. The reason matters more than the miss: **this is a paired delta
between two scorers that are nearly the same function**, because NNLS zeroed the second modality in
almost every fold. **The tiny interval bounds the gain THIS stack realised; it does NOT exclude
achievable complementarity.** ⛔ Reporting it as an "adequately powered negative" would be wrong, and
the producer flagged it rather than doing so.

⛔ **Two arms VOID, declared rather than reported**: scWAT × oWAT (n=25) had all-zero NNLS weights in
2 of 5 folds making the pooled prediction constant, so its −0.3159 is an artifact; and scWAT→oWAT
transfer is a **numerical no-op** — the coefficient vector was shrunk to zero by its own lambda
selection.

### ⭐⭐ The cross-modal transfer SURVIVES the complexity control

`executions/protein-severity-release-20260902T011748Z/` posthoc (job 21329701). The decisive question
was whether RNA-trained weights reading proteomes are just reading intensity complexity, since liver
complexity alone reaches ρ ≈ 0.48.

| quantity | value |
|---|---:|
| descriptors alone (LOO) | +0.4772 |
| protein-native axis (LOO) | +0.7288 |
| **T1 RNA-weights-on-protein** | **+0.6128** |
| **T1, partialled on Kleiner given the descriptor score** | **+0.5382** |
| protein-native, same partialling | +0.6272 |
| T1 correlation *with* the descriptor score | +0.3590 |
| T1 vs Gini | −0.3621 |
| T1 vs log dynamic range | −0.4944 |

⭐ **Controlling for the complexity score, the transfer retains +0.5382 of its association with
Kleiner** — it loses some but keeps most, and its own correlation with the descriptor score is only
+0.359. **RNA-trained weights are not a complexity proxy.**

⚠ **This is a partial correlation, not the nested increment I asked for.** Partialling controls
linearly; a nested predictive increment is the stronger form and is still worth running. But +0.5382
after control, against a raw +0.6128, is strong evidence the transfer is protein-identity-specific.

⛔ **RETRACTED 2026-09-04 — the paragraph below is wrong and was marked ✅.** Kept in place rather
than deleted, per the no-silent-overwrite rule. It survives only as a record of what was claimed.

> ~~✅ Also resolved: at low axis coverage the **deployed null is not strongly negative**, so the score
> inversion seen under heavy degradation is *a property of the restricted weight vector*, not of the CV
> geometry — the artifact that contaminated two other lanes today does **not** apply here.~~

✅ **What is established instead**: the inversion is neither a property of the weight vector nor of the
LOO geometry. It is an **evaluation-harness artifact**. `score.py` adds a **scalar** intercept, and a
constant cannot move a Spearman, so **the deployed score never inverts at any coverage tested** —
**+0.6255 at 2% coverage, 0 of 200 replicates negative**. The negative curve only ever appeared in a
LOO harness carrying a **per-fold** intercept, which the released scorer does not use.
⚠ The claim was reversed three times before it settled — weight vector (v1) → LOO geometry (v1.1) →
harness (v1.2). A v1.1 "mitigation" reading was also withdrawn once scale was matched: at 2% coverage
+0.2677 raw falls to +0.1592 scale-matched and +0.0450 with the observed scaled down.
See `release/masld-liver-protein-severity-v1.2/MODEL_CARD.md` §4d and
`executions/liver-release-correction-20260903T231124Z/verifier/joints-21412209.out`.

---

## ⭐ LIVER RNA AND LIVER PROTEIN SHARE A SEVERITY DIRECTION — but not a fibrosis-specific one

`executions/mash-masl-cross-cohort-transfer-20260902T0200Z/`, prespec `1612e435` + 3 addenda.
Liver RNA (**521 participants, 5 cohorts, F0–F4**, cohort-adjusted per-gene association) against liver
protein (PXD051911, n=58, Kleiner), **3,941 shared symbols, NO shared participants.**

**Concordance Spearman +0.4523.** Removing the top 5% of movers on either side lowers it to **0.365
without collapsing** — the agreement is **broad, not carried by an ECM/collagen handful**.

### ⛔ The specificity control is the finding: it is NOT fibrosis-specific

RNA side held fixed at fibrosis; only the protein label changes:

| protein label | concordance with RNA-fibrosis |
|---|---:|
| **Kleiner fibrosis** | **0.4523** |
| lobular inflammation | 0.4272 |
| ballooning | 0.4020 |
| NAS | 0.3913 |
| steatosis | 0.3586 |

**Fibrosis wins by 0.025.** The five labels are **not collinear** (median pairwise ρ ≈ 0.6), so the
control is informative. ⭐ **The licensed claim is "liver RNA and liver protein share a general severity
direction", NOT "the modalities agree about fibrosis".** Concordance only — no prediction language.

⚠ **DEMOTED TO EXPLORATORY by the producer's own rule.** Two sealed predictions were wrong *in the
optimistic direction* (predicted 0.05–0.35, observed **0.45**; predicted null sd < 0.03, observed
**0.056**), and its prespec demotes an arm whose predictions fail. **Being wrong optimistically is why
it added the specificity control before the full run.**

## ⛔⛔ A THIRD null-construction finding: the label-only MDE is the WRONG bar for a TRANSFER claim

A model trained on **permuted training labels** and scored on the **untouched test cohort** produces a
null with **sd 0.16–0.20**, against the label-only null's 0.105 and 0.058. **A random direction in real
expression space picks up whatever dominant structure the test cohort has**, so its AUROC swings far
more than iid noise. Confirmed by a random-direction diagnostic with **no training at all**, which
reproduces the same width — so it is the geometry of the expression matrix, not the fitting.

⭐ In the smoke, arm A2 reached **AUROC 0.82** and cleared "is this score associated with the label"
at p 0.002 — but did **NOT** clear "did the real training labels matter."
**Verdict: `SCORE_TRACKS_LABEL_but_TRAINING_NOT_SHOWN_TO_MATTER`.**
⛔ **Three distinct null-construction defects found in one day**: leave-one-group-out anti-correlation,
leave-one-out-under-shrinkage anti-correlation, and now a transfer null that must permute *training*
labels rather than compare against a label-only bar.

### The MASH/MASL unlock is real, but smaller than I briefed
- **GSE167523: n=98** (NAFL 51 / NASH 47), labels cross-checked against **two independent sources**
  with **0 disagreements**. ✅ Opening the deposited workbook **confirms** the earlier rejection rather
  than overturning it — three columns, no fibrosis stage.
- ⛔ **GSE126848 usable n is 31 (15 vs 16), NOT the 53 I gave.** The other 25 are healthy/obese — a
  *disease-presence* contrast, not MASH-vs-MASL, and correctly not pooled.
- ✅ Unusually clean substrate: both in the same frozen set, 86,369 genes, identical axis order,
  featureCounts v2.1.1, GENCODE v49, `-s 2`, single-end. **No crosswalk needed.**
- ⛔ **Cohort is perfectly confounded with instrument** (NextSeq 500 vs HiSeq 3000). Unremovable.
- MDEs from class marginals, before any fit: A1 (98→31) AUROC **0.7820**; A2 (31→98) **0.6643**.

⛔ **SEVENTH lane**: in arm A1 the 12 technical library-complexity descriptors alone (**AUROC 0.729**)
**BEAT the genes (0.704)**.

### ✅ The RNA-vs-protein agreement question was BLOCKED and correctly left blocked

`blocked_question_2.json`: *"PXD051911 contains no RNA for any participant."* Status **BLOCKED**, and
the artifact records what was **not** done: *"No crosswalk, imputation, proxy or surrogate was
constructed in place of a measurement that was never taken. Nothing was worked around."*
⭐ **That is the right response to an impossible question** — the alternative is a proxy that would have
looked like an answer. It also names where it *can* be answered (GSE267145, another lane's substrate)
rather than leaving the question open-ended.

## ⭐⭐ WHY fusion fails: both modalities make the SAME error — 2026-09-02 (final)

Job 21329712 (10:53, peak RSS 1.10 GB against a 48 G ask), **7/7 guards demonstrably_live**, 10 of 12
sealed predictions correct. Results sha256 `44f9648d…`.

### The complementarity map — error in rank percentile by stage, denominators printed

| stage | n | RNA \|e\| (signed) | H3K27ac \|e\| (signed) | closer: RNA vs H3K27ac |
|---|---:|---|---|---|
| F0 | 71 | 0.2355 (**+0.0889**) | 0.2312 (**+0.0687**) | 38 vs 31 |
| F1 | 15 | 0.2660 (−0.1865) | 0.2391 (−0.1838) | 6 vs 8 |
| F2 | 9 | 0.2222 (−0.1841) | 0.1740 (−0.1470) | 2 vs 6 |
| F3 | 4 | *not estimable* 0.4646 (−0.4646) | 0.2020 (−0.1995) | 0 vs 4 |

⭐ **Both modalities over-rank F0 and UNDER-rank every higher stage, in the same direction, on the same
participants.** There is **no A-fails-high / B-fails-low structure to exploit.** That shared,
stage-correlated error **IS** the 0.9448 residual correlation — **and it is why no combiner helps.**
**This is the mechanism, not just the outcome.**

### ⛔ A FOURTH null-centring finding: the stack-gain null is CENTRED NEGATIVE
Permutation null of the gain: **mean −0.0953**, sd 0.0921, p97.5 +0.1217, max +0.3385, p 0.141.
⛔ **"The stack is slightly worse than the best single" is exactly what a no-signal stack does.**
Reporting the observed **−0.0058 without the null mean would have read as harm** — it is in fact
*better* than the null. Stack delta MDE 0.0470, CI95 [−0.0441, +0.0260],
`NEGATIVE_adequately_powered`, and NNLS kept RNA non-zero in 2 of 5 folds — **a stack that used both
modalities and gained nothing**, not a degenerate one.

### ⭐ Partial correlation is the right statistic, and it is asymmetric
| | partial Spearman of fibrosis |
|---|---:|
| on **H3K27ac given RNA** | **+0.3194** |
| on **RNA given H3K27ac** | **+0.0130** |
| on liver given plasma | +0.7019 |
| on plasma given liver | −0.1504 |

⭐⭐ **Chromatin carries ordering that RNA does not; RNA carries essentially nothing chromatin does
not.** Reproduces the chromatin lane's D1 +0.0992 / D2 −0.0327 **from an independent statistic**.
✅ Use the partial correlation, not the residual correlation — it is scale-free, directional, and
answers the question the residual correlation buries.

### ⛔ My four-tissue overlaps were participant-level unions ACROSS VISITS
Visit-matched: **liver∩oWAT = 0**, liver∩plasma 41, liver∩scWAT 31, scWAT∩oWAT 25, scWAT∩plasma 32,
**all four = 0**. ⛔ **The 56 and 23 I quoted are unions across visits, not usable pairs.** Every pair
now reported both ways.

### Other closures
✅ **DPI chain confirmed empirically**: RNA → imputed chromatin → severity is **0.2100** against direct
RNA 0.2880, delta **−0.0780**, MDE 0.1478. Does not beat direct RNA, as predicted.
⚠ **liver→oWAT transfer rescue**: bounded null (5,715 shared symbols, n=27, delta −0.1278, MDE 0.3022).
⛔ **scWAT→oWAT is a NO-OP**: coefficient norm **3.8e-04**, lambda pinned at the grid maximum 1e6 in
**27 of 27** folds, delta exactly 0.0 with a zero-width interval. **Do not read that zero as a null.**
⚠ **And its automatic no-op flag FAILED**: it tests *value* identity, but the operative degeneracy is
*rank* identity, which Spearman sees and the flag did not. Caught by the printed norm and lambda.
⛔ **One arm VOIDED by its own guard**: scWAT × oWAT (n=25), all-zero NNLS weights in folds 3 and 4 made
the pooled prediction constant, so its −0.3159 is an artifact — `VOID_degenerate_stack`.
⚠ **ATAC instrument BLOCKED with an exact requirement**: needs the GSE296875 Data S1 workbook
(`fibrosis_categorical`, 37 of 39 donors) joined on `donor_id`. **Nothing local carries that label.**

### ✅ CROSS-LANE REPRODUCTION: two lanes, independent implementations, bitwise agreement

`results_B_20260902T013839Z` (job 21330965). Anchors recomputed by the fusion lane against the values
shipped by the lanes that produced them:

| anchor | reproduced | shipped | abs diff | tolerance |
|---|---|---|---:|---:|
| plasma tie ceiling, n=143 | 0.9284102359253253 | identical | **0.0** | 1e-09 |
| liver LOO Spearman, n=58 | 0.7290918599842632 | identical | **0.0** | 1e-06 |

⭐ **This is what makes every disagreement elsewhere in this session readable as a SPECIFICATION
difference rather than broken machinery.** Without a reproduce-the-other-lane arm, a conflicting number
is uninterpretable.

Three diagnostics worth carrying:
⚠ **The lambda-selection choice costs 0.006012** — the same recipe with LOO lambda selection gives
0.723079 against the shipped 0.729092. Small, but exactly the kind of unstated implementation detail
that becomes an irreproducible number in someone else's hands. **State the selection rule.**
⚠ **Dropping ONE participant moves the plasma ceiling to 0.9271660676** (from 0.9284102359). At n=143 a
single participant shifts it by 0.0012, so **ceilings quoted to four decimals are over-precise.**
⭐ **Reversing the prediction vector gives −0.078939, NOT −0.729092.** The negative of a good predictor
is not a correspondingly bad one under LOO shrinkage — the same anti-correlation geometry found in two
other lanes today, and a further demonstration that **ρ must never be compared to zero here.**

## ⛔⛔ PANDAS SILENTLY DESTROYS "None" AS A CATEGORY — project-wide hazard — 2026-09-02

**`pandas`' default `na_values` contains the literal string `"None"`.** In GSE296875 that is a **real
steatosis category**, and a default `read_csv` turned **18 of 38 observed donors into NaN**, merging
them with the 1 genuinely missing donor. Found and fixed by the fusion lane in its own code; its
n=38 analysis set was never affected, only the stratum map, which had silently covered 20 of 38.

⛔⛔ **THREE columns are damaged, not one** — measured on the endpoint lock, default read vs
`keep_default_na=False`:

| column | default parse | correct parse | destroyed |
|---|---|---|---|
| `fibrosis_categorical` | NaN 24, yes 8, mild 7 | **None 22**, yes 8, mild 7, NaN 2 | **ALL 22 NEGATIVES** |
| `steatosis_categorical` | NaN 19, Mild 13, Severe 4, Moderate 3 | **None 18**, … NaN 1 | 18 |
| **`steatosis_status`** | NaN 5, … | **None 4**, … | 4 |

`config/evaluation/gse296875_histopathology.toml:34` is `negative_source_values = ["None"]`, so a
default read of `fibrosis_categorical` **deletes the entire declared negative class** and merges it
with the 2 genuinely missing donors — producing a plausible-looking 24.

⭐⭐ **The rule is STRUCTURAL, and this corrects my own framing.** I called the ATAC lane's escape
"luck". It was not. **BOOLEAN columns are safe; `*_categorical` and `*_status` STRING columns are
not.** `fibrosis_any` parses to a real bool — **22 False / 15 True / 2 NaN, identically under both
settings** — so that lane was fine **by construction**. `fibrosis_status` is safe only because its free
text never happens to be the bare word.
✅ **Anyone reading `fibrosis_any` is unaffected by construction; anyone reading `fibrosis_categorical`
is destroyed by construction.** ⭐ **That tells you WHICH lanes to re-check rather than all of them.**

✅ **The fix that generalises**: assert the strata **SUM TO n**. Guard B04 now prints 38/38, shortfall 0,
on all four arms. **A map that covers 20 of 38 and says nothing is the failure mode** — the same
printed-denominator discipline that caught three other errors today.

### Two more self-found defects in the same lane
⛔ **The instrument gate used |ρ|**, so a model at ρ **−0.269** passed as "working". Now signed — and it
changes a real reading: the leave-one-well-out arm now prints *"SIGNED +0.20: False (an absolute-value
gate would have said True)"*.
⛔ **Sealed predictions Q-1..Q-10 had no coded adjudication.** A sealed prediction nobody scores is not a
prespecification. Now coded and run.

## ⛔ GSE296875 steatosis: NO working instrument — untestable, not null

| arm | ρ(RNA) | ρ(ATAC) |
|---|---:|---:|
| all-cell | +0.0505 | −0.0503 |
| hepatocyte-only *(post-hoc, declared)* | −0.0116 | −0.1422 |
| leave-one-well-out | −0.1713 | −0.2690 |

Bare-Spearman MDE at n=38 is **0.4411**; the project's own shipped instrument floor is **0.20**.
**All four arms fail both, on the signed value.** Every stack verdict there is
`UNTESTABLE_no_working_instrument` — **complementarity is untestable, not absent.**
⭐ **The post-hoc hepatocyte-only rescue did NOT rescue it**, which makes the substrate-level negative
*stronger*. Q-4 and Q-6 were wrong and voided — **both in the direction of expecting the substrate to
work.**

⭐⭐ **And it is the sharpest demonstration of the criterion defect in the whole lane**: with both
instruments at zero, the attainable residual-correlation floor is **0.997457** and the observed
**0.996224 sits BELOW its own floor.** *"Nothing about that number is a statement about
complementarity; it is a restatement of 'both residuals are the outcome.'"*

✅ **Criterion defect formally recorded** in `criterion_defect_record` — file, field, sealed value 0.8,
verified_date, the algebra, the ρ > 0.408248 equal-strength requirement, all six measured floors, three
failure directions including the C1 counterexample, and the replacement. ⛔
`config/complementarity_readiness_audit.toml` was **NOT edited**, and the record says so.
✅ Suggested repair if a threshold is wanted at all: **per pair, as floor(c=0) minus a fixed margin —
never one constant across substrates.**

---

## ⛔ FOUR-TISSUE PROTEOMICS, FINAL (10,000 draws): only liver works, a second tissue HURTS — 2026-09-02

Job 21329678, COMPLETED 01:03:19. Supersedes the smoke figures; point estimates unchanged, intervals
now final.

### Instruments
| tissue | n | ρ | ceiling | ρ/ceiling | MDE | complexity increment | verdict |
|---|---:|---:|---:|---:|---:|---|---|
| **liver** | 58 | **0.7291** | 0.9107 | **0.801** | 0.361 | **+0.2413 [+0.083, +0.420]** | works |
| scWAT | 58 | 0.0917 | 0.8975 | 0.102 | 0.361 | — | ⛔ **GATE FAILED, VACUOUS** |
| oWAT | 27 | 0.3282 | 0.8835 | 0.372 | **0.517** | +0.2083 [−0.172, +0.568] | ⛔ **withdrawn** (fails sex control) |
| plasma | 143 | 0.2120 | 0.9284 | 0.228 | 0.232 | +0.1465 [−0.055, +0.353] | under its MDE |

⛔ **Only liver clears its complexity gate.** oWAT's and plasma's increments both cover zero.

### Agreement — 10 pairs, NOT ONE clears its MDE
| pair | n | raw | CI95 | partial\|Kleiner | MDE | p |
|---|---:|---:|---|---:|---:|---:|
| liver~scWAT same-visit | 31 | +0.3133 | [+0.0012, +0.5541] | +0.3589 | 0.4850 | 0.086 |
| liver~oWAT | 25 | +0.3977 | [−0.012, +0.692] | +0.2628 | 0.5351 | 0.050 |
| liver~plasma same-visit | 41 | +0.1617 | [−0.188, +0.498] | +0.0664 | 0.4256 | 0.303 |
| scWAT~plasma same-visit | 32 | **−0.0033** | [−0.344, +0.334] | −0.0109 | 0.4779 | 0.987 |
| oWAT~plasma | 24 | −0.0974 | [−0.522, +0.364] | −0.2602 | 0.5451 | 0.647 |

⚠ liver~scWAT's CI *just* excludes zero (+0.0012) but sits **well below its MDE of 0.485** — a lucky
draw, not a detection.

### ⛔ Fusion: adding a second tissue ACTIVELY HURTS
| pair | n | liver alone | joint | increment | CI95 |
|---|---:|---:|---:|---:|---|
| **liver+scWAT (visit-mixed)** | 56 | 0.7536 | 0.6444 | **−0.1091** | **[−0.181, −0.051]** ⛔ excludes zero |
| liver+scWAT (same-visit) | 31 | 0.7126 | 0.6065 | −0.1062 | [−0.230, +0.004] |
| liver+oWAT | 25 | 0.5615 | 0.3688 | −0.1927 | [−0.396, +0.003] |
| liver+plasma | 41/48 | 0.632/0.649 | 0.648/0.643 | +0.016 / −0.006 | flat |

⛔ **The visit-mixed liver+scWAT increment is negative AND excludes zero.** Adding a vacuous second
tissue does not merely fail to help — it **degrades** a working instrument. ⚠ The cleaner same-visit
arm (n=31) covers zero, so the strongest statement rests on the visit-mixed arm.

### ⛔ Substitution: decisive no
scWAT vs liver **−0.7896 [−1.210, −0.328]** and **−0.7452 [−1.025, −0.448]**; plasma vs liver −0.6799
and −0.6496, all excluding zero. **No peripheral tissue substitutes for a liver biopsy at these n.**
⚠ oWAT vs liver (−0.1061, covers zero) looks closest but oWAT is **withdrawn** on its sex control.

### ⭐ Two controls, and only liver clears BOTH

| tissue | ρ | content-free MAX | clears content-free? | perm p | clears null? |
|---|---:|---:|:--|---:|:--|
| **liver** | **0.7291** | 0.3890 | ✅ | **0.0104** | ✅ |
| scWAT | 0.0917 | 0.3743 | ❌ | 0.650 | ❌ |
| oWAT | 0.3282 | 0.2610 | ✅ | 0.298 | ❌ |
| plasma | 0.2120 | 0.3335 | ❌ | 0.092 | ❌ |

⚠ **PLASMA IS CONTESTED and is reported as NOT ESTABLISHED.** Its bootstrap CI **excludes zero**
(+0.043, +0.376) — but perm p = 0.092 **and** it fails its content-free max. ⭐ **The bootstrap CI does
not account for the LOO shrinkage artifact; the permutation null does**, so the permutation is the
calibrated reference. ⛔ **Two controls disagreeing is a "not established", not a positive.**
⚠ oWAT is the mirror case — beats its content-free max but not its permutation null. **The permutation
null is preferred because it preserves X's covariance.**

### ⭐⭐ The liver complexity signal is COMPOSITION, not depth — this softens the confound reading

Decomposing the ρ 0.4835 that complexity alone achieves in liver:
**Gini −0.539 and entropy +0.381 lead it; library depth is only +0.240.**

⭐ **So it reads as a genuine compositional shift — a fibrotic liver is less dominated by a few
hepatocyte proteins — rather than a pure technical artifact.** That is a materially different claim
from the ATAC and MASH/MASL lanes, where depth-like descriptors carried it. **The confound family is
not homogeneous: sometimes it is technical, here it is biology the descriptors happen to capture.**
⚠ Caveat kept: the descriptors are deterministic functions of the matrix and the protein block is
~5,300 of 5,310 joint dimensions, so the +0.2413 increment measures **added predictive value, not
independence**.

### ⚠ NAS beats fibrosis in liver — and it is TIE STRUCTURE, not better prediction
Liver NAS **+0.8147** exceeds its fibrosis 0.7291. But **ρ/ceiling is 0.828 vs 0.801** — NAS has 9
levels (ceiling 0.984) against Kleiner's 4 (0.911). ⭐ **The apparent advantage is the tie structure.
This is exactly what the ceiling exists to expose**, and it is why raw ρ across endpoints is not
comparable.

⚠ **Four of nine sealed predictions wrong** (scWAT 0.15–0.45 → 0.092; partial|Kleiner ±0.20 → +0.359;
fusion −0.109 outside [−0.10, +0.05]; null shape badly). ✅ **All four run toward LESS signal than
predicted, so none can manufacture the negative A/B/C findings**, and the liver machinery control
reproduces **bitwise**. 20 guards live (18 demonstrably_live, 2 live_only_under_a_changed_constant).

### ✅ THE AGE CONFOUND IS REFUTED — liver protein is not an age proxy

I raised this because liver protein predicts **age at +0.6037** against Kleiner at +0.7291, and age and
fibrosis are usually correlated. **They are not correlated here.** Clinical marginals against Kleiner,
n=58:

| covariate | ρ vs Kleiner |
|---|---:|
| age | **−0.0388** |
| BMI | −0.1537 |
| sex | −0.1047 |

⭐ **Age barely relates to Kleiner in this cohort**, so the proteome's age signal is **orthogonal** to
its fibrosis signal. The confound I feared does not exist in this substrate.

| model | ρ vs Kleiner | vs NAS |
|---|---:|---:|
| clinical only (age+sex+BMI) | **−0.4470** | −0.1384 |
| clinical + protein | **+0.7277** | +0.8157 |
| *(protein alone, for reference)* | *+0.7291* | *+0.8147* |

⛔ **The reported increment of +1.1747 [+0.856, +1.432] must NOT be quoted bare.** It is inflated by an
artifactually **negative** baseline: a 3-feature clinical model with no real signal collapses under LOO
shrinkage to the anti-monotone leave-one-out mean — **the same −0.48-centred artifact found in this
lane and two others today.** An increment measured from −0.447 is not a 1.17-point gain in any
meaningful sense.

✅ **The honest statement is the comparison that needs no baseline**: liver protein scores **+0.7291
alone and +0.7277 with clinical covariates included** — a difference of 0.0014. **Adding age, sex and
BMI changes nothing, and the instrument is not riding on them.**

### ✅ scWAT's mapping VERIFIED by an X-linked protein — explanation 2 holds

⭐ **The elegant check**: liver's top sex discriminators are **ARSD (+0.74), GYG2 (+0.63), RPS4X
(+0.59)** — all X-linked. **scWAT's top hits include GYG2 (+0.49), the same X-linked glycogenin-2.**
**A scrambled sample-to-column mapping cannot recover an X-linked protein as a top sex
discriminator.** The join is correct. (oWAT's top hits — BTD, RPL19, CERT — are nothing sex-linked,
consistent with its failed sex control.)

⚠ Mis-mapping calibration, for completeness: scrambling liver takes sex AUC 0.955 → 0.900 (10%) →
0.796 (20%) → 0.750 (30%) → 0.541 (75%). scWAT's 0.783 *would* correspond to ~20–25% mis-mapping —
**but only if adipose and liver were equally sex-dimorphic, which there is no reason to assume**, and
GYG2 rules it out directly.

✅ **So explanation 2 holds: the matrix is sound, and subcutaneous adipose proteome simply carries
little about liver fibrosis in this cohort.** VACUOUS stands on p 0.687, and **scWAT does not need
revisiting on suspicion of an upstream defect.**

QC also clears it — scWAT is best of three on sample-sample correlation (0.8936 vs liver 0.8816,
oWAT 0.8847), dynamic range 17.14, missingness 0.211.

### The 7 cross-visit stage changes are GENUINE, not a labelling defect
ID2 F2→F4, ID15 F2→F3, ID11 F1→F2, ID31/ID42/ID47 F0→F1, **ID50 F2→F1**. **Six increases, one
decrease** — consistent with progression between visit 1 and bariatric surgery, and the single
decrease sits inside the deposit's own σ_e (0.4655–0.5092). **No arm is affected** (every instrument
predicts the Kleiner on its own row); the count appears as `n_label_mismatch` in the agreement arms.

### ⭐ The two nulls are NOT interchangeable, and oWAT proves it
The permutation null **keeps X's covariance** and is the calibrated reference for every p. The
content-free null **destroys it** and is a strictly easier bar, quoted only as max/percentiles.
⭐ **oWAT clears the easy bar (0.3282 > 0.2605) and fails the real one (p 0.298)** — which is exactly
why they cannot be swapped.

### ⛔ MY BMI REASONING WAS WRONG — retracted (it was stated verbally, never written here)

I argued that **both** adipose tissues failing the BMI positive control (scWAT −0.139, oWAT −0.160)
made an upstream data defect more likely than two independent biological nulls.

⛔ **Liver fails the BMI control too: liver → BMI = +0.1028.** Every participant in this deposit is
severely obese (**median BMI 41–42, IQR ≈ 39–45**), so **BMI is a poor control target for EVERY tissue
here** — that is range restriction, not tissue quality. Adipose failing it is **not diagnostic of
anything**, and my inference from it was invalid.

✅ **The controls that DO discriminate are sex (0.955 / 0.783 / 0.700) and age (+0.604 / +0.225 /
+0.215)**, and the oWAT withdrawal rests on those, not on BMI.
⭐ **Choose a positive control the cohort can actually vary on.** A control target that is restricted by
the cohort's own selection criteria cannot distinguish a broken matrix from a working one.

## ⭐⭐⭐ THE LIVER "COMPLEXITY" SIGNAL IS FIBROGENESIS — tested, not inferred — 2026-09-02

⚠ **The producer stopped me from building a seven-lane narrative on an inference.** Its line that the
Gini/entropy signal "reads as a compositional shift" was an **inference** — it had measured only which
descriptors lead, and a dynamic-range shift could equally come from differential sample handling.
Because I was about to give it cross-lane weight, **it tested it** (job 21331010). **It holds.**

**Top 60 liver proteins RISING with Kleiner — 22 are ECM, 5 collagens:**
TGFBI **+0.77**, COL14A1 **+0.77**, FBLN5 +0.76, NID2 +0.75, LTBP4 +0.74, EMILIN1 +0.73, DCN +0.72,
LAMC1 +0.71, NID1 +0.70, LUM +0.70. **A textbook fibrosis matrix panel.**

**Top 60 FALLING — ZERO collagens, all hepatocyte metabolic enzymes:**
AOX1 −0.67, SORD −0.62, AHCY −0.62, GSTZ1 −0.57, SELENBP1 −0.57, CBS −0.56, PCK2 −0.56, HAAO −0.56.

⭐⭐ **ECM accumulating while hepatocyte metabolic protein is lost MECHANICALLY flattens the intensity
distribution** — Gini down, entropy up, top-10 share down. **That is the same event the descriptors
measure.** So the composition reading is now supported by **protein identity**, not merely by which
descriptor ranked first.

⭐ **This substantially strengthens the liver result rather than qualifying it: what looked like a
technical confound is, in this tissue, the disease process itself.** It is why the depth-driven vs
composition-driven distinction matters and why the seven lanes must not be written as one phenomenon.

### ⛔ Multiplicity, which must travel with the decomposition
**10 descriptors x 4 tissues = 40 univariate tests, uncorrected.** At n=58 the nominal two-sided bar is
|ρ| **0.2596**; Bonferroni over 40 is **0.4275**.

| descriptor | liver ρ | clears |
|---|---:|---|
| **Gini** | **−0.539** | **nominal AND Bonferroni** |
| entropy | +0.381 | nominal only |
| top-10/50/100 share | −0.31 to −0.37 | nominal only |
| **library depth** | **+0.240** | **NEITHER** |
| missingness | −0.240 | neither |

✅ **This cuts in favour of the reading, not against it: depth does not clear even the nominal bar.**
So the supported claims are **the ordering** (dynamic range leads, depth does not) and **the Gini
result**. ⛔ **Do not lean on entropy or the top-N shares individually.**

⚠ **The producer corrected its own method**: its keyword tagger counted only 5 of 60 negatives as
"hepatocyte" because it keyed on CYP/ADH/ALDH prefixes — **AOX1, SORD, AHCY, GSTZ1, SELENBP1, PCK2,
AADAT and HAAO are all hepatocyte metabolism and match none of them.** ⛔ **Read the protein
identities, not the tag counts.**
⚠ The +0.2413 increment caveat is unchanged: **added predictive value, not independence.**

---

## ⛔ MASH/MASL CROSS-COHORT TRANSFER: does NOT clear the strict null — 2026-09-02 (final)

`executions/mash-masl-cross-cohort-transfer-20260902T0200Z/results/transfer_results.json`.
Job 21329683, COMPLETED 01:32:24. The lane that was meant to convert two rejected cohorts into external
validation on a coarser endpoint.

| arm | direction | **T** (technical) | **G** (genes) | **GT** | MDE | verdict |
|---|---|---:|---:|---:|---:|---|
| **A1** | GSE167523 (98) → GSE126848 (31) | **0.7292** | 0.7042 | 0.7042 | 0.7820 | ⛔ `UNDERPOWERED_untestable` |
| **A2** | GSE126848 (31) → GSE167523 (98) | 0.7443 | **0.8206** | 0.8202 | 0.6643 | ⚠ `SCORE_TRACKS_LABEL_but_TRAINING_NOT_SHOWN_TO_MATTER` |

⛔ **In BOTH arms `observed_exceeds_null_max: false` against N2** — the strict transfer null that asks
*"did the real training labels matter?"* A model trained on **permuted** training labels scores about
as well on the untouched test cohort. **So the score is picking up dominant expression structure that
happens to align with the label, not learned biology.**

⭐ **A2 reaching AUROC 0.8206 and still failing is the whole point of the strict null.** It clears
"is this score associated with the label" comfortably, and the transfer claim still does not stand.
Most analyses would have reported 0.82 and stopped.

⛔ **In A1, technical descriptors alone (0.7292) BEAT the genes (0.7042)** — the seventh lane, and here
the technical arm wins outright while the gene arm is untestable anyway.

### What this means for the modality inventory
⛔ **The MASH/MASL unlock does NOT deliver external validation.** GSE126848 (n=31 usable, not the 53 I
first quoted) and GSE167523 (n=98) **remain unusable for a severity claim** — now for a measured
reason rather than a missing-endpoint one. **External validation stands at three cohorts, all RNA
fibrosis: GSE268273, GSE213621, GSE276114.**
⚠ C1 concordance is labelled **`FORM 1 CONCORDANCE. This is NOT a prediction and may not be written as
one.`** — the RNA↔protein severity-direction result (0.4523) keeps that scope.

### ✅ Gate repair (job 21331024): a dead gate made live, and genes add NOTHING over technical

A previously non-firing gate is now live — `PD1_gate_now_live_corr_below_0.99`, A1 corr **0.9344**,
A2 corr **0.9224**, both under the 0.99 threshold, so the gate can now fail.

**Gene increment over the technical baseline, both arms:**

| arm | increment over T | CI95 |
|---|---:|---|
| A1 | −0.0042 | [−0.214, +0.223] |
| A2 | +0.0334 | [−0.057, +0.126] |

⛔ **Genes add nothing over technical descriptors in either direction.** Combined with A1's technical
arm (0.7292) **beating** its gene arm (0.7042), the MASH/MASL picture is internally consistent:
**the apparent transfer signal is carried by technical descriptors, the genes contribute no detectable
increment, and the transfer fails the strict permuted-training-label null.**

⭐ Three independent lines — the marginal ordering, the nested increment, and the strict transfer null —
**all point the same way.** That is what a well-constructed negative looks like.

## ⛔⛔ CROSS-TISSUE CONCORDANCE FAILS — and the null I specified was 5–12x too narrow

`executions/mash-masl-cross-cohort-transfer-20260902T0200Z/results/c3_crosstissue_concordance.json`,
addendum E `246eac5b`.

**I asked for the wrong null.** I specified *permute protein identity within one tissue*. The correct
null is *permute the Kleiner label within each tissue and recompute BOTH vectors*.

| pair | ρ | n shared | **N-identity (mine)** | **N-label (strict)** | sd ratio |
|---|---:|---:|---:|---:|---:|
| **liver~scWAT (primary)** | +0.123 | 2,384 | p 0.0005 | **p 0.254** | **5.35x** |
| scWAT~oWAT | +0.268 | 2,797 | p 0.0005 | p 0.262 | **12.40x** |
| liver~oWAT | +0.224 | 3,017 | p 0.0005 | p 0.033 ✅ | 5.58x |
| liver~plasma | +0.150 | 163 | p 0.054 | p 0.190 | 1.51x |
| scWAT~plasma | +0.281 | 175 | p 0.0005 | p 0.013 ✅ | 1.56x |

⛔ **My null declared 4 of 5 significant. The correct null clears 2 — and the primary liver~scWAT arm
is not one of them.**

⭐ **Why**: identity-permutation **destroys each vector's internal protein-correlation structure**;
label-permutation **preserves** it. **Protein correlation structure is shared across tissues** —
abundance and technical gradients behave alike everywhere — so **two independent random-label
association vectors agree substantially by chance.** That is the 0.109–0.235 null sd, and ρ 0.123 sits
inside it.

⛔ **Nothing here supports "the tissues agree about severity."** ✅ And the two pairs that DO clear are
the two pre-declared **DESCRIPTIVE** in the sealed prespec — oWAT at n=27, and plasma pairs carrying
**under 200 shared features**. The producer refused to upgrade them on a p-value.

### ⭐⭐ THE STANDING RULE, third instance in one execution
**A null that holds one thing fixed while the real analysis CHOSE it is measuring the wrong thing.**
⭐ **If a null does not regenerate the quantity end to end, it is the wrong null.** Three instances in
this single execution: the label-only MDE against a *trained* transfer; the shared-penalty gate where a
narrow block could not move the fit; and now identity-permutation against a *label-derived* vector.

⚠ Analysed features are smaller than raw accession overlap — complete-case filtering runs **within each
tissue first**: liver 4,019/7,096; scWAT 2,814/6,988; oWAT 4,032/6,988; plasma 248/427. The liver~scWAT
figure is **2,384**, not the raw 5,723 overlap.
✅ Guards all live: joins bijective in all four tissues; **follow-up quarantine 0 of 41 leaked**;
identity permutation reproduces every observed ρ to <1e-10; scWAT/oWAT digests match the four-tissue
lane's published pins exactly.
✅ **Independently reproduced from the other lane**: visit composition (liver 58/58 V1, scWAT 32 V1 +
26 V2, oWAT 27/27 V2) and the **7 of 56 discordant grades (31 same-visit, 25 cross-visit)**, recorded
with an explicit *"do not reconcile, deduplicate or average this away"* note.

## ⛔→✅ THE NESTED-GATE DEFECT: liver's gate WAS dead, and the conclusion SURVIVES the rebuild

`results_exploratory/EXPLORATORY_nested_gate_liveness.json`, job 21331199.

**`corr(joint, protein-only)` — prespecified threshold >0.99 = dead:**

| tissue | width ratio | liveness corr | gate |
|---|---:|---:|---|
| **liver** | ~5,377 : 10 = **710:1** | **0.999642** | ⛔ **DEAD** |
| oWAT | ~700:10 | 0.9990 | ⛔ dead (no consequence — arm already withdrawn) |
| **plasma** | **427 : 10 = 43:1** | **0.7703** | ✅ **LIVE** |

⛔ **Liver's +0.2413 was the non-nested difference wearing a nested label** — exactly the arithmetic I
flagged (0.2456 vs 0.2413, gap 0.0043). Ten descriptors among 5,377 proteins could not move a shared
penalty.

### ✅ But both live-gate rebuilds still exclude zero
| liver rebuild | increment | CI95 | liveness |
|---|---:|---|---:|
| original (dead gate) | +0.2413 | [+0.083, +0.420] | 0.9996 |
| **A — dimension-balanced** (10 protein PCs vs 10 descriptors) | **+0.2448** | [+0.118, +0.385] | 0.9161 |
| **B — block-scaled** (all proteins, each block ÷√dim) | **+0.2139** | [+0.076, +0.369] | 0.9582 |

⭐ **"Liver proteins add real signal over complexity" HOLDS.** ⛔ **But quote +0.2448 or +0.2139, never
+0.2413.**
⚠ The producer's own prediction was wrong in the *unflattering* direction: it expected the defect to
have **inflated** the result and for rebuild A to shrink toward zero. **It came in slightly larger.**

### ⛔ I WAS WRONG ABOUT PLASMA — the defect does NOT propagate by presumption
**Plasma's gate is LIVE at 0.7703.** ⭐ **The defect is WIDTH-RATIO dependent**: the narrow block
survives at 43:1 and dies at 710:1. **Asserting it for plasma would have been wrong**, and I had
started to.

⚠ Plasma still changes under block-scaling — increment **+0.1465 → −0.0068 [−0.214, +0.197]** — so
**plasma proteins add nothing over complexity once blocks are weighted comparably.** That *reinforces*
its "not established" verdict rather than altering it.

### ⭐⭐ The rule, and its free screen
✅ **Report `corr(joint, wide-only)` beside every nested increment.** One line, and it separates a
nested test from a non-nested difference wearing its label.
⭐ **The free screen when no liveness number is to hand: the nested-vs-non-nested GAP.** Liver's was
**0.0043**; plasma's was **0.0752**, an order of magnitude larger. **That gap alone flags a dead gate.**
⛔ **Do not assume the defect is universal.** It tracks the width ratio.

✅ The sealed `results/` was **not rewritten** — the superseded number stays in
`per_tissue_arms.json` with the correction in a separate file, **so the change is visible rather than
silent.**

### ✅ C2 (RNA↔protein concordance) DOES NOT have the defect — verified from source

The loop in `src/crossmodal_concordance.py` permutes the **outcome labels** — fibrosis within cohort on
the RNA side, Kleiner on the protein side — and **recomputes both association vectors from scratch
every draw.** That is the strict N-label null. **Feature identity is never permuted anywhere in C2.**

**Null audit run explicitly so the size of the avoided error is visible** (`results/c2_null_audit.json`):

| null | sd | max abs | p |
|---|---:|---:|---:|
| N-identity (the defective one) | 0.0160 | 0.0542 | 0.0005 |
| **N-label (what C2 actually used)** | **0.0782** | **0.3669** | 0.0005 |

⭐ **sd ratio 4.88x — the same defect ratio as C3.** The identity null *would* have been badly too
narrow here too; it simply would not have changed the answer.

✅ **Observed 0.4523 against a label-null MAXIMUM of 0.3669 — margin +0.0855, and 5.8 label-null sd from
zero.** Nothing in 2,000 label permutations came close. **The number is safe to quote.**

⭐ **Why C2 was right and C3's spec was wrong, in one line**: *"C2's null was mine and I wrote it to
regenerate the statistic end to end; C3's was specified from outside as a shortcut on the finished
vectors."* **That is the whole difference**, and it is why the standing rule generalises.

⚠ **The statistic is safe; the INTERPRETATION is unchanged and still narrow.** The specificity control
still says **general severity, not fibrosis**: fibrosis 0.4523 vs inflammation 0.4272, ballooning
0.4020, NAS 0.3913, steatosis 0.3586 — fibrosis wins by 0.025 with non-collinear labels.
✅ Sayable: *liver RNA and liver protein share a severity direction.* ⛔ Not sayable: *the modalities
agree about fibrosis.*
⚠ C2 also remains **demoted to exploratory** on its own rule (PB1 and PB3 wrong), which constrains how
hard it can be leaned on regardless of the null.

⭐ **Width-ratio confirmation across two lanes**: the gate died at **1,150:1** in the concordance lane
and **710:1** in the four-tissue lane, and **survived at 43:1** in plasma.
