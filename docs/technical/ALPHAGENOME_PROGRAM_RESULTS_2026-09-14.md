# AlphaGenome–MASLD regulatory response program — results

**Release state: candidate research layer. Nothing here is adopted. No Resource number, evidence class,
gene membership, threshold, hero gene or figure changes.** Every row of the accompanying Gene Catalog
attachment is marked `candidate_not_adopted`. The program is a scientific test of what sequence models
and the available liver measurements can and cannot say about allele effects in adult liver; it is not a
manuscript result and no panel in the paper is built from it.

Every number below is read from a deposited table under
`GWAS/finemapping/results/alphagenome_program/` and the file is named on the same line or in the section
header. Where a package's own independent checker marked a claim FAILED, that is said in section 9.
The specification these packages execute is
`scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md`, written before any analysis ran.

Integration package: `g-integration-20260915T072500Z`.

---

## 1. The three findings, stated plainly

Two of the three contradict what was written down before the work ran.

**(a) A frozen open-weight sequence model predicts signed endogenous liver caQTL allele effects far
better than the reporter task suggested.** On 32,322 Currin liver caQTL peak-lead variants in 2,643
1-Mb blocks, a frozen HyenaDNA 4,096-bp delta ridge reaches macro Spearman **0.1737
[0.1633, 0.1844]** against the measured ALT-allele effect on held-out chromosome folds, and beats the
allele-identity ridge control by **0.1301 [0.1162, 0.1440]**, with 10,000 of 10,000 paired block-bootstrap
draws positive
(`c2-endogenous-head-20260914T194000Z/tables/head_summary.tsv`). The per-seed count beside it is
**3 of 3 comparable seeds, not 5 of 5**: the Fisher z-macro of the allele-identity control is undefined at
seeds 1103 and 4721, so only three seeds admit a per-seed gain at all (`head_summary.tsv`,
`positive_gain_seeds` 3 and `comparable_seeds` 3). On the pooled out-of-fold metric, which is defined at
every seed, all five seeds beat the control. The written prediction C2.1 was that this
gain would be **below 0.082 with an interval covering zero**. It is neither. The failure does not depend on
picking the better of two backbones: the Caduceus delta ridge gain is 0.0841 [0.0709, 0.0975], also above
0.082. **The learning curve has not plateaued**: the final doubling of training loci, 12,500 to 25,858 per
fold, still delivered **+0.0173** macro Spearman with 5 of 5 seeds positive
(`.../tables/learning_curve.tsv`), and the deposit fits **+0.023 macro Spearman per doubling** over the top
three points (`c2-endogenous-head-20260914T194000Z/RESULTS.md` section 5; the raw points are in the table).

**(b) Reporter supervision adds nothing to that endogenous prediction on the route tested, and the
measured labels themselves agree only close to a peak.** Adding the out-of-fold prediction of a
reporter-trained head to a head trained on the caQTL label changes signed Spearman by **−0.0015
[−0.0070, +0.0039]**, bootstrap p 0.597, on 1,244 bridge variants in 168 1-Mb blocks
(`e2-transfer-20260914T203752Z/tables/defined_contrast_summary.tsv`). The bound that matters is
model-free: the two **measured** labels on the same variants agree at Spearman **0.2644
[0.1069, 0.4320]** when the variant is within 500 bp of the tested peak centre (209 variants, 89 blocks,
BH q 0.012 in a family of 12) and at **0.0517 [−0.0146, 0.1120]** in the prespecified 5-kb stratum, which
does not survive correction (`.../tables/measured_label_ceiling.tsv`,
`.../tables/measured_label_ceiling_bh.tsv`). No model can transfer a relation the measurements do not
share. The same package's positive control, in which both labels are generated from one shared learnable
component, tops out over its three nested arms at a gain of **+0.0142 [+0.0067, +0.0246]** —
`plus_reporter_oof_free_selection` at noise 4.0, against `plus_reporter_oof` +0.0122 [+0.0057, +0.0215] and
`plus_reporter_oof_locked_hparams` +0.0100 [+0.0054, +0.0169] at the same noise level
(`.../tables/positive_control.tsv`). The maximum point estimate stays below the 0.02 threshold E2.1 was
written against, though its upper bound 0.0246 does not, so this route could not have produced a large gain
even with perfectly shared labels.

**(c) The project-trained adult hepatocyte ChromBPNet does not beat the archived AlphaGenome Atlas
accessibility channel on either endpoint, and the gap is an upper bound.** At Currin caQTL peak leads
(1,236 matched variants, 915 blocks) ChromBPNet reaches signed Spearman 0.6565 [0.6243, 0.6857] against
the Atlas liver ATAC quantile's 0.7258 [0.6995, 0.7493]; the paired contrast is **−0.0693
[−0.0913, −0.0482]**, bootstrap p floored at 1e-4
(`c1-matched-endogenous-20260914T193517Z/tables/c1_contrasts.tsv`,
`.../tables/c1_model_metrics.tsv`). On the enlarged matched set the Borzoi fold-in allows (3,845 leads,
1,743 blocks, the archived AlphaGenome cap reached exactly) it is **−0.0910 [−0.1055, −0.0768]**, so the
enlargement widens the gap rather than closing it; the same loader reproduces −0.0693 on the 1,236-lead
subset in the same run, so that move is the enlargement and not a recipe change. At the 495
allelic-imbalance sites in 54 blocks the contrast is **−0.1481 [−0.2361, −0.0532]**, p 0.0036
precision-weighted, and −0.1476 [−0.2250, −0.0524], p 0.0046 unweighted
(`c1-endpoint2-v2-20260915T084559Z/tables/c1e2_contrasts.tsv`). Prediction C1.1 said the opposite sign.
Borzoi, now scored at both endpoints, sits between them at endpoint 1 (0.6791) and with ChromBPNet at
endpoint 2 (0.2722 weighted), not with AlphaGenome.
**Two comparator asymmetries both favour AlphaGenome**: the ChromBPNet column is a single held-out-fold
model's prediction, while the AlphaGenome column is a mean over three liver ATAC tracks of a hosted static
model whose training folds are unknown here and whose published corpus includes ENCODE adult liver ATAC and
DNase. Nothing in either package holds anything out of AlphaGenome, and track averaging removes variance
the single-fold column keeps. The supportable statement is that ChromBPNet does not beat AlphaGenome under
these conditions, not that AlphaGenome is the better sequence model.

---

## 2. The arm that never ran, and the one action that unblocks it

**AlphaGenome adaptation did not start.** The released weights `google/alphagenome-all-folds` are a gated
Hugging Face repository and this account is not granted access; the probe on 2026-09-14 returned
`GatedRepoError 403`. No local JAX or Haiku runtime and no `alphagenome_research` checkout exist. This is
recorded in `scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md` section 1 and in the rights
register `a2-rights-exposure-20260914T193626Z/tables/rights_exposure.tsv`, where the row reads
"adaptation permitted under the licence, but access not granted, 403".

**The single action that unblocks it:** accept the model terms on the `google/alphagenome-all-folds` model
page using the Hugging Face account whose token this machine holds. Everything else the arm needs is a
compute request, not a permission.

Until then every AlphaGenome quantity in this program is either an archived Atlas score or a zero-shot
hosted-API call, and the hosted API's terms prohibit using its outputs to train other models
(`a2-.../tables/rights_exposure.tsv`, `training_allowed = PROHIBITED — hosted Atlas outputs, not for
training other ML models`). No AlphaGenome output was used as a training input, target or learned
calibration anywhere in this program.

---

## 3. What ran

| Package | Deposit | What it delivers |
|---|---|---|
| A1 | `a1-mpra-reliability-v2-20260914T204500Z` | reporter measurement reliability, 4,359 elements, 239 blocks |
| A2 | `a2-rights-exposure-20260914T193626Z` | rights, exposure and defect register, 165 artifact rows |
| B1/B2 | `b1-response-layer-v2-20260914T204700Z` | locus-by-context-by-readout layer, 989 signals, 10,511 variant rows |
| B3 | `b3-saturation-20260914T193641Z` | bounded saturation sensitivity, 3,000 regions, 1,337 blocks |
| C1 | `c1-matched-endogenous-20260914T193517Z` | matched model comparison, endpoint 1 |
| C1 e2 | `c1-endpoint2-gpu-20260914T203750Z` | the GPU scoring: ChromBPNet at 500 sites, Borzoi at 3,102 variants |
| C1 e2 v2 | `c1-endpoint2-v2-20260915T084559Z` | the four-family endpoint 2 (495 sites) and the enlarged endpoint 1 (3,845 leads) |
| C2/C3 | `c2-endogenous-head-20260914T194000Z` | trainable open head on caQTL, plus learning curve |
| D1 | `d1-lineage-allelic-20260914T194213Z` | lineage-resolved allelic labels and a mapping-bias control |
| E1 | `e1-bridge-20260914T194145Z` | reporter-to-endogenous bridge census |
| E2 | `e2-transfer-20260914T203752Z` | transfer contrast, third Holm member |
| E3 | `e3-reporter-context-20260914T193515Z` | hosted-model reporter scoring, two padding arms |
| F1 | `f1-phase-20260914T193608Z` | phase substrate for local haplotypes |
| F2/F3 | `f2-haplotype-v2-20260915T063000Z` | GNMT verdict on a matched readout, and the corrected null |
| B1 refresh | `b1-dbsnp-refresh-20260915T084700Z` | the B1 indel-orientation columns recomputed against the current dbSNP deposit |

Superseded deposits that must not be read for numbers: `a1-mpra-reliability-20260914T193627Z`,
`b1-response-layer-20260914T{193909,194738,195729}Z`, `f2-haplotype-20260914T230433Z`. Each carries its own
supersession marker. The indel-orientation columns of `b1-response-layer-v2-20260914T204700Z` are
superseded on their own, and only on the 548 indel rows; the rest of that deposit stands. Read those
columns from `b1-dbsnp-refresh-20260915T084700Z/tables/research_layer_indel_orientation_refreshed.tsv`, as
`SUPERSEDED_COLUMNS_indel_orientation.md` in the layer directory says. The endpoint-1 and endpoint-2
statistics of `c1-endpoint2-gpu-20260914T203750Z` are likewise superseded on their own, by
`c1-endpoint2-v2-20260915T084559Z`, which reran them with the Borzoi column present; that deposit's GPU
scores, defect findings and identity checks stand and are the v2 package's inputs.

Resampling unit throughout: the 1-Mb block, except for the reporter arms, which resample the 239 Borzoi
long-range blocks. Donors, not cells or reads, are the biological unit wherever a donor-level statement is
made. Bootstrap seed 20260914; comparison seeds 1103, 2909, 4721, 6673, 8111.

---

## 4. Endogenous allele effects

### 4.1 Matched model comparison, endpoint 1 (Currin caQTL peak leads)

`c1-endpoint2-v2-20260915T084559Z/tables/c1e2_model_metrics.tsv`. Signed Spearman against `beta_alt`,
10,000 block-bootstrap draws, every model recomputed inside the same draw. **Tier A4 (the four spec
families plus two controls) is n = 3,845 in 1,743 blocks** after the Borzoi GPU jobs landed and Borzoi's
lead coverage rose from 1,236 to 3,878; tier A, which also requires ChromBPNet HepG2, is n = 1,601 in 1,088
blocks. The 1,236-lead readings the wave-1 C1 package published are the tier-A4 wave-1 subset column below,
reproduced by the v2 loader to four decimals.

| model | tier A4, n = 3,845 | 95% CI | tier A, n = 1,601 | 95% CI | wave-1 subset, n = 1,236 |
|---|---|---|---|---|---|
| alphagenome_atac_liver | **0.7361** | 0.7206–0.7507 | 0.7334 | 0.7110–0.7542 | 0.7258 |
| alphagenome_dnase_liver | 0.7120 | 0.6961–0.7275 | 0.7170 | 0.6944–0.7378 | 0.7108 |
| borzoi_atac | 0.6791 | 0.6630–0.6950 | 0.6920 | 0.6681–0.7143 | 0.6831 |
| chrombpnet_adult_hep_gse281367 | 0.6504 | 0.6315–0.6695 | 0.6611 | 0.6331–0.6872 | 0.6539 |
| chrombpnet_adult_hep | 0.6450 | 0.6258–0.6643 | 0.6570 | 0.6296–0.6830 | 0.6565 |
| borzoi_dnase | 0.6079 | 0.5891–0.6261 | 0.6294 | 0.6022–0.6549 | 0.6207 |
| chrombpnet_hepg2 | not in tier A4 | | 0.6109 | 0.5787–0.6416 | 0.6133 |
| control_allele_identity | −0.0008 | −0.0319–0.0305 | −0.0058 | −0.0540–0.0411 | −0.0036 |
| control_position | −0.0183 | −0.0490–0.0130 | −0.0284 | −0.0781–0.0214 | −0.0260 |

**The C1.1 contrast sign survives the enlargement and the gap widens.** ChromBPNet adult minus AlphaGenome
ATAC is **−0.0910 [−0.1055, −0.0768]** on the 3,845, against −0.0693 [−0.0913, −0.0482] on the wave-1 1,236
and −0.0764 [−0.0961, −0.0572] on tier A. Both bootstrap p values are floored at 1/10,000 draws. The
component-matched comparison is deliberate: the wave-1 subset is recomputed by the same loader in the same
run, so the move from −0.0693 to −0.0910 is the enlargement, not a recipe change.

**Matching still costs 88 percent of the label universe.** The Currin lead universe is 32,336 leads in
2,643 blocks (`.../tables/c1e2_matched_counts.tsv`); tier A4 is capped at 3,845 by the archived AlphaGenome
panel, which the Borzoi run now meets exactly. The matched set is the strong-effect tail: the same
ChromBPNet reaches 0.4495030608 on all 32,336 leads, reproducing the deposited
`currin_external_calibration.json` `point_estimates.direction.spearman` to 5.55e-17
(`c1-endpoint2-v2-.../tables/qc_addendum.json`). Absolute values here are not comparable to that; only
within-set contrasts are.

### 4.2 Matched model comparison, endpoint 2 (measured allelic imbalance)

`c1-endpoint2-v2-20260915T084559Z/tables/c1e2_model_metrics.tsv`, 495 sites in 54 1-Mb blocks,
precision-weighted by 1/se, unweighted beside it. **This is now the four-family comparison the
specification describes**: the Borzoi GPU jobs completed on 2026-09-14 at 21:35 and 21:40, after the
c1-endpoint2-gpu package had written its report, and Borzoi covers all 500 sites. The matched set is the
same 495 in 54 blocks, because the binding constraint is still the 5 sites the archived AlphaGenome route
never served, so every pre-existing row below is unchanged to four decimals.

| model | weighted Spearman | 95% CI | unweighted | 95% CI |
|---|---|---|---|---|
| alphagenome_atac_liver | 0.4507 | 0.3799–0.5393 | 0.4439 | 0.3822–0.5198 |
| alphagenome_dnase_liver | 0.4194 | 0.3365–0.5215 | 0.4121 | 0.3460–0.4984 |
| chrombpnet_adult_hep | 0.3025 | 0.2056–0.4172 | 0.2963 | 0.2107–0.4091 |
| chrombpnet_adult_hep_gse281367 | 0.2959 | 0.1886–0.3958 | 0.2424 | 0.1456–0.3368 |
| **borzoi_atac** | **0.2722** | **0.1823–0.3809** | **0.2874** | **0.2108–0.3768** |
| **borzoi_dnase** | **0.2235** | **0.1314–0.3475** | **0.2411** | **0.1648–0.3426** |
| control_allele_identity | 0.1079 | 0.0284–0.1798 | 0.0554 | −0.0041–0.1151 |
| control_position | 0.0708 | −0.0451–0.1921 | 0.0443 | −0.0577–0.1576 |

**Borzoi lands with the ChromBPNet families, not with AlphaGenome, and the two open-weight families are not
separable here.** ChromBPNet adult minus Borzoi ATAC is +0.0304 weighted [−0.0864, +0.1480], p 0.58; both
sit below the archived AlphaGenome ATAC column by an interval excluding zero (ChromBPNet −0.1481 [−0.2361,
−0.0532], p 0.0036; Borzoi ATAC −0.1785 [−0.2456, −0.1083], p 0.0002). That is the same ordering endpoint 1
gives, now with a third family in it. A non-separation on 54 blocks is not an equivalence: the ChromBPNet-
minus-Borzoi interval is wide enough to contain a difference as large as the ChromBPNet-to-AlphaGenome gap
(`c1-endpoint2-v2-.../tables/c1e2_contrasts.tsv`).

**The allele-identity control is not at chance at these sites.** On all 500 sites in 59 blocks it reaches
**0.1023 weighted with an interval 0.0206 to 0.1749 that excludes zero**, and 0.0536 unweighted with an
interval −0.0045 to 0.1117 that includes zero. The same control sits at 0.000306 (−0.0104 to 0.0108) on the
32,336 Currin leads (`c1-matched-endogenous-.../tables/c1_unmatched_per_model.tsv`). A label-free function
of the two alleles should not predict a measured allele ratio; that it does here is the expected signature
of substitution-class-dependent mapping bias in uncorrected heterozygote calls, and it caps how much of the
sequence models' 0.22 to 0.45 at these sites reads as sequence skill.

**The Borzoi column reproduces the archived panel exactly.** 200 Currin leads the archived src/79 panel had
already scored were rescored by the program's own loader: maximum absolute difference 4.44e-16 on both the
ATAC and the DNase channel, Pearson 1.000000, sign agreement 200/200, orientation-multiplier agreement
200/200. The merge rule was written before the endpoints were read (merge only if max abs diff < 1e-6 on
both channels), so the enlarged endpoint-1 Borzoi column is one column, not two provenances
(`c1-endpoint2-v2-.../qc_checks.json`).

### 4.3 Trainable open head and the learning curve

`c2-endogenous-head-20260914T194000Z/tables/head_summary.tsv`, 32,322 unique Currin peak-lead variants in
2,643 blocks, five held-out chromosome folds, Fisher z-macro, five-seed ensemble, 10,000 block bootstraps
within fold.

| arm | macro Spearman | 95% CI | gain vs allele identity | gain 95% CI |
|---|---|---|---|---|
| **hyenadna / delta ridge** | **0.1737** | 0.1633–0.1844 | **0.1301** | 0.1162–0.1440 |
| hyenadna / full ridge | 0.1600 | 0.1494–0.1709 | 0.1165 | 0.1027–0.1306 |
| caduceus / delta ridge | 0.1276 | 0.1170–0.1386 | 0.0841 | 0.0709–0.0975 |
| caduceus / full ridge | 0.1071 | 0.0958–0.1184 | 0.0636 | 0.0508–0.0763 |
| allele-identity ridge (control) | 0.0435 | 0.0325–0.0543 | — | — |

**Label orientation, and why the Currin `EA` column does not threaten it.** The label is the source `beta`
unchanged, the FastQTL slope for ALT dosage, so `beta_alt > 0` means ALT raises accessibility. The lead
file's `EA` column disagrees with that reading on 18,981 of 35,361 rows, but it is **a recode of the sign of
beta, not an independent effect-allele label**: over all 35,361 rows `EA == ALT` in exactly the 16,380 rows
with `beta > 0` and `EA == REF` in exactly the 18,981 rows with `beta < 0`, with zero allele-set mismatches
and zero `beta == 0` (`c2-endogenous-head-20260914T194000Z/RESULTS.md` section 3). A column that is a
deterministic function of the sign of beta carries no orientation information beyond beta's own, so it
cannot be evidence that beta is EA-referenced. That **strengthens** the decision to keep beta
ALT-referenced rather than weakening it. The independent check is the one that matters and it agrees: on the
3,845 leads that also appear in `direction_features_caqtl.tsv`, whose effect comes from the Currin
per-chromosome nominal files, sign agreement with `beta_alt` is 0.9943 and the rank correlation is 0.9785;
had the lead beta been EA-oriented, agreement would sit near 0.5
(`c1-matched-endogenous-20260914T193517Z/RESULTS.md`).

Dropping the 5,191 strand-ambiguous variants moves the headline by 0.0013. The reference this is measured
against is the published reporter arm: HyenaDNA 4-kb delta ridge 0.21622 [0.16278, 0.27170], gain 0.08201
[0.02743, 0.13892] on 1,033 MPRA elements
(`Analysis/MASLD_Model_Benchmark/executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv`).
**The two gains are not statistically distinguishable** — the reference interval is wide and overlaps this
one. What is established is that the gain here is large, precise and nowhere near zero, which is what the
prediction denied.

Learning curve, HyenaDNA delta ridge, macro Spearman by training loci per fold
(`.../tables/learning_curve.tsv`; the macro is undefined at the two smallest sizes for some seeds, and the
pooled curve is tabulated beside it in the deposit):

| training loci per fold | 800 | 2,000 | 5,000 | 12,500 | 25,858 (all) |
|---|---|---|---|---|---|
| macro Spearman | undefined (0 of 5 seeds) | 0.0939 (**one seed**, 6673) | 0.1123 | 0.1473 | **0.1646** |
| pooled Spearman, all 5 seeds | 0.0351 | 0.0714 | 0.1100 | 0.1469 | 0.1644 |

The 2,000-loci entry is **not** a five-seed mean: the Fisher z-macro is undefined at four of the five seeds
there, so 0.0939 is seed 6673 alone. The pooled out-of-fold metric, which is defined at all five seeds, gives
**0.0714** at 2,000. The 5,000, 12,500 and 25,858 entries are five-seed means.

Final doubling **+0.0173**, 5 of 5 seeds positive (per-seed +0.0134 to +0.0213). No flattening anywhere in
the measured range.

**Two things this curve is not.** It varies training loci and nothing else, so it says nothing about
donors directly. And the selected ridge alpha sits at the grid maximum 1e4 in 113 of 125 primary fits, so
every value here is a lower bound on what this feature set could reach with heavier shrinkage available.

---

## 5. Reporter measurement, transfer and construct scoring

### 5.1 A1 reporter reliability

`a1-mpra-reliability-v2-20260914T204500Z/tables/reliability_by_context.tsv`, 4,359 paired-SNV elements in
239 long-range blocks, raw two-versus-two Pearson averaged over the three balanced partitions, 10,000
block-bootstrap draws.

| quantity | context | raw 2v2 Pearson | 95% CI | projected 4-replicate | 95% CI |
|---|---|---|---|---|---|
| activity | HepG2_control | 0.9840 | 0.9810–0.9869 | 0.9919 | 0.9904–0.9934 |
| activity | HepG2_PAOA | 0.9812 | 0.9775–0.9846 | 0.9905 | 0.9886–0.9923 |
| activity | LX2_control | 0.9855 | 0.9824–0.9881 | 0.9927 | 0.9911–0.9940 |
| activity | LX2_TGFb | 0.9906 | 0.9888–0.9922 | 0.9953 | 0.9944–0.9961 |
| allele effect | HepG2_control | 0.7733 | 0.7393–0.8030 | 0.8722 | 0.8501–0.8907 |
| allele effect | HepG2_PAOA | **0.6938** | 0.6356–0.7383 | 0.8192 | 0.7772–0.8495 |
| allele effect | LX2_control | 0.8481 | 0.8314–0.8615 | 0.9178 | 0.9079–0.9256 |
| allele effect | LX2_TGFb | 0.8592 | 0.8400–0.8735 | 0.9243 | 0.9130–0.9325 |
| treatment minus control allele effect | HepG2 | **0.0218** | −0.0816–0.1230 | 0.0427 | −0.1777–0.2190 |
| treatment minus control allele effect | LX2 | 0.5931 | 0.5485–0.6288 | 0.7446 | 0.7084–0.7721 |

Stated as ranges over the four contexts: activity reliability runs from **0.9812** (HepG2_PAOA, the floor)
to 0.9906 (LX2_TGFb); allele-effect reliability from **0.6938** (HepG2_PAOA) to 0.8592 (LX2_TGFb) raw, and
from 0.8192 (HepG2_PAOA) to 0.9243 projected. HepG2_PAOA holds the floor of all three ranges.

**The section-3 decision rule, fixed before the run, fires on HepG2.** Its projected four-replicate
reliability of the treatment-minus-control allele effect is 0.0427 [−0.1777, 0.2190], below the 0.50
threshold, so **HepG2's allele-by-state reporter endpoint is declared not measurable on this dataset**
(`.../tables/decision_rule_section3.tsv`). That is a measurement bound, not a model failure. LX2 is
measurable at 0.7446 [0.7084, 0.7721]. Depth does not explain the HepG2 failure: its top DNA-depth tertile
reaches 0.0612 against LX2's bottom tertile 0.5137, and 5 of 6 HepG2 tertile intervals include zero.

### 5.2 E1 bridge census

`e1-bridge-20260914T194145Z/bridge_counts.tsv`, 4,359 reconstructed MPRA SNVs joined to nine endogenous
liver labels on hg38 chromosome, position and both allele orders.

| label | MPRA variants matched | 1-Mb blocks | share of the endogenous side |
|---|---|---|---|
| Currin nominal caQTL | 3,894 | 243 | 3,894 / 7,133,527 |
| GTEx v8 liver sQTL | 772 | 31 | 772 / 159,014 |
| Broadaway liver eQTL | 202 | 30 | 202 / 5,340 |
| Currin LD proxies (r2 ≥ 0.8) | 71 | 46 | 71 / 36,719 |
| Currin q<0.05 leads | 69 | 50 | 69 / 34,098 |
| allelic imbalance, union | **16** | **9** | 16 / 714 site rows |

**"Carries a caQTL test" is a weak statement.** Currin tested every variant within 1 Mb of a peak, so the
median MPRA variant is 11,336 bp from the centre of its nearest tested peak. Within 500 bp there are only
**209 variants in 89 blocks** (`.../e2_distance_strata.tsv`). Only the Currin nominal caQTL label can carry
a block-resampled transfer test at all.

### 5.3 E2 transfer

`e2-transfer-20260914T203752Z/tables/defined_contrast_summary.tsv`, nearest-peak label, primary stratum
5,000 bp, 1,244 variants in 168 blocks, macro over five held-out chromosome folds, 10,000 block bootstraps.

| stratum | n | blocks | arm A (endogenous only) | arm B (plus reporter OOF) | gain | 95% CI | p |
|---|---|---|---|---|---|---|---|
| 500 bp | 209 | 89 | +0.0978 | +0.1008 | +0.0030 | −0.0246–+0.0319 | 0.849 |
| 1 kb | 350 | 113 | +0.0523 | +0.0446 | −0.0077 | −0.0233–+0.0094 | 0.371 |
| **5 kb (primary)** | **1,244** | **168** | **−0.0191** | **−0.0207** | **−0.0015** | **−0.0070–+0.0039** | **0.597** |
| 10 kb | 1,838 | 186 | −0.0071 | −0.0089 | −0.0018 | −0.0059–+0.0030 | 0.488 |
| 100 kb | 3,692 | 239 | −0.0321 | −0.0324 | −0.0002 | −0.0028–+0.0026 | 0.901 |
| unlimited | 3,893 | 243 | −0.0295 | −0.0295 | −0.0000 | −0.0026–+0.0026 | 0.993 |

Both arms are at chance on the endogenous label in the primary stratum. The reporter head supplying arm
B's extra column is working at the published level: macro Spearman **0.2178** on its own MPRA label over
4,358 elements under chromosome-group folds, against the published 0.21622 on 1,033 elements
(`.../tables/reporter_head_self_evaluation.tsv`). Arm B fails to gain anything with a live reporter head,
not a dead one. Fold leakage is excluded by construction: zero shared elements and zero shared blocks
between the two training sets in every evaluation fold (`.../tables/fold_membership.tsv`).

### 5.4 E3 what the hosted model says about the construct

`e3-reporter-context-20260914T193515Z/tables/e3_correlations.tsv`, 296 oligos scored at 16,384 bp in two
padding arms, 1,184 API calls, 112 long-range blocks, Spearman against the measured HepG2-control allele
effect, 10,000 block bootstraps.

| stratum | arm | DNASE Spearman | 95% CI |
|---|---|---|---|
| random (n=200, 93 blocks) | native flanks | +0.144 | +0.013–+0.287 |
| random | N padding | +0.153 | +0.016–+0.302 |
| random | native minus N-pad (paired) | −0.009 | −0.070–+0.052 |
| top-effect (n=100, 57 blocks) | native flanks | +0.292 | +0.072–+0.532 |
| top-effect | native minus N-pad (paired) | +0.022 | −0.106–+0.160 |

**Native genomic flanks add nothing to the accessibility prediction of a reporter oligo**, and on the RNA
channel they actively hurt: native minus N-pad is **−0.389 [−0.585, −0.178]** in the top-effect stratum,
where the N-padded RNA readout reaches +0.346 [+0.173, +0.525] (BH q 0.0016 in the prespecified 48-cell
family) and its native counterpart is −0.042 [−0.260, +0.201]. In a native window the local RNA prediction
is dominated by whatever transcription the locus carries and the 126-bp edit is a rounding error against
it. The top-effect stratum is selected on label magnitude (allele-effect sd 1.564 against 0.385 in the
random stratum) and is not an unbiased estimate of anything.

---

## 6. Lineage labels, saturation, haplotypes and the locus layer

### 6.1 D1 lineage-resolved allelic labels

`d1-lineage-allelic-20260914T194213Z/tables/lineage_feasibility.tsv`, 30 donors across GSE281367 and
GSE244832, 1,119 target sites, reliability filter ≥10 reads, ≥3 per allele, alt fraction in [0.15, 0.85],
≥3 het donors per site.

| lineage | sites passing | 1-Mb blocks | het donor-site observations | median reads per het donor-site |
|---|---|---|---|---|
| ALL_READS | 611 | 60 | 4,480 | 37 |
| ALL_LABELLED | 398 | 55 | 3,118 | 39 |
| **Hepatocyte** | **292** | **50** | 2,173 | 39 |
| Endothelial | 64 | 23 | 505 | 23 |
| Stellate_Cell | 27 | 15 | 255 | 23 |
| Plasma_Cell | 17 | 11 | 154 | 18 |
| Macrophage | 9 | 7 | 75 | 17 |
| Cholangiocyte | 3 | 3 | 66 | 18 |
| Kupffer_Cell / LSEC | 1 / 1 | 1 / 1 | 22 / 22 | 17.5 / 14.5 |
| NK_T_Cell / B_Cell | 0 / 0 | 0 / 0 | 13 / 1 | — |

**Only hepatocyte supports a lineage-specific allele-effect label; endothelium is at the edge.** The
binding constraint is depth per lineage at these sites, not donor count: every lineage has reads in all 30
donors, but the median reads per donor-site outside hepatocyte and endothelium is 1 to 3.

The WASP-style mapping-bias control found **no allele-specific excess above a matched single-base control**:
allele-swap minus control read-failure rate is **−0.0138** in GSE281367 (bwa) and **−0.2041** in GSE244832
(bowtie2) (`.../tables/mapping_bias.json`). The prespecified primary filter removes 80 and 98 percent of
sites and measures short-read alignment fragility at 50–51 bp, not allele-specific bias; it is reported and
not used. The matched filter removes 52 of 435 and 24 of 279 reliable sites and moves every reference-bias
mean toward zero, so part of the residual reference bias is removable mapping bias and part is not
(GSE244832 labelled-cell residual after the filter is still −0.086).

### 6.2 B1/B2 locus-by-context-by-readout layer

`b1-response-layer-v2-20260914T204700Z/response_layer_signals.tsv`,
`.../response_layer_variants.tsv.gz`, `.../b1_checks.json`. 989 eligible signals (A_direct 29, B_direct 173,
C_enzyme 787), 10,511 (signal, variant) rows. **Two different denominators, kept apart:** **3,822 of the
10,511 rows** are readable, a row share of **0.364**, and **623 of the 989 signals** pass the Atlas coverage
gate, a signal share of 0.630. The two are separate counts of separate things. The signal gate is
all-or-nothing on a credible set, so the 623 gated-readable signals are exactly the 623 signals whose
`atlas_coverage_state` is `covered`; the other 366 are 353 `partial_indel` and 13 `uncovered_indel`
(`.../coverage_gate.tsv`). Every readable row sits inside one of the 623.

| measured channel | rows | signals | 1-Mb blocks |
|---|---|---|---|
| in a donor snATAC peak | 291 | 194 | 102 |
| measured allelic imbalance | 190 | 118 | 33 |
| liver eQTL direction | 7,342 | 773 | 315 |
| differential-accessibility state | 4 | 3 | 3 |
| model-API indel rescue (a prediction) | 222 | 176 | 103 |

**734 of 989 signals distinguish no explanation class from measurement.** 45 distinguish local regulatory
alone, 134 coding alone, 76 both, none three or four. No signal in this universe can distinguish RNA
processing or a candidate long-range explanation, because the deposit holds no sQTL, no contact map and no
enhancer-gene assay in adult liver at these loci. That is a bound on what the layer can test, not a
negative result about those mechanisms.

**The Atlas liver panel is still not an adult panel.** Read from the build's own scorer metadata
(`.../track_panel_composition.tsv`): ATAC 3 adult of 3; DNase 2 adult and 1 embryonic of 3; H3K27ac 2 adult
of 2; RNA and splice-site usage 3 adult, 1 child and 1 embryonic of 5. The Track 0 direction columns run on
a different 7-track panel that also pools two embryonic hepatocyte tracks. Restricting the RNA panel from 7
tracks to 5 and then to 3 adult tracks flips the predicted sign on 41,389 of 877,512 rows (4.72%) and then
32,528 (3.71%), and moves the posterior-weighted eQTL-direction concordance from 0.5623 to 0.5620 to
0.5665. **All three arms sit at chance per row** (0.5004–0.5006); no interval, null or resampling is
computed, and these are LD-correlated variants inside signals, not independent tests.

### 6.3 B3 bounded saturation sensitivity

`b3-saturation-20260914T193641Z/tables/b3_contrasts.tsv`, 3,000 archived U1 H3K27ac regions in 1,337
1-Mb blocks, 10,000 block bootstraps, BH inside four named families. **The frame is chr1–chr12 complete
plus chr13 at 63%**, not a complete chr1–chr13: at the snapshot chr1 through chr12 were finished region for
region and chr13 held 1,558 of 2,474 regions (62.97%), with chr14 onward at zero
(`.../tables/frame_by_chromosome.tsv`). Every number here and the block bootstrap under it are statements
about that frame.

B3.1 (top-skill minus bottom-skill quintile, `mean_abs`, native scale, adult liver) is met on its verdict
channel and is small: CHIP_HISTONE **+0.00166 [+0.00014, +0.00317]**, BH q 0.037, a 6.7 percent relative
difference carried by a Q1 step rather than a gradient. DNASE +0.00843 [+0.00448, +0.01240] and ATAC
+0.00656 [+0.00243, +0.01066] are larger, CAGE is null, and **AVI_SCORE runs the other way at −0.01432
[−0.02287, −0.00569]**, BH q 0.0013.

B3.2 fails in direction, not in power. In adult-liver accessibility and histone channels, promoter status
is associated with **less** predicted base sensitivity and sensitivity **rises** with distance to the
nearest TSS: DNASE promoter −0.811 [−0.881, −0.743] SD, `log10_dist_tss` +0.219 [+0.187, +0.249];
CHIP_HISTONE −0.704 [−0.777, −0.629] and +0.189 [+0.155, +0.221]. The predicted direction holds only for
CAGE (+0.760, −0.319) and the model's own aggregate score. All twelve intervals exclude zero at BH
q ≤ 3.1e-4.

A post-hoc adjustment the hashed prespecification omitted removes about 70 percent of the skill
association: adding measured H3K27ac `signal_mean` and `signal_sd` takes the standardised skill slope on
positional concentration from +0.200 to +0.054 for CHIP_HISTONE and to −0.000 for AVI_SCORE
(`.../tables/b3_posthoc_signal_adjustment.tsv`). Labelled post hoc everywhere because it was not in the
hashed file.

### 6.4 F1 phase substrate

`f1-phase-20260914T193608Z/tables/f1_summary.json`,
`.../tables/f1_block_bootstrap_intervals.tsv`. The local 1000 Genomes EUR PLINK panel retains explicit
phase on 22 of 22 chromosomes, 379 samples, 758 haplotypes, 0 unphased calls in 506,344 exported genotype
calls. It is Phase 1 low-coverage, statistically phased and imputed, so its haplotypes are model-based, not
read-backed.

**1,508 of 1,587 P5 pairs (0.9502 [0.9266, 0.9702], 250 blocks) have D prime above 0.9**, and the
consequence is the point: for 1,272 of 1,516 panel-computable pairs (0.8391 [0.7974, 0.8789]) the rarer of
the two single-variant arms has EUR haplotype frequency below 0.01, and 668 of 1,516 (0.4406
[0.3913, 0.4910]) have an arm seen in **0 of 758** EUR haplotypes. A P5 residual is therefore a legitimate
in-silico perturbation but not a statement about a combination of alleles that people carry. These 668
zeros are direct counts of phased hardcalls, the same route that gives 2 of 758 at the GNMT pair in section
6.5, so they are not the EM artefact that produced F2's zero there.

Read-backed phase from existing BAMs exists for a handful of loci at most: 9 site pairs per deposit have
both sites measured, 3 of them at ≤500 bp, and the strongest candidate is PPP1R3B at 238 bp with 4
double-het donors in GSE281367 (`.../tables/f1_read_backed_phase_summary.tsv`).

### 6.5 F2/F3 the GNMT pair and the corrected additivity null

`f2-haplotype-v2-20260915T063000Z/tables/f2v2_gnmt_verdict.json`,
`.../tables/f2v2_summary.json`. Every number is a model prediction about substitutions in one GRCh38
reference sequence scored through the hosted model API. Nothing here measures an interaction in cells or
people.

The GNMT pair rs2296805 × rs2296804 is **additive on this scale in all four channels** against a
readout-matched 100-draw chr6 null whose draws sit inside real chr6 protein-coding gene spans: residual
inside the central 95 percent with p **0.1188** (RNA), **0.2772** (ATAC), **0.4059** (DNase), **0.5743**
(H3K27ac); p floor 1/101 = 0.0099.

**Two of the four scored sequences are haplotypes almost nobody carries: each single-variant arm is on
2 of 758 European haplotypes, 0.264%.** F1 and F2 disagreed on this count and the disagreement is now
settled by recomputation from `data/1kg_eur/chr6_eur` directly
(`g-integration-20260915T072500Z/tables/gnmt_arm_frequency_resolution.tsv`, job 21766819). Two independent
phased routes — plink2 `--export vcf` and plink2 `--export haps`, both reading the `.pgen` — agree exactly:
REF1-REF2 303, **ALT1-REF2 2**, **REF1-ALT2 2**, ALT1-ALT2 451 of 758 haplotypes, panel r2 0.9782, D prime
0.9890. **F1's 2 of 758 is right and F2's 0.0000 is an artefact of the file format it read.** PLINK 1.9
`--ld` takes the `.bed`, which carries no phase, so it estimates haplotype frequencies by EM from
genotypes; of the 179 double-heterozygous individuals 177 are in coupling and 2 in repulsion, and EM
resolves all 179 to the majority phase. That moves exactly those two repulsion chromosomes into the two
common haplotypes: EM gives REF1-REF2 = 0.402375 × 758 = 305.00 = 303 + 2 and ALT1-ALT2 = 0.597625 × 758 =
453.00 = 451 + 2, and reports r2 as 1.0000 instead of 0.9782. The three candidate explanations that would
have been more benign are all excluded on the same panel: both routes use all **379 founders and 0
nonfounders**; there are **0 missing calls** at either variant and 0 haplotypes dropped; and PLINK printed a
literal `0` while printing `0.402375` on the line above, so 2/758 = 0.002639 was not lost to rounding. The
conclusion the claim rests on is unchanged and now has an exact denominator: a P5-style single-variant arm
at this locus is carried by 2 of 758 European haplotypes. One caveat survives the recomputation — this
panel is 1000 Genomes Phase 1 EUR, statistically phased and imputed, so the 2 is a phasing model's output,
not read-backed phase.

For the 199-pair family against the corrected 2,500-draw null, BH at q<0.10 rejects **3 RNA, 6 ATAC, 18
DNase and 34 H3K27ac rows**. Read those as a bound on the test, not as predicted non-additivity:

- Every rejection sits below 2 kb of separation, where the two ±1 kb chromatin readout windows overlap; in
  the ≥2 kb stratum rejections are 0 in all four channels.
- Arm-magnitude matching drops them to **0 / 0 / 4 / 0**.
- The 199 rows hold only 131 distinct variant pairs; as distinct pairs the rejections are 1 / 3 / 10 / 16.
  The three RNA rejections are **one variant pair entering the family three times**, once per COLOC signal
  naming CDK18, and both of its variants lie 12,126 bp outside the readout window.
- The RNA readout is unmatched on width pair by pair: 142 observed pairs read out over a gene span of 855
  to 687,828 bp while almost every null draw reads out over 20,000 bp.

**The existing P5 deposit's null could not have rejected below 53 tied pairs.** Its 300 draws floor the
empirical p at 0.003322, and BH at q=0.10 over m=1,572 needs k ≥ 53 pairs tied at that floor before
anything is rejected; observed ties are 2 to 16 per channel. Its reported absence of non-additivity is a
statement about the null's resolution, not about additivity
(`a2-rights-exposure-20260914T193626Z/DEFECTS.md` item e).

---

## 7. The confirmatory family

Three contrasts were named as confirmatory before any of them ran; everything else in the program is a
secondary screen with BH inside a named family. Holm at family-wise 0.05, m = 3. Recomputed after the
Borzoi fold-in enlarged endpoint 1 from 1,236 to 3,845 leads
(`c1-endpoint2-v2-20260915T084559Z/tables/c1e2_holm_family.tsv`; the E2 raw p is read from
`e2-transfer-20260914T203752Z/tables/defined_contrast_receipt.json`, key `primary_gain_bootstrap_p`).

| rank | member | raw p | Holm threshold | Holm-adjusted p | rejected |
|---|---|---|---|---|---|
| 1 | C1 endpoint 1, ChromBPNet vs Atlas ATAC, tier A4 (3,845 leads) | 1e-4 (floored at 1/10,000 draws) | 0.0167 | 3.0e-4 | yes |
| 2 | C1 endpoint 2, same models, allelic imbalance (495 sites) | 0.0036 weighted (0.0046 unweighted) | 0.0250 | 0.0072 | yes |
| 3 | E2 transfer contrast | 0.597 | 0.0500 | 0.597 | no |

**No verdict changed when endpoint 1 grew threefold.** Both C1 raw p values were already floored or small
enough that the enlargement moved the effect size (−0.0693 → −0.0910) without moving a rejection, and the
unweighted endpoint-2 p of 0.0046 gives an identical family (`c1-endpoint2-v2-.../holm_sensitivity.json`).

**What is rejected is "the two models are equal", in the direction opposite to the written prediction.** A
rejection here is not evidence for the predicted ordering. All three p values are bootstrap p values, not
permutation p values, and the 1e-4 floor means "no draw crossed zero", not "p = 0.0001".

---

## 8. Every prespecified prediction, met and not met

| id | prediction as written | observed | verdict | file |
|---|---|---|---|---|
| A1.1 | raw 2v2 allele-effect reliability in [0.70, 0.90] in all four contexts | HepG2_PAOA 0.6938 [0.6356, 0.7383]; other three inside | **not met** | a1-v2 `tables/reliability_by_context.tsv` |
| A1.2 | activity reliability exceeds allele-effect reliability in every context | holds in all four, and under every activity convention | met | a1-v2 `tables/activity_convention_sensitivity.tsv` |
| A1.3 | treatment-minus-control reliability below 0.40 in both cell lines | HepG2 0.0218; LX2 0.5931 | **not met** | a1-v2 `tables/reliability_by_context.tsv` |
| A1.4 | nine-combination convention differs from same-index by <0.05 | HepG2 0.0048; LX2 0.0021 | met | a1-v2 `tables/a1_4_convention_difference.tsv` |
| B3.1 | top-skill quintile more sequence-sensitive than bottom | CHIP_HISTONE +0.00166 [+0.00014, +0.00317], q 0.037; AVI_SCORE reversed | met, weakly and not across readouts | b3 `tables/b3_contrasts.tsv` |
| B3.2 | sensitivity rises with promoter status, falls with TSS distance | reversed in all four chromatin channels, all intervals exclude 0 | **not met** | b3 `tables/b3_contrasts.tsv` |
| C1.1 | ChromBPNet beats Atlas ATAC at Currin leads by 0.00–0.10, interval excludes 0 | −0.0910 [−0.1055, −0.0768] on the enlarged 3,845-lead tier A4; interval excludes 0, sign reversed, and further from the predicted band than the 1,236-lead −0.0693 | **not met** | c1e2-v2 `tables/c1e2_contrasts.tsv` |
| C1.2 | every model below its caQTL value at imbalance sites, and ≥1 interval includes 0 | first clause met for all **six** sequence columns (0.2963 vs 0.6450; 0.2424 vs 0.6504; 0.2874 vs 0.6791; 0.2411 vs 0.6079; 0.4439 vs 0.7361; 0.4121 vs 0.7120); second fails, lowest weighted bound 0.1314 (borzoi_dnase) | **not met** (first clause met) | c1e2-v2 `tables/c1e2_model_metrics.tsv` |
| C1.3 | allele-identity control at chance on both endpoints, abs < 0.05 | endpoint 1 tier A4 −0.0008; endpoint 2 0.1079 weighted, 0.1023 on all 500 with interval excluding 0 | **not met** at endpoint 2 | c1e2-v2 `tables/c1e2_model_metrics.tsv` |
| G1 | four-model endpoint-2 matched set reaches ≥450 of 500 | 495 sites in 54 blocks, **all four spec families present** | met (was not evaluable before the Borzoi fold-in) | c1e2-v2 `tables/c1e2_matched_counts.tsv` |
| G2 | allele-identity control keeps a positive weighted Spearman excluding 0, unweighted smaller | 0.1079 (0.0284–0.1798); unweighted 0.0554, includes 0 | met | c1e2-v2 `tables/c1e2_model_metrics.tsv` |
| G3 | endpoint-2 primary contrast keeps its negative sign | −0.1481 weighted, −0.1476 unweighted | met | c1e2-v2 `tables/c1e2_contrasts.tsv` |
| G4 | Borzoi ATAC endpoint-2 Spearman below its 0.6831 endpoint-1 value | 0.2874 at the imbalance sites vs 0.6920 on tier A (0.6831 on the wave-1 1,236) | met (was not evaluable) | c1e2-v2 `tables/c1e2_model_metrics.tsv` |
| G5 | tier A4 reaches ≥3,800 leads and the C1.1 sign survives | **3,845 leads in 1,743 blocks**, the archived AlphaGenome cap reached exactly; contrast −0.0910 (−0.1055, −0.0768) | met (was not met on size) | c1e2-v2 `tables/c1e2_matched_counts.tsv` |
| G6 | Borzoi ATAC endpoint-2 weighted interval excludes 0 | 0.2722 (0.1823–0.3809) | met (was not evaluable) | c1e2-v2 `tables/c1e2_model_metrics.tsv` |
| C2.1a | best DNA-LM delta head's gain over allele identity is **below 0.082** | HyenaDNA 0.1301 [0.1162, 0.1440]; Caduceus 0.0841 also above | **NOT MET** | c2 `tables/head_summary.tsv` |
| C2.1b | that gain's 95% interval **includes 0** | [0.1162, 0.1440]; 10,000/10,000 draws positive | **NOT MET** | c2 `tables/head_summary.tsv` |
| C2.2 | allele-identity control near chance on this label, abs < 0.05 | 0.0435 [0.0325, 0.0543] | met on the threshold; its interval excludes 0 | c2 `tables/head_summary.tsv` |
| C3.1 | curve flat at the top end, final increment below 0.01 | +0.0173, 5/5 seeds positive | **NOT MET** | c2 `tables/learning_curve.tsv` |
| D1.1 | hepatocyte ≥200 sites; stellate <50; macrophage <10; cholangiocyte <10 | 292 / 27 / 9 / 3 | met | d1 `tables/lineage_feasibility.tsv` |
| D1.2 | ALL_READS reproduces the deposit site count within 5 sites | 435 vs 435; 279 vs 279 | met (2 of 435 stored values differ; see defect D-08) | d1 `tables/lineage_feasibility.tsv` |
| D1.3 | reference bias negative in every lineage with ≥20 het sites | −0.025 to −0.063 across five lineages | met | d1 `tables/lineage_feasibility.tsv` |
| E1.1 | ≥1,000 MPRA variants carry a Currin nominal caQTL test | 3,894 | met, and weaker than it sounds (209 within 500 bp) | e1 `bridge_counts.tsv`, `e2_distance_strata.tsv` |
| E1.2 | ≤20 carry a measured allelic-imbalance value | 16 | met | e1 `bridge_counts.tsv` |
| E1.3 | ≥100 carry a Broadaway liver eQTL row | 202 | met | e1 `bridge_counts.tsv` |
| E2.1 | transfer gain below 0.02 with an interval covering 0 | −0.0015 [−0.0070, +0.0039], n 1,244 in 168 blocks | met | e2 `tables/defined_contrast_summary.tsv` |
| E2.2 | reporter head OOF macro on its own label ≥ 0.10 | 0.2178 over 4,358 elements | met | e2 `tables/reporter_head_self_evaluation.tsv` |
| E2.3 | reporter-OOF-alone abs Spearman < 0.05 on the caQTL label | +0.0105 [−0.0585, +0.0630] | met | e2 `tables/defined_contrast_summary.tsv` |
| E3.1 | native-flank correlation exceeds the N-padded one | random −0.009 [−0.070, +0.052]; top-effect +0.022 [−0.106, +0.160] | **not met** | e3 `tables/e3_correlations.tsv` |
| E3.2 | native-flank Spearman in [0.05, 0.25] | random +0.144 (met); top-effect +0.292 (above the band) | met / **not met** by stratum | e3 `tables/e3_correlations.tsv` |
| E3.3 | both arms below the 0.216 trained-head value | random +0.144 / +0.153 (met); top-effect +0.292 / +0.270 | met / **not met** by stratum | e3 `tables/e3_correlations.tsv` |
| E3.4 | the API accepts N-padded sequences | 296 of 296 returned finite values on both alleles | met | e3 `raw/api_scores.jsonl` |
| F1.1 | ≥40% of P5 pairs have D prime above 0.9 | 0.9502 (1,508/1,587) [0.9266, 0.9702] | met, by 55 points | f1 `tables/f1_summary.json` |
| F2.1 | the GNMT interaction residual falls inside the central 95% of its matched null | inside in all four channels, p 0.1188 / 0.2772 / 0.4059 / 0.5743 | met | f2-v2 `tables/f2v2_gnmt_verdict.json` |
| F2.2 | the GNMT arms are the same order as the family median, so F2.1 is uninformative | largest arms 0.0151 (rs2296805, H3K27ac) and 0.0321 (rs2296804, RNA) log2 against a family median abs single-variant effect of 0.0019–0.0060 | **NOT met** (the arms are large, so F2.1 is informative) | f2-v2 `tables/f2v2_summary.json`, `carried_forward_unchanged_from_v1/gnmt_primary_arms` and `.../per_channel/*/median_abs_single_variant_effect` |
| F3.1 | no pair survives BH q<0.10 against the corrected null | rows 3 / 6 / 18 / 34; distinct variant pairs 1 / 3 / 10 / 16; arm-magnitude matched 0 / 0 / 4 / 0 | **NOT met** | f2-v2 `tables/f2v2_summary.json` |
| F3.2 | the corrected null is not degenerate | observed/null median abs residual ratio 0.81 to 2.52 | met | f2-v2 `tables/f2v2_summary.json`, `carried_forward_unchanged_from_v1/family_level_contrast` |
| F3.3 | the corrected null's q95 differs from the deposit's by >10% in ≥1 channel | RNA +19.7% (0.007855 vs 0.006563), H3K27ac −13.3% (0.003133 vs 0.003612) | met | f2-v2 `tables/f2v2_summary.json`, `carried_forward_unchanged_from_v1/per_channel/*/null_q95_abs_residual`; the deposit-side comparator is the 300-draw P5 null recomputed in f2-v1 `tables/f3_summary_n2500.json`, `existing_p5_deposit/recomputed_null_here` |
| F2v2 P1 | the primary pair stays inside the central 95% of the repaired null in all four channels | inside in all four | met | f2-v2 `tables/f2v2_gnmt_verdict.json` |
| F2v2 P2 | the repaired RNA null's q95 is >10% smaller than the central-20-kb null's | −1.2% (0.008873 → 0.008767) | **NOT met** | f2-v2 `tables/f2v2_gnmt_verdict.json` |
| F2v2 P3 | the three chromatin q95 values agree between draw sets within 25% | ATAC +42.8%, DNase +64.3%, H3K27ac +55.7%; verdicts unchanged | **NOT met on the bound**, verdicts held | f2-v2 `tables/f2v2_gnmt_verdict.json` |
| F2v2 P4 | in-stratum recomputation reproduces 5 / 0 / 3 / 0 with RNA's five at the floor | counts reproduce; 3 of RNA's 5 at the floor, 2 at twice it | met on counts, **not on the floor** | f2-v2 `tables/f2v2_both_active_stratum.json` |
| F2v2 P5 | the measured both-active null rate exceeds 0.25% in every channel and the arms rank-correlate | 1.12 / 0.88 / 0.52 / 0.80%; Spearman +0.085 to +0.254 | met | f2-v2 `tables/f2v2_both_active_stratum.json` |

Section 4 of the specification sets no numbered predictions for B, so B1/B2 is scored against its nine
assigned corrections instead: eight were done in the layer itself, and correction 9's precondition — a
non-superseded dbSNP orientation deposit — was unmet when the layer was built (the third deposit had an
empty `tables/` at 20:47Z). **That precondition is now met.** `p0-dbsnp-orientation-20260914T210732Z`
carries a table and no supersession marker, and `b1-dbsnp-refresh-20260915T084700Z` recomputes the layer's
indel-orientation columns against it; the affected columns of the layer are named in its
`SUPERSEDED_COLUMNS_indel_orientation.md`, and the counts are in section 9 item 4.

---

## 9. Defects that change a number someone would read

The consolidated register with file, wrong value, correct value and the claim it touches is
`GWAS/finemapping/results/alphagenome_program/g-integration-20260915T072500Z/DEFECT_REGISTER.md`. The four
that most directly move a printed number:

1. **The p3a MPRA deposit swaps two contexts.** `p3a-benchmarks-20260909T190214Z/tables/mpra_dav_calls_S2_S5.tsv`
   labels TableS3 `HepG2_PAOA` (464 calls) and TableS4 `LX2_ctrl` (338); by content they are the other way
   round. The printed HepG2-PAOA stimulus-specific set of 230 should be **97** and the LX-2 TGFβ-specific
   set of 311 should be **240**. Every per-context MPRA number for those two contexts in that deposit is
   swapped, which exchanges a hepatocyte lipotoxic context for a stellate-cell baseline.
2. **The GNMT additivity prediction is recorded as holding for a pair that was never scored.**
   `docs/technical/alphagenome_atlas_execution_2026-09-09.md:618` says prediction H2 holds. The pair is
   absent from the 989-signal universe and from all 1,587 P5 pairs; four independent checks agree. It has
   now been scored, and it is additive on this scale, but the earlier record converted an untested question
   into a reported result.
3. **The allelic concordance printed in five places is superseded.** The same execution document prints
   0.655 with 33 blocks and p 5.5e-06 at lines 117, 191, 276, 437 and 556, including two headline summary
   tables, and records the correction to 0.6605 over 430 sites and 51 blocks, p 1.7e-06, only at line 580.
4. **Indel orientation in the P6f rescue and the P6d dossier the B1 layer reads came from a dbSNP deposit
   now marked SUPERSEDED.** A current, non-superseded orientation table now exists:
   `p0-dbsnp-orientation-20260914T210732Z`, which resolves **112,863 of 114,389** deferred indels (0.9867)
   against 76,780 (0.6712) for the deposit those two packages used. Measured against it by
   `b1-dbsnp-refresh-20260915T084700Z`: **38 of the 152 P6f targets** change orientation and **4 of the 25
   distinct dbSNP-oriented P6d uids** do, every one of them a same-site reciprocal, so those rows' signed
   readouts refer to the other allele. On the B1 layer, **145 of 548 indel rows** change `variant_class`
   and 157 change `orientation_source`; the 9,963 non-indel rows do not move. The newest rescue and
   dossier deposits (`p6f-indel-rescue-20260914T231617Z`, `p6d-dossier-20260915T000218Z`) are already on
   the current arbiter and show 0 changes. What did **not** change is the flip rate: the earlier deposit's
   38,438 flips against 38,342 agrees was read as an exact coin flip, and at the 1-Mb block the current
   deposit gives a block-mean flip rate of 0.5200 [0.4973, 0.5446] against 0.5004 [0.4891, 0.5122], with
   0.5 inside both. HSD17B13 rs72613567 is not affected: it resolves to `chr4:87310240:T:TA` in all four
   deposits, so the splice log2 −1.090 does not move.

Two authority conflicts are for the user to decide and were not edited:

- `docs/RESULTS.md:1825` states that weights trained on GSE267145 or GSE296875 do not leave the repository,
  which contradicts the adopted 2026-09-07 derivative-weight-release decision recorded in
  `docs/decisions/2026-09-07-derivative-weight-release.md` and in `config/resources.toml:55-56`.
- `config/datasets/gse267145_znf469_human_liver.toml:22` permits tier-2 base-resolution sequence use in
  `blockers[0]` and prohibits sequence extraction in `blockers[3]`, in the same array.

One independent-checker FAILED verdict is on the record and has been corrected: the wave-1 B1 check marked
`analysis_blocks_per_universe_NONEMPTY` FAILED because B_direct was printed as 29 named analysis blocks
when the correct count is **28**; the empty `analysis_block` string of 10 B_direct signals had been counted
as a block. The v2 deposit prints 28 and B_direct-only 18 rather than 19.

---

## 10. What is NOT established

- **That AlphaGenome is a better liver-accessibility sequence model than the project ChromBPNet.** The two
  are not held out on the same terms and their columns are not averaged over the same number of tracks.
  The measured gap is an upper bound on any genuine advantage.
- **That reporter supervision does not help endogenous prediction in general.** One route is bounded:
  one extra out-of-fold reporter column on a ridge over the same frozen 4-kb HyenaDNA delta embedding, on
  1,244 variants within 5 kb of a peak centre. Joint training, fine-tuning, a multi-task head, another
  model and another endogenous label are all untested, and this design's own ceiling, maximised over the
  three nested positive-control arms, is +0.0142 [+0.0067, +0.0246].
- **That the C2 gain exceeds the reporter task's gain.** The reference interval [0.0274, 0.1389] overlaps
  this one [0.1162, 0.1440]. What is established is that the gain here is large and excludes zero.
- **Any causal claim.** Currin peak leads are LD-selected marginal associations and a lead is not
  necessarily the causal base; the imbalance label is read-level. A macro Spearman of 0.174 does not mean
  the model identifies the causal base.
- **That ChromBPNet and Borzoi differ at the allelic-imbalance sites.** +0.0304 weighted [−0.0864, +0.1480],
  p 0.58 on 54 blocks is a non-separation, not an equivalence: the interval is wide enough to contain a
  difference as large as the ChromBPNet-to-AlphaGenome gap. (The earlier statement here, that nothing was
  known about Borzoi at these sites or about a four-model matched set, is closed: both GPU jobs completed
  on 2026-09-14 and `c1-endpoint2-v2-20260915T084559Z` delivers the four-family comparison and the
  exact 200-lead cross-run reproduction.)
- **Anything about the 28,458 Currin leads no GPU panel covers.** The Borzoi run stopped at its priority-1
  tier by design, because the archived AlphaGenome panel caps the four-family matched set at 3,845 leads
  and the remaining leads cannot enlarge it. Enlarging endpoint 2 past 495 sites would need an AlphaGenome
  call at the 5 unserved sites, a hosted-API action this program did not take.
- **That the official Calico Borzoi weights would give these numbers.** What was scored is the Hugging Face
  PyTorch port `johahi/borzoi-replicate-0`, whose weight terms are UNDECLARED, so it may be evaluated
  internally and not redistributed; the audit calls the ports independently versioned derivatives whose
  hashes and declared licences are not evidence about the official TensorFlow objects
  (`config/artifacts/models/borzoi/{exposure_audit,checkpoints}.json`). Borzoi's training corpus also
  includes GTEx v8, so its liver column is not an unseen-tissue comparison.
- **Anything about a lineage other than hepatocyte, and marginally endothelium.** The negative result is
  scoped to 1,119 Atlas-scored SNVs inside consensus peaks, this read depth, this filter and these 30
  donors. It does not say those lineages have no allelic imbalance.
- **A donor-axis claim anywhere.** The C2/C3 learning curve varies training loci only. D1's donor counts
  size a label; they do not estimate an effect.
- **Non-additivity at any local haplotype.** The F3 rejections are bounded by readout-window overlap below
  2 kb, by arm magnitude, and by the family counting the same variant pair up to six times. The GNMT
  verdict is additivity on a model's prediction scale for two arms that are not population haplotypes.
- **Any statement about phase in a MASLD donor.** F1's haplotype frequencies are 1000 Genomes Phase 1 EUR,
  statistically phased and imputed, and both LD sources are EUR only.
- **That a hosted-API log2 track ratio is comparable to an Atlas quantile.** They are separate scales and
  are never pooled in any table here.
- **Peak-level or genetics-overlapping saturation.** The archive holds one universe, U1 H3K27ac, on
  chr1–chr12 complete plus chr13 at 63% (1,558 of 2,474 regions), chr14 onward at zero; U2's 7,000 snATAC
  peak regions are 0 on disk.
- **That the layer's 0.364 readable share of (signal, variant) rows — 3,822 of 10,511, not to be read as a
  share of signals, which is 623 of 989 — is biological.** It is a property of the Atlas point-query API
  and of the posterior.
- **Anything a licence text says.** Every rights verdict in A2 is a faithful reading of what this
  repository records, cited to file and line; four sources in active use (Currin caQTL, Broadaway liver
  eQTL, GTEx v8 liver sQTL, `leafcutter_joint_v2_bulk_cohorts`) have **no registry entry at all** and have
  had no licence check of any kind.

---

## 11. Attachments

| File | What it is |
|---|---|
| `GWAS/finemapping/results/alphagenome_program/g-integration-20260915T072500Z/research_layer_sidecar.tsv` | 989-row Gene Catalog attachment, keyed (ensembl_id, credible_set_id) for the 816 signals with a gene and on signal_uid for the 173 without; every row `candidate_not_adopted`; measured and predicted columns kept apart, each with its sign convention named |
| `.../research_layer_sidecar_dictionary.tsv` | per column: measurement or model prediction, unit, sign-convention column, source file |
| `.../DATA_DEVELOPMENT_PROPOSAL.md` | what to acquire next, sized from the measured variance and block counts here |
| `.../DEFECT_REGISTER.md` | one consolidated defect list, including the two authority conflicts |
| `.../RESULTS.md` | what the integration package itself ran |
| `.../MANIFEST.tsv` | every input path with sha256 |
