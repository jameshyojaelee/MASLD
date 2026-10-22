# Current results and number ledger

**Corrected Figure 2 results adopted by the user, 2026-10-13.** Member-supported restoration of
Figure 2 positions gives a descriptive protein-altering plus canonical-splice
fraction of 23.34% for direct traits and 3.87% for enzyme traits, ratio 6.03,
pooled across ancestries. H58 (2026-10-15) shows this pooled value mixes trait
class with ancestry and fine-mapping coverage: EUR 29.88% vs 6.92% (difference
0.2296 [0.0550, 0.4236]), AMR 6.73% vs 4.54% (0.0219 [-0.0489, 0.1065]); no
stratum receives a trait-class call.
The older coordinate-dependent intervals, context comparisons and 62 robustness
refits are withdrawn. The corrected-coordinate reruns completed on
October 13: coding/splice ratio 6.026671, physical-cluster bootstrap interval
2.57459–13.51016, and all 57 comparisons in the corrected robustness menu.
These Figure 2 estimates are adopted; the corrected Figure 2D and unassembled
noncoding-context companion are installed with their corrected source tables.
Historical named releases are preserved. This October 13 decision concerned
Figure 2 coordinates; COLOC r2 is selected separately on October 14.
Figure 5B/5C and the class comparison now have completed 774-unit counterparts;
the September 844-record estimates are historical. Full evidence is in the
[October 13 adjudication](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#review-of-this-review-2026-10-13).

**Historical Figure 3F matrix, 2026-10-05, superseded October 14:** 61 of the completed September
436 primary SuSiE genes have BH FDR <0.05 in at least one existing adjacent
fibrosis contrast. Four complete 23,702-gene families yield 244 displayed
effects. Each genetic lane reproduces the strongest tier-1/2 driver's eligible
exact-pair PP.H4. No PIP estimates are displayed. Contrast denominators are
209/195/165/100 library records and 206/192/161/97 provisional genotype
components; one component enters both F0 and F1. These are historical fits,
not new donor-qualified estimates. The 61 genes, 244 effects, 23,702-gene
families and 209/195/165/100 denominators are on the 844-record basis and are
historical.

**Figure 3F on the 774 units, rendered 2026-10-12 and installed 2026-10-14:**
`figures/candidates/fig3f-no-pip-774-20261012T011128Z/`
(`source_tables/fig3f_stage_genetics_matrix.tsv`,
`fig3f_contrast_denominators.tsv`, `checks.txt`). Same selection rule: SuSiE
PP.H4 > 0.5 and BH FDR < 0.05 in at least one adjacent fibrosis contrast, each
contrast its own 24,115-gene BH family (no joint FDR over the four). At BH
q < 0.0125 per contrast, which bounds the FDR over all discoveries in the four
families at 5% but not specifically within the SuSiE subset, 15 of the 54 genes
and 16 marks remain: ACVR2B, BRI3, CDH6, CORO1C, CREB3L4, CTSH, CYP2A6, EPSTI1,
F2RL1, FAM177A1, GCAT, MTARC2, OMG, SLC17A2, TNFSF10 (review R16-02). Of the
436 primary SuSiE genes of COLOC release `susie_coloc_r2_20260924T053244Z`
(selected 2026-10-14), 426 are in the 774-unit stage table and 54 are selected, giving 216
displayed effects (61 and 244 on the 844 records). 59 effects have BH FDR < 0.05
(65 before): 30, 13, 11 and 5 for F0→F1, F1→F2, F2→F3 and F3→F4 (32, 21, 7 and 5
before). Twenty-nine of the 54 genes are selected only through an F0→F1 mark. Of the 30 F0→F1 marks, 18 keep BH q < 0.05 (same sign) with disease-only F0 and 0 without GSE130970 (smallest q 0.103). Marks kept (q < 0.05, same sign) without GSE130970 / GSE135251 / GSE162694, for F0→F1, F1→F2, F2→F3 and F3→F4: 0/6/6/0, 15/3/1/2 and 12/0/3/4; 1 of the 59 (OMG, F2→F3) is kept under all three omissions (review repair job 7550; `figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/repair/fig3f_marks_cohort_omission.tsv`). Reselecting with disease-only F0 gives 45 genes and 50 marks (21 / 13 / 11 / 5): 42 of the 54 stay, 12 leave (ACVR2B, HSPA4, MGAT2, MTARC2, MTOR, NUDT4, PABPC4, PFKFB2, PMF1, RC3H2, RPL8, WIPI1) and FCHO2, HPS6 and SEC14L2 enter. Each contrast has a complete 24,115-gene BH family and uses 205, 190,
157 and 94 units, with no identity cluster on both sides of a contrast. Ten
genes leave the panel (ADGRL2, ARMC9, ARNT, COLEC11, COPG2, ENTPD5, FLRT3, PKN3,
PPP1R3B-DT, SLC2A2) and three enter (OMG, EPHA2, GBF1); 51 are in both.

**Bulk refit adopted by the user, 2026-10-11.** The five-cohort
disease-versus-control analysis of the paper is the 774-unit refit
(`RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/`):
138 control and 636 disease units, 24,115 genes tested, 1,854 canonical calls
(1,333 up, 521 down) and 2,131 interval-null calls. Units are genetic-identity clusters,
one library each, not verified donors. The 855 held annotations are counted and
not tested. Statements below that give 844-record numbers are historical or are
labeled as not rebuilt. The values were re-derived from the count matrices in
`Analysis/Hypotheses/refit-adoption-review-20261011T151011Z/`.

**Historical Figures 3–4 identity qualification, 2026-09-29:** the genotype ledger identifies
repeated bulk libraries and conflicting source labels; GSE193066 identity calls
remain unresolved. Historical bulk roster counts and newly reconstructed fits
are conditional, not verified independent-donor inference. The correction
candidate includes an explicit roster join and a separate provisional identity
sensitivity. That restriction applied to the September candidates; the user
adopted the corrected-count bulk refit on genetic-identity clusters on October
11. Neither that decision nor the historical network analysis verifies
independent donors or qualifies GSE193066; see
[the disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).

**Correction implementation, 2026-10-05 (historical; superseded by the
2026-10-11 adoption above).** At that date no new donor-qualified bulk estimate
was available. The 1,725-DEG corrected-controls fit and 1,833-DEG name-fold
sensitivity both use historical library rosters; the latter does not establish
gene equivalence from shared names. Full annotation correspondence and supported
analysis units are required before a primary refit. The genetics-only Figure 1C
candidate partitions the completed September 436-gene set as 414 enzyme-only,
12 direct-disease/liver-fat-only and 10 both, using any eligible trait support
rather than the strongest trait driver. No synchronized result adoption follows.
See the [change and check log](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#resource-corrections-implementation-2026-10-05).

The recovered primary identity candidate contains 774 source-qualified genetic
clusters for bulk/continuum inputs, 392 endpoint records and 418 stage records.
The strict-depth sensitivity contains 738, 385 and 407 respectively. These are
candidate roster sizes, not verified independent-donor counts. GSE193066 fails
the source-person crosswalk and remains excluded from
primary replication, repeat-biopsy ICC and derived ceiling claims. The
annotation candidate retained 855 unresolved correspondences after official
reference exclusions and exact transcript-coordinate checks. On 2026-10-07 the
user decided that these 855 are held and not tested; the corrected bulk total is
1,854 canonical calls in 774 units. Annotation r8 relabels those 855 as
`held_untestable` and leaves the other 85,514 annotations and counting targets
unchanged. The 776-library recount in
`figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/`
and the October 8 corrected-count refit are complete; the user adopted the
774-unit bulk refit on October 11. This does not adopt its remaining dependent
figures or COLOC r2. Both Figure 7 review sets have six panels; the
framework exposure reconstruction retains 54 clean, 22 confounded and eight
unresolved encoder–head–study cells. These development comparisons remain
conditional on their saved models and study-specific identities.

The user-authorized annotation policy excludes seven contamination-scaffold
and 56 complete false-duplication annotations from counting, with records
retained as non-applicable. TRPM3 remains included. The revised provisional
Catalog retains all 28,017 records and marks the 14 represented reference
artifacts non-applicable for bulk RNA; source evidence classes remain historical.
The revised NMF protocol keeps fixed weights on 4,517 unchanged counting objects
(4,438 after the existing continuum signature exclusion). Historical-count
normalization diagnostics meet the fixed reconstruction criteria for eight
of ten stage axes and ten of ten continuum axes; stage k4-P2 and k6-P6 do not
meet the global Spearman ≥0.99 criterion and remain unavailable in that protocol.
These diagnostics do not establish corrected-count transport or replace original
scores. Since then the primary bulk fit and fixed-NMF scoring have run on the
corrected 774-unit counts (2026-10-08 and 2026-10-09).

Both identity-only disease-versus-control fits are complete on unchanged
historical unique-only counts, adjusting for dataset and inferred sex. These
three rows are historical diagnostics; the last two rows are the fits on the
corrected counts:

| Fit | Libraries / supported analysis units | Tested annotation rows | Canonical calls (up / down) | interval-null calls |
|---|---:|---:|---:|---:|
| Historical corrected controls (identity-only diagnostic) | 844 / unknown | 23,702 | 1,725 (1,261 / 464) | 2,036 |
| Primary identity roster (identity-only diagnostic, historical counts) | 774 / 774 | 23,535 | 1,756 (1,243 / 513) | 2,014 |
| Strict-depth identity roster (identity-only diagnostic, historical counts) | 738 / 738 | 23,940 | 1,807 (1,254 / 553) | 2,009 |
| **Adopted 774-unit refit, corrected counts** | 774 / 774 | 24,115 | 1,854 (1,333 / 521) | 2,131 |
| Strict-depth sensitivity, corrected counts (four-cohort disease contrast: GSE126848 keeps 5 controls and no disease unit) | 738 / 738 | 24,509 | 1,906 (1,344 / 562) | 2,135 |
| Library-complexity sensitivity, release sample weights fixed, + within-cohort z(log library size; total counts on the tested genes) (review R03-09) | 774 / 774 | 24,115 | 1,823 (1,351 / 472); 1,760 of 1,854 kept | not run |
| Same, + z(log normalized library size; library size x RLE factor) | 774 / 774 | 24,115 | 1,715 (1,289 / 426); 1,648 kept | not run |
| Same, + z(top-100-gene share) | 774 / 774 | 24,115 | 1,815 (1,280 / 535); 1,738 kept | not run |
| Same, + z(expression entropy) | 774 / 774 | 24,115 | 1,814 (1,287 / 527); 1,759 kept | not run |
| Same, + z(log library size) + z(top-100 share) + z(detected genes) | 774 / 774 | 24,115 | 1,605 (1,085 / 520); 1,563 kept (84.3%); all-gene Spearman 0.984 | not run |
| Same, + z(log normalized library size) + z(top-100 share) + z(detected genes) | 774 / 774 | 24,115 | 1,590 (1,100 / 490); 1,569 kept (84.6%); all-gene Spearman 0.988 | not run |

Supported units are source-linked participants or conservative genetic identity
clusters, not verified independent people. Canonical calls require eBayes
BH FDR <0.05 and absolute log2 fold change >0.5; interval-null tests the interval null
of absolute log2 fold change ≤0.25. Each BH correction covers all tested rows
separately within its roster and test family. The primary roster retains
1,572 historical canonical calls, loses 153 and gains 184; the strict roster
retains 1,539, loses 186 and gains 268. Selection also changes cohort/case
composition, filtering, normalization and quality weights, so these differences
do not isolate the effect of removing repeats. Duplicate counting targets and
unresolved identities remain among called rows; these are not corrected
physical-gene DEG totals. Saved results are in
`figures/candidates/resource-corrections-20261005-v1/identity_only_comparison_r2/`.

## Figures 3–4 correction candidate, 2026-09-29

These values are candidate analyses, not a replacement for adopted memberships.
The [technical disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md)
links complete source tables, definitions, denominators and limitations.

- GSE213621 source categories: Control 68, F0–F1 95, F2 107, F3–F4 91.
  Disease-only continuum n=402 (109 GSE162694 + 293 GSE213621); original
  categorical relabeling leaves coefficients unchanged.
- Disease-only fixed-projection coefficients, stage/sex then additionally
  composition-adjusted: Stromal ECM 0.726 [0.626, 0.825] then 0.488
  [0.344, 0.632]; Ductular injury 0.677 [0.599, 0.755] then 0.639
  [0.549, 0.729]. HC3-based fixed-effect 95% intervals; full 117-program family.
  These four coefficients are on the 844-record basis and are historical.
- The same coefficients on the 774 units, 2026-10-12 (candidate, not adopted):
  `figures/candidates/continuum-composition-774-20261012T012137Z/`
  (`comparison/meta_old844_vs_new774.tsv`, `cohort_old844_vs_new774.tsv`,
  `call_status_summary.tsv`, `cohort_n_by_arm.tsv`; script
  `scripts/analysis/fig34_corrections/01_continuum.py`). Model: GSE162694 and
  GSE213621 fitted separately; OLS of the fixed program or system score on the
  continuum axis (fixed projection, disease-only arm) plus stage and sex; the
  composition adjustment adds 15 ILR coordinates of 16 cell-type fractions; HC3
  errors; fixed-effect inverse-variance meta-analysis of the two cohorts; BH
  within family (117 programs, 116 estimable; 43 systems). Unit: one
  genetic-identity cluster per person, n=367 (109 + 258; 402 = 109 + 293
  before). The 844-record values reproduce exactly from the default inputs
  (`reproduction_844_stored_vs_rerun.tsv`, largest difference 0).
  Stromal ECM: 0.760 [0.654, 0.866] then, composition-adjusted, 0.294 [0.154,
  0.434], BH q=1.4e-4 (844 records: 0.726 then 0.488). Ductular injury: 0.707
  [0.626, 0.788] then 0.601 [0.500, 0.701] (844 records: 0.677 then 0.639). The
  direction and the BH call hold for both; the adjusted Stromal ECM effect is
  smaller. Per cohort, composition-adjusted Stromal ECM is 0.172 [-0.032, 0.377],
  BH q=0.233 in GSE162694 (was 0.432, q=6.1e-4), so it is no longer supported in
  that cohort alone, and 0.408 [0.212, 0.604], q=2.3e-4 in GSE213621 (was 0.547).
  The cause is the rebuilt cell-type fractions, not the roster: with the 774-unit
  scores and the earlier fractions the adjusted coefficients are 0.479 [0.331,
  0.627] for Stromal ECM and 0.657 [0.564, 0.749] for Ductular injury
  (`sensitivity_774_scores_earlier_fractions/continuum/continuum_meta_effects.tsv`).
  The composition-adjusted Stromal ECM estimate therefore depends on which
  deconvolution input is used.
  Calls at BH q<0.05 (disease-only, fixed projection, HC3), 844 records to 774
  units: programs unadjusted 91 to 94 (3 gained, none lost), composition-adjusted
  67 to 67 (8 gained, 8 lost); systems unadjusted 38 to 37, composition-adjusted
  11 to 24 (14 gained, 1 lost). Programs lost after adjustment: cholangiocyte
  Hypoxia (NDRG1); fibroblast ECM adhesion (TNC), Myofibroblast (EBF1/SLIT3) and
  Activated fibroblast (MMP19); hepatocyte Cell polarity (PARD3) and Leukocyte /
  immune (DOCK2); macrophage Glycolysis / hypoxia (HK2); T-cell NK / cytotoxic
  (KLRD1). Gained: cholangiocyte Complement / plasma (CFH); fibroblast Portal
  fibroblast (NR1H4) and Actomyosin / MHC-I; hepatocyte UPR / ERN1 stress and
  Sinusoidal endothelial (STAB2); macrophage Hep-like metabolic (CYP3A5); T-cell
  Translation / MYC (EEF1) and Translation / ribosome (EEF1A1). Systems gained:
  05, 08, 11, 12, 13, 16, 18, 19, 24, 27, 32, 33, 36 and 38; lost: 42.
- Why the composition-adjusted calls changed, examined 2026-10-12
  (`figures/candidates/continuum-composition-774-diagnosis-20261012T013141Z/fraction_diagnosis/`:
  `three_arm_meta_disease_only_fixed_projection_hc3.tsv`,
  `call_transitions_by_step.tsv`, `changed_systems_three_arms.tsv`,
  `axis_r2_on_ilr.tsv`, `fractions_old_vs_rebuilt_by_cell_type.tsv`,
  `integrity_fraction_files.tsv`, `integrity_units_rank.tsv`,
  `rederivation_vs_stored_cohort_effects.tsv`; and
  `fraction_diagnosis_supplement/adjusted_axis_earlier_vs_rebuilt.tsv`). Three
  arms were compared: 844 records with the earlier fractions, 774 units with the
  earlier fractions, and 774 units with the rebuilt fractions. Composition-adjusted
  calls are 11, 10 and 24 of 43 systems and 67, 65 and 67 of 117 programs;
  adjusted Stromal ECM is 0.488, 0.479 and 0.294. The change is at the fractions
  step (systems 14 gained, none lost; programs 8 gained, 6 lost), not at the
  roster step (systems 1 gained, 2 lost; programs none gained, 2 lost). The 14
  systems gained at the fractions step on the 774 units are 08, 11, 12, 13, 16,
  18, 19, 23, 24, 27, 32, 33, 36 and 38; system_05 was already called with the
  earlier fractions on the 774 units, and system_23 is not in the 844-to-774
  list of the previous bullet because it was lost at the roster step. In the
  14 systems gained overall the adjusted coefficient is 1.2 to 3.3 times larger
  with the rebuilt fractions, while the median standard-error ratio over the 43
  systems is 1.03, so the calls come from larger estimates, not smaller errors.
  Stage, sex and the 15 ILR coordinates explain the same share of the axis with
  either fraction set (GSE162694 0.8836 and 0.8845; GSE213621 0.8565 and
  0.8696); the 12 to 14% of axis variance left over correlates 0.82 and 0.87
  between the two sets, so the two adjustments leave different residual axes.
  Systems 13 and 36 have opposite signs in the two cohorts under the rebuilt
  fractions. No defect was found in the row matching, row order, zero floor,
  cell-type columns or design rank, and all 86 system-by-cohort coefficients
  re-derive to within 3.1e-13 in each arm. In the modelled disease-only units
  the mean fibroblast fraction falls (GSE162694 0.055 to 0.037; GSE213621 0.028
  to 0.019) and the macrophage fraction rises (0.029 to 0.045; 0.043 to 0.052).
  Marker check, 2026-10-12
  (`figures/candidates/continuum-fraction-markers-substitution-20261012T014841Z/marker_agreement/marker_agreement_by_cell_type.tsv`,
  `marker_agreement_tally.tsv`): Spearman rho across persons between each
  fraction and a marker-gene score from the 774-unit voom expression, in the
  modelled disease-only persons (GSE162694 n=109, GSE213621 n=258), earlier
  then rebuilt fractions. Fibroblasts 0.604 to 0.830 and 0.379 to 0.671;
  cholangiocytes 0.462 to 0.602 and 0.413 to 0.494; hepatocytes 0.287 to 0.187
  and 0.104 to -0.204. Macrophages, endothelial cells, the T/NK sum and the
  B/plasma sum differ by 0.04 or less in both cohorts. Paired bootstrap over
  persons, 2,000 draws, seed 20260929: of 24 cohort-by-fraction comparisons the
  rebuilt set is higher in 10, the earlier set in 4, and the interval of the
  difference includes zero in 10. This is a consistency check, not a
  validation, because fractions and marker scores come from the same expression
  data. Neither fraction set is better overall: the rebuilt set agrees better
  for fibroblasts and cholangiocytes and worse for hepatocytes.
  Attribution to cell types (`substitution/single_cell_type_target_features.tsv`,
  `summary_by_cell_type_14_gained_systems.tsv`,
  `cumulative_shared_order_call_counts.tsv` in the same directory): replacing
  only the fibroblast fraction with the rebuilt one moves adjusted Stromal ECM
  by -0.202 of the -0.185 total change (GSE162694 -0.253, GSE213621 -0.049).
  No single cell type accounts for the 14 gained system calls: replacing one
  type gives 10 to 19 system calls (B cells 19, cDC2s 18), not 24. B cells have
  the largest single effect on the 14 systems (median coefficient change 0.030)
  and are below the 1e-6 floor in 159 (earlier) and 120 (rebuilt) of 258
  GSE213621 persons (`figures/candidates/continuum-rare-type-arms-20261012T020132Z/arms/below_floor_16_types.tsv`).
  Replacing one type and re-closing the fractions is approximate: the 16
  single-type changes sum to 1.23 times the total for Stromal ECM.
  Alternative composition codings, prespecified before any outcome was read
  (`figures/candidates/continuum-rare-type-arms-20261012T020132Z/`: `PRESPEC.md`,
  `DEVIATIONS.md`, `arms/summary_arm_by_fraction_set.tsv`,
  `arms/agreement_between_fraction_sets_by_arm.tsv`). Same 367 persons, scores,
  axis and model; only the coding of composition changes. Calls are BH q<0.05
  within 43 systems and within 117 programs (116 estimable); coefficients are
  fixed-effect meta estimates with 95% intervals.

  | Arm | Fraction set | ILR coordinates | System calls of 43 | Program calls of 117 | Stromal ECM | Ductular injury |
  |---|---|---:|---:|---:|---|---|
  | A, 16 types, floor 1e-6 | earlier | 15 | 10 | 65 | 0.479 [0.331, 0.627] | 0.657 [0.564, 0.749] |
  | A, 16 types, floor 1e-6 | rebuilt | 15 | 24 | 67 | 0.294 [0.154, 0.434] | 0.601 [0.500, 0.701] |
  | B, rare types pooled, 13 parts | earlier | 12 | 19 | 70 | 0.418 [0.275, 0.562] | 0.730 [0.641, 0.819] |
  | B, rare types pooled, 13 parts | rebuilt | 12 | 23 | 59 | 0.203 [0.060, 0.347] | 0.643 [0.548, 0.738] |
  | C, lineage level, 7 parts | earlier | 6 | 35 | 76 | 0.464 [0.357, 0.571] | 0.691 [0.617, 0.766] |
  | C, lineage level, 7 parts | rebuilt | 6 | 34 | 62 | 0.328 [0.174, 0.482] | 0.664 [0.567, 0.761] |
  | D1, 16 types, floor 1e-4 | earlier | 15 | 25 | 60 | 0.449 [0.309, 0.589] | 0.678 [0.596, 0.760] |
  | D1, 16 types, floor 1e-4 | rebuilt | 15 | 20 | 53 | 0.233 [0.089, 0.376] | 0.654 [0.557, 0.752] |
  | No composition terms | none | 0 | 37 | 94 | 0.760 [0.654, 0.866] | 0.707 [0.626, 0.788] |

  Arm B pools plasma cells, cholangiocytes, B cells and neutrophils; arm C
  joins cholangiocytes to hepatocytes under the prespecified merge rule. Arm D2
  (floor 1e-3) is not estimable with the rebuilt fractions in GSE213621 (with
  the earlier fractions: 31 systems, 72 programs, Stromal ECM 0.479, Ductular
  injury 0.644). Arm D3 (multiplicative replacement at 1e-4) gives the same
  call counts as D1 and focal coefficients equal to three decimals. No arm makes the two
  fraction sets agree by the prespecified criteria. The lineage arm is closest:
  21 of 160 features are called in one set only (28 in arms A and D1), with
  Pearson r between sets 0.960 for systems and 0.988 for programs. It also
  explains less of the axis (R2 0.75 to 0.80 against 0.86 to 0.88 in arm A),
  so part of its agreement comes from weaker adjustment. Pooling rare types
  does not help (35 called in one set only). Ductular injury is stable in every
  estimable arm (0.60 to 0.73; unadjusted 0.707). Stromal ECM is called at BH
  q<0.05 in every estimable arm, but its adjusted size runs from 0.203 to 0.479
  by coding and fraction set (unadjusted 0.760), and the rebuilt set is below
  the earlier set by 0.14 to 0.22 in every arm.
  Reporting rule, agreed by the user on 2026-10-12: composition-adjusted
  continuum results are reported as dependent on the deconvolution input, and
  no individual adjusted system or program is named as a finding (the lists in
  the previous bullet are a record of what changed, not findings). Ductular
  injury is reported as holding after composition adjustment (0.60 to 0.73).
  Stromal ECM is reported as reduced by composition adjustment by an amount
  that depends on the deconvolution input, with the range given (0.203 to 0.479
  against 0.760 unadjusted). No adjusted system is named. The unadjusted
  results, 37 of 43 systems and 94 of 117 programs, do not use the fractions.
- Targeted corrected content-lineage MASH-versus-steatosis follow-up:
  fibroblast Stromal ECM n=26, beta=0.329 [0.098, 0.560], BH q=0.0106;
  cholangiocyte Ductular injury n=26, beta=0.0848 [−0.135, 0.304], q=0.449.
  These are projected scores on common corrected/raw cells, not native Hotspot
  units; partial correction coverage and broad-lineage substate mixture remain limits.
- Three protein correlations reproduce in 46 common participants. Corrected
  covariate-df three-pair BH q=0.00243, 0.00772, 0.00243. Paired NAS-adjustment
  attenuation intervals exclude zero, without establishing mediation.
- Fixed C3 cross-content covariance after histology and one residual component:
  discovery Fisher-averaged r=0.288 [0.236, 0.346], n=469; retrospective
  evaluation r=0.189 [0.117, 0.289], n=105. Participant-bootstrap intervals are
  conditional on fixed gene calibration; sparse categorical strata limit usable draws.
- Figure 3G matched interval-null sensitivity, superseded value from this 2026-09-29 candidate (406 pairs): OR=1.32, exact 95% CI 0.761–2.315,
  nominal McNemar P=0.358. Current 774-unit value (426 pairs): paired OR=1.321, nominal exact McNemar P=0.321
  (`figures/candidates/fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_matched_treat_summary.tsv`). Canonical overlap and interval-null remain different definitions.
- Figure 4H's historical correction candidate preserved 18 strict maxT
  conjunction outlines. The current main panel uses the 367 disease units and plots the donor-first system estimates (rho 0.670, 34/43 same direction; feature-first 0.631 and 32/43) and has no outline; the stage unit is score SD per cohort stage-code SD.
  Separate competitive tests support 0/43 connected-null and 1/43 uniform-null systems.


**Verified:** 2026-08-17 EDT  
**Continuum candidate verified:** 2026-08-21 EDT  
**Rule:** every number below is labeled by release state. Do not combine release
states silently.

## Release labels

| Label | Meaning |
|---|---|
| **Current official** | Read-only July release; usable when that release is named |
| **Validated update** | Recomputed and checked; intended for the Resource but not yet adopted with every dependent output file |
| **Validated candidate** | Complete, checked analysis retained for candidate figures or narrative integration; not an official release until synchronized adoption |
| **Preflight only** | Source/design audit completed; inferential model not run |
| **Pending rebuild** | The corrected COLOC compute is complete and adopted; these derived values stay unauthorized until they are regenerated against the adopted release |
| **Supplementary** | Valid boundary or null analysis; not a main-paper discovery |

## Continuum translation extension: validated candidate, 2026-09-06

| Quantity | Candidate result | Unit and interpretation |
|---|---:|---|
| Complete Gene Catalog sidecar | 28,234 | Unique gene identifiers; includes all 23,370 bulk-tested genes plus genetics-only genes (23,370-gene, 844-record basis; not rebuilt on the 774 units) |
| Retained replay evidence | 5,362 signal pairs; 22,356,697 variant records | All exported pairs from 816 main-trait gene–study records; lower-support pairs remain visible |
| Supported main-trait signal pairs | 925 | SuSiE PP.H4 >0.5; may include several signals per gene–study record |
| Direction-resolved supported signal pairs | 83 | Single exported signal on each side and ≥95% original shared posterior mass in one resolved marginal direction; no conditional or drug effect |
| Stage/sex-supported programs | 93/117 | Both co-primary meta-analysis BH q<0.05 and all four cohort/axis directions |
| Composition-adjusted supported programs | 68/117 | Same rule and participants, adding 15 lineage log-ratio coordinates |
| Composition-adjusted separate-cohort support | 26/117 | BH q<0.05 separately in both evaluation cohorts and both axes, with consistent direction |
| Evaluation participants | 109 + 361 = 470 | GSE162694 and GSE213621; complete denominators per testable program |
| Matched liver/plasma participants | 41 | 11 controls and 30 with MASLD; existing public cohort, unpaired with bulk RNA |
| Matched shared-protein family | 345 | Six untestable members retained in the BH denominator |
| Adjusted matched protein associations | 49/345 | BH q<0.05 after fibrosis, age, BMI, sex and acquisition-batch adjustment; not liver origin or treatment response |

Direction categories overlap at gene level. The liver-fat results include
opposing directions across two PDFF studies for the same gene; individual
signals are retained without a common directional recommendation. No
replicated direct-MASLD opposing-direction hero is asserted. Exact tables,
figure and limitations are linked in the
[technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md).
These are candidate outputs, not changes to the adopted genetic memberships,
evidence classes, or main figure selection. They were computed on the 2026-08-10
bulk fit (23,370 genes) and the 2026-08-17 COLOC release (816 gene–study
records; 767 after the corrected release) and have not been rerun
(844-record basis; not rebuilt on the 774 units). The
corrected Gene Catalog release (user-approved 2026-09-25) is described under
"MASLD Gene Catalog" below.

## Human bulk RNA-seq

| Quantity | Current official release | Validated Resource update (2026-08-10) | Corrected controls, 844 records (2026-09-23 refit; historical) | Adopted 774-unit refit (2026-10-08; adopted 2026-10-11) | Interpretation |
|---|---:|---:|---:|---:|---|
| Count-matched human samples | 1,281 | 1,281 | 1,281 | not re-ledgered | Nine retained cohorts in the corrected source registry |
| QC-passing human samples | 1,260 | 1,257 | 1,257 | not re-ledgered | Complete source substrate; three current metadata rows lack a count match |
| Cohorts in pooled disease-vs-control analysis | 5 | 5 | 5 | 5 | Same intended five-cohort comparison. In the 774-unit fit GSE126848 keeps 23 control and 7 disease units, all male |
| Pooled samples | 846 | 844 | 844 | 774 units | Updated analysis uses fragment-counted libraries and removes two failed samples. The 774 are genetic-identity clusters, one library each, not verified donors |
| Controls / disease | not re-ledgered here | 157 / 687 | 140 / 704 | 138 / 636 | Biological groups; the corrected fits use the 6 source controls of GSE130970 instead of its 23 F0 biopsies. 774-unit cohorts (control / disease): GSE126848 23 / 7, GSE130970 6 / 65, GSE135251 10 / 197, GSE162694 31 / 109, GSE213621 68 / 258 |
| Figure 3D PCA participants | 328 strict extremes | 418 | 416 | 392 | Endpoint visualization; the 774-unit roster is 113 strict controls and 279 advanced-disease units across four cohorts (844-record roster: 113 and 303); the pooled PCA is supplementary |
| Genes tested | 27,638 | 23,370 | 23,702 | 24,115 | Ensembl gene rows. Universe changes under corrected fragment counting and filters; 426 of the 23,702 rows lie on alternate contigs or patches. In the 774-unit fit none of the 855 held annotations is tested and 39 tested genes lie off chromosomes 1-22, X, Y and M |
| Canonical DEGs | 1,853 | 1,347 | 1,725 | 1,854 | Ensembl gene rows; `padj < 0.05` and `|log2FC| > 0.50` (canonical from 2026-08-12); 1,700 of the 1,725 and all 1,854 lie on chromosomes 1-22, X, Y or M |
| Up / down DEGs | 1,438 / 415 | 1,003 / 344 | 1,261 / 464 | 1,333 / 521 | Direction by disease-vs-control log fold change |
| interval-null DEGs (sensitivity arm) | 1,918 | 1,616 | 2,036 | 2,131 | `treat_fdr < 0.05` with interval null `|logFC| ≤ 0.25`; canonical 2026-06-29 to 2026-08-12 |

Adopted 774-unit source:
`RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/`
(`arm_deg_counts.tsv`, `primary_source_controls/deg_results.csv`). Model: limma
voom with quality weights, dataset + inferred sex + disease, BH over the 24,115
tested genes. The two sensitivity arms give 1,907 canonical DEGs with the 4
GSE130970 controls that have steatosis 0 and NAS 0, and 1,897 after removing the
three sex-discordant GSE130970 biopsies (771 units). Against the 2026-08-10 fit,
1,278 canonical DEGs are retained, 447 gained among the 22,978 shared genes, 128
gained among newly tested genes and 68 lost, and one retained gene changes sign;
log2FC Spearman on the shared genes is 0.958
(`primary_source_controls/old_vs_new_summary.tsv`, `old_vs_new_metrics.tsv`).
Against the 844-record corrected-controls fit, among the 23,121 genes tested in
both, 1,553 calls are shared, 126 are lost and 189 are new; over all genes 172
of the 1,725 earlier calls are absent (46 not tested) and 301 of the 1,854 are
new (112 not tested before). IGHA2 is called in both with opposite sign (-0.774
to +0.615). log2FC Pearson is 0.986 and Spearman 0.972. The sensitivity counting
arm, which credits a read shared with a held annotation to its single candidate
gene, changes no call among 24,108 commonly tested genes (log2FC Pearson 0.9999;
1,855 calls, the added one tested only in that arm). The 738-unit strict-depth
fit gives 1,906 calls; it is a four-cohort disease contrast because GSE126848
keeps 5 controls and no disease unit there
(`Analysis/Hypotheses/refit-adoption-review-20261011T151011Z/out/03_04_fit_comparisons.tsv`,
`04_call_counts_all_fits.tsv`).

Cohort dependence. GSE213621 supplies 326 of the 774 units, 68 of the 138
controls and 45.1% of the summed quality weight. Refitting without each cohort
in turn keeps 96.8% (GSE126848), 94.4% (GSE130970), 94.1% (GSE135251), 88.6%
(GSE162694) and 66.7% (GSE213621; 1,236 of 1,854) of the calls with the same
sign (`bulk-loco-summary-20261009T143029Z/loco_summary.tsv`). Without GSE213621,
618 calls are lost: 108 are no longer tested, and of the 510 still tested 505
keep their sign (251 keep BH padj < 0.05 with |log2FC| at or below 0.5, 43 keep
|log2FC| > 0.5 with padj at or above 0.05, 211 miss both) and 5 have the
opposite sign with padj above 0.48. Effects shrink: across the 1,746 tested
calls the regression slope of the leave-out log2FC on the full log2FC is 0.87
(Pearson 0.956, Spearman 0.859). 1,026 calls survive all five removals
(`Analysis/Hypotheses/refit-adoption-review-20261011T151011Z/out/05_loco_recomputed.tsv`). The spike-in calibration below (H81) connects these removals to the GSE213621
weight: an effect placed in one cohort only was called in the pooled fit in
87.45% of replicates when that cohort was GSE213621 and in at most 0.25% for any
other cohort (1.0 log2FC), and 99.8% of supported one-cohort effects (2.0
log2FC) lost their call in at least one LOCO refit. A uniform 1.0 log2FC effect
kept its call in all five refits in 90.6% of cases. Neither rate turns the 828
calls not kept in all five into source-dependent calls.

Glucocorticoid-response score as a covariate (V1 on the 774 units, 2026-10-12;
`Analysis/Hypotheses/V1-774-20261012T012943Z/results/run-20261012T013138Z/`:
`step0_reproduction.tsv`, `step2_score_by_cohort.tsv`,
`step2_score_disease_lm.tsv`, `step3_overlap_fisher.tsv`,
`step4_refit_summary.tsv`, `step4_concordant26.tsv`, `step4_gene_table.tsv.gz`,
`step5_random_summary.tsv`; job 22200894). This is a post hoc sensitivity
analysis outside the V1 Holm family; the earlier V1 prespecification said no
refit, and this run's `PRESPEC.md` states the conflict. No cohort has a
glucocorticoid exposure variable. The score is an expression proxy that may
carry disease biology: the mean, over the 70 dexamethasone-response genes V1
derived from GSE320593 (46 raised, 24 lowered), of sign times the within-cohort
z of voom log2 CPM. Model `~ dataset + inferred_sex + gc_score + group_binary`
with voomWithQualityWeights rerun, BH over 24,115 genes, padj < 0.05 and
|log2FC| > 0.5, seed 20260923, n = 774 genetic-identity clusters (138 control,
636 disease). The adopted fit was reproduced first (1,854 calls; largest log2FC
difference 4.5e-14).

| Model | Calls | Up | Down | Kept of 1,854 | Lost | Gained |
|---|---:|---:|---:|---:|---:|---:|
| Adopted fit | 1,854 | 1,333 | 521 | 1,854 | 0 | 0 |
| With 70-gene score | 1,900 | 1,416 | 484 | 1,732 (93.4%) | 122 | 168 |
| Score from the 46 raised genes only | 1,926 | 1,420 | 506 | 1,712 | 142 | 214 |
| 16-gene core list, recalled without a citation (descriptive only) | 1,628 | 1,202 | 426 | 1,552 | 302 | 76 |

With the 70-gene score, log2FC over all genes correlates with the adopted fit
at Pearson 0.9936 and Spearman 0.9896. Of the 122 lost calls, 119 fall to
|log2FC| at or below 0.5 with padj still under 0.05, 2 lose padj only, 1 loses
both, and none changes sign. All 26 genes in both the genetics and the
transcriptomics class keep their call. Control: 200 random 70-gene scores
matched on expression decile and sign pattern lose a median of 193 calls (2.5
to 97.5 percentile 29 to 753) and gain a median of 92.5 (20 to 593); the
observed 122 lost and 168 gained sit at empirical positions 0.62 and 0.26, so
the change is within what any expression-derived covariate of this kind
produces. The score is lower in disease than in control units by 0.162 (95% CI
0.084 to 0.240, nominal p 5.0e-05; cohort and sex adjusted); by cohort the
disease-minus-control difference is negative in four cohorts and +0.115 in
GSE162694. Of the 521 down calls, 9 are in the 70-gene set (odds ratio 6.78,
95% CI 2.94 to 13.8, nominal p 1.9e-05, background all 24,115 genes; 12.2 on
the 14,174-gene common background); of the 1,333 up calls, 4 are (odds ratio
1.04). The 9 are all genes that dexamethasone raises: FKBP5, KLF9, TAT, CEBPD,
SYBU, GADD45G, ERRFI1, CXCL2 and NR0B2. Eight of the 13 set genes that were
calls are lost with the score in the model (all of the 9 except NR0B2); FKBP5
goes from log2FC -0.587 to -0.340. Not done: no exposure data; surgical status
of GSE162694 and GSE213621 was not re-derived from source metadata; no
per-cohort covariate refit and none without GSE213621; no interval-null. Reading: the
disease-versus-control signal does not depend on the glucocorticoid-response
score. A small group of glucocorticoid-response genes (FKBP5, TAT, KLF9 and
others) is lower in disease biopsies than in controls, and those individual
calls should be treated as possibly reflecting differences in perioperative
exposure, not disease.

**GTEx postmortem liver fat (H55, 2026-10-15;
`Analysis/Hypotheses/H55-20261015T112558Z/RESULT_H55.md`).** Pass, same direction.
Among 1,002 disease-state calls measurable in GTEx liver, genes change with image-measured
histologic fat in 222 postmortem donors not selected for liver disease in the same
direction (rho 0.638) as between disease and control biopsies (Freedman-Lane two-sided
P 0.0002, floor 2/10,001; 10,000 permutations within six center-by-Hardy strata; seed
20260923). Also in donors who did not die on a ventilator (Hardy 1–4, 138 donors: rho
0.632, Holm P 0.0004). More than with ischemic time (joint model: rho 0.636 for fat and
-0.350 for ischemic time; difference 0.986, bootstrap 95% interval 0.625 to 1.022, Holm
P 0.0005). Exposure validity: Spearman 0.813 between the clear-vacuole fraction and the
pathologist note grade (112 donors). Controls: PLIN2 slope 1.069 log2 CPM (P 7.52e-07)
and the 197-gene steatosis composite (rho 0.717, P 0.0002) met their predictions; slide
tissue area (rho 0.354, P 0.0564) and ischemic time did not meet their negative-control
predictions, and neither stopped the analysis. 326 of the 1,002 genes have BH q < 0.05
for the fat slope, and 99.1% of these have the bulk sign. The claim is histologic
steatosis of unrecorded cause in postmortem liver; it does not assign any biopsy-cohort
call to steatosis, and attribution in the biopsy cohorts stays indeterminate. The
`gtex_postmortem_fat` annotation is a field of Catalog r9d (adopted 2026-10-16), not r8. No current
number changes.

**Adolescent cohort GSE66676 (H62, 2026-10-15;
`Analysis/Hypotheses/H62-20261015T124908Z/RESULT_H62.md`).** Null. In GSE66676, the
adult signature score did not differ detectably between 32 adolescents with NAFLD and 33
without (d 0.44, 95% CI -0.07 to 0.95); the cohort had 80% power for d of 0.71. The
median adult NAFL-versus-control shift (d 0.97) is outside this interval. Model
`score ~ nafld + sex + age + bmi` (residual df 60): coefficient 0.164 (95% CI -0.026 to
0.355), P 0.090, Freedman-Lane P 0.093 (10,000 permutations); family of one; seed
20260923. Two of 67 participants were removed for a recorded-versus-inferred sex
mismatch (S4). The IGFBP2 and IL32 positive control was not met (both moved in the adult
direction, neither at nominal P < 0.05: P 0.37 and 0.079), and the end-to-end permutation
calibration gave a false-positive rate of 0.067 against a predicted 0.035 to 0.065. The
controls are severely obese adolescents, so the result does not show
that the signature is specific to adults. GSE185051, descriptive only: in 52 children the
score correlated with NAS at rho 0.12, nominal P 0.41. No current claim changes.

**Per-gene composition sensitivity of the 1,854 calls (H56, 2026-10-15;
`Analysis/Hypotheses/H56-20261015T112556Z/RESULT_H56.md`).** Positive control failed
(S6): the 159 stromal-specific calls had a median adjusted-to-unadjusted log2FC ratio of
0.528 with the seven-group coding and the earlier fractions, where a median below 0.5 was
required in all eight flag cells. The calibration null would also have stopped the labels
(16-type coding, earlier fractions: lambda 1.254, above 1.2), and the uniform-call
negative control failed in 7 of 8 cells; the within-hepatocyte benchmark (NC2) also
failed. No `bulk_composition_sensitivity` Catalog column
is written. The Resource keeps describing the 1,854 calls as tissue-level effects. No
current number changes.

**Genetic ancestry and the 1,854 calls (H79, 2026-10-15;
`Analysis/Hypotheses/H79-20261015T125125Z/results/run-20261015T132031Z/RESULT_H79.md`).** Both tests reject (family F1 = {T1, T2}, Holm at
family-wise alpha 0.05; empirical p floors 1/201 and 1/101). Within the same
cohorts, controls include more units of African genetic similarity than disease
biopsies (31 of 138 against 37 of 636). Adding four genetic ancestry PCs from
RNA-called genotypes to the disease-versus-control model changed 233 of the
1,854 calls (lost 140, gained 93), more than within-cohort-and-group
permutations of the same PCs (median 120, 95th percentile 161; Holm p 0.0199).
64 genes met the ancestry-sensitive rule (expected false flags from the
reference: median 6, 95th percentile 17). Restricting to EUR-called units
changed 570 calls (lost 101, gained 469), more than random subsets of the same
cohort and group sizes (median 409.5, 95th percentile 501; Holm p 0.0199).

| Fit | Units | Calls (up / down) | Lost | Gained | Changed | Reference median (95th percentile) | Holm p |
|---|---:|---:|---:|---:|---:|---:|---:|
| Adopted, `~ dataset + inferred_sex + group_binary` | 774 | 1,854 (1,333 / 521) | 0 | 0 | 0 | | |
| T1: adopted + PC1 to PC4 | 774 | 1,807 (1,285 / 522) | 140 | 93 | 233 | 120 (161) | 0.0199 |
| T2: adopted model, EUR-called units only | 614 | 2,222 (1,630 / 592) | 101 | 469 | 570 | 409.5 (501) | 0.0199 |
| D2: T1 without 16 possibly mixed units (descriptive) | 758 | 1,838 (1,301 / 537) | 132 | 116 | 248 | not computed | |
| D3: adopted + PC1 only (descriptive) | 774 | 1,795 (1,288 / 507) | 147 | 88 | 235 | not computed | |

The changes are threshold crossings rather than changes in effect size. Of the
140 calls lost in T1, 127 fell to |log2FC| at or below 0.5 with padj still
below 0.05, 3 lost padj only, 10 lost both and none changed sign. Of the 625
calls with adopted |log2FC| below 0.6, 125 (0.200) were lost, against a
permuted median of 0.096 (2.5 to 97.5 percentile 0.040 to 0.187). Over the
1,854 calls, log2FC with the PCs correlates with the adopted log2FC at Pearson
0.9972 (slope 0.984, no intercept); EUR-only, Pearson 0.9926 and slope 1.069.
These agreement values are descriptive and were not tested. Of the 64 flagged
genes, 40 are adopted calls lost with the PCs and 24 are non-calls gained; the
flags carry no multiplicity correction and are a descriptive label, not a
per-gene test. GSTM1 is flagged: its PC1 coefficient is -9.19, its log2FC goes
from -1.440 (padj 0.0133) to -0.989 (padj 0.102), and mean log2 CPM is 4.496 in
controls of African and 2.667 in controls of European genetic similarity.
Among the 26 genes in both the genetics and the disease-state class, 24 keep
their call with the PCs (MLIP and ELAPOR2 are lost; MLIP is flagged) and 23 in
the EUR-only fit. A separate EUR-only rule (changed in the EUR-only fit and in
at most 4 of 100 subsets) flags another 64 genes, none shared with the PC
flags, against a median of 18 expected false flags (95th percentile 48.1).
Controls: the adopted fit was reproduced (largest log2FC difference 4.46e-14);
rs738409 G is less frequent in units of African than of European genetic
similarity within disease group (stratified exact CMH OR 0.451, 95% CI 0.270
to 0.728, p 0.000627); the true PCs carry within-cell expression signal (4,879
genes at BH q < 0.05 on the PC F-test in a cohort-by-group cell-means fit,
against a permuted median of 1.5 and maximum of 2,025); the permuted-PC median
of that count was 1.5 (NC-2, limit 50). Unit: 774 genetic-identity clusters
(138 control, 636 disease) in five cohorts; covariates dataset, inferred sex
and disease group, plus PC1 to PC4 in T1; per-gene BH over the fixed 24,115-gene
universe inside each fit is the call rule, not a test of H79; seed 20260923.
GSE135251 per-cohort membership in Figure 3C cannot be ancestry matched (only 2
EUR controls; 6 of its 10 controls are of African genetic similarity, against 171 of its 197 disease units of European similarity; `Analysis/Hypotheses/H79-20261015T125125Z/results/run-20261015T132031Z/s3_ancestry_by_cohort_group.tsv`). The PCs
come from expressed sites and may also carry library effects; the rs738409
control shows they track a known allele frequency but does not exclude a
technical component. The adopted fit remains the release and no main-figure
number changes; Catalog r9d (adopted 2026-10-16) carries the 64-gene
ancestry-sensitivity flag as a descriptive field.

**Calibration of the bulk evidence states by spike-in (H81, 2026-10-15;
`Analysis/Hypotheses/H81-20261015T125233Z/run-20261015T141242Z/RESULT_H81.md`).** Pass (T-O2 met, PC-A and PC-B pass, no stop rule).
In 50 replicates that added known log2 fold changes to the 774-unit counts
after permuting disease labels within cohort, the pooled fit called a uniform
1.0 log2FC change in 50.7 to 100.0% of genes per expression decile and called
0.0% of unchanged genes; the smallest change detected with 80% recall was
0.649 log2FC at its lowest (decile 8), 0.653 in the most highly expressed decile and
above 1.5 in the least. Changes were added by binomial thinning of control or disease counts
(Gerard 2020), each replicate was refit with the shipped pooled recipe, the
five leave-one-cohort-out (LOCO) refits and the five per-cohort fits, and 3,200
genes were spiked per replicate (320 per expression decile). There is no
hypothesis test and no family: every rate is an estimate with a 2,000-draw
replicate-cluster bootstrap interval (seed 20261923), and thresholds that
certify a low error rate are judged on the unfavourable bound.

| Expression decile (AveExpr, log2 CPM) | Recall at 0.5 | Recall at 1.0 (lower bound) | Recall at 1.5 | Detection limit at 80% recall (log2FC) | Equivalence rule usable |
|---|---:|---:|---:|---:|---|
| 1 (-4.74 to -2.74) | 0.181 | 0.507 (0.492) | 0.438 | above 1.5 | no |
| 2 (-2.74 to -2.02) | 0.372 | 0.796 (0.776) | 0.701 | 0.732 | no |
| 3 (-2.02 to -1.15) | 0.453 | 0.989 (0.984) | 0.946 | 0.666 | no |
| 4 (-1.15 to -0.06) | 0.479 | 1.000 (1.000) | 1.000 | 0.656 | no |
| 5 (-0.06 to 1.28) | 0.488 | 0.998 (0.996) | 0.999 | 0.654 | no |
| 6 (1.28 to 2.58) | 0.503 | 1.000 (0.999) | 0.999 | 0.650 | no |
| 7 (2.58 to 3.56) | 0.504 | 1.000 (1.000) | 1.000 | 0.650 | yes |
| 8 (3.56 to 4.44) | 0.508 | 0.999 (0.998) | 1.000 | 0.649 | yes |
| 9 (4.44 to 5.50) | 0.488 | 1.000 (1.000) | 1.000 | 0.653 | yes |
| 10 (5.50 to 16.72) | 0.488 | 1.000 (1.000) | 1.000 | 0.653 | no |

False `supported` calls among unchanged spiked genes: 0.0004 (bootstrap 95%
interval 0.0001 to 0.0008; largest decile upper bound 0.0040); a median of 5
unspiked genes per replicate were `supported`. Recall at 1.0 is below 0.8 in
deciles 1 and 2 by the lower bound, so a missing bulk call there does not
exclude a 2-fold disease-versus-control change; those genes are below the
detection limit rather than unchanged. They are 4,323 of the 22,261 tested
genes without an adopted call (2,165 in decile 1 and 2,158 in decile 2;
`summary/D5_detection_limit_lookup_not_supported_genes.tsv.gz`). Recall at 1.5
falls below recall at 1.0 in deciles 1 to 3 (PC-B(ii) not met) because
thinning removes more of these genes at `filterByExpr` (untested 0.559, 0.299
and 0.053 at 1.5 against 0.484, 0.200 and 0.009 at 1.0).

The unadopted L09-01 equivalence rule (random-effects TOST at margin 0.5 over
the five per-cohort fits) called a 0.75 log2FC change `tested_negative` in 0.0%
of genes in each of deciles 7, 8 and 9, the only deciles where it met all five
error limits. Elsewhere it cannot certify absence of a 0.75 log2FC change: its
size at the margin under between-cohort heterogeneity was 0.0579 overall
(bootstrap 0.0549 to 0.0610), and it called 43.1% of 1.0 log2FC effects
confined to GSE126848 `tested_negative`. The median BH-effective alpha in the
replicates was 0.0411 against 0.0177 on the adopted per-cohort fits (ratio
2.321). The Catalog `tested_negative` state remains not assessed (no negative-call
rule is defined); adopting the rule is a user decision that also needs the other L09-01 gate failures resolved.

| Spiked effect | Pooled-fit call rate by spiked cohort (GSE126848 / GSE130970 / GSE135251 / GSE162694 / GSE213621) | Calls that lose the call in at least one LOCO refit |
|---|---|---:|
| 1.0 log2FC in one cohort only | 0.0000 / 0.0005 / 0.0025 / 0.0020 / 0.8745 | 0.998 (1,755 of 1,759) |
| 2.0 log2FC in one cohort only | 0.0005 / 0.0005 / 0.0005 / 0.0315 / 0.9420 | 0.998 (1,947 of 1,950) |
| Uniform 1.0 log2FC | not applicable | 0.094 |
| Uniform 1.5 log2FC | not applicable | 0.090 |
| Uniform 0.75 log2FC | not applicable | 0.127 |
| Uniform 0.5 log2FC | not applicable | 0.540 |

In these simulations a supported effect confined to one cohort lost its call in
at least one LOCO refit in 99.8% of cases (2.0 log2FC; denominator from
GSE126848 1, GSE130970 1, GSE135251 1, GSE162694 63 and GSE213621 1,884
supported gene-replicates). This does not show that LOCO retention proves
cohort independence for real genes. A uniform effect of 1.0 or 1.5 log2FC kept
its call in all five LOCO refits in 90.6% and 91.0% of cases, so a LOCO loss for
a call of that size is uncommon under a uniform effect; LOCO losses concentrate
in calls near the 0.5 floor. No call is relabelled. Unit: spiked genes inside
replicates are calibration draws on one fixed set of 774 genetic-identity
clusters, not donors; every rate describes the recipe on these units and no p
value is reported.

Historical 844-record corrected-controls source:
`RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/`
(`arm_deg_counts.tsv`). The model is unchanged; only GSE130970 control labels
change (review item 1). The two sensitivity arms give 1,775 canonical DEGs with the
4 controls that have steatosis 0 and NAS 0, and 1,770 after removing the three
sex-discordant GSE130970 biopsies (841 participants). Against the 2026-08-10 fit,
1,326 canonical DEGs are retained, 366 gained among the 23,370 shared genes, 33
gained among newly passing genes and 21 lost; log2FC Spearman on the shared genes
is 0.983. The chromosome count joins the fit's gene rows to the chromosome column
of the corrected stage-extension table.

**Paper decision:** the five-cohort fragment-count analysis on the 774 units is
the bulk analysis of the paper (adopted by the user on 2026-10-11). Stage,
endpoint, lncRNA, evidence-class and continuum values below are quoted from
candidates built on it. Figure 3F, Figure 3G in the working layout, the Figure
4E program models and Figure 6 from Gene Catalog r8 were rebuilt as candidates
on 2026-10-12, as were the continuum composition-adjusted coefficients
(`figures/candidates/continuum-composition-774-20261012T012137Z/`). The Catalog
sidecar was rebuilt on the selected r8 inputs on October 14, as recorded
[below](#catalog-review-and-molecular-object-companion-2026-10-14); the October
12 statement that it was not rebuilt is superseded.

The validated endpoint sensitivity does not alter that estimand. Advanced
disease is `NAS >= 5` or fibrosis `>= F3`; strict controls are source-defined
controls with no observed NAS above 1 and no observed fibrosis above F0.
On the corrected labels no source label conflicts with histology. Requiring at
least five units per arm retains 392 of the 774 units (279 advanced, 113
strict controls, 6 of them from GSE130970) across four cohorts. The endpoint and
canonical log2 fold changes have Spearman rho 0.953 across all 24,115 genes and
90.9% direction agreement. The independently thresholded endpoint analysis calls
2,950 genes; 1,755 of 1,854 canonical genes are recovered (94.7%), with Jaccard
overlap 0.576 and 1,195 endpoint-only calls. The interval-null sensitivity
recovers 94.9% of canonical interval-null calls. Source:
`figures/candidates/fig3d-corrected-774-20261009T171017Z/source_tables/`
(`figs3_endpoint_deg_concordance_metrics.tsv`, `fig3d_endpoint_census.tsv`).
Historical 844-record values (2026-09-25): 416 participants (303 advanced, 113
strict controls), rho 0.940 across 23,702 genes, 89.2% direction agreement, 2,857
endpoint calls, 1,630 of 1,725 recovered (94.5%), Jaccard 0.552, 1,227
endpoint-only, interval-null recovery 92.2%. Superseded 2026-08-10 values: 418
participants (115 strict controls), rho 0.938, 1,307 of 1,347 recovered (97.0%),
1,336 endpoint-only. The endpoint analysis does not change canonical Catalog
membership.

**Figure 3C restored five-cohort results, 2026-10-14.** All five selected
cohorts now have estimable disease-versus-control models on the fixed
24,115-gene family, with BH FDR <0.05 and |log2FC| >0.5. Their union contains
13,761 distinct genes and 20,054 cohort memberships: 9,065 genes called in
one cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.
Counted by direction (review R03-01, 2026-10-15), 1,865 of the 4,696
multi-cohort genes (39.7%) are called up in one cohort and down in another
(1,529 / 286 / 43 / 7 in the two- to five-cohort tiers); 2,831 have one sign in
every calling cohort (1,933 / 636 / 218 / 44). 1,632 of the 1,865 include a
GSE126848 call; without GSE126848, 346 of 2,565 (13.5%) are mixed. Same-sign
tiers counted as the larger of up and down calls are 2,223 / 675 / 225 / 44.
Source: `figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-01_fig3c/`.
This is per-cohort membership, not the 1,854 pooled calls or a cohort-omission
analysis. The earlier four-fit 774-unit candidate (10,630 genes; 8,065 unique
to one cohort) and the historical 844-record censuses are superseded for this panel.

**Figure 3C with the GSE126848 filtered refit, 2026-10-21 (user decision D3).**
GSE126848 enters per-cohort membership through a refit after `filterByExpr` within
its 30 units (RLE factors recomputed; BH over its 15,592 genes; 3,440 calls, 1,821
up and 1,619 down). The other four cohorts keep the fixed 24,115-gene family. The
union contains 12,180 genes and 17,347 memberships: 8,341 genes called in one
cohort and 3,839 in two or more. 2,628 have one sign in every calling cohort
(1,820 / 582 / 186 / 40 in two to five cohorts) and 1,211 (31.5%) are mixed, 946 of
them with a GSE126848 call. Without GSE126848, 346 of 2,565 (13.5%) are mixed; with
the 6,147 calls, 1,865 of 4,696 (39.7%). All 3,258 refit calls that are also among
the 6,147 have the same sign in both fits, and the other 182 are not among them.
2,730 of the 6,147 calls (2,639 up) are on genes the cohort's own filter removes.
Sources:
`Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/`
(job 3121, seed 20260923) and `figures/candidates/fig3c-five-filtered-20261021T172748Z/`.
This replaces the 2026-10-14 union above for Figure 3C.

GSE126848 retains 23 control and 7 disease genetic-identity clusters, all male.
The source roster has 25 male controls and 31 disease records, including ten
women. Identity screening excluded 24 disease records (all ten women and 14 men)
and two male controls, all for `insufficient_filtered_identity_sites`.
The restored intercept-plus-disease model omits the constant sex term,
retains the established pooled expression filter and RLE normalization, and
re-estimates voom quality weights on its 30 units. It yields 6,147 calls in
its complete 24,115-gene BH family; 1,564 of 1,854 pooled calls agree in direction
(84.36%). An independent refit with separate group means reproduced all call
members; native selected fragment counts matched the adopted DGE exactly.
The other four cohort estimates are unchanged. Unequal identity eligibility
may limit representativeness; genetic clusters are not verified donors.
Library depth (review R03-02, 2026-10-15; jobs 22252035 and 22256115): retained
disease libraries have median 19.48M fragments counted on the 24,115 tested
genes against 35.30M for controls (disease deeper in 5 of 161 pairs; exact
Wilcoxon P=1.87e-5; r(group, log2 depth)=-0.709); after RLE normalization the
medians are 13.85M and 20.40M. 437 of the 6,147 calls (all up) are genes with
zero reads in all 30 libraries; each has log2FC 0.554, the precision-weighted
control-minus-disease difference in mean log2(normalized library size + 1).
With log2 library size added, the model gives 2,678 calls and keeps 2,025 of
the 6,147 (all-gene log2FC Spearman 0.799); with log2 normalized library size
it gives 3,649 and keeps 3,056 (Spearman 0.866). Neither keeps any of the 437.
A per-cohort filterByExpr refit gives 3,440 calls (1,821 up, 1,619 down) on
15,592 genes, a different BH family. Without the 437 zero-read calls the
Figure 3C union is 13,447 genes (8,843 unique to one cohort, 19,617
memberships). Median per-cohort quality weights are 0.970 (control) and 0.935
(disease); the 1.14 versus 0.72 medians are from the pooled 774-unit fit.
Depth sources: `figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-02_gse126848_depth/`
and `figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/repair/repair_libdef_log.txt`.
Sources: `figures/candidates/bulk-context-774-20261014T140904Z_22236311/`
(`cohort_disease_all_gene_results.tsv`, `native_GSE126848_checks.tsv`) and
`figures/candidates/remaining-774-panels-20261014T141620Z_22236524/`
(`fig3c_cohort_alluvial_tier_counts_from_table.tsv`).
Figure S3C plots the leave-one-cohort-out analysis on the fixed 24,115-gene set:
canonical-threshold retention ranges from 70.7% to 96.6% and all-gene effect
Spearman correlations from 0.870 to 0.995; per-cohort direction agreement among
the 1,854 canonical genes ranges from 78.3% to 100% over all five estimable
cohorts (`figures/candidates/figs3c-774-20261011T165944Z/source_tables/figs3c_cohort_robustness.tsv`).
The full refits, which repeat expression filtering after each removal, keep
66.7% to 96.8% of calls with the same sign (see Cohort dependence above); the
text quotes the full refits. Historical 844-record values: five cohorts, 15,540
genes, 23,292 memberships, 9,721 unique to one cohort, 4,326, 1,122, 302 and 69
in two to five cohorts; 77.8% to 100%; 0.846 to 0.992; 68.9% to 95.7%.
Superseded 2026-08-10 values: 15,090 genes, 9,777 unique to one cohort, 79.8% to
100%, 0.867 to 0.992, 60.7% to 98.0%. The genetic-identity cluster is the
inferential unit in every 774-unit fit.

## Recurring biological examples

These examples were chosen to explain biological questions, not because they
maximize modality count. The disease-state values below are from the
774-unit refit adopted on 2026-10-11. The pooled
RNA model uses 774 genetic-identity clusters and fixed cohort and inferred-sex terms; BH
adjustment covers all 24,115 tested genes. Adjacent-stage models use 205, 190, 157, and 94 units for F0→F1,
F1→F2, F2→F3, and F3→F4, respectively, with fixed cohort and inferred-sex
terms and a separate 24,115-gene BH family per contrast. PXD051911 protein
effects use 58 patients and adjust for acquisition batch, age, BMI, and sex;
protein q values use BH correction across all 6,071 tested proteins.

### Lower hepatocyte metabolic gene and protein abundance: GNMT–MAT1A–CYP2C19

| Gene | Five-cohort disease versus control | Adjacent fibrosis contrast(s) | Adjusted liver protein |
|---|---|---|---|
| GNMT | log2FC -0.859; canonical positive; interval-null q=6.87e-10 | F0→F1 -0.428, q=0.00913; q=0.0426 on disease-only F0 and 0.0430 with a control-status term, but q=0.207 without GSE130970. F1→F2 -0.382, q=0.0562, not supported at BH 0.05 on the 774 units | log2FC -1.136, BH q=0.00642 |
| MAT1A | log2FC -0.485; below the canonical 0.50 floor but interval-null q=3.51e-8 | F0→F1 -0.239, q=2.45e-4; q=9.5e-4 on disease-only F0 and 1.6e-3 with a control-status term, but q=0.0892 without GSE130970. No later adjacent contrast is supported | log2FC -0.473, BH q=1.22e-4 |
| CYP2C19 | log2FC -1.667; canonical and interval-null positive | F1→F2 -0.697, BH q=0.0477 (q=0.052 with a control-status term; q=0.219 without GSE135251 and 0.499 without GSE162694) | log2FC -0.459, BH q=0.0255 |

RNA and stage values: adopted 774-unit refit
(`bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv`)
and its stage tables
(`figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_all_gene_results.tsv`);
F0-arm q values from
`figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_hero_genes.tsv`.
Historical 844-record values: GNMT -0.896 (interval-null q=4.67e-12; F1→F2 -0.378,
q=0.0485; F0→F1 -0.424, q=0.0089), MAT1A -0.463 (F0→F1 -0.244, q=1.36e-4),
CYP2C19 -1.694 (F1→F2 -0.687, q=0.0414). The GNMT F1→F2 contrast was below BH
0.05 on the 844 records and is not on the 774 units.
Superseded 2026-08-10 RNA values: GNMT -0.813, MAT1A -0.433, CYP2C19 -1.438; the
earlier statement that the GNMT and MAT1A F0→F1 effects do not survive the F0
correction depended on the mislabeled GSE130970 controls and is withdrawn.

The RNA, stage, and protein directions agree for all three genes. This supports
a recurring example of reduced abundance of hepatocyte metabolic genes and proteins in
established disease. The protein model compares 46 MASLD with 12 control
participants and has no severity term. In a post hoc check that was not
prespecified (review round 3, F5bcd), protein abundance within the 46 MASLD
participants is negatively associated with NAS (rank-residual partial Spearman,
batch, age, BMI and sex fitted within the 46 participants, on quantile-normalized
log2 abundance as in Figure 5C; t approximation with 40 df: GNMT -0.431,
MAT1A -0.597, CYP2C19 -0.503; BH q 0.0066, 1.8e-4 and 0.0014 over 3 proteins x 2
histology scores) and with fibrosis stage for GNMT and MAT1A (-0.338 and -0.565;
q 0.034 and 3.0e-4) but not CYP2C19 (-0.241; q 0.12). Lower bulk abundance
can reflect a smaller hepatocyte share of the biopsy as well as lower per-cell
abundance. It does not establish one upstream mechanism, a discrete
patient subtype, longitudinal change, or an inherited regulatory program.
Their version-fixed single-cell program memberships are tested but nonsignificant,
and their spatial rows are coverage context rather than predeclared gene-level
inference.

An exploratory native DIA-MS audit in 46 MASLD participants further tests
whether the three proteins mark a coherent state. Rank-residual partial
Spearman correlations adjusted for acquisition batch, age, BMI, and sex are
0.472 for GNMT–MAT1A, 0.406 for GNMT–CYP2C19, and 0.472 for
MAT1A–CYP2C19. The 2026-09-29 participant-level correction candidate reproduces
these coefficients and accounts for fitted nuisance-model degrees of freedom:
three-pair BH q=0.00243, 0.00772 and 0.00243, respectively, replacing the earlier
0.00141, 0.00518 and 0.00141. All three remain supported after adding fibrosis stage, including
global BH correction over all 12 pair-by-adjustment tests, but none remains
supported after adding NAS or all four histology components. This candidate
audit supports a histologic-severity-linked metabolic state, not an autonomous
three-gene mechanism. It is not part of the adopted release.

## Where the three aspects and variant combinations currently stand

This subsection collects what the release already contains for the two
goal-level dimensions added to `docs/PAPER.md`, so the numbers stay in one place.
Neither dimension carries a claim.

**Aspect resolution.** Recorded fibrosis stage is the severity anchor. NAS enters
as a graded sum in the adjacent-NAS contrasts (0 to 1-2, 1-2 to 3-4, 3-4 to 5-8)
tabulated elsewhere in this document, and as an adjustment set in the candidate
proteomics audit: the GNMT-MAT1A-CYP2C19 state holds after adjusting for fibrosis
stage, including global correction over all 12 pair-by-adjustment tests, but does
not remain supported after adding NAS or all four histology components. That is a
recorded result about the components as covariates. The three aspects have not
been resolved as separate outcome axes anywhere in the release, so the evidence
state is `indeterminate` and no aspect-specific attribution is available.

**Variant combinations.** Fine-mapping partitions credible-set mass by
consequence class rather than assuming inherited variation is regulatory. One
two-variant local haplotype is characterised at GNMT, reported as an unresolved
local haplotype rather than a resolved causal variant. No systematic treatment of
variant combinations exists beyond that locus, and the colocalization design is
single-variant throughout.

**Chromatin.** Accessibility is evaluated across jointly testable peaks, with the
large majority resolving to `indeterminate` and a minority to `source_dependent`,
alongside the variant accessibility heatmap and the six-lineage accessibility
matrix. Chromatin supplies assay-native context for placing candidates; no
epigenetic effect or causal claim is available from it.

### Prespecified programs with different spatial organization: Stromal ECM and Ductular injury

The source-diagnosis correction completed October 14 excludes the two
GSE174748 disease explants, which are cirrhotic in the original [Filliol study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9949942/),
under the unchanged three-stage rule. The primary stage model
([H59](../Analysis/Hypotheses/H59-20261015T112432Z/RESULT_H59.md))
uses the 60 donors whose datasets could be ambient-corrected. Only GSE202379
and GSE244832 contain more than one stage, so they alone inform the slope
(fitted donors: hepatocytes 55, cholangiocytes 38, fibroblasts 43,
macrophages 25, T cells 53). It fits generalized least squares with a
separate residual variance per dataset and tests each program against 19,999
within-dataset stage permutations (seed 20260923), with BH over all 117
programs and the unchanged stability rule. No program is selected; the
smallest permutation q is 0.323375. Ductular injury has β=0.1365 (95% CI
0.0121 to 0.2609), permutation p=0.0367, q=0.323375; detection probability at
its earlier β=0.402086 is 0.999. Stromal ECM has β=0.0837 (−0.0672 to
0.2346), q=0.577711; detection probability at its earlier β=0.351525 is
0.794, so this analysis cannot exclude an effect of that size. Dataset slopes
differ for both (Ductular injury 0.1074 and 0.7423, Cochran Q p=0.0164;
Stromal ECM 0.0373 and 0.7556, Q p=0.00676). The 62-donor equal-variance fit
(Stromal ECM β=0.351525, q=0.085458, HC3 q=0.165240; Ductular injury
β=0.402086, q=0.037633, HC3 q=0.155826) is a sensitivity that depends on
equal variances, on its parametric null and on two GSE189600 donors without
ambient correction. It has zero calls with the same within-dataset
permutation null (Ductular injury permutation q=0.12285), without the two
GSE189600 donors (q=0.098523) or with a 50-cell lineage gate (q=0.076511).
The earlier 64-donor two-program selection (ECM q=0.0496, Ductular q=0.0187,
both HC3 q=0.105) is historical. Fixed weights and eligible donor scores are
unchanged; complete corrected tables and independent checks are in the
[source-diagnosis companion](../figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/).
IGFBP7 is not a canonical pooled DEG (log2FC 0.447, below the 0.50 floor) but increases at F2→F3
(0.416, BH q=0.00249). BICC1 is a canonical pooled DEG (log2FC 0.599) and
increases at F1→F2, F2→F3, and F3→F4 (0.426, 0.509, and 0.799; BH
q=0.0146, 0.0135, and 0.0496). These RNA values are from the adopted 774-unit
refit and its stage tables (844-record values: 0.364; 0.428, q 0.00108; 0.586;
0.420, 0.486 and 0.825, q 0.0134, 0.0142 and 0.0390). Matched spatial nulls support the Stromal ECM
program in GSE192741, whereas the Ductular injury program is indeterminate.
IGFBP7 is one of 49 fixed Stromal ECM genes and contributes 0.02474 (2.47%)
of its L1 weight. BICC1 is one of 22 Ductular injury genes and contributes
0.05436 (5.44%). These are program exemplars, not gene-level spatial tests or a
proposed IGFBP7–BICC1 molecular interaction.

The aligned ambient candidate freezes all 117 memberships, weights, and scoring
universes and uses the current source-diagnosis metadata. decontX succeeds for
GSE136103, GSE174748, GSE185477, GSE202379, and GSE244832. GSE189600 and Liver
Atlas fail because their ambient batches contain only one contamination
cluster; they are excluded from both sides of the primary common-universe
comparison and appear as `uncorrected_passthrough` only in the full-universe
sensitivity. Raw reconstructed scores pass r=0.99893–1.00000 for all 117
programs. Integer-count transport is valid for 106 programs. The 11 version-fixed
T-cell programs are `untestable` because their native non-integer scoring
substrate cannot be reproduced from complete-atlas integer counts.

In the corrected 60-donor primary common universe, with four contributing
datasets for these focal programs, Stromal ECM changes from β=0.350898
(BH q=0.127377; HC3 q=0.230021) to β=0.127891 (BH q=0.423983;
HC3 q=0.450123). Its paired correction-induced change is β=-0.223006,
a 63.6% reduction (nominal p=0.006417, the smallest ordinary p of 106
computable paired changes; HC3 p=0.024386, eighth by HC3 p, where
cholangiocyte Inflammatory / oxidative (SOD2) is first at 0.00307; BH
q=0.478071; HC3 q=0.337646), and its
median version-fixed L1-weighted ambient burden is 0.810 across corrected datasets
(range 0.689–0.929). Ductular injury changes from β=0.384613 (BH q=0.098523;
HC3 q=0.181719) to β=0.401035 (BH q=0.111792; HC3 q=0.194105).
Its paired change is β=0.016422 (BH q=1.00; HC3 q=0.999547), with median
ambient burden 0.342 (range 0.235–0.685). Neither program is
selected in the H59 primary fit, so the selection-conditioned attenuation
state is not applicable to either. Under the H59 model the ambient-corrected
betas are 0.0518 for Stromal ECM (permutation q=0.55887) and 0.1428 for
Ductular injury (q=0.26975). Neither paired change passes
the 117-program family correction. The earlier 62-donor ambient coefficients
0.401→0.130 and 0.396→0.389 are historical. The full-universe passthrough
sensitivity now has at most 62 donors, reproduces the corrected primary raw
effects, and is not described as fully ambient-corrected.

The separately validated same-atlas localization analysis pools runs to
biological donors and uses ten prespecified lineage-minus-hepatocyte contrasts.
Stromal ECM is higher in fibroblasts by β=2.271 (BH q=3.38×10^-26; HC3
q=6.88×10^-19; 53 paired donors, four datasets). Ductular injury is higher in
cholangiocytes by β=1.574 (BH q=2.04×10^-23; HC3 q=8.01×10^-21; 51 paired
donors, five datasets). Both remain the highest-scoring lineage at the 50- and
20-cell donor thresholds and retain direction in every estimable leave-one-dataset-out
fit. These results support same-atlas descriptive localization, not transcript
origin or cell autonomy.

The first test of the two fixed programs in new donors (V3,
GSE256398 single-nucleus hepatocyte pseudobulk, at least 50 nuclei per donor;
HC3; Holm over four tests) passes no test. Ductular injury: stage beta 0.694
(95% CI 0.062 to 1.326; n=13 non-cirrhotic donors; p=0.0347, Holm 0.139)
and disease versus healthy 1.080 (−0.091 to 2.250; n=17; p=0.0677, Holm
0.203). The disease group of that second contrast has 11 donors: 3 MASLD F0,
4 MASH F2-F3 and 4 MASH cirrhosis, so unlike the three-stage atlas models it
includes cirrhosis. Stromal ECM: 0.300 (−0.709 to 1.309) and 0.176 (−0.840
to 1.192), Holm 1. Power is low at these donor counts and no positive control
was prespecified or run ([V3 table](../Analysis/Hypotheses/V3-20261008T015239Z/results/v3_primary_family.tsv)).
V3's marker labels reproduce exactly only on Intel v4 compute nodes; on other
CPU types the clustering, and therefore the hepatocyte pseudobulk, differs. A
rerun on v4 nodes reproduced the V3 pseudobulk exactly, and under the H90
detection rule both programs remain testable in GSE256398 (Ductular injury 20
genes, 0.910 of weight; Stromal ECM 41 genes, 0.841)
([V0](../Analysis/Hypotheses/V0-20261021T122533Z/RESULT_V0.md)).
In an independent probe-based cohort (GSE298719; 6 normal, 7 F2-F4 MASLD
biopsies), the fixed Ductular injury program was higher in MASLD (beta 1.35 SD,
95% CI 0.29-2.42; Freedman-Lane permutation p = 0.006), robust to processing
plate, but not specific relative to matched random gene sets (86th percentile)
([V5](../Analysis/Hypotheses/V5-20261021T122723Z/RESULT_V5.md)). Stromal ECM
was untestable on that assay (0.424 of weight detected). No figure panel uses
this result.
GSE289173 (48 donors, single-nucleus RNA-seq) was opened to test whether the
Ductular injury score is higher in later-stage donors, but the test could not
be run because only 14 of 21 program genes, carrying 0.64 of the fixed
program weight, were detected in at least 1% of hepatocyte nuclei
(prespecified minimum 0.80; [H90](../Analysis/Hypotheses/H90-20261015T120337Z/RESULT_H90.md)).

The compact Figure 4–6 artwork candidate v8 contains eight individual one-page
PDFs and two manuscript-scale review proofs. It passed 99 checks covering exact
panel and source-table rosters, alternative-to-experiment routing, one-page
PDFs, Ghostscript readability, font subtypes, six-point plain typography,
`useDingbats=FALSE`, and compact title-free dimensions. The v8 record remains
candidate provenance. On 2026-08-17, the user approved a subsequent working
Figure 4 selection that adopts same-atlas localization as 4D and fixed-member
bulk transport as 4E; its independent render, checksum, PDF, font, and scope
checks passed in job 20872108. The user subsequently approved a shared-axis 4D
geometry that removes the redundant hepatocyte reference row without changing
the ten-contrast family or estimates; its byte, source, PDF, font, manifest, and
scope checks passed in job 20878079. The gradual 4E bulk line plot subsequently
removed its BH-supported-weight point-size encoding while preserving all eight
effects and its source table exactly; 20/20 adoption checks passed in job
20878757. Synchronized Resource-release adoption remains open.

### Superseded July COLOC cross-check and shared-variant anchors

The table below is from the **superseded July** Tier-1/2 COLOC release. It has
been superseded by the native r2 rederivation below. The SuSiE arm churned
substantially at the August adoption (147 tier-1/2 genes changed call), so no
July row here may be quoted as current. Coordinates are GRCh37
and denote
the top shared-variant anchor stored by COLOC; the parenthetical value is
`SNP.PP.H4`. Neither the coordinate nor `SNP.PP.H4` is a GWAS credible-set PIP
or proof that the variant regulates the named gene.

| Gene | Superseded July result | Shared-variant anchor | Interpretation |
|---|---|---|---|
| GNMT | SuSiE PP.H4=0.1356 for UKBB ALT; single-signal ABF PP.H4=0.7076 for direct EUR NAFLD | 6:42921489 (0.2125) for SuSiE; 6:42928758 (0.3462) for direct-NAFLD ABF | No primary SuSiE call; one direct-disease ABF sensitivity call |
| MAT1A | no finite SuSiE result; ABF PP.H4=0.07305 for EUR AST | 10:81143387 (0.7790) | No positive genetic call; high variant-level sharing does not rescue low locus-level PP.H4 |
| CYP2C19 | no finite SuSiE result; ABF PP.H4=0.03427 for direct EUR NAFLD | 10:96306509 (0.05238) | No positive genetic call |
| IGFBP7 | no finite SuSiE result; ABF PP.H4=0.08484 for EUR AST; direct-MASLD ABF PP.H4=0.01333 | 4:57949517 (0.7396 for AST; 0.1245 for direct MASLD) | No positive genetic call; proxy and direct phenotypes remain separate |
| BICC1 | no finite SuSiE result; ABF PP.H4=0.02493 for EAS NAFLD | 10:60583832 (0.06040) | No positive genetic call |

The current exact rows were rederived from the adopted COLOC r2 native output
in job 22236878. Coordinates are GRCh37 COLOC shared-variant anchors, not GWAS
credible-set PIP or causal variants. None of these five genes is a primary SuSiE
call; missing finite SuSiE output for the other four genes is not a negative test.

| Gene | Current strongest ABF gene–trait PP.H4 | Shared-variant anchor (SNP.PP.H4) |
|---|---|---|
| GNMT | 0.705832, direct EUR NAFLD | 6:42931261, rs2296804 (0.258990) |
| MAT1A | 0.070472, MVP EUR AST | 10:81143387 (0.765392) |
| CYP2C19 | 0.037590, direct EUR NAFLD | 10:96292988 (0.044813) |
| IGFBP7 | 0.082400, UKBB AST | 4:57949517 (0.712220) |
| BICC1 | 0.025806, MVP EAS NAFLD | 10:60583832 (0.053429) |

GNMT's strongest SuSiE pair is UKBB ALT, PP.H4 0.128219, anchor 6:42921489
(SNP.PP.H4 0.182316). Its direct-NAFLD ABF pair is sensitivity evidence and
does not establish a genetically regulated three-gene program. Exact source
rows are in `figures/main/fig2_genetics/source_tables/recurring_example_exact_pairs.tsv`.

The native GNMT allele and LD source chain is also rechecked. At nonpalindromic
rs2296805, the G allele lowers liver GNMT expression and aligns with higher
NAFLD risk. This is a proxy for the current rs2296804 anchor: their UKBB-European
r2 is 0.998001. The anchor is C/G and near allele frequency 0.5, so its reported
allele labels do not independently resolve strand. rs9471976 is nonpalindromic
and also aligns lower expression with higher risk; its r2 with rs2296805 is
0.818301. rs11752813 has r2 0.580492 but is another palindromic C/G variant;
its sign is not independently harmonized. Unsigned LD supports shared local
variation and does not establish signed phasing or identify one causal nucleotide.
Ji et al. associated rs9471976 and rs11752813 with hepatic GNMT protein in 268
human liver biopsies and reported allele-dependent reporter activity. This is
source-owned functional context, not a Resource discovery
([Ji et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3463826/)). Native
risk/eQTL/LD rows are in `figures/main/fig2_genetics/source_tables/gnmt_risk_eqtl_ld_source_chain.tsv`.

The separate target-gene biotype companion uses all 436 primary SuSiE genes:
424 protein-coding and 12 lncRNA, with every stable ID resolved from native
liver-eQTL marginal tables. These are gene-level targets; the assay does not
identify transcript isoforms. Counts describe supported genes, not enrichment
relative to an unmeasured biotype-specific assay denominator.

## Noncoding RNA source and annotation candidate

These products are validated candidates, not an adopted release. The original
requirement that every library have a reverse/same assignment ratio of at least
four was retained as a failed sensitivity. Before opening disease outcomes, it
was replaced by a cohort-level rule requiring a median ratio of at least four,
at least 90% of libraries at or above four, and every library above one.
GSE213621 additionally required no biological group to exceed a 10% failure
rate and no more than five percentage points of failure-rate spread among
groups. All five cohorts passed that amended rule without excluding a library.

| Check | Result | Interpretation |
|---|---:|---|
| Count-source rule | 5 of 5 cohorts pass | All use featureCounts 2.1.1, GENCODE v49, and reverse-stranded counting |
| GSE126848 reverse/same ratio | median 18.52; 0 of 56 below 4 | Pass |
| GSE130970 reverse/same ratio | median 16.99; 0 of 78 below 4 | Pass |
| GSE135251 reverse/same ratio | median 15.30; 0 of 216 below 4 | Pass |
| GSE162694 ten-BAM audit | median 13.97; 0 of 10 below 4 | Pass; independently counted with `-s 0`, `-s 1`, and `-s 2` |
| GSE213621 reverse/same ratio | median 10.98; 346 of 367 at or above 4; all above 1 | Passes the amended cohort rule |

The 21 lower-ratio GSE213621 libraries are distributed across control, F0/F1,
F2, and F3/F4 groups rather than confined to one stage. Group-specific failure
counts are 5/68, 5/97, 7/107, and 4/95, respectively; the failure-rate spread is
3.14 percentage points. Every library still has more reverse- than same-strand
assignment. Both the failed strict rule and the passed amended rule are stored
with their hashes and rationales.

The independent GENCODE and prespecified-program annotation lanes are valid:

| Candidate quantity | Value | Qualification |
|---|---:|---|
| GENCODE v49 gene records | 86,369 | Versioned Ensembl IDs are primary |
| GENCODE v49 lncRNA genes | 36,490 | Includes alternative loci for coverage audit |
| Canonical-chromosome lncRNA genes | 34,880 | Chromosomes 1–22, X, and Y |
| Main-text-eligible lncRNA mappings | 33,991 | Canonical chromosome and neither same-strand exonic nor complex |
| Prespecified programs annotated | 117 | Memberships, hashes, and weights unchanged |
| Uniquely mapped lncRNA memberships | 409 | Membership rows, not distinct genes |
| Programs requiring leave-all-lncRNA-out sensitivity | 49 | At least 5% of original L1 weight is uniquely mapped lncRNA |

The seven lncRNA genomic classes contain 20,085 intergenic, 7,361 antisense
exonic-overlap, 4,591 antisense gene-body-overlap, 2,037 same-strand gene-body,
1,450 opposite-strand promoter-proximal, 632 same-strand exonic-overlap, and
334 complex genes. Classification is genomic context, not evidence of lncRNA
function.

### Five-cohort pooled biotype subdivision

**Adopted 774-unit fit, Figure 3 lncRNA inset rebuilt 2026-10-09.**
Source: `figures/candidates/fig3-lncrna-corrected-controls-r8-caliper-search-20261009T173834Z/source_tables/fig3_lncrna_role/`
(`summary.tsv`, `match_balance.tsv`). lncRNAs are 26.7% of the 24,115 tested
genes and 29.6% of the 1,854 canonical positives. Of canonical positives, 960 of
1,157 protein-coding genes (83.0%) and 263 of 549 lncRNAs (47.9%) are
upregulated in MASLD; on the 844 records the lncRNA figure was 56.3%, so most
called lncRNAs are now lower in disease. The unadjusted odds of a canonical
positive call, lncRNA versus protein-coding, is 1.12 (95% CI 1.00–1.25);
adjusted for expression and detection it is **0.428** (0.376–0.488; 21,501
genes). The high-confidence disease-associated lncRNA set is **392**. In 1:1
covariate-matched lncRNA/protein-coding pairs the odds of upregulation, lncRNA
versus coding, is **0.132** (0.090–0.193; 520 pairs; discordant pairs 30
lncRNA-up against 227 coding-up). The 0.20-SD caliper failed the 0.10
standardized-difference limit (gene length −0.108; 538 pairs) and 0.15 SD was
selected (largest difference 0.082, gene length).

**Historical 844-record corrected-controls fit, lncRNA inset of 2026-09-25.**
Source: `figures/candidates/fig3-lncrna-corrected-controls-20260925T195928Z/source_tables/fig3_lncrna_role/`
(`summary.tsv`, `match_balance.tsv`), from the lncRNA producers rerun in
`RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-corrected-controls-20260923T231956Z/`.
lncRNAs are 27.0% of the 23,702 tested genes and 30.9% of the 1,725 canonical
positives. Of canonical positives, 895 of 1,090 protein-coding genes (82.1%) and
300 of 533 lncRNAs (56.3%) are upregulated in MASLD. The unadjusted odds of a
canonical positive call, lncRNA versus protein-coding, is 1.14 (95% CI
1.02–1.27); adjusted for expression and detection it is **0.507** (0.447–0.575;
21,155 genes). The high-confidence disease-associated lncRNA set is **358**. In
1:1 covariate-matched lncRNA/protein-coding pairs the odds of upregulation,
lncRNA versus coding, is **0.210** (0.150–0.294; 499 pairs; discordant pairs 41
lncRNA-up against 195 coding-up). The match caliper is chosen from the covariates
only, before the direction outcome is computed: calipers of 0.20, 0.15, 0.10 and
0.05 SD are tried in order and the first whose standardized differences all meet
the unchanged 0.10 limit is kept. On the corrected fit the 0.20-SD match failed
that limit (gene length −0.157, expression variability 0.129; 528 pairs), and
0.10 SD was selected (largest difference 0.099, gene length). The values below
are the 2026-08-10 fit (23,370 genes, 1,347 canonical positives) and are
superseded: adjusted OR 0.439, 279 high-confidence lncRNAs, and matched direction
OR 0.203 over 368 pairs.

The primary inference remains one all-gene family under the canonical threshold. The
biotype rows are subdivisions of its 23,370 tested genes and 1,347 positives,
not a separate lncRNA-only multiple-testing analysis. One effect-size floor
applies to every biotype; see the note below on why a biotype-specific floor was
rejected.

| Biotype | Tested | Canonical positive | Up | Down | interval-null comparator | Up | Down |
|---|---:|---:|---:|---:|---:|---:|---:|
| Protein coding | 14,664 | 886 | 744 | 142 | 1,137 | 918 | 219 |
| lncRNA | 6,249 | 402 | 220 | 182 | 434 | 205 | 229 |
| Other | 2,457 | 59 | 39 | 20 | 45 | 21 | 24 |

**The unadjusted lncRNA share must not be read as enrichment.** lncRNAs are
26.7% of tested genes and 29.8% of canonical positives, against 26.9% under the
interval-null comparator. That apparent rise is a property of the threshold, not of biology:
lncRNAs are lower-expressed, so their fold-change estimates are noisier, and a
floor applied to the point estimate preferentially includes them. The interval-null
comparator, which weights by standard error, holds the share at parity (ratio
1.004). Lowering the floor makes the artifact worse, not better — the lncRNA
share of positives rises to 31.7% at a 0.30 floor and 32.1% at 0.20.

A biotype-specific effect-size floor was considered and **rejected**. It would
amplify the same artifact, and a class-specific filter makes the biotype
comparison circular by construction: the class given the looser filter is
enriched automatically.

The adjusted analysis is the interpretable one, and it is unchanged in direction
and strengthened in magnitude. In a secondary descriptive gene-level logistic
model adjusted for average abundance, expression variability, gene length, and
cohort detection, lncRNA annotation was associated with **lower** odds of a
positive call than protein-coding annotation: odds ratio **0.439** (95% CI
0.381–0.505; nominal P=1.13×10⁻³⁰; 20,913 genes) under the canonical threshold,
against 0.488 (95% CI 0.427–0.557; nominal P=4.38×10⁻²⁶) under the interval-null
comparator. Once abundance and variability are accounted for, the apparent
unadjusted enrichment reverses. This describes the tested gene universe; genes
are not biological replicates, and the result is not a donor-level
disease-effect test.

Regenerated under the canonical threshold on 2026-08-12
(`bulk_disease/run_canonical_gate_2026-08-12/`), the high-confidence lncRNA set
is **279**, against 271 in the retained interval-null arm. The count rises even though
the lncRNA positives fall from 434 to 402, so the high-confidence rate goes from
62.4% to 69.4% of positives. This is a property of the threshold: a hard floor on the
point estimate selects large-effect genes, and large effects re-clear the same
floor in the leave-one-cohort-out and equal-weight refits more reliably than the
modest-but-precise genes the interval null included. Both arms are retained
side by side; neither directory overwrote the other.

Among the 434 interval-null-positive lncRNAs, **271** meet every prespecified
high-confidence rule: eligible mapping, same direction in at least four of five
cohorts without a materially opposite fifth cohort, retained direction in all
five leave-one-cohort-out fits, retained interval-null significance in at least four,
and agreement under equal library weighting. The deterministic main-text
example rule selects `ENSG00000213062.6`, an unambiguous opposite-strand
promoter-proximal lncRNA record with no separate GENCODE display symbol
(log2 fold change -0.761; pooled t=-16.55). It must be displayed by stable ID;
no familiar name may be substituted post hoc.

Three one-page supplementary candidates report the complete robustness family,
genomic-class counts, and the deterministic stable-ID example. Each PDF passed
source-table checksum, one-page, Ghostscript, and visual review. They remain
candidate supplement panels and do not add a sixth main figure.

### DEG-threshold sensitivity

The canonical threshold and its alternatives are all reported. Values below are on
the current official substrate and are reproduced by
`figures/supplementary/figS_deg_threshold_sensitivity/`.

⚠ **Every number in this section is on the superseded July COLOC release.** The
whole analysis was rerun on the adopted 2026-08-17 release
(`RNA-seq/results/genetics_state_interface/20260817T133951Z_promoted_coloc/`).
The conclusions survive; the values move:

| | July (below) | adopted 2026-08-17 |
|---|---:|---:|
| jointly testable background | 14,931 | **14,920** |
| colocalized genes in it | 447 | **437** |
| observed / expected overlap, canonical threshold | 23 / 32 | **20 / 31.2** |
| overlap odds ratio, canonical threshold | 0.690 (p = 0.094) | **0.616 (p = 0.037)** |
| matched ladder at the 0.50 floor | 23 against 24 controls | **20 against 16**, McNemar OR 1.31, p = 0.58 |
| threshold-free SuSiE rho | 0.019 (p = 0.232) | **0.016 (p = 0.34)** |

Corrected inputs (corrected COLOC and corrected-controls bulk, interface universe;
`RNA-seq/results/genetics_state_interface/20260924T061535Z_coloc_r2_corrected_controls/`):
jointly testable 14,183; genetic 406; observed / expected overlap at the canonical
threshold 23 / 30.8, OR 0.72 (nominal p 0.15); matched ladder at the 0.50 floor
23 against 17 controls, McNemar OR 1.35, p 0.43; threshold-free SuSiE rho on raw
log2FC 0.021 (p 0.20) over 3,622 genes. The matched ladder stays flat (every
floor p ≥ 0.25).

The ladder is still flat at every floor and the threshold-free estimate is still
null, so **the "do not call colocalized genes buffered or constrained" rule
stands**, and the lower nominal odds ratio does not reverse it; it is not a depletion
finding (H57: no depletion is detected among the 3,769 genes evaluated by
multi-signal COLOC). Quote
the corrected values; the July numbers below are kept because the corrections
they document are why the analysis is trusted at all.

All overlap tests below run against the **14,931-gene jointly testable
background** — genes assayed in both the bulk contrast and colocalization — and
the 447 colocalized genes within it. DEG counts are given both genome-wide and
restricted to that background; only the restricted count enters the test.

| Threshold | DEGs (genome-wide) | DEGs (in background) | Non-overlap with the genetic set | Overlap odds ratio | Fisher p |
|---|---:|---:|---:|---:|---:|
| `padj < 0.05` and `|log2FC| > 0.50` (**canonical**) | 1,853 | 1,079 | 94.85% | 0.690 | 0.094 |
| Interval-null FDR < 0.05 at `lfc = 0.25` (calibrated comparator) | 1,918 | 1,261 | 92.39% | 0.889 | 0.604 |
| `padj < 0.05`, no fold-change floor | 13,043 | 8,435 | 41.16% | 1.104 | 0.333 |

The canonical and comparator thresholds agree on 1,497 genes; 356 are
canonical-only and 421 are interval-null-only (Jaccard 0.658). The two sets are
therefore not interchangeable even though their totals are close, and any
figure must state which one it used.

**Corrected 2026-08-13: the earlier ladder used the wrong background.**
`scripts/figures/figS_deg_threshold_sensitivity_data.py` set the Fisher
background to all 27,638 bulk-tested genes rather than the 14,931 jointly
testable ones. A gene never tested for colocalization cannot count as evidence
against overlap; including such genes inflates the double-negative cell and
therefore every odds ratio. This is the asymmetric-denominator trap. The
producer now asserts against `orthogonality_audit.tsv`, which was always
correct, so the two can no longer disagree.

Two conclusions from the retired version do not survive:

- The claimed significant *enrichment* at permissive thresholds was an artifact.
  No floor: 1.612, p = 7.8e-07 → **1.104, p = 0.333**. Floor 0.10: 1.296,
  p = 0.0075 → **1.082, p = 0.412**. Neither is significant.
- The comparator threshold no longer points the other way. Its odds ratio was
  reported as 1.106 (enrichment) against the audit's 0.889 (depletion) for the
  same genes; **0.889 is correct**. The threshold choice never reversed the
  direction — the denominator did.

**The 0.50 floor is no longer load-bearing.** On the corrected background the
odds ratio decreases monotonically across all eleven floors, from 1.104 with no
floor to 0.395 at `|log2FC| > 1.00`, and is never significantly above 1 at any
floor. Depletion reaches nominal significance only at floors 0.25 and 0.30
(0.778 and 0.736; p = 0.046 and 0.033), neither surviving correction across the
eleven tests. The direction of the result is therefore stable at every
threshold, and the earlier instruction to always state the conditioning floor
is superseded — though the magnitude still depends on it and should be quoted
with the floor attached.

Across conventional floors from 0.30 to 1.00 the non-overlap fraction runs
87.0% to 99.6%. The corrected reading is that colocalized genes are
progressively **less** represented among genes with larger disease-state
responses, with no significant over-representation among weakly responding
genes at any floor.

**The ladder is an expression artifact, not a biological property.** Tested
directly: `scripts/manuscript/response_amplitude_matched.py`, results under
`RNA-seq/results/response_amplitude/`.

Colocalized genes are markedly higher-expressed than the rest of the background
(mean `AveExpr` 3.696 against 2.590, SMD 0.363), and a high `|log2FC|` floor
preferentially includes low-expression genes whose effect sizes are inflated
conditional on passing FDR. Matching each of the 447 colocalized genes to an
expression-matched control (L2 propensity, 0.2 SD logit caliper, 447 of 447
matched, max |SMD| 0.00002) removes the entire effect:

| Test | Result |
|---|---|
| Mean amplitude, `\|shrunk_logFC\|` | Δ = **+0.003**, 95% CI [−0.019, +0.025], Wilcoxon p = 0.64 |
| Mean amplitude, raw `\|logFC\|` | Δ = +0.007, 95% CI [−0.016, +0.030], p = 0.38 |
| Adding `SE` to the match (conservative) | Δ = −0.009, 95% CI [−0.036, +0.016], p = 0.71 |
| Decile-matched null, 10,000 draws | observed 0.1416 against null 0.1369 ± 0.0079, p = 0.556 |

The mean is not the estimand the ladder makes, so the ladder itself was
recomputed **within matched pairs**, scoring each colocalized gene and its
matched control DE or not at each floor and testing the discordant pairs with
exact McNemar. The matched ladder is flat: odds ratios 1.01, 1.17, 1.03, 1.11,
1.11, 1.07, 1.03, 0.95, 0.87, 1.14, 0.67 across the eleven floors, every p ≥
0.27, against the unmatched decline from 1.104 to 0.395. At the canonical 0.50
floor, 23 colocalized genes are strongly differentially expressed against
**24** of their expression-matched controls.

Conclusion: no difference was detected in these expression-matched comparisons;
equivalence was not tested. **Do not describe colocalized genes as buffered, constrained, or transcriptionally
held constant**, and do not present the ladder's decline as a property of
inherited risk. The non-overlap is what independent sampling of the
transcriptome produces for two sets of this size and expression profile —
447 genetic and 1,079 disease-state genes in a 14,931-gene background give an
expectation of 32 shared genes against 23 observed — not evidence that the two
layers occupy structurally different parts of the transcriptome.

The threshold-free analysis on the corrected COLOC release and the corrected-controls
fragment-count fit (`figures/candidates/fig3g-coloc-r2-20260924T061535Z/`, interface
tables in `RNA-seq/results/genetics_state_interface/20260924T061535Z_coloc_r2_corrected_controls/`)
uses 3,622 genes with finite SuSiE PP.H4 and raw RNA-seq log2FC estimates; the
signed association is Spearman rho=0.021 (P=0.20). Figure 3G shows these
observations without a fitted line and colours them by the canonical rule: 23
Both, 180 Transcriptomics sig. only and 383 Genetics sig. only. Superseded
2026-08-17 values: 3,617 genes, rho=0.0141 (P=0.396). These 3,622-gene values
are on the 844-record fit. Figure 3G on the 774-unit refit in the working layout
(`figures/candidates/fig3g-working-layout-774-20261012T011006Z/`, `checks.txt`
and `source_tables/fig3g_raw_treat_interface.tsv`; installed 2026-10-14)
plots 3,776 points: Neither 3,156, Transcriptomics sig. 194, Genetics sig. 400
and Both 26 (3,036, 180, 383 and 23 of 3,622 before), among 14,704 jointly
testable genes with 426 genetic. Because PP.H4 is unsigned, the test is its
Spearman rho with |log2FC|: 0.0014, nominal P=0.93 (0.0010, P=0.95 among the
3,769 genes of the H57 universe; review R14-07;
`figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R14-07_fig3g/fig3g_magnitude_statistics.tsv`
and `figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/repair/R14-07_fig3g_universes.tsv`). The signed-log2FC rho 0.022
(nominal P=0.18;
`fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_matched_treat_summary.tsv`)
is withdrawn as the wrong form for an unsigned predictor.
COLOC release `susie_coloc_r2_20260924T053244Z` is selected on 2026-10-14.

The current S3G expression-matched interval-null sensitivity uses 426 pairs:
37 calls only in the genetic member, 28 only in its matched control, 360 in
neither and one in both (paired OR=1.32143; nominal exact McNemar
p=0.321084). It is distinct from the canonical overlap and from the
hepatocyte-specificity class comparison. The native joint universe is 14,704;
annotation-adjudicated Catalog eligibility is reported separately.

Figure S3G contains the expression-matched sensitivity using the calibrated interval-null family (lfc=0.25, BH
q<0.05), not a post-FDR raw fold-change cutoff. In the historical 844-record
comparison, all 406 Genetics sig. genes
were paired one-to-one with expression-matched controls. Thirty-three pairs contain a
interval-null call only in the Genetics sig. member, 25 only in the matched control,
347 in neither, and 1 in both (paired OR=1.32; exact McNemar P=0.36).
Superseded 2026-08-17 values: 428 pairs, 30 / 21 / 377 / 0, OR 1.43, P 0.262.
This historical 844-record matched comparison detects no difference in call
frequency; it does not establish equivalence or biological complementarity.
Because PP.H4 is unsigned, it does not test risk-allele direction.

**Caveat that belongs in Methods.** Because the effect-size filter is applied
after FDR adjustment rather than folded into the test, the adjusted p-value
does not carry a calibrated error rate for the filtered set (McCarthy & Smyth
2009). The interval-null arm is retained as the calibrated comparator for exactly this
reason and its columns must not be deleted from `canonical_deg_results.csv`.

## Cross-sectional stage remodeling candidate

Historical job `19966029` completed from the 844-record, five-cohort
fragment-count substrate. All four adjacent fibrosis contrasts, three adjacent
NAS-group contrasts, and F1–F4-versus-F0 models use fixed cohort and inferred-sex
terms and complete 23,370-gene Benjamini–Hochberg families. The corrected-controls
rerun (`figures/candidates/bulk-release-corrected-controls-20260925T175215Z/`) uses
the same design with complete 23,702-gene families; the table below gives its
values. Input/output hashes,
family cardinalities, design rank, and reported FDR values passed independent
checks. Candidate v2 is validated but not adopted; job `19957749` and candidate
v1 remain failed execution provenance.

**774-unit stage tables (2026-10-08; built on the refit adopted 2026-10-11).**
`figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/`
uses the same design with complete 24,115-gene families
(`stage_extension_design_audit.tsv`). The 844-record tables further down this
section are historical; the 774-unit values are:

| Axis | Contrast | Reference / comparison units | Cohorts | Higher / lower at BH q<0.05 |
|---|---|---:|---:|---:|
| Adjacent fibrosis | F0→F1 | 103 / 102 | 3 | 645 / 693 |
| Adjacent fibrosis | F1→F2 | 102 / 88 | 3 | 261 / 137 |
| Adjacent fibrosis | F2→F3 | 88 / 69 | 3 | 381 / 131 |
| Adjacent fibrosis | F3→F4 | 69 / 25 | 3 | 184 / 97 |
| Fibrosis versus F0 | F1 vs F0 | 103 / 102 | 3 | 645 / 693 |
| Fibrosis versus F0 | F2 vs F0 | 103 / 88 | 3 | 2,155 / 2,111 |
| Fibrosis versus F0 | F3 vs F0 | 103 / 69 | 3 | 3,144 / 3,046 |
| Fibrosis versus F0 | F4 vs F0 | 103 / 25 | 3 | 3,185 / 3,858 |
| Adjacent NAS | 0→1–2 | 46 / 65 | 3 | 2,790 / 2,548 |
| Adjacent NAS | 1–2→3–4 | 65 / 117 | 3 | 2 / 1 |
| Adjacent NAS | 3–4→5–8 | 117 / 167 | 3 | 1,954 / 2,351 |
| NAS versus NAS0 | 3–4 vs 0 | 46 / 117 | 3 | 4,218 / 5,119 |
| NAS versus NAS0 | 5–8 vs 0 | 46 / 167 | 3 | 4,571 / 5,270 |

NAS0 is 45 source controls and 1 GSE162694 disease unit (GSE130970 4, GSE135251 10 and GSE162694 31 controls), so the versus-NAS0 rows and the NAS 0→1–2 step mainly contrast disease with controls. NAS 1–2 has 63 disease units and 2 GSE130970 controls; NAS 3–4 (117) and 5–8 (167) are all disease. NAS 1–2→3–4 (NAS 1–2 includes the 2 GSE130970 controls) gives 3 genes, and among disease biopsies only NAS 3–4→5–8 gives 4,305 (1,954 higher, 2,351 lower).

Counts are BH q<0.05 split by log2FC sign in `stage_extension_all_gene_results.tsv`,
with no fold-change filter. Figure 3E's F1–F4-versus-F0 totals are 1,338, 4,266,
6,190 and 7,043 under this rule, which differs from the 1,854 canonical
disease-versus-control calls requiring BH padj <0.05 and |log2FC| >0.5;
the same totals are in `figures/candidates/fig3e-s3-774-20261011T165944Z/source_tables/fig3e_deg_counts.tsv`. All fibrosis contrasts use GSE130970, GSE135251 and GSE162694 only (drawn from the 387 staged units; each contrast fits 94 to 205 of them); GSE213621 stage is binned and not used in these contrasts. Refitted without GSE130970 / GSE135251 / GSE162694 (review repair job 7550, same design; the base fits reproduce the release exactly), the BH q<0.05 counts are F0→F1 14 / 929 / 651 (155 / 115 / 140 units), F1→F2 240 / 205 / 0 (154 / 93 / 133), F2→F3 464 / 2 / 280 (138 / 54 / 122), F3→F4 102 / 6,565 / 508 (82 / 29 / 77), F2 vs F0 3,293 / 4,105 / 1,934, F3 vs F0 4,804 / 3,994 / 4,301 and F4 vs F0 4,800 / 9,546 / 4,830. All 14 F0→F1 calls without GSE130970 are among the 1,338, and 4,596 of the 7,043 F4 vs F0 calls remain without it. F0→F1 fitted in each cohort alone gives 252 / 3 / 0 (50 / 90 / 65 units). Of the 1,338 F0→F1 calls, 95.8% have the pooled sign in all three cohort coefficients of a cohort-by-stage model (95.4% in separate single-cohort fits). F4 vs F0: 4,047 of the 7,043 calls have a cohort-by-stage interaction at BH q<0.05 (2-df moderated F over 24,115 genes), and 75.4% have the pooled sign in all three cohorts (review job 1032). Tables: `figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/repair/fibrosis_counts_cohort_omission.tsv` and `figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/fig3e_cohort_heterogeneity.tsv`.
F0 on the 774 units is 103 units: 14 source controls and 89 F0 disease biopsies
(GSE130970 6 and 17, GSE135251 8 and 37, GSE162694 0 and 35;
`f0_composition_audit.tsv`). The three F0 arms
(`figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_summary.tsv`):

| contrast | all F0 (n=103) | disease-only F0 (n=89) | all F0 + control-status term (n=103) |
|---|---:|---:|---:|
| **F0→F1** | **1,338** | **783** | **814** |
| F1→F2 | 398 | 398 | 386 |
| F2→F3 | 512 | 512 | 493 |
| F3→F4 | 281 | 281 | n/a |
| F2 vs F0 | 4,266 | 3,946 | 3,679 |
| F3 vs F0 | 6,190 | 6,459 | 6,166 |
| F4 vs F0 | 7,043 | 6,929 | 6,640 |

All-gene F0→F1 effects correlate at Spearman rho 0.950 between the all-F0 and
disease-only arms (`f0_arm_effect_agreement.tsv`). The all-F0 F0→F1 count rests mainly on GSE130970: 14 of the 1,338 calls remain without it (155 units). GNMT F0→F1: q = 0.00913 (all
F0), 0.0426 (disease-only) and 0.0430 (control-status term); MAT1A: 2.45e-4,
9.5e-4 and 1.6e-3 (`f0_arm_hero_genes.tsv`). Neither F0→F1 difference holds without GSE130970 (GNMT q 0.207, MAT1A 0.0892); without GSE135251 they are q 0.0429 and 0.00251, and without GSE162694 0.0338 and 0.00631. CYP2C19 F1→F2 (q 0.0477) is q 0.0239 / 0.219 / 0.499 without GSE130970 / GSE135251 / GSE162694. The F2, F3 and F4 versus F0 differences of all three genes are negative in every fit: q < 0.001 under every F0 reference and without GSE130970, and q ≤ 0.034 without GSE135251 or GSE162694 except MAT1A F4 vs F0 without GSE135251 (q 0.0592; 11 F4 units) (`figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/repair/hero_genes_cohort_omission.tsv`). In Figure 3E and the Figure S3
stage figure (`figures/candidates/fig3e-s3-774-20261011T165944Z/`) P6 is not
shown because that axis is unavailable on the corrected counts, and the NMF
activity lines on the NAS axis are relative to NAS 1–2. A refit with NAS0 as
the reference was run on 2026-10-12 with 46 NAS0 units
(`figures/candidates/fig3e-s3-774-nas0ref-20261012T003338Z/source_tables/fig3e_prespecified_nmf_activity.tsv`).

**GSE213621 grouped-stage check (H51, 2026-10-15;
`Analysis/Hypotheses/H51-20261015T112617Z/RESULT_H51.md`, stage 2
`stage2-20261015T164829Z/result.json`).** Not blind: review job 1032 had already
fitted GSE213621 F3–F4 versus F0–F1 in 164 disease units, and its agreement with the
discovery F3-versus-F0 contrast (83.1% same sign among the 6,190 discovery calls;
all-gene Spearman 0.468; re-derived from `jobs/1032/stage_check_genes.tsv.gz`) was read
before the plan was written
(`DISCLOSURE_ADDENDUM_H51.md`). The result supports a statement about the slope and its
interval; it is not a prediction confirmed on unseen data. Units: 371 discovery disease
units refitted on bins (GSE130970 65, GSE135251 197, GSE162694 109;
`~ dataset + inferred_sex + bin`) and 258 GSE213621 disease units (F0–F1 87, F2 94,
F3–F4 77; `~ inferred_sex + bin`); the 68 GSE213621 controls are excluded. Units are
genetic-identity clusters. b is the OLS slope of GSE213621 log2FC on discovery log2FC over
discovery BH calls kept in GSE213621. Primary family: two one-sided stratified permutation
tests (1,000 bin permutations within 8 accession-quartile by sex strata), Holm,
family-wise alpha 0.05; intervals from 200 bin-by-sex stratified bootstrap draws (normal
interval primary); seed 20260923.

| Contrast | Genes | b (95% normal interval) | Percentile interval | Null mean | One-sided permutation P | Holm P | Sign concordance (null mean) | Minimum detectable slope | State |
|---|---:|---|---|---:|---:|---:|---|---:|---|
| F3–F4 vs F0–F1 | 8,449 | 0.708 (0.550 to 0.867) | 0.553 to 0.870 | 0.006 | 0.001 (floor) | 0.002 | 0.778 (0.507) | 0.243 | pass |
| F2 vs F0–F1 | 2,665 | 0.615 (0.402 to 0.828), descriptive | 0.416 to 0.852 | 0.017 | 0.001 (not interpreted) | not applicable | 0.946 (0.529) | 0.34 | untestable (S3: above 0.30) |

Reading (prespecified): grouped fibrosis-stage effects estimated in three cohorts (371
disease units) point the same way in the 258 GSE213621 disease units: replication slope
b 0.708 (95% interval 0.550 to 0.867) for F3–F4 versus F0–F1, beyond the slope that
submission order and sex structure produce. Size against held-out discovery cohorts is not
determined (ratio of mean GSE213621 slope to mean held-out-cohort slope, both on effects
from the other two discovery cohorts: 0.793, 0.488 to 1.098; F2, descriptive, 0.803,
0.440 to 1.167). F2 is untestable
under S3, and nothing about F2 transport is claimed. Positive controls: 5 of 5 kept, all 5
with the discovery sign at nominal P < 0.05 (COL1A1, EPCAM, GNMT, LUM, CYP2C19).
Negative control NC1 failed: the slope of the GSE213621 male-versus-female coefficient on
discovery F3–F4 effects is -0.155 (95% interval -0.279 to -0.031; 8,154 autosomal
genes). Two readings remain open: sex and bin are not separated in GSE213621, or sex and
fibrosis share real biology (for example sex-biased hepatocyte metabolic genes). Women
only (161 units; 49 / 55 / 57 by bin): b 0.742 for F3–F4 and 0.732 for F2.
Sensitivities, F3–F4 / F2: accession-quartile covariate 0.719 / 0.600 (not below half the
primary, so not attributed to submission order); RLE recomputed after filtering 0.709 /
0.616; up calls 0.909 / 0.739; down calls 0.577 / 0.413; all kept genes 0.674 / 0.557
(Spearman 0.467 / 0.531); discovery with a top-50-share covariate 0.680 / 0.582;
discovery without GSE162694 0.741 / 0.600. Composition-adjusted (six ILR coordinates of
seven lineage groups from `corrected_774_20261009T141331Z`): F3–F4 0.431, not below half
of 0.708; F2 0.259, below half of 0.615 (descriptive only, because F2 is untestable; no
composition reading is drawn for F2). R2 of bin on the six ILR coordinates in GSE213621 is 0.205. Both
adjusted results depend on the deconvolution input, and no gene is named from adjusted
fits. Over the printed Figure 3E calls, the GSE213621 F3–F4 slope is 0.638 (0.495 to
0.781) for F3 versus F0 and 0.378 (0.269 to 0.487) for F4 versus F0; the Figure 3E F0
includes controls. F3–F4 per-gene states over the 8,560 discovery calls (BH over the 8,449
kept in GSE213621): replicated 4,580, opposite sign 353, not replicated at adequate power
807, underpowered 2,709, not measurable 111; 525 change state with the accession
covariate. The per-gene states are a Catalog r9d field (adopted 2026-10-16). H51 is not
independent support for the 11,130 stage-supported genes, which already include
GSE213621. Figures 3E and 3F remain three-cohort results.

**Library-property check of the stage contrasts (H82, 2026-10-15;
`Analysis/Hypotheses/H82-20261015T125052Z/run-20261015T131530Z/RESULT_H82.md`).**
Untestable under rule S5: the planted-covariate positive control fell below the null
2.5th percentile in 18 of 20 draws, where all 20 were required. By the prespecified
reading no H82 number enters the paper. On 2026-10-21 the user approved a logged post hoc
deviation (decision D4; `H82-20261015T125052Z/DEVIATIONS.md`, "Post hoc reporting
deviation"): the D1 joint-adjustment counts for all 11 primary-test contrasts are reported
descriptively in Supplementary Note 7, with no F1-specific wording and no T1/T2 result
(table: `REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/h82/h82_joint_adjustment_retention_11_contrasts.tsv`).
Figure 3E counts are unchanged.

The 2026-09-08 rigor review found that the older Figure 3 stage-extension table
omitted confidence intervals in all 257,070 gene–contrast rows. A separate
[interval repair](technical/oakes_pachter_rigor_20260908/stage-interval-repair/reproduction_checks.tsv)
reproduces all 11 fits and supplies finite pointwise 95% intervals without
changing effects, SE, P values or FDR beyond numerical precision. The original
table and figures were not overwritten, and the repair is not adopted. This
does not repair the F0 reference mixture or restore the early-stage claims
qualified below.

⛔ **The F0 reference mixes source controls and F0 disease biopsies.** With
source diagnosis labels (review item 1; GSE130970 controls are the 6 biopsies
reconstructed from Hoang et al. 2019 Supplementary Table 1), F0 is 104
participants: **14 source controls and 90 F0 disease biopsies** (GSE130970 6 and
17, GSE135251 8 and 38, GSE162694 0 and 35; GSE162694's controls carry
`fibrosis_stage = NA` and are dropped upstream). F1–F4 are effectively clean: one
control in F1, one in F2, none above.

Rerun three ways on the 844-record corrected-controls fit (historical; the
774-unit arms are tabulated above) in
`figures/candidates/f0-reference-arms-corrected-controls-20260923T231939Z/analysis/f0_arms/`
(`f0_arm_summary.tsv`, `f0_composition_audit.tsv`; same model and BH family):

| contrast | all F0 (n=104) | disease-only F0 (n=90) | all F0 + control-status term (n=104) |
|---|---:|---:|---:|
| **F0→F1** | **1,426** | **929** | **975** |
| F1→F2 | 544 | 544 | 524 |
| F2→F3 | 518 | 518 | 507 |
| F3→F4 | 234 | 234 | n/a |
| F2 vs F0 | 3,904 | 3,611 | 3,402 |
| F3 vs F0 | 6,333 | 6,248 | 6,042 |
| F4 vs F0 | 6,462 | 6,325 | 6,051 |

**The 2026-09-07 control-confounding explanation is withdrawn.** The earlier arms
(`figures/candidates/f0-reference-arms-20260907T105813Z/`, job `21489343`) gave
1,471 / 7 / 258 F0→F1 calls and were read as showing that the F0→F1 signal was
mostly disease versus control. That reading depended on labelling all 23 GSE130970
F0 biopsies as controls, 16 of which have steatosis above zero. With source labels,
929 of the 1,426 all-F0 calls (about two thirds) remain with a disease-only F0
reference, and all-gene F0→F1 effects correlate at Spearman rho 0.950 between the
all-F0 and disease-only arms (`f0_arm_effect_agreement.tsv`). The all-F0 arm
remains the cross-sectional description; quote the disease-only value beside it,
because the F0 reference still contains 14 source controls.

The worked-example F0→F1 effects now remain supported in all three arms. GNMT:
q = 0.0089 (all F0), 0.0365 (disease-only) and 0.0383 (control-status term); MAT1A:
q = 1.4e-4, 6.5e-4 and 9.0e-4 (`f0_arm_hero_genes.tsv`). The earlier statement
that these two claims do not survive is withdrawn. The pooled 844-participant
disease-versus-control result for GNMT, MAT1A and CYP2C19 is a different
contrast.

Also fixed there: `scripts/figures/build_pi_stage_extensions.R:112` builds its
t critical value with an unnamed `qt()` vector and then indexes it by gene name,
so all 257,070 confidence intervals in the adopted stage-extension output are
NA. The rerun names the vector and fails closed if any interval is missing.

| Axis | Contrast | Reference / comparison participants | Cohorts | Higher / lower at BH q<0.05 |
|---|---|---:|---:|---:|
| Adjacent fibrosis | F0→F1 | 104 / 105 | 3 | 658 / 768 |
| Adjacent fibrosis | F1→F2 | 105 / 90 | 3 | 350 / 194 |
| Adjacent fibrosis | F2→F3 | 90 / 75 | 3 | 376 / 142 |
| Adjacent fibrosis | F3→F4 | 75 / 25 | 3 | 169 / 65 |
| Fibrosis versus F0 | F1 vs F0 | 104 / 105 | 3 | 658 / 768 |
| Fibrosis versus F0 | F2 vs F0 | 104 / 90 | 3 | 2,012 / 1,892 |
| Fibrosis versus F0 | F3 vs F0 | 104 / 75 | 3 | 2,996 / 3,337 |
| Fibrosis versus F0 | F4 vs F0 | 104 / 25 | 3 | 3,020 / 3,442 |
| Adjacent NAS | 0→1–2 | 46 / 65 | 3 | 2,684 / 2,437 |
| Adjacent NAS | 1–2→3–4 | 65 / 121 | 3 | 2 / 1 |
| Adjacent NAS | 3–4→5–8 | 121 / 175 | 3 | 1,860 / 2,197 |

Counts in the table directly above are from the 844-record corrected stage
release (historical)
(`analysis/stage_extensions/stage_all_gene_results.tsv` and
`stage_extension_all_gene_results.tsv`, BH q<0.05 split by log2FC sign). The
2026-08-10 values were 682 / 789, 354 / 194, 376 / 146, 174 / 70, 2,011 / 1,893,
2,986 / 3,334, 2,999 / 3,425, 2,692 / 2,435, 2 / 1 and 1,877 / 2,213.

Higher and lower are defined in the comparison group relative to the named
reference. These counts use BH q<0.05 without an added fold-change threshold.
The continuous NMF display starts at NAS1–2 because only one NAS0 donor has a
prespecified NMF score; the NAS0-versus-NAS1–2 DEG model still uses all 111
eligible participants.

The prior 661-sample, six-cohort stage candidate and every earlier read-count
stage table are provenance only. They cannot supply Figure 3 counts or be
combined with the selected five-cohort pooled result. Every stage contrast is
cross-sectional and uses participants from cohorts carrying the compared
groups; no comparison represents within-person change.

## Histology-anchored molecular continuum candidate

The validated candidate
`RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-20260818T173348Z`
adds an externally anchored molecular continuum without changing Figure 3 or
replacing recorded histology. GSE162694 and GSE213621 are modeled separately as
competitor-independent evaluation cohorts. Within each cohort, the primary
model is `outcome_z ~ continuum_z + factor(fibrosis_stage) + inferred_sex`;
fixed-effect meta-analysis then combines cohort coefficients. The cohort-PC1
and fixed external projection are co-primary, and the consensus rank is
visualization-only. All 139 observed genes from the published 145-gene
signature are removed before pathway,
Hotspot, and fixed-basis NMF scoring; ordinary transcript inference excludes
the 139 observed axis components, which remain explicit and are evaluated only
with target-specific leave-one-gene-out axes.

| Layer | Complete family / testable | Continuum-associated | Matched fibrosis-trend comparison |
|---|---:|---:|---:|
| Ordinary transcripts | 23,370 / 23,231 | 14,651 | 10,476 stage-supported; 6,836 shared; effect rho=0.314 |
| Hotspot programs | 117 / 116 | 93 | 64 stage-supported; 54 shared; effect rho=0.803 |
| Hallmark | 50 / 44 | 37 | 29 stage-supported; 27 shared; effect rho=0.693 |
| KEGG | 658 / 260 | 170 | 123 stage-supported; 100 shared; effect rho=0.622 |
| Reactome | 1,787 / 1,027 | 781 | 607 stage-supported; 525 shared; effect rho=0.643 |
| GO biological process | 7,583 / 2,811 | 2,108 | 1,786 stage-supported; 1,531 shared; effect rho=0.696 |
| GO molecular function | 1,855 / 720 | 462 | 342 stage-supported; 259 shared; effect rho=0.642 |
| GO cellular component | 1,042 / 467 | 332 | 244 stage-supported; 203 shared; effect rho=0.648 |
| Fixed-basis NMF axes | 10 / 10 | 8 | 7 stage-supported; 6 shared; effect rho=0.733 |

**774-unit rerun (2026-10-09 transcript layer; 2026-10-11 remaining layers;
panels installed 2026-10-14).** The table above is on the 23,370-gene,
844-record basis and is historical. Transcript layer
(`RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/bulk/`):
the complete family is 24,115 genes, 23,972 testable after removing the 143
observed genes of the 145-gene signature, and 15,215 are continuum-associated
(`run_summary.tsv`). Of the 1,854 canonical DEGs, 1,312 are
continuum-associated, with direction agreement 0.853 on the signature axis and
0.860 on the fixed projection (`canonical_deg_concordance.tsv`). The held-out
check (308 discovery and 466 evaluation units; the evaluation units are 140
from GSE162694 and 326 from GSE213621) gives 2,193 discovery genes, direction
agreement 0.923, 1,574 BH-supported discovery genes and all-gene Spearman 0.588
(`heldout_concordance_summary.tsv`). Stage-versus-continuum Spearman on the
signature axis is 0.348, 0.140, 0.371 and 0.403 for F0→F1 to F3→F4
(`stage_continuum_concordance.tsv`). The stage-adjusted models use 435 units
(109 GSE162694 and 326 GSE213621, the latter including 68 controls that the
stage code places at stage 0; the last column of the table below refits both
axes without them;
`logs/continuum_774/run_20261009T171017Z/comparison_844_vs_774.tsv`). The
matched fibrosis-trend comparison on the 774 units
(`figures/candidates/continuum-774-panels-20261011T165944Z/stage_continuum_comparison/comparison_summary.tsv`;
testable counts from `decision/layer_summary.tsv` under
`hac-molecular-layers-corrected774-20261011T152106Z`):

| Layer | Complete family / testable | Continuum-associated | Matched fibrosis-trend comparison | Effect rho, 367 disease units |
|---|---:|---:|---:|---:|
| Ordinary transcripts | 24,115 / 23,972 | 15,215 | 11,130 stage-supported; 7,334 shared; effect rho=0.356 | not recomputed |
| Hotspot programs | 117 / 116 | 93 | 64 stage-supported; 54 shared; effect rho=0.823 | 0.629 (81/116 same direction) |
| Hallmark | 50 / 45 | 38 | 30 stage-supported; 29 shared; effect rho=0.745 | 0.662 |
| KEGG | 658 / 302 | 206 | 145 stage-supported; 126 shared; effect rho=0.669 | 0.675 |
| Reactome | 1,787 / 1,124 | 873 | 675 stage-supported; 593 shared; effect rho=0.687 | 0.612 |
| GO biological process | 7,583 / 3,317 | 2,567 | 2,163 stage-supported; 1,933 shared; effect rho=0.740 | 0.671 |
| GO molecular function | 1,855 / 824 | 539 | 409 stage-supported; 323 shared; effect rho=0.683 | 0.523 |
| GO cellular component | 1,042 / 516 | 385 | 277 stage-supported; 242 shared; effect rho=0.700 | 0.547 |
| Fixed-basis NMF axes | 10 / 8 with a stage estimate | 8 | 7 stage-supported; 6 shared; effect rho=0.476 | not recomputed |

The first four columns use the 435 stage-complete units, including the 68
GSE213621 controls placed at stage 0. The last column refits both axes in the
367 disease units (109 GSE162694, 258 GSE213621; review round 3, R05-02,
`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/feature_layer_stage_vs_continuum_rho_by_roster.tsv`).
The stage-supported, shared and continuum-associated counts were not
recomputed on that roster.

Two NMF axes (k4-P2 and k6-P6) have no stage estimate on the corrected counts.
The paired step uses 54 GSE193066 participants whose identity calls are
unresolved (`paired/paired_layer_summary.tsv`).

Historical 844-record corrected-controls rerun of the ordinary-transcript layer
(`RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected-controls-20260925T192142Z/bulk/`):
the complete family is 23,702 genes, 23,563 testable after removing the 139 axis
components, and 14,818 are continuum-associated (`run_summary.tsv`,
`continuum_membership.tsv.gz`). Of the 1,725 canonical DEGs, 1,173 are
continuum-associated, with direction agreement 0.808 on the signature axis
(previously 923 of 1,347 and 0.832). The held-out check of a disease signature
discovered in the source-overlap cohorts (343 discovery, 501 evaluation
participants) gives 3,669 discovery genes, direction agreement 0.739 in the
evaluation cohorts, 1,897 BH-supported discovery genes and all-gene Spearman
0.323 (previously 1,146 genes, 0.873, 767 and 0.426;
`heldout_concordance_summary.tsv`). Stage-versus-continuum Spearman on the
signature axis is 0.348, 0.105, 0.379 and 0.366 for F0→F1 to F3→F4. The
matched fibrosis-trend comparison in the transcript row and the other layers in
the table were not rerun.

The two transcript-derived continuum effect vectors agree at rho=0.9991. Both
co-primary scores support positive stage-adjusted associations for the
prespecified Ductular injury and Stromal ECM programs in both evaluation
cohorts, and the additive Figure 4F adoption criteria pass. The complete
stage-versus-continuum review candidate is
`figures/candidates/histology-vs-continuum-comparison-20260818T192630Z`.
All five Resource cohorts were scored. The current descriptive continuum-window
figures show GSE126848 (n=53), GSE130970 (n=76), GSE135251 (n=214), GSE162694
(n=140), and GSE213621 (n=361). Stage displays contain the latter four
stage-available cohorts; GSE126848 has no recorded fibrosis in the selected
manifest. The first three cohorts remain source-overlap descriptive context and
do not enter the primary two-cohort meta-analysis.
Program and pathway rankings align strongly with fibrosis-based analyses, while
the weaker transcript-level concordance and significant stage-adjusted effects
show that the continuum is not merely a relabeling of stage groups.

All sliding windows are descriptive, fixed in advance, and excluded from
inference. GSE193066 paired models are descriptive: the 54 source-declared pairs have
unqualified identities (K01), so they neither upgrade nor veto cross-sectional
findings. The Ductular injury paired slope on the signature axis is 0.634
(SE 0.096; BH q 6.3e-7 across 117 programs), a paired association in
unqualified pairs, not within-person support. Stromal ECM is untestable in the
paired substrate because retained signature-excluded weight is 0.738, below 80%. The safe
claim is that the continuum complements histology by describing
cross-sectional molecular variation within recorded stage. It does not support
a biomarker, subtype, causal, therapeutic-window, or within-person
temporal-ordering claim. The synchronized Resource candidate remains unpromoted.
On 2026-08-24, hash-pinned Figure 4F–4G and selected S3/S4 panels were adopted to
the working figure authorities, then revised from two-cohort-only descriptive
views to the five-cohort context package
`histology-continuum-five-cohort-context-20260824T175647Z`. This working-figure
revision does not adopt the synchronized Resource release.

The molecular-systems layer exists to give the continuum a unit of
interpretation between single genes and whole-cohort axes. It compresses 13,092
correlated features into a few dozen auditable groupings so that within-stage
remodeling can be described as coordinated rather than gene by gene. That
purpose is organizational rather than inferential, and the stress tests below
bound it: the groupings are a navigable vocabulary and evidence that within-stage
variation is shared across pathway families and single-cell-derived programs,
not a ranking of systems by disease relevance.

The corrected outcome-blind molecular-systems atlas is
`figures/candidates/histology-continuum-molecular-systems-atlas-20260824T195148Z`.
After exclusion of the 139 observed genes from the published 145-gene
signature and GENCODE-v49 mapping, it retains 5,329 testable pathway sets and
116 testable Hotspot programs from the complete 13,092-feature ledger. Of
5,445 mapped nodes, 5,439 form one main graph component and 43 size-ranked
systems. Six nodes form three disconnected two-node microcomponents; they
remain explicit in the registry and topology view but are not named or included
in system summaries. Stage and stage-adjusted fixed-projection effects were
added only after the coordinates, component registry, and community IDs were
checksummed.

Across the 5,439 main-component nodes, stage-versus-continuum effect rho is
0.683 and 80.5% share direction. Equal-collection system summaries have
rho=0.658 and 36/43 share direction (844 records; 0.7265 and 37/43 on the 435 stage-complete
units of the 774-unit substrate, with the 68 GSE213621 controls coded as stage 0; 0.631 and 32/43 on the 367 disease units; Main Figure 4H plots the donor-first estimates, 0.670 and 34/43); feature-balanced medians have rho=0.701
and 37/43 share direction. Fixed-window direction agrees between GSE162694 and
GSE213621 for 39/43 systems under either aggregation. Across all five
descriptive cohorts, the counts are 13/43 and 12/43, respectively. These values
are descriptive summaries over nested gene sets, correlated cohorts, and
overlapping windows; no independent-sign null or inferential p value is used.
Outcome-blind topology sensitivity gives adjusted Rand indices of 0.633 at
k=20, 0.751 at k=40, 0.855 at resolution 0.7, and 0.815 at resolution 1.3
relative to k=30/resolution 1. Biological terms remain in the top-ten
representative audit rather than becoming single-medoid system names. Figure
4G, Figure 4H, and S4Q–S4U are adopted working panels; Figure 3 and the synchronized
Resource release remain unchanged.

Figure S4U pools the five descriptive system-window profiles only after
within-cohort standardization and percentile ranking. For each system and
window, collection medians are averaged within cohort and the five resulting
cohort scores receive equal one-fifth weight. The source lineage retains the
complete 43-by-5-by-9 cohort input and 43-by-9 pooled output. This is a
descriptive aligned-cohort summary, not participant-level pooling or an
additional significance test.

The separate donor-level system analysis constructs each system score within
participant, then fits `system_score_z ~ continuum_z + factor(fibrosis_stage) +
inferred_sex` in GSE162694 and GSE213621. Every stage-adjusted model is
complete-case in recorded fibrosis stage. Thirty-one of the 140 GSE162694
participants have no recorded stage in the selected manifest, so the
stage-adjusted models use 109 GSE162694 and 361 GSE213621 participants, 470 in
total, not the 501-participant descriptive roster. Fixed-effect
meta-analysis uses BH over the complete 43-system family for each co-primary
axis. Thirty-three systems meet both meta-analysis q<0.05 criteria with the same
direction in all four cohort-by-axis coefficients; 29 retain support with
feature-balanced scoring. Twenty-seven systems are fibrosis-stage-associated,
24 pass both stage and continuum criteria, and nine are continuum-only.
Within-stage/sex permutations give BH support in all four tests for 22 systems;
18 pass the stricter maxT family-wise P<0.05 criterion in all four tests.
The corrected 774-unit system stress source is
`RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-corrected774-20261009T171017Z/tables/system_stress_summary.tsv`.
It gives 33 systems passing joint maxT on both axes, zero passing the
connected-subgraph BH null on both axes and one passing the uniform
size-matched BH null on both axes. Stage-adjusted n=435 (109 +326, including the 68 GSE213621 controls, which
carry no recorded histology); these are
system counts requiring support on both axes, not counts of significant
axis-specific rows.

On the 774 units (`logs/continuum_774/run_20261009T171017Z/comparison_844_vs_774.tsv`; panels in `figures/candidates/continuum-774-panels-20261011T170854Z/` and `continuum-774-panels-20261011T165944Z/`, Figure 4G/H installed on 2026-10-14 and replaced on 2026-10-16 by the renders on the 367 disease units) the same 18 systems pass strict maxT; joint BH support goes from 36 to 37 systems, joint maxT from 35 to 33 and continuum-only maxT from 9 to 7; 0 of 43 systems exceed the connected-subgraph null and 1 of 43 (system_01) exceeds the size-matched null. The paragraph that follows describes the 844-record basis. The 43-system membership geometry (memberships, weights and coordinates) uses the adopted F_five 23,370-gene, 844-record expression-filtered universe, and the system scores in Figure 4G and 4H use the 24,115-gene expression matrix of the 774-unit fit; the corrected disease-versus-control fit keeps 23,702 genes (332 added, none lost), but the GSE213621 stage code places its 68 controls at stage 0, so the earlier Figure 4G and 4H models included those controls as their own stage category and the 4H stage axis read the control label; the installed Figure 4G and 4H (2026-10-16) are the refits on the 367 disease units, 109 GSE162694 and 258 GSE213621, without the 68 controls (review round 3, R05-02; `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/repair/fig4g/` and `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`); the other values stand on the corrected inputs (outcomes
document, "Continuum bulk transcript layer").

Two boundaries on those counts are recorded here rather than left to the reader.
First, the two evaluation cohorts do not carry the same stage scale: GSE213621
records levels 0 through 3 and GSE162694 records levels 0 through 4, yet the
secondary `stage_z` model standardizes each cohort's own scale before
fixed-effect meta-analysis, so the combined stage coefficient is not expressed in
a single common histologic unit. Second, none of the nine continuum-only systems
pass the strict criterion; all 18 maxT-supported systems are
stage-and-continuum-supported. The strict criterion is the conjunction of four
maxT-corrected axis-by-cohort tests, but the two co-primary axes give system
coefficients correlated at r=0.9993 (0.9997 and 0.9989 within cohort), and
requiring both axes removes none of the 18 (one axis in both cohorts gives 18
on either axis), so it is in effect two cohort tests. GSE162694 contributes 109
participants across ten stage-by-sex strata, so the conjunction is limited by the
smaller cohort's within-stratum permutation resolution. All nine continuum-only
systems agree in direction across both cohorts; their GSE162694 residual
correlations are smaller than their GSE213621 values. The strict arm therefore
does not by itself establish continuum association independent of recorded stage
at the system level. The
donor-first and plotted pooled-window estimands agreed at Spearman rho=0.997
on the 844 records. No like-for-like 774-unit reconciliation was run: the
774-unit donor-first windows against the fixed 844-record feature-first
windows give Spearman 0.849 (Pearson 0.922) for the equal-weight five-cohort
pool.
Main Figure 4H shows the complete 43-system fibrosis-stage versus
stage-adjusted continuum comparison. Labels remain restricted to the 15 largest
systems fixed before outcome overlay. It uses neither sliding-window values nor
membership geometry. The rho=0.658 and 36/43 same-direction count are the
844-record values. On the 435 stage-complete units of the 774-unit substrate the panel source gives Spearman rho 0.7265
(printed as 0.73 on the earlier panel) and 37/43 in the same direction; with the 68 GSE213621 controls removed from both axes (367 disease units) the values are 0.631 and 32/43 (feature-first 4H installed 2026-10-16; `figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`). Since 2026-10-17 Figure 4H plots the donor-first estimates on the same 367 units: rho 0.670 (bootstrap 95% interval 0.382 to 0.842) and 34/43 (`figures/main/fig4_singlecell_programs/panels/data/panel_deep_20261017/4H.tsv`)
(`figures/candidates/fig3-fig4-layout-774-20261013T135430Z/panels/fig4/fig4h_single_outline/fig4h_source.tsv`).
Both remain descriptive summaries of correlated systems.

The revised panel is
`figures/candidates/histology-continuum-fig4h-competitive-mark-20260827T130920Z`.
Every point, coordinate, label, and summary value is carried over unchanged; only
what the ink encodes has moved. The outline previously marked the 18 systems
clearing the four-test maxT conjunction. That criterion is set by the smaller
cohort's within-stratum permutation resolution, and the jointly powered
cross-cohort test supports 35 of 43, so neither count is scarce enough to carry
ink. The outline now marks the single system exceeding size-matched draws on both
co-primary axes. Both counts are retained in the panel source table and manifest
for the caption authority. The rebuild also repairs a defect carried by the
adopted panel: the direction-colour legend keys rendered black because the
legend-only layer contributed its fixed colour to the colour guide.

Supplementary panels S4V and S4W are
`figures/candidates/histology-continuum-system-competitive-null-20260827T131027Z`.
S4V shows the observed system_01 statistic against its own null distributions
under both references and both axes, which makes the separation between them
visible: the connected-subgraph null is broad and reaches the observed value,
while the unconstrained null is tight. S4W retains all 43 systems and plots the
deciding quantity, the worse BH q across the two co-primary axes, under each null.
Figure S4U retains all 43 systems and marks only the 18 strict maxT-supported
rows with a neutral `*` beside unchanged window tiles. In the 367 disease units,
where the unadjusted fit calls all 18 (BH q<0.05 within 43 systems), 4 of the 18
are called in all 11 estimable composition-adjusted fits, 14 lose the call in at
least one fit, 9 are called in 5 or fewer and 1 in none
(`figures/candidates/continuum-rare-type-arms-20261012T020132Z/arms/meta_effects_all_arms.tsv`). The 33 primary BH calls
remain in the complete source registry; the panel does not assign P values to
windows or `Delta`. The active Figure 4G/S4Q–S4U display revision removes
embedded explanatory footnotes without changing values, geometry, order, or
inference; the full explanations remain in the external caption authority.
Figure 4H is hash-pinned in
`fig4h_system_stage_continuum_promotion_20260825.tsv`.

The active visual revision orders both system heatmaps by the equal-cohort
pooled late-minus-early contrast, defined as the mean of windows centered at 80
and 90 minus the mean of windows centered at 10 and 20. The active S4U omits
this redundant contrast tile but retains the same row order. This is a
result-adaptive presentation order only; all underlying values, geometry,
support calls, and testing families are unchanged.

Retrospective transport checks (2026-08-28;
`RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z/PRIMARY_transfer_verdict.tsv`):
the 139-gene continuum loading built on the 23,370-gene, 844-record basis,
projected without refitting, correlates with fibrosis stage at rho 0.510 in 93
GSE174478 samples and 0.710 in 66 GSE240729 samples (FFPE). Against 1,000 matched
random-signature draws, empirical P is 0.117 and 0.032, Holm-corrected over the two
cohorts to 0.117 and 0.064. Units are samples; donor counts are not recorded. This
check was not rerun on the 774-unit basis.

### Molecular-system stress tests

The stress-test candidate is
`RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-20260827T125516Z`.
It reads only the checksummed atlas and inference candidates, verifies both
against their checksum registries, and reproduces the released collection-balanced
participant scores to within 5.3e-15 before any test runs. It changes no
geometry, community ID, Hotspot membership, or adopted panel.

A jointly powered cross-cohort permutation replaces the conjunction of four
separately corrected axis-by-cohort tests. The axis is permuted within
cohort-by-stage-by-sex strata, the two cohorts are combined by weighted Fisher z
of the residual partial correlation, and maxT runs over the complete 43-system
family. Thirty-five of 43 systems pass family-wise P<0.05 on both co-primary
axes, against 18 under the four-test conjunction, and all nine continuum-only
systems are among them. The earlier absence of continuum-only systems from the
strict arm was therefore a resolution limit of the smaller cohort's
within-stratum permutation, not evidence against within-stage association. Both
criteria are retained; the prespecified 18-system count is unchanged.

The competitive nulls ask a different question: whether a system is more
continuum-associated than an equally sized region of the same geometry rather
than more than zero. Each system is compared against 5,000 size-matched draws,
scored and fitted through the identical pipeline, under a connected snowball
sample of the prespecified mutual-30NN graph and under an unconstrained node sample.
The per-draw null statistics are retained so the supplementary panels can show
the distributions themselves. Under the conservative connected null, **no system
passes BH q<0.05 on both axes**; under the permissive unconstrained null, one
does. That system is
system_01, the largest system at 311 features, whose representative terms are
integrin binding, ECM organization, contractile myofibroblast, and cell
adhesion. Its connected-null P values are 0.0062 and 0.0078, which do not
survive BH over 43 systems, and the next-ranked system is at P=0.115. Drawn sets
are not hub-enriched relative to the tested systems: mean node degree is 22.1
under the connected null and 21.3 under the unconstrained null against 19.7 for
system_01. The 33-of-43 and 35-of-43 counts therefore describe a diffuse axis
shared across the membership map, consistent with 14,651 of 23,231 transcripts
being continuum-associated after stage adjustment (14,818 of 23,563 on the
844-record corrected-controls fit; 15,215 of 23,972 on the 774 units). Except for the ECM system
against the permissive reference, no system is shown to carry
continuum association specific to its own membership.

The 43 systems are also not 43 independent quantities. Across the 466 pooled
evaluation units of the 774-unit substrate their score correlation matrix gives a Li and Ji
effective number of 14 independent systems, with 41.1% of variance in the first
component and 74.7% in the first three; the per-cohort effective numbers are 12
and 14 (`RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-corrected774-20261009T171017Z/tables/system_redundancy.tsv`). Superseded 2026-10-09: on the 844-record basis (501 units) the values were 15 systems, 36.5% and 72.4%. BH over 43 controls multiplicity but not this redundancy.

**Historical repeat-biopsy diagnostic, qualified 2026-10-14.** The following
GSE193066 calculations retain source-declared pair labels and historical
numbers. All 164 cohort libraries remain identity-unqualified. Its 58-pair ICC,
SD 0.6433 and derived reliability/ceiling calculations are diagnostic scenarios,
not primary within-person evidence or empirically qualified label-error
calibration. They do not resolve the continuum-versus-stage comparison.

Finally, the apparent strength of the continuum relative to fibrosis stage is
largely a comparison artifact. The 2.79-fold median ratio visible in the Figure
4H source table compares a marginal stage coefficient against a
stage-adjusted continuum coefficient. Fitting both marginally in the same
participants gives a median ratio of 1.57, and binning the continuum to each
cohort's own stage coarseness gives 1.45. The reliability correction is no longer assumed. Kleiner
fibrosis staging was compared across source-declared biopsies in GSE193066, in 58 pairs
with repeat biopsies a median two years apart: exact agreement 28 of 58, within
one stage 52 of 58, and a within-participant standard deviation of 0.6433
stages with a participant-bootstrap 95% interval of 0.5170 to 0.7600. The
retained diagnostic quantity is that standard deviation rather than a primary
reliability coefficient; the same assumed variance yields reliabilities of 0.19, 0.28
or 0.47 in this deposit depending only on which stage variance is used as the
denominator. Applied as rho = 1 - 0.414 / Var(stage) in the regression-slope
convention used above, the reliability-corrected ratio spans 0.27 at the paired
subset's own restricted variance to 1.03 at a stage standard deviation of 1.2,
so a coarse, imperfectly reproducible ordinal and a continuous expression
projection remain unseparated by these data. Because identities are unqualified,
this calculation remains a sensitivity scenario rather than measured error calibration.

That measurement does not transport unchanged. The same statistic was computed
in an independent deposit, PXD051911, in fifty-seven patients with repeat
Kleiner grades: 0.5092 stages with an interval of 0.3043 to 0.6939 in the
twenty-seven initial-to-surgery-visit pairs (age differing by 0 years in 12
pairs and by 1 year in 15), and 0.4655 stages with an interval of 0.3651 to
0.5477 in the thirty follow-up pairs with no surgery-visit grade (age differing
by 2 to 4 years; 2 of the 30 lost more than 10 BMI points, and the other 28 give
0.4818, interval not computed), the 0.4655 value
formally discrepant with GSE193066 because the ratio interval of 0.5366 to
0.9555 excludes one. Neither stage composition nor real interval change accounts
for the gap: standardising to the source stage mix moves the value further away
rather than closer, and the follow-up arm shows a significant net upward
shift in stage, which can only inflate a within-pair variance. Two deposits cannot support a
pooled estimate. The historical SDs 0.4655 and 0.6433 define scenario inputs,
not an empirically qualified transportable error range. The direction of the
resulting bias is counterintuitive and matters: a larger standard deviation
implies a lower reliability, hence a lower ceiling and a larger disattenuated
correlation, so the GSE193066 value is the flattering end for any claim that a
model approaches the measurement ceiling, and 0.4655 is the conservative
choice.

This is a historical statistic across source-declared serial biopsies. It
folds reader variability, biopsy sampling variability and any real change into
one term, and no inter-reader substrate exists in any deposit in this project,
so no claim about agreement between readers is available. The disagreement is
symmetric, 15 up and 15 down with a sign test p of 1.0 and a mean delta of
-0.034 stages whose interval covers zero, which excludes a net upward drift in
stage as the explanation but does not exclude balanced bidirectional change;
the arm that could have caught the latter is null and underpowered. The intraclass
correlation in this sample is 0.1886 with a bootstrap 95% interval of -0.0918
to 0.4343 and does not clear a shuffled-pairing null, so it is reported with its
interval and never as a point reliability. The same coefficient computed on the
same label in PXD051911 is 0.62, so no statement of the form that Kleiner
staging has a reliability of 0.19 is available. The gap is not a denominator
artefact: the implied stage variances are 0.510 and 0.558, within nine per cent
of each other, while the error variances are 0.414 and 0.217, a factor of 1.9.
Because a reliability is one minus a ratio of those two quantities, a 1.9-fold
difference in error variance produces a 3.2-fold swing in the coefficient at an
essentially unchanged denominator. Across the two deposits the coefficient moves
3.2-fold while the standard deviation moves 1.38-fold, and it is that stability
under the same underlying difference, not any claim that the error itself is
constant, which is the reason the standard deviation is the quantity reported
here. Because a repeat-biopsy interval
loads real change into the within-participant term, the measured value is a
lower-bound interpretation would require qualified identities and additional
error assumptions. The 58 source-declared pairs are all
high-risk at first biopsy and their stage variance is 0.74 times that of the
full 106, so they are a selected set.

## Single-cell programs and composition

| Quantity | Value | State |
|---|---:|---|
| Prespecified Hotspot programs | 117 | Prespecified registry |
| Cell types represented | 5 | Hepatocyte, macrophage, fibroblast, cholangiocyte, T cell |
| Fitted donors with fewer than 50 cells of the tested lineage (H59 fit) | T cells 28 of 53 (17 below 20; minimum 1); cholangiocytes 6 of 38; fibroblasts 10 of 43; macrophages 5 of 25; hepatocytes 0 of 55 | Donor means enter the stage models unweighted; no T-cell program is near selection (smallest T-cell permutation p 0.1947, q 0.520) |
| Programs meeting primary association/stability criteria | 0 | H59: 60 ambient-correctable donors, per-dataset residual variance, within-dataset permutation BH over all 117 programs; smallest q=0.323375 |
| Programs meeting the 62-donor equal-variance sensitivity | 1 | BICC1 only (q=0.037633); zero without the two GSE189600 donors (q=0.098523) |
| Programs passing HC3 family correction (62-donor sensitivity) | 0 | Corrected ECM HC3 q=0.165240 and BICC1 q=0.155826 |
| Integrated libraries / biological donors / cells | 275 / 104 / 1,232,318 | Descriptive integration census; Liver Atlas contributes 42 runs, 38 assays, and 19 biological donors |
| Analyzed libraries / biological donors / cells | 273 / 102 / 1,232,285 | Two GSE189600 libraries with 13 and 20 annotated cells are excluded; donor collapse is mandatory for inference |

The corrected donor census replaces the stale 125-donor count, which treated
Liver Atlas sequencing runs as donors. The corrected disease-state refit contains at most 62
biological donors across GSE185477, GSE202379, GSE244832, GSE174748, and the
two analyzed GSE189600 donors, with no Liver Atlas rows. With dataset as a
covariate, GSE174748 (2 Healthy) and GSE185477 (3 Healthy) contain one stage
and do not inform the stage slope. In that 62-donor fit the slope is estimated
from 57 hepatocyte donors in three datasets and from 38 cholangiocyte, 43 fibroblast,
25 macrophage and 53 T-cell donors in two datasets (GSE202379, GSE244832),
with four macrophage Steatosis donors in total. The H59 primary fit also
drops GSE189600, leaving 55 hepatocyte donors. The two GSE174748
cirrhotic explants are excluded under the unchanged three-stage rule.
The descriptive 102-donor categories are 24 Healthy, 14 Steatosis,
36 Steatohepatitis, 17 Cirrhosis and 11 Unresolved. The exact Liver Atlas
[Table S6 join](../Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/source/Liver_Atlas_TableS6_donor_join.tsv)
gives four Healthy, three Steatosis, one Cirrhosis and eleven Unresolved;
all 19 remain excluded from stage inference. GSE202379 is mapped
from its source diagnoses rather than fibrosis: Healthy control to Healthy,
NAFLD to Steatosis, and NASH without cirrhosis to Steatohepatitis. All 117
program memberships and weights are unchanged and all 117 effects are
estimable. The four-selected/two-HC3 values belong to the superseded metadata
analysis and are not reader-facing results.

**The 113/114 tested and 27/26 supported pairs are not a conflict and are not
superseded; they are two arms of the bulk program map, reconciled by
recomputation on 2026-09-07.** The weighted arm requires retained source-weight
coverage of at least 0.80, tests **113** programs in a 226-test fibrosis/NAS BH
family, and supports **27** at q<0.05. The unweighted vocabulary arm requires
mapped-gene coverage of at least 0.80, which additionally admits T-cell module 9
(`hotspot_tcells_eaa9c262b58d9323`, "NK / cytotoxic (KLRD1)", weight coverage
0.783), tests **114** in its own 228-test family, and supports **26**. The 26
are a strict subset of the 27. The single program that differs is
**cholangiocyte module 18, "CFTR / ion transport"**, at weighted q = 0.0369 and
unweighted q = 0.0541 — one borderline call either side of 0.05, moved by the
scoring and the family size, not by a disagreement about the biology. The
program that enlarges the T3 family is supported in neither arm. Quote 27/113
for the weighted program map and 26/114 for the vocabulary comparison, and never
mix them in one sentence. The 27/113 and 26/114 belong to the 469-record working
roster. On the 364-unit three-cohort Figure 4E arm the unweighted vocabulary
count is 1/114 against a matched-random mean of 1.24 (p = 0.473), and NAS 0
against 0.19 (p = 1); on the 457-record four-cohort sensitivity arm it is
20/114 against 26.91 (p = 0.871), and NAS 2 against 2.565 (p = 0.687). The
774-unit weighted arms are in the Figure 4E paragraph below.

| Atlas label | Program | H59 permutation q | 62-donor OLS q | 62-donor HC3 q | State |
|---|---|---:|---:|---:|---|
| Hepatocyte | Stromal ECM | 0.577711 | 0.085458 | 0.165240 | Not selected; cannot exclude the earlier effect (detection 0.794) |
| Hepatocyte | Ductular injury | 0.323375 | 0.037633 | 0.155826 | Not selected; the 62-donor call is a sensitivity only |

The H59 column is the 60-donor primary fit. The other two columns are the
October 14 diagnosis-qualified, 62-donor maximum single-cell refit, now a
sensitivity. They do not use the bulk refit; the earlier 64-donor values
are superseded by source-diagnosis correction, not by the 774-unit bulk adoption.

A supplementary contextual screen retains the two historically selected
programs and fixed 13-pair roster: 26 planned tests, with no replacement
program or pair selection. After source qualification and the unchanged
three-stage filter, 22 are estimable and four lack 30 matched donors.
Models adjust for disease stage, dataset and log sender/receiver cell counts.
Two Ductular injury associations pass both model-based and HC3 BH families
of 26: positive hepatocyte→endothelial *CDH1*→*PTPRM* (58 source-reported
donors across four datasets; β=0.495716, q=3.627665×10^-6,
HC3 q=1.167371×10^-10) and negative endothelial→hepatocyte
*ST6GAL1*→*EGFR* (n=58; β=-0.298712, q=0.006914, HC3 q=0.048800).
All four estimable leave-one-dataset-out effects retain their respective
directions: 0.003372–0.698344 and -0.566162 to -0.069645.
No ECM association passes both families. Earlier 24-estimable/one-hit/n=60
values are historical. Same-atlas RNA covariance does not establish secretion,
receptor activation, antagonism or an intercellular mechanism.

The complete source-qualified four-stage LR screen uses 69 donor identifiers,
retains all 55,402 original LR/cell-type groups, fits 349 and marks 55,053
untestable for insufficient stage counts. Its 1,047 coefficients use native
within-cell-pair/term BH correction (27 families), then Bonferroni across nine
cell-pair families. Of the fixed 13 displayed pairs, six are testable and all
have adjusted q=1; seven are untestable. Thus this is not 0/13 tested.
The old metadata join had already dropped 38 unmatched Liver Atlas run aliases
from 117 score identifiers, yielding 79; it omitted source exclusions for ten
GSE136103 donors. The corrected 69-unit fit restores those source exclusions.
The separate lineage census's 125/42 run count was not the old fitted sample.
Complete sources and checks are in the
[qualified communication candidate](../Analysis/SingleCell/candidates/source-qualified-CCC-20261014_22237840/CHECKED.txt);
the displays remain supplementary context.

**Historical working artwork, 2026-10-05; Figure 4E superseded October 14.**
Figure 4D displays the existing 117-program discovery/content
comparison (80 matches, 37 differences). Working Figure 4E now displays the
existing weighted linear bulk coefficients for all 113 eligible programs on
both fibrosis and NAS axes, with BH across 226 tests. Its 27 fibrosis and four
NAS linear calls remain conditional on historical bulk identity; they do not
establish a difference in total remodeling when nonlinear responses are allowed.
No models or program definitions were changed by this 2026-10-05 artwork update.

**Figure 4E program models on the 774 units, run 2026-10-12 (four-cohort arm
selected and installed 2026-10-14).** Source:
`figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/`
(`arm_summary.tsv`, `old_vs_new_agreement.tsv`,
`stromal_ecm_ductular_by_arm.tsv`). Figure 4E has never used the five-cohort
disease-versus-control fit. It uses the four cohorts that record both fibrosis
stage and NAS: GSE130970, GSE135251, GSE162694 and GSE174478. GSE174478 is
outside the 774 units and is not identity-resolved. Model: within each cohort,
OLS of the program score on sex, fibrosis stage and NAS with HC3 standard
errors; random-effects meta-analysis across cohorts (REML, Knapp-Hartung); BH
over eligible programs times two axes. The composition adjustment adds three
within-cohort CLR principal components.

| Arm | Records | Cohorts | Eligible programs / BH family | Fibrosis / NAS calls (q<0.05) | After composition adjustment |
|---|---:|---:|---:|---:|---:|
| Historical working panel (844-record roster) | 469 | 4 | 113 / 226 | 27 / 4 | 0 / 0 |
| 774 units only | 364 | 3 | 112 / 224 | 0 / 0 (smallest q 0.079 / 0.159) | 0 / 0 |
| 774 units plus GSE174478 | 457 | 4 | 112 / 224 | 16 / 1 | 3 / 0 |

The 774-unit rosters are GSE135251 207, GSE130970 71 and GSE162694 86 units,
plus 93 GSE174478 records in the four-cohort arm. The three composition-adjusted
fibrosis calls in the four-cohort arm are Ductular injury (0.177 [0.149,
0.205]), cholangiocyte Rho-GTPase / adhesion (0.212 [0.174, 0.250]) and
Hepatocyte stress (SLC38A1) (0.267 [0.213, 0.321]), each at BH q=0.042. Effects
agree with the working panel: Pearson 0.984 (fibrosis) and 0.973 (NAS) across
112 programs for the three-cohort arm, 0.997 and 0.998 for the four-cohort arm,
and all 27 fibrosis and 4 NAS working-panel calls keep their sign in both arms.
The four-cohort calls are not explained by the cohort count alone. All 20
(16 fibrosis, 1 NAS, 3 composition-adjusted fibrosis with the rebuilt
fractions) have tau2 = 0 and between-cohort Q of 0.09 to 1.61, below the
expected value of 3, so their untruncated Knapp-Hartung SE is 0.17 to 0.73 of
the REML SE. With the Knapp-Hartung variance floored at the REML value
(metafor test="adhoc"), both arms give 0 calls in every model (smallest q
0.193 in the four-cohort arm). A REML z-test gives 47/28/21/10 calls in the
four-cohort arm and 37/14/13/11 in the three-cohort arm; with three or four
cohorts it ignores the uncertainty in tau2, and the Figure 4E code uses
Knapp-Hartung (`scripts/manuscript/program_observability_map/analysis_lib.R`)
(`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/source_tables/fig4e_calls_by_arm_and_hk_variant.tsv`).
With the earlier cell-type fractions the four-cohort composition-adjusted
count is 1 fibrosis and 0 NAS instead of 3 and 0 (H85), and the three-cohort
arm stays at 0 and 0 under both fraction sets. Across the 112 programs, adjusted
effects correlate between the two fraction sets at Pearson r 0.935 (fibrosis) and
0.988 (NAS), and the within-cohort shuffle negative control passed (100 shuffles,
seed 20260923; `Analysis/Hypotheses/H85-20261015T124643Z/run-20261015T130425Z/P2_between_set_agreement.tsv`).
The prespecified positive control, Ductular injury, keeps a positive adjusted fibrosis
effect with a nominal 95% CI above zero under both sets (earlier 0.229 [0.135, 0.322],
rebuilt 0.177 [0.149, 0.205]), but it reaches BH q < 0.05 only with the rebuilt fractions
(q 0.161 with the earlier set), and it is subject to the Knapp-Hartung qualification
above. The old records without
GSE174478 (376 records, three cohorts) give 1 fibrosis and 0 NAS calls. Fibrosis effects for the two prespecified programs (per fibrosis
stage, 95% CI, BH q):

| Program | Arm | Unadjusted | Composition-adjusted |
|---|---|---|---|
| Stromal ECM | Working, 469 | 0.339 [0.258, 0.420], q=0.029 | 0.189 [0.024, 0.354], q=0.226 |
| Stromal ECM | 774 only, 364 | 0.320 [0.103, 0.538], q=0.214 | 0.170 [-0.056, 0.396], q=0.421 |
| Stromal ECM | 774 plus GSE174478, 457 | 0.326 [0.211, 0.440], q=0.045 | 0.146 [0.039, 0.253], q=0.230 |
| Ductular injury | Working, 469 | 0.389 [0.340, 0.437], q=0.010 | 0.230 [0.138, 0.321], q=0.118 |
| Ductular injury | 774 only, 364 | 0.376 [0.331, 0.422], q=0.079 | 0.186 [0.159, 0.212], q=0.126 |
| Ductular injury | 774 plus GSE174478, 457 | 0.375 [0.351, 0.399], q=0.004 | 0.177 [0.149, 0.205], q=0.042 |

The 2026-10-13 Figure 4 candidate used the four-cohort arm (774 units plus
GSE174478, named in the legend), which the main session recommended and the
user told it to proceed with; both arms had been computed on 2026-10-12 before
that choice. The review-round-3 candidate Figure 4E (R13-03, R12-01) uses the
three-cohort 774-unit arm (0 fibrosis and 0 NAS calls), and the four-cohort
arm is a sensitivity whose calls need GSE174478 and the untruncated
Knapp-Hartung variance (candidate
`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4e/`). The installed Figure 4E uses the three-cohort 364-unit arm: no fibrosis or NAS program reaches BH q<0.05 (smallest q 0.079 and 0.159); the four-cohort arm is a sensitivity.

The previous Figure 4E transported 46 mapped Stromal ECM and 21 Ductular injury
member genes into bulk fibrosis-stage contrasts. At F4 versus F0, 36 and 18
genes meet BH q<0.05, representing 79.3% and 86.4% of mapped program weight;
weighted directional agreement is 98.2% and 100%, respectively. These are
member-level summaries, not program-level significance tests. On the
corrected-controls fit these F4 summaries are unchanged and all eight program
stage effects change by less than 0.001
(`figures/candidates/fig4e-corrected-controls-20260925T175215Z/fig4e_old_vs_new.tsv`). The TF analysis
remains supplementary: 268 TFs were tested within each F-versus-F0 contrast,
with 35, 50, 44, and 85 meeting BH q<0.05 in F1, F2, F3, and F4, respectively.
The source table retains all 268 TFs and complete contrast-specific BH families.

The corrected candidate is recorded at
`Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/`.
Do not reconstruct its calls from older donor-collapse or four-program tables.
It predates the ambient-RNA sensitivity summarized above and therefore cannot
by itself establish hepatocyte origin.

The 117 memberships and weights are fixed; no external outcome may trigger
module rediscovery, reweighting, or renaming. Final cell-composition numbers
must be recomputed against the 1,257-sample corrected fragment substrate.

Forty-nine programs carry at least 5% lncRNA L1 weight and therefore required a
leave-all-lncRNA-out sensitivity in the original single-cell scoring framework.
The exact stored Hotspot score was retained, and only the gene-removal delta
from the source-reconstructed graph was added. Reconstructed all-gene scores
correlated at least 0.99896 with stored scores. This 0.995 technical floor was
finalized during source reconstruction and is not the biological acceptance
threshold. At the biological-donor level, **48 of 49** programs retained both
Pearson and Spearman correlation of at least 0.90 and preserved the stage-effect
direction. The minimum donor Pearson and Spearman correlations were 0.9688 and
0.9442.

The sole failure was cholangiocyte module 28: removing four lncRNA members
(11.8% of L1 weight) changed its small stage coefficient from +0.0164 to
-0.00325. It remains in the complete registry but cannot be described as
lncRNA-independent. Both prespecified programs (Stromal ECM and Ductular injury) passed this
lncRNA-removal sensitivity: Stromal ECM
had donor Pearson 0.9987 and Ductular injury had donor Pearson 0.9979,
with direction preserved for both. This sensitivity tests dependence of the
program score on its lncRNA members; it does not establish that any lncRNA is
functional.

### Validated noncoding plot package

The read-only candidate
`figures/candidates/noncoding-resource-2026-08-12/plot-upgrade-v3/` contains
nine individual PDFs and nine source tables. Independent checksum, page,
Ghostscript, source-cardinality, and scientific-census validation passed, as
did manual visual review. Repository-wide tests and the Resource scope
validator also passed in SLURM job `19779385`. The earlier `v1` and `v2` candidates are explicitly
rejected for source-table serialization and visual-label defects, respectively.

The added composition panel makes the correct null visible: lncRNAs constitute
26.74% of all 23,370 tested genes and 26.86% of the 1,616 interval-null-positive genes
on the 2026-08-10 fit. On the corrected-controls fit they are 27.0% of 23,702
tested genes and 30.9% of the 1,725 canonical positives
(`fig3-lncrna-corrected-controls-20260925T195928Z/source_tables/fig3_lncrna_role/summary.tsv`).
On the adopted 774-unit fit they are 26.7% of 24,115 tested genes and 29.6% of
the 1,854 canonical positives
(`fig3-lncrna-corrected-controls-r8-caliper-search-20261009T173834Z/source_tables/fig3_lncrna_role/summary.tsv`);
the interval-null-positive share is quoted on the 23,370-gene basis only (not rebuilt on
the 774 units).
This is comprehensive inclusion, not preferential lncRNA enrichment. The
program panels show all 117 original lncRNA-weight fractions and the complete
49-program sensitivity family; 48 preserve the donor-level stage direction.
The Figure 6 panel visualizes five deterministic molecular-object routes, four
of which are assigned to records (r8 holds no regulatory DNA element object), and
contains no corrected genetics outcome or experimental follow-up screen result.

### Validated noncoding-DNA candidate and plots

**Coordinate qualification, 2026-10-13.** Credible-set instance summaries below
remain descriptive. Locus-balanced and context-family summaries use historical
sentinel parsing and must not be treated as corrected-coordinate results. The
restored coding/splice point estimate is given under "Protein-altering fraction
of fine-mapped risk" below; other coordinate-dependent comparisons and their
uncertainty were unavailable at the initial recheck. Corrected-coordinate
reruns are now adopted below; the table immediately following
retains its historical locus-balanced and context values.

The read-only DNA-only candidate
`GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/`
uses the read-only fine-mapping substrate and does not read partial corrected
COLOC outcomes. Of 6,467 reliable Tier-1/2 study-by-locus credible-set
instances, 6,415 (99.20%) have raw SuSiE PIP sums within 0.90–1.10. Their
42,480 study–credible-set member rows collapse to 18,539 unique variants.

The complete consequence table was generated through the official Ensembl
GRCh37 REST VEP endpoint after orienting every stored allele pair against the
GRCh37 reference sequence. The first pass that assumed a fixed upstream allele
order is explicitly rejected; it misclassified the PNPLA3 I148M sentinel. The
accepted pass correctly classifies both PNPLA3 I148M and TM6SF2 E167K as
protein-altering and assigns 100% of normalized PIP mass to a mutually
exclusive consequence category.

| Candidate quantity | Direct MASLD/liver fat | Liver enzymes |
|---|---:|---:|
| Reliable credible-set instances | 207 | 6,208 |
| Mean protein-altering PIP mass | 15.77% | 3.26% |
| Mean canonical-splice PIP mass | 0.04% | 0.00% |
| Mean synonymous/UTR PIP mass | 8.46% | 5.11% |
| Mean other-noncoding PIP mass | 75.73% | 91.63% |
| Locus-balanced protein-altering PIP mass | 22.37% | 3.39% |
| Locus-balanced canonical-splice PIP mass | 0.31% | 0.00% |
| Locus-balanced synonymous/UTR PIP mass | 4.71% | 3.57% |
| Locus-balanced other-noncoding PIP mass | 72.61% | 93.04% |
| Noncoding PIP mass within 2 kb of a GENCODE TSS (pre-repair; H58 repaired 13.84% / 14.67%) | 7.87% | 4.46% |
| Noncoding PIP mass in a liver ABC enhancer (pre-repair; H58 repaired 0.43% / 3.70%) | 0.23% | 1.57% |
| Noncoding PIP mass open in at least one liver lineage (pre-repair; H58 repaired 9.26% / 8.44%) | 7.32% | 4.15% |
| Noncoding PIP mass with unresolved deposited context (withdrawn 2026-10-15: failed lookup; H58 repaired 0% / 0%) | 52.11% | 60.43% |

Regulatory-context percentages overlap and must not be summed. Credible-set
instances are descriptive study-by-locus units, not independent biological
replicates. The unresolved fractions came from failed lookup of sites absent
from the context table, almost all in non-European studies (H58, 2026-10-15).
After re-annotation by position (promoter and open-chromatin flags at GRCh38
positions after liftOver of 16,324 sites, all lifted 1:1; ABC overlap at GRCh37
positions), no
noncoding posterior is unresolved, and 77.90% (direct) and 77.69% (enzyme) of
noncoding posterior lies outside the deposited promoter, open-chromatin and ABC
annotations (summed over all credible sets of each class). The DNA result
therefore does not assign a localized regulatory mechanism to most noncoding
signals.

#### Family multiplicity correction (2026-08-17)

⛔ **Historical; inference withdrawn.** The August table and September
physical-cluster intervals both used faulty sentinel positions. The October 13
coordinate recheck below initially established a restored coding/splice point
estimate only. Neither the p, q and "Survives" columns here nor the September
intervals support the restored-coordinate analysis. Corrected bootstrap and
context reruns are now adopted in the corrected-coordinate section below.

Six direct-versus-enzyme comparisons are made on the same groups. They are one
family and are corrected together. Source:
`GWAS/finemapping/results/fig2_multiplicity/20260817T161128Z/`. Percentages here
are cluster-weighted and so differ from the instance-weighted table above.

| Test | Direct | Enzymes | Ratio | Permutation p | **BH q** | Survives |
|---|---:|---:|---:|---:|---:|:--|
| Coding PIP mass (protein-altering + splice) | 22.54% | 3.44% | 6.55 | <1e-4 | **6.0e-4** | yes |
| Unresolved deposited context | 39.13% | 73.55% | 0.53 | 2.0e-4 | **6.0e-4** | yes |
| Promoter-proximal (within 2 kb of a TSS) | 14.57% | 3.72% | 3.92 | 0.0021 | **0.0042** | yes |
| Resolved, no context | 38.52% | 18.70% | 2.06 | 0.0057 | **0.0085** | yes |
| Open in at least one liver lineage | 9.80% | 5.43% | 1.81 | 0.185 | 0.222 | no |
| **Liver ABC enhancer** | 0.10% | 0.93% | 0.11 | **0.398** | **0.398** | **no** |

⛔ **The ABC enhancer depletion is retired.** Historical bootstrap p=2e-4 and
permutation p=0.398 do not establish a depletion, and both used the faulty
cluster positions. ABC fractions are tiny and heavily zero-inflated; the
restored-coordinate context comparison below is an adopted descriptive
estimate with sparse ABC overlap. Direct-trait loci may not be described as
biologically ABC-depleted.

Two methodological points this correction established, both binding on future
work in this section.

1. **Which null.** The historical p is a bootstrap p,
   `2*min((diff<=0).mean(), (diff>=0).mean())`, which asks whether a difference
   clears zero given within-arm sampling variability. A claim about trait class
   needs a permutation null that shuffles the class label. The two disagree in
   both directions across this family: ABC passes the bootstrap and fails the
   permutation; coding does the reverse under a mis-specified statistic. Both are
   reported in the historical table with BH applied to permutation p. These
   tests are not current; the September review found no defensible trait-label
   exchange when physical clusters are shared between classes.
2. **Which statistic.** A permutation statistic of `|log(ratio of means)|` is
   invalid on these zero-inflated masses. 74.1% of enzyme clusters carry zero
   coding mass, so a random subset lands 7× below the pooled mean easily and never
   7× above; over 20,000 shuffles, 100% of exceedances were low-side. Use the
   difference of fractions.

The groups are **operational 1-Mb clusters** of credible-set sentinel positions,
not independent loci: no LD is consulted in forming them, a 1-Mb window can span
more than one independent signal, and long-range LD can link clusters.

**Historical ascertainment analysis, 2026-08-17; not a restored-coordinate
comparison.** Direct clusters in that analysis sit at twice the
standardized effect size of enzyme clusters (median β_std 0.0368 versus 0.0183)
and come from studies with roughly one seventh the effective sample size (21,893
versus 160,000), because the direct arm is mostly case-control with few cases.
Matching on β_std (standardized mean difference −0.002) moves the ratio only from
6.87 to 5.32, and adjusting for it attenuates the trait-class coefficient by 2.2%
(0.1258 to 0.1230, cluster-robust over 271 clusters); adding effective sample
size raises it to 0.1464. Source:
`GWAS/finemapping/results/ascertainment_control/20260817T152705Z/`. This controls
discovery power at the lead variant only, not fine-mapping resolution, LD
reference panel, or ancestry composition.

The accepted five-panel render is
`figures/candidates/noncoding-resource-2026-08-12/noncoding-dna-v5/`.
It contains a Figure 2 consequence panel, a Figure 2 nonexclusive context
matrix, a Figure 5 six-lineage accessibility matrix, a supplementary
30-variant accessibility heatmap, and the complete credible-set noncoding-mass
distribution. Each PDF has a read-only source table. Structural validation,
sentinel checks, Ghostscript parsing, unit tests, Resource-scope validation,
and manual visual review passed in SLURM job `19783326`.

Historical August state: these DNA-only values were validated candidates,
not synchronized manuscript numbers. At that time, target biotype, final
examples and genetics-dependent Catalog entries still required a rebuild
against the August 17 COLOC release. This is not the current r2 or companion
state described below.

The title-free locus-balanced consequence composition was first adopted as the
working main Figure 2D under explicit user approval on 2026-08-13. Its source
tables retain the instance-weighted estimates, cluster-bootstrap intervals,
and 1-Mb grouping provenance. This panel-level adoption does not adopt the
corrected COLOC release or convert these values into synchronized manuscript
numbers. Its historical coordinates, fractions and caption are superseded by
the adopted October 13 correction and installed Figure 2D described below.

## Inherited genetic variation and genetic–state overlap

**COLOC r2 adopted, 2026-10-14.** The selected genetic input is
`GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/`; the August directory
is retained as history. Independent checks confirm all 50 studies, 948,503
records unique under the study/chromosome/gene-ID key and all five recorded output hashes.
The correction reran 1,307 previously positive gene-study rows across 32
studies; other negative or unavailable results and the other 18 studies were
carried unchanged. Of the corrected rows, 1,233 remain positive, 49 are eligible at PP.H4 <=0.5 and 25 have insufficient shared posterior and are untestable: in every SuSiE signal pair, at least one signal kept less than 0.5 of its posterior on the variants the GWAS and eQTL fits share (the coloc 5.2.3 `overlap.min` default, recorded per pair by `GWAS/finemapping/src/coloc_signal_eligibility.R`). The primary tier-1/2 set contains 436 genes, ABF contains 822,
ABF-only contains 560, their union contains 996, and the full-portfolio SuSiE
set contains 708. Seven tier-1/2 genes retain an explicit priority SuSiE-untestable state; these genes also have finite lower-support results in other studies.
The current method-availability recount (job 22239724) finds finite ABF for
all 948,503 master records, of which 947,903 have an identified gene; 600
chromosome-specific anonymous `abf_only` records are excluded from gene counts.
Finite SuSiE is present in 21,799 rows (2.298% of all master rows), spanning
4,869 stable Ensembl genes. The 35 primary studies have 664,022 master rows
(663,602 identified), 12,892 finite SuSiE rows (1.942%) and 3,918 genes with a
finite result. Tier 1 has 771 finite SuSiE rows across 392 genes; Tier 2 has
12,121 across 3,888. These availability counts include 20,517 numerically
finite carried negative rows, which were not freshly shared-posterior-qualified.
The freshly qualified finite subset has 1,282 full-portfolio and 805 primary
rows. The 25 failed pairs touch 21 genes, and 11 primary failed pairs touch
nine genes; other studies give finite results for each. After remaining
positive calls are prioritized, 17 full-portfolio and seven primary genes
retain the untestable state. The full portfolio has 1,181 ABF-positive genes
and a 1,494-gene SuSiE-or-ABF union. No current memberships change in this
availability recount. The exact master-key and missing-ID reconciliation is
in `GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/method_availability_22239724/`.

**SuSiE-COLOC simulation calibration (H80, 2026-10-15): stopped, no number.**
The prespecified simulation stopped at stop rule S2 before any simulated value
was computed. Its reproduction control refit 20 stored liver-eQTL SuSiE fits
from observed z; 19 matched within the 1e-4 PIP bound and one
(ENSG00000090612) differed by up to 0.0269. Production eQTL fits were
two-stage: genes not converged at 500 iterations were refit at 5,000 (the
2026-08-15/16 refit), which reproduces that gene to 6.3e-9, while the
simulation fitted only the first stage. The statement that common-variant
colocalization has limited power (`fig2.md`) stays qualitative; no simulated
power or false-positive rate may be quoted. Source:
`Analysis/Hypotheses/H80-20261015T130802Z/RESULT_H80.md`.

Figure 1C, Figure 2A/B/C/E and Figure 5D now use their checked r2 PDFs and
source tables. Figure 2F/G/H/I are render-identical to the r2 candidates and
retain their artwork with synchronized source tables. Figure 2D retains the
independent October 13 coordinate correction. The broader Resource release
remains a separate decision. Checks:
`GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/review.json`
and `figures/candidates/r2-panel-adoption-check-20261014/panel_comparison.json`.

**774-unit values (2026-10-09).** These use the adopted 774-unit bulk refit and
COLOC release `susie_coloc_r2_20260924T053244Z`, adopted on 2026-10-14. Source:
`RNA-seq/results/manuscript_release/candidates/evidence-classes-corrected-r8-20261009T142027Z/orthogonality_audit.tsv`
(SuSiE PP.H4 > 0.5, all traits).

| Quantity | 774-unit value | Qualification |
|---|---:|---|
| Primary SuSiE-nominated genes | 436 | COLOC r2; the 7 SuSiE-untestable genes are not positive and are excluded from the overlap background below |
| Background evaluated by multi-signal COLOC | 3,769 | Native stable-Ensembl interface genes with SuSiE state positive or tested at PP.H4 ≤ 0.5 (H57 U_eval). Excludes 10,928 genes with only an ABF result (935 canonical DEGs) and 7 SuSiE-untestable genes. Distinct from 14,704 (native interface), 14,720 (symbol-grouped) and 14,682 (Catalog) |
| Genetic genes | 426 | SuSiE PP.H4 > 0.5 |
| Disease-state genes in that background | 220 | Canonical threshold; 1,155 in the 14,704 interface |
| Overlap with disease-state set | 26 | Against 21.21 expected within 10 AveExpr quantile bins; MH OR 1.295 (95% CI 0.835 to 2.007), exact two-sided P 0.2755 (H57, one test). Expression not fully controlled: a list of 426 non-colocalized genes matched on AveExpr carried 13 DEGs against 21.3 expected (O/E 0.610); no direction stated. The CI includes 1; ORs 0.46 to 1.77 had less than 80% power to be detected. Unstratified: 24.87 expected, Fisher OR 1.055, nominal P 0.83 |
| Non-overlap | 93.9% | `(426 − 26) / 426` |
| Non-overlap expected without stratification | 94.2% | `1 − 220 / 3,769`. Two sets of these sizes are mostly disjoint whether or not the biology differs; do not present non-overlap as a finding. The withdrawn 14,720-row values were 33.4 expected, OR 0.76, 92.2% |

Rows above are from H57
(`Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/results.json`);
the audit file named above holds only the withdrawn 14,720-row values.
SuSiE-or-ABF in those 14,720 rows, without expression stratification, gave 959
genetic, 57 overlap, 94.1% non-overlap, OR 0.73, nominal p 0.022 (one of several
definitions; no multiplicity correction applied); these values do not use the
3,769-gene basis of the table. As a diagnostic, the single-signal ABF map (ABF
PP.H4 ≥ 0.5, 788 genes) in the 14,704-gene interface within 10 AveExpr quantile bins gives 49 shared against
43.6 expected, MH OR 1.15 (95% CI 0.85 to 1.57), nominal exact P 0.41 (H57 D8). The Figure 3G interface
(`figures/candidates/fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_matched_treat_summary.tsv`)
uses 14,704 jointly testable genes, 426 genetically anchored, 1,155 canonical
disease-state and 26 in both; under interval-null 38 of the 426 are calls, and the
expression-matched paired comparison gives OR 1.32, nominal p 0.32 over 426
pairs. The canonical overlap odds ratio in the interface universe is not in that
table and is not quoted.

**Signed TWAS-disease agreement (H60, 2026-10-15; candidate, Reading 3).**
ALT TWAS z came from OTTERS (P0.05 model, Broadaway liver-eQTL weights, UKBB
ALT GCST90019492); this OTTERS run is dated 2026-03-24 and is not part of any
adopted release. Disease t came from the adopted 774-unit refit (774
genetic-identity clusters, 138 control and 636 disease; not verified donors).
Across 14,703 genes the AveExpr-decile-stratified Spearman rho was 0.0107
(200-block jackknife 95% CI −0.0059 to 0.0272; within-chromosome circular-shift
P 0.1994, 10,000 shifts; family of one, two-sided alpha 0.05; seed 20260923; 0
of 200 jackknife blocks change the sign). Genome-wide signed rank agreement
between ALT-TWAS z and the disease t is bounded inside ±0.03 without
composition adjustment. With hepatocyte fraction in the design, rho 0.0167
(CI 0.0001 to 0.0332; shift P 0.0449); this interval is not inside ±0.03, and
because the primary test did not reject, no direction is stated. The bound
covers only genome-wide rank agreement: it excludes sign agreement of about 0.60 or
more among the 1,527 genes with |z| > 3, but not full sign agreement among the
103 |z| > 3 canonical genes (expected rho about 0.0225). D5, |z| > 3: rho
0.0421 (95% CI −0.0075 to 0.0916), sign agreement 0.5141 against 0.4994
expected; D6: 54 of 103 same sign (0.524). Cohort-omission fits range 0.0094
to 0.0123, none with the opposite sign. The COLOC-excluded step was not tested
(rho 0.0057, n 9,089), so agreement beyond colocalized loci is not shown. All
predicted controls passed and stop rules S1 to S8 did not trigger; no control
tests the real TWAS-to-disease pairing. The Figure 1B and Figure 3G numbers
are unchanged. Source:
`Analysis/Hypotheses/H60-20261015T112434Z/out-final-20261015T120955Z/RESULT_H60.md`.

**Kumak gene lists (H71, 2026-10-15; null).** Of 944 NAS and 1,013 fibrosis
Kumak genes matched to the Catalog, 258 and 301 are bulk-supported; in the
bulk/COLOC interface 661 of 902 and 714 of 943 have no multi-signal COLOC
result. Among the 3,769 evaluated genes (241 NAS and 228 fibrosis list genes),
NAS: 36 colocalized against 28.905 expected within 10 AveExpr bins, MH OR
1.322 (95% CI 0.910 to 1.920), exact P 0.1474, Holm-adjusted P 0.2949;
fibrosis: 28 against 24.885, MH OR 1.159 (0.764 to 1.759), exact P 0.5056,
Holm-adjusted P 0.5056 (family of two, Holm). The expression-matched
random-list negative control (NC-1, 1,000 lists per Kumak list) missed its
prediction: mean stratified observed/expected 0.931 for NAS (outside the
predicted 0.95 to 1.05) and 1.000 for fibrosis (inside it), and the fraction of random lists with P < 0.05
was 0.000 for both (predicted 0.02 to 0.08). Both values lie inside the S6 stop
limits, so S6 was not triggered. The colocalization rate did not
differ detectably from that of genes of the same expression level; planted
shifts of OR about 1.74 and 0.48 (NAS) and 1.75 and 0.48 (fibrosis) were
detected in 82% and 65%, and 77% and 62%, of simulations. Do not write
"similar", "the same", "consistent with" or "orthogonal". Source:
`Analysis/Hypotheses/H71-20261015T124239Z/run-20261015T125815Z/RESULT_H71.md`.

The next table is historical. Its values are measured on the **corrected COLOC release**
(`GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/`, review items 2 and 5)
and the **844-record corrected-controls bulk fit**, in the evidence-class universe (bulk and
atlas genes;
`RNA-seq/results/manuscript_release/candidates/evidence-classes-corrected-20260925T175215Z/orthogonality_audit.tsv`).
They were candidates and are superseded by the 774-unit table above:

| Quantity | Corrected value | Qualification |
|---|---:|---|
| Primary SuSiE-nominated genes | **436** | coloc's shared-posterior eligibility check applied; 26 of the adopted 462 lost, 0 gained |
| SuSiE-untestable tier-1/2 genes | 7 | ADAM15, ASAP3, CKB, GSTM2, LILRB3, PSMC5, RITA1; kept in denominators as not positive, not as negatives |
| Jointly testable background | **14,243** | Genes assayed in both the bulk contrast and the colocalization universe. The Figure 3G interface uses a different universe (14,183 on corrected inputs) |
| Jointly testable genetic genes | **409** | SuSiE PP.H4 > 0.5 within the background |
| Disease-state genes in that background | **1,082** | Canonical threshold, restricted to the jointly testable background |
| Overlap with disease-state set | **24** | Against **31.1 expected** (409 × 1,082 / 14,243); Fisher OR 0.75, nominal p = 0.22 |
| Non-overlap | **94.1%** | `(409 − 24) / 409` |
| Non-overlap expected by chance | **92.4%** | `1 − 1,082 / 14,243`. **Two sets of these sizes are mostly disjoint whether or not the biology differs; do not present non-overlap as a finding** |
| Union of named genetic evidence | **996** | Tier-1/2 SuSiE-or-ABF union |
| Tier-1/2 ABF-nominated genes | **822** | ABF values reproduce the adopted master |
| Full-portfolio SuSiE genes, all 50 studies | 708 | Adopted 745 |

The following values are historical, measured on the **COLOC release adopted
on 2026-08-17** (`GWAS/finemapping/results/susie_coloc/`) with the 2026-08-10
bulk fit. Both the July and August genetic memberships are superseded by r2:

| Quantity | Adopted value | Qualification |
|---|---:|---|
| Tier-1/2 GWAS strata | 35 | Liver-specific/direct and accepted proxy strata; the registry was not changed by the rerun |
| Primary SuSiE-nominated genes | **462** | **Adopted release (2026-08-17); supersedes the July value of 473** |
| Jointly testable background | **14,920** | **Recomputed 2026-08-17**; genes assayed in both the bulk contrast and the adopted colocalization universe. Rebuilt from the adopted release, not carried over from July's 14,931 |
| Jointly testable genetic genes | **437** | **Recomputed 2026-08-17** from the adopted 462; supersedes 447 |
| Disease-state genes in that background | **1,066** | **Recomputed 2026-08-17**; canonical threshold, restricted to the 14,920-gene jointly testable background. Supersedes the July value of 1,079. Only this restricted count enters the overlap test |
| Overlap with disease-state set | **20** | **Recomputed 2026-08-17**; canonical threshold (`padj<0.05`, `\|log2FC\|>0.50`), against **31.2 expected**; Fisher OR 0.616, nominal p = 0.037 |
| Non-overlap | **95.4%** | `(437 − 20) / 437`; version- and threshold-specific, not a universal biological constant |
| Non-overlap expected by chance | **92.9%** | `(437 − 31.2) / 437`, i.e. what independent sampling gives for a 437-gene set against 1,066 of 14,920. The observed 95.4% is 2.6 points above it. **Two sets of these sizes are mostly disjoint whether or not the biology differs; do not present non-overlap as a finding** |
| Overlap / non-overlap under the interval-null comparator | 33 / 92.4% | `(437 − 33) / 437 = 92.4485%`; the retired canonical, retained as the calibrated arm. **Producer located 2026-09-07**: `RNA-seq/results/genetics_state_interface/20260817T162108Z_promoted_coloc_raw_treat/matched_treat_by_arm.tsv`, `n_case_treat = 33` at `treat_lfc = 0.25` on the `susie_pp4_0.5`/`official_27638` arm. Its `input_manifest.json` carries the same three sha256 inputs as the canonical run, so the two arms are on one release. The earlier "92.5%" was a rounding of 92.4485% |
| Union of named genetic evidence | **1,013** | **Adopted release (2026-08-17); supersedes the July value of 1,030** |
| Tier-1/2 ABF-nominated genes | **822** | **Adopted release; supersedes 862** |

The corrected union and direct-disease scopes: SuSiE-or-ABF 930 genetic, 56
overlap, 94.0% non-overlap, OR 0.77 (nominal p 0.06); direct-disease SuSiE-or-ABF
150 genetic, 10 overlap, 93.3%, OR 0.87 (nominal p 0.88). The interval-null-comparator row
was not recomputed in the evidence-class universe; in the Figure 3G interface
universe 34 of 406 genetic genes are interval-null calls
(`20260924T061535Z_coloc_r2_corrected_controls/matched_treat_by_arm.tsv`).

Adopted-value source: `RNA-seq/results/genetics_state_interface/20260817T133951Z_promoted_coloc/`
(job `20871359`), which rebuilds the background from the adopted release rather
than reusing July's, and reports the official 27,638-gene read-count substrate
and the 23,370-gene five-cohort fragment substrate side by side without pooling
(23,370-gene, 844-record basis; not rebuilt on the 774 units).

**Do not report a depletion.** On the 774 units with COLOC r2 (adopted
2026-10-14), restricted to the 3,769 genes evaluated by multi-signal COLOC and
stratified by 10 AveExpr quantile bins, no depletion is detected: 26 shared
against 21.21 expected, MH OR 1.295 (95% CI 0.835 to 2.007), exact P 0.2755
(H57; expression not fully controlled, so no direction is stated). The earlier
unstratified test gave OR 0.76, nominal p = 0.20 in the 14,720-row
evidence-class universe (OR 0.757, 33.46 expected, in the 14,704-gene
interface); both counted genes without a multi-signal COLOC result as
negatives (10,941 and 10,928). On the 844-record corrected inputs it gave OR 0.75,
nominal p = 0.22 in the evidence-class universe and OR 0.72,
nominal p = 0.15 in the Figure 3G interface universe; on the adopted inputs it
was OR 0.616, nominal p = 0.037. Nothing here establishes that colocalized genes
are buffered or constrained: the expression-matched ladder is flat, and the
threshold-free SuSiE estimate on corrected inputs is rho = 0.021 (p = 0.20). The
following ABF check is on the adopted inputs and was not rerun. The ABF
threshold-free correlation looks strong (rho = −0.066, p = 8.6e-16) and is an
expression artifact — partialling out `AveExpr` leaves rho = −0.016 (p = 0.051),
and on absolute effect it vanishes (p = 0.34). See `EXPRESSION_CONFOUND_CHECK.md`
in the run directory; report the SuSiE arm and the matched ladder instead.

**Post-rerun correction (2026-08-17).** The corrected 50-study × 22-chromosome
COLOC rerun is COMPLETE and PROMOTED to canonical. Measured directly from
`GWAS/finemapping/results/susie_coloc/gene_level_coloc_tier12.csv`, named genes
only (one blank-symbol row excluded, as before):

| | July | post-rerun | Δ |
|---|---:|---:|---:|
| Tier-1/2 SuSiE > 0.5 | 473 | **462** | −11 |
| Tier-1/2 ABF > 0.5 | 862 | **822** | −40 |
| Tier-1/2 SuSiE-or-ABF union | 1,030 | **1,013** | −17 |
| full-portfolio SuSiE > 0.5 | 736 | **745** | +9 |

Superseded by the 2026-09-24 corrected release (review item 2,
`susie_coloc_r2_20260924T053244Z/r2_summary.json`): tier-1/2 SuSiE 436, ABF 822,
union 996, full-portfolio SuSiE 708.

The tier-restricted counts FALL while the full portfolio RISES — the recovered
eQTL fits are spread across all 50 studies, whereas tier-1/2 draws on 35. The
SuSiE arm churns substantially: **147 tier-1/2 genes change their call (79 out,
68 in)**, driven by retaining ~14.1% of variants that were previously dropped at
every locus as strand-ambiguous. All 7 QA guards pass against the adopted
canonical, including C8 positive/negative control recovery.

⚠ **Downstream status is now split.** Main Figure 2 A-I has been rederived
from the adopted release and explicitly adopted as canonical. Still stale
against the adopted release as of 2026-08-17: the atlas rebuild chain
(`RNA-seq/run_atlas_rebuild_chain.sh`),
`scripts/figures/fig4a_overview_candidates.py`, evidence-class sizes,
genetics-dependent Gene Catalog fields, other worked-example consumers, and
the ~30 single-cell/spatial/ATAC/perturbation consumers enumerated in
`docs/archive/progress/2026-08-17_coloc_rerun_promoted_to_canonical.md`. The
`relayreg`/`orientsub` products are stale for a different reason: they ran on
2026-08-15 against the rerun before the 145 recovered eQTL fits were merged, so
they completed silently on a portfolio headline of 712. Rebuilt on 2026-08-17:
the portal ancestry export
(`masld-atlas-v2/public/data/coloc_*`), the Figure 1 GWAS panels, and canonical
Figure 2 A-I.

COLOC adoption authorizes the genetic counts and overlap arithmetic above; the
separate round-2 approval authorizes the canonical Figure 2 A-I examples. It
does not authorize evidence-class sizes or Gene Catalog entries. On the corrected
inputs, evidence classes and a corrected Figure 2 A-I candidate
(`figures/candidates/fig2-coloc-r2-20260924T053244Z/`) now exist, and the corrected
Gene Catalog release, whose adoption the user approved on 2026-09-25, is built
(see "MASLD Gene Catalog").

### Program lineage labels versus program content, all 117, 2026-09-07

Source: `Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/`;
reasoning in that directory's `FINDINGS.md`. The registry is unchanged — no
membership, weight or name was altered, and no disease outcome was read.

A program's lineage name records the cell-type subset Hotspot was run on, not a
content test. Scoring every program's fixed member genes against donor-level
per-lineage pseudobulk (runs collapsed to donors; Hepatocytes 70, Endothelial
66, Fibroblasts 56, Cholangiocytes 51, T cells 39, Macrophages 35), and
z-scoring each gene across lineages so the index reads specificity rather than
abundance:

**80 of 117 programs (68%) have content matching their label; 37 (32%) do not.**

The index was validated against the donor-paired same-atlas contrasts before
anything else was read: it independently returns Fibroblasts for Stromal
ECM/IGFBP7 and Cholangiocytes for Ductular injury/BICC1, matching beta 2.271 and
1.574. The producing script exits non-zero if it does not.

The mismatches are systematic, not scattered. **17 of the 37 are programs
discovered in a non-hepatocyte subset whose content is hepatocyte**, and their
own names say so: "Hepatocyte xenobiotic (CYP2C)" labelled fibroblast, "Plasma
protein (APOA1)" labelled macrophage and separately fibroblast, "Bile-acid
transport (ABCB11)" and "Plasma lipoprotein (APOE)" labelled cholangiocyte.

The mechanism is ambient RNA and it is measured. Against the existing all-117
decontX burden table, the 17 non-hepatocyte-labelled programs with hepatocyte
content have median ambient burden **0.432** (IQR 0.361-0.589) against **0.061**
(IQR 0.041-0.084) for the 60 non-hepatocyte programs whose label matches;
Mann-Whitney one-sided p < 1e-5, rank-biserial **+0.998**. Across all 117,
ambient burden and hepatocyte content correlate at rho = +0.19 (p = 0.040).

Consequences: Figures 4C-4E must label programs by content with the discovery
subset retained in the legend as provenance; the two featured programs are the
general case rather than exceptions; and the vocabulary null (469-record roster: 26 of 114
fibrosis-supported against 18.4 expected, p = 0.134; 364-unit Figure 4E arm:
1 against 1.24, p = 0.473; 457-record sensitivity arm: 20 against 26.9,
p = 0.871) now has a cause, because a
third of the vocabulary does not measure what its name says. This does not
rescue the registry on fibrosis.

### Are the two gene classes distinguishable? Corrected 774-unit comparison, 2026-10-14

The corrected comparison uses the adopted 774-unit RNA fit and COLOC r2,
with the original expression-matching rule, outcomes and seeds unchanged.
Its stable-ID joint universe contains 14,704 genes: 426 genetic positives,
1,155 canonical RNA calls and 26 shared, hence 400 genetic-only and 1,129
state-only genes. Of the 1,129 state-only genes, 935 have no multi-signal COLOC
result and 194 were tested at PP.H4 ≤ 0.5 (H57 D6). A post hoc rerun on the 193
tested state-only genes that are single-cell observable gave 129 pairs, mean
difference −0.342 (paired sign-flip nominal p 0.0154, 100,000 draws, seed
20260923; sign split 57/72, p 0.218; minimum detectable paired d 0.247); among
the 281 pairs below, the 56 with a tested state-only partner have mean −0.748
(p 0.00148) and the 225 with an untested partner −0.145 (p 0.233)
(`Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave3/c7/20261015T185127Z/c7_results.tsv`,
arms C and D). Single-cell observability retains 399 genetic-only and
1,124 state-only genes; 281 pairs satisfy the unchanged 0.05-SD expression
caliper. The protein comparison has 73 matched pairs. Genes and matched pairs
are the units of this comparison, not independent patients or independent
genetic loci.

| Outcome | Matched pairs | Mean genetic-minus-state difference | Post hoc paired sign-flip nominal p |
|---|---:|---:|---:|
| Hepatocyte specificity | 281 | −0.265515 log2 units | 0.01482985 |
| RNA-effect-adjusted protein stage response | 73 | −0.015601 | 0.487465 |

Both expression-balance checks pass (SMD 0.00371 and −0.00062). The prespecified
direction, greater hepatocyte specificity among genetic genes, fails: the
observed direction is opposite. The hepatocyte difference attenuates from the
historical 272-pair mean −0.435 and nominal sign-flip p=2×10⁻⁵. On the
corrected pairs the positive/negative split is 126/155 (nominal sign-test
p=0.094677), and nominal Wilcoxon p=0.008919. The sign-flip diagnostic uses
100,000 draws, seed 20260907, conditional on the realized matching. It is
post hoc and does not correct dependence among genes, shared loci or gene
sets; it is not confirmatory evidence of biological separation or opposition. Flipping signs jointly within the 216 H10 loci that hold the 281 matched genetic genes gives nominal p=0.01606 (100,000 draws, seed 20260907), and 99 of 216 locus means are positive (sign test p=0.247); both are post hoc (`Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Claims/R14-10/r14_10_recompute.log`).
No difference is detected in the residualized protein arm; this does not
establish equivalence. The raw protein comparison remains circular and is
excluded from the scientific claim.

The prespecified separate negative control is now recorded: 1,000 within-pair
label swaps per declared outcome, seed 20260907, without rematching or new
fits (job 22236840). The simulated means are −0.001071 and −0.000441,
with Monte Carlo SE 0.003442 and 0.000721; analytic expected means are zero.
No new acceptance tolerance was added. This diagnostic does not calibrate the
different relabel/rematch null or change the failed directional prediction.

The original relabel-and-rematch procedure is retained for numerical
reproduction only. Its two-outcome BH values are 0.00159984 and 0.459654
(nominal p=0.00079992 and 0.459654). The previously identified covariate-null
mismatch remains: rematching after relabeling removes the actual expression
imbalance. These p/q values must not be quoted as calibrated significance.
The primary-arm rematch-null SD is 0.076056 versus paired sign-flip SD
0.108878. The hepatocyte-versus-fibroblast positive control passes at 2.6749
log2 units.

Source: [corrected class comparison](../figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/class_difference_results.tsv),
[paired diagnostic](../figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/corrected_null_diagnostic.tsv),
and [matching eligibility](../figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/matching_feasibility.tsv).
Compute job 22236311 completed; the prespecification hash is unchanged.

### Historical class comparison, 2026-09-25, superseded by the 774-unit comparison

**Historical 844-record analysis.** The 272-pair hepatocyte-specificity test and
71-pair protein test below used the September corrected-controls bulk fit,
not the adopted 774-unit refit. Its estimates are superseded by the
2026-10-14 comparison above and retained only as historical results.

Source: `RNA-seq/results/class_difference_test/cdt-corrected-20260925T203546Z/`
(`VERDICT.md`, `results/class_difference_results.tsv`,
`results/corrected_null_diagnostic.tsv`, `results/matching_feasibility.tsv`), job
`21905239`. Same sealed prespecification as the 2026-09-07 run
(`prespec/PRESPEC.md`, sha256 `f5213ae6...`, verified at runtime); only the
substrate changed: the corrected COLOC release and the corrected-controls bulk fit.
Class counts match the corrected Figure 3G interface (14,183 genes; 406
colocalized; 1,077 disease-state; 23 both).

Membership overlap cannot carry a claim, because it cannot fail: 92.4% of a
409-gene set is expected to miss a 1,082-gene set in a 14,243-gene background by
arithmetic alone, against 94.1% observed (844-record basis; on the 774 units
with COLOC r2 (adopted 2026-10-14), among the 3,769 genes evaluated by
multi-signal COLOC, 94.2% expected without stratification against 93.9%
observed for 426 and 220). This test asks instead whether the
classes differ in something measurable.

Each colocalized gene was matched 1:1 to a disease-state gene on bulk `AveExpr`
(deterministic nearest neighbour, hardest-first, 0.05 SD caliper). Genetic-only
genes average `AveExpr` 4.16 against 0.47 for state-only genes (difference 3.69
log2 units; SMD 1.40), and 272 of 382 genetic-only genes (71%) have a partner
inside the caliper.

| | hepatocyte specificity (primary) | protein stage response, residualised (secondary) |
|---|---:|---:|
| matched pairs | 272 | 71 |
| SMD after matching | +0.0034 | -0.0020 |
| mean paired difference | **-0.435** | -0.031 |
| median paired difference | -0.441 | -0.008 |
| sign split (positive / negative) | 119 / 153, p = 0.045 | 35 / 36, p = 1.0 |
| Wilcoxon p | 3.7e-5 | 0.33 |
| paired sign-flip null p (100,000 draws) | **2.0e-5** | 0.18 |
| prespecified permutation p | 1.0e-4 (floor of 10,000 draws) | 0.14 |
| minimum detectable paired d | 0.170 | 0.333 |

**Verdict: the prespecified prediction is not supported, and the classes differ in
the opposite direction.** The prediction recorded in advance was that colocalized
genes are *more* hepatocyte-specific than their expression-matched disease-state
partners. Among 272 matched pairs they are *less* hepatocyte-specific (mean paired
difference -0.435; paired sign-flip p = 2.0e-5; Wilcoxon p = 3.7e-5). Because this
direction was not prespecified, it is a post hoc result. The protein-level arm is
null (71 pairs; prespecified permutation p = 0.14, sign-flip p = 0.18). The raw,
un-residualised protein comparison is circular for the reason given below (d =
-0.43, p = 2.0e-4) and is not reported as a result. No mechanism is inferred from
the hepatocyte-specificity difference.

**Superseded 2026-09-07 run.** Source:
`RNA-seq/results/class_difference_test/cdt-20260907T105813Z/`, job `21489259`, on
the 2026-08-17 COLOC release and the 2026-08-10 bulk fit. Its values below are
kept as the earlier record; the same negative direction was present but not
distinguishable from zero under the paired null (p = 0.064).

Each colocalized gene was matched 1:1 to a disease-state gene on bulk `AveExpr`
(deterministic nearest neighbour, hardest-first, 0.05 SD caliper). Matching is
mandatory, not a sensitivity: the classes differ in expression by SMD **1.576**
(mean `AveExpr` 4.048 against -0.285), because colocalization requires a
well-powered eGene. Only **242 of 415** colocalized genes have a partner inside
the caliper at all.

| | hepatocyte specificity (primary) | protein stage response, residualised (secondary) |
|---|---:|---:|
| matched pairs | 242 | 54 |
| SMD after matching | +0.0053 | -0.0006 |
| mean paired difference | **-0.2035** | -0.0371 |
| median paired difference | +0.0096 | -0.0242 |
| sign split | 123 / 119, p = 0.85 | 24 / 30, p = 0.50 |
| Wilcoxon p | 0.088 | 0.288 |
| paired sign-flip null p | **0.064** | 0.154 |
| minimum detectable paired d | 0.180 | 0.381 |

**2026-09-07 verdict (superseded): the classes were not distinguishable on either
axis at that sample size.** The prespecified rule additionally required the direction to match the
prediction recorded in advance (colocalized genes *more* hepatocyte-specific);
the observed sign is the opposite, so "separable" is not claimed on any reading.

Two things must travel with this result.

- **The prespecified permutation null is anti-conservative and its p = 0.0036
  must not be quoted.** Relabelling genes and rematching destroys the covariate
  imbalance that defines the real comparison: under relabelling both
  pseudo-classes have the same expression distribution, matching becomes easy,
  and differences shrink. On the 2026-09-07 pairs that null had sd 0.0684 against the
  statistic's own standard error of 0.1091 — **1.59x too narrow**. The correct
  paired null randomises which member of each realised pair carries the label,
  preserving the matching, the matched-share selection and the expression
  asymmetry; it had sd 0.1097 and gave p = 0.064 on the 2026-09-07 pairs. On the
  corrected pairs the relabelling null sd is 0.077 and the paired sign-flip null
  sd is 0.107; the sign-flip p is the one to quote.
- **The un-residualised protein comparison is circular.** Raw, it gives
  d = -0.554, p = 0.0002. Disease-state genes were *defined* by their bulk RNA
  response and protein tracks RNA, so that difference confirms the selection.
  Residualising on `|bulk logFC|` removes it (p = 0.154). Do not report the raw
  form.

Positive control passed: hepatocyte markers separate from fibroblast markers by
2.67 log2 units on the same index (predicted > 2.0 in writing beforehand), so
the null is a real null rather than a broken measurement.

**Language consequence.** Do not describe the two classes as "separable",
"biologically opposed", or "complementary" as a finding. Report the prespecified
directional test as failed and the opposite-direction hepatocyte-specificity
difference as post hoc. The largest difference between the classes is liver
expression level: 3.69 log2 units of mean `AveExpr` between genetic-only and
state-only genes on the corrected inputs (4.3 on the 2026-09-07 inputs).

### Program-pair covariance in the GSE193066 holdout, 2026-09-08 (historical; identities unqualified)

Source: `RNA-seq/results/gse193066_holdout/holdout-20260907T105813Z/`, prespecification
`PRESPEC_A_coordination.md` (sha256 `fd00249637a1aa90...`) hashed before the cohort was
read and verified at runtime; reasoning and caveats in that directory's `VERDICT.md`.

GSE193066 was sealed for this test and had never been opened. The four
communities, the 2,426 coordinated edges, the 113 program definitions and the
scoring function were taken unchanged from discovery; **nothing was refit**.
Analysis set 105 participants, first biopsies only — the holdout is 160 samples
but 106 participants, and 58 contributed a repeat biopsy, so a sample-level
analysis would have pseudoreplicated 58 people.

| statistic | observed | null mean | p |
|---|---:|---:|---:|
| edge sign concordance, 2,426 coordinated edges | **0.8974** | 0.7095 | 0.001 |
| fixed-partition modularity vs label permutation | **0.2519** | -0.0064 | 0.005 |
| fixed-partition modularity vs degree-preserving rewiring | **0.2519** | -0.0054 | 0.005 |
| ductular-portal community density, lineage-matched null | **0.8762** | 0.6809 | 0.001 |

All three predictions recorded before opening were correct, including that
holdout modularity would sit below the discovery value of 0.358.

⚠ **Quote the effect, not the p, for the concordance statistic.** Its null is at
**0.709, not 0.5** — non-coordinated program pairs also agree in sign 71% of the
time — so the result is a 19-point excess over that baseline. The null is also
narrower than an independent one, because each of the 1,000 draws takes 2,426
edges from only 3,769 available non-coordinated edges and draws therefore share
about 64% of their content. A binomial test at 0.5 would have declared success on
the baseline alone and was forbidden by the prespecification for that reason.

GSE193066 participant identities are unqualified (K01), so this comparison is
historical and does not count as replication. In GSE213621, H52's program-group
permutation control, which removes covariance between programs, still gave a
mean excess of 0.0251 over the non-coordinated baseline (95th percentile
0.0809), above the prespecified limits; the concordance excess quoted here is
the same statistic. n = 105, one cohort, one platform, and GSE193066 holds
essentially no F4 (one sample).

### Program-pair covariance in further cohorts and in tissue (H54, H52, H88), 2026-10-15

Title clause. "The coordination of multicellular remodeling" is removed from the
title under the joint H52/H54 rule, row H52 untestable / H54 supports
(`Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/JOINT_RULE_H52_H54.md`;
`Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/TITLE_CLAUSE_DECISION.md`). Program-pair covariance is
reported in Results as a discovery-cohort result with the replications below.
H52 and H54 use different edge sets (2,328 pairs rediscovered on 371 refit
units; the 2,426 fixed T4 pairs) and different nulls; their numbers are never
pooled. GSE193066 counts toward neither. Neither test addresses "multicellular":
bulk program scores mix cell abundance, biopsy content and program activity.

Discovery (fixed T4, `interaction_summary.tsv` and `cohort_blocks.tsv` under
`BULK-PROGRAM-MAP-v9/systems/T4_interaction/`): 469 source records (GSE130970
76, GSE135251 214, GSE162694 86, GSE174478 93), 113 programs, 6,195 pairs in the
testing family, 2,426 coordinated beyond the leading residual component (650
same-lineage, 1,776 cross-lineage). Modularity of the 4 fixed communities is 0.358 against null means of 0.138 (200 rank-matched label permutations) and 0.089 (200 degree-preserving rewirings); one community holds 15 programs, 9 fibroblast and 6 cholangiocyte (`interaction_summary.tsv`, `modularity_null.tsv`, `program_communities.tsv` in the same directory).

**H54, two further bulk cohorts**
([RESULT](../Analysis/Hypotheses/H54-20261015T111917Z/RESULT_H54.md)). Verdict:
replicated in both primary sets at k = 1 and k = 3; the F0-F1 secondary sets
also pass. Unit: participant, one biopsy each. Design: sex + fibrosis stage +
NAS bin. Statistic: fraction of coordinated pairs with a positive partial
correlation; excess = that fraction minus the mean of 200 within-lineage
program-label permutations; interval = 1,000 participant-bootstrap draws within
fibrosis stage. p* is the largest of the four permutation and bootstrap p
values; Holm over A91 and B87 (primary) and separately over A66 and B53
(secondary).

| Set | n | Fraction positive, k = 1 (random-placement mean) | Excess k = 1 [95%] | Excess k = 3 [95%] | p* | Holm |
|---|---:|---|---|---|---:|---|
| A91, GSE281797 | 91 | 0.8478 (0.6901) | 0.1577 [0.1323, 0.1818] | 0.1867 [0.1177, 0.1884] | 0.000999 | rejected, primary |
| B87, GSE225740 | 87 | 0.8336 (0.5220) | 0.3116 [0.2698, 0.3235] | 0.2407 [0.1864, 0.2419] | 0.000999 | rejected, primary |
| A66, GSE281797 F0-F1 | 66 | 0.8478 (0.6821) | 0.1658 [0.1320, 0.1891] | 0.1844 [0.1039, 0.1938] | 0.000999 | rejected, secondary |
| B53, GSE225740 F0-F1 | 53 | 0.8293 (0.5257) | 0.3036 [0.2511, 0.3106] | 0.2236 [0.1594, 0.2347] | 0.000999 | rejected, secondary |
| A37, A66 without "No Pathology" | 37 | 0.8590 (0.7143) | 0.1446 [0.1048, 0.1791] | 0.1716 [0.1060, 0.1796] | descriptive | outside family |

A91 and A66 print the same fraction positive because their sign flips cancel
exactly (104 flips at k = 1, 52 each way; 112 at k = 3, 56 each way); the
underlying z values differ (k = 1 z correlation 0.988 over 6,194 finite pairs;
recomputed from `stage2/edge_z_A91.tsv.gz` and `stage2/edge_z_A66.tsv.gz`). Limits stated with every reading: participant overlap
with the discovery cohorts was not checked by genotype in either cohort;
GSE225740 identities cannot be checked (no reads, no deposited sex, age or
batch); the GSE281797 composition arm failed its acceptance rule (median
hepatocyte fraction 0.4387, accepted 0.60 to 0.95) and was not read; in
GSE225740 the result held after adding three composition components, log library
size and top-100 count share (k = 1 excess 0.2237, descriptive; estimated
composition is not measured composition); a size-matched
label null also rejected in both sets; within-cohort TMM moved the GSE193066
fraction by 0.0202 (0.8974 to 0.9176); 29 of the 66 A66 participants have
diagnosis "No Pathology". Controls: the code reproduced the GSE193066 fraction
0.897361912613355 and comparator 0.709472008490316 (differences below 1e-15);
the label null detected the GSE193066 signal (p 0.000999, null mean 0.6061);
gene permutation gave mean fractions 0.4851 (A91) and 0.4827 (B87);
informative-null pairs were not above the null (p 0.530 and 0.978); the label
test rejected in 0 and 1 of 200 gene-permuted datasets (limit 16). Power gate:
MDE at k = 1 0.0831 (A91) and 0.0817 (B87) against the GSE193066 excess 0.2915.
No numbered stop rule fired.

**H52, GSE213621**
([RESULT](../Analysis/Hypotheses/H52-20261015T111454Z/out/RESULT_H52.md)).
Verdict: untestable, stop rule 5. The 371-unit coarse rediscovery reproduced the
fixed 2,328-pair set exactly, and the leave-GSE162694-out positive control
passed (sign concordance 0.9361 against 0.7487, 109 units). Negative controls on
permuted GSE213621 data, 100 permutations each, fail if the mean excess exceeds
0.02, its 95th percentile reaches 0.05, or the mean sign-matched excess exceeds
0.02: gene permutation passed (mean −0.0284,
95th percentile 0.0099); program-group permutation failed (mean 0.0251, 95th
percentile 0.0809). No GSE213621 outcome statistic was computed.

**H88, HMSMA Visium**
([RESULT](../Analysis/Hypotheses/H88-20261015T130436Z/out/RESULT_H88.md)).
Verdict: Bulk-sign only. Unit: donor (32 arrays, 9 female and 23 male); 44
programs measurable in all 32. T1 (160 coordinated against 372 non-coordinated
cross-lineage pairs sharing no genes): excess 0.0257 Fisher z, label-null 95%
range of D −0.0107 to 0.0147, Holm p 0.0002, MDE 0.0182, donor-bootstrap
interval 0.0197 to 0.0326, positive in 31 of 32 donors; program relabeling
within lineage p 0.0024 (holds). T2 (134 against 80 pairs with positive bulk
covariance): excess 0.0116, p 0.1240, MDE 0.0280; program relabeling p 0.1373
(fails). Of the 372 non-coordinated pairs in T1, 292 have bulk r_disjoint of
zero or below. The four secondary variants kept the T1 excess (ratios 0.87 to 1.34;
Holm p 0.0004 each). Controls: zonation reproduced (mean z −0.1705, negative in
31 of 32), planted difference 0.03 detected in 0.998 of draws, T1 false-positive
rate 0.038, spot-shuffle |D| at most 0.0006, portal-tract co-location mean z
0.0680 (positive in 28 of 32). Reading: "Program pairs with positive bulk
covariance co-vary across spots, but coordinated pairs do not exceed other
positively covarying pairs (T2 0.0116, p 0.1240; an excess of 0.0280 or more
would have been detected with 80% power)." This is not spatial evidence for
fixed condition 1 of the joint rule.

### The severity instrument does not yet demonstrate within-person tracking, 2026-09-08

**Historical source-declared-pair diagnostic, qualified 2026-10-14.** GSE193066
identities remain unqualified. The fitted numbers below are retained for
historical description only, not primary within-person tracking, ICC or
label-error calibration.

Same directory, prespecification `PRESPEC_B_instrument_change.md`
(sha256 `fee244640b27afee...`). The label-free diagnostic was computed and
reported before any Kleiner value was joined.

Across the 58 paired biopsies, within-person change in latent severity and in
recorded Kleiner stage correlate at **Spearman 0.236** (P = 0.074; bootstrap 95%
CI **-0.002 to 0.457**), against a minimum detectable rho of **0.361** fixed
before the estimate was read. Sign concordance among the 30 participants whose
stage actually changed is 19 of 30 (P = 0.20). **Verdict: inconclusive, reported
as a bound** — not evidence the instrument fails and not evidence it works.

In the 28 pairs whose recorded stage did not change, the score still moved with
SD 0.383, against 0.340 across all pairs. At this scale score movement does not
carry stage movement. Three constraints were stated in advance: the release axis
is kallisto and GSE193066 is STAR; all 58 source-declared pairs were high-risk
at first biopsy; and serial biopsy variability includes real change and
sampling. The historical SD inputs 0.4655 and 0.6433 do not establish the
label's error or make a one-stage change interpretable as true change.
**The instrument remains a cross-sectional measure.**

### Which trait class drives the colocalized map, corrected 2026-09-24

Source: `GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/`
(`r2_summary.json` and `gene_level_coloc_tier12.csv`). The SuSiE split uses the
SuSiE driver (`driving_trait`) and the ABF split uses each gene's own ABF driver
(`abf_driving_trait`, review item 5).

| set | genes | strongest evidence from a liver-enzyme trait | from direct disease / liver fat |
|---|---:|---:|---:|
| SuSiE PP.H4 > 0.5 (primary) | 436 | **423 (97.0%)** | **13 (3.0%)** |
| ABF PP.H4 > 0.5 | 822 | 705 (85.8%) | 117 (14.2%) |

Primary set by trait: ALT 137, AST 158, GGT 128; NAFLD 8, PDFF 4, NASH 1. ABF set
by trait: ALT 246, AST 253, GGT 206; NAFLD 88, PDFF 20, NASH 9.

**PDFF input qualification, verified 2026-10-13.** SuSiE and COLOC r2 consumed
the coordinate-correct AF-augmented Pazoki PDFF copy. Correcting reversed
allele labels left its non-allele columns unchanged. The two double-lifted
hg19 copies were not inputs to these pipelines; the uniform sign reversal
does not change their unsigned posteriors. C2orf16 and RABGEF1 are the two of
436 primary genes with SuSiE support involving this study, and both also pass
the primary threshold in another study. The blanket coordinate uncertainty
for these inputs is withdrawn. [Unverified] Direction-sensitive and model
outputs beyond the inspected pipelines remain untraced. Execution and copy
checks are linked in the
[verified review](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#review-of-this-review-2026-10-13).

**Eight previously missing PDFF genetic correlations completed, 2026-10-14.**
The corrected allele-orientation input now gives native LDSC estimates for
all eight pairs. The installed LDSC code reused an unfiltered allele vector
after filtering incompatible variants; prefiltering by the exact native
allele rule avoids that error without changing Z scores or allele labels.
Duplicated rsIDs (12,123 in each input) are excluded rather than arbitrarily
choosing one record; 118,338–125,501 incompatible or ambiguous shared IDs
are then excluded. Each pair retains approximately 1.15 million unique
allele-compatible IDs before the LD-reference intersection.

| PDFF pair | rg | SE | Nominal p |
|---|---:|---:|---:|
| MVP NAFLD | 0.7389 | 0.0687 | 5.7694e-27 |
| MVP cirrhosis | 0.6155 | 0.0891 | 5.0011e-12 |
| FinnGen cirrhosis | 0.5890 | 0.1719 | 0.0006 |
| UKBB ALT | 0.4979 | 0.0466 | 1.2033e-26 |
| MVP ALT | 0.3790 | 0.0512 | 1.4121e-13 |
| UKBB AST | 0.1860 | 0.0493 | 0.0002 |
| MVP AST | 0.2219 | 0.0522 | 2.1106e-5 |
| FinnGen obesity | 0.2556 | 0.0498 | 2.9104e-7 |

FinnGen cirrhosis has observed-scale h2 0.0026 (SE 0.0012); its correlation
is imprecise and does not establish a precise bounded relationship. The other
seven traits also have positive estimated heritability. These summary-GWAS
estimates use the native LD-block jackknife; SNP counts are not biological
sample sizes. No correction of the original matrix hypothesis family is
claimed. Jobs 22236252 and 22236390 completed; all 256 allele-rule combinations
matched installed LDSC, no LDSC exceptions remained and source hashes were
unchanged. Source:
`GWAS/ldsc/rg_pdff_compatible_20261014T140601Z_22236252/pairwise_estimates.tsv`.

**Fixed 12-trait genetic-correlation family completed, 2026-10-14.**
Job 22237000 repaired the remaining 27 non-PDFF allele-compatibility failures;
job 22237084 checked all native outcomes, unchanged input hashes and positive
estimated observed-scale heritability in both traits for all 66 pairs. No pair
remains computationally missing. The versioned terminal matrix preserves 28
successful July non-PDFF fits and the three October 12 corrected-sign PDFF fits,
adds eight October 14 PDFF and 27 non-PDFF repaired fits, and records each source.
The repaired fits use native allele-compatible joins and exclude all duplicated
rsIDs; retained historical fits preserve their original preprocessing. This is
a candidate matrix, not a uniform 66-pair refit or a new installed figure.
P values remain nominal in the fixed 66-pair family. Six unconstrained estimates
exceed one, all involving Ghodsian NAFLD or FinnGen cirrhosis. For example,
PDFF–Ghodsian NAFLD 1.1925 (SE 0.2360) and MVP NAFLD–FinnGen cirrhosis 1.0263
(SE 0.2543) are retained with uncertainty, not clipped or interpreted as true
correlations above one. Weak Ghodsian NAFLD and FinnGen cirrhosis heritability limits precision.
Source: `GWAS/ldsc/rg_missing_compatible_20261014T150911Z_22237000/rg_12trait_66pair_terminal.tsv`.

Superseded values (`GWAS/finemapping/results/coloc_trait_scope/cts-20260907T105813Z/`,
adopted 2026-08-17 release): SuSiE 448 / 14 of 462 (ALT 150, AST 164, GGT 134;
NAFLD 8, PDFF 5, NASH 1); ABF 758 / 65 of 823. The ABF split used the SuSiE driver
whenever a finite SuSiE result existed and included one blank gene identifier, so
it did not describe each gene's strongest ABF evidence.

⚠ **Scope.** `gene_level_coloc_tier12.csv` holds one row per gene (18,976 rows in
the corrected release), and `driving_trait` names the trait behind that
gene's **strongest** colocalization. This therefore measures which trait class
supplies each gene's best evidence. It does **not** show that a gene has no
weaker colocalization with the other class; that needs a per-gene-by-trait
posterior matrix, which this file does not contain. Do not write "exclusively"
or "only".

Read with the protein-altering result below, these numbers are descriptive.
Pooled across ancestries, direct-trait credible sets carry a larger share of
posterior probability in coding variants than liver-enzyme credible sets; H58
(2026-10-15) found this pooled contrast confounded with ancestry and
fine-mapping coverage, with no trait-class call in any stratum. The strongest
colocalization evidence for most genes comes from a liver-enzyme trait, and far
more enzyme than direct-trait gene-by-GWAS pairs were tested (12,121 against
771). Some coding effects may
act through protein function without changing steady-state RNA; others change
RNA abundance through nonsense-mediated decay or splicing
([Rivas et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25954003/)). The overlap
between the genetic and disease-state gene sets is a separate gene-level
estimate. Any role of the coding difference in that overlap is a hypothesis:
this comparison does not hold gene eligibility, phenotype, study power or
expression observability fixed, and did not test it. The retired "complementary,
largely non-overlapping maps" framing is not restored.

**Sex difference at enzyme signals (H72, 2026-10-15;
`Analysis/Hypotheses/H72-20261015T134833Z/run-20261015T142309Z/RESULT_H72.md`).**
Pass with heterogeneous genes. At the shared liver-enzyme signal of 13 of 301
genetically anchored or concordant genes, the per-allele effect differed between women
and men on both the raw and the standardized scale (BH q < 0.05 over 281 SNP tests per
scale; UK Biobank, 184,990 women and 159,146 men for ALT). No gene was bounded within a
female-to-male ratio of 0.67 to 1.5, and 271 could not be bounded. The 13 genes share 9
SNPs; 12 have the larger effect in men and THBS3 (ALT) in women. Thirteen further genes
were different on one scale (BH q < 0.05) and not detected on the other; this design
cannot tell a scale effect from lower power on the second scale. Four genes (SAR1A,
TYSND1, HKDC1, DYRK2) are scale-dependent: significant on both scales, but the raw
per-allele effects have opposite signs in women and men, with one sex near zero, while the
standardized effects have the same sign. The prespecified sentence for this state does not
describe them, so they are reported only by state. Of the 426 Catalog r8 genes (400
anchored, 26 concordant), 123 had no enzyme colocalization above PP.H4 0.5 (not
applicable) and 2 had no usable UK Biobank variant (FCGR2B, CFHR1). The median per-SNP
detectable difference was 0.0242 SD per allele on the standardized scale. Positive
controls (urate) and genomic-control lambda (0.884 to 1.061) passed. Averaged over GGT
loci, enzyme-locus effects were proportional to each sex's trait SD; this re-analyses
public UK Biobank sex-stratified GWAS already analysed by Flynn et al. (2020) and
Bernabeu et al. (2021) and is not a finding of this paper. Nothing here says that a gene
acts differently in MASLD. Gene membership comes from COLOC r2 and Catalog r8; the
proposed `enzyme_sex_effect_state` column is an annotation that must be rebuilt for any
later Catalog release.

### Protein-altering fraction of fine-mapped risk: coordinate recheck, 2026-10-13

The October 13 review independently restored sentinel positions against native
credible-set member coordinates. The chromosome.position decimal encoding
truncates trailing zeros: GCKR 2.2773094 was read as chromosome 2, position
2,773,094 rather than 27,730,940. This mislocates 715 of the 6,415 gated sets.
Every gated set has exactly one trailing-zero restoration within 1 Mb of all
its native members. Inputs:
`GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_pip_architecture.tsv`
and `credible_set_members_annotated.tsv` in the same directory. The PIP-sum gate
and within-cluster weighting are unchanged; the independent calculation did
not import either affected Figure 2 script.

| Quantity | Restored-coordinate point estimate | Historical September calculation; faulty coordinates |
|---|---:|---:|
| Gated credible-set instances | 207 direct / 6,208 enzyme | Same |
| Physical 1-Mb clusters | 233 total; 19 direct / 228 enzyme / 14 shared | 264 total; 24 / 259 / 19 |
| Protein-altering plus canonical-splice posterior fraction | 23.3428% direct / 3.8732% enzyme | 22.68% / 3.39% |
| Ratio | **6.026671** | 6.69 |
| Difference of fractions | **0.194695** | 0.192906 |

The descriptive direction persists. The initial recheck established the point
estimate without new uncertainty. The shared parser repair and completed
reruns below now supply adopted corrected bootstrap and robustness outputs.
The historical 62-comparison menu becomes 57 because there are 19 rather than
24 direct clusters to omit. No restored-coordinate permutation p is reported;
named historical releases are preserved. Corrected Figure 2D and its
noncoding-context companion were regenerated from the adopted tables and
installed on October 13.

#### Adopted corrected-coordinate results, 2026-10-13

Job `22218816` completed the six-comparison bootstrap family, full robustness
menu, descriptive trait-directness tables and separate ascertainment sensitivity
in `figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/`.
The user adopted these corrected Figure 2 results on October 13. The original
output directory is retained without rewriting historical analyses. All 21 focused checks passed,
including an independent clustering/weighting regression over 6,415 gated
credible sets and 715 restored positions. All 57 recorded input/snapshot
checksum entries were unchanged before and after the run.

| Quantity | Adopted corrected result | Unit and qualification |
|---|---|---|
| Coding/splice posterior fraction | 23.3428% direct / 3.8732% enzyme | 207 / 6,208 gated credible-set instances, weighted within operational 1-Mb clusters |
| Coding clusters | 233 physical; 19 direct / 228 enzyme / 14 shared | Pooled across both trait classes; no claim of LD independence |
| Ratio and 95% bootstrap interval | **6.026671 [2.574593, 13.510165]** | Multiplicity producer, 10,000 physical-cluster draws; a sampled shared cluster carries both classes; pooled across ancestries, descriptive only, no trait-class call (H58) |
| Difference and 95% bootstrap interval | **0.194695 [0.069229, 0.339940]** | Difference of posterior fractions, same draws |
| Noncoding-context clusters | 227 physical; 18 direct / 223 enzyme / 14 shared | Context subset clustered separately using coordinates recovered from all members |
| Corrected robustness menu | 57 comparisons, ratios **4.368132–8.275691** | 1 reference +13 direct-study +18 enzyme-study +19 direct-cluster omissions +5 windows +1 LD-block comparison |
| Smallest difference-interval lower bound | **0.022373** | All 57 difference intervals above zero; each has 10,000 finite draws; these are correlated sensitivities, not 57 independent replications |

The [primary bootstrap table](../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/fig2_family_multiplicity.tsv)
and [complete robustness comparisons](../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/robustness/coding_mass_robustness.tsv)
retain all estimates. Seed 20260817 is combined with a stable test/fit identifier.
The robustness reference uses a different identifier and gives ratio interval
2.58449–13.55838; this is a separate bootstrap stream, not a conflicting point
estimate. Permutation p/q are absent in the primary family. Bootstrap tail
diagnostics and their BH columns are not confirmatory hypothesis tests.
The legacy trait-directness table also carries bootstrap tail diagnostics;
these must not be reported as calibrated trait-class p-values.

**Ancestry and subclass stratification (H58 and review round 3, 2026-10-15).** The ancestry stratification (H58 family C) was adopted on 2026-10-17 as the Figure 2D ancestry companion (`figures/main/fig2_genetics/panels/Fig2D_coding_by_ancestry.pdf`); the direct-subclass split is supplementary (Figure S2S, `figures/main/fig2_genetics/panels/FigS2S_direct_subclass_split.pdf`). Neither is part of the adopted ratio table above, and no number changes. Ancestry composition: 106 of 207
direct and 82 of 6,208 enzyme credible sets are European, and 15 of 19 direct
clusters contain a European direct set. H58 family C (10,000 shared draws over
81 pooled clusters, seed 20260923; 98.333% intervals, Bonferroni over three):
EUR+AMR 0.1497 [0.0388, 0.2718] (direct 20.98% vs enzyme 6.00%), EUR 0.2296 [0.0550, 0.4236] (29.88% vs 6.92%), AMR 0.0219 [-0.0489, 0.1065] (6.73% vs 4.54%) (95% intervals shown); no call for any member. All non-European studies pooled (10 direct and 224 enzyme clusters) give 6.99% against 3.68%, difference 0.0331 [-0.0383, 0.1274], a descriptive value with no calibration control (`Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009/descriptive_refits.tsv`, row nonEUR_pooled). EUR failed
coverage (enzyme 8 of 767 lead windows by locus identifier, 23 of 767 within 500 kb by position (`Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009/G7_coverage.tsv`); direct 47 of 47), power (MDE 0.304 >
0.195) and calibration (N1a: 181 and 188 of 200 at 95% and 98.333%); AMR failed
calibration (N1b: 141 and 150 of 200); EUR+AMR carries the EUR coverage failure
and failed calibration (N1c: 173 and 184 of 200). Direct subclasses against
all enzyme clusters (3.87%): PDFF 40.14% (37 sets, 12 clusters; difference
0.3626 [0.1835, 0.5455]) and NAFLD+NASH 12.44% (170 sets, 13 clusters; 0.0857
[-0.0164, 0.2096]); 10,000 physical-cluster draws, seed 20260923. Direct-cluster
concentration: the top 5 of 19 direct clusters carry 76.0% of the summed direct
cluster fractions and one PDFF credible set at chr11:823,586 carries 22.5%;
ratios after dropping the direct rows of the top 1 to 5 clusters are 4.93,
4.15, 3.45, 2.69 and 1.96 (a biased diagnostic, not a test). The ungated
ascertainment medians below (effective n 21,893 vs 160,000; standardized effect
0.0368 vs 0.0183) cover 208 and 6,259 credible sets, including the 52 that fail
the PIP-sum gate. Comparing the trait classes within one ancestry requires SuSiE
fine-mapping of all 767 UKBB ALT, AST and GGT lead windows and of the MVP
European ALT, AST and NAFLD GWAS; this has not been done. Figure 2A support
count (H58): of the 436 multi-signal genes, 100 have a colocalized lead within
250 kb of a same-study fine-mapped window (89 unpadded; 190 within 250 kb of
any window), and 94 are supported only by the three MVP European strata,
which have no credible sets.

The five adopted context differences (promoter 0.138659, ABC -0.005811,
accessible 0.059563, resolved-no-context 0.226098, unresolved -0.413643) are
withdrawn: H58 (2026-10-15) showed the unresolved mass came entirely from failed
lookup, almost all of it in non-European studies (summed noncoding PIP 3,604.4
from sites absent from the lookup table and 47.8 from other allele order, of
which 0.27 came from European rows). After position re-annotation (16,324
sites, 1:1 liftOver) the pooled cluster-weighted differences are promoter
0.0995 [-0.028, 0.248], ABC -0.0338 [-0.052, -0.017], accessible -0.0312
[-0.102, 0.046], resolved-no-context -0.0896 [-0.229, 0.038], and unresolved 0.
Pooled values mix ancestries; within AMR (99% intervals, Bonferroni over five)
only resolved-no-context passed calibration (-0.0957 [-0.344, 0.126], no
difference detected); promoter, ABC and accessible received no call. Context categories overlap and must not be summed.

The [ascertainment sensitivity](../figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/ascertainment/ascertainment_control_summary.tsv)
preserves its historical **ungated** design: 6,467 instances, 19 direct and
234 enzyme clusters, baseline ratio 6.178788 and effect-matched ratio 5.782896
(400 matched instances, matching SMD -0.001927). All three adjustment models
ran on 6,466 effect-observed instances and 239 pooled clusters. Native lead
position lookup was already correct; this repair changes sentinel clustering.
Its old permutation diagnostics are not adopted trait-class tests and its
ungated fractions are not the primary gated estimand. Adoption preserves these
inference limits. Implementation and checks are in the
[rerun record](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#coordinate-parser-repair-and-reruns-2026-10-13).

The installed [Figure 2D](../figures/main/fig2_genetics/panels/Fig2D_pip_architecture_by_trait_directness.pdf)
prints protein-altering and other-noncoding fractions of 23.1% / 73.5% for
direct traits and 3.9% / 92.5% for enzyme traits, with 3.6% synonymous/UTR
mass labeled for enzymes. The full four-class estimates above use the
[corrected source table](../figures/main/fig2_genetics/panels/Fig2D_pip_architecture_by_trait_directness_source.tsv).
The combined coding/splice ratio and its interval are in the
[panel interval table](../figures/main/fig2_genetics/panels/Fig2D_corrected_coding_interval.tsv).
The installed [noncoding-context companion](../figures/main/fig2_genetics/panels/Fig2_companion_noncoding_annotation_gap.pdf)
remains outside A-I and displays fractions summed over all credible sets of a
class (direct / enzyme). The previous companion printed unresolved 52.11% / 60.43%,
promoter 7.87% / 4.46%, accessible 7.32% / 4.15% and ABC 0.23% / 1.57%; H58
(2026-10-15) withdrew these (failed lookup). Repaired values, same estimator,
as printed on the installed panel:
unresolved 0% / 0%, promoter 13.84% / 14.67%, accessible 9.26% / 8.44%, ABC
0.43% / 3.70%, resolved with none of these contexts 77.90% / 77.69%. Its region-balanced intervals are separate
estimates. Job 22219710 regenerated both PDFs; source values, exact page sizes,
embedded fonts, PDF integrity and rendered text passed checks, and both panels
were visually inspected. Previous Figure 2D files were preserved; other
existing main-panel PDFs are unchanged. The
[adoption and artwork record](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#figure-2-adoption-and-regenerated-panels-2026-10-13)
contains the execution details.

**Historical September uncertainty, withdrawn for current inference.** Sources:
`GWAS/finemapping/results/fig2_multiplicity/20260923T215748Z/fig2_family_multiplicity.tsv`
and `GWAS/finemapping/results/fig2_robustness/20260923T215800Z/`
(`robustness_summary.json`, `coding_mass_robustness.tsv`). They used seed
20260817 and 10,000 physical-cluster bootstrap draws, with no defensible
trait-label exchange and no reported permutation p. Their ratio interval
3.27–12.99 and difference interval 0.084–0.315, and the separate robustness
reference intervals 3.28–13.06 and 0.087–0.316, used faulty coordinates.
The `nominal_p_bootstrap_diff` and `nominal_q_bootstrap_BH` columns remain
unreported diagnostics. The following historical values are retained solely
to identify the outputs requiring replacement:

| Historical check; not valid restored-coordinate inference | fits | historical ratio range | historical lowest difference-interval lower bound |
|---|---:|---|---:|
| reference, 1-Mb clusters | 1 | 6.69 | 0.087 |
| leave-one-study-out, direct | 13 | **5.13 - 9.76** | 0.043 |
| leave-one-study-out, enzyme | 18 | 5.63 - 7.25 | 0.077 |
| leave-one-cluster-out, direct | 24 | 5.70 - 6.98 | 0.066 |
| window size 0.25 to 5 Mb | 5 | 5.30 - 6.69 | 0.068 |
| ancestry-matched LD blocks | 1 | **7.07** (27 vs 491 clusters) | 0.078 |

The historical observation that all 62 difference intervals stayed above zero
does not establish robustness after coordinate restoration. MVP_NAFLD_AMR
contributes 95 of the 207 direct credible-set instances; the historical drop
gave ratio 9.76. The historical `2022_36402844_PDFF_EUR` drop (14 instances)
gave ratio 5.13, difference interval 0.043–0.246; the LD-block comparison gave
7.07. These historical values are superseded by the corrected comparisons
adopted above.

Superseded (2026-09-07, `fig2_robustness/20260907T105813Z/`, ungated credible sets, row-level
permutation): ratio 6.549, refit range 4.99 to 9.58, worst permutation p 8.0e-4, family
BH q 6.0e-4, 95 of 208 direct instances.

### Adopted Figure 2 A-I round-2 release

**Historical August 17 panel set.** Figure 2D's coordinates, estimates, source
tables and caption were replaced by the adopted October 13 correction above;
the remaining panel sequence and release states are unchanged by that decision.

The adopted A-I set was produced in
`figures/candidates/fig2_round2_promotion_20260817T211348Z/`. Job `20879541`
rerendered all nine panels together from adopted inputs and passed
source-cardinality, exact fractional-point page-box, one-page PDF,
embedded-font, Ghostscript, prohibited-label, raster-equivalence, visual, and
Resource-scope validation. The resulting PDFs and sidecars are canonical under
`figures/main/fig2_genetics/panels/`. Panel B uses the full 0--1 axes and a
portrait plotting area measuring 1.488 by 1.970 inches; only 0.018 inch remains
between the x-axis-title baseline and the page bottom.

The variant consequence annotation was rerun directly from the adopted COLOC
aggregate in job `20871989`: 19,851 unique variants were classified, and all
five prespecified coding, intronic, promoter, and intergenic sentinels passed.
For the 1,013 named Tier-1/2 multi-signal-or-single-signal COLOC union, the
strongest lead-variant annotation per gene is:

| Consequence | Genes |
|---|---:|
| Coding | 36 |
| UTR | 73 |
| Promoter | 107 |
| Intron | 455 |
| Intergenic | 342 |

The multi-signal ancestry partition is 352 EUR-only, 70 shared with EUR, 26
non-European-only but below the lead-variant genome-wide-significance threshold, and
14 non-European-only genes passing that threshold (13 EAS, 1 AFR). These are Figure
2 source-table values, not final evidence-class or Gene Catalog memberships.

The cross-ancestry matrix is exhaustive under its displayed rule, not
hand-selected: 13 genes have PP.H4 > 0.5 in four of the five ancestry panels,
and no gene has support in all five. The roster is GGT1, PANX1, EPHA2, C2orf16,
ACTG1, SHROOM3, MLIP, MSL2, EFHD1, GOT2, GLDC, HSCB, and CHEK2.

**Corrected Figure 2A–C/E–I adopted, 2026-10-14.** Rebuilt on the
corrected COLOC release in `figures/candidates/fig2-coloc-r2-20260924T053244Z/`
(validation passes); values read from the rendered PDFs. Panel A: 436
multi-signal, 560 single-signal only, union 996; SuSiE PP.H4 bands >0.9,
0.8–0.9 and 0.5–0.8 contain 244, 89 and 103 genes. Three genes drawn as
single-signal only could not be tested by multi-signal COLOC and are not
negatives. Panel B (installed 2026-10-15 from the round-3 review repair): x is the
SuSiE PIP of the same GWAS and lead variant that gives each gene's maximum
PP.H4, taking the largest SuSiE PIP when that variant has several fine-mapping
records. Of 3,918 genes, 883 have a fine-mapping record at that variant and 866
a SuSiE PIP; 327 are shown, 18 with PIP ≥ 0.9 and PP.H4 < 0.5 and 18 with
PIP ≥ 0.9 and PP.H4 > 0.5. The three MVP European strata have no fine-mapping,
so the 1,348 genes whose maximum PP.H4 comes from them have no x. Of the 327, 62
have several records and for 50 the SuSiE PIPs differ; taking the smallest SuSiE PIP gives 326 shown and 17
with PIP ≥ 0.9 and PP.H4 < 0.5 (POLL falls below 0.99; RORA x 0.433 instead of 0.440).
THBS3 has its maximum PP.H4 (0.963) at two leads (UKBB_ALT 1:155129562 and
UKBB_GGT 1:155086571), and only the UKBB_ALT lead has a SuSiE PIP (0.051); taking
the other lead gives 326 shown with the same corner and top-right counts. The previous panel
showed 695 points, 152 and 56, with x taken as the maximum over four methods, all
strata and all leads. Panel C consequence counts: coding 35, UTR 70, promoter 104, intron 450, intergenic 337 (996; the 2026-10-17 redraw assigns intronic and intergenic by direct overlap with hg19 knownGene transcripts and prints intron 524, intergenic 263); of the 35 coding genes, 14 have the lead variant in
their own coding sequence or splice site and 21 have it only in another gene's
(29 distinct lead variants; round-3 review candidate panel
`figures/candidates/review-r3-fixes-20261015T112148Z/F2abce/Fig2C_coding_noncoding_split_ownership.pdf`). Panel E ancestry partition: 340 EUR-only (127 SuSiE-tested in at least one non-EUR study, 213 not tested in any), 64 shared,
23 non-European-only below the genome-wide-significance threshold, and 9
non-European-only passing it (8 EAS, 1 AFR). These are gene-level categories:
19 of the 32 non-European-only genes have their lead variant within 1 Mb of the
lead variant of a European multi-signal-positive gene (the same 19 with each
gene's minimum-p lead or maximum-PP.H4 lead), and the 9 passing genes form 8
loci, 4 of which (5 genes) have no European multi-signal-positive neighbour
within 1 Mb. Panel F keeps the same 13 genes;
RORA and FABP1 keep PP.H4 0.983 and 0.998. These supersede the August genetic
counts and sources above. F/G/H/I retain their render-identical artwork;
D uses the separately adopted October 13 coordinate-corrected panel.

Late supplementary Figure 2 panels J, L, M, and O--R were independently
rederived in candidate
`figures/candidates/figS2_JR_promoted_coloc_20260817T204537Z_20878109/`.
Job `20878109` passed source, PDF, embedded-font, 6-point-text, prohibited-label,
and Resource-scope checks. Panels K and N are withheld rather than carried
forward: K has 31,036 historical rows without adopted COLOC gene/posterior
linkage, and the 35 adopted GWAS strata contain no stage-associated GWAS stratum
for N. Canonical supplementary PDFs remain unchanged pending separate adoption.

Early supplementary panels A--D and F--I were rederived and independently
validated through job `20879572`; E is withheld because the adopted BBJ_GGT--ACADS
has no multi-signal support. The principal changes are EFHD1--BBJ_AST PP.H4
0.84 to 0.650043; method-concordance rho 0.5453 to 0.5134; full-portfolio
multi-signal genes 736 to 745 and single-signal genes 1,234 to 1,181; and
tier-1/2 multi-signal threshold counts 473/369/267 to 462/345/252 at PP.H4
0.5/0.8/0.9. These supplementary files remain candidates; main A--I alone was
authorized for round-2 canonical adoption.

**Nearest gene at colocalized loci (H74, 2026-10-15; supplement only).** At
114 of 259 tier-1/2 liver-trait loci with a colocalized liver eQTL gene
(SuSiE-COLOC PP.H4 > 0.5), the protein-coding gene with the TSS nearest the
GWAS lead variant of the colocalizing signal was not a colocalized gene
(E1 = 0.440, 95% CI 0.381 to 0.501). At 27 loci the nearest gene had been
tested and showed a distinct liver eQTL signal (PP.H3 > 0.8); at 86 loci it had
no multi-signal colocalization test. A colocalized nearest gene counted as
agreement even when it is not the gene with functional evidence: at 1p13 the
nearest protein-coding gene to the anchor is PSRC1 (colocalized), not SORT1;
this positive control was not blind.
Colocalized genes lay nearer the lead variant than tested genes placed at
random in the same window (observed rate 0.766 against null 0.241;
Holm-adjusted p T_A 4.8e-77, T_B 0.0001, the smallest value 10,000
permutations allow; calibration, not a finding). The
nearest gene was untested at the TM6SF2 and APOE loci; at the GCKR locus the
nearest gene (MPV17) was untested and GCKR was colocalized. Descriptive: E1
was 0.27 in the 88 loci with the fewest protein-coding genes within 500 kb and
0.61 in the 88 with the most; by the tier of the
anchor row, 7 of 7 tier-1 and 107 of 252 tier-2 loci had a nearest gene that was
not colocalized;
E1 was 0.405 with the SuSiE top SNP as anchor. Do not call the nearest gene wrong or
the colocalized gene causal. Source:
`Analysis/Hypotheses/H74-20261015T124604Z/run-20261015T125815Z/RESULT_H74.md`.

COLOC is shared-signal evidence. It nominates genetically anchored candidates;
it does not by itself identify a causal gene or prove mediation.

**Accessible-chromatin colocalization (H53, 2026-10-15; PASS; Catalog r9d
fields, adopted 2026-10-16).** Unit of description: the GWAS signal
(study, lead). Of 581 European tier-1/2 study-lead signals with a liver-eQTL
colocalization in r2, 191 share a signal with a bulk-liver accessible region
(Currin 2025 caQTL, 138 donors; coloc.bf_bf PP.H4 > 0.8 and lead LD r² ≥ 0.5),
at least 45 of them already reported by Currin 2025 for the same trait and
peak (a lower bound, because the proxy-lead match was not implemented; Broadaway
2024 not checked); 315 map to a different variant, 60 were not
tested in this source, 9 colocalized without confirmed lead LD and 6 were
indeterminate. Of the 436 primary SuSiE genes, 143 have a
colocalized accessible region; at 108 the same region also colocalizes with
the gene's eQTL (three-way, PP.H4 > 0.8) and 35 are reached only through the
r2 GWAS-eQTL pair (11 with an ABC-only predicted link) and are not called
linked to the region. Tier-1 direct signals: 6 colocalized, 16 different
variant, 2 indeterminate, 1 not tested. Primary estimate: 80 of 238
confirmatory liver loci (chr1 and chr19 excluded), 0.336 (Wilson 95% CI 0.279 to 0.398),
conditional on prior liver-eQTL colocalization; not the share of liver GWAS
loci with a caQTL. T1 (family of one): liver-only 172 against platelet-only 78
loci, adjusted OR 3.350 (profile 95% CI 1.622 to 7.397), Wald P 0.001645;
unadjusted 70/172 against 11/78 (descriptive). This agrees with Currin's
trait-level pattern and is not presented as new. Stop rules S1 to S5 passed;
shifted-caQTL negative control 0.0407 and 0.0513 without the LD rule, which
misses the prediction of 0.03 or below per arm but meets the S4 stop bound of
0.10, and 0.0233 and 0.0256 with it; NC2 0 of 292. The chromosome 1 and 19
positive controls had been examined in the pilot and are not blind; only the
chromosome 5 controls were blind. Sensitivities: 223 signals at
p12 = 1e-5; 253 at PP.H4 > 0.5 with the LD rule. The caQTL and eQTL donor sets
may overlap, so three-way support is not independent. Source:
`Analysis/Hypotheses/H53-20261015T112655Z/RESULT_H53.md`.

## Assay observability and evidence-state accounting

These results define the interpretive boundary around the inherited-variation
and multicellular-remodeling maps. Apparent negatives are assigned states that
reflect what each assay could actually establish rather than being counted as
evidence against a gene or program.

| Boundary | Result | Evidence interpretation |
|---|---:|---|
| CosMx complete-program observability | 12 / 117 programs | The other programs cannot be treated as negative whole-program tests |
| GeoMx complete-program testability | 0 / 117 programs | `untestable`: no gene axis loaded; direction and normalization unresolved; dropped from inference |
| GSE287826 complete-program testability | 0 / 117 programs | `untestable`: no gene axis loaded and no authoritative AOI-to-donor key; observability not yet determined |
| Visium/ST complete-program observability (GSE192741, Vu, HMSMA, Yakubovsky) | 114 / 117 programs each | 3 partial in each; 0 untestable |
| ATAC jointly testable peaks | 39,914 | Denominator after cohort and lineage filters |
| ATAC `indeterminate` peaks | 38,066 / 39,914 | Available cohorts do not resolve a replicated accessibility state |
| ATAC `source_dependent` peaks | 1,844 / 39,914 | BH q < 0.05 in one cohort only; 1,505 of 1,844 have the same logFC sign in the other cohort. No peak reaches q < 0.05 in both cohorts with opposite signs (0 `discordant`) |
| Catalog `tested_negative` calls | not assessed | No negative-call rule was defined for either domain; the negative-rule fields are empty in all 56,466 evidence rows, so a count of zero would be fixed by construction |
| Bulk L09-01 equivalence rule, candidate (H81 spike-in) | usable in 3 of 10 expression deciles (7, 8, 9) | False `tested_negative` at 0.75 log2FC 0.0% in those deciles; size at the margin above 0.05 elsewhere (0.0579 overall with between-cohort heterogeneity). Not adopted; the Catalog `tested_negative` state stays not assessed |

These counts are not validation failures. They show that modality-vote counting
would silently treat targeted-panel absence, missing donor joins, source
dependence, and unresolved replication as equivalent negative evidence. The
Resource instead exports the state, its reason, the unresolved alternative,
and the deterministic next-experiment rule in machine-readable form.

## MASLD Gene Catalog

### Catalog review and molecular-object companion, 2026-10-14

The paper's selected r8 Catalog uses the adopted 774-unit bulk fit and COLOC r2.
The [reviewed copy](../figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/)
preserves the exact October 9 scientific tables and all original class
memberships. Its local review completed in job `22236566`, followed by a
read-only verification. The reviewer is explicitly **Codex remaining_catalog
scientific review**, not a substituted user signature. The historical source
bundle is unchanged; no external portal was published.

| Quantity | Selected r8 value |
|---|---:|
| Gene identities | 28,233 |
| Gene evidence rows | 56,466 |
| Genetically anchored | 400 |
| Concordant (both domains supported; direction not compared) | 26 |
| Established-state associated | 1,129 |
| Unresolved | 13,127 |
| Untested | 13,551 |
| Jointly eligible Catalog identities | 14,682 |
| Prespecified programs / memberships | 117 / 6,633 |
| `tested_negative` calls | not assessed (no negative-call rule defined) |

The [local interface](../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/access_v2/index.html)
was rebuilt from these same gene calls with visible restrictions: 774
conservative genetic-identity clusters, not verified independent donors; 855
annotations counted but held untestable; no gene-level causal inference or
program-to-member-gene support. All five review images were visually inspected,
and all 20 stored examples independently reproduce from native source
thresholds. All 55 recorded input files retained their checksums.

The [molecular-object companion](../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/index.html)
contains 28,233 records from the same r8 classes. Its assay_applicability field
reads `indeterminate` for bulk RNA in all 1,853 RNA-supported records and for
chromatin and spatial context in every record; proteomics reads `indeterminate`
in 21,305 records and `not_applicable` in 6,928, and no record carries a protein
measurement. A candidate correction, not adopted, sets bulk RNA to the gene's
RNA call, proteomics to quantified or not quantified in PXD051911 (6,042 genes
quantified: 5,952 by direct symbol match and 90 after symbol remapping), and
chromatin and spatial to not joined at gene level. Protein significance is
reported in table T1 of
`Analysis/Hypotheses/H86-20261015T125018Z/run-20261015T125926Z/RESULT_H86.md` and the Figure 5B legend, not in the Catalog. Its count matrix, annotation
summary, annotation correspondence and unit roster match the qualified input
certificate. Zero annotations remain undisposed; the 855 held annotations
remain explicitly untestable. Observed object routes are protein-state/context
14,176, decisive missing assay 13,576, lncRNA RNA product 413, and lncRNA locus
disambiguation 68. Every record retains `no_authoritative_link` with target-link
basis `none`; no credible-set PIP mass or cross-assay eligibility is inferred
from static overlap. Two annotation fields of Catalog r9d (adopted 2026-10-16), not part of r8, add
to these records without changing the link field or any class. A separate
r9d field records fine-mapped protein-altering or canonical-splice
credible-set members per gene (31 genes at normalized PIP ≥ 0.5; H87); it does
not create a variant-to-gene target link for noncoding members, does not
change any class and misses non-canonical splice effects. A
regulatory-region field records whether a bulk-liver accessible region shares
the colocalizing GWAS signal (143 of the 436 primary SuSiE genes; at 108 the region
also colocalizes with the gene's eQTL, and the other 35 are not called linked to
the region; H53). These routes propose experiments and do not establish
mechanism or therapeutic benefit.

The three overlap universes are reconciled in the
[identity tables](../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/identity_reconciliation_v3/summary.json):
14,720 pre-adjudication symbol-grouped rows; 14,704 native bulk/COLOC stable-ID
join rows; 14,682 jointly eligible Catalog identities. Identity adjudication
removes 31 jointly eligible source rows among the 59 quarantines, leaving
14,689 rows; seven genetically untestable records then remain outside the
joint Catalog class universe. The 31 jointly eligible quarantines contain no
positive genetic or disease-state class. One of the other 28 quarantines has
unilateral RNA support, so quarantine is not equivalent to absence of evidence.
Over the bulk fit, 114 tested Ensembl IDs have no Catalog record: 109 are
multi-Ensembl quarantines and 5 are second IDs of a symbol kept under another
ID (two POLR1HASP, two ABCC6P2 and one CES1P1), with no quarantine record. Counted by symbol instead, 103 of the 114 share their symbol with another Ensembl ID tested in the bulk fit (the count in the Figure 6B legend). 35
of the 114 have BH q < 0.05 and one, ENSG00000278732 (Metazoa_SRP, log2FC
+0.535), is a canonical call, so the Catalog holds 1,853 of the 1,854 canonical
calls. Of the 61 Tier-1/2 COLOC IDs without a record, 32 are multi-Ensembl
quarantines, 22 are quarantined as absent from GENCODE v49, six are the genetic IDs of the six
genes described below, and one (ENSG00000223797) is a second ENTPD3-AS1 ID whose
symbol is kept under ENSG00000280739. The 76 program-member IDs without a
record (82 of 6,633 memberships) are outside the bulk fit and every COLOC
table.
The seven nonpositive genetics records are PSMC5, ASAP3, RITA1, ADAM15, CKB,
LILRB3 and GSTM2. The native stable-ID interface keeps them as nonpositive
denominator rows; the Catalog marks their interface class untested.

Direct source-symbol comparison gives 16 categorical-only and zero
native-interface-only symbols, each explained by differing bulk/genetic IDs.
Catalog-versus-native jointly eligible IDs differ by six Catalog-only unique
GENCODE-v49 identity resolutions, all unresolved, and 28 native-only IDs (21
quarantined and seven genetically untestable). The six (CTSLP2, ITFG2-AS1,
MATR3, PAXX, POLR2J4 and SUGT1P3) are keyed on their bulk Ensembl ID, but
their genetics row is the COLOC record of a different ID that is absent from
GENCODE v49; the bulk IDs have no COLOC record. Under a same-ID rule they would
be untested (unresolved 13,121, untested 13,557). The earlier reported 24/8 split
is not reproduced from these keys and is superseded by these explicit source
tables. Semicolon-delimited multi-ID text is not a stable gene identifier.

The 28 historical Myojin class-assay rows still test earlier classes and are
retained as historical, nonconfirmatory records. Rerun with the r8 classes and
the fixed Plan 40 code on the 7,678 screen genes in the four prespecified
classes (84 untested genes excluded; OLS with covariates log1p HLF TPM, minimum
guide count and HLF Chronos), the prespecified omnibus partial F for the three
class indicators is 0.395, Freedman-Lane permutation p = 0.761 (100,000 draws
within 25 guide-count and expression strata). Standardized effects against
unresolved genes are -0.044 for genetically anchored (n = 207), 0.035 for
established-state associated (n = 348) and 0.181 for concordant (n = 11; HC3
95% CI -0.218 to 0.365 LFC). Counting the 84 untested genes in the reference
gives F 0.399, p = 0.760. No contrast reaches the prespecified standardized
effect of 0.20, so the nonconfirmatory verdict is unchanged. The unit is the
gene in a pooled CRISPR screen in one cell line, not a donor. GEN terminal restrictions on
context rescue and powered-negative claims remain in force; its old class
counts do not qualify current r8 classes. Neither source is used to create a
gene-level negative or to change a current Catalog class.

The [current Figure 1–6 source-context companion](../figures/candidates/figure-claim-context-20261014T185701Z_22240302/checks.json)
preserves all 320,281 rows from 147 working source files, with their assay units,
hypothesis families, applicability, unresolved alternatives and experiment
rules. Of 12,119 gene-keyed source rows, 11,727 have a unique r8 gene match;
392 remain explicitly unmatched or ambiguous. All 11,023 program-keyed rows
join to the fixed program identities and membership hashes. These are source
row counts, not numbers of distinct genes or biological replicates. The current
117-program stage context has no selected program in the H59 primary fit
(smallest permutation q=0.323375); the 62-donor equal-variance sensitivity has
one (BICC1/Ductular injury) and zero HC3-selected programs. There are zero
strict negative calls. It remains separate from r8's fourteen historical
multimodal program contexts, which are the only program-level evidence rows in the Catalog:
two programs (Stromal ECM and Ductular injury) in seven datasets from the 2026-08-07 candidate. Liver protein (PXD051911) is untestable for both programs under the Figure 5E scoring rule (at least eight measured members and at least 20% of L1 weight retained): Stromal ECM (IGFBP7) has 8 measured members carrying 16.6% of its weight, and Ductular injury (BICC1) has 2 carrying 9.0% (`figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5E.tsv`). No 117-program stage result or HMSMA spatial row is in
the Catalog. No current gene call, native assay value, program
membership or weight changes through these joins.

Figure 6C worked entries (r8 classes, Catalog r9d values; `figures/main/fig6_gene_catalog/panels/data/panel_deep_20261017/6C.tsv`). GLP1R is untested: no Tier-1/2 SuSiE result (best ABF PP.H4 0.022, MVP NAFLD EAS) and not detected in the bulk fit. HKDC1 is supported in both domains: SuSiE PP.H4 0.988 (UK Biobank AST; > 0.5 in 3 of 5 Tier-1/2 studies with a SuSiE result) and log2FC +1.07 (95% CI 0.68 to 1.47; BH q 8.0e-7). MTARC1 is unresolved: no Tier-1/2 SuSiE result (best ABF PP.H4 0.046, FinnGen NASH) and log2FC −0.15 (−0.22 to −0.08; q 1.2e-4); its fine-mapped missense variant (PIP 0.76 in a PDFF GWAS; PIP ≥ 0.5 in 3 studies) does not enter the genetic call. THRB is genetically anchored: SuSiE PP.H4 > 0.999 (UK Biobank GGT; > 0.5 in 1 of 4 studies) and log2FC −0.20 (−0.26 to −0.14; q 4.2e-10).

### Catalog r9d fields, adopted 2026-10-16

Catalog r9d was adopted on 2026-10-16 with the user's approval. It keeps the 28,233 r8 gene
identities and all r8 class memberships and adds annotation fields from
prespecified analyses. No field below changes an evidence class or a
next-experiment route. Class-dependent counts use r8 classes, COLOC r2 and the
adopted 774-unit bulk fit. All counts in this section are read from the r9d build (seed 20260923):
the r9c build, which corrected the states that r9 wrote for the caQTL link and
single-cell lineage fields, plus the H79 and H81 fields below.

Assay observability (H66). For each of the 28,233 Catalog genes, the Resource
records whether its protein was detected in liver (58 participants) and plasma
(143 participants) by the PXD051911 proteomes, whether it is targeted by the
Govaere CosMx panel, and in which of nine liver cell lineages its RNA is
detected across donors. Liver protein was detected in at least one sample for
7,137 genes and in at least half of the 58 participants for 6,161. Plasma
protein was detected in at least one of 184 samples for 579 genes and in at
least half of the 143 baseline participants for 456. The CosMx panel targets
982 Catalog genes: 942 through a single-gene probe and 40 through one of 19
probes shared between genes. Of 23,593 genes with a single-cell pseudobulk row
(104 donors), 13,029 are detected in at least one of nine lineages (raw count
of at least 5 and CPM of at least 10 in at least half of that lineage's
eligible donors); 4,640
Catalog genes have no single-cell row. Where an assay did not observe a gene,
the entry reads "not observed in [assay]; absence is untestable, not
negative", and no gene receives a `tested_negative` call. Plasma detection
rose little with liver RNA abundance (0.123 between the highest and lowest
abundance deciles, below the prespecified 0.30) and ships with that failure
noted. The lineage field names the lineage with the highest fraction of donors
in which the gene is detected. It reads `not_detected` for 3,643 genes that
are detected in no donor of any lineage and `tied` for 518 genes whose nonzero
maximum is shared by two or more lineages. Ambient RNA makes hepatocyte transcripts appear
in other lineages, so this field measures detection across donors, not
expression level, lineage specificity or cell of origin.

Liver protein coverage by class (H86). Among r8 classed genes, a protein row was
quantified in the PXD051911 liver differential table (6,071 proteins) for 173 of
1,129 established-state genes (15.3%), 203 of 400 genetically anchored genes
(50.8%) and 11 of 26 concordant genes (42.3%); 387 of 1,555 classed genes in
total (24.9%). A gene without a protein row was not quantified and has no
protein evidence in either direction.

Follow-up annotation (H78). A companion table gives gnomAD v4.1 LOEUF and pLI,
the DepMap Public 24Q4 common-essential call and median Chronos gene effect, and
the Human Protein Atlas secretome location and blood concentrations. Coverage
over the 28,233 genes is 17,152 (60.8%) for LOEUF, 16,761 (59.4%) for DepMap and
18,675 (66.1%) for HPA. These fields describe how a gene might behave in a
follow-up experiment. They do not enter any evidence class or route.

Healthy-liver zonation fields (H63). The mRNA field (five non-steatotic
live-donor livers, Visium; Yakubovsky et al. 2026) calls 432 genes central and
177 portal; 3,429 are detected but not zonated, 15,769 have detection too low
to be measured in Visium and 8,426 are not in the assay. The protein field
(hepatocyte protein zonation, Weiss et al. 2026) calls 393 genes central and
426 portal; 897 are measured and not zonated and 26,517 are not in the source
table. Both are healthy-liver reference values owned by their sources. The
test of agreement between the two layers is reported with the spatial results.

Fine-mapped coding variants (H87). 31 Catalog genes carry a fine-mapped
protein-altering or canonical-splice credible-set member at normalized SuSiE
PIP >= 0.5 (20 unresolved, 7 genetically anchored, 3 established-state
associated, 1 untested, 0 concordant); 15 more have one at PIP 0.1 to 0.5 and
50 more below 0.1. This field reports fine-mapped coding variants; it does not
change the eQTL-colocalization class and does not identify a causal gene. The
top variant is missense for PNPLA3 (PIP 0.706), MTARC1 (0.757), TM6SF2 (1.0),
APOE (1.0), SERPINA1 (0.598) and GCKR (0.968). MBOAT7 and TRIB1 have no entry.
The field reads the credible sets of the fine-mapping run `GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/` (release label `noncoding-dna-precoloc-candidate-2026-08-12`). That run has no adoption record of its own; its credible sets enter the paper through two adopted analyses that read the same files, the corrected Figure 2D (installed 2026-10-13) and this Catalog r9d field (adopted 2026-10-16).

Liver caQTL link field (H53). 191 of 581 European tier-1/2 liver GWAS signals
(study leads) with a liver-eQTL colocalization share a signal with a bulk-liver
accessible region (138 donors; at least 45 already reported by Currin 2025); the test is reported with the genetics results. For the 436
genes with a COLOC r2 Tier-1/2 SuSiE signal, `caqtl_gene_link` reads:
three-way supported 108 (the region also colocalizes with the gene's eQTL),
GWAS-eQTL pair only 35 (not called linked to the region), `not_linked` 221 (the
accessible region maps to a different variant), `not_tested_in_source` 60,
`indeterminate` 6 and `colocalized_lead_ld_unconfirmed` 6. Not tested and
indeterminate are untested states, not negatives. The eQTL and caQTL donor sets
may overlap, so three-way support is not independent replication.

Postmortem histologic fat field (H55). The `gtex_postmortem_fat_state` field
covers the 1,853 Catalog genes among the 1,854 disease-state calls (one call has
no Catalog row). At BH q < 0.05 over the 1,002 calls measurable in GTEx liver
(222 postmortem donors), 229 genes rise and 94 fall with histologic fat in the
same direction as between disease and control biopsies, and 3 fall with fat
although they are higher in disease; 676 are not detected, 584 are below the
GTEx expression filter and 267 have no GENCODE v26 match. The slide tissue-area
negative control missed its prediction (rho 0.354, p 0.0564); its stop rule
needed both p < 0.05 and a tissue-area rho at least half the fat rho, and only
the second held. The field concerns histologic steatosis of unrecorded cause in
postmortem liver and does not assign any biopsy call to steatosis.

Sex difference at enzyme signals field (H72). For the 426 genetically anchored
or concordant genes, `enzyme_sex_effect_state` reads heterogeneous 13, one
scale only 13, scale-dependent 4, underpowered 271, untestable 2 and not
applicable 123 (no enzyme colocalization above PP.H4 0.5); no gene is bounded
within a female-to-male ratio of 0.67 to 1.5. The 13 heterogeneous genes sit at
9 distinct SNPs (DDT, GSTT2, GGT1 and MIF share one), so they are not 13
independent observations. The test is reported with the genetics results. The
field describes UK Biobank enzyme GWAS and says nothing about sex differences
in MASLD.

Nearest-gene check (H74). Supplementary only; no Catalog field is added unless
the user adopts one.

Composition sensitivity (H56). No Catalog field was added because the
prespecified stromal positive control failed: the median ratio of
composition-adjusted to unadjusted log2FC for 159 stromal-specific calls was
0.528 in one of eight adjustment settings, where a value below 0.5 was required
in all eight. The 1,854 bulk calls remain described as tissue-level effects.

Ancestry sensitivity (H79). Fields `rna_ancestry_sensitive_h79` and
`rna_ancestry_change_h79`: 64 genes meet the
H79 ancestry-sensitive rule (changed status when four genetic ancestry PCs
enter the model and in at most 9 of 200 permuted-PC fits), against a median of
6 expected false flags. Of these, 40 are adopted calls (20 established-state
associated, 1 concordant, MLIP, and 19 untested; post hoc count against the r8 classes, which r9d keeps) and 24 are not. The label is descriptive and changes no class.

Bulk detection limit (H81). Fields `rna_expression_decile_h81`,
`rna_detection_limit_log2fc_h81` and `rna_below_detection_limit_h81`: for each of the 22,261
tested genes without an adopted call, the detection limit of its expression
decile (log2FC reached with 80% recall in the spike-in). In deciles 1 and 2
(4,323 genes) a missing call does not exclude a 2-fold change, and the entry
reads below the detection limit rather than implying no change. No call or
class changes. Of the 22,261 genes, 113 have no Catalog record (identities
quarantined in r8 as irreducible to one Ensembl gene; 73 of them in deciles 1
and 2) and carry these fields in the sidecar of IDs without a Catalog record.

**Historical September Catalog; superseded by r8 selected for the paper on
2026-10-11.** `RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-20260925T175215Z/`
(corrected bulk primary arm, corrected COLOC release `susie_coloc_r2_20260924T053244Z`,
canonical disease-state rule). The user approved its adoption on 2026-09-25. It
passed PASS-00 to PASS-05 (`passport_gate_status.tsv`); its manual review and
Figure 6 reconstruction were then pending. Those September actions are not
the current manuscript requirements: r8 now supplies the paper's Catalog
counts, and the October 11 decision does not require a portal seal for the
manuscript (see current r8 status below and `docs/STATUS.md`).

| Quantity | Historical September 844-record Catalog | Earlier candidate (2026-08-07 bundle) | State |
|---|---:|---:|---|
| Genes | 28,017 | 30,754 | `portal_export/passport_bundle_status.json` |
| Evidence rows | 56,034 | 61,508 | Same |
| Prespecified programs | 117 | 117 | Consistent with program registry |
| Program memberships | 6,633 | 6,633 | `passport_program_membership.parquet` |
| Contexts | 14 | 14 | Same bundle status |
| `tested_negative` calls | 0 | 0 | Expected under current conservative source/coverage rules |
| Primary class: genetically anchored | 383 | 398 | `passport_gene_index.parquet` |
| Primary class: concordant | 23 | 33 | Same |
| Primary class: established-state associated | 1,054 | 1,214 | Same |
| Primary class: unresolved | 12,700 | 12,998 | Same |
| Primary class: untested | 13,857 | 16,111 | Same |
| Catalog v2 documented fields | 19 | Independently validated candidate requirements |
| Molecular-object follow-up routes | 5 | Deterministic, outcome-blind rules |

**Figure 6 from Gene Catalog r8, 2026-10-12 (historical candidate record; its
28,233-gene counts are selected and local review completed on October 14).**
The selected r8 Catalog is
`RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/`,
with the adopted 774-unit bulk fit and COLOC r2. Its class counts
replace the September Catalog table for the paper; synchronized Resource
adoption remains separate.
`figures/candidates/fig6-catalog-r8-20261012T011202Z/` (`review/verify_counts.txt`,
`source_tables/`). r8 has 28,233 genes and no tested-negative row. Figure 6B
cells, genetics call by RNA call: both supported 26; genetics supported and RNA
indeterminate 400; genetics supported and RNA untestable 10; RNA supported and
genetics indeterminate 1,129; both indeterminate 13,127; genetics indeterminate
and RNA untestable 4,222; RNA supported and genetics untestable 698; RNA
indeterminate and genetics untestable 8,621. Totals: RNA untestable 4,232; genetics untestable 9,319, of which 9,312 have no Tier-1/2 COLOC record and 7 are SuSiE-untestable for insufficient shared posterior (`figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/panels_20261015T173321Z/f6_repair_summary.json`, `fig6b_call_matrix_3x3_source.tsv`). Figure 6D bars: untested 13,551,
unresolved 13,127, established-state associated 1,129, genetically anchored 400,
concordant 26 (both supported; direction not compared). The working 6B and 6D were still on the 2026-08-07 catalog of
30,754 genes (class totals 16,111 / 12,998 / 1,214 / 398 / 33, table above).
Three points a reader needs: (1) 708 genes classed "untested" have one supported
call (698 RNA, 10 genetics), because the other domain is untestable; (2) the RNA
"supported" call requires BH padj < 0.05 and |log2FC| above 0.5, so 11,301
bulk-tested genes with BH padj < 0.05 and a smaller change are "indeterminate"
(counted from `passport_evidence_long.parquet` of r8; 11,335 in the fit itself);
(3) the r8 PP.H4 column shows the best SuSiE PP.H4 for the 436 supported genes
and the best ABF PP.H4 for 15,231 of the 18,478 indeterminate genes, with no
field naming the method. 555 indeterminate rows show a value above 0.5 (152 at
or above 0.8); 222 of these have a finite SuSiE PP.H4 of at most 0.5 that the
ABF value replaces (GNMT shows ABF 0.706; its SuSiE value is 0.128). The
genetics call uses SuSiE only. A candidate per-gene table, not adopted, carries
the best Tier-1/2 SuSiE and ABF PP.H4 in separate columns, each with its GWAS
study and that study's ancestry, and for RNA the log2FC, standard error, t, p
and BH q from the 774-unit fit (138 control and 636 disease genetic-identity
clusters).

The earlier bundle predated the 2026-08-17 COLOC release and used the interval-null
disease-state rule, so the class changes combine review items 1, 2 and 5 with the
canonical rule. The six SuSiE-untestable genetics records carry "untestable /
insufficient shared posterior" and no statistic. Corrected Catalog inputs are in
`Analysis/Multimodal_Program_Projection/candidates/gen-corrected-20260925T175215Z/`.

The noncoding-aware v2 requirements are read-only at
`Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-2026-08-11/`
and passed an independent checksum, field, vocabulary, route, and separation
review in `gene-catalog-v2-contract-review-2026-08-11/`. It contains the public
schema, a human-readable data dictionary, and five follow-up rules, but no gene
outcomes. Regulatory DNA, transcription through a lncRNA locus, and the mature
lncRNA molecule remain separate candidate objects. The first regulatory-DNA
outcomes come from H53 (2026-10-15): 191 GWAS signals share a signal with a
bulk-liver accessible region, and at 108 Catalog genes that region also
colocalizes with the gene's eQTL; these are candidate annotations and change no
class. Genetics-dependent entries
remain on hold until they are rebuilt against the COLOC release adopted on
2026-08-17; the compute dependency is cleared, the rebuild is not done.

Gene Catalog vocabulary:

| State | Meaning |
|---|---|
| `supported` | The prespecified assay supports the evidence claim |
| `discordant` | An adequately observed assay points in a conflicting direction |
| `tested_negative` | The claim was adequately tested and unsupported |
| `untestable` | Coverage or design could not test the claim |
| `not_applicable` | The assay does not address that biological question; not used in the r8 gene calls (0 of 56,466 rows) |
| `source_dependent` | Interpretation depends materially on one source or mapping assumption; the r8 field is true in all 56,466 rows and does not separate genes |
| `indeterminate` | No call under the prespecified rule. For RNA this includes 11,301 genes with BH q < 0.05 whose absolute log2FC is at most 0.5: the direction is resolved but the change is below the effect floor |

## Molecular and physical context

### Validated single-cell ATAC Context v3 candidate (not adopted)

The isolated candidate
`Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/`
uses GSE244832 and GSE281367 only. `NON_GENETIC_READY`, `INTEGRATION_READY`,
and `FIG4F_READY` passed independent validation. Canonical figures, the Gene
Catalog, and release IDs were not changed. The genetics check wrote
`blocked_upstream_release` and no biological genetics rows because at the time
of that run no adopted corrected COLOC manifest existed. Adoption happened on
2026-08-17, so the variant-posterior replay is now runnable and has not been run.

| Candidate quantity | Value | Qualification |
|---|---:|---|
| Consensus peaks | 676,658 | Condition-blind, fixed-width 500-bp, non-overlapping hg38 intervals |
| Hepatocyte / stellate / macrophage peaks | 206,874 / 173,191 / 142,456 | Broad macrophage lineage |
| Cholangiocyte / T-NK peaks | 113,584 / 40,553 | GSE281367 T/NK is a combined source label |
| Peaks supported by both cohorts | 347,509 | 131,663 hepatocyte; 110,048 stellate; 55,300 macrophage; 45,901 cholangiocyte; 4,597 T/NK |
| Sparse count matrices | 16 | Ten consensus matrices plus six old-peak-space sensitivity matrices |
| Direct fragment recounts | 10 of 10 exact | One seed-fixed donor-peak combination per cohort and lineage; counting unit is one deduplicated fragment-file record |
| Official chromVAR rows / tested contrasts | 8,790 / 7,032 | JASPAR2024 CORE vertebrate motifs; eight contrast-testable lineage-by-cohort families |

Matched MASH-versus-normal DA used donors as replicates, separate cohort fits,
TMM normalization, `voomWithQualityWeights`, and BH correction within each
lineage-by-cohort family. GSE244832 MASL donors were excluded. The prespecified
default condition-blind `filterByExpr` rule left no GSE244832 macrophage peaks,
so macrophage DA is `untestable`; the filter was not relaxed after observing
that result.

| Lineage | Donors, GSE244832 normal/MASH | Donors, GSE281367 normal/MASH | Joint peaks | Significant, GSE244832/GSE281367 | Strict replicated | Directional replication: same sign / GSE281367 discoveries; donor-relabeling p; BH q over three lineage rows | Main dynamic criterion |
|---|---:|---:|---:|---:|---:|---:|---|
| Hepatocyte | 5 / 9 | 6 / 6 | 30,650 | 1 / 1,016 | 0 | 849 / 1,016; p = 0.0495; q = 0.074 | Did not pass the sensitivity-stable replicated-peak criterion |
| Stellate | 5 / 9 | 6 / 6 | 9,264 | 820 / 15 | 4 | 15 / 15; p = 0.0340; q = 0.074 | Did not pass: directional q = 0.074 is above 0.05 |
| Macrophage | 4 / 8 | 6 / 6 | 0 | 0 / 0 | 0 | not testable; p = 1 | Untestable after default filtering |

Directional replication is calibrated against all 2,002 relabelings of the 14
GSE244832 donors (9 MASH), refitting the release model each time
(condition-blind `filterByExpr`, TMM, `voomWithQualityWeights`, robust eBayes);
the GSE281367 fits are fixed. The identity relabeling reproduces the saved
GSE244832 log fold changes (maximum absolute difference 5e-15). Exact p is the
fraction of relabelings, identity included, with at least the observed same-sign
count: hepatocyte 99 / 2,002, stellate 68 / 2,002. BH adjustment uses the
prespecified family of three lineage rows in `04_run_da.R`, with untestable
macrophage at p = 1. Over the two testable lineages only, q would be
0.0495 for both; that family was not prespecified and is not used for the
eligibility call. The earlier binomial q values (2.64×10⁻¹¹⁰ and 4.58×10⁻⁵)
treated peaks as independent trials and are withdrawn: under donor relabeling
the hepatocyte same-sign count has mean 505 and SD 211 peaks, against SD 16
under peak independence.

Across the 39,914 jointly testable peaks, evidence states were 4 `supported`,
1,844 `source_dependent` (q < 0.05 in one cohort only; same sign in the other
cohort for 1,505), and 38,066 `indeterminate`, with no `discordant`
peak. All four supported peaks were stellate-lineage peaks, retained direction
under both leave-one-donor-out analyses, and did not reverse under the
GSE244832 batch or broad MASLD sensitivities. Under the same 2,002 GSE244832
donor relabelings, 45 gave at least four supported peaks (p = 0.0225, a single
test added in review; null mean 0.17). The prespecified candidate dynamic-main
eligibility rule also requires
directional q < 0.05; with the donor-relabeling q of 0.074 (three lineage rows),
stellate no longer meets it. The four peaks are candidates; they are not a
main-panel claim and do not prove a causal gene or causal cell type.

The prespecified-program projection contains exactly 234 rows, 117 programs per
cohort. Program scores were testable for 116/117 programs in GSE244832 and
117/117 in GSE281367. MASH contrasts were testable for 78/117 and 117/117,
respectively; GSE244832 cholangiocyte and T/NK contrasts remain explicitly
untestable. Across the paired cohort rows, 156 were `indeterminate` and 78
`untestable`; no program met the strict two-cohort replicated-program rule.
No official chromVAR motif met the strict q<0.05, same-direction replication
rule. chromVAR remains supplementary, and these results do not define a
regulatory network or a tested-negative motif family.

**ATAC statistics correction, 2026-10-22 (candidate tables only).** Corrected
tables are in
`Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-statistics-correction-20261022/`;
the named August candidate remains unchanged. Moderated standard errors replace
79,828 primary peak errors and 80,878 legacy peak-space sensitivity errors.
For the 195 testable cohort-program rows, removing the variance trend on
donor-centered scores changes all 195 standard errors, p-values and BH q-values;
the largest absolute p/q changes are 0.1623094354 / 0.1983822922. BH retains the
78 testable programs in GSE244832 and 117 in GSE281367. No program has q < 0.05
before or after correction. Saved disease effects, peak p/q-values, all 39
untestable program inferences and every evidence state are unchanged. The dependent
program source table carries the corrected p/q-values and an accurate
testable-program BH label; its preceding spatial/protein rows are byte-identical.
`numerical_changes.tsv.gz` records every changed numerical cell, including
duplicated exports. Figures, genetic context and official releases are untouched.
Detailed calculations and checks are in
[the technical correction report](technical/code_logic_performance_followup_20261021/atac_result_update_20261022/REPORT.txt).

The Figure 5F source, after filtering, contains 30 exact cells: two cohorts × five
lineages × promoter-gene coverage, program-score testability, and contrast
testability. Its one-page PDF passed source rederivation, checksum, font,
structural, and visual review. Candidate-only Gene Catalog adapters preserve
30,089 gene-lineage-state promoter rows and 13,266 gene-program-cohort rows;
the optional integration maps at least one ATAC v3 field to 14,510 of 27,187
catalog genes without changing any protected L8 value or `layers_active` vote.

Static peak observability and disease-associated DA are distinct. Distal peaks
are not assigned to genes, and the GSE281367 T/NK context remains
source-dependent. The adopted-COLOC replay was executed and recorded on
2026-08-20: 5,362 signal pairs and 22,356,697 variant-posterior rows reproduce
816 historical colocalized gene-study pairs over 462 genes. The adopted r2
Figure 5D uses one deterministically selected primary signal pair per gene-study pair and
intersects normalized SNP.PP.H4 posterior mass with static peaks from both ATAC
sources across five shared lineages. The installed display counts genes, each at its best state over its study pairs:
22 direct MASLD/liver-fat genes and 424 liver-enzyme genes (from 28 and 739 pairs,
767 total; ten genes in both classes). For each lineage it plots the percent of genes
whose best pair places at least 50% of the primary-signal posterior in peaks of at
least one source or of both sources (gene bootstrap, 10,000 resamples within trait
class, seed 20261017): 0% to 18.2% of direct genes and 3.5% to 13.2% of liver-enzyme
genes for at least one source, and no gene under the uniform-posterior expectation.
Its source table is
`figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5D.tsv`,
computed from
`figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv`; the historical panel had
29 direct and 787 enzyme pairs.
All primary pairs retain at
least 0.95 mapped posterior mass, so none is untestable; most lineage contexts
are partial rather than both-source accessible, and more than half of the 3,409 partial pair-lineage rows (1,943) place under 1% of the posterior in peaks of either source (239 of 767 pairs, 13 of 28 direct, are under 1% in every lineage). Against a locus-matched expectation (the fraction of each pair's tested variants in peaks), liver-enzyme genes carry 2.5 to 3.1 times the expected mass on average in every lineage (n=424 genes, gene-bootstrap intervals above 1), but in each lineage at most half of them (94 to 212 of 424) exceed their own expectation; for the 22 direct-trait genes the intervals include 1 (Figure 5D locus-matched background companion; 2,000 gene resamples, seed 20260923; table
`figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_locus_background_summary.tsv`). Peak overlap is regulatory-DNA
context, not target-gene assignment or a functional test. Pair-level PP.H4
reproduced the adopted release within the documented tolerance, but no
adopted per-variant posterior vector exists for an external equality check.

### Validated spatial candidate (not adopted)

**V4, October 8 candidate; qualifications checked 2026-10-13.** In
`Analysis/Hypotheses/V4-20261008T015029Z/summary_signflip_exact.tsv`, HMSMA has
32 analysis units, mean centered residual Moran's I 0.0982723, all positive.
The prespecified one-sided exact sign-flip nominal p is 2.3283e-10; Holm
adjusted p is 4.6566e-10 over the two-test HMSMA/Vu family. [Inference] The
description of these units as 28 disease participants and four controls rests
on inferred exclusion of individuals 27–29 as controls, without an
authoritative person-to-clinical-label join (`PRESPEC_V4.md`, roster section).
This result is spatial consistency under proxy abundance adjustment; it does
not establish a disease/control or stage contrast. Its adjustment is weaker than
Figure 5F's: HMSMA scores are residualized for a 15-marker hepatocyte proxy, and
the GSE192741 and Vu rows for the cell2location hepatocyte abundance only, not all
16 abundances. Its 0.0983 therefore cannot be compared with the 5F excess (0.040
and 0.050) and does not replicate it. The Vu units are arrays holding several
biopsies each, not donors. It is distinct from the
older matched-gene calibration and array-level summaries below and is not
adopted by this documentation correction.

The read-only candidate `spatial-resource-candidate-2026-08-11` passed 32
analysis/requirement checks and the six-panel Figure 5 source/PDF checksum check.
It is suitable for adjudication and downstream candidate builds, but it is not
the current official Resource release.

| Candidate quantity | Value | Interpretation |
|---|---:|---|
| Registered dataset or analysis decisions | 14 | Every source is retained, supplemental, skipped, or dropped explicitly |
| Assays requiring complete coverage | 7 | GSE192741, Vu, HMSMA, Yakubovsky, CosMx, GeoMx, and GSE287826 |
| Prespecified-program coverage rows | 819 | Exactly 117 programs × 7 assays; no outcome statistics or winner selection |
| Gene-context rows | 23,576 | Assay-facing measurement/testability records |
| Confirmatory effect rows | 20 | Exactly 2 predeclared programs × 10 registered source contexts, including untestable and not-applicable rows |
| HMSMA technical spots / arrays | 130,097 / 35 | Technical coverage only; neither is a biological donor count |
| Candidate Figure 5 output files | 6 PDFs / 6 source tables | All PDFs are one page and checksum-linked to their exact source table |
| Downstream spatial figure output files | 4 PDFs / 4 source tables | Figure 1 capability, Figure 6 Gene Catalog, all-program coverage, and HMSMA descriptive panels validated in candidate `spatial-figure-integration-candidate-2026-08-11-r2` |

Coverage status across the seven required assays:

| Assay | Observable | Partial | Untestable | Biological-unit qualification |
|---|---:|---:|---:|---|
| GSE192741 Visium | 114 | 3 | 0 | 5 sections collapsed to 4 donors for inference |
| Vu et al. Visium | 114 | 3 | 0 | 10 arrays; donor mapping unresolved, so results remain source-dependent |
| HMSMA Visium | 114 | 3 | 0 | 35 arrays; donor/clinical mapping unresolved |
| Yakubovsky normal-liver spatial | 114 | 3 | 0 | Normal-liver reference, not MASLD validation |
| Govaere CosMx | 12 | 97 | 8 | Targeted-panel observability; no whole-program inference when coverage fails |
| Govaere GeoMx | 0 | 0 | 117 | Dropped from inference because direction/normalization remain unresolved |
| GSE287826 | 0 | 0 | 117 | No authoritative 30-AOI-to-15-donor key |

The confirmatory family was recorded before execution and contains only the Stromal
ECM and Ductular injury programs, fixed in advance from the superseded August
stage fit; neither passes the corrected 117-program rule. The prespecified
confirmatory design residualized spot scores for the hepatocyte abundance only;
the all-16 model plotted in Figure 5F was added after that analysis as a
supplementary falsification of composition dependence. The values below are
those plotted in Figure 5F: spot scores residualized for all 16 cell2location q05
abundances plus log library size and log detected genes, expressed as excess
over the matched-gene null mean, with 9,999 matched-gene draws (seed 42 for
GSE192741, 1042 for Vu) and BH over the two programs within each source. In
GSE192741, Stromal ECM had donor-first excess 0.04050 (matched-null SD 0.01098;
4 of 9,999 draws at or above it, P=0.0005, q=0.001) and is `supported`.
Ductular injury had excess 0.01532 (SD 0.01060, P=q=0.0793) and is
`indeterminate`, not negative.

In Vu, Stromal ECM had excess 0.05034 (SD 0.01033; 2 of 9,999 draws, P=0.0003,
q=0.0006). Its evidence state is `source_dependent` because the ten arrays cannot
be joined to biological donors. The r8 Catalog program-context row for Vu
Stromal ECM reads `supported` and should read `source_dependent`. Ductular
injury had excess 0.002654 (SD 0.008228, P=q=0.3158) and is `indeterminate`.

In the prespecified model, with the cell2location hepatocyte abundance as the only
composition covariate (the discovery lineage of both programs), the excesses were larger: Stromal ECM
0.1390 in GSE192741 and 0.1471 in Vu, each with no null draw at or above it
(P ≤ 1e-4, the resolution floor of 9,999 draws; q ≤ 0.0002), and Ductular injury
0.02587 (P=q=0.0755) and 0.01922 (P=q=0.1158). The two models give the same
evidence states. The four r8 Catalog spatial program-context rows carry these
hepatocyte-only estimates, not the plotted all-cell values.

The matched-gene null resamples genes; donors and arrays are held fixed. Under
the plotted all-cell model, residual Moran's I was positive for Stromal ECM in
3/4 GSE192741 donors and 10/10 Vu arrays, and for Ductular injury in 4/4 donors
and 9/10 arrays, so neither positive sign nor cross-unit consistency
distinguishes the programs. Stromal ECM excess over each donor's own
matched-gene null was positive in 3 of 4 donors (0.0156, 0.0657, 0.1023 and
-0.0216 for H35 to H38; exact one-sided sign-flip P=0.19, descriptive; the floor
with four donors is 0.0625). Stromal ECM exceeded its gene-level matched null in both sources and Ductular injury did not; neither source gives a donor-level test (4 donors in GSE192741, arrays without a donor key in Vu); no prespecified test compared the two
programs. Descriptively, the per-unit Stromal ECM minus Ductular injury excess
was positive in 3 of 4 GSE192741 donors and 10 of 10 Vu arrays (arrays are not
donors). **Generic spatial autocorrelation can mimic program transportability,
and an assay-matched null is required to separate the two.** This is not a
negative Ductular injury result, and the Vu contrast remains source-dependent
rather than donor-level replication.

Canonical Figure 5F pairs this calibration with a registered H&E overlay for
Guilliams et al. section JBO019 (steatotic donor H37), selected outcome-blind
because its 1,164 graph-eligible spots equal the dataset median. The displayed
Stromal ECM score uses 26 measured fixed members, retains 0.5475 of the original
L1 weight, and is residualized for all 16 cell2location abundance factors plus
library size and detected genes. Re-derivation from the spot matrix reproduced
the recorded dataset-level residual Moran's I exactly (0.0582734890); this is the
four-donor collapsed value, not the overlay section's. Section JBO019 itself has
residual Moran's I 0.1213 and excess 0.1023 over its matched-gene null, the
largest of the five GSE192741
sections (excess -0.0151, 0.0462, 0.0657, 0.1023 and -0.0216 for JBO014, JBO015,
JBO018, JBO019 and JBO022): it was chosen by spot count, not by score, but it
shows the strongest section. The overlay adds no test:
it is a representative physical view without pathologist region masks and does
not authorize localization to a named histologic structure.

HMSMA is deliberately descriptive: across 35 arrays, median per-array Moran's
I was 0.1619 for Stromal ECM (range 0.03844–0.4849) and 0.01237 for Ductular injury (range
-0.02244–0.08640). No population P value, disease contrast, fibrosis/NAS model,
histologic localization, or MALDI colocalization is authorized without the
authoritative joins.

The candidate portal exports the registry, gene context, 117-program coverage,
and two-program effects as release-linked products. Legacy numeric
`s5_spatial` fields remain compatibility-only and are excluded from rankings,
active-layer counts, evidence-strength bars, and default decisions. Candidate
gene and program pages consume the categorical tables and expose testability,
source dependence, and biological versus technical units; the files remain
unstaged from the public data directory until synchronized adoption.

| Layer | Current conclusion | Qualification |
|---|---|---|
| Tissue protein | Useful orthogonal transportability evidence | Must remain assay-native and abundance-aware |
| Secreted protein | Supplementary independent context | Cannot be collapsed with tissue proteome into one score |
| Chromatin/regulatory context | Supports cell-context observability | Does not convert a colocalized candidate into a causal gene |
| Single-cell localization | Donor-aware cell-labelled context | Cells and runs are never biological replicates; ambient-sensitive programs do not establish lineage origin |
| Spatial organization | Physical-context support in accepted cohorts | Existing cohorts do not provide a complete MASLD control/donor/histology/metabolite design |
| Yakubovsky | Zonation and local lipid-context reference | Normal liver reference, not MASLD disease validation |
| GSE287826 | Could not be tested inferentially | No authoritative public donor key |

The current canonical Figure 5 contains all six adopted panels and the vector
assembly `fig5_molecular_context.pdf`. Panel 5D carries the variant-consistent
ATAC context; 5E carries assay-native program context; and 5F carries the
all-cell-adjusted matched-gene calibration plus one registered Stromal ECM H&E
overlay. The calibration carries the spatial inference; the map only shows the
score's physical distribution in an outcome-blind representative section.

Figure 5E (installed table
`figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5E.tsv`)
covers all 117 programs. Thirty-one contain at least one lineage-open colocalized
member gene; for three of them the gene is RORA, a direct-trait gene. Protein
scores are testable for 95 programs (n=46; at least eight measured member proteins
and 20% retained weight), and Freedman-Lane permutation P (9,999 permutations) is
BH-corrected over those 95. Twenty-two programs have a robust protein-histology
association (14 negative, 8 positive); 8 of the 22 are among the 31 lineage-open
programs. The panel shows 45 programs: the 31 lineage-open programs and the 14
robust programs with no lineage-open gene.
**Figure 5B/5C corrected on 774 units, 2026-10-14.** Installed panel 5B
contains 3,054 measured genes from the fixed historical July 10,044-gene
universe, with 211 protein-significant genes (188 concordant, 23 discordant)
and 2,843 not protein-significant. Fourteen protein-measured universe genes
lack an accepted bulk row. Spearman rho=0.457835 and direction agreement
65.9136%, rounded to 66% on the panel (rho 0.346 across all 5,836 measured proteins, 0.487 in the 2,390 transcriptomic-arm genes and 0.305 in the 664 genetic-only genes; the transcriptomic arm was selected from RNA-based analyses that include the retired 846-sample bulk contrasts). Among protein-significant genes, sign agreement is 170/177 when RNA BH padj <0.05 and 18/34 otherwise. OLS slope=0.455025,
intercept=−0.023131 and R²=0.183900. The earlier 211/185/26 and 64% working
panel and September 844-record 210/183/27 and 65% estimates are superseded.
Panel 5C retains the fixed 25 proteins, with nine canonical RNA calls and
23 protein calls. CD74 enters and HLA-DQA1 leaves the RNA call set relative
to September. HLA-DQA1 uses primary ENSG00000196735.14 (RNA log2FC=0.379801,
BH padj=0.024969), fails |log2FC| >0.5, and has an open RNA marker. Its
protein log2FC=1.393, padj=0.07964, so that marker is also open. The historical
filled RNA call used ENSG00000206305.12. Selection and estimation share the
protein cohort; the fixed rows do not provide independent validation.
Job 22236311 completed, focused source and panel checks passed, and prior
PDFs are preserved. Source: [774-unit molecular context](../figures/candidates/bulk-context-774-20261014T140904Z_22236311/figure5/).
Coordinate-only maps of both programs are retired from the active panel tree.
Candidate `fig5-completion-2026-08-24-v6` passed 80 panel/source/assembly checks
in job `21068950`; post-copy multimodal, Ghostscript, and Resource-scope checks
passed in job `21068987`.

The release-linked downstream candidate adds four one-page panels without
changing any spatial call: Figure 1 separates biological from technical units
and displays source/join/coverage checks; Figure 6 contrasts the supported
Stromal ECM program, indeterminate Ductular injury program, and untestable HMSMA clinical
branch; the supplement reports all 819 program-by-assay coverage rows and all
70 HMSMA array-by-program descriptive rows. Each panel has an exact source
table and checksum manifest. The candidate is validated and visually reviewed,
but remains unpromoted.

The visually accepted publication-layout candidate is
`spatial-publication-figures-candidate-2026-08-11-r4`. It contains 13 individual
one-page PDFs: 11 were simplified and two already-clear panels, Figure 5E and
the all-117-program coverage heatmap—were retained byte-for-byte. Across the 11
redesigned panels, plot text fell to 37.7% of the prior total. This is a visual
change only: exact estimates, evidence states, unit definitions, and complete
families rederive from the recorded source tables. The focused cell-context panel
shows the Stromal ECM–fibroblast and Ductular injury–cholangiocyte comparisons needed for the
paper message; the complete 64-row factor summary remains in the validated
cell-context candidate.

The result-led successor is
`spatial-impact-figures-candidate-2026-08-11`. It retains 12 r4 panels
byte-for-byte, replaces the main 5A and 5F with the compact calibration logic
and matched-null discrimination panels, and carries the former source-state
matrix into the supplement, for 14 one-page PDFs total. Its four-row
transportability profile independently rederives the prespecified program hashes,
9,999-draw nulls, unit-sign counts, accepted evidence states, and all-cell
composition sensitivity. It adds no inferential family or spatial score and
writes no canonical output. This candidate passed source, checksum, page,
font, claim-boundary, and visual review and is the preferred spatial assembly
family pending scientific adjudication and synchronized adoption.

### Validated spatial diagnostics (not adopted)

Four bounded candidate lanes now test the spatial result without expanding the
two-program confirmatory family or changing Figure 5 evidence states:

| Diagnostic candidate | Validated content | Paper role |
|---|---|---|
| `spatial-full-composition-sensitivity-candidate-2026-08-11-r2` | 82 checks; all 16 shared cell2location q05 factors plus the same two QC covariates; four dataset-by-program rows; full-rank 19-column designs | Supplementary falsification of composition dependence |
| `spatial-unit-null-diagnostics-candidate-2026-08-11-r2` | Three PDFs and three exact source tables covering physical-unit heterogeneity, matched-null calibration, and weight sensitivity | Makes donor collapse, source dependence, and indeterminacy auditable |
| `spatial-cell-context-attribution-candidate-2026-08-11-r2` | 480 within-section/array correlations, 448 donor/array reporting rows, 64 complete factor summaries, and one full-family heatmap | Descriptive context attribution only; no spot-level P values |
| `hmsma-mask-sensitivity-candidate-2026-08-11-r2` | 210 rows = 35 arrays × 2 programs × 3 UMI thresholds; all 35 input H5ADs checksum-recorded | Tests stricter proxy masks while donor and clinical joins remain unresolved |

`spatial-diagnostic-bundle-candidate-2026-08-11` is the single validated index
for these products: five accepted components, six one-page PDFs, 58 candidate
output files, and 27 executable code dependencies are checksum-linked, with
canonical writes disabled. Superseded visual candidates and failed validators
are not registered in the bundle.

With all 16 abundance factors included, the Stromal ECM centered spatial effect
fell from 0.1390 to 0.04050 in GSE192741 (`q=0.001`) and from 0.1471 to
0.05034 in Vu (`q=0.0006`). It therefore remains supported within GSE192741
and within-source robust in Vu, but most of the original magnitude is
composition-linked and Vu remains source-dependent. Ductular injury remained
indeterminate: 0.01532 (`q=0.0793`) in GSE192741 and 0.002654 (`q=0.3158`)
in Vu after full adjustment. The accepted supplementary attenuation panel is
`spatial-full-composition-figure-candidate-2026-08-11-r4`; its source table
retains all four rows and the noncausal claim boundary.

The full 16-factor descriptive heatmap localizes the Stromal ECM program most
strongly to fibroblast-rich context: median within-unit Spearman rho was 0.522
across four GSE192741 donors and 0.398 across ten Vu arrays, with all reporting
units positive. Ductular injury–cholangiocyte covariation was weaker (rho 0.101, 3/4
GSE192741 donors positive; rho 0.174, 10/10 Vu arrays positive). These are not
independent validations because program scores and cell2location abundances
come from the same Visium RNA matrix.

The physical-unit panel shows why aggregation matters: H35's two Stromal ECM
sections had residual Moran's I of 0.080 and 0.328 with hepatocyte abundance as
the only composition covariate (0.003 and 0.077 under the all-cell adjustment of
Figure 5F), but they
contribute one donor-collapsed value. No section- or array-level P value is
reported. In HMSMA, Stromal ECM remained positive in all 35 arrays under 1,000- and
2,000-UMI proxy masks (rank rho 0.977 and 0.937 against the 500-UMI baseline;
sign agreement 1.00). Ductular injury was less stable at 2,000 UMIs (rho 0.750; sign
agreement 0.80), reinforcing its indeterminate status. This sensitivity cannot
replace an authoritative tissue mask because the available H5ADs were already
filtered at 500 UMIs.

The CosMx and Yakubovsky audit supports no expansion: CosMx retains only its
source-native, array-level IL32–macrophage vignette, and Yakubovsky remains a
three-donor normal-liver boundary with both programs indeterminate. No new
panel or inferential claim is warranted for either source.

**Gene-level healthy-liver zonation annotation (H63, 2026-10-15)**
([RESULT](../Analysis/Hypotheses/H63-20261015T125909Z/run-20261015T131038Z/RESULT_H63.md)).
Verdict: Pass; no stop rule fired. Over 1,337 non-landmark Gene Catalog r8 genes
measured in both layers, mRNA position along the central-portal axis in five
non-steatotic live-donor Yakubovsky livers (M4 to M8) agreed in direction with
Weiss et al. 2026 hepatocyte protein zonation: Spearman rho 0.3914 (95%
gene-bootstrap interval 0.3429 to 0.4392), above the 0.30 floor, and positive in
each of the five donors (0.288 to 0.409). The gene-permutation p (9.999e-05, the
floor) is too small because genes covary; the donor condition and the 0.30 floor
govern the reading. Agreement is higher in hepatocyte genes (rho 0.586, 562
genes) than in other genes (0.267, 775 genes), as the prespecification expected
because Visium spots also hold portal-tract and central-vein non-parenchymal
cells while Weiss et al. measured single hepatocytes; this split is descriptive. Among T1 genes called in both layers, 20 of
165 disagree, 17 of them mRNA central with protein portal. CPS1 is not zonated
at the mRNA level but portal at the protein level; this disagreement was known
before the plan was written. Weiss et al. reported an
RNA-protein zonation comparison against single-nucleus RNA (their Supplementary
Fig. 9), so this is a check for the Catalog, not a new finding. Healthy liver
only; it says nothing about zonation in MASLD. H63 uses the five non-steatotic
donors M4 to M8; the three-donor description above belongs to the earlier
two-program lipid-zone analysis. The candidate columns
(`h63_zonation_columns.parquet`) are not in `catalog_reviewed_r8`.

## Model framework and transfer limits (Figure 7): validated candidate, working figure

Release state: **Validated candidate** for every number in the first table (each re-derived from
its producing file on 2026-09-07); **development evidence, not adopted** for the second table.
Figure 7 asks what current predictive models learn from the Resource maps and how one would know
it is not prior exposure. Every score is reported beside a per-cohort audit of the model's
training exposure and beside its task-native baseline. Registry: `Analysis/MASLD_Model_Benchmark/MODELS.md`.

| Quantity | Value | Unit and source |
|---|---:|---|
| Released bulk-RNA severity instrument, out-of-cohort Spearman with Kleiner stage | **0.6631** (tie ceiling 0.9523; in-cohort ridge 0.6094) | 109 participants, GSE268273; `release/masld-severity-v1.2/expected_ood.json` (bitwise = v1.1). GSE268273 also chose the shipped transform: the highest Spearman of four coherent, deployable arms (0.6506 to 0.6631), under the prespec 2c rule that `release-instrument-20260901T211017Z/prespec/AMENDMENT_03.md:53-54` kept; the amendment was written before the quantile arms were run and after the first three transform families' results were read |
| Same instrument, two additional external cohorts, Spearman / tie ceiling in the Figure 7D units | **0.709 / 0.866 / 0.680** (record-level 0.677 / 0.774 / 0.794) | GSE213621 263 genetic identities, GSE276114 81 MASLD and 82 chronic viral hepatitis participants (Spearman 0.6680, 0.7860, 0.5443); record-level: GSE213621 299, GSE213621 4-level 367, GSE276114 177 participant records; `executions/external-validation-expansion-20260901T190000Z/results.json`. No demonstrated advantage over the recorded transductive in-cohort ridge; intervals including zero do not establish equivalence. In GSE276114, 9 or 10 of the 81 MASLD and 24 of the 82 chronic viral hepatitis participants have HCC (peritumoral tissue at least 1 cm from the tumor; source-study counts, no per-sample label; `MODEL_RESULTS_SUMMARY.md:3135-3176`) |
| Historical repeat-biopsy ICC diagnostic, out-of-contract scoring | ICC **0.7658** vs **0.1886** | 58 source-declared GSE193066 pairs, identity-unqualified; no primary within-person reliability, ICC or label-error calibration claim. Updated `release/masld-severity-v1.2/MODEL_CARD.md` §4.2 reports only 3 pairs with both biopsies inside shipped acceptance bounds. Not accepted-input deployment performance; source-reported, not independently recalculated in this review |
| Activity head (v1.2), Spearman with NASH-CRN activity sum, four cohorts never trained on | 0.4367 / 0.6054 / 0.4900 (+ one cohort under source terms) | 106 / 93 / 91 participants; sealed gate PASS at the minimum 2 of 4; no two-axis claim; `MODEL_CARD.md` §4.4 |
| Pretraining exposure of the external validation cohort | **109 of 109** GSE268273 samples are BulkFormer pretraining rows (one-to-one, constant offset); 1,141 of 1,142 project samples in its training split | `executions/rerun-floor-20260905T172141Z/EXPOSURE_MEASURED.md:28-34`, `SESSION_FINDINGS.md:631`. An "unseen-sample" claim is unavailable for that model |
| Pretrained bulk model under realistic gene absence | Spearman **0.6631 → 0.4318** with 2,472 genes removed; matched ridge does not fall | `SESSION_FINDINGS.md:720-725` |
| Generic nonlinearity on the severity task (TabICLv2 vs ridge on the same 30 or 100 training-fold PCs) | Prespecified GSE268273 test (109 participants): four comparisons, **−0.014 to +0.001**. Eight comparisons against ridge on raw stage (two in GSE268273, six post hoc in GSE213621 under 3-level and 4-level coding and in GSE276114): mean **+0.0031**, range −0.014 to +0.027. All below MDE 0.0667 | GSE268273 109, GSE213621 299 and 367 records (two codings of one cohort), GSE276114 177 participants; `executions/tabicl-nonlinear-comparator-20260905T000000Z/results/results.json` (deltas) and `secondary_deltas.json` (post hoc) |
| Bulk RNA ↔ H3K27ac coupling, log2 CPM squared-error skill, rank-4 reduced-rank regression | **0.1221** [0.0827, 0.1592] forward; **0.1348** [0.1024, 0.1653] reverse | 99 participants, GSE267145, five pre-registered participant folds, participant bootstrap; `executions/chromatin-state-thesis-20260906T124724Z/out/a2_results.json`. Rank 4 is a floor: skill rises to 0.1505 / 0.1626 at rank 10; nested rank 11–12 beats rank 4 by +0.0206 [0.0136, 0.0276] |
| RNA + H3K27ac late fusion vs the better single assay, three-state macro-F1 | **+0.0030** [−0.0754, +0.0806] | 99 participants; `chromatin-state-thesis…/RESULTS.md:64`. Mechanism: residual correlation of the two assays' errors 0.945 (`MODEL_RESULTS_SUMMARY.md:4993`); the interval does not exclude a five-point gain |
| DNA language-model heads, signed variant-effect direction | **0.5016 / 0.5015** vs allele-feature ridge 0.5033 vs measured chance 0.4998 | `MODEL_RESULTS_SUMMARY.md:86`. No positive margin over the allele baseline is established |
| Joint-space refusal bound for the severity instrument | measured, **not shipped**: admits GSE281797 (0 of 94), refuses GSE225740 (73 of 93) for coverage 0.605 | `executions/severity-v12-joint-contract-20260907T162640Z/RESULTS.md`; shipped as a label-free diagnostic only |

The 2026-09-08 [focused rigor assessment](technical/OAKES_PACHTER_RIGOR_INTEGRATION_ASSESSMENT_2026-09-08.md)
independently reproduced the four severity-transfer correlations, not every
comparison in this table. Raw Spearman values in the additional three contrasts
are 0.638198, 0.747294 and 0.673865; the displayed 0.677/0.774/0.794 values are
ratios to tie ceilings. The native comparator uses full-evaluation-cohort
expression scaling, and its saved aggregate paired intervals were not
independently reproduced from native participant predictions. Repeat-biopsy ICC
does not establish accuracy or absolute-stage calibration and remains a
historical diagnostic on unqualified GSE193066 identities, not primary
within-person reliability. Other model rows
retain their prior evidence state, not a new independent reproduction claim.

A **candidate rebuild of panel 7B alone** (2026-09-21, not promoted, not adopted) positions
each gain against **its own minimum detectable effect** rather than on a shared raw axis, so
that unlike outputs can be compared without being collapsed; the native gain, interval and
MDE are printed beside every row. MDE is derived per row as **2.802 × SE**, the 80%-power
convention this repository already uses, with SE recovered from that row's own 95% interval.
It is a power bar, not a significance bar, which matters because three rows have intervals
excluding zero while sitting below their own MDE. Read this way the panel separates
increments that clear their design's floor — RRR rank 4 transferred **+2.80×**, graded
fibrosis **+1.28×**, MPRA hyenadna **+1.03×** — from those that do not, and shows the two
endogenous cis-eQTL embedding rows at **−1.34×** and **−1.59×**, measurably worse than a
distance ridge. The rebuild also adds the campaign's input-length contrast as eight per-fold
rows: the same adapted head clears its MDE in **4 of 4** folds against native ATAC at a
matched 2-kb input (**+1.7× to +4.2×**) and in **0 of 4** at the native 1-Mb input. Those
eight rows are development evidence from the regulatory campaign and are labelled as such.
Two rows with no definable MDE were retired from the gain axis and are listed with the panel
that still carries their content. Source:
`Analysis/MASLD_Model_Benchmark/executions/fig7_panels_rebuild-20260921T191140Z/`.

A **five-panel candidate redesign of Figure 7** (2026-09-22, not promoted, not adopted;
`Analysis/MASLD_Model_Benchmark/executions/fig7_redesign-20260922T143859Z/`) replaces the
six method-name forest lists with panels that plot observations on real axes, in the style
of Figures 2-4. 7A plots one point per participant: the released severity score by fibrosis
stage in the three cohorts whose records it was not fitted on (109 / 299 / 177 records; GSE268273 also chose the shipped transform; Spearman
**0.6631 / 0.6382 / 0.6739**, tie ceilings 0.952 / 0.942 / 0.848, each re-derived from the
per-participant scores through the released path and asserted equal to the recorded value).
7B is the two training-set-size curves with their per-seed points. 7C is the masking
contrast: pretrained BulkFormer **0.6631 -> 0.4318**, its random-initialisation twin
**0.5193 -> 0.2881**, ridge on the same genes **0.6445 -> 0.6564**, read from
`bulkformer-feasibility-20260905T000000Z/results/mask_ablation{,_ridge}.json` rather than
from the session summary. 7D shows each sequence task in its own units with the task-native
baseline as a point and dashed line (MPRA: allele-identity ridge 0.1342, HyenaDNA 0.2162;
endogenous |z|: distance ridge 0.2344, HyenaDNA 0.2126; signed direction: chance 0.4998,
allele ridge 0.5033, best head 0.5016) and, beside it, the adapted head's per-fold MSE
improvement over native ATAC at 2 kb and 1 Mb (development evidence). 7E is a tile matrix of
per-study gain over the liver-specific encoder for six encoders and two heads, with 22 of 84
cells grey (exposure unresolved, not drawn) and 18 marked as confounded. No number changed;
the framework schematic, bounds ladder and verdict list leave the main figure for caption
and supplementary material. `figures/main/fig7_model_framework/` is unchanged.

Development evidence, GSE267145 (99 participants, one sample each; not adopted, not an
external validation; `Analysis/MASLD_Model_Benchmark/executions/chromatin-beyond-histology-20260907T161521Z/RESULTS.md`):

| Quantity | Value |
|---|---:|
| Activity head fitted on four other cohorts (weights fixed before this analysis, no chromatin) vs H3K27ac steatosis-selective component, partial Spearman given steatosis, ballooning, lobular inflammation, fibrosis and sex | **0.434** [0.181, 0.629]; Freedman–Lane 0/5,000; MDE(80 %) 0.293 |
| Same, additionally given 14 RNA and 7 H3K27ac library descriptors and BayesPrism hepatocyte fraction | 0.426 [0.163, 0.672] |
| Nested-CV: component re-selected in every training fold (index 3 after concentration residualisation), loading stability mean \|cos\| | **0.711** (fixed index 5: 0.353) |
| Bulk-RNA ridge (42,163 genes, λ by inner CV) predicts that component out of fold, partial given the grades, sex and RNA descriptors | **0.552** [0.351, 0.704]; regenerated 500-draw null q97.5 0.229, 0/500; MDE 0.343 |
| Within steatosis grade 1 (n = 43), top-8 vs bottom-8 by RNA prediction, measured chromatin component | 1.36 SD (random 8-vs-8 q97.5 0.74; 0/2,000) |

Boundaries. Cross-sectional ordering only; a within-batch ranker whose absolute level is not
recoverable; no diagnostic, prognostic, treatment-response or clinical claim. The shared
within-histology RNA–chromatin axis is one cohort, discovery-motivated, unadjusted for age and BMI
(absent from the deposit), and has no second H3K27ac cohort; what the axis is remains open. It
does not change the NAS-component separability state, now `indeterminate`: its per-region
association statistic takes the same values under one shared program and under three programs
that each follow their recorded grade. That statistic concerns per-region association patterns
rather than a donor-level factor conditional on the grades.
Derivative weights trained on GSE267145 or GSE296875 are releasable under the
[adopted 2026-09-07 owner decision](decisions/2026-09-07-derivative-weight-release.md).
Raw and processed source data, including per-participant values, are not rehosted;
GSE126848, used as an external RNA-transfer cohort in the benchmark, is one of the five
Resource bulk cohorts and is not contact-naive at the paper level.

Training-set size, **validated candidate as of 2026-09-21** (Figure 7 supplement;
`Analysis/MASLD_Model_Benchmark/executions/learning-curve-20260913T143818Z/RESULTS.md`).
Its state was raised from development evidence by an independent reconstruction that
recomputed every curve statistic from the deposited per-seed points through a second
code path: **10 quantities agree to 1.1e-16**, the sealed prespecification's sha256
matches both results files, and all four written predictions re-score to the same
outcome. No model was refit and no fold rescored, so this reviews the curve statistics,
not the fits. Source: `executions/learning-curve-review-20260921T221859Z/out/checks.json`.
Each shipped head was refit on nested random subsets of its own training set, five seeds per
size, and scored on the held-out data of the shipped number. The full-data refit reproduces
both shipped values (absolute difference 0 and ≤3.3e-16).

| Quantity | Value |
|---|---:|
| Severity instrument, out-of-cohort Spearman (109 GSE268273 participants), by training **samples** (pool has no donor key) | 130: 0.6357 · 208: 0.6479 · 287: 0.6498 · 364: 0.6589 · 442: 0.6584 · 521: **0.6631** |
| RNA→H3K27ac cis head, out-of-fold skill (99 GSE267145 participants), by training participants per fold | 20: 0.0320 · 31: 0.0797 · 43: 0.0973 · 55: 0.1187 · 66: 0.1312 · 78: **0.1426** |
| Same subsets, cis head minus the same recipe without cis genes | 20: 0.0200 · 43: 0.0661 · 78: **0.0946** |

Quadrupling the severity pool adds 0.027. The severity curve is **flattening but not
flat**: its first increment is **+0.0122** and its last **+0.0047**, both reproduced in
the reconstruction. The cis head has not plateaued at 78 participants.
The cis-head recipe uses min(81, training participants − 1) expression components, and no
subset reaches 81, so participant count and component count rise together along that curve.
Only the difference row compares arms with the same component count; the reconstruction
confirmed independently that npc equals training participants minus one at every point and
that the 81-component cap is never reached, so the **nested-minus-globals** series
(**+0.0200 → +0.0946**) is the component-matched statistic and the raw nested series is not.

### Ceilings, matched baselines and the fold-0 adapter confirmation, 2026-09-22 (development evidence, not adopted)

Five reference measurements were made so that each Figure 7 lane is read against what its own
data allow, and the closed fold 0 of the Currin caQTL lead set was then spent on a single
prespecified confirmation. Nothing below changes an adopted number, weight, threshold, gene
class, hero gene or figure. Sources:
`Analysis/MASLD_Model_Benchmark/executions/fig7-step0-ceilings-20260922T175014Z/RESULTS.md` and
`GWAS/finemapping/results/alphagenome_campaign/fold0-confirmation-20260922T191005Z/RESULTS.md`.

**Fold-0 confirmation.** The selected `adapter_r16_last5_symmetric` five-seed ensemble, trained
on folds 1-4 and scored on fold 0 (6,574 variants, 524 blocks), beats AlphaGenome native ATAC
at a matched 2,048-bp input by **+0.03303** [0.02768, 0.03865] and the matched fixed head
by **+0.02459** [0.01937, 0.03028], both clearing Holm at
p 0.0010 against minimum detectable effects of 0.0114 and 0.0066 fixed before the fold was
opened. Reading, fixed in advance: adapting the trunk adds skill at matched input, not merely
that a trained head helps. **Fold 0 is spent and no further analysis may be run against it.**
Fold 0 was closed to adapter training and selection, not untouched: fixed-recipe predictions
existed for 5,822 of its rows. Guard P.1's reproducibility budget was set after fold-1 drift
was seen; every other element was sealed beforehand. The 1-Mb contrast remains unresolved and
is not part of this confirmation.
The recipe was chosen on folds 1-4 by a leave-one-fold-out selection that was not nested:
each outer fold's recipe was chosen from validation fits on the other development folds, and
those fits had been trained on the outer fold (review item 12). The fold 1-4 selected-recipe
results are therefore exploratory development estimates, and the prespecified confirmation
claim rests on this fold-0 contrast only.

**The chromatin cis head's increment has plateaued; its absolute skill has not.** Against the
globals-only baseline the gain rises 0.0200 to 0.0946 across 20 to 78 training participants per
fold and looks unsaturated. Against reduced-rank regression, the strongest global model, the
cis-specific gain is flat from 31 participants per fold onward: 0.0242, 0.0225, 0.0234, 0.0215
and 0.0232. RRR alone accounts for 83% of the nested arm's rise. The recorded statement that
this head has not plateaued describes absolute skill and the globals contrast, and should not
be read as a donor-limitation claim for the cis-specific increment without the RRR row beside
it. RRR's own rank grid binds in 3 of 5 folds at full size, so +0.0232 is itself an upper bound.

**caQTL: a restricted causal-lead oracle sensitivity, not a task-wide bound.** Measurement
error alone leaves a simulated ceiling of 0.973 and a selection-aware empirical-Bayes ceiling
of 0.953, both uninformative. The restricted causal-lead oracle sensitivity at the 1-kb tested
window assumes that each lead is causal with probability equal to its PIP, that a lead drawn
as causal is predicted with its oracle effect, and that every lead not drawn as causal is
predicted by one constant per fold (the fold mean of the causal draws). Real non-causal leads
can carry predictable marginal association through LD with nearby causal variants, so this is
a sensitivity conditional on that restriction, not a limit for every sequence model. On all
32,322 leads it gives **0.807** [0.804, 0.810]. On the 3,627 leads where both models were
scored on the same folds, the same simulation gives **0.8774** [0.8699, 0.8846], and
AlphaGenome local zero-shot ATAC reaches 0.7441 against HyenaDNA's 0.2344. Intervals are
2.5-97.5% ranges over 2,000 simulation draws (`0b_caqtl_ceiling/out_21845878/run/ceilings.tsv`).

**Severity: label-ceiling simulations are historical sensitivity scenarios, not calibrated ceilings.**
The simulation assumes an independent label-error term using source-reported
repeat-biopsy diagnostics separated in time: 58 identity-unqualified GSE193066 pairs, a
median two years apart (within-participant SD 0.6433 stages), and 30 PXD051911 arm A2 pairs
(0.4655). Serial biopsies include biological change and biopsy sampling as well as reread
error, so neither value is a direct estimate of single-occasion reader error.
GSE193066 identity remains unresolved, excluding its SD 0.6433 and derived
ceilings from primary calibrated claims. The technical
report names 0.4655 as the conservative scenario for a claim that the instrument approaches its
ceiling. For the three cohorts where v1.2 is scored (`latent_s_from_source` transport; primary
GSE268273 split for the two grouped cohorts; ceiling intervals are 2.5-97.5% ranges over 4,000
draws; `0a_severity_ceiling/out/label_ceilings.tsv`):

| cohort | n | achieved Spearman | recorded-stage tie ceiling | simulated ceiling, SD 0.6433 | achieved / ceiling | simulated ceiling, SD 0.4655 | achieved / ceiling |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE268273 | 109 | 0.6631 | 0.9523 | 0.799 [0.728, 0.858] | 0.830 | 0.876 [0.830, 0.913] | 0.757 |
| GSE213621 | 299 | 0.6382 | 0.9423 | 0.779 [0.728, 0.822] | 0.820 | 0.873 [0.847, 0.895] | 0.731 |
| GSE276114 | 177 | 0.6739 | 0.8483 | 0.710 [0.651, 0.763] | 0.949 | 0.761 [0.716, 0.800] | 0.885 |

In these historical scenarios, ceiling-normalised spread within v1.2 is 0.130 at SD 0.6433 and 0.154 at SD 0.4655, against
0.331 and 0.337 when v1 and v1.2 are mixed, so under either scenario the cross-cohort spread is
substantially instrument version rather than cohort difficulty. The ceiling covers the label
only; the RNA sample's own biopsy-to-biopsy variation cannot be measured here, because only 3
of 58 GSE193066 repeat-biopsy pairs are inside the acceptance contract.

**Two heads whose penalty sat at the grid ceiling behaved oppositely.** The LX-2
allele-by-treatment head was under-shrunk: at the shipped ceiling it was worse than predicting
no interaction (-0.00247) and with the grid extended it collapses onto the null (-0.0000138)
with a sevenfold tighter interval, so the deposited value was a grid artefact and the
defensible statement is that the interaction is not learnable from these fixed-weight embeddings. The
HyenaDNA caQTL delta ridge did not move at all: the wide grid selects 1e4 in all 25 fits and
the macro Spearman is identical to 16 digits.

**Current Figure 7C/F source-exposure correction, 2026-10-14.** The exact
evaluated TranscriptFormer TF-Sapiens checkpoint matches all 30 prediction
receipts and the 50,000-cell embedding receipt. Its complete 644-dataset
inventory identifies five clean studies and two training-seen studies
(GSE185477 and Guilliams Liver_Atlas), rather than seven unresolved studies.
Current 7C/F use these states. Geneformer V2 retains zero qualified clean
studies; the historical statement that two displayed models have none is
superseded. Pooled TranscriptFormer macro F1 remains 0.9388747 over all seven
studies and five broad cell classes, so it is not clean-only or disease-stage
performance. All existing scores, gains and intervals are unchanged.

The native clean-subset comparison is reaggregated from the original saved
10,000 paired donor-bootstrap draws (seed 20260824), after reproducing the
existing UCE subset to 1e-12. On five clean studies and 80 donor identifiers,
TranscriptFormer minus the liver encoder macro-F1 difference is -0.004170
(95% conditional interval -0.008732 to -0.000436) for the linear head and
-0.003490 (-0.008311 to -0.000181) for the two-layer head. These are descriptive
fixed-model comparisons, conditional on the five studies and evaluable cell
classes, not population-of-studies generalization or new significance claims.
The liver encoder's fifteen seed/fold receipts confirm disjoint training/query
studies, no held-query adaptation and training-only feature selection.
[Source tables and checks](../Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/README.md).

**Figure 7B, 7E and 7F legend values.** 7B: across the 7 displayed reporter-assay rows (6 models), HyenaDNA gains 0.082008 Spearman over the allele-identity ridge (95% interval 0.027427 to 0.138918; Bonferroni-corrected percentile lower bound 0.008947 from the same 10,000 draws); it is the only row whose corrected lower bound is above zero (`Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7a/exec/mpra/mpra_bonferroni.tsv`). 7E: in GSE83452 the diet arm has 6 resolved and 16 persisted participants and the surgery arm 14 and 1 (`Analysis/MASLD_Model_Benchmark/executions/gse83452-intervention-response-20260831/results_20260831/results.json`, `histology.P2_by_arm`). 7F: the smallest held-out study, GSE189600, contributes 2 donors and 149 of 50,000 cells (0.3%) (`Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/study_metrics.tsv`).

**Historical Figure 7E exposure candidate (earlier callout).** Resolving TranscriptFormer's declared corpus by CELLxGENE
dataset UUID moves the panel from 44 clean / 18 confounded / 22 unresolved to **54 / 22 / 8**.
The eight remaining unresolved cells are Geneformer V2, kept unresolved because its corpus
table carries DOIs rather than dataset ids and absence at DOI level does not establish
cleanliness. Source:
`Analysis/MASLD_Model_Benchmark/executions/encoder-exposure-step4-20260922T175715Z/`.

**Two defects found, neither corrected here.** The Resource study `Liver_Atlas` is Guilliams
2022 (GSE192742/GSE192740, CELLxGENE e84f2780, donors H02-H38), not MacParland 2018
(GSE115469); six model exposure crosswalks and `CONTAMINATION_AUDIT_2026-08-23.md` treat it as
MacParland and record its public date as 2018-10-22 instead of 2022-01-11. All 42 of its
libraries are human runs in the Guilliams SRA study and all 19 donors match the Guilliams
deposit and none matches MacParland. No Figure 7E cell changes, but the recorded evidence is
about the wrong study and public date is load-bearing for exposure. Separately, the caQTL gain
is printed as `0.1301 [0.1162, 0.1440]` in
`docs/technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md:155` and in the c2 execution's own
results, where its producing table records `0.130152 [0.116147, 0.143989]`; correctly rounded
that is `0.1302 [0.1161, 0.1440]`. Neither defect changes a claim.

### Regulatory-model evidence corrections, 2026-09-15

These corrections retain the original comparisons and their candidate status.
Local AlphaGenome macro Spearman **0.7441** refers to **3,627** matched Currin
variants, not the full 32,322-variant population. The reporter variant-bin gain
is secondary development evidence; the original pooled-region comparison did
not establish superiority. HepG2 PAOA allele-effect reliability **0.6938** has
interval **[0.6356, 0.7383]**, which contains 0.70.

Signed-splice direction **0.5777**, interval **[0.5405, 0.6152]**, on **637 pairs
in 247 blocks** is weak directional information. It fails both inconsistent
original criteria: P1 required >0.65, excess >0.05 and block P<0.01;
the D1 wording used >0.60 and P<0.05. Its magnitude correlation is **0.0357**;
the failed ranking/event comparison remains unchanged. The original claim of
“no direction” is too strong.

An exact-duplicate correction retains **621 unique scored pairs** after removing
16 additional copies from those 637 rows; the **247 genomic blocks** are
unchanged. Concordance is **0.5765** versus marginal sign expectation **0.5002**,
with block-bootstrap interval **[0.5383, 0.6135]**, original nominal
block-majority P=**0.000301**, and magnitude Spearman **0.0284**. The original
2,000-draw, seed-123 procedure and inconsistent thresholds are retained; this
is a weighting correction, not new biological replication. Event detection
still falls below proximity on 603 matched pairs. Independent reconstruction
passed for direction, its interval and nominal P, and magnitude in
`GWAS/finemapping/results/alphagenome_campaign/week1-20260915/coordinator-splice-21773779/`;
the complete comparison is in `review/splice-dedup-21773602/` under the same
campaign root. Historical outputs remain unchanged. Exact source-cluster
denominators, gained/lost junctions and new reference-allele checks remain
unfinished, so this does not evaluate an improved event definition.

The historical Corgi differential skills used scaling fitted on held outcomes.
Its RNA comparator also trained on every output region. Those reconstruction
values are not training-calibrated donor × region prediction errors. The new
candidate comparison separates those evaluation populations and fits scaling
and features within training partitions. Exact source records and independent
reconstruction checks are in
`GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/`
and the campaign Catalog build's `independent_prediction_checks/`.

### Initial regulatory context comparisons, 2026-09-15 (not adopted)

The full-source reporter comparison uses **4,358 of 4,359 source constructs**, **239
long-range locus blocks**, and **four paired experimental replicates** per
condition. Replicates are cell-line experiments, not participants. Compression,
scaling and penalty selection use nested training folds. No fitted interaction
arm has lower point-estimate error than zero interaction.

| Interaction model | LX-2 RMSE | HepG2 RMSE |
|---|---:|---:|
| Zero interaction | 0.424805 | 0.238265 |
| Direct interaction | 0.426949 | 0.239220 |
| Separate condition heads | 0.426901 | 0.239220 |
| Shared head plus condition residual | 0.425267 | 0.238565 |

Units are log2 reporter-activity treatment contrasts. Direct and separate-head
MSE-improvement intervals versus zero are below zero in both cell lines;
the shared-residual intervals include zero: LX-2 **[−0.00085043, 0.00003659]**,
HepG2 **[−0.00031214, 0.00000494]**. These are descriptive intervals from 2,000 paired locus and
joint experimental-replicate resamples, preserving shared-control covariance;
they omit fitting/selection uncertainty and support no confirmatory claim.
The source is `GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/lx2_full_21773115/`.
The initial 1,033-construct result remains in `model/lx2_21772727/`. Its zero
RMSEs were 0.418884 (LX-2) and 0.223411 (HepG2), also better than every fitted
arm. The new and archived representations of these 1,033 constructs differ
slightly (relative L2 difference 0.004003; maximum absolute difference 0.078125).
Thus changes between the two comparisons cannot be attributed solely to the
larger population; every within-comparison arm uses identical representations.
After allele subtraction and forward/reverse-complement averaging, the fitted
effect features differ by relative L2 **0.030797** (correlation **0.999526**).
All input sequences, variant bins and declared checkpoint paths match; the
numerical cause is unresolved.

The donor pilot uses **99 GSE267145 participants**, **128 regions in 64 tiles**,
and simple linear RNA and counted-interval sequence features. It is not a
full nonlinear context-model comparison. All 85 native-scale error summaries
were independently reconstructed from source counts and saved predictions.

| Donor-difference model and evaluation | Skill versus zero difference | MSE difference, 95% interval |
|---|---:|---:|
| Shared pairwise, held donors at trained regions | 0.011212 | −0.003243 [−0.009975, 0.003470] |
| Same RNA projection, region-specific model, held donors at trained regions | 0.037220 | −0.010766 [−0.027491, 0.004918] |
| Shared pairwise, held donors at 64 held regions | −0.007252 | 0.002343 [−0.006394, 0.009721] |

Negative MSE differences favor the model. Units are squared log2(1+H3K27ac
CPM), using full source-library totals and training-only fitting. Intervals
use 1,000 donor-by-tile resamples and are descriptive, without multiplicity
claims. The region-specific RNA model has no matched held-region result.
The exclusion ladder includes overlapping genes, ±1 Mb, chromosome, matched
random exclusions, technical inputs and shuffled RNA. Exact target exclusion
and batch-stratified shuffling lack their respective source annotations.
These results do not establish trans regulation or a patient treatment response.
Sources: `GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/b2-21771990/`
and `coordinator-check-21772825/` under the same campaign.

The stronger **99-participant** comparison is complete and independently
reconstructed. At **128 previously trained regions**, nested-selected RNA has
RMSE **0.352445**, versus **0.382406** for training means, in log2(1+H3K27ac
CPM). The paired donor-by-region-component MSE reduction is **0.022017**,
95% interval **[0.013319, 0.032763]**, BH **q=0.001399** in the complete
56-contrast trained-region family. Pairwise FiLM does not surpass this RNA
baseline: MSE improvement **−0.017042 [−0.023892, −0.010485]**, q=0.001399.
Across held donors and **64 held regions**, sequence-only RMSE is **1.291878**,
unrestricted FiLM-profile **1.293698**, and product-kernel **1.312455**. FiLM's
MSE improvement over product-kernel is **0.048884 [−0.224284, 0.311642]**,
q=0.946422 in the complete 40-contrast held-region family: no established
context benefit. Chromosome-excluded selected RNA retains trained-region
prediction (RMSE **0.352545**; dimension-matched random exclusion **0.352674**).
All **30 trained-region and 22 held-region exclusion contrasts have BH q=1**;
this does not establish equivalence or trans regulation. Intervals condition
on fitted predictions and omit fitting and selection variability. Independent
review reconstructed 104 point summaries and 174 contrasts; neural fitting
and inner selection were not repeated. Sources: `model/donor_context_21773708/`
and `review/donor-context-21774201/` under the campaign.

**A three-point training curve shows the 128-region null was support-limited,
and that a non-neural baseline leads throughout.** The same comparison was refit
at **128**, **1,024** and **8,192** training regions against one byte-identical
**1,024-region** held panel, hash-verified across all three fixtures, with 99
donors and 2,760 neural fits at each point. The leading arm's mean region-donor
Spearman rises **0.0907 → 0.1249 → 0.1547**, RMSE falls **1.1914 → 1.0825 →
1.0553**, and the donor-difference MSE falls **0.3005 → 0.2931 → 0.2919**. The
gain per eightfold increase is **+0.034** then **+0.030**, so the curve is nearly
linear in log training size and has not plateaued; the earlier 128-region null
was at least partly a shortage of training support rather than an absence of
signal.

The arm that leads is `product_kernel`, a non-neural baseline, at every point. At
8,192 regions **all 32** architecture contrasts run against it are negative, with
q as low as **0.001665**: the profile-objective arms lose most (**−0.063** to
**−0.028** in donor-difference MSE) and the residual and pairwise neural arms by
**−0.011**. No neural arm overtakes the baseline at any training size tested,
though the gap narrows from **0.0345** to **0.026**.

**The target-chromosome exclusion remains a null, and the bound tightens roughly
twelvefold.** For the leading arm the contrast moves from **+0.00015
[−0.00026, +0.00063]** at 128 regions to **−0.00002 [−0.00005, +0.00000]** at
8,192: more training support bought a much tighter bound on the same negative,
not a signal. Of the 22 exclusion contrasts, 3 reach BH q below 0.05 at 128
regions, **0** at 1,024 and **1** at 8,192, and none of these is coherent. The
three at 128 disagree in direction between the annotated-gene and
random-chromosome comparators of the same arm, and the single hit at 8,192
(`additive_nn_profile`, **−0.0317**) is contradicted by that arm's own
random-chromosome control (**+0.0031**, q 0.275). Nothing here establishes
equivalence or trans regulation.

One dissociation belongs on the record because it inverts the ranking. At 8,192
regions `sequence_only` has the **best RMSE of any arm** (**1.0353** against the
leading arm's **1.0553**) and the **worst** region-donor Spearman (**−0.0784**).
RMSE is dominated by how well a model places a region's overall level, the donor
quantity by variation between donors within a region; reporting RMSE alone would
have ranked the sequence-only control first. All intervals carry the deposited
unit, donors within fold and genomic components, and exclude fitting and
selection variability. Sources: `model/donor_curve_fixtures_21784771/` and
`model/donor_curve_fit_{128_21784884,1024_21784885,8192_21784886}/` under the
campaign.

The corrected local Catalog specimen contains **61,316 evidence records**, with
**7,626 variant**, **6,613 gene**, **186 region**, **340 event**, **117 program**
and **10 context** identities (plus 989 source signals). These are attachment
coverage counts, not newly supported genes. Unknown splice identities do not
form shared unnamed entities. **878** allelic-accessibility records affected
by corrected allele identities remain masked; **120** available records passed
the deposited-allele identity check. Single-SNV and synthetic phased two-SNV
inference examples each return seven native ATAC/DNase tracks, with repeated
reference predictions identical. Synthetic phase is not natural-haplotype
validation. This initial specimen is `catalog-21772728/`, with inference
examples in `inference-21772546/`, under the campaign directory.

The independently checked extension `catalog-currin-21773574/` preserves every
initial record and adds **32,322 measured Currin caQTL associations** and
**343,200 out-of-fold development predictions** from all twelve recipes on
28,600 matched variants. It contains **436,838 evidence records** at **39,766
variant**, **6,613 gene**, **32,508 region**, **340 event**, **117 program** and
**11 context** identities (plus the same 989 signals). Every measured SE joins
the same variant and peak and agrees with the source nominal variance. No
measurement SE is assigned to a prediction. Fold-specific training exposure,
source-beta units and missing individual prediction uncertainty remain visible.
These are stored benchmark predictions: fitted head weights were not saved or
deployed. No new gene links or genetic memberships are inferred. Independent
checks: `review/catalog-currin-21773580/` under the campaign.

The twelve fixed-representation recipes completed on **28,600 identical
variants in 2,628 genomic blocks**, with all 60 fold fits successful. Both input
lengths retain the same comparison rows. Length-specific coverage remains
separate: 2,048-bp target pooling covers 28,616 variants; 16,384-bp pooling
covers 32,306, with 16 sequence exclusions. Independent reconstruction matched
all twelve error and correlation summaries. Macro Spearman uses the established
Fisher-z mean of the five fold correlations.
The coordinate audit found zero cross-fold input overlaps at all three lengths.
Historical 1-Mb bootstrap bins are not verified independent LD units: 252,442
pairs of overlapping native windows cross those bin boundaries. The completed
native comparison therefore also reports whole-chromosome uncertainty, with
identical predictions and separate, explicitly named exploratory families.
Only 22 chromosomes are available, fewer in the single held adaptation fold;
this sensitivity does not establish pretraining or source-cohort independence.
Exact variant-and-peak matching excludes 94 of 1,696 specialist-source records
that refer to another measured peak, leaving 1,602 before score-coverage checks.
The excluded records retain their own source targets and effects.

Native 1-Mb attempt accounting is complete for all **32,322 eligible rows**:
**30,611 scored**, **1,711 excluded**, and **zero unresolved**. This is coverage
accounting, not a full-population superiority result. The downstream matched
comparison initially stopped on an unnamed variant-index column after scoring
completed. The repaired rerun and precision strata are complete, retaining the
failed output and reusing saved predictions. Independent numerical reconstruction
of this latest comparison has not run. The deposited table reports native ATAC
macro Spearman **0.584373** over the **30,611 scored variants**, with
training-calibrated RMSE **0.544389** in source caQTL-beta units. The earlier
0.7441 remains a different, smaller-population result.

On **6,441 shared validation variants**, the matched table gives native ATAC
RMSE **0.536015** and signed Spearman **0.579456**, versus **0.545425** and
**0.564337** for the lowest-RMSE adapter (rank 16, last three blocks, symmetric
pooling). That table scores native at 1 Mb against 2-kb adaptation, so it does
not test adaptation at equal context length. Sources:
`model/native/completion.json`, `model/final_21775091/` and
`review/precision-21775092/` under the campaign.

The equal-context-length test has now run, and it separates the two effects.
Native scoring was repeated at 2,048 and 16,384 bp with the same 501-bp
ATAC/DNase readout, hardware and training calibration, on the same 25,748
permitted rows. Native ATAC signed Spearman rises with input length: **0.5396**
at 2 kb, **0.5689** at 16 kb, **0.5795** at 1 Mb. Against native at its own 2-kb
input the adapter improves MSE by **+0.018150** [+0.009541, +0.027674], BH
q **0.0017**; against native at 1 Mb the same adapter loses **0.010176**
[0.003534, 0.017302], q **0.0089**. Both signs survive whole-chromosome
resampling over five units, and the equal-length contrast repeats on the
unrestricted 2-kb population of **6,834 variants in 524 blocks**
(**+0.016877** [+0.008332, +0.024932], q **0.013**).

Adaptation therefore does add skill at equal input length, and does not recover
what 1 Mb of native context supplies.

Those numbers were one seed on one validation fold. The sweep has now been run at
**five seeds across the four admissible validation folds**, 240 fits in all, with
**fold 0 never trained on and never validated on** so that one fold stays
untouched for a single later confirmation. Three things separate that a
single-seed table could not. **These four-fold selected-recipe results are
exploratory development estimates** (review item 12, 2026-09-23): the
leave-one-fold-out selection chose each outer fold's recipe from validation fits
that had been trained on that outer fold, so the selection is not nested. The
prespecified confirmation claim rests on the fold-0 contrast only.

First, the matched-length result is robust. The leave-one-fold-out selected arm,
as a five-run ensemble, beats native ATAC at its own 2-kb input in **all four
folds**: **+0.02579**, **+0.02193**, **+0.03810** and **+0.02921** in MSE, mean
**+0.02875**, every interval excluding zero and every p at the 1/2,001 resampling
floor.

Second, the comparison against 1-Mb native does not resolve. Across the four
folds it gives **−0.00254** (p 0.36), **−0.01016** (p 0.010), **+0.00551**
(p 0.027) and **−0.00055** (p 0.85), mean **−0.00193**. Two folds have intervals
excluding zero and they point in **opposite directions**, so which fold
validates, not the model, sets the sign. The single-seed value reported above is
confirmed rather than superseded: the same arm on the same fold gives
**−0.01021** [−0.01661, −0.00403] as a five-seed ensemble against **−0.010176**
[−0.017302, −0.003534] from the single fit. It was one fold of four, and the
least favourable one for adaptation.

Third, the arm ranking was never resolvable at one seed. The across-seed standard
deviation of MSE is **0.0089 to 0.0135** depending on the recipe, while the gap
between the two best recipes is **0.0030**. Over four folds and five seeds the
lowest-MSE recipe is **rank 16, last five blocks, symmetric pooling**
(**0.31391**), not the rank-16/last-three/symmetric arm the single-seed table
named (**0.31688**). Leave-one-fold-out selection does pick the last-five arm in
all four folds, so the selection is stable even though the per-fold winner is
not: fold 2 prefers last-three. Averaging five runs is worth **−0.0115** MSE on
its own, more than the entire gap between the top two recipes, so a five-run
ensemble and a single run are reported as separate quantities and never
interchanged.

What survives at the class level is the ordering, not the winner: all eight
adapters beat both partial query/value updates and both matched fixed-head
controls in the four-fold macro table. Sources:
`model/nested_f1_21783075/` through `model/nested_f4_21783078/` and
`model/nested_selection_21784704/` under the campaign. No finalist is declared,
the fold-0 confirmation has not been run, and nothing is adopted. The earlier reading, that the table did
not support adapter superiority to native ATAC, described the unequal-length
comparison only. The limits are unchanged: one adaptation seed, one held
validation fold, five validation chromosomes, historical 1-Mb bins not verified
as LD units, a significance-selected source population, and a native track
readout that remains a different object from a learned scalar head even at
matched input. Independent numerical reconstruction reproduced all **94 MSE
contrasts** and **64 model points** to a maximum absolute difference of
**3.6e-15**, importing no producer inferential function. Sources:
`model/native_lengths_21775189/`, `model/native_short_compare_21775217/` and
`review/native-lengths-21783036/` under the campaign.

The bounded adaptation sweep completed **12 arms × 5,000 steps**, one seed,
on **6,834 validation variants in 524 genomic blocks**. The eight adapter
Spearman values range from **0.541389 to 0.558563**; matched fixed variant and
symmetric heads give **0.528753** and **0.526776**, respectively. Adapter RMSEs
range from **0.546808 to 0.560796** in source caQTL-beta units, versus
**0.565534** and **0.567125** for those controls. The two partial query/value
updates give Spearman **0.525718** and **0.527409**. A post-pilot paired analysis
now supports **six of eight adapter MSE improvements** under both historical-bin
and whole-chromosome resampling: their 95% intervals exclude zero and BH q is
below 0.05 in each complete ten-contrast exploratory family. Neither partial
update supports improvement; the 3e-6 update increases error. The largest
observed adapter reduction is **0.020317 RMSE units**, from **0.567125 to
0.546808**, with chromosome-bootstrap interval **[0.012041, 0.024484]**.
Only **five validation chromosomes** support that sensitivity. The 2,000-draw
intervals condition on measured labels and the fitted single-split, single-seed
models; source effect-estimation error, refitting and selection uncertainty
are absent. Independent reconstruction reproduced all twelve RMSE/correlation
points and twenty MSE/RMSE interval and BH calculations. Five-seed nested
selection and superiority to the native model remain unresolved. Sources are
the representation-comparison files indexed by the campaign README and
`model/sweep_summary_21773124/`, `review/adaptation-precision-21773558/` and
`coordinator-adaptation-21773575/` under the campaign.
The twenty Gaussian likelihood arms are also complete, with exact source-SE
variant-and-peak joins independently checked. On the same **6,834 validation
rows**, source-SE Gaussian versus squared error gives **one improvement, five
deteriorations and four unresolved contrasts**; versus constant-SE Gaussian,
**two improvements, two deteriorations and six unresolved**. These counts
require the 95% MSE interval and BH q<0.05 and agree across historical-bin
and chromosome resampling, each with its complete ten-contrast family.
Rank-4/last-three/symmetric improves RMSE from **0.560796 to 0.553774** versus
its squared-error counterpart, but signed Spearman declines from **0.543506
to 0.535913**. The best observed RMSE among all thirty compared fits remains
**0.546808**, the prior rank-16/last-three/symmetric squared-error adapter.
This one-split, one-seed comparison conditions on fitted predictions; the
likelihood variance is not validated individual prediction uncertainty.
Source: `review/effect-likelihood-21774090/` under the campaign.

The separate compiled-route zero-adapter output check has passed: all **56**
comparisons are exactly equal on **32 development sequences**, across four
adapter configurations and two pooling definitions. Repeated calls and fresh
compilations also agree exactly. Independent review reconstructed all 32
scalar contrasts. This supports the tested forward implementations; it does
not establish general numerical equivalence, explain the separate archived/new
reporter-feature difference, or replace the pending full native comparison
and five-seed nested selection. Sources: `model/zero_outputs_21773146/` and
`review/zero-outputs-21773161/` under the campaign.

Measured native ATAC performance has now been assessed separately from output
drift. The exact source matrix contains **349,685 counted peaks × 138 donor
columns**. A coordinate-hash sample selected **512 calibration and 512 validation
regions** before reading activity values, within source-ascertained consensus
peaks no wider than 2,048 bp. All eleven original/adapted backbones completed
the same reference windows and unchanged adult liver ATAC heads. Each model's
affine calibration uses training regions only. The target is donor-mean
log2(1+CPM), using every deposited peak for each donor's library denominator;
off-peak reads and source-wide VST/GC corrections are not used.

| Regional ATAC predictor | Validation RMSE, log2(1+peak CPM) | Spearman |
|---|---:|---:|
| Training-mean target | 1.073990 | Undefined for a constant prediction |
| Original backbone, training-calibrated | 0.741495 | 0.676903 |
| Rank-16, last-three-block symmetric adapter | 0.726547 | 0.673187 |

The displayed adapter had the largest allele-effect RMSE gain in the prior
pilot; all eleven backbones remain in the complete ATAC table. Its ATAC MSE
improvement is **0.021944**, with intervals **[−0.035837, 0.091658]** under
region-component resampling and **[−0.035184, 0.066360]** under chromosome
resampling. This does not establish either improvement or noninferiority.
The rank-4/last-three/symmetric adapter and both partial-update pilots increase
MSE with BH q<0.05 under both schemes. The partial-update RMSEs are **0.745696**
and **0.759436**. The original model clearly exceeds the training-mean baseline
in this comparison; it does not establish superiority to a width/GC baseline,
which was not run.

Each scheme uses 1,000 paired donor draws shared across calibration and
validation, with separately resampled training/validation regions and affine
refitting. There are 512 validation window components and only five validation
chromosomes. BH is applied separately to each complete eleven-contrast
exploratory endpoint/resampling family. These intervals condition on the
fitted adaptations and omit adaptation-selection uncertainty. Independent
weighted reconstruction passed all twelve point summaries and twenty-two MSE
interval/P/BH contrasts; ranking intervals were not independently reconstructed.
Sources: `model/native_retention_21774009/` and
`review/native-retention-21774089/` under the campaign. This is same-source,
2-kb regional adult-liver ATAC assessment, not independent-cohort validation,
donor prediction, or whole-model native-function accuracy.

Natural-background source checks now resolve the **138-column author-deposited
Currin accession map to GSE26105**. Independent reconstruction from the count
object and public accession lists passed. No final corrected phased-call or
phase-quality asset was listed in the inspected Currin GEO/Zenodo deposits;
the accession map does not independently verify genotypes or phase.

A bounded ATAC information census, independently recounted from six indexed
BAM intervals, found the following direct two-site coverage. Denominators are
distinct barcode/fragment geometries with a callable allele at either site,
conditional on the fixed fetched interval and MAPQ≥30/base quality≥20.

| Source and donor | GRCh38 site pair, 1-based positions | Direct two-site molecules / callable molecules |
|---|---|---:|
| GSE281367 Z01 | chr8:9325848 A>G; chr8:9326086 A>G | 0/37 |
| GSE281367 Z04 | Same chr8 pair | 0/56 |
| GSE281367 Z05 | Same chr8 pair | 0/37 |
| GSE281367 Z08 | Same chr8 pair | 3/93 |
| GSE244832 D15 | chr19:45326124 G>A; chr19:45326536 G>C | 0/11 |
| GSE244832 D17 | chr22:43928847 C>G; chr22:43928850 C>T | 2/3 |

Z08 has three mate-supported ALT/ALT observations; D17 has one REF/REF and one
ALT/ALT observation. D17's 47 alignment records collapse to three fragment
geometries, so reads cannot be treated as independent phase observations.
None of the seven inspected Z04 templates excluded for inconsistent template
length supplied calls at both sites. Absent BAM sample tags in D15/D17 leave
donor identity based on curated paths and source metadata. These are selected
donor-by-pair information summaries, without independent genotype truth,
mapping-bias correction, a natural interaction test, or a change to the GNMT
null. Sources: `data/read-phase-21773997/`, `data/currin-donor-crosswalk-21773997/`,
`review/phase-counts-21774168/` and `review/currin-crosswalk-21774168/` under the
campaign. Further threshold changes do not supply the missing phase evidence.

The target-aware transfer comparison retains **2,765 variant–peak pairs at
1,636 variants in 172 long-range blocks**, without choosing targets by
association significance. The physical-overlap subset contains **112 variants
in 63 blocks**. The initial procedure selected zero predictions for both
methods in all five folds; its identical-prediction intervals cannot establish
an information limit. An explicitly recorded development amendment compares
the best learned candidates from the same training grid, retaining the original
zero-selection results. All-pair specialist/shared RMSEs are **0.474934/0.473595**;
their MSE improvement is **0.001271 [−0.000319, 0.002499]**. Physical-overlap
improvement is **0.002504 [−0.004673, 0.009555]**. All nine estimable stratum
intervals include zero; none passes BH correction across the complete
12-stratum exploratory family (three strata are untestable). These 2,000
paired-block intervals condition on measured source effect estimates and
fitted predictions; they omit measurement and refitting uncertainty. They do
not exclude useful transfer from other models or data. Source:
`model/transfer_r2_21773128/`, with the original in `model/transfer_21773114/`.

Candidate Figure 7 comparisons are in
`figures/candidates/regulatory-week1-20260915-v4/`; the earlier subset panels
remain in v3. The adaptation uncertainty panel is in
`figures/candidates/regulatory-adaptation-week1-20260915-v2/`; its printed
reductions and rendered intervals were inspected against the checked table.
The Figure 6 query specimen is in the campaign's
`figures/Figure6-21753521/`. The additional measured/predicted variant example
is in `figures/candidates/regulatory-query-week1-20260915-v1/`; all thirteen
printed effects match the independently checked Catalog query, and the PDF
was rendered and inspected. The example follows source order, with all twelve
recipes and unavailable prediction uncertainty shown explicitly. It is not
selected as a successful prediction or established gene mechanism.
Native scoring and the matched model uncertainty are complete and
independently reconstructed. Protected outcomes remain closed; none of these
candidates is adopted.

The target-gene branch is complete as an inventory and is underpowered as a
comparison. On the pinned DC-TAP source, the geometry and uniform-power rules
admit **1,302 element-gene pairs over 372 elements and 14 chromosomes**, and
those rules leave **19 positive pairs** for the expression-decrease endpoint and
**27** for any measured change. Fixed distance alone reaches average precision
**0.1966** against a constant-score reference of **0.0146**, but the
chromosome-bootstrap interval is **[0.0452, 0.6511]**: an order of magnitude
wide, so no model comparison placed against this baseline can resolve anything.
Target discrimination is weaker still, with **11 of 240** multi-target elements
carrying a detected target and a tie-averaged top-ranked precision of **0.0167**.
The assayed context is K562, not liver. The branch is therefore reported as a
power bound and no model is fitted to it; adding one would produce an interval
that contains the baseline by construction. Sources:
`data/target-admission-v1/` and its `followup-21775222/` under the campaign.

The splice branch is bounded by its population, not by its method. Exact source
identity resolved **221,828 variant-event pairs** over **1,485 unique events**
in **1,414 clusters**, with zero clusters carrying conflicting gene-annotation
strands, and the historical table joined **630 of 657** rows to an exact event
identity. But `new_effect_evaluation_ready_rows` is **0**, because the only GTEx
v8 sQTL data on disk is the significant-pair table: a population selected on the
outcome, with no all-tested denominator. The per-SNP all-associations files are
not published in the public GTEx bucket, whose
`sQTL-Tissue-Specific-All-SNP-Gene-Associations/` prefix holds a zero-byte
placeholder. The tested **phenotype** universe is published, as the LeafCutter
intron clusters actually quantified in liver plus the cluster-to-gene grouping,
and is being acquired; until it is in place, no splice effect comparison is run
and the deposited identity work stands on its own. Sources:
`data/splicing-identity-21773601/`, `review/splice-dedup-21773602/` and
`coordinator-splice-21773779/` under the campaign.

### Two-model molecular-value development, 2026-09-22 (not adopted)

The new comparison retains the existing 99-participant GSE267145 RNA-to-H3K27ac
profile model and treats Currin caQTL allele effects as a separate, assay-native
task. It does not reopen spent Currin fold 0 or turn the existing severity
instrument into a new model. The fixed comparison specification is
`scripts/analysis/masld_two_models/design.json`; the candidate analyses are under
`GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/`.

[Observed pattern, not a guarantee] The findings below concern inspected
development sources and do not establish accuracy in a new cohort.

In the inspected five-participant-fold H3K27ac development set, the released
global reduced-rank RNA model has squared-error skill **0.11919** and the
global-plus-local model **0.14133** against fold-training region means. A
histology, sex and RNA-library-descriptor ridge, with its penalty chosen inside
each outer training set, has skill **0.01971**. The local RNA model's paired
skill increment over that comparator is **0.12162** [0.09542, 0.14682]
(participant bootstrap; Holm P=0.00150 across three declared contrasts). The
global-to-local increment is **0.02214** [0.01734, 0.02662]. Adding estimated
hepatocyte mRNA fraction to the comparator changes its skill to **0.0202**.
Those final BayesPrism fractions were jointly estimated across all 99 RNA
profiles, so this is a cohort-adapted composition sensitivity. Improvement
over the specified linear, numeric-grade comparator does not establish general
conditional independence from histology. These are same-cohort prediction
comparisons, not evidence of chromatin mechanism or
external H3K27ac transport. The exact candidate values and fitted nuisance
choices are in `chromatin_conditional_21849192/results.json`; independent
numerical reconstruction of the target, participant axis and overall skills
agreed to numerical precision.

An additional metadata-defined comparison matched 22 disjoint pairs (44
participants) on all four histology grades, sex and evaluation fold. Skill
for their measured chromatin differences was **0.11070** [0.05349, 0.16267]
for global-plus-local RNA and **0.08859** [0.03661, 0.13617] for global RNA.
The paired increment was **0.02211** [0.01423, 0.02871]; the increment over
the technical comparator was **0.10699** [0.05460, 0.15297]. Both centered
bootstrap contrasts had Holm P=0.00020 in this two-contrast sensitivity.
Twenty pairs were women with fibrosis 0; one male pair had fibrosis 0 and
one fibrosis 1. No fibrosis 2–4 pairs were represented. Intervals resample
disjoint pairs conditional on the shared fitted models and omit retraining
uncertainty. Crossed pair/chromosome resampling gave a similar local-model
interval [0.05249, 0.16264]. This supports within-histology prediction only
in the matched low-fibrosis population. Source:
`within_histology_21849320/results.json`.

A bounded retrospective acquisition test used 1,000 outcome-blind,
nearby-gene-annotated regions, five disjoint evaluation folds, 45 base-training
participants per fold and additional-measurement budgets of 10%, 20% and 30%
of the 99-person cohort. The primary budget adds 20 measurements to the 45
already available; 20 represents 54.1–60.6% of the 33–37-person selectable
pool, not 20% of that pool. At this primary budget,
random selection captured **0.5846** of candidate-pool prediction error and
reduced held-participant MSE by **0.00632** after refitting. RNA diversity
captured **0.6481** and reduced MSE by **0.00779**; global/local disagreement
captured **0.6290** and reduced MSE by **0.00768**. Relative to random, the
five-fold paired capture differences are +0.06349 and +0.04438, both with
enumerated two-sided fold sign-flip P=0.0625 and Holm P=0.1875 in their three-policy
family. Their refit advantages are +0.00147 and +0.00136 MSE units; neither
passes the parallel Holm comparison (both P=0.375). The five training sets
overlap, so fold-resampling intervals are a sensitivity summary, not a claim
of independent replications or exact biological-unit inference. The original
histology policy balances counts across exact grade strata with deterministic
ties; it is not proportional stratified random sampling. Global rank-12 RRR
is fitted across all 96,460 responses; local correction and evaluation use
the fixed 1,000-region panel. This bounded recipe differs from the full
released fitted model. Sources are
`chromatin_acquisition_21849197/acquisition_summary.tsv` and
`chromatin_acquisition_21849197/primary20_fold_comparisons.tsv`.

A second calculation gives every evaluation participant equal weight and
averages the three random-selection losses within participant. The primary
20-measurement refit gains over random are **0.001486** [−0.000325, 0.003911]
for RNA diversity, **0.001323** [−0.000615, 0.003844] for disagreement, and
**−0.000497** [−0.001735, 0.000744] for histology balancing. All three
centered participant-bootstrap contrasts have Holm P=0.4962. These 20,000
draws resample 99 participants conditional on fitted predictions; they omit
retraining and random-policy Monte Carlo uncertainty. They do not establish
that any selection policy improves refitting more than random. Source:
`chromatin_secondary_review/acquisition_primary20_paired_participant.tsv`.

The single permitted sequence-by-RNA residual experiment used 1,024 fixed
regions, 99 held participants, fixed sequence features and rank 4/8/16
selection within each outer training fold. Sequence features use the permitted
`counted_interval` interpretation with its recorded 1-bp uncertainty; correspondence
to called-peak boundaries remains unresolved, as specified in the
[sequence-use decision](decisions/2026-09-15-gse267145-sequence-use.md).
The selected interaction's residual
skill is **−0.000563**, versus **−0.000545** for the additive control,
**−0.000841** for sequence-only and **−0.024011** for RNA-only. Interaction
minus additive is **−0.000018** with paired-participant bootstrap interval
[−0.001704, 0.001574]. This increment uses baseline-model residual SSE as
denominator. Independent reconstruction on the planned original chromatin-skill
scale gives **−0.0000159** [−0.001520, 0.001375], below the +0.01 margin.
The negative increment cannot qualify on either scale; the simpler
chromatin model remains selected. This is a held-participant comparison at
modeled regions, not an unseen-region or external-source result. Its base fit
uses rank-12 global RRR across 96,460 responses, with local correction and
evaluation restricted to 1,024 regions.
Inputs, choices and participant-level errors are in
`sequence_rna_residual_21849236/`.

The released profile scorer's numerical input check used 20 already-exposed
GSE267145 participants and therefore makes **no prediction-accuracy claim**.
Its fixed raw transformation gave the same individual's profile alone or
within a ten-person batch (maximum absolute difference <2×10⁻¹⁴), and
multiplying that RNA library by 0.5, 2 or 4 left the raw profile unchanged.
The batch-marginal route changed the same ten predictions when ten other
participants were added (profile RMS difference **0.0481** log2-CPM units;
maximum **2.0852**). With random genes retained at 90%, 75% and 60%, raw
profile RMS shifts were **0.0503**, **0.1101** and **0.1840**. At 90% coverage,
retaining the lowest-abundance genes while omitting the highest changed the
profile much more (RMS **1.1928**) than omitting the lowest-abundance genes
(RMS **0.00174**). Coverage alone therefore does not characterize numerical
stability. These are perturbations of source-training inputs, not measured
accuracy under missing genes. Source:
`profile_input_stability_21849262/input_stability.tsv`.

A held-participant perturbation test used the bounded rank-12 recipe with
1,000 local-correction/evaluation regions, keeping continuous RNA estimates.
Here and below, “RNA mass” means summed input abundance among the 42,163
modeled genes, not physical RNA mass or a whole-transcriptome fraction.
Complete-input raw global-plus-local skill was **0.14743**, versus **0.13952**
with batch-marginal mapping; their paired difference was −0.00791
[−0.01591, 0.00021]. These are not the released tuned-recipe values, whose
skill on this same panel is 0.15782. Random retention of 60% of genes reduced
raw skill to **−0.14621**, while marginal skill was **0.12848**. Removing the
highest-abundance 10% of training-defined genes retained only **9.91%** of
RNA input mass on average and produced raw skill **−7.54401**, versus
marginal **0.13112**. Removing the lowest-abundance 10% retained over
99.9999% of input mass. This deliberate compositional stress does not describe
ordinary 10% random dropout; a gene-count coverage threshold alone is
insufficient to describe robustness. On the same 50 participants, scoring
ten-person batches instead of their full 17–21-person held folds decreased
marginal skill by **0.01470** [0.00589, 0.02381]; raw predictions were unchanged.
Multiplication of library values by 0.5 or 2 was invariant within numerical
precision, but is not a read-thinning or empirical read-depth test. All
intervals are unadjusted, paired fixed-prediction sensitivities conditional
on query-batch composition, with no retraining or external-source uncertainty.
The participant SSEs and contrast points were independently reconstructed.
Sources: `held_input_robustness_21849340/summary.tsv` and
`held_input_robustness_21849340/paired_sensitivity_contrasts.tsv`.

A subsequent adaptive-development check corrected observed-gene CPM by the
expected observed RNA-mass fraction, estimated only from each fold's training
RNA. Missing genes remained at training log-CPM means. This kept individual
scoring independent of a query batch and left complete-input predictions
unchanged. With random 90%, 75% and 60% gene retention, corrected local-model
skills were **0.14651**, **0.13600** and **0.12177**, respectively. Paired gains
over raw normalization were **0.02922** [0.01843, 0.04083], **0.09434**
[0.07520, 0.11415] and **0.26799** [0.23128, 0.30831]. Omitting the
training-highest-abundance 10% gave corrected skill **0.13263**, instead of
raw −7.54401. This correction assumes the remaining genes carry a similar
fraction of RNA mass to the training reference. It supports an optional
individual-sample sensitivity, not a universal replacement for marginal
mapping: most corrected-versus-marginal intervals included zero. The
intervals remain descriptive and unadjusted, and the analysis followed the
observed failure of raw normalization. Only the bounded rank-12, 1,000-region
recipe has held-participant evaluation; extension to the released full model
has no separate accuracy estimate. Source:
`reference_mass_robustness_21849357/`.

The requested proportional histology-stratified random policy was evaluated
separately with three fixed seeds, largest-remainder stratum quotas and the
same acquisition splits, target transformations and model recipe. At 20
additional measurements, its mean captured pool-error fraction was **0.59602**
and its participant-weighted refit MSE reduction was **0.00721**, versus
**0.00629** for random sampling. These are descriptive comparisons; the
original deterministic balancing result remains separately identified.
Source: `histology_sampling_21849335/`.

The native-plus-adapter comparison completed on 2026-09-23 using all 60
planned two-fold, five-seed fits. Evaluation includes **24,430/25,748**
development leads (**94.88%** coverage), across 17 chromosomes and 2,073
historical 1-Mb bins. The 1-Mb eligibility restriction removes 1,318 leads;
all five models are compared on the same retained rows. Currin fold 0 was
not re-evaluated. MSE in squared source ALT-dosage beta units is:

| Model | Matched development MSE |
|---|---:|
| Native 2 kb, training calibrated | 0.330824 |
| Matched fixed head 2 kb | 0.318440 |
| Five-seed adapter 2 kb | 0.302252 |
| Native 1 Mb, training calibrated | 0.300161 |
| Regularized native-plus-adapter combination | 0.290705 |

The combination reduces MSE by **0.009457** [0.008305, 0.010500] versus
native 1 Mb and **0.011548** [0.008407, 0.014633] versus the adapter.
Both chromosome-bootstrap contrasts have Holm P=0.000666 across the two
prespecified comparisons. These are conditional intervals for saved
predictions, excluding refitting, recipe-selection, fitted-seed and runtime
numerical variation. The improvement over the better simple model is
**3.15%**, below the fixed **5%** development complexity margin. Therefore
**retain calibrated native 1 Mb as the selected simple allele model**;
the combination is a positive incremental-information result that does not
meet the escalation criterion. Statistical significance does not change
that decision, and these inspected development results are not external
confirmation. The adapter recipe came from the non-nested fold 1-4 selection
(review item 12), so the adapter and combination values are exploratory
development estimates. Source: `allele_combination_21849200/results.json`.

A descriptive experiment-prioritization check selects the highest absolute
held-fold predictions within each development fold. At the primary 20% budget
it selects 4,887 leads. Native 1 Mb captures **26.43%** of total absolute
measured beta, versus **20.00%** expected under uniform sampling, with
**98.51%** predicted/measured direction concordance among selected leads.
The adapter and combination capture **26.16%** and **26.49%**, respectively.
This is enrichment within an association-significance-selected source,
not an estimate of experimental success, causal-variant detection or
new-variant accuracy. Pointwise chromosome intervals condition on fixed
selected sets and predictions; no superiority tests were performed.
Source: `allele_prioritization_21849355/`. A separate reconstruction of all
nine model/budget selections, estimates and pointwise intervals agreed within
2.8×10⁻¹⁶; details are in `allele_independent_review.json`.

The exact-pair example `examples/allele_pair.json` returns the prespecified
chr6:42,911,813 C>T source association near GNMT: measured beta −0.338619,
native-1-Mb prediction −0.045810 and adapter prediction −0.144173.
These agree in sign but underestimate the measured magnitude. They are
association coefficients, not isolated allele-editing effects. The source
peak is an experimental readout to investigate; its target-gene link and
any MASLD effect remain unsupported.

The 32-variant outcome-blind inference-cost check used the same L40S GPU and
unchanged checkpoint files, but fresh native raw ATAC effects differed from
saved predictions by up to **0.00561**, and fresh adapter beta predictions by
up to **0.00887**. Changing the native requested outputs to ATAC alone did
not resolve the mismatch. An earlier same-card, separate-process native
check had shown similar variation. A completed unchanged-scorer diagnostic
confirmed exact within-process repeats but maximum between-process differences
of **0.01306** raw native ATAC units and **0.00572/0.00347** adapter ensemble
beta units at batch one/four. Original batch-four geometry did not reproduce
the saved adapter predictions. Native direction differed from saved output at
1/32 variants; adapter ensemble direction did not change in this small check.
This is runtime numerical variation with unchanged weights, distinct from
fitted-seed variation. Its exact backend cause remains [Unverified]. These
observed differences do not define a numerical tolerance bound. Fresh sequence
inference remains experimental; the exact-table interface reproduces the saved
evaluation values. Measured warm latencies were **1.627 s** for the evaluated
native ATAC+DNase request and **0.0253 s** for all five adapter seeds per
variant, excluding initialization; they are implementation diagnostics, without
an admitted cost/performance ratio or noninferiority claim. Stored predictions
used for matched accuracy comparisons are unchanged. Sources:
`inference_cost_21849351/` and `reproduction_discriminator_21849377/`.

A public-source inspection found 17 exact donor identifiers shared by processed
RNA and H3K27ac matrices in GSE128072. They are normal liver tissue, the RNA
matrix is TPM, and the H3K27ac intervals require hg19-to-GRCh38 reconciliation.
The source is therefore a candidate for a separate transport design, **not**
an external validation of the released count-input profile model. Wenz
GSE277774 remains unadmitted for allele-effect transfer. Its public
FDR-selected 175-sample file has **14,076** significant peak leads and RASQUAL
Pi, not Currin FastQTL beta or the all-tested denominator. At the stated SNP
positions, **13,688/14,076** reported REF alleles match the local GRCh38
reference, **238** match ALT and **150** match neither; **305** exact variant
tuples overlap Currin leads. This supports GRCh38 interpretation for most
reported rows but does not prove every allele is harmonized. The final 175
donor crosswalk, row-level mismatches and checkpoint exposure remain unresolved,
so no external model score is reported. Source checks are in
`source_metadata/admission.json`, `source_metadata/wenz_fdr5_identity.json` and
`source_metadata/wenz_admission.json`. [Inference] Published collection histories
support donor non-overlap: [Wenz](https://pmc.ncbi.nlm.nih.gov/articles/PMC13012615/)
used BioTIP specimens collected in 2012–2020, whereas
[Currin](https://pmc.ncbi.nlm.nih.gov/articles/PMC12212085/) matched deceased-donor
samples to genotypes documented in 2011. This does not resolve final-175 sample
membership. The potentially informative Wenz Table S1 supplement could not be
retrieved during this review; its identifiers remain unverified. Details are in
`source_metadata/wenz_donor_provenance_update.json`. The completed four-fold
comparison remains single-source candidate evidence; no external superiority
claim is made from these results.

Six individual 6-point vector panels, their source tables and the full caption
are in `figures/candidates/figure7_two_models_20260922_21849207/`.
The independent numerical reconstruction passed for chromatin skills,
matched-pair results, acquisition errors, all five allele MSEs, the 60
stacking-input split records, four separately refitted stacking coefficients,
paired allele intervals and Holm values, and all main-panel estimates.
The six rendered panels also passed text and visual inspection. Records:
`audit_21849207.json` under the campaign root and
`figures/candidates/two-models-partial-proof-20260922/final_21849207/six_panel_review.json`.
The 60 nested GPU fits took 3:12:21 on one L40S with 4 CPUs; peak host RSS
was 5,222,284 KiB. Named model releases and the working Figure 7 are unchanged.

### H&E to molecular programs: external test in non-postmortem liver, 2026-10-03 (not adopted)

The GTEx-fitted image model (Model C, `scripts/analysis/he_molecular/`, PRESPEC_C
Amendment 6) was applied once, without refitting, to 50 TCGA-LIHC donors with a
cryosection H&E slide and RNA-seq from the same non-tumor liver sample. The
prespecified primary endpoint was not met. Unit: donor (n = 50). Covariates in the
partial Spearman: fat-vacuole fraction, age, sex. Null: 1,000 shuffles of the image
prediction within tissue source site. Multiple testing: BH over the 16-program family.

| Quantity | Primary scores | Residualized on RNA negative controls |
|---|---|---|
| Median partial Spearman over 14 testable programs | 0.204 (p = 0.103) | 0.231 (p = 0.001) |
| Programs with BH q < 0.05 and partial Spearman ≥ 0.30 | 4 of 16 | 3 of 16 |

- One program meets the rule in both versions: Myeloid / immune (DOCK2),
  `hotspot_fibroblasts_81a8061005f0bc22` (0.420, q = 0.031; residualized 0.317, q = 0.033).
- Primary only: Macrophage / immune (CD163) 0.383, Hepatocyte stress (SLC38A1) 0.411,
  Resident macrophage (LIPA) 0.485. Residualized only: Complement / plasma (CFH) 0.406,
  Amino-acid metabolism (GOT1) 0.355.
- Two of the 16 are untestable in TCGA (gene coverage 0.74 and 0.79, below 0.80) and
  count as not transferring: Cell motility / filopodium and hepatocyte Complement (CFH).
- Negative controls: hypoxia 0.313 (p = 0.038), immediate-early 0.104 (p = 0.22). Hypoxia
  exceeds the family median, so under Amendment 6 item 7 the result is not separated from
  tissue-handling state.
- Guard: the 45-donor GTEx sealed gains reproduced to 5e-15 before TCGA was read. The 28
  observed statistics were recomputed with a second implementation (max difference 2e-16;
  `c6_verify/`).

Reading: an immune-content signal carries from postmortem paraffin slides to surgical
cryosections; the hepatocyte metabolic programs that passed in GTEx do not carry
consistently. This is reported as not demonstrated in this cohort, not as
`tested_negative`: n = 50, 28 donors from one site, cryosections, and tissue beside
a tumor. Phikon-v2 pretraining included TCGA slides (images only). Non-tumor histology
was not retrieved, so no MASLD statement follows. The Visium within-slide check
(PRESPEC_C section 7) was not run: only ~4 µm/px preview images exist locally and the
model needs 0.5 µm/px. Source:
`Analysis/MASLD_Model_Benchmark/executions/he-molecular-c6-tcga-20261003T133755Z/`
(jobs 1600-1605).

### Mouse diet models against the human program stage trend, 2026-10-03 (not adopted)

Candidate analysis, prespecified in `scripts/analysis/preclinical_fidelity/PRESPEC_F.md`
(two amendments, both recorded there). Public mouse data only; the in-house MCD
experiment and all pooled-screen data are excluded. Output: `Analysis/MASLD_Model_Benchmark/executions/preclinical-fidelity-20261003T153626Z/`
(`f2_fidelity_r2` is the corrected statistics run; `f2_fidelity` is kept).

Human reference: OLS of each fixed program score on fibrosis stage and sex in GSE130970
(76), GSE135251 (214), GSE162694 (109) and GSE213621 (361) participant libraries,
fixed-effect meta-analysis. 92 of 117 programs are testable in all four cohorts and 62
are stage-associated at BH q<0.05 over those 92. Leave-one-cohort-out Spearman over
programs is 0.797 to 0.828 (median 0.801). Some libraries are from the same person, so
human standard errors are slightly too small.

Mouse: 17 diet-versus-lean-control contrasts from 15 GEO series; animals are the unit.
GSE196908 (GAN, RPKM table of 9,809 genes) has no testable program and is untestable.
In the other 16 contrasts, 37 to 56 programs are testable on both sides.

| Result | Value |
|---|---|
| Primary endpoint: same-diet-class minus different-class similarity of mouse program profiles | −0.012 (0.381 vs 0.393; 33 and 85 series-independent pairs); class-label permutation p=0.479; not met |
| Class-balanced version | −0.030; p=0.582 |
| Spearman with the human stage trend, median over 16 contrasts | 0.726 (human ceiling 0.801) |
| Contrasts passing the animal-label shuffle null at BH q<0.05 (family of 16) | 2: FPC 20 wk, GSE162876 (0.785, q=0.032); LIDPAD, GSE159911 (0.806, q=0.032) |
| Contrasts at q=0.0502 | 5: MCD GSE156918 (0.729), CDAA GSE294630 (0.764), WD GSE220575 (0.713), FFC GSE292565 (0.753), GAN 24 wk GSE294629 (0.802) |
| Contrasts with no agreement | HFD GSE274914 (−0.450, 95% interval −0.662 to 0.088); WD 8 wk GSE305484 (−0.389, −0.632 to 0.054) |
| Stage-associated human programs testable in at least one mouse contrast | 35 of 62; the other 27 are untestable in mouse, not negative |
| Programs significant in the opposite direction in at least one contrast | 6 |

Reading: mouse diet models of different classes change the testable human programs in
a similar pattern, and that pattern agrees with the human fibrosis-stage trend in most
series. These data do not show that diet class determines which programs are
reproduced. Agreement is a similarity of direction in bulk liver; composition and
within-cell change are not separated in either species, and diet class is confounded
with laboratory. An independent recomputation (statsmodels OLS, pandas rank
correlation) reproduces the human coefficients, the 16 agreement values and the class
difference to 1e-15.

### In-silico saturation over measured liver regions (P4), 2026-09-17 (not adopted)

Every possible substitution was scored across a region universe defined on
2026-09-09, before any score was retrieved: **96,460** H3K27ac regions of the
released chromatin predictor (U1) and **6,930** snATAC consensus peaks selected
by disease-change evidence, cohort significance, genetic overlap or program
linkage (U2). The reduction over all **103,390** regions completed with
**zero failures**. The analysis prespecification was written on 2026-09-14,
with retrieval at 71,000 of 103,390 regions and no region's scores opened
outside a 200-region mechanics pilot. One amendment was registered before any
statistic: the primary sensitivity scale is the raw predicted effect relative to
the track's own mean level, not the Atlas calibrated quantile, because inside a
measured liver peak the quantile saturates near its ceiling. Both scales are
reported throughout, and the quantile version is labelled saturating.

Two of the five written predictions held.

**F1, the primary question, returned a null, as predicted.** The H3K27ac regions
the released predictor can reliably predict (**4,335** matched of 4,364) are not
more sequence-sensitive than matched regions it cannot: median relative
sensitivity **0.086844** against a matched-null median of **0.086815**, a
difference of **+0.0000296**, p **0.97**, BH q **0.97**. Prediction P2 also held,
with the correlation between a region's per-fold predictor skill and its relative
sensitivity at **0.027** within strata, against a registered bar of 0.15. The
secondary DNase channels do move (**+0.005026** for adult-verified liver DNase,
q **0.0059**), which is reported as it stands and is not read back onto the
predictor's mechanism. Sixteen of the 24 F1 cells are computed; the remaining
eight are life-stage groups with no eligible track and are reported as
indeterminate, not as null results.

**F2 failed its prediction on the primary scale, and the two scales disagree in
sign.** Promoter regions (**13,217**) were predicted to be more localised. On the
relative scale their share of sensitivity in the best 50 bp is **0.1700** against
a matched null of **0.2087**, a difference of **−0.0387**: they are less
localised, not more. On the saturating quantile scale the same comparison gives
**+0.00095**, the predicted direction. Every F2 p sits at the resolution floor of
**1/(10,000+1)**, so the family is saturated and BH cannot separate its members;
the correct statement is that these differences are at the floor of 10,000
matched draws, not a p value. The prespecification anticipated the scale
question, which is why the disagreement is reportable rather than a choice made
after the fact.

**F3 held on which factors dominate and failed on whether composition differs.**
HNF4A, RXRA, SP1, FOXA2 and the cohesin/CTCF tracks do lead the top-5 lists of
liver regions, as predicted. But **20 of 218** tested TF tracks survive BH at
0.05, against a prediction that none would. Thirteen are less frequent in
predictable regions and seven more frequent. The largest is adult liver HNF4A,
in the top 5 of **91.9%** of predictable regions against **86.7%** of matched
controls (q **2.8e-21**); POLR2A (**3.0%** against **4.6%**), MAX and GABPA run
the other way. A reading consistent with that split, offered as interpretation
rather than result, is that predictable regions lean toward lineage-factor
occupancy and away from general promoter machinery [Inference]. A further 120
tracks did not reach the registered 50-occurrence rule and were not tested.

**F4 is computed but is one comparison reported twice.** Thirty-two of 48
contrasts are computed and 16 are indeterminate on the minimum-group rule. All
32 are positive and at the draw floor: snATAC peaks carrying a disease-change
evidence state are more sequence-sensitive than matched indeterminate peaks,
most in GSE244832 hepatocyte ATAC (**+0.0852**) and least in H3K27ac
(**+0.0212**). The `original_nonindeterminate` and `source_dependent` groups
return near-identical statistics because **1,844 of 1,848** such peaks are
`source_dependent`; they are the same peaks under two labels and must not be
counted as two findings. The prespecification named this underpowering before
the family ran.

**P5 failed in the opposite direction.** Relative sensitivity was predicted to
correlate positively with a region's measured H3K27ac signal mean at Spearman
above 0.2, on the reasoning that stronger regions carry more motif content. The
observed correlation is **−0.136** for H3K27ac, **−0.070** for ATAC and
**+0.016** for DNase over all 96,460 regions. The stated mechanism is not
supported.

**F5 now has its background, and the family is reported.** Its size-matched
random gene-set null needs **18,570** promoter regions that were never in the
2026-09-09 universe. Those regions were defined from source promoter annotation
alone, with the universe receipt written and checksummed before any score was
requested, and retrieved from the hosted Atlas over **17 h 16 m** into a
**204.8 GB** archive; the retrieved set equals the originally missing set exactly,
with no region on either side only. Every program's coverage gap is closed: the
missing-region count is **0** for all 117 programs, against 18,570 before. The
number of testable programs is unchanged at **11 of 117**, because the other 106
fall below the registered 200-region minimum by size rather than by coverage; the
largest indeterminate program has **193** regions and the smallest computed one
**204**. Retrieval bought an admissible null, not more programs.

Across those 11 programs, **113 of 385** TF contrasts reach BH q below 0.05
within their own program, 65 positive and 48 negative, with **72** of the nominal
p values at the resolution floor of **1/(1,000+1)**.

**A registered sensitivity shows the result is not promoter base composition,
against the written prediction.** F5's prespecified null matches gene sets on
size alone, while F1 through F4 match on GC, width, promoter status and measured
signal. The unstratified hits depleted promoter-general factors and enriched
liver-identity factors, which is what a CpG content difference produces, so a
GC-matched variant was registered on 2026-09-18 with its predictions written
before any matched score was computed. The confound is real: **6 of 11** programs
have a promoter-region mean GC outside the central 95% of the size-matched null,
five of them at the draw floor, in both directions (the two hepatocyte programs
at **0.656** and **0.667** against a null of **0.631**, the cholangiocyte,
macrophage and T-cell programs as low as **0.584**), and matching brought the
null's GC closer to each program's in every case. But the written prediction that
fewer than half the contrasts would survive **failed**: **94 of 113** survive
under a GC-decile-matched gene-set null with **zero sign flips**, and only 19 are
lost, of which **16** are promoter-general factors and **1** is liver-identity.
F5 is therefore reported as robust to promoter GC composition.

One bound belongs beside it. The surviving contrasts concentrate in the
non-hepatocyte programs (**47 of 140** macrophage tests, **27 of 70**
cholangiocyte, **21 of 70** T-cell) while the two hepatocyte programs give
**4 of 70** and the single fibroblast program **0 of 35**, against a TF panel
drawn from liver and hepatocyte ChIP. Much of what this family finds is that
a non-hepatocyte program's promoters carry a different liver-TF sensitivity
profile than random liver-expressed genes, which is closer to restating the
programs' cell identity than to a statement about what the sequence model learned
[Inference].

**F6, the genetic-overlap family, ran as a new registered family and failed three
of its four predictions in the same direction.** The 2026-09-09 universe left this
family empty because the variant-peak overlap table did not exist when it was
defined, eight hours earlier. That universe was not reopened. F6 is a separate
family over separate regions, with its own prespecification written on 2026-09-17
and checksummed before any of its scores was requested, over the **6,920** peak
coordinates carrying a queried variant, of which **5,686** were new and retrieved.
Admission used lineage and peak coordinate only; the overlap table's hosted
variant-effect columns were never read, because selecting regions on a model score
and then testing a model score is circular. Controls are promoter background and
original snATAC peaks in the same lineage with no overlap, matched on GC decile,
width decile, promoter status and signal mean and sd quintiles, with 10,000 draws.

**Within those strata, peaks that overlap a queried variant are less
sequence-sensitive than matched peaks that do not.** All 48 contrasts computed
and none were indeterminate, but 12 of the 48 are the `liver_atac_adult_verified`
track group returning values identical to `liver_atac` in every row, so they are
one comparison reported twice and are counted once here, as in F4. Of the **36**
distinct contrasts, **24** reach BH q below 0.05 and **every one of the 24 is
negative**; none is positive. The largest is
GSE244832 hepatocyte ATAC at **0.2292** against a matched null of **0.2487**, a
difference of **−0.0195** at the resolution floor of 1/(10,000+1), and it
replicates in GSE281367 at **−0.0132** and in stellate cells at **−0.0098** and
**−0.0042**.

The qualifier "within those strata" is load-bearing, and an independent
reconstruction is what established it. Before matching, the same comparison runs
the other way: **+0.0474** for GSE244832 hepatocyte ATAC, and matching reverses
the sign in **7 of 12** reconstructed contrasts, including all six hepatocyte
ones. The reason is that F6 borrowed F5's promoter-defined background, so the raw
control pool is **96.2%** promoter against the group's **64.0%**; after matching
those stand at 0.640 against 0.640 and GC at 0.5700 against 0.5693. A naive
comparison of these two sets would have reported the opposite result. The
reconstruction re-derived GC, width, promoter status and donor signal from their
sources and rebuilt the statistic on its own region pool and its own decile
edges; it agrees in sign on every contrast, with a largest gap of **5.7e-4**,
about a tenth of the median effect.

P6a predicted no difference. Its first half holds, since overlap peaks are not
more sensitive, and its second half fails: the difference sits well outside the
matched null's central 95% interval, in the less-sensitive direction. Decision
rule D2 named only the more-sensitive failure, so this direction was not
anticipated. P6b holds as an ordering and fails as an exclusion: ATAC moves most
and DNase next, but H3K27ac does move in hepatocytes (**−0.0029**, **−0.0019**)
while staying null in stellate cells (**−0.0003**, q 0.50). P6c fails; the
queried subset is not indistinguishable from the unqueried subset, with 6 of 12
contrasts below q 0.05 and every significant one negative. P6d fails, and not on
the factor it named: HNF4A did not move, and **11 of 126** TF contrasts clear BH,
led by SP1 depleted in all four cohort-lineage cells (**−0.071** stellate,
**−0.037** and **−0.016** hepatocyte) with GABPA, CTCF and NR2F2 enriched.

A reading consistent with the direction, offered as interpretation and not tested
here, is that a region whose predicted output is highly sensitive to single-base
substitution is a region where variation is less tolerated, so the peaks carrying
a common variant are the less sensitive ones [Inference]. Two limits bound it.
The variant list is liver-trait fine-mapping output, not a random sample of common
variation, so ascertainment follows association power and gene density rather than
sequence. And whether a peak was tested at all depends on that same list, which
the matching strata do not control.

One limit travels with every number above: each value is a predicted effect for a
substitution that mostly exists in no person, and is never presented as measured
chromatin or as a variant effect in a donor. The earlier qualifier that the
archives were present but not integrity-validated no longer applies. All
**103,390** original regions, the **5,686** genetic-overlap regions and the
**18,570** background regions were opened and checked: every payload carries its
quantiles layer, the track count its own request recorded, and a finite
in-range quantile sample. Exactly two of the 18,570 carry a reference gap base,
and both archive 1,497 substitution rows rather than 1,500, which is the correct
count: they are overlapping 500-bp windows covering a single N at
chr10:58267943, and an N position carries no substitution. The check's original
3 x interval-bp expectation was wrong for those two and now reads the
substitutable base count from the reference instead. Nothing here
nominates a gene, a factor or a region for follow-up. Sources:
`data/p4-full-21772842/families/`, `data/f5-f6-families-21797329/`,
`data/f5g-gc-sensitivity-21797337/`, `data/f5bg-integrity-21797404/` and
`review/f6-overlap-21797374/` under the campaign; prespecifications
`scripts/analysis/alphagenome_atlas/43_saturation_prespec.json`,
`scripts/analysis/alphagenome_campaign/data_p4_u2c_prespec.json` and
`scripts/analysis/alphagenome_campaign/data_p4_f5_gc_prespec.json`.

## Cross-assay modules: v2 validated candidate, 2026-09-12 (not adopted)

A second, smaller module family built so that every assay has a real test to
fail. The gene universe is fixed first to symbols measured in **all four** of
five-cohort bulk RNA, the single-nucleus atlas, both Visium sources and the
PXD051911 liver proteome (**4,468 genes**); modules are built inside it from
co-expression alone, fixed under a recorded checksum, and only afterwards
labelled with the existing five-cohort stage results. The 117 prespecified
Hotspot programs, their weights and their testability rules are unchanged;
the largest Jaccard between any new module and any of the 117 is 0.15 (median
0.03), which is an overlap statistic and not evidence of a new coordinate
system. Producer `scripts/analysis/cross_assay_modules/`; candidate panels
`figures/candidates/cross-assay-modules-2026-09-12-v2/`; run root
`RNA-seq/results/cross_assay_modules/cam-v2-20260912T130258Z` (figures
`figure_v4/`); 63 validation checks passed. A second review the same day
found seven labelling and presentation defects (`GRAPH_K_AMENDMENT.md`,
amendment 5); their corrections are in the numbers below and change no fitted
value.

**This supersedes the 2026-09-11 candidate, whose numbers are retracted.**
That panel applied the discovery direction both at scoring and in a one-sided
test, so for the 101 of 160 modules whose direction is negative, agreement was
tested in the wrong tail: a reversal could print as support and an agreement as
a negative. It also defined negatives by power rather than by a bound, left the
competitive test uncorrected, tested paired couplings for modules failing the
chromatin coverage rule, and reported an ATAC estimate unadjusted for the
processing well. Module construction was not implicated and is unchanged: the
v2 membership reproduces the v1 checksum byte for byte. Corrections are listed
in `scripts/analysis/cross_assay_modules/GRAPH_K_AMENDMENT.md`, amendment 4,
and are not outcome-blind, since v1 had already read every outcome.

| Quantity | Value | Interpretation |
|---|---:|---|
| Shared universe | 4,468 genes | Measured in all four reference assays |
| Modules | 160 | Mean 22.6 genes; 3,611 genes assigned; 94 discovery-supported |
| Discovery participants | 760 fibrosis, 407 activity | Four and three of the five cohorts carry the grade |
| Modules meeting the coverage rule, liver proteome | 160 / 160 | Compare 26 / 117 for the prespecified programs |
| Modules meeting the coverage rule, CosMx and plasma | 0 / 160 | Coverage-only columns |
| Module-by-assay cells | 1,600 | 160 modules x 10 assays |
| `replicates` cells | 144, over 87 modules | Two-sided BH q<0.05 with the oriented slope positive; 6 of them are in the single-nucleus column, a coarse-diagnosis proxy rather than a fibrosis stage |
| `discordant` cells | 8, over 8 modules | Reversals; five in H3K27ac. v1 could not print this state |
| `tested_negative` cells | 112 | BH-corrected bound excludes a standardized slope of 0.2 (230 at 0.3) |
| Cells exceeding co-expressed sets | 1 | `cam_m70` in GSE268273, after BH over the family |

**Associations replicate; one test exceeds co-expressed sets.** On the fibrosis axis 102 of
160 modules pass joint Benjamini-Hochberg correction over the 320 module-by-axis
tests and 105 beat a permuted label; on the activity axis 93 and 97. Against
10,000 size-matched connected subgraphs of the same co-expression graph,
**zero** modules survive correction on either discovery axis, and one cell
survives in the whole assay panel. Per assay, 2 to 10 modules clear an
uncorrected competitive p<0.05 against 5.5 to 8 expected. The competitive
machinery is calibrated: run under 20 within-stratum label permutations of the
discovery fibrosis axis it rejects at 0.056; that shows the reference does not
reject spuriously and says nothing about its power to detect specificity. The
supportable statement is that one module-by-assay test exceeds the
connected-set reference after the per-assay correction; a module can reproduce
its histologic association in an independent cohort while being no better than
other co-expressed sets of its size, and the two are reported as separate
states. Proteome cells replicate only if they also hold under ten
intensity-complexity descriptors at BH q<0.05 over the family (31 of 34
nominal survivors).

**Negatives now carry a bound.** 112 cells exclude a standardized slope of 0.2
in the discovery direction after correction, with medians of the one-sided 95
percent upper bound from -0.04 to 0.06 by assay. The ATAC column resolves
nothing at all: with 39 donors in eight wells, a well-adjusted estimate and a
within-well permutation, the median upper bound is 0.37, so all 110 testable
cells are unresolved rather than negative.

**Within-participant coupling is the strongest column, and its primary
adjustment is the well.** In GSE296875, where RNA and ATAC come from the same
nuclei, 33 of 150 modules testable in both assays couple at BH q<0.05 given the
processing well, all positive, largest rho 0.81; severity is not in that
adjustment. Two sensitivities carry their own within-well permutation p and BH
correction: given well and the binary fibrosis label, all 33 remain at q<0.05;
given well and the donor's hepatocyte nuclear fraction, the composition common
cause of both scores, 19 of 33 remain at q<0.05 although a median 96 percent of
the effect size is retained and 31 of 33 keep at least half. In GSE267145, 8 of
135 couple given steatosis, ballooning, lobular inflammation, fibrosis and sex,
largest |rho| 0.49, two of them inverse; that is the arm with histology in the
primary adjustment. Neither arm compares a module's coupling with that of
matched co-expressed sets, so neither establishes module-specific coupling or
shared regulation.

**Bulk severity signal is substantially cell composition.** Composition-adjusted
slopes are a median 41 percent of unadjusted on the fibrosis axis, where 73 of
102 remain significant, and 17 percent on the activity axis, where 4 of 93
remain. No bulk or proteome cell in this panel is independent of composition.
The single-nucleus per-lineage arm is within-lineage donor-level evidence, not
composition-free: 27 of 1,280 module-by-lineage tests pass correction.

**Boundaries.** Module membership is not reproducible under cohort removal
(adjusted Rand 0.28 to 0.52), while the module score evaluated **on the omitted
cohort** is (median Spearman 0.939 to 0.965, against 0.21 to 0.46 for random
same-size sets and 0.42 to 0.80 for the held-out first principal component). The
family is a set of scoring coordinates with stable values and unstable rosters:
no module's gene list is a finding, and no member gene inherits its module's
state. Fourteen modules are marked score-unstable and are excluded from the main
panel rows. Modules without discovery support were not selected to be null and
are an unselected comparison class. The two Visium columns answer spatial
coherence rather than stage association: 18 of 160 modules show spatial excess
over matched genes in GSE192741, built within physical section and collapsed to
four donors; the other 142 are unresolved, not negative, since the test carries
no bound on spatial excess and 90 of them sit above their null mean. The Vu
source cannot resolve donors and is source-dependent throughout. The
single-nucleus column tests a coarse diagnosis category and is a labelled proxy
for the fibrosis family, not a fibrosis stage. Eight modules include one member
whose proteome measurement is a multi-gene protein group; none of them
replicates in the proteome, and the 160-of-160 proteome coverage count includes
them. Two of the six prespecified Figure 5K examples (`cam_m147`, `cam_m06`)
fail the score-stability rule and are marked as such rather than replaced. This
is a validated candidate, not an adopted release, and it may not redefine any
evidence class, hero gene, threshold or figure.

## Public functional stress tests

| Workstream | Result for this paper | Disposition |
|---|---|---|
| Myojin palmitate survival | Did not confirm the proposed Resource-level enrichment | Complete supplement; survival is not interchangeable with lipid storage |
| Acute induction/localization/reversal | Strict multi-assay adoption criteria did not pass | Complete supplement with valid nulls and source checks |
| Risk/state double dissociation | Full two-arm criteria did not pass | Complete supplement; no functional-partition headline |
| Chronic reversal/risk bridge | Full risk-to-state architecture criteria did not pass | Complete supplement; no causal relay claim |
| Multicellular response topology | Did not produce a source-independent main-paper mechanism | Complete supplement; no replacement search |

These analyses strengthen the Resource by defining what acute public models can
and cannot validate. They do not justify reopening mechanism hunting or moving
the paper away from its Resource identity.

## N4 within-person agreement eligibility correction, 2026-10-14

Job 22236285 completed the retrospective correction on historical fragment
counts with the 117 fixed programs unchanged. The inferential grouping remains
the source-declared pair/genetic cluster, not a verified biological person.
All 164 GSE193066 libraries remain identity-unqualified; N1 non-mixture status
does not establish person identity. Primary replication, ICC and derived ceiling
claims using that cohort remain excluded.

The 22 different-visit pairs (3 GSE130970, 19 GSE193066) retain stage ICC
0.073529 (95% person-group bootstrap interval -0.367460 to 0.480198), median
program ICC 0.705732 (0.384154 to 0.840625), and best-program activated stellate
(PDGFRA) ICC 0.905771 (0.769980 to 0.944250). Of 108 testable programs, 102
point estimates exceed stage ICC. These are retrospective descriptions;
intervals use 20,000 grouped bootstrap draws, seed 20261007. The historical
October 7 claim of **69 stage-BH calls is withdrawn**: raw bootstrap tail
fractions were descriptive diagnostics, not calibrated null-test p-values.

The historical **87 random-null calls are withdrawn**. The corrected random
pool applies the same N1 eligibility as observed GSE193066 libraries, requires
called genotype components, excludes observed source/component units and
conflicting source/component duplicates, and matches cohort and exact stage
without replacement or stage substitutions. Only 27 eligible GSE193066 null
units remain. Stage 2 requires 12 endpoints but has 8, stage 3 requires 19 but
has 2, and stage 4 requires 1 but has none. Consequently the different-visit
random-null comparison is **untestable**, with no null p-values or BH calls.

The separate 28 untyped pairs (3 GSE135251, 25 GSE213621) retain stage ICC
0.636545 (0.324318 to 0.809406) and median program ICC 0.108216 (-0.141212 to
0.310574). Their exact-stage eligible null is estimable: **zero programs pass
BH q < 0.05**, with all 117 fixed programs retained in the hypothesis family
and untestable programs assigned p = 1. The 2,000 null draws use seed 20261008;
GSE213621 alone (25 pairs) also gives zero null-BH calls. Pair type remains
unestablished, so this comparison does not resolve whether the repeated
libraries are repeat biopsies or duplicated/mislabeled samples.

Sources: [corrected summary](../Analysis/Hypotheses/N4-eligibility-corrected-20261014T140730Z/summary_by_group.tsv),
[exact-stage capacity](../Analysis/Hypotheses/N4-eligibility-corrected-20261014T140730Z/null_capacity_1.tsv),
and [identity disposition](../Analysis/Hypotheses/N4-eligibility-corrected-20261014T140730Z/GSE193066_identity_disposition.tsv).
Focused exact-stage, no-replacement, deterministic-seed and insufficient-pool
checks passed on compute; observed ICCs and intervals reproduce the historical
outputs, while invalid stage-BH and different-visit null calls are absent.
No additional N4 job is pending. These outputs do not qualify GSE193066 or
substitute historical counts for the adopted primary bulk refit.

## Historical corrections that must not return

| Retired statement | Correct interpretation |
|---|---|
| 1,277/1,259 human samples as current | Pre-recount documentation; corrected source registry is 1,281 count-matched and 1,257–1,260 QC-passing depending on release |
| 54 disease-significant Hotspot modules | Run-level pseudoreplication; donor-level result was 22 in that historical analysis, while the present prespecified selection uses the Plan-20 criteria above |
| 58 Andrews donors with F-stage | 58 sequencing runs mapping to 46 biological donors |
| cross-cohort inferred F-stage accuracy 0.638 | Run-level/leaked result; not authorized for cross-cohort fibrosis inference |
| Fixed 13 historically nominated ligand–receptor pairs | Source-qualified screen: six testable (all adjusted q=1), seven untestable; heatmap is descriptive |
| 290 sex-dimorphic genes | Retired random-slope result; fixed-cohort model yields 8 and the former female-progressor narrative is retired |
| k=6 discrete patient subtypes | Factor subspace is stable, but hard partitions are not; use continuous axes |

## Where the full history lives

## N3: source-qualified hepatocyte target-expression check

The October 14 refit supersedes the October 9/13 N3 estimates. The original
four-stage model retains 70 biological donor records from four cohorts,
adjusting for dataset and inferred sex. Only GSE174748 SRR14602057/58 change
from Steatosis to source-documented Cirrhosis; stage counts are 14 Healthy,
11 Steatosis, 34 Steatohepatitis and 11 Cirrhosis. All native count matrices,
the fixed 51-target roster, controls and 1,000 within-dataset permutations per
fit (seed 20260923) remain unchanged. Expression eligibility gives separate
45-target raw and 40-target corrected BH families.

Ordinary-test calibration fails in both fits: mean nominal-positive fractions
are 0.078889 raw and 0.078125 corrected, above the prespecified 0.07 limit.
The prespecified empirical-permutation BH fallback gives nine raw calls and
**zero corrected calls**; the ordinary-BH counts of 12 and six are diagnostic.
Corrected LEPR beta=0.544008, ordinary BH q=0.014856 and primary empirical BH
q=0.074212, so the earlier “LEPR robust” conclusion is withdrawn. No detected
corrected calls does not establish absent effects, powered negatives or
clinical target effects. Both fixed positive-control composites pass.

All 208 included libraries were corrected with native empty-droplet profiles;
none of the 70 included donor records has passthrough. The historical
nine-target ambient diagnostic roster remains fixed separately. Its native
attenuation labels use ordinary BH and remain diagnostic, not primary calls.
The nominal contamination-stage coefficient is 0.177292 logit per step
(95% CI 0.049382–0.305203; p=0.007348). Incomplete ambient removal, no doublet
analysis and no second correction method remain limitations. Both complete
fits and the repaired summary pass the thirteen checks in
[verification.json](../Analysis/Hypotheses/N3-source-diagnosis-20261014/verification.json).

The previous 2,220-line number ledger is preserved verbatim in the documentation
archive. Use [technical/NUMBERS_HISTORY.md](technical/NUMBERS_HISTORY.md) to find
historical thresholds, superseded releases, and source files. New prose and
figures should cite this compact ledger unless a technical audit explicitly
requires the archive.

## Review triage 2026-10-21 (release 2026-10-21-r9)

Run: `Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/` (PLAN.md, DECISIONS.md, analyses/AN1-AN7/RESULT.md, installed.tsv). All genetic values are on COLOC r2 (`susie_coloc_r2_20260924T053244Z`).

- **Shared-posterior check of carried negatives (AN1; sensitivity, not a new COLOC release).** The reimplementation reproduced all 816 rechecked rows exactly. Of 12,087 carried primary negatives, 343 fail on every signal pair, 1,744 on some and 10,000 on none. Under the user's rule D1 (tested if any tier-1/2 row is eligible): 3,776 → 3,736 genes, 26 shared vs 21.1 expected, MH OR 1.30 (0.84-2.02), p 0.27. Under the adopted priority rule: 3,589 genes, OR 1.28 (0.83-1.99), p 0.28. 3,769 stays prespecified. PP.H3 counts: 3,918 → 3,877; 3,335 → 3,306 (3,101-3,416).
- **Figure 2D at genome-wide leads (AN2; user decision D2).** Direct 29.65% vs enzyme 6.67% coding/splice posterior; difference 0.2298 (0.0808-0.3958); ratio 4.45; 16 direct and 115 enzyme clusters. Previous all-lead values 23.34% vs 3.87% (ratio 6.03) are superseded. Non-European reliable sets: 6,227.
- **Figure 3C (AN3; user decision D3).** Five cohorts with the GSE126848 filterByExpr refit (3,440 calls, BH over 15,592 genes): union 12,180, mixed sign 1,211 of 3,839 (31.5%).
- **H82 joint library adjustment (AN3; user decision D4, post hoc descriptive deviation).** F1 99 of 1,338; F2 2,113 of 4,266; NAS1-2 2,595 of 5,338.
- **Spatial family (AN4; descriptive, not prespecified).** With the all-16-adjusted matched-gene null, 68 of 114 (GSE192741) and 74 of 114 (Vu) programs exceed it at BH q < 0.05; Stromal ECM ranks 42 and 29.
- **p12 prior (AN5).** At p12 = 1e-6, 122 of 436 supported genes fall to PP.H4 ≤ 0.5 (111 anchored, 8 both, 3 untested); 29 gain at 1e-5. Three-way genes: 108, or 100 with r² ≥ 0.8 between caQTL and eQTL leads. Of 14 genes at PP.H4 = 1.0, 11 fall to ≤ 0.5 after removing GWAS variants inconsistent with the LD panel (post hoc filter; no class change).
- **GTEx area adjustment (AN6).** Fat rho 0.631 vs 0.638; paired difference -0.0065 (-0.033 to 0.004). The ancestry Holm p 0.0199 is not at its null's floor (2/201 and 2/101 before Holm).
- **H59 sensitivities (AN7).** QC-passing-cell donor scores: 0 of 117 selected (smallest q 0.304). T-cell 20-cell gate: 0 of 11.
