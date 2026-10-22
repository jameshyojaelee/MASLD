# Continuum and sliding-window analysis: code and scientific review

Reviewed 6 September 2026. This is an independent technical review, not a change to manuscript claims or an adopted analysis.

## Recommendation

Retain the continuum as a descriptive molecular coordinate and retain the sample-level, stage-adjusted association analyses with their existing boundaries. Do not use the current coherence N1/N2 results as confirmatory evidence of changing gene coordination. Correct N1’s axis reconstruction first, then calibrate the composition-conditional test. The recent matched-gene results do not establish that this particular signature outperforms comparable gene sets.

The main sliding-window implementation does not treat overlapping windows as independent replicates. The most consequential defects are in the newer coherence calibration and a transfer-support flag, rather than in the basic window means.

## Scope and evidence

I reviewed recent working files regardless of Git state, following the current analysis READMEs and paper authorities. “Not including the model” is interpreted as excluding the predictive-model benchmark, model architectures/training, random-forest benchmark internals, and model-performance claims. Statistical regressions, null construction, and the PCA/projection operations defining the continuum are included.

The review traced score construction in 02/03, transcript methods in 10, pathway scoring in 20, Hotspot analysis in 31, paired contrasts in 40, windows in 63/64, system aggregation/inference in 67/70, stress tests in 73, the recent coherence scripts 80–85, and matched-null/control-free-transfer scripts. This is not an exhaustive review of legacy code, TF/NMF internals, or every figure.

Evidence includes source inspection, stored tables, both existing test suites run on compute nodes, a fresh reconstruction of the frozen and N1 axes, and a seeded null simulation. S4U was inspected by PDF text extraction and rendered PNG. Large production analyses and their 2,000–10,000-draw nulls were not rerun.

## Continuum versus sliding windows

| Object | Calculation | Question |
|---|---|---|
| Continuum | One molecular score per sample, from cohort-fitted signature PC1 or externally fixed discovery loadings | Where does this sample sit along the expression direction? |
| Sliding windows | Mean molecular scores among samples within a fixed continuum-percentile interval | How does the average profile vary along the ordering? |
| Coherence | Member-gene correlation summarized by leading eigenvalue fraction or mean Fisher-transformed correlation | Do genes vary together differently at different positions? |

The published signature contains 145 genes; 139 occur in the Resource substrate, with additional cohort-level coverage checks. The fixed projection uses normalized log counts and discovery-oriented loadings. The cohort-PC1 sign is aligned to that external direction, rather than selected using target histology. The mean of the two percentile scores is presentation-only.

Nine windows are centered at 10%, 20%, …, 90%, with width 20 percentile points. The first spans approximately 0–20%, the second 10–30%, and the last 80–100%. Neighboring windows overlap by half their nominal width. Actual membership depends on ranks, ties, and endpoint rules.

These are within-cohort relative ranks, not a common fibrosis unit or time. A 90th-percentile sample in one cohort need not be more severe than a 70th-percentile sample in another. Smoothness can result partly from overlapping samples.

The window mean is an unadjusted expression summary. The main continuum coefficient instead comes from a sample-level regression adjusted for categorical fibrosis stage and inferred sex. A rising curve alone therefore does not establish within-stage association, a transition, gene-coordination change, or a treatment window.

Sources: [projection][projection], [window intervals][windows], [window means][window_means].

## Findings requiring correction or narrower interpretation

### 1. High: N1 simulates a different axis from the observed coherence contrast

**Verified implementation defect.** N1 uses the frozen registry for the observed ordering, but applies discovery weights directly to gene-standardized log-CPM for simulated orderings. The frozen continuum instead uses quantile-normalized log2(count+1), without that gene-wise standardization. These are different measurement scales and transforms.

Fresh reconstruction on real data:

| Cohort | Samples | Frozen projection vs registry, Spearman | N1 reconstructed axis vs registry, Spearman | Different tercile assignments |
|---|---:|---:|---:|---:|
| GSE126848 | 53 | 1.000000 | 0.978874 | 8 |
| GSE130970 | 76 | 1.000000 | 0.988954 | 6 |
| GSE135251 | 214 | 1.000000 | 0.987891 | 14 |
| GSE162694 | 140 | 1.000000 | 0.979970 | 12 |
| GSE213621 | 361 | 1.000000 | 0.979693 | 38 |

High correlation does not repair a null intended to reproduce the exact assignment to coherence groups. Sample-space rotation preserves the covariance of the matrix being rotated, but this does not reproduce the joint relationship between the actual frozen score and expression outcomes.

**Consequence:** current N1 artifact-floor P values cannot be treated as calibrated for the frozen continuum. This does not prove the biological contrast false or predict how much corrected P values would change.

**Correction:** construct a joint matrix containing the actual score-generating signature measurements and the expression outcomes; rotate both together and apply the original projection to its own measurement block. Require exact reproduction of observed scores/ranks before null generation. If changing scales, transform loadings consistently. Also state the rotational null’s distributional assumptions: covariance preservation alone is not a general null for arbitrary non-Gaussian expression.

Sources: [N1 transform and score reconstruction][N1], [actual normalization][normalization], [fresh diagnostic results][fresh].

### 2. High: N2 does not establish composition-independent coherence

**Verified implementation; calibration unverified.** N2 fits a linear Gaussian model for the rank-normalized continuum conditional on selected lineage log ratios, using the same samples. It draws independent Gaussian residuals with one fitted SD. The script supplies no conditional-distribution calibration, heteroscedasticity checks, or independently fitted conditional law. In-sample R-squared measures fit, not validity of the resampling distribution.

The selected GSE162694 test includes 140 samples and eight lineage variables, with recorded R-squared 0.766. In the residualized arm, the leading-eigenvalue contrast gives nominal P=0.8996 for Ductular injury and P=0.01949 for Stromal ECM. These are conditional on the resampling assumptions, not evidence of cell-intrinsic ECM coordination or proof that composition causes the ductular result.

Cholangiocytes are omitted for poor observability, correctly disclosed by the source. Their omission limits how strongly a ductular program can be cleared of composition effects. The selected subcomposition also does not condition on every possible lineage change.

**Correction:** describe this as sensitivity to measured lineage proportions. Validate the conditional law and empirical type-I error for this exact coherence statistic before interpreting its tail probability inferentially. Independent estimation/cross-fitting may help but does not itself solve misspecification. Null compatibility is not evidence of a specific compositional mechanism.

Sources: [N2 fitting and draws][N2], [stored N2 results][n2_results]. Conditional-randomization validity depends on the resampling law; robustness results require additional assumptions: [Niu et al.](https://arxiv.org/abs/2211.14698).

### 3. Moderate: the coherence pilot is narrower than the recorded specification

**Verified mismatch, not proof of an unauthorized change.** The specification describes a 50-neighbor local-linear estimator across 117 programs, leave-one-cohort-out analysis, and a particular hierarchy of nulls. Scripts 81–85 instead examine two focal programs through bottom-versus-top terciles. N1/N2 calibrate Pearson statistics; only the initial descriptive script also calculates Spearman statistics. The split-permutation script permutes all samples within cohort, rather than stage-by-sex strata or fitted uniform residuals.

The N2 source explicitly selects GSE162694 because it cleared earlier tests. That is an exploratory follow-up on the same observations, not independent replication. Stored N1 residualized leading-eigenvalue contrasts in GSE213621 are -0.03887 for Ductular injury and +0.02278 for ECM, with BH q=0.459 and 0.692. These results do not show replicated positive coherence increases across both independent evaluation cohorts.

Multiplicity differs between experiments. Scripts 82/83 apply BH across ten program-by-cohort tests separately for each arm/statistic. N2 exports nominal P values without a corrected family. “Two focal programs, no correction” in the specification does not justify selecting favorable results across programs, arms, statistics, and cohorts as confirmatory.

**Correction:** label this as a two-program tercile pilot, record the changed estimand, and choose a primary statistic/arm and family before a confirmatory extension. Do not conflate it with the nine descriptive windows or describe the 117-program local analysis as completed.

Sources: [specification][spec], [split permutation][split], [N2 selected-cohort rationale][N2], [stored N1 results][n1_results].

### 4. Moderate: the Fisher-z statistic has an incorrect universal bias correction

**Verified mathematical defect.** The coherence scripts add 1/(n-1) to mean Fisher-transformed off-diagonal correlation and describe this as subtracting a negative null bias. For independent symmetric Gaussian variables, changing one variable’s sign maps r to -r with equal probability. Since atanh is odd, its null expectation is zero. Adding a positive constant creates bias.

A seeded 10,000-draw check at n=30 gave uncorrected mean 0.001881, Monte Carlo SE 0.001929, and “corrected” mean 0.036363.

**Consequence:** absolute coherence values and unequal-size comparisons are shifted. For equal-size top and bottom groups, the constant cancels exactly from their difference. This does not by itself invalidate the current equal-tercile contrasts or their P values, and it does not affect the leading-eigenvalue statistic.

**Correction:** remove the universal adjustment. Any bias treatment must match the actual statistic and sampling assumptions.

Sources: N1 statistic at line 76 of [the N1 script][N1]; [fresh null check][fisher].

### 5. Moderate: a transfer-support flag ignores its accompanying Holm correction

**Verified output inconsistency.** The transfer script calculates Holm P values for its two-cohort primary family but defines F_b_beats_matched_null solely by whether observed rho exceeds the empirical 95th percentile.

The stored GSE240729 row has rho=0.71003, nominal P=0.031968, Holm P=0.063936, and F_b_beats_matched_null=TRUE. GSE174478 has rho=0.51039 and Holm P=0.116883. Neither clears 0.05 in that corrected family.

**Correction:** rename the existing flag “above pointwise null 95th percentile” or require the corrected threshold for a family-level support flag. Keep the FFPE boundary label and competitor comparisons separate. The correlations and stored Holm P values themselves are not wrong.

Sources: [decision code][transfer], [actual transfer rows][transfer_results].

### 6. Scientific boundary: axis association is not signature or system specificity

The recent original-cohort matched-gene test returns nominal P=0.042957 and 0.040959, but Holm P=0.081918 for both. Within-disease sample counts are 109 and 293. This fails the specified superiority criterion; it is not an equivalence test establishing that the signature and random sets are identical.

The system competitive analysis asks whether a system exceeds a same-size comparison set, rather than whether it is associated with the continuum at all. Its recorded result is no system supported on both axes under the connected-graph reference and one under the unconstrained reference.

The connected reference preserves size and connectivity, not uniform sampling over every possible connected set, identical compactness, or exact collection composition. It is a defined sensitivity reference. Calling it “conservative” does not establish a universal ordering of the two nulls.

The two co-primary axes are highly redundant: stored agreement is 0.99805 and 0.99851 in the two evaluation cohorts. Requiring both is a stability condition, not two independent replications. Independent evidence comes from separate cohorts.

**Recommendation:** retain the narrower cross-sectional expression-variation claim. Do not infer a uniquely informative signature, unique biological systems, subtypes, or mechanistic specificity from these counts.

Sources: [matched-gene verdict][matched_results], [competitive null construction][stress], [recorded system stress results][results].

### 7. Scientific boundary: repeat biopsies years apart do not isolate measurement error

The current results discussion uses a within-person SD from biopsies a median two years apart to construct a reliability correction for the stage-versus-continuum comparison. It also acknowledges that this statistic includes real interval change and biopsy sampling variation.

These statements need to remain connected. Variance of paired differences includes both measurement error and genuine change, with covariance terms if those components are dependent. It identifies measurement-error variance only under additional assumptions. Transporting that SD across cohorts adds another assumption, and ordinal staging limits a classical continuous additive-error correction.

**Recommendation:** retain it as a sensitivity analysis under stated assumptions, not a measured correction proving continuum-stage equivalence or superiority. Predictive-model reliability calculations are outside this review.

Source: [current interpretation and caveats][reliability].

## Strengths and limits of the established analysis

The main system model is:

    system_score_z ~ continuum_z + factor(fibrosis_stage) + inferred_sex

It is fitted separately in GSE162694/GSE213621 and combined by fixed-effect meta-analysis. Stored model tables confirm n=109 and n=361, totaling 470 complete cases, rather than the 501 scored evaluation samples. Windows can include unstaged samples; their denominator is not the regression denominator. Cells and windows are not used as independent replicates.

The code checks unique sample/analysis-unit keys. These do not independently prove one biological participant per deposited sample; this review did not reconstruct all original donor crosswalks. The coherence specification itself distinguishes the 844-row sample roster from a verified donor roster. The paired extension uses an explicit 54-person/108-biopsy crosswalk and within-person differences, adjusted for fibrosis change.

The transcript analysis uses RNA-seq-appropriate voom quality weights and limma. Signature components are excluded from ordinary transcript claims and evaluated through leave-one-gene-out axes. Signature genes are removed from downstream program/pathway scores. This removes direct feature reuse, but not biological correlation with the signature or shared cell-composition effects.

The 117-program and 43-system families remain explicit in multiplicity correction, including untestability. System inference retains complete-family BH, seeded within-stage/sex permutations, and maxT results. The newer joint permutation and earlier four-test conjunction have different power and hypotheses; their 35-versus-18 supported-system counts must retain their criterion labels.

The three source-overlap cohorts supply descriptive context rather than primary molecular-association inference. Independent evaluation is two cohorts, not five. The paired molecular analysis asks a separate within-person question, but the source publication already used the Japanese paired cohort; it is not entirely new validation of the axis. [Kamzolas et al.](https://www.nature.com/articles/s42255-026-01543-7) explicitly describes that use.

## Sliding-window and figure assessment

Window membership is left-closed/right-open except for the final upper endpoint. SEs are descriptive within-window summaries, not simultaneous uncertainty bands or uncertainty on a longitudinal trajectory.

The pooled heatmap first summarizes each feature within a window, then takes collection medians, averages collections equally, and finally averages the five cohorts equally. Participant-level inference instead aggregates features within participant before regression. These operations do not generally commute. Script 70 records their reconciliation, making this a documented estimand difference rather than an undisclosed defect.

Equal cohort weight is defensible for an “average cohort” display. It is not a pooled participant mean: the 53-sample cohort and 361-sample cohort receive equal weight. Percentile alignment also aligns relative extremes despite different disease compositions.

I extracted and visually inspected S4U: 43 system rows, nine columns, and 18 row-level asterisks, without a late-minus-early contrast tile. No erroneous printed numerical claim was identified. Its external caption remains essential: asterisks refer to separate sample-level maxT tests, not significant windows or curves. The late-minus-early row ordering is outcome-adaptive presentation, distinct from the outcome-blind system geometry.

Smooth blue-to-red profiles do not establish temporal order or a localized transition. Keep cohort-specific profiles and sample counts alongside the pooled display.

Sources: [window aggregation][system_windows], [participant aggregation and reconciliation][systems], [reviewed PDF][panel].

## Verification and saved evidence

| Check | Result |
|---|---|
| Existing core suite | PASS, six test files; job 21468813 |
| Existing molecular-layer suite | PASS, including system-stress tests; job 21468813 |
| Frozen projection reconstruction | Spearman 1 against registry in all five cohorts |
| N1 reconstruction before rotation | Different ranks/terciles in every cohort |
| Stored system-model denominator | 109 + 361 = 470 |
| Seeded Fisher-z null check | Claimed correction introduces positive bias |
| Window endpoint fixture | Nine windows retained; counts saved |
| S4U | PDF text extracted and PNG visually inspected |
| Full production analysis/nulls | Not rerun |

The first diagnostic attempt stopped on a cohort-level signature row absent from GSE126848’s normalized matrix. The rerun applies the cohort-intersection rule used by production projection. Both existing test suites had already passed and were not repeated. The first batch’s overall success status is not proof every subcommand passed: exit_status.tsv explicitly records the failed diagnostic. Corrected diagnostic job 21468899 completed successfully.

Commands submitted:

    sbatch docs/technical/continuum_review_20260906T121513Z/checks.sbatch
    sbatch docs/technical/continuum_review_20260906T121513Z/recheck.sbatch

The wrappers invoke the existing rnaseq Rscript, both existing test runners, and checks.R. Each requested one CPU, 8 GB, cpu/nslab, with a 90-hour limit. Elapsed times were 31 seconds and 8 seconds. First-job batch MaxRSS was 356,772 KB; the second was reported as zero, which is not interpretable as measured zero usage. No downstream analysis jobs were submitted.

Source inspection used rg, sed/nl and file reads. PDF inspection used pdftotext and pdftoppm. No production code, existing result, figure, or paper authority was modified. The Resource scope validator was not run because this review changes none of those authorities.

New files are confined to this timestamped review directory: REVIEW.md, checks.R, checks.sbatch, recheck.sbatch, their logs, diagnostic TSVs and sessionInfo, and PDF text/PNG inspection artifacts. These support the review rather than replace an analysis.

See [diagnostic script][checks], [core-test log][core], [molecular-test log][molecular], [initial exit statuses][jobs], and the linked diagnostic result tables.

The next concrete correction is N1’s score-generating measurement block and exact observed-score reproduction check. Until the corrected null runs, preserve coherence results as exploratory and retain the narrower cross-sectional continuum interpretation.

[N1]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/83_focal_program_coherence_n1.R:106
[N2]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/85_focal_program_coherence_n2_composition.R:127
[spec]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/80_program_coherence_prespecification.json
[split]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/82_focal_program_coherence_calibration.R:88
[normalization]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/02_build_unsupervised_axes.R:185
[projection]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/03_build_projection_rf_axes.R:253
[windows]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/lib_molecular_layers.R:337
[window_means]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/64_render_five_cohort_context.R:107
[system_windows]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/67_overlay_molecular_systems_atlas.R:538
[systems]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/70_infer_molecular_systems.R:173
[stress]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/molecular_layers/73_stress_test_molecular_systems.R:489
[transfer]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/histology_anchored_continuum/control_free_transfer/01_project_control_free_cohorts.R:581
[transfer_results]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/control_free_transfer/cft-v1-20260828T141625Z/PRIMARY_transfer_verdict.tsv
[matched_results]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/matched_null/mnc-fib-v1-run-20260828T135030Z/VERDICT.txt
[n1_results]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/program_coherence/hac-coherence-n1-20260828T022956Z/focal_coherence_n1_floor.tsv
[n2_results]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/program_coherence/hac-coherence-n2-20260828T024128Z/focal_coherence_n2_composition.tsv
[results]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/RESULTS.md:693
[reliability]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/RESULTS.md:739
[fresh]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/continuum_review_20260906T121513Z/recheck/n1_axis_reconstruction.tsv
[fisher]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/continuum_review_20260906T121513Z/recheck/fisher_z_null_check.tsv
[panel]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/histology-continuum-molecular-systems-s4u-caption-free-20260825T130145Z/panels/s4_molecular_systems_pooled_continuum_windows_with_inference.pdf
[checks]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/continuum_review_20260906T121513Z/checks.R
[core]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/continuum_review_20260906T121513Z/core_tests.log
[molecular]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/continuum_review_20260906T121513Z/molecular_tests.log
[jobs]: /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/continuum_review_20260906T121513Z/exit_status.tsv

