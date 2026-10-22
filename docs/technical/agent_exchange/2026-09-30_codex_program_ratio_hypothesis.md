# Separate hypothesis: count-weighted missing-to-observed program ratios

2026-09-30. Mathematical proposal and metadata evidence only. This is a separate candidate from the earlier standardized missing-contribution target. No original 117-program definition, weight, assay observability call, biological claim or failed model decision changes.

**Assessment:** a count-weighted ratio cancels an unknown common library scale exactly. It does not recover the original Hotspot score, identify missing expression from the panel, or eliminate gene-specific assay capture. The current data permit a bounded synthetic check, not an independently validated reusable MASLD model. Fitting remains unjustified until a concrete receiving use and compatible target measurement are specified.

## Fixed objects and evidence

The existing [program admission result](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/program_measurement_admission_21995500/summary.json) reports eight complete program targets, 714 admitted observed symbols and the 69/20/20 development-person allocation. Its target roster and one-to-one feature requirements remain unchanged. The parent [measurement report](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/docs/technical/agent_exchange/2026-09-30_codex_program_model_options.md) documents RNA units and external-source limits.

An additional small-file metadata census inspected 7,093 frozen membership rows across 117 programs: **zero negative and zero zero-valued weights** in each of `source_weight`, `canonical_weight_text` and `original_l1_weight`. This is an executed metadata check, not a count/score analysis. The membership SHA256 is `feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b`; admission coverage SHA256 is `e571982d2d522d543492fcc5e3455305cf4ca3a14c4b03e012bf66e59fbf721b`.

## Exact estimand

For one biological person's fixed tissue aggregate i, let `c_ig` be nonnegative, assay-native gene quantities on the admitted complete reference. In GSE268273 these are RSEM expected counts, not absolute transcript molecules. Let `w_pg` be the existing frozen original-L1 weight for program p. Partition its original members into the fixed observed set `O_p` and missing set `U_p`, using only the frozen feature mask and metadata admission.

```
A_ip = sum_(g in O_p) w_pg c_ig
B_ip = sum_(g in U_p) w_pg c_ig
R_ip = B_ip / A_ip                      if A_ip > 0
T_ip = ln(1 + R_ip) = ln((A_ip+B_ip)/A_ip)
```

`A` and `B` have original-weight native count units. `R` and `T` are dimensionless, and `T` is a convenient nonnegative regression target. The logarithm is natural. There is no gene centering, gene logCPM transform, gene SD division, pseudocount or subset weight renormalization in these sums. Both sums use the original weights unchanged. A common multiplication of all weights cancels from R, but the proposed calculation retains the original weights to preserve the count-unit interpretation of A/B.

For any positive person-specific scale d, replacing every `c_ig` by `d*c_ig` multiplies A and B by d and leaves R/T unchanged. In particular,

`[sum_U w*(1e6*c/full_library_total)] / [sum_O w*(1e6*c/full_library_total)] = B/A`.

Thus a complete-library denominator is unnecessary **for this count-weighted ratio**. A deployable panel need not reconstruct the full-library total to provide panel-relative predictor features. The hypothetical corrected count-weighted total is `A*(1+predicted R)` in the observed assay's weighted-count units, only when predicted R refers to the same complete-assay counting rule. Its missing fraction is `R/(1+R)`. Neither quantity restores absolute CPM, cohort-standardized atlas units or lineage activity in mixed tissue.

**Identification limit:** observed counts `[4,2]` can coexist with missing count 10 or 100; with weights `[0.2,0.3,0.5]`, both give identical panel input but R is respectively 25/7 or 250/7. Cancellation removes nuisance scale. Predicting R still needs independently learned covariance and can fail under biological or assay shift.

## Why this cannot recover a standardized weighted projection

The earlier proposed standardized weighted projection, which is distinct from the original Hotspot module score, uses gene transformations that do not factor as a common positive scale:

`sum w*(log2(1+1e6*c/library_total)-training_gene_mean)/training_gene_SD`.

The nonlinear log1p and gene-specific centers/scales prevent denominator cancellation from a ratio of its missing and observed sums. Even before standardization, fixing observed/missing counts at 2/8 while changing complete-library total from 100 to 1000 changes the ratio of their logCPM values, although their raw-count ratio remains 4. After centering, either sum can be negative or cross zero. A ratio then lacks a stable missing-fraction or correction interpretation.

This proposal therefore changes the estimand openly. It cannot relax the previous target contract, reinterpret its negative result, claim restoration of 117 original atlas scores or make unmeasured biological programs observed. A ratio of standardized scores is rejected. If a future frozen membership ever contains signed weights, this specification stops rather than silently taking absolute values or splitting signs. The current metadata are positive.

**Original Hotspot score distinction, verified source review:** `Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py:124–129` specifies raw counts, the depth-adjusted negative-binomial (`danb`) model and the per-cell `n_counts` depth key. Line 277 calls `hs.calculate_module_scores()`. The exported gene weights are separately calculated as mean within-module local-correlation Z at lines 280–288; these are the weights retained by this ratio, not verified module-PC loadings. The installed `hotspotsc` version 1.1.3 implements module scoring by depth/model-centered gene residuals, neighbor smoothing (lambda 0.9), and a newly fitted, sign-oriented PC1 (`hotspot/modules.py:13–43`; called by `hotspot/hotspot.py:534–576` in `/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/lib/python3.12/site-packages`). [Unverified] The historical execution's exact installed package bytes/version were not established by this review. The original script call and separate exported-weight calculation are verified; the logCPM formula above must not be attributed to original Hotspot scoring. This clarification changes no ratio algebra or frozen membership.

**Current scope restriction:** `Analysis/MASLD_Model_Benchmark/OVERALL_PLAN.md:159–160,175–176` permits the 117 programs only as read-only post hoc projections, not training targets or model-selection endpoints. The hypothetical fitting paragraphs below are mathematical proposals, not authorized execution under that contract. No ratio fit or program-input accessibility fit is prepared or approved by this report. Assay-support metadata and separately labeled observed contributions remain descriptive; using the 117 programs to fit or select a learned model needs an explicit scope amendment before execution.

## Zero cases, offsets and eight-target roster

`B=0,A>0` gives R=T=0 and is valid. `A=0` makes R undefined, including `A=B=0`; abstain, do not create a finite pseudo-observation. Nonfinite inputs, negative counts, ambiguous features and incomplete targets are also refusals. Very small positive A can yield large/noisy ratios; scale invariance does not provide a measurement precision guarantee. No lower-count cutoff or practical precision rule is chosen from held outcomes.

The eight complete-target roster remains explicit:

| Frozen program UID | Original members | Observed members | Ratio state from metadata |
|---|---:|---:|---|
| hotspot_cholangiocytes_98aa639a71794d01 |28|3|Structurally possible; requires positive A |
| hotspot_hepatocytes_22bce5d40b6799df |24|3|Structurally possible; requires positive A |
| hotspot_hepatocytes_8a9b435cec18211d |31|5|Structurally possible; requires positive A |
| hotspot_macrophages_31188513a8def476 |31|2|Structurally possible; requires positive A |
| hotspot_macrophages_a4bdecbfa23b855a |20|1|Structurally possible; requires positive A |
| hotspot_macrophages_df2c8f49c9e7ab89 |21|3|Structurally possible; requires positive A |
| hotspot_tcells_05f8b3c103f336ef |20|0|Structural abstention, no program-specific observed denominator |
| hotspot_tcells_3004bf769feb366a |27|8|Structurally possible; requires positive A |

There are seven structurally possible ratios within the unchanged eight-target pilot. None is yet proven positive or variable in real data. The zero-observed T-cell target cannot be rescued by borrowing another program's denominator, using all-panel mass, or assigning missing fraction 1: those are different targets and the fraction would be a trivial constant. Do not enlarge the target roster using relaxed mappings or filtered-source absence as biological zero.

Adding an arbitrary offset a in native count units, `(B+a)/(A+a)` or `ln(B+a)-ln(A+a)`, destroys common-scale invariance. For synthetic A=1.4,B=5,a=1, the offset ratio is 2.5; after multiplying counts by 100 with the same a it is 501/141. `ln(1+B/A)` avoids that offset problem and accepts B=0, but still requires A>0. The 1 is dimensionless, not an added count. Target- or depth-dependent smoothing is not selected here.

## One concrete conditional prediction hypothesis

**Hypothesis:** panel-relative RNA composition predicts T for the seven structurally possible targets with better held-person error than a fit-set mean or scalar observed-program fraction, and the benefit remains in a new source with a declared compatible measurement target.

If subsequently authorized, use only observed 714-gene input counts: `x_ig=log2(1+1e6*c_ig/sum_(h in admitted panel)c_ih)`, with centers/scales learned on fit persons only. Zero panel total abstains. These inputs are invariant to common depth; full-library totals, complete-RNA QC variables, stage, held missing genes and external-cohort centering cannot enter prediction. Preserve the existing 69 fit/20 calibration/20 internal-test allocation and technical-run grouping. This is development, because GSE268273 is already project-exposed.

A sufficient bounded candidate is multivariate ridge on T, with each predicted T truncated at 0 before `R_hat=expm1(T_hat)`. This does not change biological weights. A single shared penalty could be chosen within the existing five donor folds over the already proposed Gram-scaled grid `0,0.01,0.1,1,10,100,infinity`, with larger-penalty tie breaking. Refit all feature/target standardization within each inner fold. No new model architecture is proposed, and no fit or tuning was run.

Necessary comparators are the fit mean of T; a scalar model using the within-program observed contribution divided by panel total; and identical-feature ridge imputing missing counts or relative abundances before constructing R, with its assay units stated. Common-penalty linear imputation followed by a weighted linear sum is not a new independent algorithm; ratio/nonlinear transforms change that equivalence and must be checked explicitly. An observed-only total assumes R=0 and is a useful practical baseline. Full-total agreement alone remains insufficient because shared observed genes supply accuracy mechanically.

The missing-component endpoint is donor-averaged T squared error across the fixed seven possible targets, optionally standardized by fit-only T variance. Report every program, structural and person-level abstention, error conditional on successful prediction, and deployment coverage; do not hide zero-denominator failures by reporting only successful donors. Training-constant T is untestable. An improvement in T is an improvement in a relative count measurement, not new biological regulation or disease prediction. Application tolerances and an external endpoint remain unspecified, so this is a mathematical candidate rather than a ready claim of practical utility.

Biological people, not cells/sections/regions/runs, are inferential units. If uncertainty is later estimated, a paired donor bootstrap must keep each person's full target vector together. The existing 2000-draw seed 20260930 is a reproducible candidate convention, not a coverage guarantee under assay shift. A single aggregate contrast would be primary; any seven program-specific p-values would be one BH family. No inferential test, CI, conformal interval or nominal p-value was calculated here.

## Panel capture and truth-assay mismatch remain

With gene-specific assay capture `c_ig^(a) = d_i^(a) * b_g^(a) * molecule_abundance_ig`, common d cancels but b does not. The ratio becomes `sum_U w*b_g*m / sum_O w*b_g*m`. A panel's observed composition likewise depends on its gene-specific b. A uniform calibration constant is insufficient unless those efficiencies are proportional across genes, which is unverified. Probe effects, isoform/length dependence, segmentation and differing recovered cell mixtures remain.

RSEM expected-read allocations and complete-scRNA UMI counts represent different native measurements. Their dimensionless ratios can differ because relative gene capture/length factors differ. Calling both dimensionless does not make them the same biological estimand. The known seven-person CosMx/complete-scRNA source would require a separately justified assay adapter or a development target native to the same complete-RNA reference; a bulk computational mask alone cannot supply either. This ratio does remove the need to compare unknown total depths, but it does not remove these gene-specific or tissue sampling differences.

## Independent validation and actual scientific use

The existing source metadata resolve seven people with targeted and complete RNA in liver cancer T/B regions, not matched tissue pieces or cells. Same-person regions stay together; there are at most seven biological replicates, not 15 sections or 34 RNA samples. Adjacent non-tumor paired targeted data remain one person. The exact Panel980p probe/control axis is unresolved; the current Govaere 968-name/fixture 714-name mask cannot be reused. Complete-scRNA-assisted source cell labels cannot establish panel-only lineage recovery. Frozen discovery ancestry through `Liver_Atlas` and person disjointness remain unresolved. No source expression object was opened for this proposal.

If those barriers are resolved and a model/adapter is frozen using other people, the source could test **qualified cancer-region technical transport of a separately declared RNA ratio**. It cannot currently validate independent MASLD deployment or a bulk-RSEM biological ratio. Fitting an adapter on the same seven receiving people would consume the proposed confirmation source; retrospective region splits would not restore independence. GSE276114 may instead permit qualified bulk computational masking transport after its count semantics/person joins are established, which would still not validate actual targeted RNA.

There is a plausible narrow use: estimate how much count-weighted program signal an observed subset misses relative to its own observed contribution, report personalized panel completeness and, under a same-assay target, a total/observed correction multiplier. This could identify panels that fail when unmeasured members vary independently. It is not the original atlas score and does not itself resolve metabolic, inflammatory or fibrotic mechanisms. Ratio formation is elementary; ridge imputation and uncertainty alone do not establish methodological novelty. Current evidence supplies no practical error tolerance, independent performance, biological advantage over existing scoring or reason researchers need this output instead of a directly trained assay-native score.

**Disposition:** preserve this as an explicit separate candidate and execute only the synthetic invariants next. Do not fit merely because exposure cancellation is convenient. Advance only if a stated researcher measurement use and independent compatible target/adapter make its error interpretable. It presently does not meet the user's goal of a new independently validated model.

## Checkable artifact and checks

[Synthetic helper](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/program_ratio_invariants.py) verifies frozen-weight metadata, unchanged eight-target roster, seven structurally possible ratios, hand-derived 25/7, common count/weight scaling, `log1p` invariance, count-unit correction, zero-missing acceptance, invalid-case refusal, and counterexamples for offsets, logCPM background, capture bias and unidentifiable missing expression. It reads only membership/coverage metadata and synthetic numbers. Runtime is restricted to a SLURM allocation and saves its executed source hash/environment inventory with results.

[Launcher 1588](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_program_ratio_invariants.sbatch) requests one nslab CPU job, 1 GB, 10 min, no array. Exact submission, for root after inspection:

```
sbatch --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1588.out --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1588.err /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_program_ratio_invariants.sbatch
```

Executed before submission: metadata sign census, Python compilation, launcher `bash -n` and diff whitespace checks. Results: all passed. Root subsequently submitted job 21998523 (1588); scheduler accounting reports COMPLETED, exit 0, elapsed 1 second, and actual stderr [1588_21998523.err](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1588_21998523.err) is empty. The command above is the prepared example; the actual root submission used job-specific log names [1588_21998523.out](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1588_21998523.out) and `1588_21998523.err`.

The [executed synthetic summary](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/program_ratio_invariants_21998523/summary.json) independently reports all checks passing, source SHA256 `6790fa3a76ae663bac9f06abc776168be6469a2e9236e5270b9d8f0657ba45da`, the exact frozen metadata hashes, Python/environment inventory, unchanged eight-target roster, seven structurally possible ratios and one structural abstention. Maximum common-scale ratio discrepancy is `4.440892098500626e-16`. The identical observed-input counterexample yields `25/7` versus `250/7`; offset, logCPM-background and capture-bias counterexamples also passed. Invalid zero-denominator, signed-weight, negative/nonfinite-count and shape cases were refused. This is computational confirmation of the stated algebra and failure cases, **not predictive fitting, observed recovery performance or independent validation**. Memory accounting is reported separately by root; a zero RSS record would not establish actual peak memory.

No expression/count read, fit, protected-data access, dependency installation or authority-file change was performed. Initial worktree had extensive unrelated modifications; they were preserved. This agent submitted no job and updated only this report after inspecting the completed result.
