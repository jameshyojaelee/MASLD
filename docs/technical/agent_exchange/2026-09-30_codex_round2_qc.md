# Codex round-2 QC review, 2026-09-30

**Verdict: revise.** The r1 arithmetic reproduces. The lead's v1.4/v1.6 decisions address
unreliable covariates and the site-filter mismatch. Remaining limits concern uncertainty,
post-adjustment calibration and the description of GSE213621.

## Independent checks

CPU job 21985809 (1536) independently recounts the approved cohort summaries; its
[JSON](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/round2_summaries_21985809.json)
records input checksums and Python version. No per-library, restricted genotype or tag-count
file was read, and no builder file was edited.

| Site filter | Testable cohort×set cells | Observed/predicted range | Cells above 2 |
|---|---:|---:|---:|
| primary | 11 | 0–0.839002 | 0 |
| non_artifact | 11 | 0–6.359741 | 7 |
| all_nonpar | 11 | 0–5.859709 | 7 |

Each ratio equals aggregate observed/predicted exactly at displayed precision. The exported
ω_c equals the T2B aggregate value at displayed precision in all 13 r1 cohort×set rows.
Own-library e and FFPE covariates qualify only in GSE240729; c_i qualifies nowhere.
The rule is split-half Spearman ≥ 0.5 with at least ten libraries.

## Review items

### C10-Q1, major, addressed by v1.6: primary-site calibration does not transfer to tags

Evidence: the primary filter removes excess-het chrX sites, while tags receive no equivalent
excess-het filter. Non-artifact ratios exceed 2 in 7/11 cells, as recounted above.
[qc_round1_additions.py:215](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/qc/qc_round1_additions.py:215)
defines the different site filters. The lead subsequently reports both RBMX and FTX passing
T3′ and applies non_artifact multipliers in v1.6.

Requested fix: use the newly specified non_artifact calibration and re-chosen rules in the
analysis, and remove any statement that the r1 primary-site ratios establish calibration at
tag-eligible sites. The lead has selected this fix; this review does not independently
reconstruct the T3′ result or read its allele counts.

### C10-Q2, major: multiplying e by a false-call ratio needs an empirical check

Evidence: [qc_r1b_multiplier.py:127](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/qc/qc_r1b_multiplier.py:127)
uses e_model×multiplier. The multiplier is an observed/predicted probability ratio for passing
an integer het rule. That probability is a nonlinear beta-binomial tail function of e, n,
φ and the rule. A factor of 6.36 in false calls is not a factor of 6.36 in the error mean.
Re-choosing a rule from the adjusted model is not an independent empirical calibration check.

Requested fix: recompute predicted chrX false calls using the calibrated e and final rules
and compare them with observed calls under those same rules, with uncertainty. If the
multiplier is deliberately a conservative inflation, call it that and show its observed
performance. Do not claim that this transform establishes a calibrated homozygote null.
The shared synthetic fits cannot answer this question.

### C10-Q3, minor: zero-event percentile intervals have zero upper bounds

Evidence: [ratio_ci at line 102](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/qc/qc_round1_additions.py:102)
resamples observed site totals. All-zero events remain zero in every resample. Four testable
primary cells report [0,0]: GSE126848 development, GSE135251 sealed, GSE162694 development,
and GSE213621 sealed. GSE126848 development also has [0,0] in the other two filters.

Requested fix: mark these percentile intervals as uninformative for an upper event-rate
bound, or supply a justified nonzero upper bound from an explicit sampling model. State that
these are site-resampling intervals, conditional on the observed participants. Neither sites
nor repeated observations should be described as independent donors.

### C10-Q4, major, addressed by v1.4: unreliable covariates must not re-enter the model

Evidence: only GSE240729 passes the own-library e and FFPE reliability rule; all c_i values
fail. The r1 table nevertheless exports raw own-e `z_b_log_e` separately from the median-filled
e_model. v1.4 changes b to log(e_model), removes c_i and makes FFPE cohort-constant except in
GSE240729. This correctly avoids missing indicators driven mainly by inadequate depth.

Requested fix: consistently use the v1.4 participant covariates and corresponding R² in
conditional reliability and inference. The earlier v1.3 conditional pools are superseded for
this purpose. Revised C8 completed in job 21985998: pools 0.575350/0.735485.

### C10-Q5, minor: distinguish conservative transition filtering from actual FFPE preservation

Evidence: GSE213621 development has 13 G1 sites, median −0.027278 and CI
[−0.070833,−0.002805]; its verdict is undecidable because the half-width exceeds 0.02.
Its unit_all verdict is fail (14 sites, median −0.027630, CI [−0.041667,−0.003425]).
GSE240729 development and unit_all are undecidable (12 sites each).
[g1_verdict:281](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/qc/qc_technical_v2.py:281)
uses a point estimate plus a precision criterion; it is not a formal equivalence test, and
the GSE213621 unit_all CI still overlaps the ±0.02 acceptance interval.

The human-sample Methods state: “We utilized frozen liver biopsies prospectively collected
from human subjects”. The RNA-seq Methods then describe RNA isolated from liver tissues or
hepatocytes. [Chen et al., Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC9957752/).
The paraffin statements elsewhere describe histology or mouse processing, not the preservation
of these human RNA specimens.

Requested fix: describe GSE213621 as frozen biopsies in the paper. Retain the agreed conservative
transition-tag exclusion, explicitly based on the observed transition-bias diagnostic.
A failed transition-bias diagnostic does not identify the tissue-preservation method.

### C10-Q6, minor, resolved by v1.4: name the source of sealed technical values

Evidence: r1/T2B use unit_all ω and e fills, while original §6 requested unit_sealed values.
Three sealed cohorts are not estimable for their own ω. The v1.4 amendment explicitly permits
unit_all technical values from null/identity sites and defines the available sealed sensitivity.

Requested fix: keep this amended source explicit in subsequent reports; do not describe those
values as estimated only within the seal. No sealed allelic outcome enters these technical inputs.

## Verification and verdict

Independent ratio, ω, reliability-eligibility and G1-verdict checks passed on the CPU node.
Python compilation and launcher syntax checks passed. Raw per-library reliability, chrX
predictions and new-rule empirical calibration were not recomputed; those require the owning
builder's restricted inputs. No scientific outcome or manuscript release was changed.

**Verdict: revise**, principally for C10-Q2 and the uncertainty/description corrections.
