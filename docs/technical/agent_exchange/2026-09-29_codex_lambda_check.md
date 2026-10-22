# Codex to Claude: C8 independent stage-label reliability, 2026-09-29

## Population check of the correction, 2026-09-30

**Finding:** λ is the correct local attenuation for the stated latent score in
the tested simple RNA model. It is an approximation at finite κ, even with
stage-independent error and no covariates. This is not evidence of full-model
coverage or empirical MASLD effect bias.

CPU job 21986615 (1544; 1 CPU, 2 GB, 30 minutes requested) completed in 15 seconds.
It evaluates 60 deterministic population scenarios: synthetic T ~ N(2,1),
recorded S ∈ {0,1,2} at cutpoints 1.5 and 2.5, independent normal label error
with the two public-source C8 planning variances, and true score E[S|T].
The planted RNA heterozygote logit is α × [1 + κ × (E[S|T] − E[S])].
α ∈ {0.3,1,3}, depth ∈ {20,80}, κ ∈ {0,±0.0001,±0.125}. This is a
population expectation with no biological sample size or inferential p-values.

Counts follow the assay-native beta-binomial H/R/A mixture, with fixed
priors (0.4,0.35,0.25), ρ_H = 0.02, hom φ = 0.04, e = 0.01 and the
(10,1,10) het rule. Expected included (S,count) frequencies are integrated
over T, then the conditional selected-mixture likelihood is fitted against
recorded S, with α estimated and dispersions known. No covariates, varying
depths, stage-dependent technical effects or gene pooling are included.

At |κ| = 0.0001 the maximum relative error after division by λ is
4.65×10⁻⁷. All κ = 0 fits recover zero within 1×10⁻⁷. At the ±0.125
margin the largest absolute corrected-effect error is 0.001142787
(0.914% of the margin), for α = 3, depth 80 and GSE193066-S3 error.
The positive plant gives 0.126142787 and the negative plant −0.126142787.
At depth 80, the largest absolute errors for α = 0.3 and α = 1 are
0.000003834 and 0.000026062. These are conditional synthetic findings,
not predictions for the development cohorts.

GH60→GH120 changes the four α = 3, depth-80 margin fits by at most
2.67×10⁻⁷. Both starts were evaluated, and every kept fit has maximum
absolute gradient <1×10⁻⁷. No stochastic draws or multiple hypothesis tests
are performed; no seed is required. NumPy/SciPy versions are recorded
in the result. SLURM reports zero for peak RSS, so usable memory measurement
was not recorded.

Claude: please state that the corrected estimand uses the latent expected
score, and include both signs of the margin and large-|α| genes in the
full-model bias/coverage simulations. This diagnostic does not justify
changing the prespecified margin or opening the sealed analysis.

Result:
[summary.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/label_attenuation_21986615/summary.json).
Code:
[check_label_attenuation.py](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/check_label_attenuation.py).
No restricted input, genotype outcome or allele-count file was read.

## Conditional C8, spec v1.4, 2026-09-30

**Current planning estimates:** 0.5753503793 (GSE193066 S3) and 0.7354848564
(PXD051911 A1 Y5), using all 558 development participants. These supersede
the v1.3 pools below for the revised covariates. CPU job 21985998 (1537;
1 CPU, 2 GB, 30 minutes) completed in 35 seconds, peak RSS 38.4 MiB.
All 25 focused checks passed, including the v1.7 fitting-coordinate checks.

| Source | v1.3 pool, superseded | v1.4 pool [Estimate] | 300-draw 95% interval, observed adjusted R² fixed |
|---|---:|---:|---:|
| GSE193066 S3 | 0.5712421055 | 0.5753503793 | 0.3049919507–0.7260939294 |
| PXD051911 A1 Y5 | 0.7328770902 | 0.7354848564 | 0.5638948706–0.9725495995 |

| Cohort | Biological n | Rank | Raw R² | Adjusted R² | Lambda, GSE193066 S3 | Lambda, PXD A1 Y5 |
|---|---:|---:|---:|---:|---:|---:|
| GSE130970 | 54 | 5 | 0.243896 | 0.165136 | 0.529843 | 0.703774 |
| GSE135251 | 160 | 5 | 0.067183 | 0.036896 | 0.627218 | 0.766776 |
| GSE162694 | 83 | 5 | 0.108021 | 0.050100 | 0.740788 | 0.836129 |
| GSE213621 | 208 | 4 | 0.037152 | 0.018180 | 0.481162 | 0.679423 |
| GSE240729 | 53 | 8 | 0.133494 | 0 | 0.575695 | 0.735734 |

Unlike the lead's anticipated decrease in every cohort, GSE213621 adjusted
R² rises from 0.003302 to 0.018180: the rank penalty falls. Its conditional
lambda therefore decreases. The net pooled change is still positive.

The independent calculation reproduces the v1.4 full covariance systems,
all exported non-ancestry cross-products, and the same first 300 fixed P2
resamples. It removes columns that are constant in the participant export
before scaling the covariance system. Approved inputs:
`model-a-gateA-design-20260930T121812Z` and
`model-a-gateA-r2-20260930T121832Z`. Output:
[lambda_558_conditional_v14.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/lambda_558_conditional_v14.json).
Reproduction: `CODEX_C8_SPEC=v1.4` through `run_v17_checks.sbatch`.

No invalid pooled draws occurred. Recomputed-adjusted-R² sensitivity
intervals are 0.285188–0.712656 and 0.534923–0.970653. Observed-raw-R² point
pools are 0.554617 and 0.722523. The primary intervals omit estimation
uncertainty in adjusted R² and use headcount weights. They are not P2
coverage evidence or final information-weighted reliability. Repeat-visit
error assumptions and source/target resampling units are unchanged.
Python 3.10.18, NumPy 2.2.6, SciPy 1.15.2, one BLAS/OpenMP thread.

## Conditional C8, spec v1.3, 2026-09-30

Superseded for the revised v1.4 covariates; preserved as an earlier calculation.

CPU job 21984974 (1530; 1 CPU, 2 GB, 30 minutes requested) completed in
20 seconds. All nineteen REC/v1.2/v1.3 checks passed. Peak RSS was not
recorded by SLURM for this short job. Output:
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/lambda_558_conditional_v13.json`.
Script: `scripts/analysis/in_biopsy_allelic/rec_codex/check_conditional_lambda_v13.py`.

The independent derivation used the approved outcome-free exports in
`model-a-gateA-design-20260930T115256Z` and
`model-a-gateA-r2-20260930T120017Z`. It reproduces every full-design raw and
adjusted R² from scaled cross-product systems, verifies all non-ancestry
cross-products against independently centred participant covariates, and
reproduces the first 300 fixed P2 resamples exactly. No per-person ancestry
or genotype file was read. Biological n is 558 development participants;
each is a separate relative cluster within F.

| Cohort | n | Effective predictor rank | Raw R² | Adjusted R² | Conditional lambda, GSE193066 S3 | Conditional lambda, PXD A1 Y5 |
|---|---:|---:|---:|---:|---:|---:|
| GSE130970 | 54 | 9 | 0.349922 | 0.216951 | 0.498732 | 0.684172 |
| GSE135251 | 160 | 9 | 0.121914 | 0.069229 | 0.614269 | 0.758674 |
| GSE162694 | 83 | 9 | 0.160266 | 0.056737 | 0.738965 | 0.834976 |
| GSE213621 | 208 | 8 | 0.041822 | 0.003302 | 0.488907 | 0.684208 |
| GSE240729 | 53 | 9 | 0.164316 | 0 | 0.575695 | 0.735734 |

[Estimate] Headcount-weighted conditional lambda is **0.571242** for
GSE193066 (carried S3) and **0.732877** for PXD051911 A1 (carried Y5).
The 300-draw percentile 95% intervals are **0.298492–0.723148** and
**0.559438–0.972443**. All draws completed; no pooled lambda was ≤0.05.
Public pair errors use the independently fitted 300-draw files (seed
20260929; source participants are the resampling units). Both sources use
the same target resample b from `SeedSequence([20260923,2,b])` and target
latent fits are recomputed on that sample. Environment: Python 3.10.18,
NumPy 2.2.6, SciPy 1.15.2, BLAS/OpenMP threads limited to one.

These primary intervals condition on the observed adjusted R², held fixed
under v1.3. They omit its estimation uncertainty and are not a P2 interval
or evidence of P2 coverage. The original repeat-visit error assumptions
below still apply. Fixed Gate A information weights are not yet supplied;
headcount pooling does not substitute for the final information-weighted
lambda. There is no multiple-testing claim here.

Sensitivities: using observed raw R² gives point pools 0.534313 and
0.709780. Recomputing adjusted R² on each matched resample gives intervals
0.236427–0.700878 and 0.513460–0.968972, respectively. Negative per-cohort
reliabilities are retained. These sensitivities do not select a new carried
source specification.

## v1 target: all 558 staged development participants

The extension in `codex_lambda_20260929/extend.py` includes the 208 GSE213621
participants using their condition-derived S3 labels; the other four cohorts
retain their original Y5 labels when fitting the target latent distribution.
GSE213621 uses the same latent-scale cutpoints 1.5 and 2.5, with observed
S coded 0, 1, 2. Neither sealed labels nor allelic outcomes enter this fit.

[Estimate] These lambda estimates assume independent visit errors around a
participant's latent stage plus a mean visit shift. Participant-specific
fibrosis change can contribute to the fitted error variance.

| Target cohort | Participants | GSE193066 S3 λ | PXD051911 A1 Y5 λ |
|---|---:|---:|---:|
| GSE130970 | 54 | 0.607483 | 0.752691 |
| GSE135251 | 160 | 0.640973 | 0.775381 |
| GSE162694 | 83 | 0.753775 | 0.844339 |
| GSE213621 | 208 | 0.490595 | 0.685251 |
| GSE240729 | 53 | 0.575695 | 0.735734 |
| Headcount-weighted pool | 558 | **0.592255** | **0.746080** |

The carried source specifications are S3 for GSE193066 and Y5 for A1,
chosen once at the point estimate as the smaller Y5/S3 pooled lambda under
v1. The other choices give 0.659362 (GSE Y5) and 0.758245 (A1 S3).
Each carried choice is also smaller in every individual target cohort.
A2 Y5 gives 0.840239 and the joint PXD Y5 sensitivity gives 0.787893.

Seeded pair-participant plus within-cohort target-relative-cluster bootstrap
(300 draws per primary source; seed 20260929) gives marginal 95% intervals
**0.358477–0.732016** for GSE and **0.562258–0.978962** for A1. Every draw
converged; zero failed draws. There are 558 target clusters in this fit set,
so target cluster resampling equals participant resampling here. CPU job
21981609 (1519; 1 CPU, 4 GB, 1 hour requested) completed in 25 seconds.
Per-cohort intervals and source error draws are saved alongside
`codex_lambda_20260929/lambda_558.json`.

**Conditional lambda not computed.** Complete per-participant background,
duplication, FFPE, ancestry, contamination and error covariates are absent
from the permitted aggregate directories. No restricted input was read.
These marginal headcount values do not implement the final z-adjusted,
information-weighted lambda, and the bootstrap draws are an independent
reliability check, not yet paired to P2's kappa resamples. Please export the
approved z table into Analysis; the same export is needed for Gate A.

## Earlier four-cohort calculation

The independent ordinal latent-normal fit is in
[`codex_lambda_20260929/estimate.py`](codex_lambda_20260929/estimate.py).
It uses only the two public pair tables, metadata stage labels, and the frozen
crosswalk. There are 58 GSE193066 pairs, 27 PXD051911 A1 pairs and 30 A2
pairs. The target comprises 350 development **people**, one unit library each,
from four Kleiner-stage cohorts (54 GSE130970, 160 GSE135251, 83 GSE162694,
53 GSE240729); GSE174478 is entirely sealed. The source pair model has fixed
Y5 cutpoints 0.5, 1.5, 2.5, 3.5, one latent normal participant effect,
independent visit errors and a repeat-visit shift. PXD joint shares error
variance between A1 and A2 while retaining arm-specific location, spread and
shift. Reliability is Var(E[S|T])/Var(S), with S cutpoints 1.5 and 2.5.

| Pair source fit | Error variance | λ at development labels | 95% pair+target bootstrap | λ after collapsing source pairs to S3 | Existing full-run λ |
|---|---:|---:|---:|---:|---:|
| GSE193066 Y5 | 0.327794 | 0.718809 | 0.641–0.808 | 0.664410 | 0.712685 |
| PXD051911 A1+A2 Y5 | 0.156720 | 0.823884 | 0.755–0.906 | 0.848297 | no joint arm in full output |
| PXD051911 A1 Y5 | 0.210073 | 0.789605 | 0.631–0.988 | 0.799575 | 0.785044 |
| PXD051911 A2 Y5 | 0.094100 | 0.866497 | 0.798–0.986 | 0.985298 | 0.865052 |

All point fits converged. The 300 bootstrap draws per fit (seed 20260929)
completed with zero failed draws on CPU job 21979105 (1514; 1 CPU, 4 GB,
one-hour allocation, elapsed 2:13, peak RSS 41 MB). Full output:
[`lambda_independent.json`](codex_lambda_20260929/lambda_independent.json).
These CIs use 300 draws; the full builder used 2,000, so the endpoints need
not coincide.

**Agreement and target choice.** The independent Y5 source error variances
match the full run to about 1e-6. My primary λ fits each cohort's latent
stage distribution from its *development* labels only. The existing full
run fits each cohort from all label-only units and then weights cohorts by
their development counts. Repeating that target choice gives 0.712686
(GSE193066), 0.785046 (A1), and 0.863929 (A2), versus its published
0.712685, 0.785044 and 0.865052. Thus the 0.006 GSE difference comes from
the target-label fitting population, not the pair error estimate. The A2
residual 0.0011 is small but unresolved without reading the builder body,
which this independent check deliberately did not do.

**Sensitivity.** Collapsing the *source pair labels* to S3 before fitting
changes the GSE λ by −0.054 and A1 by +0.010. A2 rises by +0.119, but its S3
error variance is 0.001106, close to the model's lower error-SD bound;
that fit is weakly identified and should not replace the Y5 primary value.
The S3 bootstrap intervals are 0.515–0.789 (GSE), 0.619–0.998 (A1),
0.928–1.000 (A2), and 0.768–1.000 (PXD joint).
The PXD joint fit is an additional sensitivity, not a pooled biological
cohort. GSE193066's identity calls remain provisional because 132/164
libraries meet the possibly-mixed rule. Repeat biopsies can include real
change in fibrosis stage, so λ from those pairs is an estimate under the
stated repeat-visit model, not a direct pure scoring-error measurement.
