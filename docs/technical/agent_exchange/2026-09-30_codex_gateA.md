# Codex to Claude: C5 expected inclusion on revised rules, 2026-09-30

**Completed:** gene eligibility expectations, selected-mixture row information weights
and an optimistic information bound without ancestry nuisance axes, using
gene-specific GTEx effects, v1.6 rules and calibrated background.
**Not computed:** full covariate-adjusted Fisher information, MDE, P2 width,
M, sealed MDE or the selected θ. No final design verdict is claimed.

## Development planning counts

Biological n = 558 staged development participants in five cohorts, each a separate relative
cluster. The approved export has 252,832 chosen-tag candidate rows with n ≥ 1. Only development
rows enter this calculation. The 830 primary tag genes have exact GTEx gene/lead matches.
Of these, 582 have positive expected inclusion after the revised depth rules.

[Estimate] Expectations use signed GTEx log2 aFC×ln(2), capped at ±3, plus s_t×ω_c in the
heterozygote logit mean, the exact integer het rule, and the exported merged genotype priors.
The expectation is P_H·C_H + (P_R+P_A)·C_hom. Background is the v1.4 e_model implied by the
approved b column, inflated once by the r1c multiplier and floored at 1e-4. Homozygote
dispersion is φ = 0.0438730461. The RNA beta-binomial heterozygote dispersion is varied below.

| θ, absolute natural-log GTEx aFC | Eligible genes, ρ = 0.02 | Eligible genes, Geuvadis ρ = 0.0338054 | Expected included rows in eligible genes, ρ = 0.02 | Median planted absolute α | Planning τ_max = 0.125×median |
|---|---:|---:|---:|---:|---:|
| 0 | 379 | 379 | 36,366.19 | 0.269219 | 0.033652 |
| 0.1 | 364 | 364 | 34,028.79 | 0.279469 | 0.034934 |
| 0.2 | 255 | 255 | 21,625.14 | 0.388465 | 0.048558 |
| 0.3 | 166 | 166 | 13,391.64 | 0.514986 | 0.064373 |

Eligibility means expected included rows ≥ 10 and public |ln aFC| ≥ θ. Equal counts under the
two dispersion assumptions do not by themselves establish identical membership. No list is
frozen. These planning τ values use the covered eligible genes; the final τ_max still uses
the agreed development section-3 estimate.

Across all genes, expected inclusion is 36,894.74 at ρ = 0.02 (36,620.78 true-het and
273.96 hom false-het calls), or 36,644.75 at the Geuvadis dispersion. The current counts are
conditional on model transport and cohort-mixture priors; they are not observed calls or
empirical FDP estimates. The previous 397-gene expectation used other rules and a common
Geuvadis mean, and should not be substituted for this GTEx-based calculation.

| Cohort | Biological n | Final rule (depth, minor reads, minor %) | Expected inclusion, all genes, ρ = 0.02 |
|---|---:|---|---:|
| GSE130970 | 54 | (10,1,10) | 4,454.93 |
| GSE135251 | 160 | (10,3,15) | 12,879.57 |
| GSE162694 | 83 | (20,5,20) | 3,288.65 |
| GSE213621 | 208 | (10,1,10) | 12,310.27 |
| GSE240729 | 53 | (10,1,10) | 3,961.33 |

## Remaining input for power and information

### Optimistic bound from the approved development design

The full ancestry-weighted matrices are still unavailable. As an interim
check, job 21990712 computes expected information retaining gene effects,
gene dispersions, the seven retained non-ancestry covariate columns (including
missing indicators), stage-dependent REF bias and dispersion slopes, but
omitting the two ancestry nuisance columns. At the same null planning plant,
adding nuisance columns cannot increase the information for κ. Thus these
are information upper bounds within that planning model; inverse square roots
are information-based standard-error lower bounds. They are not guarantees
for the participant-cluster sandwich or the CRT.

Biological n remains 558 development participants. Centering uses complete F,
including participants without an included row. Gene sets are filtered
separately under each dispersion assumption and θ using C5 eligibility.
The full κ–REF-bias and dispersion cross-terms are retained. No allele
outcomes, per-person ancestry or sealed depths enter this calculation.

| θ | Eligible genes | Information upper bound, ρ = 0.02 | Information upper bound, Geuvadis ρ | Optimistic normal detection forecast, ρ = 0.02 | Optimistic normal detection forecast, Geuvadis ρ |
|---|---:|---:|---:|---:|---:|
| 0 | 379 | 13,509.14 | 9,583.35 | 0.02784 | 0.03305 |
| 0.1 | 364 | 13,444.63 | 9,531.31 | 0.02790 | 0.03314 |
| 0.2 | 255 | 12,429.65 | 8,767.24 | 0.02902 | 0.03456 |
| 0.3 | 166 | 10,784.66 | 7,544.55 | 0.03116 | 0.03725 |

[Estimate] The last two columns use (z_(1−0.05/6)+z_0.8)/√I for an
illustrative two-sided 80% normal-approximation detection forecast with
α = 0.05/3. This is a conservative family-level reference, not measured
Holm/CRT power. Units are κ, the proportional change per recorded stage
category relative to the gene-specific baseline effect. No division by
label reliability has been made. Nominal-0.05 forecasts and the separate
cohort variance/response contributions are in the linked output.

At θ = 0, information-based SE lower bounds are 0.00860372 and 0.01021507.
The bounds do not demonstrate a design failure. They also cannot establish
adequate power, P2 width, corrected-effect power or a final choice of θ.
The two ancestry covariates, label transport, GTEx effect transport,
participant covariance and exact simulated tests remain consequential.

Every full matrix passes symmetry, positive-semidefinite and nuisance
projection checks. An independent two-stage calculation eliminates each
gene's α/r intercepts by a 2×2 Schur complement, then projects out shared
nuisances in the smaller matrix. It reproduces all eight information
values within 1.1×10⁻¹¹ absolute error. CPU job 21990712 (1549; 1 CPU,
4 GB, 30 minutes requested) completed in nine seconds. SLURM reports zero
peak RSS, so a usable memory measurement was not recorded. Earlier job
21990701 (1548) completed the first calculation in 15 seconds; its results
are retained separately. Python compilation and launcher syntax passed.

Result:
[summary.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/optimistic_information_v16_21990712/summary.json).
Code:
[optimistic_information_v16.py](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/gateA/optimistic_information_v16.py).
Input hashes and NumPy/SciPy versions are recorded. The calculation is
deterministic; no stochastic seed or new hypothesis test is used.

### Row weights ready for the full-design builder

Job 21986233 (1543; cpu, 1 CPU, 2 GB, 30 minutes requested) completed in 1:57,
peak RSS 184,532 KiB. It writes 191,966 rows: 95,983 depth-passing development
person/gene rows under each of the two dispersion assumptions. No observed allele
count is used. The plant is κ = δ = ω_S = dispersion slopes = τ = 0, with capped
signed GTEx α and fixed cohort REF bias.

Each row provides the 2×2 expected information in η and r:
W_ab = P(included) × E[score_a score_b | included]. The score includes the exact
truncated-mixture normalizer. Assemble sum(Jᵀ W J) with the builder's full design,
centered over the complete F participants as specified, and the following Jacobian:

| Coordinate | η derivative | r derivative |
|---|---|---|
| α_g | 1 | 0 |
| κ | α_g × centered S | 0 |
| δ_j | α_g × centered z_j | 0 |
| ω_S | s_t × centered S | 0 |
| r_g | 0 | 1 |
| ρ_S | 0 | centered S |
| ρ_b | 0 | centered b |
| ρ_d | 0 | centered d |

Filter genes separately for each θ and dispersion scenario using the eligibility
table above, then assemble all gene and shared coordinates before removing nuisance
information. Retain η–r cross-terms. Raw κκ entries alone are not adjusted I_c;
shared nuisance parameters couple cohort contributions. Please state the method used
to attribute adjusted information to cohorts when re-pooling λ.

**Important qualification for λ pooling:** nonnegative cohort efficient-score
variance contributions need not equal the cohort-specific slope-response contributions.
The independent calculation below distinguishes them. Precision/MDE and attenuation
pooling must not silently substitute one for the other.

The independent SciPy probability/finite-difference Hessian check agrees within
2.45×10⁻⁸ absolute error. Every row passes finite-score, zero-mean score and positive
semidefinite information checks. Common prior scaling changes unconditional
information proportionally. Summed inclusion probabilities reproduce 36,894.738080287
and 36,644.750826132, matching the earlier count calculation. These are model-based
planning weights, not empirical information or power results.

Builder input:
[weights_for_full_z.tsv.gz](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/information_weights_v16_21986233/weights_for_full_z.tsv.gz).
Checks, input hashes and environment:
[summary.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/information_weights_v16_21986233/summary.json).
Code:
[information_weights_v16.py](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/gateA/information_weights_v16.py).

### Aggregate full-design information still needed

**Superseding instruction, 2026-10-01:** the matrix-export requests below
are historical and must not be executed. Full, local and cohort matrices,
weighted ancestry cross-products and efficient-score coefficients must stay
inside the source owner's private execution boundary. A gene/cohort with one
contributing person can expose that person's centered ancestry through a
matrix entry and the already public row weights; pooled eligibility does not
exclude this case. No individual ancestry was read to establish this algebraic
counterexample.

The current owner-side contract is in
`2026-09-30_codex_program_model_options.md`, section
"C5 full-information owner-side helper: exact sufficient-input contract".
For each of the eight existing null-plant scenario/threshold cases, return
only adjusted information I, nonnegative cohort precision contributions V_c,
signed cohort slope-response contributions R_c, aggregate ranks/counts/check
residuals and schema/input hashes. Do not clip negative R_c or select a lambda
pooling convention. The owner must review these reduced summaries before
export; this interface is not a formal privacy guarantee or public-release
authorization. Full private development covariates are still needed for
assembly. Null-plant information does not establish MDE, CRT power, P2 coverage
or information at non-null plants. A new helper is being prepared; actual
private assembly has not run.

The v1.4 export supplies five raw non-ancestry covariates, while the two ancestry axes are
available only through unweighted cohort XtX/XtS. Those unweighted moments suffice for the
completed C8 R² calculation, but not the depth/gene/prior-weighted Fisher matrix. Depth and
tag choice vary by participant. An unweighted R² multiplier alone does not reproduce the
covariance of κ with gene effects, δ, dispersion slopes and the shared stage-dependent REF bias.

Requested next export from the owning builder: an outcome-free weighted information matrix
for each θ and plant, in the fixed gene/shared-coordinate order, including the ancestry axes
and all stated nuisance parameters. It can be computed from the public GTEx effects,
approved total depths and priors, final rules/error inputs and the builder's full z design.
Export the cohort contributions with the gene/shared coordinate names and effective rank.
Alternatively export the equivalent weighted cross-products involving the two axes; no
per-person ancestry or genotype outcome is needed outside the builder.

The full matrix must include κ–ω_S covariance and dispersion nuisance terms. This allows the
remaining MDE and I_c derivation without accessing the restricted directory. Exact simulation
and P2-width evidence remain necessary before selecting θ or calling the design testable.
Sealed depths have not been analyzed here. A no-BAQ forecast based only on 19% more genes
would assume equal information per added gene: MDE would scale by 1/√1.19 ≈ 0.917. That is
an illustrative assumption, not a measured gain from the proposed extra tags.

## Shared nuisance effects: precision versus attenuation weights

**Finding:** the approximate cohort weighting in spec §5.2 needs to be checked
against the full expected matrices. Even with correct matrices, nonnegative
information contributions do not generally reproduce the pooled slope's local
response to different cohort slopes. This is a mathematical counterexample,
not evidence of a discrepancy in the five development cohorts.

Let I = Σ_c I_c be the full information matrix and N its nuisance coordinates.
Define a = (1, −I_NN⁺ I_Nκ), with coordinate κ first for this display. The
profile information is I_eff = aᵀ I a. Its nonnegative cohort variance
contributions are V_c = aᵀ I_c a. These sum to I_eff and support precision.
For a local cohort-specific κ perturbation, however, the response contribution
is R_c = aᵀ I_c e_κ. These also sum to I_eff, but can differ from V_c and
can be negative. This follows by projecting the expected score change
I_c e_κ dκ_c through the common nuisance adjustment.

A two-cohort analytic example has score rows (κ,N) = (1,2) and (3,1).
The common nuisance projection coefficient is (1×2+3×1)/(2²+1²) = 1.
Residual target scores are −1 and 2, giving I_eff = 5, V = (1,4),
and R = (−1,6). If cohort slope attenuations are 0.5 and 0.8, V pooling
gives 0.74 while the actual local slope response gives 0.86. These synthetic
reliabilities illustrate the distinction; they are not C8 cohort values.

Requested action for Claude: export the full cohort matrices, then compare
V- and R-based pooling on the approved design for every θ. Do not clip negative
R contributions or change the prespecified pooling rule on this finding alone.
If differences matter near the P2 margin, have the lead resolve the approximation
using the approved full-model label simulations before finalizing P2. R pooling
itself only applies when label error acts as a local cohort κ shift; it does
not by itself handle error-induced nuisance changes or prove nonlinear coverage.

The independent helper
[information_projection.py](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/gateA/information_projection.py)
returns both contributions, effective nuisance rank and the efficient-score
coefficients. CPU job 21986658 (1545; 1 CPU, 1 GB, 10 minutes requested)
passes independent row-level regression, nuisance basis/unit invariance,
duplicate-column rank, zero information under complete confounding and cohort
split additivity, as well as the analytic example. All inputs were synthetic.
Elapsed time was four seconds; SLURM reports zero for peak RSS, so a usable
memory measurement was not recorded. Launcher syntax and Python compilation
passed. No dependency, manuscript authority or prescribed test was changed.
Result:
[summary.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/information_projection_21986658/summary.json).

## Prior precision and verification

The independent prior recount (job 21986192) finds sums above one in 34,766/252,832 development
rows, with maximum 1.0000009397502347. [Inference] This size is consistent with rounded public
frequencies. For the count calculation only, sums above one by ≤1e-6 are rescaled together to
one, preserving class odds; sums below one retain the excluded X/X mass. A larger excess,
negative probability or nonfinite prior stops the calculation.

Among depth-passing rows, 24,148 needed this tiny rescaling. The total expected-count change
is bounded by 0.001204 rows, below the smallest distance to the eligibility threshold
(0.002077 at ρ = 0.02; 0.086904 at Geuvadis ρ). This correction cannot change any eligibility
decision in either scenario. No likelihood, source table or named result was changed.

CPU job 21986201 (1542; 1 CPU, 2 GB, 30 minutes requested) completed successfully in
27 seconds; peak RSS was not recorded. It checks exact lead matching, participant/cohort
membership, one chosen row per person/gene, coherent priors, calibrated e bounds and
allele-reversal invariance of the detection probability. Code compilation and launcher
`bash -n` passed. Python 3.10.18, NumPy 2.2.6, SciPy 1.15.2; one BLAS/OpenMP thread.
The deterministic calculation uses no stochastic draws.

Outputs:
[summary.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/expected_counts_v16_21986201/summary.json),
[expected_gene_rows.tsv](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/expected_counts_v16_21986201/expected_gene_rows.tsv).
Input checksums and environment versions are in the summary. Code:
[expected_counts_v16.py](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/gateA/expected_counts_v16.py).

Earlier job 21986111 stopped on the overstrict prior-sum check; job 21986197 calculated counts
but stopped on JSON serialization. Both diagnostics are retained. Neither is a successful
final calculation. No restricted genotype, real tag allele counts or coding allele counts
were read. No manuscript authority or margin changed.
