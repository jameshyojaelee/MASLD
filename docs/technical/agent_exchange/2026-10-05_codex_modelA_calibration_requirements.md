# Model A numerical repair and prospective calibration requirements

2026-10-05. User authorized proceeding with numerical repair and independent
validation-contract work. This records the existing specification's requirements;
it does not replace its thresholds, choose a new model from outcomes, or open
protected biological data. Only synthetic inputs are authorized for the numerical
repair tests. Original and completed parallel P3 results remain unchanged.

## Binding scientific target

The target is a common proportional association of gene-specific cis allelic
effects with cross-sectional fibrosis categories, with gene-specific departures
and a prespecified persistence bound. These are bulk RNA allelic counts, with
participant/relative clusters as biological units. Mixture genotype uncertainty,
allele orientation, read-based inclusion and the inclusion normalizer belong to
the likelihood. Neither genes nor Monte Carlo draws add biological replication.

The binding source is
[the likelihood specification through v1.7](2026-09-29_model_A_likelihood_spec_v1.md),
including its Appendix B user decisions and
[the design/fixture answers](2026-09-30_claude_to_codex_design_fixture_answers.md).
Current fixed-P3 accounting is 25 converged and 25 failed/nonconverged draws;
its p=1 at observed LR=0 does not establish calibration.

## Numerical checks before another inferential experiment

Repair candidates must preserve the mathematical selected-count likelihood,
orientation, parameter bounds, two prescribed starts and original convergence
criterion. Preserve the original failed-fit accounting; no dropping failed draws,
shrinking variance bounds or guessing likelihood values after an exception.
Local nonconcavity alone does not invalidate the finite Gaussian integral.
The existing conditional likelihood is bounded above by zero, providing a useful
Gaussian-tail bound for independent integration.

The original prescribed approximation is adaptive Gauss–Hermite with 20 nodes,
with the mode/scale held fixed in its score calculation. A candidate that needs
a different approximation must label that departure and demonstrate accuracy;
do not silently call it the frozen v1.7 implementation. Start with the actual
saved failure points, successful positive-variance points and the exact variance
boundary. Compare likelihoods and gradients with independent scalar integration,
higher-node refinements and finite differences at declared steps. A convergence
flag without those comparisons is insufficient.

Existing specification accuracy requirements (§4): five fixture genes spanning
5–500 rows; compare independent scalar integration at epsabs=1e−12 within
1e−8 relative at tau in {0.25,1,4} times tau_max; compare adaptive20/adaptive40
within 1e−8 on the fixture. Fixed-vector agreement (§11) uses absolute likelihood
tolerance1e−6 and gradient atol1e−8 plus rtol1e−6. The final projected-gradient
criterion remains1e−6, including the right variance-boundary condition.
Report failed numerical comparisons and approximation disagreements explicitly.

## Genuine calibration is a separate experiment

Generate independent datasets, not additional resamples of the one B-FIX
dataset. Known generating parameters are available to the simulation evaluator,
not supplied as fit labels or warm starts. Draw gene random effects before
read counts and het selection so inclusion-induced distortion is actually tested.
Freeze seed streams and scenario assignments before inspecting results; retain
every failed dataset/refit in the specified denominator and report failure rates.

Appendix B3 already requires at least4,000 warp-speed Monte Carlo replicates
for size and coverage bands and200 full-CRT replicates. The latter fails if
the one-sided binomial lower95% bound for size exceeds0.05. Preserve the
approved P2 boundary plants kappa=±Delta*lambda and tau=lambda*tau_max,
Delta=0.125, simulated stage-label error and technical plants: reference-bias
slope, doubled false-het rate at S=2, non-BB tails with depth–stage coupling,
tag switching, LD mixture and omitted depth-dependent bias. Exact scenario
amplitudes and size/coverage bands must be taken from the controlling
specification/owner, not invented here. These numbers are requirements, not
submitted workloads.

Preserve the §4/v1.5 selection fallback: report tau bias and tau-bound coverage
before and after the n* restriction. Absolute bias above0.1*tau_max triggers
the prescribed depth restriction; persistent bias retains the explicit flag
"inclusion-selection not corrected". The revised n* rule uses the99% fraction
of the attainable infinite-depth detection probability and a10,000-depth
ceiling, not an unattainable absolute99% detection target. Failure to find
a qualifying n* must not produce an invented restriction. Preserve §5.2's
tau-component fallback: rejection above0.06 at tau=lambda*tau_max triggers
the conditional parametric bootstrap of the signed root at tau=tau0 using
the §5.3 machinery and1,000 draws.

For every scenario report null rejection rate with binomial uncertainty, bias,
interval coverage, power at the declared alternatives, failed fits and runtime.
As supplemental diagnostics recommended in the prior contribution review,
use matched common-slope reduced models to isolate selection/genotype handling;
an observed-genotype simulation comparator is an oracle with additional input,
not a directly deployable RNA-only comparator. These are not additional
binding requirements in the likelihood specification. Retain the Holm family
{P1,P2,P3} at0.05, with p2=1 when P2 is untestable. Do not use a smaller p-value
as evidence of utility when size is wrong.

The ordinary test draw counts are P1=1,000, P2=2,000 and P3=1,000 (§5), distinct
from simulation replicate counts. Monte Carlo draws must not be extended after
seeing a p-value. Report conditional versus unconditional calibration separately.

## Scale, missing inputs and independent human replication

A small synthetic numerical pilot can check calculations and execution cost;
the60-person/eight-gene B-FIX design cannot certify the actual558-person,
five-cohort development design. Full covariate/genotype-prior adjustment and
information-weighted label attenuation still need the owner-side C5 information
summaries and their exact design bindings. No restricted molecular data should
be exported to the unrestricted tree to resolve those inputs.

The specified independent participant evaluation already exists as protected
sealed F237, including the whole GSE174478 cohort and participant hash parts
of the development sources (§6). Its existence is not permission to open it.
Identity/eligibility, pre-opening power/MDE, numerical agreement and calibration
must satisfy their existing requirements before an authorized opening. Alpha is
refitted rather than transported; replication requires the declared direction
and test, not just a favorable point estimate. Preserve the separate Resource
and experimental follow-up paper scopes.

For a different RNA-to-chromatin or sequence prediction model, that sealed RNA
allelic roster is not a compatible validation endpoint. Those models need their
own source-native measured targets, input/output contract and verified independent
donor joins. This requirement must not be collapsed into the Model A workflow.

## Current execution boundary

No calibration jobs or biological fits have been submitted from this document.
An isolated numerical repair and ready-to-review owner data requests are being
prepared in parallel. Root will inspect the implementation and show the exact
job resources before submitting bounded synthetic checks. Wider simulations
follow demonstrated numerical accuracy and an adequate measured runtime.

## Initial numerical checks submitted

The isolated candidate uses composite adaptive integration of the normalized
prior mass and derivative moments, with log-shape beta-binomial arithmetic.
This is an explicit approximation revision, not silent equivalence to frozen
20-node quadrature. Original numerical modules and completed results are
preserved. Two independent source reviews found no material error in the
Gaussian normalization, variance scores, Hessian covariance terms or tail
allowances. Numerical correctness remains unverified at submission.

Executed copies are frozen in
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/repair_p3_launch_20261005T121740Z/`.
Candidate `repair_p3_integral.py` SHA256:
`461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3`.
Checker `repair_p3_fixed_point_check.py` SHA256:
`89434fd7e32b06cb93d9d31f18b646eaa646736ec329ab48cc00fa4409561f72`.
The initial synthetic case protocol SHA256 is
`8e7ddf8b0a7bb18afd9692b53c3560f4258fdc373c4b822aee36aa94395eabe9`.

| Actual job | Fixed diagnostic case | CPU | Memory | Time |
|---|---|---:|---:|---:|
|22039102|Original failed draw0/gene2 positive variance|1|4GB|90h|
|22039103|Successful draw18/gene2 positive variance|1|4GB|90h|
|22039104|Successful draw19 exact variance-zero boundary|1|4GB|90h|

All are individual cpu/nslab jobs, indexed by actual ID, with32MiB output
caps. Kernel/checker AST, metadata-only admission and launcher syntax checks
passed before submission. Root caught and corrected an unimported test-module
assertion before launching. All three were confirmed RUNNING with empty
startup stderr; this is not a completed numerical check.

Positive-variance checks are initially restricted to the one24-row synthetic
gene. They compare with an independently written full-support mixture
probability and scalar integration, tighter/domain/grid refinements and
declared finite-difference steps. The zero-variance branch compares against
a separately loaded unmodified reference module. Full five-gene/5–500-row
accuracy, complete fit agreement and genuine calibration remain separate
requirements. Candidate refusal is retained, not replaced by an invented fit.

Source inspection corrected the earlier proposed tail-domain failure mechanism:
domain expansion already uses the quadrature-only log lower mass estimate.
[Inference] Intermediate exponentiated relative tail allowances can overflow
for very small masses, but overflow alone does not establish a false refusal;
the subsequent expansion uses the finite log estimate. Deep-row scale checks
remain required before500-row use. No such failure has been observed in the
initial24-row cases. Executing snapshots remain unchanged.

Actual Claude CLI review was attempted with tools disabled and failed with
the weekly-limit message, reported to reset at14:00 America/New_York. Its
response SHA256 is
`361e5b1aac86b74596f2f2707cf11a607e3a434c29562def64bb44ae1128447b`.
This is not a completed Claude scientific review. Local independent review
continues; the ready source-owner requests are in
[the request packet](2026-10-05_codex_independent_data_requests.md).


## Independent reference arithmetic cost comparison

Job22039109 completed in13seconds with exit0. The isolated cumulative
integer-rising-product calculation was compared with the retained literal
full-support formula over210 fixed synthetic depth/eta/rho cases. All were
finite. Maximum absolute log-PMF difference was6.184563972055912e-11;
selected-support conditional log-PMF difference was4.729372449219227e-11.
At depths<=20 the maximum log-PMF difference was3.410605131648481e-13.
Candidate full-support mass differed from1 by at most2.7847057992858026e-11.
Arithmetic timings totaled4.47097472473979seconds literal versus
0.015357142314314842seconds cumulative. These are fixed-formula arithmetic
timings, not an end-to-end integration, fit or calibration speed estimate.

Output: `repair_p3_reference_cost_20261005T124000Z/summary.json` in REC.
Output SHA256: `ad53929ecb40f94add47666f96f333daefd3403d06da3d695bfe011ebe6d190d`.
The working checker now has an optional hash-guarded fast reference; its
literal reference remains available and the executing initial copies remain
unchanged. Any partial derivative shard remains incomplete until the entire
fixed seven-axis/three-step family is inspected. No new calibration or
biological fit has been run.


## Accelerated boundary variance check

Job22039112 completed in25seconds with exit0. Frozen checker SHA256
`b7a6c4c364876e4df71eb69afd57a2fd029fcaf5bce7e43ed696a29e811f8dde`
uses the independently reviewed fast-reference option and only variance axis28.
All three original steps were retained. Its case, draw, fixed gene/rows/p/v,
count hash, environment, kernel settings, outer methods and steps match the
original22039104. Maximum boundary-baseline numerical difference is
1.1951328815484885e-11; full logmass agrees exactly.

For gene2 the independently rederived exact right variance score is
-362.18636151336807. Three-point one-sided scalar differences at steps
.001/.0003/.0001 are -321.0069847614889/-354.74402755416656/
-361.1036035918147, respectively. Differences from the exact score are
41.17937675187915/7.442333959201505/1.082757921553366. These do not
establish derivative accuracy at the required tolerance. The error decreases
with step; truncation and floating arithmetic remain unbounded. Do not label
this completed execution a passed numerical check or replace its steps.
The full seven-axis family and positive-variance checks remain required.
Original22039102–22039104 remain live and unchanged. No fits or calibration.


## High-precision diagnostic environment failure

Job22039117 failed in17seconds with exit1, before numerical evaluation:
`ModuleNotFoundError: No module named mpmath`. Login-node distribution metadata
reported1.3.0, but that did not establish import availability under the job's
`PYTHONNOUSERSITE=1` isolation. Its frozen source and failure record are
preserved. No high-precision derivative result exists from that job.
A standard-library Decimal implementation is being prepared instead; no
dependency installation, environment changes or statistical retuning.


## Independent high-precision boundary score result

Job22039129 completed in1minute44seconds with exit0. All18 literal versus
shared-prefix full/selected-support arithmetic comparisons agree exactly
at80/110digits. At the saved24-row gene2 point, the independent110-digit
right variance score is -362.18636151336717613469575298096352441315355.
It agrees with the candidate analytical score within9e-13. The80/110-digit
finest-step scores differ by1.3527352441315355049011905004698329e-25;
within110digits the1e-20/1e-25 scores differ by about6.56e-36.
Root matched all11 fixed-input/environment/settings fields against22039104.
This independently verifies the boundary score at this synthetic point; it
does not verify positive-variance integration, fit convergence, calibration
or biological validation. Original declared-step discrepancies remain
recorded. Job22039121 remains live and preserved, not replaced or canceled.

Result: `repair_p3_fixed_point_boundary19_22039129/independent_high_precision_boundary.json`
in REC; SHA256 `130a4e0477bcf8367faab885fd591ecaf830218d8bf5b83c6b5cc2ddb3875795`.
The frozen launch is `repair_p3_boundary_decimal_fast_launch_20261005T130300Z`.
Its OpenMP environment variable was spelledOMP_THREADS rather than
OMP_NUM_THREADS; OpenBLAS/MKL were explicitly limited to1. No heavy
OpenMP computation occurs in the standard-library Decimal calculation.
Use the correctOMP_NUM_THREADS spelling for future launchers.


## Predictor cache point comparison

Job22039132 completed in9seconds, exit0. All144 row/offset/Hessian-mode
predictor comparisons and six complete likelihood/gradient/Hessian point
comparisons matched exactly. Changed-parameter refusal and restoration after
normal/exception paths passed. This is one saved failed synthetic point,
not full-fit or numerical-integration equivalence.

Three-point evaluator timings (uncached versus populated cache) were
.038462501019239426/.03590993583202362seconds without Hessian and
.04279948212206364/.041953857988119125seconds with Hessian. Savings
are modest in this single-run measurement; no end-to-end speed claim.
The cache is not incorporated into executing kernels or fits. Initial
integral results remain necessary before another numerical-method revision.

Output: `repair_p3_fixed_point_failed0gene2_22039132/predictor_cache_comparison.json`
in REC, SHA256 `b00393fb46eb7a06e80d5ad8af61a1d4aa5b2d2e671d7b1a2a804e79643cc8a3`.


## Positive-variance adaptive integral refusal

The initial candidate in job22039103 (saved draw18/gene2 positive variance)
refused with `Candidate integrand evaluation cap exhausted` inside
quad_vec subdivision. The100000 integrand-evaluation limit was actually
reached. This candidate has not passed numerical accuracy at this point.
The job remains live for the separate independent scalar/finite-difference
calculations; refusal is not a completed job or a replacement failed CRT draw.

Refusal SHA256: `e547a00f2c9358cb9d1dfef4cb2c6121ae8c6efbbf4948d317c6727f81f05f26`.
A separately declared fixed composite Gauss-Legendre candidate is being
prepared at unchanged mass/gradient/Hessian tolerances, retaining every
finite-grid interval and empirical order/domain/grid refinement comparisons.
Those comparisons will not be described as certified error bounds. No
likelihood/model change is intended; a successful numerical comparison alone
does not establish fit agreement or statistical calibration.


## Composite GL fixed-point submissions

Two independent source reviews found no blocking probability/derivative
formula defect for the bounded checks. Root moved the compute guard before
numerical imports and corrected the evaluation-count label (repeated nodes
are recomputed, not called unique). Frozen GL SHA256:
`45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c`.
Frozen checker SHA256:
`3a7e42b13529f41b0f2bd8eeab17c8f8c5bc66410ba2ccdc734274a59ec29a2c`.
Copies are in `repair_p3_composite_gl_launch_20261005T132600Z` in REC.

Jobs22039200 (failed0gene2) and22039201 (positive18) each request1CPU,
4GB,90hours oncpu/nslab, with32MiB output caps. They use the guarded
fast independent probability reference and baseline-only execution. They
retain empirical16/32/64 order/domain/grid/tail diagnostics; agreement is
not a certified error bound. Independent scalar comparison remains required.
The16→32 comparison is retained but not an acceptance condition.
Full seven-axis derivatives, full five-gene/depth accuracy, fit agreement,
calibration and biological validation are not established by these jobs.
Both initial adaptive positive checks hit the100000-evaluation cap and
remain preserved with their independent calculations still live.


## Original boundary completion and positive18 scalar reference

Original boundary job22039104 completed in1hour25minutes13seconds, exit0.
All21 declared finite-difference records (seven axes, three steps) exist.
Completion does not erase the recorded coarse-step truncation discrepancies;
the separate high-precision boundary score remains the stronger derivative
check at this saved point.

The independent literal scalar calculation in22039103 returned
-69.97130792295839 at z bounds24 and26, with exactlyzero saved
refinement difference. Its reported quadrature log-error estimate is
2.502093446291813e-14, not a certified floating-error bound. The adaptive
candidate had refused, so its candidate-difference field remainsnull.
This saved scalar value can independently assess22039201 when it returns.
Scalar reference SHA256: `d14d4c1e1ef45a06075d7d3e6ea87251704fd3b98715051b8eec5ec58b6c712d`.

## Completed composite likelihood comparisons

Jobs22039200 and22039201 completed with exit0 in38minutes37seconds and
30minutes21seconds, respectively. Maximum resident memory was58724K and
54596K. These are two saved synthetic gene2 points, not model fits.

For failed0gene2, candidate logmass is -87.25979655580069 versus the
independent expanded scalar reference -87.25979655583522, a difference
of3.453237695794087e-11. The previously saved normalized scalar result
is -87.25979655580508. For positive18, candidate logmass is
-69.97130792302882 versus independent -69.97130792295943, a difference
of-6.939160357433138e-11. Both absolute differences are below1e-9.

The positive18 empirical combined refinement allowances are5.2154e-12
for logmass,2.4212e-10 maximum gradient and5.1794e-8 maximum Hessian,
below their unchanged1e-9/1e-8/1e-5 budgets. It used51072 actual node
evaluations; repeated nodes were recomputed. These empirical allowances
are not certified error bounds. Independent gradient/Hessian comparisons,
full five-gene/depth accuracy, fit agreement, calibration and biological
validation remain required. Hash-guarded baseline reuse is being prepared
to avoid repeating these completed baseline calculations in derivative shards.

Independent receipt arithmetic also gives failed0gene2 empirical combined
allowances4.732214620162267e-12 /6.496686522083905e-11 /
1.0312797282807651e-8, with53760 node evaluations. Its full29-coordinate
score differs from the previous direct integration by at most
3.514824542527606e-11. Neither saved point has an independently checked
Hessian. The observed scalar discrepancies exceed the empirical logmass
allowances, reinforcing that these allowances omit floating arithmetic error;
the scalar discrepancies remain below the declared1e-9 point tolerance.

## Parallel Hessian and finite-difference checks with baseline reuse

Independent source review, AST parsing, both case metadata inspections,
14 launcher syntax checks and git diff --check passed. An initial frozen
metadata inspection failed because the stable-BB probe copy was missing;
the exact pinned probe was added before submission and both inspections
then passed. No failed numerical point was substituted.

Jobs22039656–22039669 are14 individual cpu/nslab jobs, each1CPU/4GB/90h,
one per axis2/10/16/17/24/25/28 for each saved positive-variance point.
All14 were observed RUNNING, and all14 logged successful guarded parent
baseline reuse within8seconds. Every perturbed scalar and score integral
is computed afresh. Each shard retains all three declared scalar FD steps
and both finer score-FD Hessian column checks. No single shard satisfies
the full seven-axis family. Fits and calibration remain downstream.

Frozen launch: `repair_p3_hessian_reuse_launch_20261005T142226Z` in REC.
Checker SHA2569ce9ae5ab4fd72312ffdd7e1f88cfd5917bbd25da2ab93e1366575d9bd200dcf;
reuse helper SHA2569de3fc38f954f5b101d1e03a057ce5501943910b3c32fbccbb6dc9ff23a96327.
Four parent receipt hashes per case and actual job IDs are saved there.
The requested1500–1599 name range was historically exhausted; names1520–1533
were reused only after checking for no active scheduler collision. Actual job
IDs uniquely identify these new submissions in the locked job index.

Concurrent complex-step jobs22039439–22039452 are separate independent
first-gradient diagnostics, not Hessian checks. At inspection13/14 had
completed; positive18 retained42/42 records with maximum gradient discrepancy
1.2169821e-10 versus GL across seven axes/three complex steps/two domains.
Available failed0 comparisons differed by at most1.63036e-10. Their flags
do not declare accuracy accepted. Prepared conditional jobs22039592/3
completed144/144 exact point comparisons each, including Hessian modes;
these timings do not establish end-to-end fit speed. Executing frozen copies
were not modified by this baseline reuse change.

## Completed independent first-gradient comparison

All14 complex-step jobs are now completed. Root matched eight fixed input,
parameter and environment fields against each GL parent, verified all seven
axes and the complete three-step/two-domain Cartesian family, zero refusals
and42 successful records per case. Independently rederived maximum absolute
gradient discrepancies versus GL are1.630356682991163e-10 for failed0gene2
and1.2169820706731116e-10 for positive18. Both are below1e-8 at these
two synthetic points. This checks seven selected gradient coordinates;
it does not establish all-coordinate or all-point accuracy, Hessian accuracy,
fit agreement or calibration. Complex truncation and floating errors remain
uncertified. The84 comparisons and input file hashes are saved in
`repair_p3_hessian_reuse_launch_20261005T142226Z/independent_gradient_comparison.json`
under REC. The14 Hessian shards remained live at the scheduler check.

## Next accuracy checks and collection readiness

The bounded stdlib collector `repair_p3_hessian_reuse_collect.py` passed
primary source review and metadata inspection. Its readiness execution found
all14 shard summaries still missing and correctly wrote no completed-results
file. SHA25668b4735aeee377b83c666f902f26a6779f96b1c554ceec7959fd2e87c39a7e0a.
It requires the exact axis/step family and retains 29-coordinate Hessian
differences or explicit refusals, separately from accuracy acceptance.

The original `check_quadrature_v1.py` constructs five synthetic genes with
5/25/75/200/500 rows, seed20260923 and tau_max0.05, evaluated at
tau multipliers0.25/1/4. It tests the original GH integration, not this GL
candidate. Its scalar reference truncates at12 prior SD without a relative
tail check, and its tolerance uses1e-8*(1+abs(loglikelihood)) rather than
strict relative1e-8. Those historical results cannot certify the candidate.
An isolated candidate checker is being implemented on the exact same fixture
and all15 gene/variance combinations, with an independently normalized scalar
reference, explicit tail/refinement diagnostics and strict relative tolerance.
Original GH20/GH40 discrepancies and refusals will be preserved.

Full shared-fixture agreement remains a separate eight-gene/60-person/173-row
requirement. Legacy fixed-vector and fit runners evaluate the original
approximation, so they will not be submitted unchanged as candidate checks.
No full-fit or calibration claim follows from collection readiness.

## Fifteen parallel five-gene accuracy submissions

Final checker SHA2564aec076799a157551d86c1c2fc3da8901e436f8779cca190643c3398fd0f0b21
passed primary and independent source review, AST parsing, all15 launcher
syntax checks and all15 metadata inspections. Nonfinite integration diagnostic
allowances are preserved as explicit tagged values; accepted numerical
likelihood/gradient/Hessian outputs must be finite. This fixes a source-inferred
serialization risk, not an observed accuracy failure. All fixture construction
and RNG calls are identical across shards; only the tested gene/tau point differs.

Submitted jobs22039932–22039944,22039946,22039947 cover exactly fivegenes
times three tau multipliers. Each individual cpu/nslab job requests1CPU,
4GB,90hours; total15CPUs/60GB. All15 were observed RUNNING with empty
stderr immediately after submission. No numerical result is asserted yet.
Frozen launch: `repair_p3_five_gene_accuracy_launch_20261005T143525Z` in REC,
including source hashes, actual submission receipts and separate job logs.
Each single-point summary declares the combined15-point family incomplete;
the full family and all numerical refusals must be considered together.
These likelihood checks do not establish independent Hessian accuracy, whole
fit agreement, statistical calibration or independent biological validation.

## Live accuracy warning retained

Job22039943 (gene3, tau multiplier4) logged NumPy `overflow encountered
in multiply` in an outer-product calculation. At inspection the scheduler
still reported RUNNING, with no final point record or summary. This is an
observed warning, not a completed job failure or an accepted numerical result.
Its origin and effect on the final comparison remain unverified. Preserve
the stderr and evaluate finite-output guards, any explicit refusal and the
independent scalar comparison when the point finishes. No restart or
replacement point was submitted.

## First completed accuracy points and summary Boolean correction

Jobs22039933/34 completed with exit0 in5minutes42seconds and5minutes28seconds.
They are gene0 (five rows) at tau0.05/0.2. Candidate-minus-expanded-scalar
loglikelihood differences are7.105427357601002e-15 and8.881784197001252e-15,
respectively; scalar domain refinements are zero. Original GH20/GH40 also
agree at these two points. This does not establish the other13 cases.

Both raw per-point `candidate_scalar_relative_LL_pass` values are true,
while their requested-point summary values are false. Root identified the
source cause: NumPy Boolean comparison values were tested with Python
`is True` before serialization. The working checker now explicitly converts
each of the three comparison Boolean values to Python bool; AST and
git diff --check passed. New working SHA256:
2bb3993247ab2c276a8ef9feeb71a4729be2fbd331f80052ce242a1a41778b70.
Frozen executing copies and result files remain unchanged. Collection will
retain the raw summaries and distinguish this diagnosed Boolean-summary
discrepancy from independently rederived numerical agreement. No point was
rerun or substituted to correct a summary flag.

## Shared fixed-vector likelihood checks

Checker `repair_p3_bfix_vectors.py`, SHA256
08021479792640cac64ba32770498697d02ea4d6769e18b1285f9aab098e6ba5,
passed primary and independent source review, AST parsing,32 launcher syntax
checks and32 metadata inspections. It preserves the synthetic60-person,
eight-gene,173-row fixture and its four original vectors. Each gene's random
effect is integrated separately. Boundary variance Hessian entries remain
undefined; only the28-coordinate nuisance block is retained at v0.

Jobs22040192–22040223 are32 individual cpu/nslab jobs,1CPU/4GB/90h each,
total32CPUs/128GB. Frozen launch/source hashes/actual receipts are in
`repair_p3_bfix_vectors_launch_20261005T144634Z` under REC. All16 v0 jobs
completed in6seconds with exit0; the16 positive-variance jobs were observed
RUNNING. No completed-job stderr errors were observed.

Root independently summed the eight completed gene contributions for each
boundary vector. Theta1 candidate/scalar loglikelihoods are
-777.262094102323 / -777.2620941021423, difference-1.8064838513964787e-10.
Theta2 values are -550.9440251149389 / -550.944025114964, difference
2.5124791136477143e-11. Both are below the1e-6 fixed-vector likelihood
tolerance. Raw input hashes and sums are saved as
`boundary_vector_likelihoods.json` in the launch. Independent29-coordinate
gradient accuracy, positive-variance vectors, Hessian accuracy and full fits
remain unverified. These are fixed synthetic likelihood checks, not
calibration or biological validation.

## Third completed five-gene accuracy point and compute feasibility

Gene0 at tau0.0125 (job22039932) completed in11minutes58seconds, exit0.
Candidate/scalar likelihoods are -11.6205018788385 / -11.620501878838505,
difference5.329070518200751e-15, with zero scalar domain refinement.
Root rederived the difference and criterion. All three five-row gene0 points
now agree within9e-15; the12 larger-row cases remain required.

The completed saved GL baselines used53760 and51072 evaluations in
2317 and1821 seconds including scalar checks and other work. Every completed
five-row point used37632 GL evaluations. These checker times are not pure
integrator benchmarks. Repeating this cost inside full fits and required
4000 Monte Carlo/200 full-CRT calibration is not yet practical.
An isolated mixture adaptive Gauss-Hermite numerical candidate is being
implemented alongside the preserved GL checks. It retains the same probability
model, normalized prior and fixed-u moments, including a positive prior
proposal and all located peaks to avoid overlapping-peak double counting.
No speed, accuracy or fitted-model result is claimed. It requires independent
source review, full likelihood/derivative/Hessian comparisons and multimodal
synthetic checks before use in fitting. The calibration requirement is unchanged.

## First positive-variance shared-vector points

Jobs22040212/22040220 completed with exit0 in2minutes55seconds and
3minutes44seconds. Both test gene4 with three rows, at theta3 variance0.0064
and theta4 variance0.0004. Candidate-minus-independent-scalar likelihood
differences are1.3322676295501878e-12 and6.359357485052897e-13.
Root independently rederived both differences and the descriptive1e-6/8
per-gene budget comparisons. Neither record contains a numerical refusal.
Only18/32 fixed-vector points are complete at this check; these two small
gene contributions do not verify either full positive-variance vector,
independent gradients/Hessians, optimizer agreement or calibration.

## Bounded mixture GH numerical experiments

Two independent source reviews plus primary inspection found no blocking
mass/moment formula defect for bounded saved-point experiments. Mixture GH
SHA256a00c80741b21f681f268ad24550cca7ea09949c5a7107ef8457b792943885143;
checker SHA2561dfa4d7dd2763e3c898b381fb028b592df6ee6a01734c027e6bbeea5d5e43a7f.
The optional method is pinned, archived and mutually exclusive with GL;
baseline reuse remains GL-only. It preserves all derivative axes/steps and
the same likelihood, prior normalization and fixed-u moments. Nearby peak
proposals are retained; their complete mixture denominator avoids double
counting. Node/search refinement and tails remain empirical, not certified.

Jobs22040408/22040409 are two independent cpu/nslab jobs,1CPU/4GB/90h each,
for failed0gene2/positive18, baseline-only with the separate scalar reference.
Both launcher syntax and metadata inspections passed before submission.
Frozen launch: `repair_p3_mixture_gh_launch_20261005T145831Z` under REC.
No performance, numerical accuracy, fit or calibration result is asserted.
Existing GL jobs and named results remain preserved. Future deep-row checks
must retain nonfinite intermediate diagnostic allowances explicitly.

## Mixture GH peak-search refusals and fixed-vector GL refusal

Both initial mixture GH jobs completed their diagnostic execution with exit0
in55/40seconds, but neither candidate integral returned. Each refused at
`Located mixture peak has unresolved stationary score` before quadrature;
their independent scalar checks completed. These are candidate failures,
not successful numerical agreement or demonstrated speed improvements.
A source-only correction is being prepared to refine sampled likelihood
maxima using the analytic score within the same bounded search window.
The stationary-score threshold is unchanged; unresolved brackets remain
explicit refusals. Preserve22040408/9 and compare any revised candidate
separately rather than replacing the failed records.

Shared-vector job22040217 (theta4/gene1,22rows,v0.0004) completed with a
candidate GL refusal: order refinement exceeded its unchanged precision
allowance. Its independent scalar loglikelihood is -71.86135031599478.
At the last saved bound15/order32-to64 comparison, mass/maximum-score/
maximum-Hessian allowances are3.83693196105488e-13 /
9.821149390193837e-10 /7.822545438379777e-6. The Hessian allowance exceeds
the required quarter-budget2.5e-6 intermediate check. Root rederived these
saved values; this point does not establish candidate accuracy, and all-four
fixed-vector agreement remains unproven. No threshold was loosened and no
replacement point was selected.

## Additional fixed-vector precision refusals

Jobs22040218/19 (theta4 genes2/3,24rows each) completed with the same
GL order-refinement refusal. Their independent scalar likelihoods are
-81.89087065605297 and -78.98638404438505. Root rederived the last saved
order32-to64 mass/maximum-gradient/maximum-Hessian allowances:
5.400124796467048e-13 /1.880438179358893e-9 /8.799589420810567e-6,
and1.9042545318627604e-12 /1.998710538733618e-9 /3.727368045500559e-6.
Each Hessian allowance exceeds the unchanged2.5e-6 intermediate check.
At this inspection21/32 fixed-vector jobs are complete, including three
candidate refusals; a completed diagnostic job is not a passed point.

## Revised mixture GH geometry submissions

Helper SHA2565a44eb71e1cd2adaac4aefa61b71b125542bee053fccc617e9335e83dd9cb390
adds analytic-score root refinement within each sampled maximum's original
window, requiring positive-to-negative endpoint scores. Stationarity and
precision thresholds remain unchanged; missing brackets and unresolved
curvature produce explicit diagnostic refusals. Raw and refined locations,
scores, curvatures and bracket endpoints are retained. Primary and independent
source review, AST, launcher syntax and metadata inspections passed.

Jobs22040508/22040509 test the revised candidate at the same two saved points,
each individual cpu/nslab with1CPU/4GB/90h, baseline-only. Frozen launch:
`repair_p3_mixture_gh_refined_launch_20261005T150445Z` in REC. The previous
two peak-search refusals remain preserved. The proposed explanation of those
refusals is insufficient value-optimizer stopping precision [Inference]; the
correction's effect remains numerically unverified until these jobs return.

## Revised positive18 mixture GH order refusal

Job22040509 completed in2minutes52seconds, exit0, with a candidate refusal
after148.99seconds: `Mixture GH order refinement exceeds unchanged budgets`.
Peak refinement passed this point; the candidate still did not pass integration
accuracy. Its3754 evaluations retained eight order records across bounds8/15.
At bound15, GH40/80/160/320 loglikelihoods are -69.97101165494794,
-69.97130223668464, -69.9713080940409 and -69.97130792376944.
The final160-to320 mass/maximum-score/maximum-Hessian allowances are
1.7027144849636873e-7 /1.2471998662565502e-5 /0.006350641960268671,
all above the unchanged intermediate budgets. Independent scalar checking
completed. This is not evidence of an accurate faster integrator.
Failed0 revised job22040508 remained live at inspection.

The working checker now also preserves `geometry_reports` on candidate
refusal, for future proposal diagnostics. Existing frozen jobs/receipts remain
unchanged; geometry details absent from22040509 cannot be reconstructed as
observed evidence without another numerical check.

## Additional positive-variance vector likelihood comparisons

Jobs22040222/23 (theta4 genes6/7,21/24rows) completed with candidate and
independent scalar values present, without numerical refusals. Their
candidate/scalar loglikelihoods are -57.38247182720901 /
-57.38247182720831 and -74.7596141308477 / -74.75961413084731.
Root rederived differences-6.963318810448982e-13 and
-3.836930773104541e-13, both below the descriptive1.25e-7 per-gene budget.
These additions do not repair theta4's three retained candidate precision
refusals or establish full-vector, gradient or Hessian agreement.

## First completed25-row accuracy case

Job22039936 (gene1,25rows,tau0.05) completed with no candidate or scalar
refusal. Candidate-minus-independent-scalar loglikelihood difference is
1.4210854715202004e-14; root independently rederived it and the original
relative-loglikelihood criterion. Original GH20-minus-GH40 is
-1.4210854715202004e-14. This brings completed accuracy cases to4/15;
the remaining row/variance combinations are not inferred to pass.

The revised positive18 mixture refusal has negligible final mass-tail
allowance (about1.8e-20), and each order's likelihood is identical at bounds8
and15. Thus its saved order disagreement is not resolved by domain expansion.
The broad prior proposal is a plausible convergence limiter [Inference];
the current receipts do not isolate individual component contributions.
Future diagnostic comparisons will preserve proposal widths/weights and
component-specific moments. Two fixed numerical weight policies are being
prepared for bounded tests; these modify the quadrature proposal, not the
likelihood, derivative targets, calibration requirements or acceptance budgets.

## Four numerical proposal-policy checks

Helper SHA25668f8311749795c145ce836bc56876255ba1ad967af04df2342a19eb121cae24b
and checker SHA25697687ec12978a66a4dc0ded0299a5ea0b326696716a8a1b896a5ddb59ef7e027
passed primary/independent source review, AST, four launcher syntax and four
metadata inspections. They preserve positive prior support, full mixture
normalization, all located peak proposals, node orders and precision budgets.
Invalid or underflowed weights/contributions refuse explicitly. Component
mass, score, conditional Hessian and centered outer-product contributions
are recorded separately. No component-specific convergence claim exists yet.

With J peak proposals and s=min(1,minimum standardized peak SD), scale_sqrt
sets prior weight min(1/(J+1),s*sqrt(logmass_tolerance)). Scale_precision sets
min(1/(J+1),s*epsilon), epsilon=.01*min(logmass_tolerance,
gradient_tolerance/(1+max(cp)),hessian_tolerance/(1+max(hp)+max(cp)^2)).
The Hessian term is omitted for score-only evaluations. Thus proposal weights
may differ by Hessian mode; independent score/FD checks remain required.
Each peak receives the remaining weight divided by J; with no located peaks
the prior alone has weight1. These are prespecified numerical hypotheses,
not model parameters or certified accuracy bounds.

Jobs22040801–22040804 test both policies at the same two saved points,
each individual cpu/nslab job with1CPU/4GB/90h. Frozen launch:
`repair_p3_mixture_gh_policies_launch_20261005T151726Z` under REC.
All previous equal-weight and peak-search refusals remain preserved.
Neither policy is adopted for fitting, calibration or biological analysis.

## First precision-policy candidate return and new completed GL checks

Positive18 scale_precision job22040803 returned its candidate integral at
235.76seconds with5738 likelihood evaluations. Root matched all11 saved
input/environment/settings fields against the GL baseline22039201. Candidate
loglikelihood is -69.97130792301192. Candidate-minus-GL is
1.6896706256375182e-11; maximum29-score and29x29-Hessian differences are
4.044249379830944e-10 and1.252683432539925e-7. The GL baseline used51072
evaluations. These are one-point comparison/count observations, not a
full-fit speed or independent Hessian-accuracy claim. The separate scalar
check was still running at this inspection. Final prior proposal weight is
8.526815342805169e-17, positive; all components remain represented.

Positive18 scale_sqrt job22040804 completed with a retained order-refinement
refusal after3754 evaluations. Its final prior weight is1.2721020023811323e-5.
Mass/maximum-score/maximum-Hessian allowances are1.568025709293335e-9 /
2.5237701185461765e-7 /4.3277829861447856e-5, above intermediate budgets.
Candidate LL at order320 is -69.97130792282297, but this alone does not
establish gradient/Hessian accuracy or justify accepting the refused point.

GL accuracy job22039935 (gene1,25rows,tau0.0125) returned candidate and
scalar LL -66.2473467243169, saved difference0. Theta3/gene1/6/7 fixed-vector
jobs22040209/14/15 have candidate-minus-scalar differences
-4.547473508864641e-13 /7.460698725481052e-13 /-1.2221335055073723e-12.
Root rederived these four comparisons. At this check5/15 accuracy cases
and26/32 fixed-vector diagnostics are complete; the latter include three
retained candidate refusals. Independent full-family derivative/Hessian,
fit and calibration requirements remain unchanged.

Positive18 precision-policy job22040803 subsequently completed in4minutes
22seconds, exit0. Independent scalar LL is -69.97130792295943 with zero
domain refinement; candidate-minus-scalar is -5.24948973179562e-11,
below the1e-9 fixed-point likelihood tolerance. This strengthens the saved
likelihood comparison at one synthetic point; gradients/Hessians still
require the declared independent comparisons and broader fixtures.

## Complete four-policy comparison

All four policy diagnostics are now complete. One candidate integral returned
(positive18/scale_precision); the other three refused. Root reviewed and ran
the bounded stdlib collector, adding reconciliation for retained quadrature
orders in refused points as well as successful ones. The complete result is
`repair_p3_mixture_gh_policies_launch_20261005T151726Z/mixture_policy_collection.json`
under REC. It preserves all failures and explicitly declines automatic policy
selection, derivative accuracy, full-fit agreement and calibration claims.

Across36 saved quadrature orders, component mass fractions reconstruct1
within2.22e-16, scores within5.69e-14 and Hessians within1.82e-12;
proposal-weight recipe residuals are zero. These are arithmetic checks, not
independent numerical-accuracy bounds. Scale_precision epsilon uses saved
metadata and has not been independently rederived from cp/hp arrays.

At failed0/scale_precision, final160-to320 allowances are2.5170265881632655e-10,
2.10132451558467e-8 and1.6639627634164806e-6. The score allowance exceeds
its quarter-budget. The order320 LL -87.25979656172582 also differs from
the independent scalar by about5.89e-9, exceeding1e-9. This makes a passing
refinement comparison alone insufficient. Three distinct peak locations were
retained, with duplicate proposals; the prior component's normalized mass
contribution is2.615380368785254e-10 despite proposal weight8.01944122297595e-18.
Both failed-point policies remain unvalidated; no threshold was relaxed.

## Completed fixed-vector likelihood collection

All32 BFIX diagnostics completed, including the three retained theta4 candidate
refusals. The reviewed stdlib collector wrote
`repair_p3_bfix_vectors_launch_20261005T144634Z/bfix_vectors_collection.json`
under REC. Eight-gene candidate-minus-independent likelihood differences are
-1.8064838513964787e-10 for theta1, 2.5124791136477143e-11 for theta2,
and3.467448550509289e-11 for theta3, each below1e-6. Theta4 has no valid
eight-gene candidate sum because genes1/2/3 refused Hessian refinement.
All-four-vector likelihood agreement is therefore false. Completion does not
certify gradients, Hessians, fitted parameters or statistical calibration.

The failed0 precision-mixture component review found a central-pair logmass
change of only6.04e-12 between GH160/320, while the prior component oscillated
across orders. Its independent total likelihood error remains about5.90e-9.
[Inference] Wider local proposal components could better resolve the central
shoulders sampled poorly by the broad prior. The next bounded numerical option
adds a Gaussian of exactly four times the local width at every located mode,
retaining the prior and original proposals. All policies include every component
in the complete denominator and share nonprior weight equally. Original
orders, caps and accuracy thresholds remain unchanged. This option requires
fresh independent comparisons and is not adopted for fits or calibration.

Wider-component jobs22041642–45 test failed0/positive18 with each of
scale_precision/scale_sqrt. Frozen launch:
`repair_p3_mixture_gh_wide_launch_20261005T154158Z` under REC; checker SHA
c4edba988955b8b922e28a227c4b936061b073df3b2995c60f14785d5565413c;
helper SHA8a21be3a5d4f7c731b13cdfe5a01165f640be36556fc8227ec6a8c67441828a0.
Two source reviews found no blocker for bounded testing. AST/compilation,
checker metadata inspection, shell syntax and diff checks passed. Each
individual job uses cpu/nslab,1CPU/4GB/90h. Numerical results remain pending.

Read-only follow-up localized theta4 genes1/2/3 refusals to H[v,v]:
order32-to64 differences7.812727e-6,8.799540e-6 and3.727365e-6 exceed
the2.5e-6 initial budget. Likelihood and gradient refinements pass and expanded
tail allowances are small. [Inference] Cancellation in small-variance moment
subtraction could contribute; this is not confirmed. Concurrent local-GL
jobs22041594/95 and synthetic observed-fit job22040873 were observed RUNNING;
they are separate experiments and were not duplicated here. Await their
results before targeting only these three failed genes. No tolerance change
or model adoption follows from the localization.

## Wider-proposal results and reconciliation

Root reviewed `repair_p3_mixture_wide_collect.py`, SHA
bdedc1d261f284f7c66836f38df1038aad9b9baee6fdeb4f8b5180d8ffd6bdc4.
Its metadata inspection passed; readiness correctly withheld output while
failed0 jobs22041642/43 remained live. It verifies all four jobs, original
and fourfold-width components, component moments and frozen GL/narrow
comparisons, preserving every refusal. Result writing requires its own SHA.

Positive18 wider jobs22041644/45 completed with candidate precision refusals.
At order320/domain15, likelihood errors versus the independent scalar are
-9.95191115293892e-7 and-9.950317689799704e-7, both above1e-9.
Final order-refinement LL/maximum-gradient/maximum-Hessian allowances are
6.470124381507959e-5 /0.0044444689872398906 /0.5551742433648396 and
6.469682524823385e-5 /0.004444666838100918 /0.5551671153374481.
Each used6154 evaluations. Wider proposals are not validated at this point;
no passing refinement or speed result is inferred from job completion.

Concurrent local-GL jobs22041594/95 completed, exit0, in5:58/4:20.
Acceptance, source identity and numerical agreement remain under review here.
Their completion alone does not establish accuracy or fit agreement.

The read-only local-GL comparison now confirms both integrals returned.
All eight fixed identity fields/count hashes/environment and fifteen guarded
source files matched. Final LL equals original GL bitwise at both points;
maximum gradient differences are0/3.46e-17 and Hessian differences
5.11e-15/3.55e-15. Evaluation counts remain53,760/51,072. Reported integration
times are176.887/128.387seconds; original completion elapsed includes startup
and is a different timing scope. Dimension-expanded runs still use the same
24 synthetic rows, so these do not establish full-scale or full-fit speed.

Overall comparison flags remain false: huge tail allowances at rejected
initial domains differ in absolute floating arithmetic beyond the absolute
accuracy budgets. Accepted final moments and refinements agree. Preserve the
false flags and the rejected-domain records; this is not independent Hessian
validation or a completed full-family comparison.

At theta4 v=.0004, the raw variance second-moment polynomial has a1/v^2
factor of6.25e6. [Inference] Likelihood arithmetic errors near1e-12 could
therefore contribute to H[v,v] differences near1e-5. Existing reports lack
separate conditional-H and score-outer terms and cannot prove this cause.
A targeted matched-domain order128 diagnostic for only the three failed
genes is being prepared. It retains16/32/64 comparisons and unchanged
budgets/caps; it does not claim a completed extra full-domain ladder, which
could exceed the100,000-evaluation cap.

All four wider-proposal checks have now completed and refused. Root ran the
reviewed collector, producing
`repair_p3_mixture_gh_wide_launch_20261005T154158Z/mixture_wide_collection.json`.
Failed0 precision/square-root order320 likelihood errors are
5.3253966768807004e-5 /5.3178265218889464e-5; both exceed1e-9.
Each used15,817 evaluations, elapsed773.24/739.47seconds. All four outcomes,
preceding narrow refusals and component arithmetic remain preserved. This
wider-proposal experiment is negative and is not selected for fitting.

Targeted order diagnostics22041792–94 now run theta4 genes1/2/3, respectively.
Frozen launch `repair_p3_theta4_order_launch_20261005T155906Z` under REC.
Checker SHA549cd4a46901eb9473331ebbac7ebd989052c6c660ec5604a84f672670f146d4;
helper SHA5458e862cdd7871e4c392451bf5e41f1218784f7b5529a37d6ac8f67006110a5.
Root reviewed the copied point/Kahan/tail formula against the original GL
source; AST, metadata inspections and shell syntax passed. Domains15/16/16,
spacing.25 and orders16/32/64/128 require28,800/30,720/30,720 evaluations.
Each individual job uses cpu/nslab,1CPU/4GB/90h. Original refusals and all
precision budgets remain; this is a diagnostic, not full integral acceptance.

The live synthetic observed fit22040873 has completed its null stage;
both starts converged and warm was retained under the original tie rule.
The marginal stage is live; two positive objective calls took617.8/557.8s.
Warm/cold sharding is being prepared for future fits only. It requires the
same pinned shared null, original starts, all failures, original winner/tie
rules and subsequent sequential boundary steps. Current work is not restarted.

The future synthetic observed-start prototype is now source-reviewed:
`repair_p3_observed_start_shard.py` SHA
7ad2577a9fa74a4ff83aea4b71ede577c67809b9dd88abc0ac49d890278ab232 and
`repair_p3_observed_start_collect.py` SHA
86d4e123d683f11b33e65fa85f1101db96eaccfb14933e638274b0eb0c2dccc9.
It uses ordinary imports and explicit setup, not runtime source rewriting.
Root verified the original warm/cold definitions, kappa attempt order,
unfiltered likelihood-based selection and1e-6 tie rule against the archived
fitter. Both inspection modes passed with the saved completed shared null.
No shard jobs or numerical tests were run. This prototype is restricted to
the pinned60-person/173-row fixture and converged shared null; it does not
implement general calibration datasets. Boundary restart, null snap and final
diagnostics remain explicitly deferred, so its collection cannot constitute
a completed fit. Current observed fit22040873 continues unchanged.

Root independently rederived null cold-minus-warm LL7.23048287909478e-11,
below1e-6, confirming warm retention and the saved convergence flags.

## Independent positive-variance arithmetic reference

Three jobs22042086–88 now check only theta4 genes1/2/3 variance scores and
H[v,v] using the native selected genotype-mixture probability in Decimal.
Frozen launch `repair_p3_decimal_positive_launch_20261005T161736Z` under REC.
Checker SHAafa3ef54af7f66ea3c651c85bbc1234b91a78fb8593c3f271573c7a2d15e559f;
helper SHA6ff2857ce42e408478f9c5dc7cc13ed6c70bcfdd77bd8f2169fff8de9656ac40.
The previous Decimal boundary helper remains unchanged and hash-pinned.

Both80/110 digits use all20/40/80 normalized Gaussian rules. Double nodes
initialize Newton iteration only; accepted roots, weights, probabilities and
moments use Decimal. Probabilists' Hermite weights are
(n-1)!/[n He_(n-1)(z)^2]. Fixed-u prior score is H2(z)/(2v), and raw variance
second derivative is H4(z)/(4v^2); log-Hessian subtracts the squared mean
score. Positivity, symmetry, root correction and Gaussian moments0–8 are
checked. All order/digit results and tail diagnostics remain visible.

Two source reviews, AST/compilation, all frozen metadata inspections and
shell syntax passed. Root added the explicit immutable fixture seed to the
protocol. Each job1CPU/4GB/90h, estimated0.5–2h and280 probability evaluations.
There is no new dependency, actual genotype read, fit, or accuracy acceptance.
Tail diagnostics divide by an estimated mass and do not bound interior GH
error. Full Hessian, full-model fit and statistical requirements remain intact.

The metadata-only owner request now explicitly requires development F558
identity qualification as well as F237 readiness. Current authorities do not
verify genetic clusters as independent people. Historical roster totals cannot
support actual power or independent validation until their eligibility and
person/library/kinship grouping are reconciled. No protected roster was read.

All three matched-domain order diagnostics22041792–94 completed. Root checked
the exact jobs, source/input hashes, p/v/row identity, environment, all29-vector/
29-matrix dimensions, and independently rederived every adjacent-order
difference and allowance. Collection:
`repair_p3_theta4_order_launch_20261005T155906Z/order_diagnostic_collection.json`.
Order64-to128 maximum Hessian allowances are4.673729716436969e-6,
3.2099634498407755e-6 and4.128591178440335e-6, all above2.5e-6.
Likelihood and score allowances remain below their corresponding budgets.
Increasing to128 alone therefore does not resolve these three refusals.
The high-precision references remain necessary to distinguish arithmetic
amplification from unresolved quadrature error; neither cause is yet proven.

The three Decimal jobs started with the exact BFIX seed20260930 and common
full-count SHA7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad.
Rows22/24/24 and variance.0004 match the three targeted points. Initial root
rules and point files are present without failure; full order/digit agreement
is still pending. This is distinct from the five-gene accuracy seed20260923.

A separate, unadopted native BB arithmetic prototype,
`repair_p3_bb_mode_anchor.py` SHA
a31258aaac0fb39e91dda7aa634f0e810a6d0da144511b34d9c1b73a7028c266,
computes full-support probabilities from adjacent ratios about a mode,
normalizes the native law, and subtracts normalization score/Hessian moments.
Two source reviews found no mathematical blocker. [Inference] It may reduce
large-log subtraction error, but no numerical check or accuracy/speed benefit
is established. Its standalone scope is n<=100000; estimated peak allocation
is48MiB excluding interpreter overhead. A bounded independent Decimal
probability/derivative regression is being prepared. Existing kernels,
conditionals, integrators and running fits remain unchanged.

## First completed positive Decimal references and measured gene parallelism

Decimal genes1/3 jobs22042086/88 completed in7:45/10:34; gene2 remains live.
Gene1 final H[v,v] is -669.1628061535175323778677314021881818912986305...
Order40-to80 differs by2.32e-65 and80-to110 digits by9.07e-70 in absolute
H[v,v]. Gene3 final value is -3196.75235366888095684002380198879361661888896...
Order40-to80 differs by7.17e-68 and final digit comparison by5.62e-69.
These are empirical refinements, not certified interior error bounds.

Root independently subtracted gene1 reference from every saved double GL
order: errors16/32/64/128 are -2.923881e-6,+5.056349e-6,-2.756378e-6,
+1.907533e-6. [Inference] This strengthens the arithmetic-amplification
explanation for the nonmonotone order differences. Full29-coordinate Hessian,
fit and statistical checks still remain required; no reference result alone
accepts the integration method.

Concurrent eight-gene comparison22041812 completed. Root rederived the timing
ratio727.9817257486284/160.83139715343714=4.526365738488958 and verified exact
LL/gradient/Hessian array equality from both raw outputs. This is a measured
speedup at one synthetic fixed point, not a whole-fit or biological-scale
performance guarantee. Parallel observed-fit job22042114 is now separately
listed in the shared ledger; the existing sequential fit22040873 remains live.

Mode-anchored BB regression22042221 is submitted,1CPU/4GB/90h. Frozen launch
`repair_p3_bb_mode_anchor_launch_20261005T163221Z` under REC; checker SHA
559eb29393f1165dc77608b4d54845ca3bad624ae4697f48dedcd75c2f12f684.
Root added end-of-run source/input hash checks, then independently reviewed
the finite-difference formulas and scoped criterion. Two source reviews,
AST/compilation, frozen metadata inspection and shell syntax passed. The
18-point Decimal logPMF/g/H family is the independent accuracy comparison;
28 broader cases and210 conditional-row comparisons remain descriptive.
Maximum representative depths2426/2466/3000 require an upper bound2,597,508
Decimal prefix logs. The method remains unadopted pending these results.

## Completed selected BB regression and positive Decimal references

Mode-anchored BB job22042221 completed with exit0 in6:16. All18 selected
Decimal logPMF/score/Hessian comparisons passed their prescribed criteria.
Root and an independent reviewer rederived every saved difference and ladder
stability value from raw records. Maximum absolute errors were5.4905e-15 for
logPMF,1.2510e-14 for scores and2.5154e-14 for Hessian entries; the maximum
Decimal step/precision difference was5.6731e-31. All17 source/input pins
match. This verifies the18 selected BB points only, not arbitrary tails,
conditional mixtures, integrals, full model fits or statistical calibration.

All three positive Decimal jobs22042086–88 now have completed six-record
summaries. Gene2 completed with exit0 in20:21. Its final110-digit/order80
H[v,v] is -1862.1566154713177483064059233192545395626186654...
The order40-to80 H[v,v] difference is1.2478e-64 in absolute value, and the
80-to110 digit difference is1.6578e-69. These are empirical refinement
agreements; the interior quadrature error is not certified. All four new
jobs have empty stderr. The previous three theta4 GL refusals remain intact.

Matched positive-variance calls3–8 from the live sequential22040873 and
parallel22042114 fits have exactly equal saved variances, likelihoods and
maximum score summaries. Root rederived speed ratios4.5363,4.2043,4.1745,
4.4323,4.5729 and4.5940. These are observed evaluation timings; neither fit
has a completed marginal result. They do not establish equality of all final
parameters or statistical failure behavior.

Three isolated anchored-BB/original-GL checks are being prepared. Their
admission requires the completed selected BB regression and all three
Decimal summaries, with byte pins for those records and the original
theta4 refusals. The unchanged full GL16/32/64 domain/grid/refinement rules
must pass before this calculation can be considered for further fitting.
Live SLURM accounting temporarily failed to connect to its database;
`scontrol show job` independently confirmed the completed job states.

Root reviewed and ran the stdlib Decimal collector, SHA
1b58b79814711d27c25326550cbb9e9fb2ee9b62e51695d5cb6922031c1887a1.
Metadata inspection and collection passed; all21 saved order/digit refinements
were independently recomputed. Output:
`repair_p3_decimal_positive_launch_20261005T161736Z/decimal_positive_collection.json`.
Every double GL order16/32/64/128 is compared with110-digit/order80 Decimal.
H[v,v] errors across those orders are respectively:

| Gene |16|32|64|128|
|---|---:|---:|---:|---:|
|1|-2.923881e-6|+5.056349e-6|-2.756378e-6|+1.907533e-6|
|2|+5.528580e-6|-4.744249e-6|+4.055291e-6|+8.453771e-7|
|3|-1.997671e-6|-2.317331e-6|+1.410035e-6|-2.718554e-6|

All accuracy-acceptance flags remain false, and original refused calculations
are included. The comparison supports testing the native BB arithmetic
revision; it does not certify the full integral or calibrate any hypothesis.

Full-integral checks22042346–48 are now submitted for genes1/2/3,
1CPU/4GB/90h each, cpu/nslab. Frozen launch:
`repair_p3_theta4_anchored_gl_launch_20261005T164218Z` under REC.
Checker SHA041400eeaf1706d939f3f01caa8bd48f568fc55c23c349de8b780628ce105aee;
14-receipt manifest SHA
cbb2feffcaf9f7775209c6161f879df8bc00e000c3998400b2cbd6b2caa39528.
Root reviewed the isolated original conditional-helper patch, unchanged
integration settings, completed-job checks and preservation of original
refusal bytes. AST, all three frozen metadata inspections (59 pins each),
shell syntax and live prerequisite exit-zero checks passed. Submissions are
recorded in the shared locked/fsynced job ledger. There are no new dependencies,
actual-outcome reads, model adoption or figure/release changes.

All three new jobs started and wrote their compute protocols. Root checked
rows22/24/24, variance.0004, the exact synthetic count SHA, all four live
COMPLETED exit-zero admission records and empty stderr. Their current measured
MaxRSS values are54032/54036/54064K, below each4GB allocation. The full-integral
results remain pending; no threshold has been relaxed or result accepted.

## Full29-coordinate score comparison and observed-fit refusal

The completed full-gradient family22041901–29 is already collected at
`repair_p3_full_gradient_collection_20261005T164404177840Z/collection.json`
under REC. Root preserved the earlier partial collection and used the
complete existing result without executing a duplicate collector. Source
review of `full_gradient_collect.py` found the saved-ratio reconstruction
missing from its arithmetic checks. Root separately rederived all720
per-gene probability derivative/mass ratios and174 model-score sums from raw
saved numerators and pinned masses, verified24 source/metadata pins and every
record hash, and checked the nominal1e-8+1e-6*abs(reference) comparison.

All174 nominal comparisons passed. Maximum absolute score difference was
1.7053380929610285e-10, at most0.01705325 of its nominal tolerance. The
720-record family includes all29 coordinates, three complex steps and two
domains at the same saved positive18 point,173 synthetic rows/eight genes.
The model score is the sum of per-gene derivative/mass ratios, not a ratio
of pooled numerators and masses, because gene random effects are independent.
Independent arithmetic report:
`repair_p3_full_gradient_ratio_recheck_20261005T165339Z/summary.json`.
Complete collection SHA
ba7d2e6b7d41cfdcb7bdb362b808158a1bd66a7bee30fc0a8505529a6a519c5a.
All accuracy/fit/calibration flags remain false:69 comparisons lie outside
their uncertified error estimates and256 step/domain pairs lie outside
their combined estimates. These estimates do not bound floating arithmetic
or complex-step truncation. Whole-Hessian and kept-fit accuracy remain open.

Parallel observed fit22042114 FAILED with exit1 after28:11, confirmed by
live `scontrol show job`. It refused gene2 at marginal call12,
variance0.0009636462009520182, without a Hessian request. The expanded-domain
32-to64 variance-score disagreement2.7571399244333766e-9 exceeded the
unchanged initial-quarter allowance2.5e-9; the maximum corresponding tail
allowance was6.7360e-15. This was a score-precision refusal, not a pool,
wall-time, memory or Hessian failure. Seven other genes returned successfully;
no full likelihood sum was formed. The inherited warm-fit exception stopped
before the outer cold fit/restart. The exact28 parameters and all eight
worker records are preserved, so an isolated point test needs no fit rerun.

Failed packet SHA
d3e7946dc977906052a34a36a12fc9898ba0d82fd77497a52df41f83f8a7f305;
protocol SHA05c69cdb3aa293a7b1150bfc16c0a28487d76d731307236341e4e2c0de5b0604.
The point uses the original synthetic observed count SHA7009be4... and
fixture seed20260930, not the positive18 draw's count vector. Existing
sequential22040873 remains live. No replacement fit was submitted, failed
draw omitted, precision relaxed or accepted fit claimed. A source-only
isolated test of the exact refused point is being prepared alongside the
already running theta4 arithmetic tests.

The existing sequential fit22040873 also FAILED, exit1, after1:33:30.
Root checked that it refused exactly the same call12/gene2/p28/variance as
parallel22042114. All six refused order/domain likelihoods and all29 scores
were identical, not merely within a tolerance. Comparison:
`repair_p3_observed_failure_comparison_20261005T170102Z/summary.json` under
REC. Sequential/parallel time-to-this-failure was3.31858; this is not an
accepted-fit speed measurement. It establishes matching behavior for this
one numerical refusal, not all exception classes or calibration draws.

New exact-point diagnostic22042612 is running,1CPU/4GB/90h, cpu/nslab.
Launch `repair_p3_observed_failure_anchor_launch_20261005T165931Z` under
REC; checker SHA558dd4ffe99924349f1b0523878a6d576411c75052c7d519607f50fc5d66740a.
Two source reviews, AST, frozen metadata inspection (55 pins) and shell syntax
passed. Compute protocol confirms the saved variance, original count SHA,
seed20260930 and24 target rows. It first checks three selected maximum-depth
BB points at u=0,+sqrt(v),-sqrt(v) against native probability-only Decimal
finite differences. Only if these pass does it run fresh unchanged original
score-only GL plus independent80/110-digit,20/40/80-order Decimal quadrature.
The previous full failure packet is copied without modification. This tests
one saved point, not the full parallel callback or a new fit; results remain
diagnostic. A separate source-only batch implementation of the same anchored
BB law is being prepared to retain node batching if arithmetic tests succeed.

## Exact alternative variance-score identity, source derivation only

For fixed parameters and one gene's native selected conditional likelihood
L(u)=exp(ell(u)), Gaussian differentiation gives
d(phi_v)/dv=(1/2)d2(phi_v)/du2. Two integrations by parts give the exact
whole-line identity

    d(log Z)/dv = (1/2) E_posterior[ell_uu + ell_u^2].

The native predictor is affine in u with slope S_i. Thus ell_u is the sum
of S_i*s_eta,i, and ell_uu the sum of S_i^2*h_etaeta,i. The original
conditional derivatives include genotype mixtures, selection normalizers,
orientation and reference bias; no additional orientation sign is introduced.
Counts, priors, supports, covariates, dispersion and other parameters are held
fixed. The right boundary score is(ell_uu(0)+ell_u(0)^2)/2, with inverse-variance
units. This matches the prior-score form used by the current integrator.

For fixed physical finite endpoints a,b the prior-score identity additionally
has [L*phi_v' - L'*phi_v]_a^b/(2*Z_[a,b]). Standardized quadrature limits move
in physical u as v changes, so their endpoint-motion terms cannot be ignored.
Smooth finite native supports and positive dispersion make the whole-line
boundary terms vanish under Gaussian decay. Parameter scores stay
E_posterior[ell_p]. A heat-equation representation for a cancellation-resistant
full Hessian would additionally need fourth u derivatives for H[v,v] and mixed
third derivatives for H[p,v]; the current prior-moment representation computes
those entries with fewer conditional derivatives but retains subtraction risk.

Root and an independent reviewer derived this identity directly from the
frozen native predictor and conditional law. No implementation, numerical
acceptance, new scientific hypothesis, tolerance change or method adoption
follows. It is a mathematically equivalent score fallback only if the current
arithmetic revision still fails near v=0; it does not complete the Hessian.

## Two completed anchored integrals and node-batch comparisons

Theta4 genes1/2 jobs22042346/47 completed without refusal. Gene2 completed
in26:31 and gene1's reported computation took2449.84s. Root checked all59
source/input pins per point and independently rederived all six Decimal
comparisons. Final110-digit/order80 comparisons are:

| Gene | GL minus Decimal logmass | Variance score | H[v,v] |
|---|---:|---:|---:|
|1|+5.16264e-14|-4.67470e-13|+5.15580e-10|
|2|+1.19954e-16|+2.65616e-13|-3.24850e-9|

Gene2's full GL refinement allowances independently rederive to mass7.15e-37,
maximum score2.13e-12 and maximum Hessian8.30e-9, within unchanged thresholds.
This resolves that point's previous refinement refusal; independent full
Hessian accuracy is not established. Reviews are saved as review_gene1.json
and review_gene2.json in the anchored launch. Gene3 and the exact refused-fit
point22042612 remain running. The latter's three selected BB points passed;
its GL/Decimal integral results remain pending.

The source-only node-batch helper was completed before its delegated agent
hit a usage limit. Root finished the source review locally and wrote
`repair_p3_bb_mode_anchor_batch_check.py`. New helper
`repair_p3_bb_mode_anchor_batch.py` retains per-node anchors, full-support
normalization and compensated moments. Directional cumulative sums start at
each anchor; no global-prefix subtraction is introduced. It caps32 nodes and
100000 full-support/requested entries, with an explicit conservative array
estimate below256MiB. These are synthetic diagnostic compute bounds, not
scientific eligibility filters.

First batch comparison22043044 completed with21/21 selected Decimal criteria
met and bitwise-equal batch/scalar timing-case arrays. It gave no clear speed
benefit: scalar/batch ratios1.0085,0.9331,1.0107. Root preserved those results.
The only following arithmetic-path edit changed per-entry Python conversion
to ndarray.tolist before fsum, retaining the float64 values and summation order.
Follow-up22043255 also passed21/21, with bitwise-identical timing-case LP/g/H.
Ratios1.0841,1.0884,1.0827 indicate a modest component improvement; this is
not an integration or whole-fit speed measurement. Measured peak RSS67272KiB
was below4GB allocation and the declared array-capacity limit.

Root independently rebuilt all21 saved differences and Decimal ladder
stabilities for both jobs. Maximum LP/score/Hessian errors were1.31724e-14,
7.53995e-13 and1.71639e-12; maximum ladder difference2.24338e-29. The family
comprises the original18 points and three at the actual saved failed synthetic
fit point, not21 biological participants. Input cases were identical across
both runs. Detailed independent review is at
`repair_p3_bb_mode_anchor_batch_fast_launch_20261005T172617Z/independent_selected_review.json`.

Current batch helper SHA
ea6e3c2627100e6a044732cf666d07096b4972e78cbeec08036052399ca25248;
checker SHA20a54d4addbbdd813c24cd2af5d5ad26f3f69562b945a4e311a63c16e343a5f6;
unchanged21-case manifest SHA
c94cab92b987df76298d87539e0479c9fad9f52acd51d7f5e5d555ffb28aa5f3.
Both jobs used1CPU/4GB/90h, cpu/nslab, with AST, frozen metadata and shell
syntax checks before submission. No existing integrator/fitter was changed.
Batch conditional-mixture integration, fitting, B3 statistical calibration
and independently qualified biological validation remain incomplete.

All three anchored theta4 integrals22042346–48 are now complete without
refusal. Gene3 completed in51:56; its H[v,v] difference from independent
110-digit/order80 Decimal is-1.25923e-10, variance-score difference-2.93760e-13
and logmass difference+2.19160e-14. Root collected the exact three-gene family
at `repair_p3_theta4_anchored_gl_launch_20261005T164218Z/anchored_gl_collection.json`.
All parameter/count/seed/row identities, source/input pins, finite29-vector/
29-matrix dimensions and six reference comparisons per gene were checked.
Full refinement allowances were independently rederived for every gene:
maximum Hessian allowances2.97018e-9,8.30323e-9,3.92765e-9, within1e-5.
The original refused points are preserved. This is successful numerical
refinement at these three points; a full independent Hessian certificate,
accepted observed fit, B3 calibration and biological validation are absent.

The extremely small reported mass-refinement allowances are not certified
floating-arithmetic bounds: observed GL-minus-Decimal logmass differences are
larger, although well within the nominal1e-9 threshold. No tail/refinement
estimate is being substituted for the independent comparison.

## Exact failed-fit point completed, 2026-10-05 14:00 EDT

Job22042612 completed with exit0 in56:44. At the saved call12/gene2 point,
the anchored native BB arithmetic returned the unchanged score-only GL
integral without refusal. Its logmass is-71.12453532557393 and variance score
171.28390376546417. Differences from independent110-digit/order80 Decimal
are+7.55646e-15 and-8.55710e-13, respectively. Decimal40-to80 differences
at110digits are2.00599e-24 in logmass and-1.89973e-19 in variance score.
The original sequential/parallel fit failures remain preserved.

Root checked all55 source/input pins, exact p/v/count/call identities and
original failure bytes, and independently rederived all six GL/Decimal
comparisons, seven Decimal refinements and the complete GL refinement
allowances. Review: REC
`repair_p3_observed_failure_anchor_launch_20261005T165931Z/independent_review.json`.
An optional subsequent filename inspection failed because it used the wrong
Gaussian-rule glob; the completed arithmetic review was already saved and
was separately confirmed. No numerical job failed in that inspection.
Maximum combined score-refinement allowance is2.84217e-13, smaller than the
observed independent variance-score difference. Such allowances therefore
remain empirical refinements, not certified floating-arithmetic bounds.
Full independent29-score/Hessian accuracy, an accepted observed fit,
statistical calibration and qualified human validation remain unfinished.

New isolated `repair_p3_bb_mode_anchor_cache.py` retains one normalized
full support for consecutive tail/central/observed-count requests at exactly
the same n/eta/r/Hessian mode. Exact float keys, invalid-input checks,
read-only retained arrays and independent returned copies prevent approximate
reuse or caller mutation. The scalar arithmetic and native conditional law
are unchanged. Retained arrays are capped at8MiB; this is not a total-process
memory bound. No existing fitter or integrator was edited.

The companion `repair_p3_bb_mode_anchor_cache_check.py` requires all21 saved
independent BB criteria, exact cache-invalidation/mutation checks and bitwise
conditional LP/score/H agreement for173 synthetic rows at theta4 and failed
call12, each at three u values with/without Hessians. Only then does it time
the24-row conditional calculation and evaluate the exact saved score-only GL
point with unchanged orders, precision thresholds and evaluation cap. The
complete GL refinement ladder must agree bitwise with22042612. No whole-fit
speed claim or calibration follows from this diagnostic.

## Claude source review and revised reuse check

Claude CLI source review completed successfully with tools disabled after
the reported quota-reset time. Its source-only response is saved at REC
`claude_cache_review_20261005T180824Z/response.txt` (SHA
7bd1e255693b1074216ac5a28ac7ebe9992a6f4c30f454f41886efdc92cc592c).
Claude confirmed full-support indexing preserves the scalar law, but identified
unenforced thread ownership, potentially lossy non-float64 key inputs and
comparison/test-coverage limits. It did not execute checks or approve a fit.

First reuse diagnostic22044584 passed21 saved independent BB criteria and
2076 conditional LP/score/H comparisons at the two synthetic vectors. Its
`np.array_equal` comparisons establish exact numerical equality, not strict
dtype/signed-zero byte equality; earlier prospective descriptions of this
check as bitwise were too strong. Initial conditional timing ratios are
1.3907,1.3931,1.3926 without Hessians and1.5991,1.6519,1.6489 with Hessians.
These are component timings; the full score-only GL comparison is still live.

Separate `repair_p3_bb_mode_anchor_cache_v2.py` enforces constructing
process/thread ownership and admits only float64 predictors or exactly
represented bounded Python integers. The companion v2 checker compares
dtype, shape and raw bytes; checks39 calls/12 misses/27 hits; exercises
score-only-to-Hessian invalidation, endpoint/U-shaped cases, unsupported input
types, caller mutation and rejection from another thread. It cross-checks
the scalar SHA during metadata inspection and balances scalar/reuse timing
order. It retains the same scalar computation; the executing v1 sources
remain unchanged.

V2 component-only job22044727 and v1 complete-integral job22044584 each
use1CPU/4GB/90h oncpu/nslab, individual non-array submissions. Executable
checks remain on compute nodes. AST, metadata-only inspections and shell
syntax passed before submission. Two earlier static source-location/pin
failures were fixed before any job submission and their launch directories
were preserved. V2 launch:
`repair_p3_bb_mode_anchor_cache_v2_launch_20261005T181527Z`;
checker SHA89cd3f462b41370741bdbe79df5bbeea848828a2a8283994b299c2b725bea479,
helper SHAfcd05d2f3a919bfd6257a402bfd8e869c0941d6d03797f3212ad34e34196802f.
No observed fit replacement, statistical calibration or biological adoption
was submitted. Technical-only changes do not require the Resource-scope
validator; it was not run for this work.

V2 component check22044727 completed exit0 in57seconds. All21 selected
independent BB criteria passed. All2076 native conditional LP/score/H
comparisons passed strict dtype/shape/raw-byte equality. Reuse counts were
exactly39 calls,12 misses and27 hits; unsupported predictor types and use
from another thread were rejected. Measured process peak77744KiB is75.92MiB,
below4GB allocation. Balanced-order conditional timing ratios were
1.4218,1.3923,1.3684 without Hessians and1.6401,1.6312,1.6417 with Hessians.
These remain component measurements, not complete-fit acceleration.

Root checked all92 source/input pins,21 point/reference identities and
criteria, reconstructed the independent Decimal stability ladders and all
six timing ratios, and reviewed the byte-comparison assertions. Saved review:
`repair_p3_bb_mode_anchor_cache_v2_launch_20261005T181527Z/independent_review.json`.
22044584 remains running its full score-only integral comparison, with no
stderr or failure receipt at the last inspection. Required Hessian-family
checks remain13/14 complete and five-gene scalar comparisons7/15 complete;
their collectors declined to write a complete collection. The14 independent
Hessian jobs remain live; actual numerical work is in srun step0, so the
shell-only batch-step RSS is not the numerical process memory. Both original
observed fit failures remain unchanged. B3 calibration and qualified human
validation remain incomplete.

Files added in this bounded reuse work:

- `/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_bb_mode_anchor_cache.py`
- `/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_bb_mode_anchor_cache_check.py`
- `/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_bb_mode_anchor_cache_v2.py`
- `/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_bb_mode_anchor_cache_v2_check.py`

The two Codex-owned exchange documents and shared locked job index were
updated. Commands executed: stdlib source/receipt checks and AST parsing,
metadata-only `--inspect`, `bash -n`, two individual `sbatch` submissions,
`scontrol`/`squeue`/`sstat`, the two existing family collectors and
`git diff --check` (passed). No commit, production dependency or human-authority
change was made by this work.

## Same-process reuse diagnosis and complete reference families, 2026-10-06

Prior goal turn: progress, through collection of the14 independent Hessian
jobs and identification of the cache comparison failure. This turn preserved
that failure and tested its first original GL rule in one process. New script
`repair_p3_cache_rule_diagnosis.py`, SHA
2bad5dfb6fa9fd966c6aa970578403caf7807c481914e675b1324a23b6039e28,
was AST/shell/metadata checked, archived with94 input/source pins, and submitted
as individual22075506 oncpu/nslab,1CPU/4GB/90h. It completed exit0 in1:28.
All24576 scalar/cache conditional pairs at1024 nodes matched dtype/shape/raw
bytes. Full-model versus gene-only predictors and node-level aggregate scores
also matched, and the first-rule logmass and29-score vector matched bytes.
Measured process peak75052KiB. The other eight GL rules were not rerun; this
does not erase22044584's exact cross-job equality failure or prove its cause.
Output: REC`repair_p3_cache_rule_diagnosis_22075506/summary.json`.

The five-gene experiment is now15/15 terminal. Both collectors retain three
candidate integral refusals and the original saved-summary numpy.bool_ bug.
Complete collection:
`repair_p3_five_gene_accuracy_launch_20261005T143525Z/five_gene_accuracy_collection.json`;
independent arithmetic collection:
`repair_p3_five_gene_collection_20261006T231455825881Z/collection.json`.
Original GH20/GH40 satisfy the recorded likelihood-only comparisons in this
five-gene family. This does not establish derivative, full-fit, likelihood-ratio
or calibration accuracy. Candidate likelihood accuracy across all15 is false.
The14 finite-difference reuse shards are also complete and collected at
`repair_p3_hessian_reuse_launch_20261005T142226Z/hessian_reuse_collection.json`;
their accuracy flag remains false.

New stdlib-only `repair_p3_hessian_step_diagnosis.py` retains2352 adjacent
Richardson extrapolations from all2940 saved independent curvature records.
This is explicitly a postdiagnostic synthetic numerical analysis. It retains
every pair, not a best step. Fourth-order central and fourth-order forward
one-sided stencils are both checked, and incompatible stencil pairs would be
retained without extrapolation. None were incompatible in this family. An
initial source check correctly stopped on the forward stencils; the frozen
failed attempt was preserved before the script was corrected to recognize
both native stencils. Output:
`repair_p3_hessian_step_diagnosis_20261006T231708Z/diagnosis.json`.

Extrapolation does not resolve the entire family. For the1e-4/3e-5 pair,
three comparisons per case remain above1e-5. They localize to dispersion
column10: score coordinate24 atfailed0gene2/domain26 and variance coordinate28
atpositive18/domain24, across the three complex steps. Maximum differences
are1.05756e-5 and1.81934e-5. Coarser/finer pairs and all failures remain visible;
no step was chosen to declare success. Reference estimates are not certified
error bounds. The next independent curvature check should address these
dispersion-bound mixed derivatives on their exact original synthetic counts,
rather than reusing the distinct B-FIX count fixture. No observed fit,
calibration run, actual allelic outcome access or scientific adoption occurred.

### Independent Decimal derivative checks, 2026-10-06

New `repair_p3_decimal_mixed.py` derives native beta-binomial eta/r mixed
curvature using Decimal rising products and reciprocal sums, separately from
the candidate float derivatives and the prior complex-step reference. Source
SHA256: `73cd7ecbeb5df22d8398aefcec9cb987eabc5d02456a4f8d20012efe6abc0980`.
It reconstructs the original P3 draw counts, not the B-FIX observed fixture.

Jobs 22075524/22075525 completed in 4:40/4:39. Each retained 12 row/u/precision
records with two probability-only finite-difference steps (1e-16, 1e-20),
covering the deepest target row and deepest opposite-orientation row in the
full draw, at u=0 and +/-sqrt(v), with 80/110 digits. Root independently
recomputed all 192 reported differences from the saved analytic and FD values;
all meet the unchanged absolute 1e-18 criterion. Maximum differences are
3.50080e-29 (failed0gene2) and 2.22383e-31 (positive18). All source/input hashes
remain identical. Process-reported peaks were 84,244/84,224 KiB. Results:
REC `repair_p3_decimal_mixed_check_failed0gene2_22075524/summary.json` and
`repair_p3_decimal_mixed_check_positive18_22075525/summary.json`.

The initial launches 22075515/22075516 failed before derivative computation
because gene2 has only one orientation. Those failures remain saved. The
corrected pilot selects its opposite orientation from the complete original
draw; it does not change any counts or target-gene membership.

Claude independently rederived the source formulas and found no algebra error,
but correctly identified that this pilot does not test global omega_S chaining
or agreement with the native selected likelihood. Its source-only review is
REC `claude_decimal_mixed_review_20261006T232737Z/response.txt`, SHA256
`e89721eb5463aa336b962ad0ff219f4cab8e97dd2f8e9e372d07b638fdbd3a31`.
This review is not an executed numerical check.

New `repair_p3_decimal_global_chain_check.py` checks all 24 target-gene rows
at all three u values and both precisions. It independently recomputes the
complete native parameter map for probability-only FD of p[10]/p[24], compares
the summed derivatives with `conditional_moments`, and ties each row to native
predictors and selected likelihood using the separately checked mode-ratio BB
arithmetic. It records observed/support posterior heterozygote weights for
every row and requires nonzero stage somewhere in the target gene. Native
absolute criteria are LP 1e-9, scores 1e-8, mixed curvature 1e-5; Decimal FD
remains 1e-18 and 80/110-digit stability 1e-50. No tolerance was widened.

AST, bash syntax and both source-only inspections passed. Individual jobs
22075537/22075538 were submitted with cpu/nslab, 1 CPU, 4 GB and 90 hours each;
both were RUNNING at the submission check. Live status must be read from
SLURM. Frozen launch:
REC `repair_p3_decimal_global_chain_launch_20261006T233701Z/`;
checker SHA256 `78de18780af4f62c956cfefd9d390b14fa5f4bb683ac58ce40133693ae415fb6`.
These jobs do not integrate the Gaussian random effect, fit a model, establish
curvature accuracy of the previous candidate, or run calibration. The planned
parallel integral family has not been submitted. Human outcomes remain closed.

Claude's second source-only review also found no algebra/index error, but
requested native Jacobian comparison and per-row FD criteria to prevent
cancellation in the 24-row sum. Review:
REC `claude_decimal_global_chain_review_20261006T233742Z/response.txt`,
SHA256 `7a41b23469fa9d04f79ede6fa557d0b2bcfb9b27b75cf6d3b1e8ac5fc5aeebd2`.
It did not execute code or inspect the separately reviewed probability/anchor
dependencies. Its minor observations about zero coefficients and native global
assembly are addressed in the next checker as well.

The first saved u=0/80-digit record passed in each running job: all 24 native
row comparisons and the global FD criterion. Minimum observed posterior
heterozygote weights were approximately 0.864 and 0.983, so the rows do exercise
non-negligible homozygote contributions. These partial records are not full
80/110-digit family acceptance.

New `repair_p3_decimal_global_chain_v2_check.py` retains each row's FD
differences from the same probability evaluations, requires every row and the
sum to meet 1e-18, verifies native Jacobian entries for coordinates 10/24, and
compares the full native `Model.evaluate` result for LL, both global scores and
their mixed curvature. It requires the saved stage/dispersion coefficients
16,24,25,26,27 to be nonzero; both original saved points satisfy that condition.
No parameter or count was changed. AST, shell syntax and both source-only
protocol inspections passed before submitting individual jobs 22075551/52,
each cpu/nslab, 1 CPU/4 GB/90h. Frozen launch:
REC `repair_p3_decimal_global_chain_v2_launch_20261006T234117Z/`;
checker SHA256 `8f80b5303d7e119e83e4fd486a25e9a2523aa0e75ab09e469557bfc3c110cdc0`.
The earlier live checker and its saved results are preserved. Both versions
remain diagnostic synthetic checks; neither establishes Gaussian integral
accuracy, a constrained optimum, calibration, or human validation.

### Decimal Gaussian mixed-moment family, declared before integration

The preceding all-row checker jobs 22075537/38 completed in 8:16/8:17 and
passed their six records per case. The strengthened per-row/native-Jacobian
jobs 22075551/52 remain the required prerequisite for the integral family.

New `repair_p3_decimal_gl_moments.py` constructs Decimal GL16/GL32 roots with
Newton refinement of float seeds, checks every polynomial moment through
degree 2*order-1, and integrates six raw moments. For z=u/sqrt(v), the moments
are 1, ell_r, ell_omega, gv, ell_r,omega+ell_r*ell_omega, ell_r*gv,
where gv=(z^2-1)/(2v) is the normal-prior derivative at fixed u. The integrand
weight is exp(ell-z^2/2); the normal constant cancels in posterior means.
No parameter-dependent rescaling is used. The collector forms both mixed
curvatures only after summing all chunks over the complete declared domain.

The entire family is fixed before evaluating these native integrals:

| Rule | z bound | Panel spacing | GL order | Panels | Chunks per case |
|---|---:|---:|---:|---:|---:|
| 0 | 24 | 0.25 | 16 | 192 | 16 |
| 1 | 26 | 0.125 | 16 | 416 | 16 |
| 2 | 26 | 0.125 | 32 | 416 | 16 |

Both original P3 cases, all three rules, and both 80/110-digit precisions are
required. Each of 96 individual jobs covers one case/rule/chunk and computes
both precisions. This is not an array and does not generate new Monte Carlo
draws. The prospective allocation is 1 CPU, 4 GB, 90h per job on cpu/nslab,
384 GB total if all 96 run together. The live nslab limit was checked as
850 CPUs and 7 TB per user; other users' or sessions' jobs are preserved.

The whole-line tail calculation uses deliberately loose analytic bounds on
conditional scores and mixed curvature, plus Gaussian Mills bounds. It is
reported separately from quadrature error. Decimal evaluation is not interval
arithmetic, and domain/order/precision convergence does not certify an error
bound. No complete-Hessian or fit acceptance follows from this selected family.

New `repair_p3_decimal_gl_chunk.py` requires both completed strengthened
admission results with their original counts/parameters pinned, checks all
source/input hashes, and preserves raw moments and node/panel coverage. New
`repair_p3_decimal_gl_collect.py` requires all 96 terminal successful jobs,
complete panel coverage and identical inputs. It retains every rule and
precision, comparing both H[10,24] and H[10,28] with the prior candidate in
both cases. Declared diagnostics: 80/110 precision differences <=1e-30,
successive refinement differences <=1e-8, and the unchanged candidate mixed
entry nominal tolerance 1e-5. It does not select a best refinement or average
chunk covariances; refusals and incomplete families remain visible.

Compute job 22075560 checked polynomial identities and a closed Gaussian
example independently of the native likelihood. At v=0.5 and ell=-u^2/2,
scores ell_r=ell_omega=u and ell_r,omega=0 imply integrated mixed curvature
1/3 and r/variance curvature zero. All six rule/precision records passed
against 1e-40; maximum raw-moment difference was 7.64e-77. Pi for this check
was independently computed using Machin arctangent series. Job completed in
5s; process peak 14,136 KiB. Result:
REC `repair_p3_decimal_GL_rule_22075560/summary.json`.
This proves the checked quadrature identities, not native integral accuracy.

The strengthened native checks 22075551/52 subsequently completed in 7:22/
7:23 with all six records per case passing. Root checked all 2,304 saved
per-row FD differences and independently rederived all 96 summed differences.
Maximum per-row differences were 2.24597e-28 and 1.43035e-31, against 1e-18;
all 80/110-digit stability checks passed. Reported process peaks were
89,840/89,988 KiB. These checks include native Jacobians and the native full
LL/global-score/mixed-curvature assembly, on unchanged synthetic points.

Claude rederived the GL moments, roots, partition and analytic tail algebra
without finding an error. Review:
REC `claude_decimal_GL_moments_review_20261006T234504Z/response.txt`, SHA256
`4c4737b6edf15a4a2315c42120d0711d9eb434e6f9d83b0115bacc192f275d27`.
This was source review only. It correctly noted that the first Gaussian mock
has zero odd moments, so it cannot test their ordering or the product
subtractions. The revised checker uses ell=-(u-0.75)^2/2, scores u and u^2,
and nonzero conditional cross-curvature 2/7. All six moments are nonzero and
their exact expectations follow from the shifted Gaussian posterior. Explicit
`require` checks replace Python assertions. Compute job 22075613 passed all
six rule/precision records against 1e-40 in 5s; maximum raw-moment difference
6.33e-77. Result: REC `repair_p3_decimal_GL_shifted_rule_22075613/summary.json`.
The first symmetric mock and its result remain preserved.

The review also notes that absolute raw-moment tail bounds do not certify
covariance accuracy without a lower mass bound and a numerical integration
error bound. This limitation remains explicit. Rule0/rule1 jointly vary domain
and spacing, so their difference does not attribute error to either factor.
Rule1/rule2 separately vary order. If the declared family does not converge,
its complete results must remain visible before declaring further refinement;
no favorable rule may be selected to clear a discrepancy.

All 96 shell syntax checks and all 96 source-only protocol inspections passed.
The individual integral jobs were then submitted under the displayed resource
header. Every job ID and output directory is recorded in
REC `repair_p3_decimal_GL_parallel_launch_20261006T235003Z/jobs.json`.
All 96 were confirmed RUNNING at the first live scheduler check, with zero
pending jobs and 384 GB total requested memory. Initial collection:
REC `repair_p3_decimal_GL_parallel_collection_initial_20261006T2351/collection.json`,
which correctly reports 0/96 complete and zero refusals. This is a dated
observation, not live status. The shared input manifest has 60 checked pins;
the runner SHA256 is
`7b716bd4a425c87a77ffbab47b10e544981656fde314908b9a2bec1345d24f4d`.
No new fit, calibration, actual/sealed outcome access or scientific adoption.

### Parallel arithmetic alternatives, 2026-10-06 evening (NY)

While the 96 Decimal moment jobs run, a separate numerical route checks reuse
on all nine saved rules of the original failed observed-point integral. New
`repair_p3_cache_rule_family.py` launches eight independent subprocesses of
`repair_p3_cache_rule_family_shard.py`, longest rules first, in one allocation.
Each subprocess reconstructs the original 60-synthetic-participant, 173-row
B-FIX fixture and compares scalar versus v2 support reuse in the same process.
It checks every selected row probability/score, predictor, node aggregate and
rule sum by dtype, shape and raw bytes. The declared family totals 51,072 nodes
and 1,225,728 conditional pairs over gene2's 24 rows. It uses the original
saved rule geometry, rather than re-executing the original GL callable, and
does not claim independent accuracy or a fitted model. All cross-job deltas
against the previous saved ladder are retained separately.

Job 22075903 was submitted as one cpu/nslab job, 8 CPUs, 8 GB, 90h. Both AST
checks, nine source-only inspections and shell syntax passed first. Frozen
launch: REC `repair_p3_cache_rule_family_launch_20261006T235622Z/`;
runner SHA256 `736e435c637495bcdafc017eab70528031e417d8fd1e1b5bb7545eca035456a7`.
Its final submission-helper print raised a duplicate-key TypeError after
the job receipt and ledger were written; SLURM confirmed the existing job
RUNNING, and it was not resubmitted. The earlier cross-job failure 22044584
remains preserved, even if the new same-process family passes.

A second arithmetic hypothesis replaces the O(n) full-support mode recurrence
with requested-count gamma differences. New helpers
`repair_p3_bb_gamma_difference.py`, `_v2.py` and `_v3.py` preserve the native
beta-binomial rho, eta/r derivatives, integer selection and mixture law.
For shapes >=32, stable `log1p`/`expm1` Stirling differences replace differences
of large special-function values; smaller shapes separate j=0 to support
underflowing shapes without discarding their log mass. The expansion and
polygamma identities follow [NIST DLMF 5.11](https://dlmf.nist.gov/5.11) and
[5.15](https://dlmf.nist.gov/5.15). This is a numerical representation of the
same count model, not a new biological model or a different assay method.

| Synthetic component job | Result | Balanced component timing |
|---|---|---|
| 22076024, eight-term v1 | Passed 21 independent Decimal cases, eight direct rising cases and 2,076 native conditional comparisons | Scores ~17% slower than mode ratios; Hessians ~2.06x faster |
| 22076036, v2 skips unused second derivatives | Same checks passed; all 2,076 native results byte-identical to v1 | Scores 1.139–1.140x faster; Hessians 2.115–2.119x faster |
| 22076101, four-term v3 and expanded checks | All component checks passed, including explicit helper call counts | Scores 1.727–1.819x faster; Hessians 4.211–4.234x faster |

Each job used cpu/nslab, 1 CPU, 4 GB, 90h, and completed in 48s, 47s and 30s,
respectively. The timings cover the fixed gene2 rows at 32 u values in three
alternating-order repeats, separately for scores and Hessians. They are not
whole-integral or whole-fit speed measurements. Process-reported peaks were
76,064, 76,136 and 76,116 KiB. Results are under REC
`repair_p3_bb_gamma_difference_check_22076024/`,
`repair_p3_bb_gamma_difference_v2_check_22076036/` and
`repair_p3_bb_gamma_difference_v3_check_22076101/`.

Claude independently rederived the eight-term formulas, small-x series and
chain without finding an algebra error. Review:
REC `claude_gamma_difference_review_20261007T000229Z/response.txt`, SHA256
`9cb43bd34059b2ccc6610a7971a16c03c6ccd396d5dd01c100b8fe39976a3b18`.
It was source review, not execution. It identified weak small-x coverage and
an unverified dynamic-helper lookup. V3 records the native helper call counts
for all twelve 173-row comparisons, verifies the original functions' globals
and BB lookup, and adds counts 999/1000/1001 around x=0.001. Each helper was
actually called in every comparison. Relative derivative error <=1e-10 is
additionally required where the expected value exceeds 1e-250; extreme
underflow uses the declared absolute limiting-value checks. Shape 31.999 and
32.001 bracket the expansion boundary; exp(log(32)) itself need not land on
the same branch on every platform.

Four terms were selected from the analytic positive-shape remainder before
the v3 numerical test, rather than from biological outcomes. At shape32,
the first omitted-term bounds for the log-Gamma difference and its first/
second log-shape derivatives are approximately 4.785e-17, 4.307e-16 and
4.737e-15. These bound the Stirling truncation, not floating roundoff or the
complete likelihood/integral error. Eight-term versions and all timing results
remain preserved. Root checked all v3 saved criterion flags, references,
nonzero helper counts and unchanged input/source hashes. The largest v3
independent selected logmass difference was 1.32154e-11; all unchanged
LP/score/Hessian criteria passed. Native conditional differences versus mode
ratios reached 1.32143e-11, 8.06466e-13 and 2.44427e-12, respectively.

The next check is the complete original score integral at the exact previously
failed observed-fit call12, using the distinct B-FIX count hash 7009be4a...,
not either P3 resampled count vector. New
`repair_p3_gamma_observed_point_gl_check.py` requires the completed v3 component
checks and identical admitted helper source. Job 22076126 uses 1 CPU, 4 GB,
90h on cpu/nslab. Frozen launch:
REC `repair_p3_gamma_observed_point_GL_launch_20261007T001154Z/`.
It runs the unchanged GL callable/settings, retains all nine rule comparisons,
compares all 29 scores against the previous saved integral at 1e-8 and LL at
1e-9, and retains all six Decimal references. The final independent reference
is prespecified as 110 digits/order80. The known order20 variance-score
disagreement remains visible; it is not converted into a passing order20
claim. Helper source SHA256:
`546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3`.
AST, source-only inspection and shell syntax passed before submission.
This is one gene's score integral, not a whole observed fit or a full
independent Hessian check. All fit, calibration, identity, C5/lambda and sealed
replication requirements remain. No actual/sealed allelic outcomes were read.

At the 2026-10-07 00:15 UTC collection, 30/96 Decimal integral chunks were
terminal successful and passed the collector's actual result checks: source
pins, original count/parameter identity, both precisions, node counts and
polynomial criteria. Sixty-six were still running; no refusal was present.
Collection: REC `repair_p3_decimal_GL_parallel_collection_20261007T001516Z/`.
Failed0gene2 had16/48 chunks done and21,904/46,080 panel-complete nodes recorded;
positive18 had14/48 and22,320/46,080. These are numerical workload fractions,
not progress percentages for model validation. The complete 96-chunk family
is still required before forming the declared reference comparisons. The reuse
family and the saved-point gamma GL job remained live. Current long-check
allocations total276GiB (66x4GiB +8GiB +4GiB). Live status remains SLURM's.
No new C5/lambda, identity or calibration-setting answer was present in the
existing source-owner request and design-answer files inspected this turn.

### Completed arithmetic and boundary checks, 2026-10-06 late evening (NY)

Whole nine-rule support-reuse job 22075903 completed in 24:07. All 51,072
nodes and 1,225,728 conditional pairs matched byte for byte between cached and
scalar calculations in the same processes. All nine new scalar rule outputs
also exactly matched the earlier cached original-GL outputs. The earlier
scalar job still differs by at most 1.98952e-13 in scores, with identical
likelihoods. The cache is not necessary to produce the cross-job difference;
its cause remains [Unverified]. Earlier job 22044584 remains failed. Root
review: REC `repair_p3_cache_rule_family_review_20261007T002143Z/review.json`.

Gamma v3 whole saved-point job 22076126 completed in 13:10. This is the
original B-FIX observed synthetic count vector (7009be4a...), call12, gene2,
v=0.0009636462009520182, not a P3 resampled count vector. All nine original GL
rules and all 29 scores passed the declared comparisons. All six Decimal
references remain reported, including the order20 variance-score disagreement.
Against the prespecified 110-digit/order80 reference, likelihood and variance
score differences were -1.51981e-11 and -7.04605e-11. Root independently
recomputed all saved differences and checked source pins and terminal status:
REC `repair_p3_gamma_observed_point_GL_review_20261007T003022Z/review.json`.
Pure GL time was 781.30 seconds. [Estimate] Earlier stdout markers place the
corresponding scalar integration near 2,278 seconds, about 2.92 times longer.
This is cross-job timing, not a controlled benchmark or whole-fit speed.

Prepared `repair_p3_gamma_parallel_model.py` keeps the original parent BB,
zero-variance path, fit starts, bounds, ties, and restarts; only positive-
variance fork workers install Gamma v3. Claude's source review found a missing
accumulated Hessian criterion, snapshots that could alias arrays, incomplete
failure metadata, and overly broad numerical-failure classification. These
were corrected before execution. Current source SHA256:
`fb869f667efd2fecbb5be3b378aa7a88c617e193ce1aea70cbe7f0c084f61013`.
The positive eight-worker path has not run. Workspace estimates are explicitly
not measured peak memory or general biological-scale claims.

Boundary job 22077790 completed in 44 seconds. It retained all nine prespecified
vectors: saved native null, native cold start, and saved call12, each with kappa
shifts 0/-0.1/+0.1. The new parent's v=0 likelihood, 29 scores, and native
Hessian results were byte-identical to original calls. An isolated Gamma/native
comparison passed LL absolute 1e-8, score absolute 1e-8 plus relative 1e-6,
nuisance-Hessian absolute 1e-5, and unchanged restart decisions at all points.
Maximum absolute differences were LL 2.42153e-11, score 7.90335e-8, and nuisance
Hessian 1.62984e-9; the score maximum passes the combined absolute/relative
criterion, not an absolute-only 1e-8 claim. Native undefined variance-Hessian
entries remain undefined. Root recomputed every comparison from saved arrays
and checked all source/input pins:
REC `repair_p3_gamma_boundary_review_20261007T004500Z/review.json`.
This is not an independent full Hessian, positive-variance worker test,
observed model fit, calibration, or human replication.

At the 2026-10-07 00:42 UTC collection, 62/96 Decimal reference chunks were
terminal successful, 34 were running, and none had refused. Actual completed
outputs were checked; no curvature was formed before the whole family
completed. Collection:
REC `repair_p3_decimal_GL_parallel_collection_20261007T004200Z/collection.json`.
All 96 chunks and prespecified refinements remain required. No actual or sealed
allelic outcomes were read, and no model or paper numbers were adopted.

Claude's follow-up source review completed under REC
`claude_gamma_parallel_model_followup_20261007T004311Z/`; response SHA256
`d661fd7e5be630151259e895b34a59c86b2e174bc3137b260847e5eb8b1236d3`.
Root confirmed `kernel.require` raises a ValueError subclass, which the original
resampling functions catch as a fit failure. Allocation, scope and callable
checks now raise `ModelIntegrityError(RuntimeError)` instead. Integration state
is reset before preflight and at v=0, and worker refusal metadata now includes
its phase/type/message. Original fitter methods do not compare node counts;
positive GL keeps `nodes` as a compatibility argument, as the earlier integral
implementation does. Plain worker errors without a recorded integration
ladder conservatively stop; their numerical classification has not been broadened.
Global gates remain strict, including possible accumulated numerical refusal.

The corrected model SHA256 is
`06b384d6a6b6185dc5c0eac747f2529a8491c8ac07a7a14167ea876b05cf8fc2`.
The preceding boundary pass applies to its earlier pinned source, not silently
to this revision. New boundary job 22077922 repeats the full nine-point family
and additionally attempts a positive call in its declared one-CPU allocation:
it must stop as RuntimeError before any worker starts, with cleared state.
Launch: REC `repair_p3_gamma_boundary_v2_launch_20261007T004620Z/`.
AST, source/input inspection, and shell syntax passed before submission.

Job 22077922 then completed in 42 seconds with exit0. All nine boundary
comparisons passed again, with byte-identical original-parent results. The
one-CPU positive attempt stopped as RuntimeError before any worker and left
no stale integration state. Root rederived all nine numerical comparisons,
checked the preflight record and rechecked source/input pins:
REC `repair_p3_gamma_boundary_v2_review_20261007T004900Z/review.json`.
At 00:48 UTC the result-reading collector accepted64/96 Decimal chunks,
32 running, zero refusals; both rule0 families were complete. Rule1 and rule2
still require their complete families. Collection:
REC `repair_p3_decimal_GL_parallel_collection_20261007T004800Z/collection.json`.
The next numerical comparison is the complete positive-variance G8 evaluator,
sequential versus eight workers. It has not been submitted. Full independent
curvature, fit and calibration requirements remain unresolved.

### Complete positive-variance evaluator comparison, 2026-10-06 evening (NY)

The next comparison is now submitted as two independent jobs, rather than
serially waiting for both implementations in one allocation. Job22077927 uses
8CPUs/8GiB for the eight-worker evaluator; job22077928 uses1CPU/4GiB for a
separate direct scalar-GL loop over genes. Both use cpu/nslab and90h, with all
four numerical thread settings fixed to1 before imports. Live scheduler limits
were rechecked: cpu maximum14d; nslab per-user850CPUs/7TB. These are individual
jobs, not an array. Launch:
REC `repair_p3_gamma_full_point_launch_20261007T004948Z/`.

Both use the exact saved failed observed-fit call12 parameters and
v=0.0009636462009520182, all173 synthetic B-FIX rows,8genes,60 synthetic
participants, and the original count hash7009be4a.... This deliberately tests
the point the previous implementation could not evaluate, rather than a newly
chosen easier vector. All29 scores and the29-by-29 Hessian are requested.
Each gene's GL tolerances are LL1e-9, score1e-8, Hessian1e-5/8; the full sums
must meet LL1e-8, score1e-7 and Hessian1e-5. The direct sequential loop does
not call the parallel worker function. Both retain every gene result/refusal,
and non-numerical integrity failures stop the job. No optimizer is called.

New checker `repair_p3_gamma_full_point.py`, SHA256
`928afdf60d5154b7f358c73a6628aceab9a7b87e18d374b4c4a58c2fc8edca19`,
requires the corrected nine-point boundary admission, identical sources,
exact fixture and numerical environment. Source AST, both metadata inspections
and both shell syntax checks passed before submission. Source/model files are
pinned while these jobs run. Both started and wrote their checked protocols;
this does not imply a successful positive-variance evaluation.

New result collector `repair_p3_gamma_full_collect.py`, SHA256
`71f435216937443c85f8b1af0951d5de3739c65cf3c847f14d6fdbce7ba7e2b6`,
requires both jobs terminal COMPLETED/exit0, checks the actual executed sources,
protocols, eight gene results, full dimensions, and original input pins. It
rederives all totals from ascending gene contributions, checks accumulated
allowances, compares per-gene and total LL/scores/Hessians at1e-9/1e-8/1e-5,
and checks structurally zero cross-gene alpha/r derivatives. Exact report
agreement is recorded separately; cross-job floating arithmetic need not be
byte-identical. All numerical refusals remain visible, without a success claim.
The waiting branch ran successfully with both live jobs:
REC `repair_p3_gamma_full_collection_20261007T005300Z/collection.json`.
The completed-pair branch has not run. This comparison can establish agreement
of two implementations; it cannot establish independent Hessian accuracy,
model convergence, calibration, human replication or paper adoption.

### Full-point source review and refusal handling, 2026-10-06 evening (NY)

Claude reviewed the submitted runner and first collector without tools or
numerical execution. Review: REC
`claude_gamma_full_point_review_20261007T005250Z/response.txt`, SHA256
`2cb19f1ed55466f0289b33fc4f77641cb16ee4afe3a6927ba0922908301cda5e`.
Root confirmed that scalar GL creates a separate kernel module on every call.
Its `IntegrationFailure` class therefore differs from the runner's separately
loaded class. The frozen sequential runner can misclassify a recorded numerical
refusal as an integrity failure. A resulting job failure is not evidence of a
different likelihood or a failed scientific refit. Preserve the raw records;
do not call it successful or automatically restart a live job.

The working runner now recognizes numerical refusals through their ValueError
base plus explicit integration ladder and evaluation-count attributes, rather
than class identity. Failed-domain nonfinite diagnostics are preserved with
explicit `nonfinite_float` tags; successful results still require finite JSON.
Future runner SHA256:
`da3e220f0dbe82b336ddd03b648c9467904143d7005a8ed9988e2cbb0da49390`.
A lightweight check used two distinct ValueError subclasses with the recorded
contract; both were recognized, while an unrecorded ValueError and RuntimeError
were rejected. AST passed. This is not execution of the revised numerical path.
The running jobs execute the frozen launch copy, not the working runner. Root
rechecked their executable hash and every live input pin after the edit; all
remain unchanged. The earlier source and results remain preserved.

The collector was also strengthened before a completed-pair comparison. It
verifies the executable Settings expressions through AST and derives their
values from the pinned kernel, rather than treating a prose settings string as
evidence. Both paths use identical LL1e-9/score1e-8/H1.25e-6 per-gene criteria
and all other original Settings defaults. It now ties output directories to
job IDs, requires every executed source, checks failed worker dispositions,
treats cross-gene derivative defects as integrity errors, and reads global
allowance limits from the protocol. Original report files contain no explicit
tolerance fields, so the source check is necessary. Exact ascending gene-sum
checks remain: source inspection confirms both paths use that same float64
order. Process RSS remains labeled as process RSS, excluding worker totals.
New collector SHA256:
`253c15c0b2949154cc04968fb1dfc0d18f4f46bb6894cbfb9bfb65f59efee4aa`.
The original launch collector is preserved; the new frozen copy is
`repair_p3_gamma_full_collect_v2.py` in the same launch directory.

The new waiting branch and source-settings check passed at 00:58 UTC with
both jobs still RUNNING. Actual completed-pair comparison remains not run:
REC `repair_p3_gamma_full_collection_v2_20261007T005814Z/collection.json`.
The Decimal collector again checked64/96 completed chunks,32 running and
zero refusals at the same time:
REC `repair_p3_decimal_GL_parallel_collection_20261007T005814Z/collection.json`.
Current scheduler CPU use confirmed eight active CPUs for the forked evaluator
and one for the sequential calculation. No jobs were cancelled or duplicated.

### Gamma conditional batches and gene-local moments, 2026-10-06 evening (NY)

A separate acceleration experiment now combines the unchanged existing
gene-local GL/moment implementation with scalar native conditional calls using
the checked Gamma v3 BB. Its batch callback loops the original conditional at
each eta; it does not use the earlier floating prefix BB arithmetic. This
can reduce repeated predictor and dense-matrix work, but no whole-integral
speed improvement is yet established. The two full G8 evaluator jobs remain
running with unchanged sources and inputs.

New checker `repair_p3_gamma_local_check.py` first requires byte identity at
every original conditional output for all173 B-FIX rows, theta4 and saved
call12, u=0/+sqrt(v)/-sqrt(v), both score-only and Hessian calls. It then attempts
the saved call12 gene2 score integral and Hessian integral separately, with
unchanged LL1e-9/score1e-8 and tightened per-gene Hessian1e-5/8, batch32,
original GL orders, domains and evaluation cap. It retains numerical refusals
and saves accepted integrals before reference comparisons. Accepted numerical
arrays must be finite and have the original29-coordinate dimensions; rejected
ladder nonfinite diagnostics receive explicit tags. The original Gamma whole
score reference is checked for exact point/count identity and length29.

Claude's read-only review found no blocking numerical defect for bounded
execution. Review:
REC `claude_gamma_local_review_20261007T010125Z/response.txt`, SHA256
`80d6e3b02f62e5d0bbb9419bec7e9088e02f77f10208e36d6416a31015c1a8aa`.
Root addressed its result-preservation findings before submission and confirmed
the predictor dispersion formulas, three-field environment, Settings dataclass,
kernel SHA and reference schema from the actual sources/JSON. Source AST,
input/reference inspection and launcher shell syntax passed. Submitted checker
SHA256: `edbcb159b5aaf396e9309c01632d5767e78f659cce4326882c5900c73ba15f66`.
Existing local GL/moments SHA256:
`254a1efce7d8a873f8caeda1f0630cb48e35286e15bd952352246ab26b01e4de` /
`95ec9f270d18bfbbc8f31f0857308f5f2acad053661a30b208049b7a3f59f31f`.

Job22077951 uses1CPU/4GiB/90h on cpu/nslab, individual job name1597.
Launch: REC `repair_p3_gamma_local_launch_20261007T010143Z/`.
It passed all2076 component byte comparisons and entered whole-integral
evaluation. Root independently reconciled all four173-row families against
the actual synthetic rows, checked every point flag and nonzero Gamma helper
call count, and rechecked executable/input hashes:
REC `repair_p3_gamma_local_component_review_20261007T010625Z/review.json`.
This proves batching/stacking identity at those points, not new independent BB
accuracy, whole-integral accuracy, fit agreement, calibration or human replication.

At01:06UTC the Decimal result collector accepted74/96 completed chunks,
22 running and0 refusals:
REC `repair_p3_decimal_GL_parallel_collection_20261007T010627Z/collection.json`.
It still forms no curvature before the full family is complete. No actual or
sealed outcomes were read, and no scientific release numbers were changed.

### Batched Gamma BB arithmetic, 2026-10-06 evening (NY)

New `repair_p3_bb_gamma_batched.py` evaluates the same Gamma v3 rising-product
formulas across a bounded node-by-count matrix. Every shape chooses its own
>=32 Stirling or j=0-separated small-shape branch; count0/1 overrides remain.
The native rho transformation, four Bernoulli terms, score/Hessian chain and
scalar total-shape term are unchanged. It is limited to64 nodes,100000
node-count entries and n<=5000, with explicit refusal before matrix allocation.
It does not implement the conditional genotype mixture or an integral yet.
Helper SHA256:
`5ae40297221af8158de9dad04d32294e163ef2740f344b055d2acad60e9c0b4a`.

Checker `repair_p3_bb_gamma_batched_check.py`, SHA256
`3f184a2fbff25d1e4adf5e6b501851363b063c706d8d0277d8c26274b4128d7e`,
requires the completed Gamma v3 checks and the original independent BB
receipts. Job22078120 (individual name1500,1CPU/4GiB/90h,cpu/nslab) completed
in15 seconds with exit0. Launch:
REC `repair_p3_bb_gamma_batched_launch_20261007T011135Z/`.
It passed64 independently rederived80-digit rising-product triples across
eight mixed log-shapes (including -1200 underflow and both sides of32),
21 independent110-digit selected BB points, and30 scalar comparison settings:
n0/5/100/3000/5000 by r-15/-3/5 by score/Hessian, each with64 eta nodes from
-30 to30. Original absolute and relative criteria were unchanged. Mixed-row
permutation and score-only branch checks passed; oversized64-by5001 support
was refused. Process-reported peak RSS was52884KiB, not total cluster memory.

Three alternating-order component timings each for scores and Hessians used
32 eta nodes, n3000, r-15, and300 tail counts. Timed scalar and batched outputs
were exactly equal. Scalar/batch score times were18.26–18.50ms /3.48–3.74ms
(4.95–5.25 times); Hessian times27.77–28.87ms /6.39–6.73ms (4.13–4.50 times).
These are short component measurements in one job, not conditional-mixture,
whole-integral or whole-fit speed, and not a runtime guarantee. All raw repeats
remain saved. Root rederived every saved nominal criterion, reconciled all
30 settings and64-node evaluations, recomputed timing ratios, and checked
source/reference hashes:
REC `repair_p3_gamma_batched_review_20261007T011414Z/review.json`.

Claude independently rederived the formulas and checked the mixed-row
indexing, finding no blocking defect. This was source review, not execution:
REC `claude_gamma_batched_review_20261007T011258Z/response.txt`, SHA256
`a5bb6cd8d3fa43b6581ca1455a879ec5ccdd5ea4ad60e2b55d8c7bae941b33a7`.
Its future requirements remain relevant: full support must split node batches
above the entry cap; scalar-source hashing should move outside repeated calls;
native conditional-mixture assembly and complete integrals need separate
checks. Scalar/batch agreement cannot detect shared cancellation error; the
independent rising and selected BB references are therefore retained.

At01:13UTC the Decimal collector checked84/96 completed chunks,12 running
and0 refusals. Full family curvature remained withheld:
REC `repair_p3_decimal_GL_parallel_collection_20261007T011308Z/collection.json`.
Both full G8 jobs and the scalar-Gamma local-integral job remain running.
No optimizer, calibration, actual/sealed outcome or paper adoption is implied.

### Native conditional-mixture batching, 2026-10-06 evening (NY)

New `repair_p3_gamma_batched_conditional.py` inherits the existing conditional
genotype-mixture and integer-selection assembly. It replaces node BB evaluations
with the checked batched Gamma arithmetic. The inherited node splitting keeps
each full count-axis batch at<=100000 entries and returns every requested node.
The older prefix arrays made during preparation are unused at nodes. Raw genotype
priors, orientation, homozygote components and selection/complement branches
are retained. Detached row/parameter/derived-constant snapshots reject stale
contexts; only owned arrays are made read-only. Caller row/signature arrays are
explicitly excluded from freezing. No model or optimizer method is modified.
Context SHA256:
`007e46bc796810e1adaf2c7cec0c273706e9e47aa9770acbe636b334600923a1`.

Claude found no blocking dispatch, orientation or slicing error in source:
REC `claude_gamma_batched_conditional_review_20261007T012526Z/response.txt`,
SHA256 `7fc99913ee7967373f412fb7f0f155dd38433e8fcb4d273222f236ac81feafc2`.
Root verified the prepared `.p` value comparison, original Model.G and float
clone arrays, and matching entry limits from the actual sources. Before
submission, state probes were strengthened to require successful calls before
mutation and after every restore, match the expected error messages, and change
the cached homozygote inclusion constant rather than re-test the inherited r
guard. Only a cloned synthetic row and parameter copies are mutated. Input
array flags are recorded and checked, including an actual array in the clone.

Job22078218 (individual1500,1CPU/4GiB/90h,cpu/nslab) runs new
`repair_p3_gamma_batched_conditional_check.py`, SHA256
`fdb377ddf2ec09cd0901e9550105f0be2a7ccc0555c76ed676c6f8962ea54564`.
Launch: REC `repair_p3_gamma_batched_conditional_launch_20261007T012522Z/`.
Source ASTs, all prerequisite/dependency hashes, input inspection and shell
syntax passed before submission. The primitive Gamma admission is required;
neither thresholds nor saved call12 parameters/counts are changed.

All24220 conditional comparisons passed:2076 at both fixed vectors and three
u values, plus22144 at64 wide eta values across all173 synthetic rows, both
score-only and Hessian modes. Maximum LP/score/Hessian differences versus
native scalar Gamma were0 /9.76996e-15 /5.68434e-14. Both selection branches
ran; observed/tail/selected BB support call counts were360 each in the wide
family. Seven n3000 rows required splitting at64 nodes. Every chunk returned
all requested nodes. Parameter, cached-constant and row-array mutations were
rejected, with accepted calls before and after restores. Root reconciled all
six173-row families against the actual fixture, rederived every saved nominal
criterion, checked node counts and state probes, and rechecked executable and
input hashes:
REC `repair_p3_gamma_batched_conditional_review_20261007T013327Z/review.json`.

The job has entered the complete saved gene2 score/Hessian integrals, using the
unchanged existing local GL, LL1e-9/score1e-8, H1e-5/8, original rule/domain
family and evaluation cap. Whole-integral accuracy/speed remain pending. A
completed diagnostic job will not count as a passing integral unless its
actual score-comparison and Hessian-return criteria pass. Independent full
Hessian, optimizer convergence, statistical calibration and human replication
remain separate requirements. Scalar-source hashing still occurs at each BB
batch and remains a possible later performance improvement.

At01:31UTC,90/96 Decimal chunks passed actual-result collection,6 running and
0 refusals: REC
`repair_p3_decimal_GL_parallel_collection_20261007T013104Z/collection.json`.
The slow workers continue completing panels. Their runtime differences are
not diagnosed; the Decimal implementation evaluates full selected supports.
No existing job was restarted, cancelled or relabeled as a scientific success.

### Batched whole gene2 integrals completed, 2026-10-06 evening (NY)

Job22078218 completed7:15/exit0. Score and Hessian integrals both returned at
the original failed call12 point. Each retained all9 rules and51072 nodes.
Recorded score/Hessian times were142.68977s/217.83606s. Both likelihoods
were-71.12453532558924; both29-score vectors agreed within1e-8. Returned
arrays were finite and inactive-gene alpha/r derivatives exactly zero. Local
accumulated allowances met LL1e-9/score1e-8/H1.25e-6.

All18 rule-wise LL/score comparisons passed against previous scalar Gamma GL.
Final score-only differences were LL-9.94760e-14 and maximum score6.19593e-12.
Root embedded each local score vector using its saved global indices before
comparison. Report/error arrays remain local dimension15; any future G8 sum
must embed those allowances in29 coordinates rather than adding15-wide arrays.

All six saved Decimal references (20/40/80 orders at80/110 digits) remain
reported. Root rederived every Hvv from raw variance second moment minus score
squared. Against the already prespecified110-digit/order80 reference, new
LL/variance-score/Hvv differences were-1.52975e-11/-7.66564e-11/-1.86838e-8,
within nominal1e-9/1e-8/1e-5. New Hvv=-238940.93632921067 is in variance-v
units. Order20 still fails nominal variance-score/Hvv checks. Original
`independent_Hvv_accuracy_accepted=False` flags remain false: selected nominal
agreement is not a certified quadrature bound or full independent Hessian
validation. This is B-FIX7009.../v0.0009636462009520182, not a P3 resample.
Root review: REC `repair_p3_gamma_batched_whole_review_20261007T014502Z/review.json`.

The earlier scalar-Gamma local score returned in1672.61350s. Root checked
identical counts/parameters/settings, all9 rule geometries and51072 nodes.
Recorded scalar-local/batch score ratio=11.72; earlier global-scalar/batch
ratio=5.48. These are different jobs/nodes, not controlled whole-fit speed.
The scalar-local Hessian and both full G8 jobs remain running. At01:39UTC the
Decimal collector still accepted90/96 chunks,6 running,0 refusals:
REC `repair_p3_decimal_GL_parallel_collection_20261007T013958Z/collection.json`.
Next is full G8 evaluation with batched conditional arithmetic and correctly
embedded error arrays. Full independent curvature, fitting, calibration,
qualified identity/C5/lambda and sealed human replication remain unresolved.

### Full G8 reference and two acceleration checks, 2026-10-06 evening (NY)

The older eight-worker scalar-Gamma evaluation22077927 completed54:06/exit0.
All eight genes returned, LL=-531.675870331499. Root independently rederived
the ascending-gene LL/29-score/29x29-Hessian sums, exact inactive-gene zeros,
nonnegative accumulated allowances and source/input hashes. All passed.
Reference: REC `repair_p3_gamma_full_parallel_22077927/`;
review: REC `repair_p3_gamma_full_parallel_review_20261007T020139Z/`.
This is one fixed B-FIX call12 point, not a converged fit or independent H test.

New `repair_p3_gamma_batched_model.py` retains original v=0 calls and inherited
fit methods, uses checked batched Gamma conditional arithmetic with local GL
in eight positive-v workers, and explicitly embeds each15-coordinate error
array into29 coordinates before summation. Parent BB remains original.
Claude's two read-only reviews found no coordinate/summation defect and
identified failure-classification issues. The final source scans every failed
worker for integrity errors before numeric refusals, preserves only the
successful prefix before the earliest failed gene, checks state on failure,
and recognizes only explicit numerical messages from checked helper origins.
State errors with attached ladders remain integrity errors. Conditional shapes
are checked separately from arithmetic. Missing helper counts are reconstructed
from first-row GL batches entered, with raw missing counts retained. Every
listed arithmetic message was checked against its defining source. All helpers
and wrapper must share one archive. Python3.10.18 was verified; the new checker
also uses an explicit ascending LL loop to avoid a later-version sum change.
Reviews: REC `claude_gamma_batched_model_review_20261007T015918Z/` and
`claude_gamma_batched_model_v2_review_20261007T020358Z/`.

Full comparison22078827 is running as one individual1501,cpu/nslab,
8CPU/8GiB/90h job. All14 boundary, inherited-method, coordinate/error-array
and failure-priority checks passed before entering the positive-v integral.
Every gene must match the complete scalar reference at LL1e-9/g1e-8/H1e-5,
with identical quadrature geometry, exact inactive zeros, independently
rederived embedded error sums and unchanged inputs. No fit is running here.
Launch: REC `repair_p3_gamma_batched_full_v3_launch_20261007T020734Z/`.
Model SHA256: `73ea34cbf52f12927d36a05b6ef288a72b2020bf63c3163d2b0c0ad9d8b114ab`.
Checker SHA256: `ca735174f15308c9a7e776692bcfe56ef62247b094da3caf0288d23fc1b26448`.
Earlier prepared launch directories were not submitted and remain preserved.

A separate speed hypothesis now has positive fixed-point evidence. Job22078807
completed2:07/exit0,1CPU/4GiB/90h,cpu/nslab, peak process RSS82324KiB. It uses
the original adaptive GH mode search, derivative formulas and native conditional
assembly, replacing only BB arithmetic with the checked scalar Gamma primitive.
All prespecified20/40/80/160 orders returned full G8 LL/29-score/29x29-Hessian.
All four passed strict absolute LL1e-9/g1e-8/H1e-5 and owner LL1e-6,
g-atol1e-8/rtol1e-6 comparisons with22077927. Every gene's LL also agreed
within2.93e-12. All orders are retained; none is admitted for fitting.

| GH order | Seconds | Absolute LL difference | Maximum score difference | Maximum H difference |
|---|---:|---:|---:|---:|
|20|11.1013|9.09e-13|2.09e-9|4.75e-6|
|40|17.3879|3.07e-12|1.57e-9|2.19e-6|
|80|31.2334|4.55e-13|9.99e-10|1.74e-6|
|160|58.8997|2.27e-12|1.23e-9|4.19e-7|

Root rederived every saved difference and nominal pass flag from actual arrays,
checked dimensions/finiteness/counts/source hashes, and retained all outputs:
REC `repair_p3_gamma_native_GH_review_20261007T020912Z/`.
Launch/output: REC `repair_p3_gamma_native_GH_launch_20261007T020554Z/` /
`repair_p3_gamma_native_GH_check_22078807/`.
Checker SHA256: `d3ee6ce7ba47d46e37316a5bf3d10f34a11f1f9b75fdd709666005c5592355cf`.
These are different jobs/nodes and one saved point, not a controlled whole-fit
speed claim. Adaptive GH here is distinct from the previous fixed-prior Decimal
GH20 reference and the negative mixture-GH experiments. Next is checking the
other frozen points and per-gene derivatives before any corrected observed fit.

At02:01UTC,92/96 Decimal chunks passed collection,4 running,0 refusals:
REC `repair_p3_decimal_GL_parallel_collection_20261007T020113Z/`.
The four workers have progressed into110-digit evaluation;80-digit work is
retained. No mixed-curvature result is formed before all96 complete. Older
sequential full G8 and scalar-local Hessian jobs remain active. No jobs were
cancelled/restarted, actual/sealed outcomes read, or human replication claimed.

Full batched update02:11UTC:22078827 completed3:36/exit0. Positive G8
evaluation211.06543s; job process peak RSS309836KiB. All eight genes returned
and passed every comparison, with identical per-gene quadrature geometry.
Total LL is exactly the scalar reference; maximum full score/H differences
are1.1102230246251565e-16/2.2737367544323206e-13. Worst per-gene score/H
differences are8.881784197001252e-16/2.2737367544323206e-13. Root independently
rederived all eight local-to-global error mappings, every embedded allowance,
all ascending LL/g/H/error sums, inactive zeros, all14 preflight probes and
source/input hashes. All passed. Full allowances remain empirical estimates,
not certified error bounds. Review:
REC `repair_p3_gamma_batched_full_review_20261007T021141Z/`.
The earlier eight-worker scalar positive evaluation took3245.876s: about15.4x
recorded ratio across different jobs/nodes, not a controlled whole-fit claim.
This completes the fixed-point batched G8 comparison. Native adaptive GH's
11.1s result still requires other frozen points and per-gene derivative checks;
no observed-fit, B3 calibration or independently validated human result exists.

### Native adaptive GH across frozen vectors, 2026-10-06 late (NY)

The previous turn made progress: full G8 batched GL and all four native GH
orders completed the saved call12 comparison. This turn extends that evidence
to every gene at all four prespecified B-FIX vectors, then rechecks call12
per gene. Job22079062 completed5:15/exit0,1CPU/4GiB/90h,cpu/nslab. All40
gene/parameter points and112 integral evaluations returned. The16 zero-variance
points use original BB arithmetic and nuisance-only28x28 Hessians; the24
positive points each retain20/40/80/160 orders and29x29 variance-v Hessians.
No new order is selected and no optimizer is run by this check.

Root independently rederived80 scalar likelihood comparisons,32 call12
per-gene LL/g/H comparisons,72 per-gene adjacent-order refinements, all14
ascending full-model sums and all9 full-model refinements. All nominal criteria
passed; every array is finite in its applicable coordinates, inactive-gene
derivatives are exactly zero and all source/fixture hashes match. Maximum
independent scalar LL difference=1.54330e-11. Maximum adjacent-order per-gene
score/H difference=2.06382e-9/9.10319e-6; maximum call12 per-gene H difference
versus scalar Gamma GL=2.58195e-6. These checks do not certify full independent
Hessian accuracy or floating/quadrature error bounds. All older candidate
differences and the three theta4 scalar-GL refusals remain preserved.
Review: REC `repair_p3_gamma_native_GH_family_review_20261007T022502Z/`.
Launch/output: REC `repair_p3_gamma_native_GH_family_v2_launch_20261007T021744Z/`
and `repair_p3_gamma_native_GH_family_22079062/`.
Checker SHA256: `2ed40e20d1b30239a8784de8ecd28578a0f720353ea803e60af5ef8341e6078c`.

Claude's read-only review identified summary flag typing and missing
unavailable/full comparisons. Before submission, flags were explicitly cast to
Python bool, reference derivative dimensions/finiteness were checked before
numerical comparisons, missing full refinements were explicitly counted as
failures, and each call12 full sum was compared with the saved full answer.
Review: REC `claude_gamma_native_GH_family_review_20261007T021447Z/`.
The mode-search algorithm remains original, but its probability/derivative
values use corrected Gamma BB arithmetic at positive variance. This is not a
claim that old rounded mode values remain byte identical.

New `repair_p3_gamma_native_model.py` scopes that positive-variance arithmetic
to the original native marginal call and restores original BB in a finally
block. The original null/v=0/negative-v path and all12 fitter/design methods
are inherited. The positive default remains20 nodes; only the checked
20/40/80/160 family is exposed by this synthetic helper. It captures and
checks the Gamma callable identity, parameters, caller writeability, rows and
original evaluator callables. Event-output failures are integrity errors.
An initially proposed nonfinite-result exception was removed after review to
preserve the original returned numerical values and optimizer behavior.

Source review: REC `claude_gamma_native_model_review_20261007T021956Z/`.
Root inspected the actual original likelihood: score-only marginal returns
two items (not H=None); inherited fit methods build no sub-models and contain
no broad RuntimeError/Exception catch; the only cached inclusion term uses
unchanged fixed homozygote BB, not the patched heterozygote BB. Original
mode/node helpers are local to the pinned native marginal source. Older
22077922 comparisons already check native/Gamma agreement at v=0; the new
helper preserves original boundary arithmetic exactly.

Job22079139 is running,one1500,cpu/nslab,1CPU/4GiB/90h. It requires all112
scoped-model results to match the completed native family byte for byte,
checks full original zero-boundary bytes including unavailable variance-H
entries, and tests numerical, integrity and output-failure restoration,
Gamma-callable mutation rejection, score-only behavior and recovery.
Launch: REC `repair_p3_gamma_native_model_launch_20261007T022416Z/`.
Model/checker SHA256:
`b7da2dc1b5c3b928d89f1a3626f8cb5be95c6365d26da98173d0449d956b9da8` /
`b4c42ddf3ea2b1a380ffef04303e1e655f27779b815ede06ecff5c3e60edd8ef`.

Older scalar-local22077951 completed1:06:56/exit0. Its2076 conditional
comparisons remained byte identical. Score/H integrations returned in
1672.61350s/2327.88473s. Against the batched local H result at the same point
and all9 identical geometries, LL/g/H differences are9.95e-14/6.25e-12/1.80e-8,
within nominal criteria. Review:
REC `repair_p3_gamma_local_completed_review_20261007T021528Z/`.
The sequential G8 comparison22077928 remains live with genes0 and1 saved.
Decimal collection02:21UTC remains92/96 complete,4 running,0 refusals:
REC `repair_p3_decimal_GL_parallel_collection_20261007T022112Z/`.
No old job was restarted. P3 resample points, the native-model state checks,
observed fitting, B3 calibration and qualified human replication remain separate
requirements. No actual/sealed outcome, authority or release was changed.

Scoped evaluator update02:31UTC:22079139 completed5:44/exit0, peak batch
RSS62056KiB. All112 native-model results matched the saved family byte for
byte, all16 original full zero-boundary checks passed, and all16 state probes
passed. Probes cover numerical refusal propagation, BB restoration after
numerical and integrity errors, integrity taking priority over a numeric
exception, recovery and cleared failure metadata, changed Gamma-callable
refusal, exact original boundary calls while Gamma is disabled, event-output
error classification/restoration, score-only two-item interface and byte
agreement with the H call, inherited methods and unchanged final rows/BB.
Root checked the exact112 prescribed case/gene/order records, all16 probe
results, all source/input hashes and the executed byte-comparison code.
This equivalence check retains comparison records rather than duplicating raw
subject arrays; independent numerical claims use the separately retained
family arrays and root rederivations. Review:
REC `repair_p3_gamma_native_model_review_20261007T023146Z/`.
The scoped evaluator is now checked for this frozen synthetic family. Next are
the two saved P3 resample points before the original synthetic observed-fit
recipe. This is not observed-fit convergence, full independent curvature,
B3 calibration, human replication or model/paper adoption.

### P3 native-mode limitation, 2026-10-06 late (NY)

The previous goal turn made progress: the frozen-vector GH family and scoped
evaluator checks completed. The next frozen P3 check22079693 completed all8
requested evaluations, exit0. Four returned and passed; four retained the
original ValueError `nonconcave mode in adaptive quadrature`. All refusals are
at failed0gene2, v=0.9234818604998617, before native GH integration. Increasing
the node order cannot resolve this unchanged mode-search step. The diagnostic
job completing is not a passing numerical evaluation at that point.

All20/40/80/160 orders returned at positive18, v=0.00711885901713953. Maximum
LL/score/H differences against the frozen GL candidate were6.92e-11/6.93e-11/
2.58e-8. All mixed r2–omegaS and r2–v entries passed, as did the three adjacent
order comparisons. Root rederived the actual173-count SHA256 for both draws,
all eight dispositions, all positive18 array differences and mixed entries,
the three refinements and source/input hashes:
REC `repair_p3_gamma_native_P3_review_20261007T024625Z/`.

Both draws were generated with original BB and verified before any Gamma helper
was imported or constructed. Their counts and full row keys are retained in
REC `repair_p3_gamma_native_P3_check_22079693/draw_reconstruction.json`.
The original five-case family and observed call12 results remain valid within
their tested scope; they do not cover the refused P3 point. No fast native
method is admitted for every P3 point and no observed-fit job was submitted.
Launch/checker:
REC `repair_p3_gamma_native_P3_v2_launch_20261007T023737Z/`,
SHA256 `265f9fd5b3907c173c75ad74e40069d9a2ed442e2f37c2cd6de7275bdd79c1b7`.

Claude's pre-run source review led to earlier original-BB identity capture,
generation of both count vectors before Gamma setup, explicit parameter and
variance checks, actual model-row snapshots, ArithmeticError retention for all
orders, strict True admission flags and preserved nonfinite diagnostics:
REC `claude_gamma_native_P3_review_20261007T023412Z/`.

Next are two individual batched Gamma local-GL evaluations at these same points,
using the completed exact count vectors, unchanged grid/domain/order rules and
caps, LL1e-9/g1e-8/H1.25e-6. This tests the integration route that already passed
full G8 call12, without requiring a concave native search result. Both points
retain comparisons to the existing full GL arrays and both mixed entries;
local15-coordinate error arrays are embedded into29 coordinates before checking
allowances. Fitting remains pending these numerical results and a suitable
integration policy. The failed native outcomes will remain visible.

The owner reply still ends at design-answer §11 and supplies no frozen B3 plant
allocation/amplitudes or coverage-target clarification. An asynchronous user
question now requests those settings; it does not pause independent numerical
work. C5/lambda, qualified participant metadata and untouched human validation
remain unresolved. No calibration assumptions or actual/sealed data are added.

### Both P3 batched integrals passed, 2026-10-06 late (NY)

Jobs22079995/22079996 completed5:23/3:46, both exit0, as individual1500/1501
cpu/nslab,1CPU/4GiB/90h allocations. They evaluate the complete unchanged
local-GL rule family with checked Gamma conditional batches, at failed0gene2
and positive18 respectively. The large-variance point returned despite the
native mode-search refusal. Original BB remains untouched; input rows are the
exact cached synthetic draws with the full173-row key mapping checked.

| Frozen point | Integral seconds | Rules / nodes | LL difference | Maximum g difference | Maximum H difference |
|---|---:|---:|---:|---:|---:|
|failed0gene2|320.0566|9 / 53760|3.88e-11|1.88e-11|4.80e-10|
|positive18|222.2799|9 / 51072|6.56e-11|5.05e-11|1.49e-8|

Every nominal comparison and LL1e-9/g1e-8/H1.25e-6 accumulated allowance passed.
Root independently rederived every array difference, mixed entry, each15-to29
error mapping and all allowances; checked exact inactive zeros, finite29-score
and29x29-H arrays, settings and source/input hashes. Both settings records show
unchanged default LL/g rules and all caps, with only the H allowance tightened
from1e-5 to1e-5/8. Independent full-Hessian/rigorous-error acceptance is still
not claimed.

Mixed r2–omegaS/r2–v differences are2.22785e-12/-3.98986e-16 for failed0gene2,
and-6.82766e-13/2.34891e-12 for positive18. These are in the original variance-v
coordinates. All four native refusals remain in the earlier record; the GL
results do not turn those native attempts into successes.
Root review: REC `repair_p3_gamma_batched_P3_review_20261007T030158Z/`.
Launch: REC `repair_p3_gamma_batched_P3_v2_launch_20261007T025355Z/`.
Checker SHA256: `f140752fa8fe11766f4a4a2cd0f6740c5aaedb4a36ce03a9b0eab9fde7ae0d48`.
Outputs: REC `repair_p3_gamma_batched_P3_failed0gene2_22079995/` and
`repair_p3_gamma_batched_P3_positive18_22079996/`.

Claude's review preceded the final submitted checker. The changes explicitly
check both target-gene identities, reference shapes/finiteness, all four
fixture-file pins and the kernel pin; record actual/default Settings; check
both imported-original BB identities; validate count support; and retain any
nonfinite report-only diagnostics as explicit tags while requiring finite
returned numerical arrays. Actual local-GL source contains no native marginal
or adaptive-mode call and accepts the shared settings by attributes. The
1e-9/1e-8 comparison thresholds are retained, not relaxed to obtain a pass.
Review: REC `claude_gamma_batched_P3_review_20261007T024829Z/`.

Next is the original synthetic observed-fit recipe, retaining both original
warm/cold starts and numerical failure dispositions. The checked native Gamma
and batched-GL backends can be compared in separate jobs; an automatic hybrid
or altered mode search has not been admitted. A native nonconcavity refusal is
retained as a refusal, not interpreted as an inaccurate returned likelihood.
Any accepted fit still needs its native convergence flag, original projected
score cutoff and final/restart/start comparisons. Numerical success is separate
from independent B-FIX draw agreement, B3 calibration and qualified human
replication. The five-gene scalar-quadrature requirement also remains binding.

Decimal collection02:49UTC is92/96 complete,4 running,0 refusals:
REC `repair_p3_decimal_GL_parallel_collection_20261007T024916Z/`.
No mixed-curvature truth is formed from partial chunks. The old sequential G8
job remains live. The owner B3 question remains pending; no missing plant
settings were inferred from fixture values. No actual/sealed outcomes read,
statistical calibration, model validation or Resource adoption claimed.

### Separate original observed fits launched, 2026-10-06 late (NY)

Native Gamma job22080548 and batched-GL job22080549 are RUNNING, as separate
individual cpu/nslab jobs1500/1501. Headers: native1CPU/4GiB/90h; GL8CPU/8GiB/90h.
Both use the exact original synthetic173-row/60-participant B-FIX counts,
original converged null warm vector and inherited fit/fit_marginal methods.
All original starts, bounds, tie, null-snap, restart and stopping rules remain.
No automatic backend fallback is used. Original native nonconcavity refusals
remain possible; a returned job does not by itself establish fit convergence.

Runner: scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_observed_fit.py,
SHA256910351726be9279d92bbd484dffce4934f908ed6c26b7957493dba6cf2a4b7a7.
Launch: REC `repair_p3_gamma_observed_fit_v2_launch_20261007T031426Z/`.
Claude review: REC `claude_gamma_observed_fit_review_20261007T030950Z/`.
AST, both shell syntax checks and both static input inspections passed before
submission. Static inspection also checked the actual baseline JSON shape and
complete source coverage across the completed numerical experiments.

The review led to explicit source-coverage enforcement and a three-state final
40-node comparison: returned, unavailable after numerical refusal, or zero
boundary using the same original evaluator. Zero-boundary agreement is not an
order-refinement check. Unknown ValueErrors are not silently labeled arithmetic
refusals. Actual GL source uses a context-managed process pool that joins before
returning; unexpected state errors in final diagnostics still stop the job.
Each fit is saved before diagnostics, native convergence flags are unchanged,
and projected scores are independently recomputed with the original active-set
rules and strict1e-6 cutoff. Starts, final LL/kappa/LR/tau differences and all
refusals are retained. Missing terminal files after scheduler kill remain
incomplete. Neither fit has returned yet at this entry.

Decimal collection03:15UTC remains92/96 complete,4 running,0 refusals:
REC `repair_p3_decimal_GL_parallel_collection_20261007T031500Z/`.
No partial mixed-curvature result is interpreted. The sequential G8 scalar
job22077928 is still running; no restart or cancellation. A new five-gene Gamma
likelihood-only checker is being reviewed against all15 unchanged independent
scalar references. Its first static inspection caught archived-source path
keys; basename uniqueness plus exact digests corrected that check, with the
failed pre-submission attempt preserved. No five-gene Gamma result yet.

B3 plant allocation/amplitudes and coverage target remain requested from the
specification owner; no response or scientific settings are inferred. C5/lambda,
qualified human identities,50-draw agreement, B3 and untouched human replication
remain open. No actual/sealed outcomes read, statistical calibration, biological
validation or Resource adoption claimed.

### Both null fits checked; five-gene Gamma check launched, 2026-10-06 late (NY)

Both22080548/22080549 returned the same observed null LL-529.9352754997254.
Difference from the original baseline is1.13687e-12, kappa difference-9.62563e-14.
Root independently rederived the28-coordinate bound/gradient active mask and
projected score from each saved diagnostic: maximum1.8007784444162844e-7<1e-6.
Both native convergence flags are true, retained start is warm, and cold-start
convergence is false in both records. The unsuccessful cold start is retained.
These checks accept only the null; both marginal fits remain running.
Root review: REC `repair_p3_gamma_observed_null_review_20261007T031952Z/`.
Native positive score calls have returned; GL's first positive full8-gene
score call returned in96.7120s. These are evaluations, not completed fits.

Five-gene Gamma checker job22081217 was submitted as individual1502,
cpu/nslab1CPU/4GiB/90h after Claude review, AST/shell and15-point static checks.
It uses the exact original seed20260923, row counts5/25/75/200/500 and all three
tau multipliers. Gamma replaces only heterozygous BB arithmetic during native
GH20/GH40 calls and is restored in finally; original fixed homozygote masses,
selection, predictors and native mode algorithm remain. No G8-wrapper claim.
All independent direct-scalar references are reused unchanged, with both
original/expanded tail estimates and refinement checked before allocation.
All raw fixtures, parameters and original source digests are checked. Saved
absolute differences accompany the prespecified relativeLL1e-8 and GH20/GH40
absolute1e-8 criteria. Scalar arithmetic/quadrature errors are not certified;
no independent score/Hessian, fit or calibration claim follows from a pass.

Runner: scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_five_gene_check.py,
SHA25699585404c43650a3619842cc907a9639195fe9f5b973f8242f83a0745fbda31d.
Launch: REC `repair_p3_gamma_five_gene_v3_launch_20261007T032034Z/`.
Review: REC `claude_gamma_five_gene_review_20261007T031701Z/`.
No five-gene result yet. B3 scientific settings remain pending. No actual/sealed
outcomes read, biological validation or Resource adoption claimed.

### Native observed fit and all15 Gamma likelihood checks complete, 2026-10-06 late (NY)

Previous goal turn made progress: two original observed fits and the Gamma
five-gene likelihood check were submitted after Claude reviews. This turn
collected new terminal evidence. Native observed job22080548 completed0:0 in
8:55, peak57660KiB. Both marginal warm/cold starts converged; their LLs are
-529.9352754997254/-529.9352754997226, within the original1e-6 tie tolerance.
Warm was retained and snapped to the exact null: v0, LR0. Original baseline
LL difference1.1368683772161603e-12, kappa difference-9.625633623500107e-14,
tau difference0. Null and marginal projected-score maxima are both
1.8007784444162844e-7<1e-6. Independent root recomputation checked active masks,
all28/29 scores, finite28x28 nuisance Hessians, exact null snap, both starts,
comparison arithmetic and unchanged source/input digests. The 40-node comparison
uses the same original zero-boundary evaluator and is explicitly NOT positive
order refinement. Native refusal at the earlier failed0 fixed point remains.
Root: REC `repair_p3_gamma_observed_native_review_20261007T032510Z/`.
Output: REC `repair_p3_gamma_observed_native_22080548/`.

Five-gene Gamma job22081217 completed0:0 in5:00, peak75900KiB. All15 requested
points returned both GH20/GH40 integrals (30 total), with no refusals. Maximum
absolute LL difference from the unchanged direct scalar reference is
9.094947017729282e-13; maximum GH20/GH40 difference is the same. All prespecified
relativeLL1e-8 and GH20/GH40 absolute1e-8 criteria pass. Root independently
rederived all30 differences, all15 refinements, exact Cartesian coverage,
finite17-score arrays and all source/input digests. This completes this Gamma
native likelihood-only family, not independent gradient/Hessian accuracy,
full50-draw agreement, calibration or human validation. It does not make the
old non-Gamma GL candidate's three refusals disappear.
Root: REC `repair_p3_gamma_five_gene_review_20261007T032730Z/`.
Output: REC `repair_p3_gamma_five_gene_check_22081217/`.

The batched-GL observed fit22080549 remains running. Its first seven positive
score evaluations returned; observed times64–104s motivate profiling before
any further optimization. Profiling job22081450, individual1500cpu/nslab,
1CPU/4GiB/90h, runs the unchanged already checked failed0gene2 P3 checker under
stdlib cProfile. No numerical candidate is changed, no live fit restarted.
Launch: REC `repair_p3_gamma_batched_profile_launch_20261007T032359Z/`.
The five-gene and native fit collectors do not depend on this performance job.

The next native-Gamma P3 runner extracts the complete original draw loop body
by AST, preserving original conditional draws, section4 starts, exception
classes and failure-as-exceedance rule. All50 individual scripts passed static
input inspection and shell syntax checks; they are not submitted yet, pending
Claude review. No held independent comparison outcomes were read. Each returned
fit gets its original projected-score and native40 comparison, kept separate
from original statistics; refusals and integrity failures remain distinct.
Staged: REC `repair_p3_gamma_draw_launch_20261007T032825Z/`.
B3 settings, C5/lambda and qualified human replication remain open; no Resource
adoption or biological model validation is claimed.

### Draw-review corrections and measured runtime costs, 2026-10-06 late (NY)

Claude's first draw-runner review identified missing refusal traces outside
positive integration, post-fit diagnostics potentially obscuring a completed
original draw, and an implicit loader-fixture tie. No draw jobs were submitted.
The revised runner logs conditional_draw and full inherited fit_marginal
exceptions through recording-only proxies and rethrows unchanged. It saves
recipe_summary before diagnostics and retains diagnostic failures separately;
unknown caught exceptions are explicitly unclassified rather than admitted
numerical refusals. The native finite-point Runtime guard stays unchanged:
invalid fitting state remains an incomplete job, not a failed scientific refit.
Original AST body and exception classes remain unchanged. Each launch claims
its draw once, and exact drawn counts are checked against reconstruction from
the original baseline null before fitting. The original baseline null is now
pinned; no held owner draw outcomes are read. Loader.FIXTURE is explicitly tied
to the pinned fixture. Fitting-refusal details are captured before diagnostics
can reset the native helper's last-error field.

Revised runner SHA256cb6fdc744e12d020bb4dbd4f629c6543b4330ac010ebc75724ffc7bd1288212f.
Staged launch: REC `repair_p3_gamma_draw_v2_launch_20261007T033424Z/`.
All50 static input checks and all50 shell checks pass; still not submitted,
pending focused follow-up review. The per-draw GH40 and projected-score results
will remain separate from unchanged original LR/exceedance/failure records.
With observedLR0, all original P3 exceedances are already mathematically true;
these50 fits test the required draw-LR/failure agreement, not calibration/power.

Profile22081450 completed0:0 in6:52, peak67696KiB. It reran the unchanged
failed0gene2 checker; nominal comparisons, allowances and fixed-input checks
still passed. Stdlib cProfile recorded411.0096s/99,726,969 calls. Batched rising
products took215.9388s cumulative (165.4759s self); BB batch assembly290.0868s
cumulative, and local integration403.3528s. These nested times must not be
summed. File read_bytes7.4559s cumulative and SHA2561.8203s self account for
about9.3s; source-file I/O is not the main measured cost. No source-hash removal
or numerical-candidate optimization was made. Profiled timings are not an
unprofiled speed estimate. Source/math changes would require new comparisons.
Profile and machine-readable top costs: REC
`repair_p3_gamma_batched_profile_launch_20261007T032359Z/profile_summary.json`.
Native draw parallelism remains the immediate acceleration. GL observed fit,
sequential scalar G8 and four Decimal chunks are still running. No biological
validation, calibration or Resource adoption claimed.

### Independent mixed curvature complete;50 native draws running, 2026-10-06 late (NY)

All96 Decimal GL chunk jobs completed0:0 with no refusals. Final four elapsed
times3:38:33/3:38:53/3:40:04/3:43:25; peakRSS64776/64808/65840/64708KiB.
Collection: REC `repair_p3_decimal_GL_parallel_collection_20261007T033640Z/`.
Root independently rehashed all192 raw80/110-digit moment files, checked all
96 actual scheduler terminal states, reconstructed every panel coverage and
six raw-moment sum, then rederived both covariance formulas at all12 combined
rule/precision points. No chunk covariance was summed and whole-line analytic
tail bounds were not multiplied by16. Root:
REC `repair_p3_decimal_GL_full_review_20261007T033840Z/`.

All fixed precision1e-30, refinement1e-8 and old-candidate comparison1e-5
criteria pass. Maximum precision differences8.37038e-76/1.14963e-76 for
failed0gene2/positive18; maximum refinement differences4.93137e-13/4.22919e-50.
Finest independent mixed r2–omegaS/r2–v values are:

| Frozen synthetic point | H(r2,omegaS) | H(r2,v) |
|---|---:|---:|
| failed0gene2 | -7.7770234981502963e-6 | -8.9727340818244200e-7 |
| positive18 | -3.2466128815343160e-4 | 2.0218169108601960e-3 |

Root additionally compared the actual checked Gamma-batched answers against
ALL six independent rule/precision results per point. Largest absolute mixed
entry differences3.6177396e-13/3.7104743e-13; all pass the fixed1e-5 criterion.
Root: REC `repair_p3_gamma_Decimal_mixed_review_20261007T034236Z/`.
These are the four disputed entries at two fixed synthetic points. They do
not certify all Hessian entries, all parameter values, floating/quadrature
error bounds, full fits or statistical calibration. Prior failed finite-
difference comparisons and native-mode refusals remain recorded.

Both Claude draw reviews were addressed before submission. The follow-up found
two remaining unlogged operations: original draw predictors and model
construction. Recording-only proxies now cover them too, and an invariant
requires every original fit=None record to have a captured exception. Known
native nonconcavity requires the final raising frame in pinned original source;
unknown caught exceptions remain labeled unknown. The unique draw claim is
written immediately before execution, after setup; source/schema/header checks
precede it. Actual original Model.__init__ accepts gene_count, and all headers
explicitly declare1CPU. Original fit/draw/restart/error rules are unchanged.
Reviews: REC `claude_gamma_draw_review_20261007T032744Z/` and
`claude_gamma_draw_followup_20261007T033405Z/`.

Submitted all50 original supplied P3 draws as individual22081569–22081618
jobs, cpu/nslab1CPU/4GiB/90h each (50CPUs/200GiB total). No array. All50 were
observed RUNNING, no failure.json or job failure, after submission. Ledger
appends and individual submission receipts were locked/fsynced; draw assignment
is fixed0–49 and no replacement/duplicate draw was submitted.
Launch: REC `repair_p3_gamma_draw_v3_launch_20261007T034152Z/`.
Runner: scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_draw.py,
SHA2569a3c978dccc685cf648afa35e30c419a4b4b2eb9944a6391c4fb1480b2d16b81.
All50 static checks, source/fixture inspections and shell checks passed before
submission. The complete original AST draw body is retained. Each draw checks
exact reconstructed counts against original-baseline-null counts before fit,
records exceptions before original catch, preserves conservative exceedances,
and saves original recipe statistics before separate projected-score/GH40
checks. Nonfinite fitting state remains an incomplete integrity failure;
diagnostic failure does not erase a completed original draw. No held owner
comparison results read. Returned50-draw results still need independent LR,
failure and applicable final numerical agreement; no such agreement claimed.

Native observed fit and15-point Gamma likelihood check are complete as above.
GL observed fit22080549 and scalar sequential22077928 remain live. B3 owner
settings, C5/lambda, qualified participant identities, independent human
replication and Resource adoption remain unresolved. No actual/sealed outcomes,
calibration, biological validation or adoption claimed.

### Both observed fits complete and44/50 draws collected, 2026-10-07 early (NY)

Previous goal turn made progress: native observed fit and all15 Gamma likelihood
points were checked, all96 independent Decimal chunks completed and were
independently recomputed, and50 individual native-Gamma draws were launched.
This turn adds a bounded collector, collects actual terminal draw evidence and
completes the GL observed-fit comparison. No held owner draw outcomes read.

GL observed job22080549 completed0:0 in48:39, peak304324KiB. Its two marginal
starts converge; warm/cold LLs-529.9352754997254/-529.9352754997187 are within
the original1e-6 tie tolerance, warm retained. Final v0/LR0 and the retained
28 parameters/LL/29 scores exactly equal the native job22080548. Null and
marginal projected-score maxima1.8007784444162844e-7<1e-6. Root independently
rederived both bound masks, projected scores, finite nuisance Hessians, exact
null snap and unchanged source/input digests in BOTH jobs; 40-node comparison
atv0 remains the same boundary evaluator, not positive order refinement.
Root: REC `repair_p3_gamma_observed_backend_review_20261007T040931Z/`.
Wall times8:55(native1CPU) and48:39(GL8CPU) are observations under different
allocations, not a controlled fitspeed benchmark. All prior refusals stand.

New scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_draw_collect.py,
SHA256ba4d57cc7a773056b01defedbedad1c56646eeaa0111811c41a20f497a2ed0a8,
uses only bounded JSON/source reads and actual scheduler state. Fixed roster,
executed sources, exact count-vector bytes, original LR/failed/exceedance rules,
fitting exception records and final diagnostics are checked independently.
It requires every draw job COMPLETED0:0 and complete output before emitting
an aggregate. Known failed scientific refits remain NaN-LR/exceed=True;
missing/integrity-failed jobs remain incomplete. Unknown caught exceptions are
reported separately and are not admitted as numerical refusals.

Claude review: REC `claude_gamma_draw_collector_review_20261007T034903Z/`.
The corrections retain nonfinite native40 output as a failed diagnostic rather
than aborting collection, mark order refinement not applicable atv0, bind the
diagnostic parameters/LL/score to the actual fitted LR point, and independently
project BOTH null and marginal saved scores using their original1e-8/1e-10
bound tolerances. The executed-source list is exactly the runner's eight
sources; terminal jobs missing output remain incomplete. Runtime/environment
and variance-v units are explicit. No independent fifty-draw agreement claimed.
AST and collections of actual saved synthetic outputs passed; these are
fixture results, not biological data or calibration.

Latest collection: REC `repair_p3_gamma_draw_collection_20261007T041042Z/`.
44/50 COMPLETED0:0,6 running,0 pending,0 unexpected collection/diagnostic issues.
23 recipe-successful fits are allv0/LR0 and pass the saved null/marginal score
checks; native40 order refinement is inapplicable at this boundary. The21
original failed refits include20 captured native nonconcavity refusals and
one returned nonconverged fit(draw11; null projected maximum2.1271003368761657e-6).
All are retained as original conservative exceedances. An interim failure
proportion is not interpreted: failures can finish earlier. No aggregate,
held-reference comparison or calibration claim is formed until the family is
complete. Six specific draw processes remain live; no duplicate/refit restart.

### Unchanged mixture-GH family with Gamma arithmetic, 2026-10-07

Native nonconcavity and the observed GL cost motivate a separate fixed-point
recheck of the earlier numerical family. All six equal/scale_sqrt/scale_precision
by wideFalse/True policies are retained at BOTH frozen P3 gene2 points:12 total.
This changes the old mixture experiment's BB arithmetic premise, not the frozen
statistical estimand. It remains a numerical revision of native20 quadrature,
with no automatic backend fallback, no policy selection and no fit admission.

Claude review: REC `claude_gamma_mixture_review_20261007T040005Z/`.
The review correctly found that tightening Htol also changes precision-policy
proposal weights. Before any submission, defaults were restored exactly;
G8 Htol/8 is now only a postcheck on returned allowances, not a rule/weight/
refinement/tail change. The first unsubmitted attempt remains preserved.
Further checks explicitly tie gene2/24rows and all28 parameters to the Gamma
reference; bind independent raw Decimal protocols/counts/parameters and moment
hashes to this exact point; select the unique rule2/digits110 reference by its
metadata; assert mixed coordinates(10,24,28); and prove the conditional global
BB lookup. Instrumented Gamma call counts verify actual arithmetic use.
Private BB substitution is restored in finally. Known kernel IntegrationFailure
or explicit Gamma nonfinite arithmetic refusals are retained; unknown errors
fail the diagnostic job instead of being labeled numerical. No cap relaxed.

Runner: scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_mixture_check.py,
SHA256ebef52fa5dea3b8833a3971d2b2e20b070699fea49e2ce05e7afd11d3541d1a3.
Launch: REC `repair_p3_gamma_mixture_v2_launch_20261007T040713Z/`.
All12 AST/static input checks and shell checks passed before submission.
Submitted12 individual cpu/nslab jobs22081797–22081808,1CPU/4GiB/90h each
(12CPUs/48GiB total); no array. Latest snapshot1completed/11running,0 job errors;
full12-point family must be retained before interpretation. Prior old-arithmetic
failures remain; no numerical policy is chosen from the new results.

Sequential scalar reference22077928 remains live and has saved genes0–4.
B3 scenario allocation/amplitudes and coverage target still await the owner;
C5/lambda and qualified participant identities remain open. No actual/sealed
outcomes, biological validation, calibration, Resource adoption or paper claim.

Latest draw collection: REC `repair_p3_gamma_draw_collection_20261007T041734Z/`:
47/50 completed0:0,3 running,0 unexpected record/diagnostic issues. Mixture-GH
snapshot7/12completed,5running,0 job failures; no subset policy interpretation.
Both observed fits remain checked and complete; the scalar sequential job is
still running. No held owner comparison outcomes or biological data read.

### Complete Gamma mixture-GH family, 2026-10-07

All12 jobs22081797–22081808 completed0:0. Independent root inspection is
REC `repair_p3_gamma_mixture_review_20261007T042724Z/review.json`.
Eleven policies/points explicitly refused the unchanged order-refinement
criterion. Only positive18/scale_precision/wideFalse returned; its LL/score/
Hessian differences from the checked Gamma GL reference were1.10845e-12,
4.46576e-11 and9.78253e-9. Its nominal and posthoc G8 allowance checks passed.
No policy returned at both points. This family does not qualify as a general
replacement; all12 outcomes are retained and no passing policy is selected.

The native50 collection at REC
`repair_p3_gamma_draw_collection_20261007T043515Z/` remains47completed/3running,
with0 unexpected record/diagnostic issues. The21 original failed refits remain
conservative exceedances; the26 returned successful fits are v0/LR0. These are
interim fixture results, not a final failure proportion or calibration.
Draws18,25,30 and the scalar sequential reference remain live.

A separate pure-GL full-draw feasibility pilot is being prepared for original
draws0 and18. It preserves original starts, bounds, tie/null rules, draw body,
counts, and numerical-refusal semantics. It will neither replace native draws
nor support fifty-draw agreement. Native Gamma five-gene checks have passed;
a corresponding Gamma GL five-gene check remains necessary before generalizing
that backend. No missing B3 plant or coverage setting is inferred from this
pilot, and no actual or sealed outcomes are read.

Additional root check: REC
`repair_p3_gamma_mixture_mixed_review_20261007T044850Z/` independently
rederives the passing mixture result's two cross-derivatives from the saved
110-digit rule2 raw moments. Differences are3.884376544374679e-13 and
-9.619787605386065e-13, both within the fixed1e-5 criterion. This supports
only those two entries, not full-Hessian accuracy or family acceptance.

Seven scalar sequential gene outputs are now saved. Root inspection at REC
`repair_p3_gamma_sequential_partial_review_20261007T045233Z/` compares their
actual LL,29-score and29-by29-Hessian arrays with the corresponding parallel
workers. All seven meet the fixed comparison tolerances; they are not byte
identical. Maximum absolute LL/score/Hessian differences across these seven
genes are3.48166e-12/1.44880e-11/1.80444e-8. The eighth gene and terminal total
remain pending; no final pair agreement is claimed.

Claude's first GL-draw review identified useful corrections before submission:
save GL final scores before any native40 comparison; reject unknown caught
exceptions as integrity failures; save completed-evaluation progress; and
check the eight-CPU/eight-GiB allocation before fitting. These are implemented
in the separate GL pilot runner. Its pins/claims use their own launch directory.
The GL observed null parameters must exactly equal the checked NEW native Gamma
observed null. They differ slightly from the older original baseline; instead
all draw counts must still exactly equal the older original construction.
No original native runner or numerical helper was edited. The modified pilot
passes AST, both input inspections and bash-n; Claude follow-up is pending.

### Two pure-GL draw pilots submitted, 2026-10-07

Claude follow-up at REC `claude_gamma_GL_draw_followup_20261007T045203Z/`
confirms the six substantive corrections and the inherited fit methods' refusal
propagation. Its final small suggestions are implemented: original draw-body
compilation in inspect; exact original-baseline count agreement before claiming
a draw; enforced separate GL pins/claims directory; and completed-job positive
marginal-call count. The original loop was independently inspected: prescribed
ValueError/FloatingPointError/LinAlgError catches, raw fit/LR save, conservative
failed-draw exceedance and supplied context names are unchanged. Unknown caught
errors remain explicit job failures. No helper or live native runner changed.

Runner: scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_GL_draw.py,
SHA25666342d0caf36d3284ab515384adf3b43a063c6140d146b91112d39379b3edc38.
Final launch: REC `repair_p3_gamma_GL_draw_pilot_v3_launch_20261007T045626Z/`.
AST, unchanged original draw-body compile, both pinned input inspections,
source hashes, exact new native/GL null parameter identity and bash-n passed.
The first two launch copies were not submitted. Submitted only two individual
cpu/nslab jobs: draw0=22082078, draw18=22082079, each8CPU/8GiB/90h; total16CPU/
16GiB, no array. These are a separate backend feasibility test of two existing
synthetic draws. They do not replace any native draw or admit a numerical
fallback, fifty-draw agreement, calibration or human validation. Original
native draws18/25/30 and scalar reference22077928 continue unchanged.

Startup verified from live SLURM: both pilots RUNNING,0pending/0 job errors.
Both exact-count prechecks completed and fit_started progress records exist;
claims, protocols and all17 archived source hashes were independently checked.
No pilot fit has completed yet. The native family remains47/50complete with
3running; scalar reference has7/8gene outputs saved. Technical documentation
diff checks pass. Resource-scope validation was not run: no paper/figure/
release/Catalog/portal authority was changed.

### Separate GL-pilot collection, 2026-10-07

Added scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_GL_draw_collect.py,
SHA2564d6a979e32862176a00c0eba7f505c39d32b0dd7515495bf75a9ca11b37c2c39.
It reads only the fixed two-pilot launch and its pinned synthetic inputs,
checks archived sources, claims and actual scheduler state, and rederives
original conservative-refit/LR/count arithmetic when completed records exist.
GL final score checks are separate from the optional non-gating native40
comparison. The null score is explicitly saved, not independently recomputed;
its original1e-8 projection tolerance is preserved, with1e-10 for marginal.
Neither pilot is pooled with native50 or converted into a calibration p-value.

Claude supplied-text review: REC
`claude_gamma_GL_collector_review_20261007T050121Z/`. Primary corrections
include explicit scheduler/terminal-output errors, separate GL recipe failure
counts and fitting-integration-failure details, a lower-bound label on capped
progress counts, diagnostic-inclusive positive-call counts, zero-boundary
allowance non-applicability, and reportable allowance exceedances. Native40
agreement arithmetic is rederived when finite; refusals/nonfinite comparisons
remain separate from GL score checks. Constructor refusals may retain absent
counts. Live progress uses one bounded byte snapshot and ignores an unfinished
last line. The per-file32MiB cap matches the runner; total reads remain128MiB.
The corrected review-source receipt distinguishes the actually supplied source
from subsequent primary corrections. No numerical helper or live job edited.

AST and actual live collection passed at REC
`repair_p3_gamma_GL_draw_collection_20261007T050612Z/`:0/2completed,
0 prescribed GL recipe failures observed,0 issues. Completed-fit branches are
NOT RUN because both jobs remain live. Latest live scheduler state: native
3running, GL2running, scalar reference1running;0pending/0job errors. Native
completion remains47/50; scalar outputs7/8genes. Measured batch peak RSS is
55820/55828KiB for the GL parents so far,57944/57872/57860KiB for the three
native jobs and60916KiB for the scalar job; the GL allocation remains8GiB per
job and worker peak memory is not yet measured. Downstream calibration still
awaits the B3 owner settings; no newer owner reply was found. No statistical
calibration, biological validation or paper adoption is claimed.

### Complete scalar/parallel Gamma reference pair, 2026-10-07

Scalar job22077928 completed0:0 in4:23:33; all8genes and the full result
returned. The parallel reference22077927 also completed0:0. Collection:
REC `repair_p3_gamma_full_pair_collection_20261007T051539Z/`.
Root independently recomputed each full likelihood,29-score vector and29x29
Hessian from all8 saved gene contributions:
REC `repair_p3_gamma_full_pair_root_review_20261007T051742Z/`.
Full absolute differences are3.410605131648481e-12 for LL,
9.379164112033322e-13 for the largest score entry, and
1.6065314412117004e-8 for the largest Hessian entry. Every gene and the full
sum meet the fixed comparison tolerances; accumulated numerical allowances
also meet the unchanged global limits. Arrays/reports are not byte identical.
This completes same-law evaluator agreement at one full synthetic point,
not independent full-Hessian accuracy, a fitted-model comparison, calibration
or biological validation. Scalar peak RSS:115340KiB from process inventory,
74108KiB from sampled SLURM batch statistics.

A remaining numerical check is being prepared independently of the live fits:
the fixed15-point Gamma GL family for G5/Z1/P16, with5/25/75/200/500 rows and
tau=.05*(.25,1,4). Existing independent selected-probability scalar quad
references (epsabs1e-12) and their uncertainty are reused, with exact raw
record/fixture/count/parameter/seed/source/environment checks. Each point
retains native Gamma GH20/GH40 and adds dense17-coordinate composite GL using
unchanged reviewed Gamma3 arithmetic and GL settings. All15 outcomes must be
retained. This checks the dense numerical family; it does not broaden the
G8-only batched conditional scope or supply independent gradient/Hessian
accuracy. No resampling, fits, human data or missing B3 plants are introduced.

Claude review: REC `claude_gamma_GL_five_gene_review_20261007T051244Z/`.
Primary corrections add pre-call G5/Z1/P16 and row-order checks; actual Gamma
BB invocation counts; explicit executed Settings; narrow known numerical
refusal messages; raw returned GL arrays before structural checks; separate
native and GL return flags and a GL-only relative-likelihood criterion. Row
is a regular frozen dataclass in the pinned original source. The original
fixture construction and native GH loop remain unchanged by AST comparison.

Submitted all15 individual Gamma GL accuracy points as jobs22082310–22082324,
1CPU/4GiB/90h each on cpu/nslab (15CPUs/60GiB total; no array). Runner:
scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_GL_five_gene_check.py,
SHA256ac62db5d5a6bae283ddd73bad09f6b6fa72ce5362282f6e6082bdcba26a6803b.
Launch: REC `repair_p3_gamma_GL_five_gene_v3_launch_20261007T051812Z/`.
All15 final pinned input inspections and bash-n checks passed before submission.
All15 actual protocols and six archived sources/job passed root hash/identity
checks after startup. Latest snapshot15running/0pending/0job errors;12points
have saved both native GH results, passed their fixed native/scalar criteria,
and entered GL. No completed Gamma GL accuracy point has been interpreted yet.
The first two preparation copies were not submitted. Existing numerical helper
sources, native draws and GL pilots were not edited or restarted.

Latest native collection: REC `repair_p3_gamma_draw_collection_20261007T051914Z/`:
47/50completed,3running,0issues. GL fit pilots22082078/22082079 remain live.
Scalar reference is terminal and no longer requires monitoring. Downstream:
finish these numerical/fixture comparisons, then address B3 calibration with
the owner's missing scenario allocation, plant amplitudes and coverage target;
qualified participant replication inputs remain pending. No calibration or
biological claim is inferred from numerical agreement. Technical documentation
checks pass; Resource-scope validation was not run because no paper/figure/
release/Catalog/portal authority changed.

### Three Gamma GL accuracy points complete, 2026-10-07

Added scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_GL_five_gene_collect.py,
SHA25653aa42086efa51672fca3961d9b71f8043f4d79ece3cadc7c09b03c1862f6e4a.
It requires the exact15-job roster and source/fixture/settings hashes, actual
COMPLETED0:0 state, raw scalar references, native GH20/GH40 checks, returned GL
arrays and per-point likelihood arithmetic. Refusals remain distinct from
accuracy failures. Missing/failed jobs cannot complete the family and produce
an explicit error status; collection returns nonzero on recorded errors.
No result is pooled into calibration, and G8 batched scope is not broadened.

Claude supplied-text review at REC
`claude_gamma_GL_five_collect_review_20261007T052757Z/` found no false-pass path.
Primary additions expose terminal errors, missing scheduler states, GL-only
family criteria, absolute per-point tolerances and native scalar differences.
Reported GL arrays must match the actual final integration report, empirical
refinement-plus-tail arrays must have full17/17x17 dimensions and meet the
unchanged executed Settings, and Gamma call counts must lie within the2–3
BB-support calls per row/node allowed by the pinned observed-plus-selection
code. These are numerical consistency checks, not certified error bounds or
independent derivative accuracy. Hessian asymmetry is descriptive only.

Actual final collection: REC
`repair_p3_gamma_GL_five_gene_collection_20261007T053310Z/`:3/15completed,
12running,0pending/0issues. Completed jobs22082310/22082315/22082320 are gene0
with5rows at tau.0125/.05/.2. All three pass the GL-to-scalar relative1e-8
likelihood rule and the native GH20/GH40 requirements. GL-minus-independent
scalar LL differences are8.881784197001252e-15,1.0658141036401503e-14,
1.0658141036401503e-14. Gamma BB call counts are376320/376320/408329.
Root separately rederived these likelihood comparisons from actual point
outputs and the original scalar records at REC
`repair_p3_gamma_GL_five_gene_partial_root_review_20261007T052951Z/`.
The larger25/75/200/500-row points remain pending; no full-family pass is claimed.
AST and actual completed-point collection pass. No numerical helper or live
runner changed; no new jobs submitted in this collection turn.

The native family remains47/50complete with3running; both separate GL fit
pilots remain live and have returned positive-variance evaluations, but no
pilot fit is complete. No job errors observed. B3 owner settings, calibrated
inference and qualified participant replication remain open. No Resource
claim/adoption changed, so Resource-scope validation was not run.

### Supported64-node batch comparison, 2026-10-07

Added scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_batched_p3_batch64_check.py,
current SHA2567d500ed9337f3e3637995cb9e2ed51c8d378ef6b4889704464bb26f21ac7a677.
This changes only the already supported requested GL node batch from32 to64
at BOTH fixed P3 gene2 points. All14 numerical sources, GL orders/grid/domain/
limits, parameters, count vector and tolerance settings remain unchanged;
current numerical helpers and live fits were not edited. The comparison pins
and checks the actual32-node protocol/settings/results, captures raw returned
arrays before diagnostic assertions, checks the exact numerical node path,
records raw differences and exact numerical equality, and preserves workspace
and count-axis limits. CPU/host information is recorded for timing context.
Claude review: REC `claude_gamma_batch64_review_20261007T054008Z/`.

Three checker attempts22082562/22082563/22082610 failed after returning the
integral because my assertion incorrectly required full64-node conditional
chunks on every row. The synthetic gene2 fixture has three depths2024/2158/
2466; the existing100000-entry limit correctly splits64-node requests into
[49,15]/[46,18]/[40,24]. No numerical precision criterion failed or was relaxed.
The corrected assertion checks the existing min(requested,100000//(n+1)) rule
and complete return of all requested nodes. Failed attempts remain intact and
are not counted as successful scheduler jobs or P3 resampling failures.

Positive18 raw arrays were recovered from22082610 without recomputing its
integral: REC `repair_p3_gamma_batch64_positive_raw_review_20261007T055545Z/`.
The job retains FAILED1:0 from its checker assertion. The other point had no
saved raw arrays in its first attempt, so only that terminal failed job was
rerun as22082668 (cpu/nslab,1CPU/4GiB/90h, no array). It completed0:0 and passed
the corrected count-capacity check and all numerical comparisons. Initial two
jobs and the positive diagnostic rerun used the same per-job resource header.
Original live fits and the fifteen-point family were not restarted.

Both numerical returns independently verified at REC
`repair_p3_gamma_batch64_pair_review_20261007T060216Z/`. Positive18 is exactly
numerically equal to32-node LL/score/Hessian arrays. Failed0gene2 differs by
2.842170943040401e-14 in LL,9.318379401435095e-13 in the largest score entry,
and5.820766091346741e-11 in the largest Hessian entry. Both have identical
settings and all integration-node/refinement paths to their32-node references,
and obey the count-axis chunk rule for all24rows. Mixed entries against the
independent finest Decimal reference are [-6.43641e-14,4.97992e-17] and
[1.73788e-13,-3.71047e-13], respectively, within the unchanged1e-5 criterion.
This is two-point numerical evidence, not full-eight-gene admission or
independent full-Hessian accuracy. It does not replace the running GL fits.

Observed single-job failed0gene2 elapsed times are320.0566265s at32 and
214.8816415s at64, both1CPU/4GiB. These are not a controlled hardware benchmark
or a general speed guarantee. Full-eight-gene checks remain required before
using the larger batch in new fits. The initial unsubmitted launch remains
preserved; executed launches are REC
`repair_p3_gamma_batch64_v2_launch_20261007T054337Z/`,
`repair_p3_gamma_batch64_positive_repair_launch_20261007T054942Z/`, and
`repair_p3_gamma_batch64_failed_repair_launch_20261007T055627Z/`.

The independent15-point accuracy collection is now4/15complete,11running,
0issues: REC `repair_p3_gamma_GL_five_gene_collection_20261007T060311Z/`.
Newly completed22082311 is gene1/25rows/tau.0125; it passes the native and GL
likelihood criteria, with GL-minus-scalar4.263256414560601e-14. Full-family
acceptance remains pending. Native collection at REC
`repair_p3_gamma_draw_collection_20261007T060312Z/` remains47/50complete,
3running,0issues. Both GL fit pilots remain running. B3 owner settings,
statistical calibration and qualified biological replication are still open.
No actual/sealed data, Resource claim or adoption changed. Documentation/syntax
checks pass; Resource-scope validation was not run for technical-only changes.

### Full-eight-gene64-node result and draw11 diagnostic,2026-10-07

Job22082978 completed0:0 (cpu/nslab,8CPUs/8GiB/90h,one individual job).
All8genes and the full summed likelihood,29-score vector and29x29 Hessian
are exactly numerically equal to the checked32-node result22078827.
Integration rules, node paths and logical node counts are unchanged.
All14 inherited boundary/error-mapping preflights passed. Per-row requested64
batches obey the unchanged100000-count-entry chunk limit; deep rows split.
Original14 numerical sources are unchanged. No fit or calibration was run.

Positive evaluation took150.14199898578227s versus211.06543478555977s for32.
Different hardware prevents treating this as a controlled speed estimate.
Total job elapsed2:41; sampled SLURM batch MaxRSS284316KiB. Parent RSS is
separately labeled; neither measurement is asserted to sum concurrent workers.
Executed inputs: REC `repair_p3_gamma_full64_launch_20261007T061453Z/`.
Independent sum, all8gene arrays, node counts, nonnegative mapped errors,
source hashes and exact14-preflight review:
REC `repair_p3_gamma_full64_root_review_20261007T061643Z/`.

Separate helper `repair_p3_gamma_batched_model64.py` changes requested batch32
to64, its own numerical-refusal filename, and adds a parent256MiB workspace
guard. Its source SHA is
`6f2d170724a73d29b09c6e3283c5edea1d9679dc23de766877dfd0baf8818a7b`.
The exact guard AST and pinned integer workspace function were checked on
admissible and oversized inputs without numerical imports; the real filename
refusal classifier and unrelated-origin control passed. This does not test an
actual64-node integration failure. Claude source review and primary responses:
REC `claude_gamma_full64_review_20261007T061402Z/`.

One separate diagnostic fit is running as22082995,draw11,8CPUs/8GiB/90h,
cpu/nslab,no array. Draw11 previously returned a nonconverged native null.
Its original supplied draw body, counts, starts, ties, convergence/failure rules
and observed null are preserved. It is not a replacement draw or a new
50-draw comparison. Original32-node draw0/draw18 fits remain untouched.
Runner `repair_p3_gamma_GL64_draw.py` SHA
`c021b8faed60529e2cab7e2f69e405a28296dab1e54442a77c6549518bab91b4`;
launch REC `repair_p3_gamma_GL64_draw_pilot_20261007T062020Z/`.
Static AST, original source pins, inspect and bash syntax checks passed;
startup reconstructed exact original counts and began fitting.
Separate collector `repair_p3_gamma_GL64_draw_collect.py` reports0/1complete,
1running,0issues; its terminal-fit checks are not run yet. It keeps this
diagnostic separate from all original draws. Claude review:
REC `claude_gamma_GL64_draw_review_20261007T061956Z/`. Its speculative
failure-branch concerns were resolved against the pinned original loop;
source-only review is not a numerical result.

Latest existing collectors report native47/50complete,3running,0issues;
GL32 pilots0/2complete,2running,0issues; independent likelihood family
6/15complete,9running,0issues. The six completed points cover all three
prespecified variances for genes0/1 and pass their likelihood criteria.
No all15 acceptance, converged GL pilot, independent full-Hessian accuracy,
statistical calibration or qualified human replication follows.
Missing owner B3 settings remain requested in the existing independent-data
request. No actual/sealed data, Resource figure, scientific authority or adopted
claim changed. Resource-scope validation was not run for technical-only edits.

### Corrected-Gamma adaptive quadrature hypothesis,2026-10-07

Two individual jobs22083040/22083041 (failed0gene2/positive18) were submitted,
each cpu/nslab,1CPU/4GiB/90h; no array. Both now retain FAILED1:0 from the
unchanged adaptive rounding-error refusal described below. The previous completed eight-gene64-node
comparison justifies retaining the batched GL calculation. This separate
hypothesis tests whether corrected Gamma3 arithmetic lets the unchanged
adaptive Gauss-Kronrod kernel return accurate results with fewer evaluations.
No speed or successful return is established yet. Existing native50draw,
GL32draw0/18,GL64draw11 and15-point accuracy jobs remain untouched.

New checker `scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_adaptive_p3_check.py`
SHA `7f35becb54430a08af5be31c4d2332d89ffea9982cb859789905a01b7770bafa`.
It reconstructs exactly the same173 synthetic counts and24 gene2 rows at
both already-frozen points. It uses original conditional mixture/selection
functions in an isolated module, with only the BB arithmetic replaced by the
already checked scalar Gamma3 function and restored in finally. No probability
law, optimizer, start, quadrature cap or fallback changes. All adaptive grid
intervals, mode breakpoints, domain/tail rules and refinement settings are
retained; only the empirical Hessian allowance is tightened to1e-5/8.
The field no_native_mode_search refers to the original GH approximation;
the adaptive normalized-prior kernel still locates modes for breakpoints.

Returned LL/29score/29x29 Hessian must match checked GL within
1e-9/1e-8/1e-5. The1.25e-6 allowance is an empirical per-gene error estimate,
not a claim of actual1.25e-6 or certified full-Hessian accuracy. Raw finite
returns are saved before comparison assertions. Known original adaptive
refusals retain their message, traceback and counted BB/conditional/first-row
calls; absent integration ladders remain explicit null. Neither a refusal
nor a scheduler failure may be counted as a passing comparison.

Executed launch REC `repair_p3_gamma_adaptive_launch_20261007T062706Z/`;
both pinned-input inspections, AST and bash syntax checks passed, and both
jobs wrote protocol/settings before numerical work. An earlier unsubmitted
preparation with a wrong assumed manifest filename is preserved separately.
Claude source review and primary response:
REC `claude_gamma_adaptive_review_20261007T062615Z/`.
Source inspection confirmed the direct original BB lookup and stdlib/third-party
imports; original fixed homozygote inclusion remains unchanged. From saved
GL arrays, max|H+ggT| is77904.0093 and10003.0068 for the two points,
below the scale at which the review estimated a precision-floor concern.
This does not establish adaptive convergence or certify subtraction error.

Both jobs have now terminated with the exact original-kernel refusal:
"Composite integral failed: Target precision could not be reached due to rounding error."
Elapsed3:30/3:03, sampled batch MaxRSS79136/80420KiB. They evaluated3410/3655
unique first-row nodes,81840/87720 conditional calls and204555/175470 Gamma BB
calls before refusing. No accepted integral, raw returned arrays, answer or
success summary exists. All sources, reconstruction/settings and failure
origin were independently verified at REC
`repair_p3_gamma_adaptive_root_review_20261007T063231Z/`.
This rejects the unchanged adaptive method at both diagnostic points, even
with corrected Gamma arithmetic. No cap or tolerance was relaxed; neither job
was restarted. Retain the checked batched GL calculation for continuing fits.
These numerical-method refusals are not added to the original50 P3 draws.

Owner B3 scenario/plant/coverage settings are still missing; the inspected
owner reply still ends at section11. No held comparison, actual/sealed outcome,
Resource scientific authority, figure or adopted claim changed.
Resource-scope validation was not run for these technical-only changes.


### First converged GL draw and supplementary score disposition,2026-10-07

Draw0 job22082078 completed0:0 after1:41:16,8CPUs/8GiB,cpu/nslab.
Its original supplied recipe returned converged null and marginal fits:
LL=-511.5784700491991,v=0,LR=0; all28 nuisance parameters exactly equal to
the null. Both warm/cold starts converged in both fits, with the original warm
tie retained. The largest independently rederived projected gradient is
4.189013858013624e-9 for both fits, below1e-6. The negative right variance
score is-541.6573914312444, correctly inactive at the lower boundary.
All26 positive-variance GL evaluations returned with eight genes before the
final zero-variance solution. Sampled SLURM batch MaxRSS302200KiB.

This replaces no original draw. The native GH method previously refused
draw0 at the saved nonconcave gene2 point; the separate GL diagnostic now
completes this draw's original starts/failure rules. Its173 synthetic counts
exactly match the original reconstruction, SHA
`5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691`.
Gene-label row keys were checked through the frozen gene_order; they must not
be compared directly with the dense integer-gene keys in the older receipt.
An initial root check assumed those representations were identical; no
numerical calculation failed. The corrected independent result is REC
`repair_p3_gamma_GL_draw0_root_review_20261007T064057Z/`.
Source/archived hashes, count reconstruction, LR/exceedance arithmetic,
bounds, projections, start ties and positive-evaluation count all passed.

The separately recomputed native40 boundary LL/score are exactly equal, but
v=0 has no quadrature-order comparison. Positive integration error allowances
are correctly non-applicable at the final boundary; no certified full-Hessian
or positive-fit accuracy claim follows. The existing collector now reports
one of two GL32 pilots complete, zero recipe failures and zero issues; the
other pilot remains running. The GL64 draw11 diagnostic has begun positive
evaluations. All50 agreement draws, the15-point likelihood family, genuine
calibration and qualified human replication remain incomplete.

A distinct prospective boundary-score hypothesis was reviewed with Claude:
REC `claude_boundary_score_hypothesis_20261007T063806Z/`.
For F_g(u)=exp(h_g(u)), Gaussian expansion gives
d_v log E[F_g(sqrt(v)Z)] at0 = (h_g''(0)+h_g'(0)^2)/2.
That derivative is already in the native source. A clipped upper-tail score
with full null refits could avoid positive-v fits in this supplementary test,
but supplies no tau estimate/bound, P2 equivalence or sealed posterior
projection. Estimated nuisances, active parameter bounds, selection and
unrepresented within-participant cross-gene dependence still require
independent calibration. Fixed-design parametric resampling is not itself
evidence of invalid donor inference; the dependence assumptions matter.
No score implementation, new family membership, RNG/tie tolerance or
resampling job was adopted. Original P3 remains the frozen LR procedure;
this hypothesis does not meet the requested complete-model objective.

No held comparison or actual/sealed outcome was read. Resource scientific
authorities and figures are unchanged. Syntax/technical documentation checks
pass; Resource-scope validation was not run for these technical-only edits.

### Adaptive GL order: two points and full eight genes passed,2026-10-07

Saved ladders at both P3 points and all8 observed-fixture genes already showed
16/32 agreement within unchanged quarter and combined error allowances.
This offline diagnostic is REC
`repair_p3_gamma_order32_saved_review_20261007T064641Z/`.
Potential repeated-node reuse was only13.33%/14.04% at the two P3 points,
so no node cache was implemented. Requested batch size64 is distinct from
quadrature order16/32/64.

A separate numerical helper now tries orders16/32. If tails fail, it expands
the domain using the existing formula/caps. If resolved tails but coarse
16/32 allowances fail the quarter budget, it computes64 and rechecks both
tails and the quarter budget. On the final finer grid it computes64 when
16/32 quarter or combined full-budget agreement fails. The original final64
branch requires the full combined budget, as before. Probability law,
row functions, moment/covariance formula, Kahan accumulation, grid/domain/tail
limits and precision budgets are unchanged. These are empirical precision
checks, not certified quadrature/rounding bounds.

New sources under scripts/analysis/in_biopsy_allelic/rec_codex:
- `repair_p3_local_gl_adaptive_order.py`, SHA
  `43d3cd87909d1fafde7d8a78b6e330fbdfa8a325ef32b8802f010e09c3a81d47`;
- `repair_p3_gamma_order_p3_check.py`, SHA
  `9fe45b5bf02d8ff39af58d1107a607ee1e033932aea7bdcf38e9a75012789f91`;
- `repair_p3_gamma_order_model.py`, SHA
  `04542c4685ab7cee72d0d1dd5120f18a60dc8570adc35c5b363776528afad3c2`;
- `repair_p3_gamma_order_full_check.py`, SHA
  `cbf797908cfad35f90cdcd0688742810712d734079b456d169ead851b2d19ab8`.

Two fixed-point jobs22083499/22083500 completed0:0, each cpu/nslab,
1CPU/4GiB/90h. Elapsed53s/36s, runner times50.1809s/32.6178s and sampled
batch MaxRSS62224/62296KiB. Both selected32 at both accepted stages.
Nodes23040/53760 and21888/51072, exactly57.1429% fewer.
Differences versus checked GL:
failed0gene2 LL6.82121e-13/score7.77522e-12/Hessian3.66078e-11;
positive18 LL6.53699e-13/score2.75602e-11/Hessian5.68434e-9.
Both actual Hessian differences also lie below the1.25e-6 per-gene budget.
Mixed entries minus independent finest Decimal references:
[-2.85929e-13,1.01841e-16] and[4.51093e-14,-1.28602e-13].
All executed-rung geometries/arrays, counts, source hashes and empirical
mapped errors were independently checked at REC
`repair_p3_gamma_order_pair_review_20261007T065322Z/`.
Launch REC `repair_p3_gamma_order_launch_20261007T065005Z/`.

Six exact decision-AST probes passed for no escalation, coarse+fine64,
fine-only64, combined-budget64, a still-refused64 and tail recheck after64:
REC `repair_p3_gamma_order_control_review_20261007T065148Z/`.
These used mock rule values and stdlib only, not scientific numerical tests.
Claude source review and primary response are in REC
`claude_gamma_order_review_20261007T064946Z/`.
Fine64 full-budget semantics were checked against the actual original source;
this does not relax an original quarter requirement.

Full8gene job22083751 completed0:0 in1:15, cpu/nslab,8CPUs/8GiB/90h.
All14 original boundary/error-mapping/inheritance checks and all8gene/full
comparisons passed. The positive evaluation took62.7362196s versus150.1420s
for the earlier fixed-order64-node-batch run and211.0654s for batch32.
Different hardware prevents a controlled timing ratio or full-fit speed claim.
Nodes172800 versus403200, again57.1429% fewer; all8 accepted orders are32.
Full differences versus checked reference are LL2.27374e-13,
score1.40744e-10 and Hessian9.42382e-8. Sampled batch MaxRSS283252KiB.
Independent sums, errors, inactive zeros and executed-rule subset checks:
REC `repair_p3_gamma_order_full_root_review_v2_20261007T065745Z/`.
The earlier root receipt's per-gene same_node_path label was corrected to
executed_rule_subset in v2; its numerical values did not change.
Launch REC `repair_p3_gamma_order_full_launch_20261007T065454Z/`.
No existing fit was restarted or switched. New-order full-fit agreement
remains to be tested before using it for a full draw family.

A runtime scenario calculation at REC
`repair_p3_gamma_runtime_scenarios_20261007T065044Z/` uses the single selected
draw0 cost, not a representative or calibrated forecast. If every resample
cost the same1:41 runner time on8allocated CPUs,50draws would reserve674.85
CPU-hours;1000draws13497.10;200x1000 P3 refits2.699million. Even ideal106
concurrent8CPU jobs would give132.65days for that last scenario. Actual costs,
full-scale dimensions and owner allocations are unknown; P1/P2 and observed
fits are excluded. No warp-speed cost was invented without its construction.
This motivates reducing numerical work rather than assuming parallelism alone
makes the prescribed calibration practical.

No original P3 statistic, B-FIX draw, actual/sealed outcome,117program
definition, Resource figure or adopted scientific claim changed.
Owner B3 settings, full15-point likelihood accuracy, all50 agreement draws,
genuine calibration and qualified human replication remain incomplete.
AST/inspect/bash/documentation checks passed; Resource-scope validation was
not run for these technical-only changes.

### Adaptive-order full-fit comparison submitted,2026-10-07

Job22083807 is running draw0 with the new adaptive-order GL, compared with
the already completed fixed-order GL32 result22082078. One individual job,
cpu/nslab,8CPUs/8GiB/90h; no array. Existing draw18/draw11 fits and native50
draws remain untouched. Startup verified the exact original counts and began
the original null/marginal fit. No final fit, comparison or speed result exists
yet. The old draw0 first positive evaluation began after about20minutes
including its first225.84s positive evaluation, so a null-fit phase without
positive progress events is not by itself a stalled process.

New runner `scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_order_draw.py`,
SHA `60ed9337dabe380d092b21780ba6db8ce4133bd3c971a340981f2eb237ad2465`;
launch REC `repair_p3_gamma_order_draw_pilot_20261007T070225Z/`,
pins SHA `f3d086646cabf0f69577d744b0ff49afcf43705a60f2b32c32fbafd5248997ab`.
All original16 sources remain unchanged; two checked numerical helper sources
are added. Full8gene proof and the completed reference draw hashes are required.
The original supplied draw body, observed null, starts, failures, restart/tie
rules and LR statistic are retained. Reference fit labels are only used after
the new fit, never as optimizer inputs. Progress records selected coarse/fine
orders for all8genes and actual quadrature-node counts.

The paired-fit comparison is saved before success summary and must pass:
LL/nullLL1e-6; raw gradient/null-gradient atol1e-8+rtol1e-6; kappa1e-4;
tau1e-3 if referenceLR>1, otherwise LR1e-6; independently projected null and
marginal scores strictly below1e-6. The raw paired-fit gradient condition is an
additional conservative diagnostic, distinct from fixed-point derivative
accuracy. If it fails at separately optimized points, preserve the returned
fits and diagnose at shared points rather than claiming scientific failure or
restarting blindly. All28 nuisance differences, actual variance and tau
difference are reported; they are not substituted for the stated criteria.

New collector `scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_gamma_order_draw_collect.py`,
SHA `94d0eac15526bb50056b7ddf81e5ef50a983acab1d97b7002c7c05ad69a9c0a4`.
It independently rederives the paired-fit arithmetic, fixed source/claim
identities, counts, LR/exceedances, projected scores and order-progress metadata.
Its live-input branch reports0/1complete,1running,0issues; terminal-fit
comparison checks are not run yet. Original native and earlier GL collections
remain separate. AST, inspection and bash checks passed. Claude source-only
review plus original-loop/primary responses:
REC `claude_gamma_order_draw_review_20261007T070210Z/`.
The actual original loop flags both null/marginal nonconvergence and uses
LR>=observedLR-1e-6; speculative strict-exceedance concerns do not apply.

A separate static resource inspection of the slow five-gene accuracy fixtures
is REC `repair_p3_gamma_G5_batch_feasibility_20261007T071234Z/`.
Those fixtures are G5/Z1/P16,17augmented coordinates,9local coordinates,
and depths20–36. The500-row point has14465 count-axis entries and a declared
GL workspace estimate9339648bytes at batch64, within unchanged100000-entry
and256MiB limits. This is a resource calculation, not numerical validation.
The G8/256-row context is not widened. Next useful bounded check is a separate
G5 context at the already completed25-row/tau.05 point against scalar
conditional and dense GL results before testing the full15-point family.
No G5 batch helper, new G5 job or scientific acceptance exists yet.

No calibration workload, held comparison or human/sealed data was read.
No Resource authority or figure changed. Full50 agreement, full15 likelihood
family, B3 settings/calibration and qualified human replication remain open.
Resource-scope validation was not run for these technical-only changes.

### Separate G5 batching: probe passed; fifteen-point family running (2026-10-07)

The new G5/Z1/P16 fixture context is separate from the existing G8 context.
Its scope is at most500 rows with depths20–36, unchanged100000 count-axis
entries and requested batch64. The generic base differs only in MAX_ROWS
256→500; immutable preparation, Gamma arithmetic and original likelihood
sources are unchanged. This is synthetic numerical work, not human inference.

The 25-row/tau.05 probe22083924 completed0:0 in23seconds on1CPU/4GiB.
All150 conditional calls (including repeated stage-zero predictor values)
passed against scalar Gamma: max differences LP0, gradient1.11022e-16,
Hessian6.93889e-18. The GL calculation took17.4508707seconds, used21888
quadrature nodes and547200 row/node evaluations, and selected32/32.
Against completed dense GL22082316, full17-score/17x17-Hessian differences
were LL4.26326e-14, gradient9.37916e-13, Hessian5.83896e-10. Independent
scalar LL difference was0. Native Gamma GH20/GH40 criteria also passed.
Root review rederived these comparisons, node coverage, exact inactive zeros,
local/full embedding and empirical allowances; environment, settings, counts
and selected point agree with the completed dense run.
Evidence: REC `repair_p3_gamma_G5_probe_root_review_20261007T072503Z/review.json`.

New files under scripts/analysis/in_biopsy_allelic/rec_codex:
- repair_p3_batched_conditional_g5.py:
  bf8c1354ea6f4b463370dbe102f080bf18635ebab60106f4adc39e802239948b
- repair_p3_gamma_batched_conditional_g5.py:
  94bdfae3654b16f6735c96d513945bf43075fe6325fcf2fe4a766fc22ef369f2
- repair_p3_gamma_G5_batch_check.py:
  15deb995d09c808416db7959730313e8b72d1addd230337b718e43aaffeb5f02
- repair_p3_gamma_G5_batch_family_check.py:
  64b130f6ba6d569c97dda84b063b7e2c0b0b2db7dbe34bc28a399261733c935c
- repair_p3_gamma_G5_batch_collect.py:
  886558824c921ccc2e1359073046f5d0a549981b8da5a9bfe545252be90f78d9

Probe launch REC `repair_p3_gamma_G5_probe_launch_20261007T072347Z/`,
pins d273ec741b2e14b3e77c1edba7fd6c8bd290602b0b693b6ac2f307c2907d7cbb.
After that completed proof, fifteen individual jobs22083938–22083952
were submitted, each cpu/nslab,1CPU/4GiB/90h, no array; total15CPU/60GiB.
Family launch REC `repair_p3_gamma_G5_batch_family_launch_20261007T072543Z/`,
pins ea630bba40d99d6aaf773a59d3b900f91fc8103fe2fff232edda4df7db507e36.
All15 static inspections and bash syntax checks passed. Each job pins the
completed probe and independently checked scalar family, reconstructs all805
original counts, compares every selected row at u=0,+tau,-tau in both Hessian
modes, retains Gamma GH20/GH40, and saves raw GL arrays before assertions.
The scalar references and original GH outputs are unchanged; GH20/GH40 are
recomputed with already checked scalar Gamma arithmetic.

At the latest collection,9/15 points completed and passed,6running,0issues;
no full15 acceptance is claimed. REC
`repair_p3_gamma_G5_batch_collection_20261007T073055Z/collection.json`.
The new stdlib collector rederives mass criteria, all17 dimensions,
9-to17 embedding, exact inactive zeros, empirical refinement-plus-tail
allowances, row/node coverage and support-batch counters. Gamma counters here
count support batches; they must not use the dense per-node counter rule.
Completed branches have executed on9 points. Claude source-only review is in
REC `claude_gamma_G5_collector_review_20261007T072803Z/`, with primary response.
The revised collector additionally recomputes first/final raw-array differences,
checks node coverage per row, reports support calls per batch, and exits nonzero
if a completed family fails a criterion. The speculative Settings concern was
resolved against the actual literal scalar defaults. These strengthened checks
passed on9 completed points; full15 remains pending.

Claude probe review REC `claude_gamma_G5_probe_review_20261007T072333Z/`
identified a fail-closed refusal-recording issue: the context wraps ValueError
as RuntimeError, so the runner's scalar ValueError whitelist cannot classify
that path. It propagates as a failed job, never a false success. Successful
outputs are unaffected; submitted sources remain unchanged. Primary response
records verified dense environment/settings/count equality and limits:
conditional probes do not cover all integration nodes;150 calls are not150
distinct predictor values; these shallow fixtures never exercise chunk
splitting; memory estimates are not measured peak. Dense/batched agreement
shares Gamma arithmetic and does not independently establish full derivative
accuracy. Error allowances are empirical, not certified.

Native50 is48/50complete with2running and0unexpected diagnostic issues.
Adaptive-order draw0 job22083807 remains running; by positive evaluation8
it used32/32 for all8 genes,173952 nodes and61.7260seconds for that evaluation.
This is not a completed paired-fit comparison or a controlled runtime ratio.
Earlier dense jobs and fits remain live and unchanged. No B3 calibration,
held comparison, restricted/human/sealed outcome read or paper adoption
occurred. Five authority documents and figures are unchanged.
Resource-scope validation was not run for technical-only changes.

### Draw11 retains the original convergence failure; analytic mass hypothesis checked (2026-10-07)

GL64 draw11 job22082995 completed0:0 in1:11:22, sampled batch MaxRSS313120KiB.
All25 positive-variance GL evaluations returned; no fitting integration failure,
caught exception or unknown exception occurred. The original prescribed
recipe nevertheless fails: null/marginal projected score
2.1271003368761657e-6 exceeds the strict1e-6 criterion. Returned v=0,
LL=-512.1899058535445. The null, final parameters, LL and29-score vector
are exactly equal to native Gamma draw11 job22081580; counts also match.
Native40 is identical at this boundary, where quadrature-order refinement
does not apply. The LR remains NaN and exceedance remains true under the
unchanged original failure rule. This is not an accepted fit or LR.

The saved null and marginal warm starts are nonconverged and cold starts
are marked converged. Cold-minus-warm LL gaps are-3.944933268940076e-11
and-1.3642420526593924e-11, respectively. The original rule keeps warm unless
cold is at least1e-6 higher; the boundary/null snapshot preserves the null
failure. No different start was substituted. Root rederived projection,
source pins, native equality, tie arithmetic and failure accounting:
REC `repair_p3_gamma_GL64_draw11_root_review_20261007T073350Z/`.
Collector REC `repair_p3_gamma_GL64_draw_collection_20261007T073940Z/`
completed1/1,1prescribed recipe failure,0unexpected issues.
The existing independent-data request document now asks the specification
owner to confirm this disposition or supply a separately frozen future
convergence-aware tie rule and its required calibration. No rule was changed.

A separate numerical acceleration hypothesis was reviewed with Claude in
REC `claude_analytic_initial_domain_review_20261007T073507Z/`.
For the selected BB genotype mixture, the numerator and inclusion denominator
eta scores lie in the BB score interval[-m*n,(1-m)*n], since mixture weights
and hom components are eta-independent. Their difference is bounded by n;
using2n is conservative. Predictors have eta(u)=eta(0)+stage*u with
dispersion fixed in u. Hence L=2*sum(n*abs(stage)) bounds the gene log
likelihood slope in u. With c=L*sqrt(v), integrating the pointwise bound
ell(sqrt(v)*z)>=ell(0)-c*abs(z) against the standard Gaussian gives
log mass >=ell(0)+log(erfcx(c/sqrt(2))).

This is only a proposed initializer: computed ell(0) and arithmetic are not
certified. It must fall back to the old initial domain if unusable and may not
introduce a new cap refusal. Actual LL/g/H tail and grid/order checks remain
mandatory. A different domain changes nodes and may change returned values
within tolerance; no bitwise identity or speed claim is made. The reviewed
r-score/Hessian moment extensions were not implemented or accepted.

New source:
`scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_analytic_mass_bound_check.py`,
SHA945dd46bd1082ae7a81126a4e3c4922278189d38d7cb165a0ad49c8bfc799660.
It reconstructs the exact805 synthetic counts before patching private BB
arithmetic, compares all15 scalar references, checks pointwise slopes at
u={-8tau,-tau,+tau,+8tau}, and tests120 eta-score/finite-difference cases:
depths20/36/3000, both integer selection-window edges, eta{-40,-10,0,10,40}
and four valid source-probability vectors including zero components.
No integrator, fit, threshold, resampling or human data changes.

Initial checker22084053 FAILED1:0 after20seconds because its adversarial
source vector(0,0,1) assigned all probability to the uncalled source category,
yielding zero weight for every modeled genotype. The conditional distribution
is undefined there. This was a checker-fixture error, not a bound counterexample
or scientific P3 failure. Its archived source and log are preserved. The
corrected checker explicitly requires nonzero genotype support, replaces that
invalid case with(.5,.5,0), and saves mass results before adversarial checks.

Corrected check22084074 COMPLETED0:0 in8seconds,1CPU/4GiB/90h,cpu/nslab.
All15 mass comparisons and120 score probes passed; no fit/integrator changed.
Launch REC `repair_p3_analytic_mass_bound_v2_launch_20261007T074245Z/`,
pins5cf41c9a39af8be0b708c25c54116d7bf61b149b610f33de710379bef7bb90aa.
Root independently evaluated erfc/asymptotic arithmetic, L and scalar gaps:
max log-bound difference2.842170943040401e-14; minimum gap below the
scalar reference after its error estimate0.9412676541234521;
max score finite-difference difference2.5132309122710694e-8.
REC `repair_p3_analytic_mass_bound_root_review_20261007T074453Z/review.json`.
This admits a bounded integrator comparison next, not a new production method.

Current G5 family is12/15completed/passed,3running,0unexpected issues:
all200-row points passed; only500-row points remain. Intermediate tail-envelope
overflow warnings occurred in larger cases; final finite arrays and unchanged
empirical precision checks still govern acceptance. No all15 claim exists yet.
Native50 remains48/50complete,27returned recipe successes and21failures,
2running,0unexpected diagnostic issues. Adaptive-order draw0 is still live
and had reached positive evaluation24, all8genes32/32,172800nodes,
41.4837seconds for that evaluation; no final comparison yet.
Live SLURM and fresh collectors, rather than these dated counts, control status.

No B3 calibration, held comparison or qualified human replication exists yet.
No Resource authority, figure, frozen program or scientific readiness changed.
Scope validation was not run for these technical-only changes.

### Adaptive-order draw0 full fit completed and matched (2026-10-07)

Job22083807 completed0:0 in42:48, sampled batch MaxRSS309128KiB.
The completed paired-fit collector reports1/1complete,0prescribed failures,
0issues. The final null fit is exactly equal to the completed fixed-order
GL32 draw0 reference22082078. Final LL=-511.5784700491991, all28 nuisance
parameters, all29 gradient entries, kappa, tau and LR exactly agree.
Returned v=0,LR=0; null/marginal projected score4.189013858013624e-9
passes the strict1e-6 rule. Both fits preserve original warm/cold convergence
and selection. All26 positive-variance evaluations returned before the final
boundary solution. Counts, all19 archived sources, pins, projection and paired
comparison were checked independently:
REC `repair_p3_gamma_order_draw0_root_review_20261007T074750Z/review.json`.

The actual run exercised adaptive order64: at evaluation14,
v=3.3445937817022095e-5, gene2 required refined order64 and gene5 required
coarse order64; the other206/208 gene evaluations used32/32. This is an
executed numerical fallback, beyond the earlier decision-mock checks.
Total saved positive quadrature nodes4570240. An initial root assumption that
every gene used32/32 was rejected against the progress records and corrected
before saving the root report; no returned numerical result changed.
The positive variance-evaluation sequence is not exactly equal to the older
fit's sequence, even though final parameters and statistics agree exactly.

Measured runner times are2564.5560seconds current versus6073.6945seconds
previous. These are separate hardware runs with both batch/order differences,
not a controlled hardware benchmark or a general runtime forecast. This is
one successful paired full fit, not fifty-draw agreement or calibration.
The separate draw11 remains a correctly retained original recipe failure.
A broader frozen synthetic comparison can now use the completed draw0 proof;
the original native50 and earlier draw18 diagnostic remain separate.

Latest G5 collection REC
`repair_p3_gamma_G5_batch_collection_20261007T074828Z/collection.json`
is12/15completed/passed,3running,0unexpected issues. The three500-row jobs
are confirmed live with sampled MaxRSS94260/94444/94736KiB. No all15
acceptance yet. Native50 still has2live jobs. Analytic mass-initializer
hypothesis has a completed compute/root check and source-assumption response
in `claude_analytic_initial_domain_review_20261007T073507Z/primary_response.json`;
no integration change has been implemented.

No B3 workload, protected outcomes, independent human validation or paper
adoption. Technical checks passed; Resource-scope validation not run because
authorities and figures remain unchanged.

### Fixed fifty adaptive-order draws submitted in parallel (2026-10-07)

The completed adaptive-order draw0 proof now supports a separate fixed50
synthetic P3 family. Jobs22084113–22084162 are50 individual jobs,
cpu/nslab,8CPUs/8GiB/90h each, no array. Total request400CPUs/400GiB.
Before submission, live usage was63CPUs;463 with these jobs is below the
live nslab850CPU cap. cpu max wall time was verified as14days and nslab
memory cap as7TB. No QOS/partition switch, cancellation or restart occurred.
All50 jobs are RUNNING; all50 have passed exact original-count checks and
written their claims/protocols. Largest sampled batch MaxRSS56396KiB during
the null-fit phase. No terminal new-family result exists yet.

Launch REC `repair_p3_gamma_order_family_launch_20261007T075115Z/`:
- runner repair_p3_gamma_order_family_draw.py SHA
  f27a3203641584f26637cfd8e2a05d532498159968cf0a4ef4386c2b1df40944;
- pins SHA a2bfd7657cbab245332603a6a0490ee4ae1dae99b5eb69c984a29f5d2dcdcc41;
- jobs.json SHA1606bfaa5ca6a1bae51dc38700790b695d15cadd4a55aaa99b6267318c0d4684.
All50 static inspections and shell syntax checks passed before submission.
Each receipt and job-index append was flushed/fsynced.

New source under scripts/analysis/in_biopsy_allelic/rec_codex:
`repair_p3_gamma_order_family_draw.py` broadens the passed draw0 runner to
the fixed draws0–49. All18 numerical/design helper sources are unchanged
from the completed fit proof. It pins that proof's summary, protocol, raw fit,
returned diagnostic, paired comparison and independent root review. Actual
converged fit vectors and helper hashes are checked directly. The original
draw AST,173 original rows,60 synthetic participants,50x173x2 supplied uniforms,
observed null, starts, ties, null snapshots, failures and LR/exceedance rule
are unchanged. Count equality is checked before any original draw/statistic.
This adds no replacement draws and does not mix native and GL results.

New source `repair_p3_gamma_order_family_collect.py` SHA
c87c10463ded5d1a013139cac14af98543833075a98a02dd06bcb82c49e6d207
collects this family separately. It checks all50 job identities, executed
archives, protocol/claim identities, saved counts, original failure/LR
arithmetic, projected scores, empirical allowances and native40 diagnostics.
The draw0 output is compared independently with the completed adaptive-order
proof using the prescribed LL/g/kappa/tau-or-LR and strict score criteria;
full fit-dictionary exact equality is additional reported evidence.
Completion requires every original recipe to finish and no unexpected issues;
all50 convergence is a separate explicit field. A known failed recipe is
never silently treated as a converged numerical fit.

The initial live collection reports0/50complete,50running,0issues:
REC `repair_p3_gamma_order_family_collection_20261007T075552Z/`.
Its new-family terminal branch is not executed yet. The pre-comparison
terminal checks were inherited unchanged from the completed GL64 collector;
the source equivalence record is
REC `repair_p3_gamma_order_family_static_review_20261007T075406Z/`.
Six saved-data control checks accept the actual matching proof and reject
LL, kappa, gradient, strict-score and failed-recipe mismatches:
REC `repair_p3_gamma_order_family_comparison_control_20261007T075552Z/`.
These are stdlib arithmetic/control checks, not numerical fits or calibration.

Claude source-only runner review and primary response:
REC `claude_gamma_order_family_review_20261007T075115Z/`.
No blocking changed-line issue was found. Pre-claim qualification failures
intentionally remain integrity/job failures, never statistical draws;
all50 actual prechecks passed. The proof/reference files are separately pinned
and actual fit vectors/helper hashes are checked, rather than trusting only
a root-review flag. The collector now supplies the requested draw0 output
comparison with the original tolerances. Unknown caught errors, unsuccessful
scheduler states or missing terminal summaries cannot establish completion.
The actual proof ends at v=0, so the boundary/native40 path is already tested.
The positive progress hook follows successful complete8-gene aggregation.

Latest G5 collection REC
`repair_p3_gamma_G5_batch_collection_20261007T075619Z/collection.json`
remains12/15complete and passed,3confirmed live500-row jobs,0unexpected issues.
Earlier dense jobs, native50 remaining draws and the older draw18 pilot
remain live and separate. No duplicate/replacement of their results was made.
No full50 independent agreement, B3 calibration, protected outcome read,
qualified human replication or Resource adoption is claimed.
Resource-scope validation was not run for these technical-only changes.

### Base generation checked; repeated draw0 exposes null-fit instability (2026-10-07)

The fixed50 adaptive-order family has a genuine repeat-agreement failure.
Draw0 job22084113 COMPLETED0:0 in15:31, but its original recipe fails.
Counts and18 helper hashes match the completed proof22083807. Returned v=0,
LL=-511.5784700491861 differs from the proof by1.3017142919125035e-11;
kappa differs by-3.951253724210346e-8, all28 parameters by at most
1.3513185057822952e-5. Full gradient difference is0.00016957607908807404.
The selected null/marginal warm starts are nonconverged; cold starts are
converged, with LL differences inside the original1e-6 warm tie. Projected
gradient1.1843350636815686e-5 fails the unchanged1e-6 criterion. LR remains
NaN and exceedance true, exactly as prescribed. No cold fit was substituted.
Root diagnosis:
REC `repair_p3_gamma_order_repeat0_root_review_20261007T082110Z/`.

Latest new-family collection
`repair_p3_gamma_order_family_collection_20261007T081929Z/`
has9/50completed recipes,3converged and6prescribed failures,41remaining,
and1agreement issue: draw0 versus its completed proof. Collector exit1 is
correct. This supersedes the startup-only0issue report; the family has not
passed repeat agreement or calibration. Existing outputs and live jobs remain
unchanged. No job was restarted or canceled because of elapsed time or failure.

The proof ran on Intel ne1dc3-007 and the new repeat on AMD ne1dc6-007.
To distinguish evaluator accuracy from optimizer-path variation, two separate
fixed-point jobs22084270/22084271 were pinned to those nodes,1CPU/4GiB/90h,
cpu/nslab. Both COMPLETED0:0 in4/2seconds. They evaluate both saved parameter
vectors using the original BB and scalar Gamma3, with no optimization.
New checker `scripts/analysis/in_biopsy_allelic/rec_codex/check_null_fixed_points_arch_v1.py`,
SHAb29045d92a988f4313ca8fd48d0554c2306a9a2bab6c34eef4c63032bc998e52.
Launch REC `null_fixed_points_arch_launch_20261007T081500Z/`,
pins3832d14f2d36fc068f51deaca0fc68ef528b3af2c3586296ca43346b5f26471f.

At the proof point, original projected gradients are4.189013858e-9 Intel and
4.199771475e-9 AMD; both converge. At the repeat point, both are about1.18433e-5
and fail. Across architectures, original fixed-point LL differences are at
most2.27374e-12, gradients7.01390e-10 and Hessians4.50697e-10.
Gamma3 cross-architecture gradients differ by at most2.54362e-13; replacing
original arithmetic by Gamma3 at either fixed point does not repair the
repeat's projected score. Root comparison:
REC `null_fixed_points_arch_review_20261007T081605Z/review.json`.
The CPUs are Xeon E5-2697v3 and EPYC9254; actual NumPy/runtime inventories
were saved on both nodes. Full-fit hardware causality is not established.
The different convergence outcomes concern the points reached by optimization,
not a large change in evaluating the same point.

V1.7 owner answers item8 leaves Newton-polish settings/scaling to each
implementation. A separate bounded polish trial can therefore retain the
scientific likelihood,1e-6 convergence cutoff, starts and warm/cold tie rule.
[Inference] Strict floating-LL-increase acceptance may stall near a flat optimum;
the fixed-point checks do not prove that mechanism. Any numerical repair needs
a new saved comparison and must not relabel these existing failed draws.

In parallel, a necessary piece of prospective calibration was implemented:
`scripts/analysis/in_biopsy_allelic/rec_codex/calibration_generate_base_v1.py`,
SHA176c4f18d85fa42951c0f59420be5b2299fb50b981267c8aae74f8a82e07d2ef,
and `check_calibration_generate_base_v1.py`,
SHAbc2e0b78a3c0ccca6252fe654524162be20e051cae24b555a6069eb19981eb91.
The generator draws one normal stage slope per gene before genotype classes
and unconditional BB counts, then applies the integer read-inclusion rule.
It never redraws to force inclusion and preserves empty datasets. Genotype
weights condition on the nonzero modeled support, matching the selected
likelihood's normalization; the uncalled/uncalled state is outside that support.
No B3 label-error law, plant amplitude, scenario allocation or fit is supplied.
Parameter truth is returned separately for the simulation evaluator.

Each dataset uses SeedSequence(seed,spawn_key=(replicate_index,)).spawn(3)
for effects/classes/counts. Explicit replicate IDs reproduce independently of
execution/job order. Constructor checks integer depth/minor/percent with
percent<=50. Current scope is at most8genes,7covariates,805template rows,
depth<=5000; this is synthetic scope, not the actual development design.
The scalar generator keeps only bounded zero-tau count CDFs; no persistent
positive-tau node/row cache is introduced.

Initial job22084238 passed a weaker checker; Claude correctly found that its
positive-tau inclusion moment alone could not distinguish a disconnected
random effect. Those initial source/output snapshots are preserved and are
not treated as effect-path validation. Review and primary response:
REC `claude_calibration_base_generator_review_20261007T080726Z/`.
The strengthened checker covers stage signs-1/+1, two genes, orientation-1,
nonzero omega_S/dispersion slopes, genotype-specific preselection count means,
negative-orientation class probabilities and cross-gene effect independence.

Strengthened job22084262 COMPLETED0:0 in7seconds (runner6.33559seconds),
1CPU/4GiB/90h,cpu/nslab; AST, inspection and bash checks passed.
Launch REC `calibration_generate_base_v2_launch_20261007T081313Z/`,
pinsdf702a67159e66f43b967c5484b67100dd861019225c99b648e8d0fc7dfd6139.
It checked20,000zero-tau datasets,6,000positive-tau datasets, seeded/indexed
replay and10empty datasets. Of the20,000zero-tau rows for the declared probe,
7178 passed read inclusion. Independent original-likelihood selected counts
and scalar quadrature of inclusion/u-times-inclusion agree within declared
Monte Carlo checks. Three signed heterozygote preselection effect/count
covariances are separated from zero by14.1675/12.6847/14.0983 measured MC SE,
exceeding the prespecified6SE discrimination check; their signed expectations
and empirical values agree. Genotype-specific means distinguish swapped hom
classes. Root rederived source hashes, recorded arithmetic and comparison
limits:
REC `calibration_generate_base_root_review_20261007T081745Z/review.json`.
These are generative-law checks, not fitted-model rejection-rate calibration,
power, coverage, biological replication or a B3 dataset allocation.

The existing owner request now reflects the successful numerical proof,
failed repeat and checked base generator. Owner B3 amplitudes/allocations,
coverage target, C5 summaries and F558/F237 identity/eligibility answers remain
absent; no values were invented. Actual development G/Z/depth/resource bounds
also remain unknown. No protected molecular/clinical/sealed outcomes were read.
Latest G5 collection081929 is12/15passed with3live500-row jobs and0issues;
native50 is48/50complete with2live jobs and0unexpected diagnostic issues.
No full50 independent agreement, statistical calibration, qualified human
replication or paper adoption. Authorities/figures unchanged; scope validation
not run for technical-only changes.

Subsequent live update: G5 collection
`repair_p3_gamma_G5_batch_collection_20261007T082239Z/collection.json`
checks14/15 completed and passed,0issues. Jobs22083942/22083947 completed0:0
in54:49/54:41 with sampled MaxRSS94292/94508KiB; these are the500-row
tau.0125/.05 points. Job22083952 (500rows,tau.2) is still confirmed live.
This does not resolve the separate full-fit draw0 convergence disagreement.

### Score-integral null polish checked; separate full-fit trials running

The separate saved-point test repaired both known null convergence failures
on Intel and AMD without changing the likelihood, parameter bounds or 1e-6
scientific cutoff. The already converged control was unchanged. This is a
numerical diagnostic, not calibration, fifty-draw agreement or human validation.

Final checker jobs22084403/22084404 completed0:0 in1:17/0:27, each
1CPU/4GiB/90h,cpu/nslab. Intel sampled MaxRSS52172KiB; AMD's reported zero
does not measure its actual peak. Each checked proof0, repeat0 and failed11.
Repeat0 projected gradients fell from1.18433e-5 to1.51230e-9/1.20549e-9;
failed11 fell from2.12710e-6 to2.04845e-9/2.05890e-9 (Intel/AMD).
Rerunning the original constrained polish alone left these above1e-6.
Each failed point needed one accepted new step, the control zero.

The proposal uses a regularized Newton direction and the original bounds.
It accepts a step only if both original/Gamma projected gradients decrease
by at least half the damping fraction and both 16/32 directional-score
integral gain estimates exceed32times their observed arithmetic difference.
The checker separately recomputes accepted increments with original BB scores.
Positive gains are about9e-12 while directly subtracted LL differences can
be negative at rounding scale. Across CPUs, final nuisance parameters differ
by at most7.06902e-12; LL differences are at most6.20731e-11.
All three final boundary variance scores remain negative.
The gain checks are empirical, not certified error bounds; shared likelihood,
mixture and predictor errors are not excluded.

Root review REC
`score_integral_polish_v3_root_review_20261007T084056Z/review.json`.
Final check launch `score_integral_polish_v3_launch_20261007T083901Z/`.
Earlier preliminary check jobs22084378/79 and22084397/98 and their snapshots
remain separate. Claude source reviews:
`claude_score_integral_polish_review_20261007T083314Z/` and
`claude_full_polish_trial_review_20261007T083919Z/`.
The first review prompted the sufficient-decrease and arithmetic-noise checks.
The second found no definite blocking source error; its primary response
records pool lifecycle/function checks, actual scheduler node identity and
remaining numerical limits. Sources were reviewed, not executed by Claude.

New sources under `scripts/analysis/in_biopsy_allelic/rec_codex/`:
- `repair_p3_score_integral_polish.py`,
  SHA9920b8bfb7a2e91ffd448e045086c4c2beee3aedcfc234950905bb300c8c1f0a.
- `check_score_integral_polish_v1.py`,
  SHA5bafa25a61b920b0395d5b30e4e4d984a16ecc8922c2049014a62d7ac7fa3f4d.
- `repair_p3_gamma_order_polish_model.py`,
  SHA6a661d41ac8faa3d1672554af8b2612c056f0a7af51d7b808c9a453a9f0c22b9.
- `repair_p3_gamma_order_polish_draw.py`,
  SHA1a9f1e307d4fff6c2e28f3148b8d531e25d2586e63444e8abca415901a80e38c.

Only constrained_polish is overridden in the separate full-fit model.
Original starts, warm/cold tie, count generation, P3 draw body, variance
boundary handling, integration and scientific cutoff remain inherited.
V1.7 item8 permits implementation-specific Newton-polish settings.
No existing fit or failure is relabeled and no start is selected differently
outside the original tie rule.

Three individual8CPU/8GiB/90h,cpu/nslab trials were submitted:
22084419 draw0 on ne1dc3-007,22084420 draw0 on ne1dc6-007,
22084421 draw11 on ne1dc4-045. All are running and passed the actual
original-count prechecks. Launch REC
`repair_p3_gamma_order_polish_launch_20261007T084132Z/`,
pins66e8d8f261591314ccda9de7384c842bc3fec0de460a6688ccc76383786266f0.
Actual scheduler nodes and same-node earlier compute CPU inventories are
saved in `scheduler_node_identity.json`; there was no fresh CPU probe in
these three runners. Full-fit reproducibility remains unverified.

The existing unchanged50-draw collection at084255 has37terminal recipes,
30converged and7prescribed failures,13running,0pending, and the known draw0
agreement issue. Collector exits1 for that issue. G5 is14/15passed,0issues,
with22083952 still running. Native50 remains48/50terminal,2live,0unexpected
issues at the083420 collection. The new trials do not replace these results.
B3 scenario/plant settings, coverage target, C5 information and qualified
development/replication identity inputs remain absent. No protected data read,
scientific authority/figure change or paper adoption; scope validation not run.

Subsequent collection084607:42/50terminal,33converged,9prescribed failures,
8running,0pending; the known draw0 agreement issue remains. New full trial
22084420 has actually exercised the wrapper: at177.285seconds its null polish
returned original projected gradient1.20549e-9 and Gamma8.94068e-12.
This is an intermediate null-polish result, not full-fit convergence.
All three new trials remain running; sampled batch MaxRSS is about55MiB each
before positive integration. G5's last job sampled MaxRSS94880KiB.

### G5 complete; full-fit polish still fails a positive integral

The completed G5 family now has15/15 passing selected-likelihood checks,
0issues: REC
`repair_p3_gamma_G5_batch_collection_20261007T084735Z/collection.json`.
Root review:
`repair_p3_gamma_G5_batch_root_review_20261007T084823Z/review.json`.
Across5/25/75/200/500rows at tau.0125/.05/.2, GL versus independent scalar
LL differs by at most9.66339e-13; native GH20/40 differs by at most9.09495e-13.
The last job22083952 completed0:0 in1:20:24,MaxRSS94880KiB.
Its actual GL time was4733.842seconds. The500-row tau.0125/.05 times were
3214.393/3219.695seconds, not14/19seconds. These are fixed synthetic integrals,
not fit/calibration checks, and full derivatives remain without an independent
reference. The G8 context was not widened.

Separate full trial22084420 completed the recipe in10:36,exit0, but the
prescribed P3 draw failed. The null polish succeeded; positive call14 then
refused in genes2/5 at v=3.344593757321353e-5 with
“GL order refinement exceeds unchanged precision allowance.”
The returned fit is absent, LR is NaN and exceedance remains true under the
unchanged recipe. Sampled MaxRSS294264KiB.
Collector `repair_p3_gamma_order_polish_collect.py` records the scientific
failure separately from scheduler failure, saves the numerical reason and
batch memory, and compares the two draw0 trials with the prior converged proof
under the prescribed tolerances. Only its live/failed-recipe branches have
executed; successful terminal comparison remains untested.
Latest collection:
`repair_p3_gamma_order_polish_collection_20261007T090839Z/`.
Intel draw0=22084419 and draw11=22084421 remain confirmed live.
No full-fit convergence/agreement claim.

The saved refused initial64 variance-score allowances were2.64834e-9 and
4.50965e-9, above the initial quarter budget2.5e-9. A separate helper allows
initial64 to proceed under the full budget to the unchanged final expanded
domain/halved-grid/order combined test. Initial32 retains the quarter check.
Both fixed gene checks exercised initial64 and returned within that final
empirical test. Gene5 job22084481 completed0:0 in1:46; its native40/80 score
differences are8.58589e-9/4.25720e-9. Gene2 job22084480 exited1 after saving
a completed numerical return: native40 differs by4.62231e-8, native80 by
5.61841e-9. Its LL differences are at most7.10543e-13.
This is unresolved variance-score reference precision, not qualification of
the new helper. No Hessian/full-fit use or adoption.
Evidence:
`GL_initial64_refusal_launch_20261007T085933Z/` and
`GL_initial64_refusal_root_review_20261007T090950Z/review.json`.

[Inference] The variance score's1/v amplification can turn tiny conditional
LL rounding into the observed order disagreements. The16/32 disagreement
does not shrink by orders of magnitude at64. The next independent numerical
check should compare the Gaussian integration-by-parts identity
score_v = .5 E_posterior[ell_uu + ell_u^2] with the saved normal-moment score,
using scalar integration before changing any fitted model. This is an
equivalent derivative hypothesis, not a new estimand or calibrated model.

In parallel, bounded profiling found BB arithmetic dominated selected
conditional cost. Initial job22084454 failed before profiling because saved
counts were applied in gene-block order rather than Model's individual/gene
sort. Corrected22084463 completed0:0 in7seconds. The original fixture creator
was inspected to establish its stored count order.
Paired22084482 completed0:0 in16seconds on ne1dc3-012:
500rows x3grids x64nodes =96,000 conditional values per backend.
Existing count-prefix and Gamma outputs differ by at most
1.36780e-13 LL,1.71419e-13 gradient,3.74812e-13 Hessian.
The profiled passes took8.85580seconds Gamma and3.47573seconds prefix.
These one-pass timings include profiler and unequal small bookkeeping costs;
their2.54790 ratio is descriptive, not an end-to-end speed estimate.
Raw paired arrays are retained in
`profile_G5_conditional_22084482/outputs.npz`; root metadata inspection:
`G5_prefix_pair_root_review_20261007T090240Z/review.json`.
No independent full-integral/derivative correctness follows from those arrays.

A separate500-row/tau.2 prefix integral now runs as22084526,
1CPU/4GiB/90h,cpu/nslab, using the original adaptive-order integration source
and exact seeded fixture/counts. It will compare full17-coordinate arrays with
the completed identical Gamma point and retain the scalar LL checks.
Launch `repair_p3_prefix_G5_fixed_v2_launch_20261007T090839Z/`.
Initial launcher22084512 failed before numerical evaluation from an extra
load argument; corrected source passed AST/load-signature, metadata and bash
checks. No current G8/source/fit was modified.

New sources under `scripts/analysis/in_biopsy_allelic/rec_codex/`:
- `repair_p3_gamma_order_polish_collect.py`,
  SHAbc652d0fe233b2eab912ca290dfa871fa48880a7f849e39e1cc66202b4acc788.
- `profile_G5_conditional_v1.py`,
  SHAacdff342822730803d8213edc6af4cce88f683ab72a5db2b04fc250e58482967.
- `repair_p3_local_gl_initial64.py`,
  SHAcea90f5189767b99de7171d8506cad2f583e4cd6c7eac533b8f31cf6b8bbde29.
- `check_GL_initial64_refusal_v1.py`,
  SHAa0bcae3a789e3be6c26c6d3c676f53d60bc5329e5f60c9daeba3d165c7a2c1ac.
- `repair_p3_prefix_conditional_g5.py`,
  SHAf0c5784f05bfa2168264da64dc1c5b32768c30adad5a2754bf538032af556839.
- `repair_p3_prefix_G5_fixed_check.py`,
  SHAcbe7c9c87b4cb9dda16f193b5cd456879c8ad340b1c818568643451ea75296b7.

Claude review and source-grounded primary response:
`claude_initial64_and_comparisons_20261007T085957Z/`.
The original50 collection090839 has43terminal recipes,34converged,
9prescribed failures,7running,0pending and the known draw0 agreement issue.
Native50 remains48terminal/2live/0unexpected issues at090839.
B3 scenarios/plants/coverage/C5 and qualified human inputs remain missing.
No sealed/protected data, calibration, authority/figure change or paper
adoption; scope validation not run for these technical-only changes.

### Stable variance-score references and the original constrained solver

The two heat-identity checks completed0:0:22084571 in1:00 and22084572 in1:06,
sampled MaxRSS54608/54520KiB. They evaluate the same saved genes2/5 and
v=3.344593757321353e-5 using
score_v = .5 E_posterior[ell_uu + ell_u^2].
This follows Gaussian integration by parts, with the original selected
RNA-count likelihood and unchanged stage/u units. Each backend used scalar
quadrature on standardized domains24/26, epsabs/epsrel1e-12, and conservative
native-BB derivative tail envelopes. Error estimates are empirical.

Gene2 Gamma heat score is603.5539492528011; count-prefix is603.5539492544757.
Their difference is1.67461e-9; the largest prefix domain change is1.90358e-9.
Gamma's estimated score allowance is1.32350e-10, prefix2.27046e-10.
Native40 differs from Gamma heat by-4.90043e-8; native80 by-8.39953e-9.
The new GL result differs by-2.78113e-9 from Gamma heat and-4.45573e-9
from prefix heat. Gene5 Gamma/prefix scores are-8.483609169389252 and
-8.483609169360351; their GL differences are about-6.1/-6.4e-10.
All declared domain/backend checks passed without relaxing a threshold.
This resolves the two diagnostic reference comparisons more directly than
the noisy Gaussian second-moment subtraction. It does not establish uniform
variance-score precision, mixed/variance Hessians, full fits or calibration.

Root review independently rederived score/allowance arithmetic and checked
Gaussian tails via stdlib erfc:
REC `heat_variance_score_root_review_20261007T092923Z/review.json`.
Launch `heat_variance_score_v2_launch_20261007T091815Z/`.
Initial22084559/60 failed at import from an omitted tests_v1 dependency;
the complete local import closure was checked before the corrected jobs.
Claude independently rederived the identity, chain rule, units, envelopes,
tail and ratio-allowance algebra:
`claude_heat_variance_score_review_20261007T091619Z/`.
The backends share priors, predictors, homozygous terms and selection/mixture
logic; they are independent numerical derivative/arithmetic representations,
not independent whole-model validation.

After the numerical runs, source guards were added for the28 nuisance
coordinates,29 gradients, exact gene-local embedding, count hash and reviewed
GL d/dv source hashes; original BB/import/context identity checks were added.
The actual two saved inputs passed the new metadata inspection. Numerical
kernels are unchanged. These additional runtime identity checks were not
rerun. Executed checker SHAba0c8c99ee2a789749f213b78ad8529c82bd95771d05ff3a4058f59a0cb957d4;
current `check_heat_variance_score_v1.py`
SHA01c5fb298224e3e79fd7c90bbbcdbad84c35479aaeab35eaca5e045cf59a7623.

All three first-wrapper full trials are now terminal. Intel0=22084419
completed0:0 in49:02,MaxRSS310568KiB; final LL/parameters/gradient/null are
exactly equal to the prior converged draw0 proof. It never called the new
polish. Draw11=22084421 completed0:0 in46:40,MaxRSS312300KiB; selected null
and marginal converge, v=LR=0 and final projected gradient9.99438e-11.
Both marginal starts converged. Its null warm start remained nonconverged,
with LL-512.1937008376393; null cold converged at-512.1899058535839 and was
selected under the unchanged rule. AMD0=22084420 remains the known positive
integration failure. Thus2/3 selected full fits converge, not3/3 agreement.
Root review:
`repair_p3_gamma_order_polish_full_root_review_20261007T094110Z/review.json`;
collection `repair_p3_gamma_order_polish_collection_20261007T093840Z/`.

The direct clipped-Newton replacement made less progress from the draw11
warm point than the original constrained solver. It should not replace that
solver. A separate wrapper now calls the original constrained solver first,
then score-integral refinement. Saved warm replay22084645 completed0:0 in22s
(the solver/refinement took11.29977s). Original projected gradient fell
.00188503 ->1.18369e-6 under the original solver ->2.04627e-9 after one
accepted refinement step. Final LL is within4.68390e-11 of the earlier
near-optimal warm value. The scientific cutoff and tie remain unchanged.
Root review:
`original_first_warm_polish_root_review_20261007T093723Z/review.json`;
launch `original_first_warm_polish_launch_20261007T093355Z/`.
This is a saved warm-point replay, not a new full-fit qualification.

A separate Gamma-only null-fit check retains the original fitter, starts and
tie, then re-evaluates in original arithmetic. Draw0=22084585 completed0:0
in11:03 (fit662.93824s), both starts converged and warm kept. Original
endpoint projected gradient2.74362e-8, Gamma/original LL difference-5.11591e-13,
maximum gradient difference1.80648e-9. This demonstrates a valid endpoint,
not a speed gain. Draw11=22084590 remains live.
Launch `Gamma_null_refit_launch_20261007T092056Z/`;
root draw0 review `Gamma_null_refit_root_review_20261007T093805Z/review.json`.

New additional sources under `scripts/analysis/in_biopsy_allelic/rec_codex/`:
- `check_Gamma_null_refit_v1.py`,
  SHAac4ac10d8a204cfdbc66b6ebc2ffbf126bc9073e862a555c3edc1e5a1456f467.
- `repair_p3_gamma_order_polish_model_v2.py`,
  SHA75d7b5bc4d4ef3cbf31482151338deb3aabf0030fa538d9db04663ae6f83b45e.
- `check_original_first_warm_polish_v1.py`,
  SHA82681df79a67bd4c68cb1251138dbad0cad9b66ade4c77231b3ca6907033d8e2.

The next variance-derivative candidate should use u=sqrt(v)z under a fixed
standard-normal measure, deriving both score and Hessian by that chain rule.
It must be checked against these heat references and the existing full
nuisance/mixed Hessian cases. Its derivative-tail envelopes must be derived
for that representation; the old Gaussian-moment tail bounds cannot simply
be reused. Retain the original constrained solver before refinement.
No positive full fit or fifty-draw run using either new candidate is admitted
by these fixed-point checks alone.

The unchanged original50 collection093840 is46terminal,36converged,
10prescribed failures,4running,0pending, with the known draw0 disagreement.
G5 prefix integral22084526 was live at that collection. B3 scenario/plant/coverage,
C5 and qualified human inputs remain absent; no sealed/protected data read,
calibration, paper adoption or scientific authority/figure change. Scope
validation not run for these technical-only changes.

### Completed prefix integral and pathwise derivative check, 2026-10-07

Job22084526 completed0:0 in37:36, MaxRSS89548KiB. At the exact saved
G5/gene4/500-row/tau.2 point, prefix GL took2176.27863s. Its scalar-reference
LL difference is2.04636e-12. Full-array differences against the identical
Gamma point are LL1.59162e-12, gradient2.81114e-12 and Hessian1.19094e-10.
Root independently checked the saved arrays, exact inactive zeros, pins,
executed sources and coverage of all59904 nodes for each of500 rows.
Maximum empirical mass/gradient/Hessian allowances are1.97898e-64,
6.65140e-12 and8.02984e-11. Review under REC:
`prefix_G5_integral_root_review_20261007T095240Z/review.json`.
The earlier Gamma integral took4733.84208s on a different node. This is not
a controlled speed ratio, qualification of all15 prefix points, independent
full-derivative validation, a full fit or human validation.

Gamma-only null draw11=22084590 failed1:0 after31:58, MaxRSS51904KiB.
Both starts ran; warm remained nonconverged and cold converged. The unchanged
LL tie selected warm. Original endpoint projected gradient1.04983169e-6
exceeds the strict1e-6 cutoff. Preserve the failure; do not substitute cold
or relabel earlier results. Root:
`Gamma_null_refit11_root_review_20261007T095643Z/review.json`.
Gamma-only arithmetic therefore does not resolve this convergence problem.

New diagnostic sources are `repair_p3_pathwise_moments_v1.py` and
`check_pathwise_moments_v1.py`. They differentiate under u=sqrt(v)z, including
mixed/variance Hessians, and derive separate selected-count derivative tails
using normal moments of1, |z| and z^2. They do not enter a fitter.
Initial jobs22084740/41 refused quadrature precision. Their source snapshots
and errors remain preserved. A separate retry scales derivative components
using z=0 and both one-SD nodes, requests quadrature1e-12, and saves numerical
refusal diagnostics. Scientific LL/gradient/Hessian acceptance remains
1e-9/1e-8/1e-5. Full-Hessian domain agreement would not establish independent
Hessian correctness. Claude source-only review could not run because its CLI
reported a session limit; request and response are retained under REC
`claude_pathwise_moments_review_20261007T095533Z/`.

The unchanged50 collection `repair_p3_gamma_order_family_collection_20261007T095400Z/`
still records46 terminal,36 converged,10 prescribed failures and4 running,
with one known draw0 disagreement. No calibration or human evidence is added.

The scaled pathwise retry is
`pathwise_moments_scaled_launch_20261007T095700Z/`, jobs22084749(gene2) and
22084750(gene5), each1CPU/4GiB/90h cpu/nslab. Gene5 completed its numerical
comparisons but correctly exited1:0 after52s, MaxRSS53888KiB: variance score
agrees with the Gamma/prefix heat references within5.88006e-10/6.16907e-10,
while variance-variance Hessian changes2.06185e-5 between domains24/26,
exceeding1e-5 even before adding empirical allowances. Gradient domain change
is1.33878e-9. The score result is useful; the joint derivative check fails.
Root rederived the moment ratios, covariance subtraction and propagated
allowances from the saved arrays:
`pathwise_gene5_root_review_20261007T095924Z/review.json`.
Gene2 remains live at this update. No full fits using this representation.
Current moment source SHA313af3fba14b71fe208ef69329a66caf48df694d3749b1c1c0a70467e9ffc7c7;
scaled checker SHAfa76f7a2a6f19ba5222386593dfad5504fc702061251dc948e083057ffcf82ba.
Initial checker SHA97b7ee9a8cf1e47db215c70369fa5fb69fa2a7db5f949a81e3b30f134ceffcb0.
AST checks, launcher `bash -n`, both pinned metadata inspections and the
three-document `git diff --check` passed. Scope validation not run because
these changes only record synthetic diagnostics in technical documents.

### Extended precision and independent variance curvature, 2026-10-07

The first pathwise gene5 failure is in coordinate(28,28), the variance
curvature. A separate native selected-count conditional helper now evaluates
integer rising products, oriented genotype mixtures and direct inclusion sums
in long double. The probability model, integer support and assay remain the
same. Its scope is one saved G8/Z7/P28 synthetic gene, at most32 rows and
100,000 count-axis entries; it is not installed in a fitter.
Source `repair_p3_longdouble_conditional_v1.py`,
SHA66c269cf1679ded3461564537429988925670a05ff1beb7c3b94e8625009031d.

Frozen checker `check_pathwise_longdouble_v1.py`,
SHA81e45a6b1fae6b48b6e17e37fd80e34a9077feb2a5b007f92b378cf594bfa15d,
ran in `pathwise_longdouble_launch_20261007T100452Z/`.
Gene5=22084793 completed0:0 in1:03, MaxRSS59744KiB. Domain24/26 differences
are LL0, gradient2.95231e-12 and Hessian4.32993e-8; all propagated domain
allowances pass unchanged thresholds. Variance-score differences from the
two heat references are7.45360e-12/3.63550e-11. All25 rows at three nodes
agree with Gamma within LL3.76588e-12, score2.24266e-13, curvature1.49117e-13.
One deepest row(n2994) was also checked against80-digit, probability-only
Decimal finite differences at step1e-10, with maximum score/curvature
differences about1.1e-14/1.2e-14. There is no opposite orientation in this
gene, so this is one independent row check, not two. Raw FD values were not
retained; root checked the reported comparison and rederived integral ratios,
covariance and allowance arithmetic, not the raw FD differences:
`pathwise_longdouble_gene5_root_review_20261007T100747Z/review.json`.

Gene2=22084806 completed its comparisons but correctly failed1:0 in42s,
MaxRSS58764KiB. Domain differences are LL0, gradient1.47793e-12 and
Hessian4.61005e-8, but its two Hessian allowances8.77183e-6/8.26641e-6
sum above1e-5. Both heat comparisons pass. One deepest row(n2466) agrees
with Decimal FD curvature within1.02319e-12. The separate tighter-request
check22084835 still fails the combined allowance, despite stable arrays;
do not call either a joint pass. The original Gamma check22084749 remains
live and has not been restarted or cancelled.

An independent whole-line variance reference uses80-digit probabilities and
Gaussian density derivatives, with no candidate conditional derivatives:
I'/I=E[(z^2-1)/(2v)] and I''/I=E[(z^4-6z^2+3)/(4v^2)].
Decimal-refined probabilists' Hermite roots/normalized weights at fixed
orders32/64/128 are checked against their polynomial moments. All224 nodes
are evaluated with eight workers; covariance is formed after summing complete
raw moments. Checker `check_decimal_variance_curvature_v1.py`,
SHA456c155be7c74510c173605afab51c4d4940a40e3e8f5cf9039025048e48ef8d.
Jobs22084836(gene2) and22084837(gene5) completed0:0 in2:00/1:43,
MaxRSS101700/108880KiB. Candidate-minus-reference magnitude is5.23245e-8
for gene2 variance curvature and1.29128e-8 for gene5; score differences
are1.92441e-12/9.89991e-13. Both pass unchanged LL/score/curvature thresholds
1e-9/1e-8/1e-5. Order64/128 differences are below4.5e-67 in these saved
curvatures; this is empirical order agreement, not a certified error bound.
Root independently summed all node moments, rederived LL/score/curvature,
checked normalized Gaussian moments0–6, and verified pins/archived sources:
`decimal_variance_curvature_root_review_20261007T101621Z/review.json`.
This independently checks the variance coordinate at two saved gene points,
not the complete nuisance/mixed Hessian, an optimizer, calibration or humans.

The variance second numerator is now being integrated separately in its own
units to avoid multiplying a common vector max-norm error by its large scale.
First diagnostic22084861 refused its overly strict scalar quadrature request;
the saved scalar estimated error is3.18581e-6, versus8.24099e-6 from the
vector error. A separate retry changes the scalar requested absolute error
from1e-8 to1e-5 while retaining the final combined LL/g/H acceptance
1e-9/1e-8/1e-5. Passing that final combined check remains required. No
existing result is reclassified. Launches and complete failure diagnostics
are retained in `pathwise_component_validation_launch_20261007T101529Z/`
and the subsequent `pathwise_component_budget_launch_` directory.

Original50 collection `repair_p3_gamma_order_family_collection_20261007T100600Z/`
is47 terminal,36 converged,11 prescribed failures,3 running and0 pending,
with the known draw0 disagreement. Draw49 retained selected null/marginal
nonconvergence and NaN LR. B3 settings, C5, qualified human inputs and final
model validation remain outstanding. No paper/figure authority changed.

The separate-numerator retry22084870 completed0:0 in1:11,
MaxRSS58840KiB and passes the unchanged combined limits. Gene2 domain
Hessian difference is5.54137e-8; the maximum combined Hessian allowance is
7.66527e-6, below1e-5. Variance curvature differs from the independent
Decimal reference by5.55842e-8. Root verified the scalar numerator replacement,
per-component errors, covariance propagation, input pins and archived sources:
`pathwise_component_budget_root_review_20261007T101933Z/review.json`.
Thus both saved gene points have passing pathwise domain checks and independent
variance-coordinate checks. This does not qualify a complete full-fit backend.
Successful checker `check_pathwise_longdouble_component_budget_v1.py`,
SHAe3aec3a0a92bfbb483ceee75249b56704ff7220784f7e93dccf70ac6bd355e74.

A separate acceleration candidate now evaluates the same pathwise moments
using fixed-standard-normal Hermite orders32/64/128, requesting256 only if
64/128 disagree. At least128 is required. It retains long-double conditional
arithmetic and moment summation; no optimizer or backend fallback is supplied.
This is explicitly different from the original adaptive20-node approximation.
Fixed-measure differentiation includes the eta/u chain; it does not reuse the
original frozen-mode score calculation. Order agreement is empirical, not
a certified whole-line quadrature error bound.
Sources: `repair_p3_standard_hermite_pathwise_v1.py` and
`check_standard_hermite_pathwise_v1.py`.
Three individual fixed-point comparisons were submitted: one saved full
eight-gene LL/gradient/Hessian comparison(8CPU/8GiB/90h), and the two saved
Decimal variance comparisons(each1CPU/4GiB/90h), all cpu/nslab.
The full comparison requires original173 counts, all eight genes, reference
LL/g/H agreement1e-9/1e-8/1e-5 and summed empirical order allowances
1e-8/1e-7/1e-5. It adds no fits or biological evidence.
After Claude's stated session reset, the source-only review was retried under
`claude_extended_variance_review_20261007T102159Z/`. The CLI subprocess
timed out after240s; no review response was saved. Source and request remain
available for the owner. No Codex subagents were retried.

All three standard-Hermite fixed-point checks passed. Full8=22084914 used
1792 nodes, took6.04129s, and differs from the saved GL LL/g/H by
1.90994e-11/1.62459e-10/2.68920e-8. Summed empirical order allowances are
0/1.40244e-12/2.07217e-9. Every gene selected128 after32/64/128.
The earlier full GL call took211.06543s on a different node; no controlled
speed ratio or full-fit speed claim follows. Independent variance checks
22084915/16 took6.09578/5.59946s and pass. Their variance-curvature errors
against Decimal are1.20777e-7/4.45860e-9. Root independently rederived all
global sums, order differences, reference comparisons and node counts:
`standard_hermite_pathwise_root_review_20261007T102819Z/review.json`.
Helper SHA292a6dad676ad2a461ea5b32435e85b5091759017932a4204002315244830a6d;
checker SHA5bd464430b6e3b863a570120675496d295c79afe03224a5fb97b5898ef277a2c.

Separate full-fit pilots22084966(draw0) and22084967(draw11) are now running
in `standard_hermite_pilot_launch_20261007T103340Z/`, each8CPU/8GiB/90h
cpu/nslab. Both passed preflight exactly: full8 results equal the checked
Hermite output, score-only and Hessian modes have identical LL/gradient,
and the original zero-variance boundary is byte-exact. Both draws reproduce
the prior173 counts before fitting. They use the exact supplied original
generating-null vector and unchanged original P3 draw AST, two starts,
section4 initialization, tie and strict1e-6 convergence. The original
constrained null solver runs before the separately checked refinement.
The only overrides are positive marginal integration and constrained polish.
Score-only calls avoid unused Hessian work. Every numerical refusal stays
failed; unknown exceptions are not accepted as numerical refusals.

New pilot sources:
- `repair_p3_standard_hermite_pathwise_v2.py`,
  SHA6f7c0af1c39652a881393801a74722ab3a5639d28abad3706854e34cbebe69ba.
- `repair_p3_standard_hermite_model_v1.py`,
  SHA8adb65ccafd2a05937b0874ea7dc7795f559b1a86093e2b6f806c83499ede5e9.
- `run_standard_hermite_pilot_v1.py`,
  SHAc5810b4e264f912b114343eb02a2cfd1a3b986a30f7af767015c29681fa7a03a.

These pilots do not replace any old draw or establish fifty-draw agreement,
calibration, human validation or a useful deployable biological model. The
synthetic scope remains G8/Z7/P28 with bounded counts. Larger real design
and identity inputs are still absent. AST, pinned metadata inspections,
launcher syntax and technical-document whitespace checks passed; Resource
scope validation not run because no paper, figure or release authority changed.

### Full-pilot refusal and exact scheduling comparison, 2026-10-07 11:05 UTC

All three fixed-normal Hermite full pilots are terminal prescribed failures:
22084966(draw0),22084967(draw11),22085062(draw18). Their runtimes are
910.46551/1180.21607/406.97994s; batch MaxRSS53472/53488/73080KiB.
Each raises HermiteRefusal on the first positive-variance evaluation, retaining
NaN LR and true exceedance. Earlier saved-point passes remain valid only for
those points. The draw18 pilot adds the separately checked parallel null shifts;
it does not qualify a fifty-draw family or repair positive integration.
The first two pilots were collected independently in
`standard_hermite_pilot_collection_20261007T110300Z/`:2 terminal,0 converged,
2 prescribed failures,qualification false. Collector exit1 is expected.

The original runners omitted the refusal point and completed null snapshot
on this caught-exception path. Current source runners now save both, while
returning the same null result and retaining the existing exception rule.
Archived executed runners and numerical helpers were not modified. Therefore
the three old failure points cannot be reconstructed exactly from their output.

New single-start diagnostic22085142 uses the older saved converged draw0 null
and exact saved draw0 counts, without repeating a null fit. It completed0:0
in16s(sampling MaxRSS0 is unusable), reproducing refusal at variance
0.9234818606206378 in all eight genes. Orders128/256 disagree by0.00301939
to1.64710332 in per-gene LL and0.06825770 to26.08242255 in maximum score.
This is not the actual lost pilot endpoint and not a full P3 refit. Complete
parameters/orders are saved in `standard_hermite_refusal_check_22085142/`.
Runner SHA d35d09304f976d55ddc21e091a35a1215168f21e6a374bf91c81962da892a743;
pins SHA83cec4cde41696f00d1f58bdc7496c095652887bc29f88e64a73c17259cd8360.

The exact original three-shift null loop was independently parallelized via
ordered Pool.map. Checks22085026/27 preserve every endpoint array byte-for-byte,
including convergence. Serial/parallel times are28.21685/15.82997s(draw0)
and35.86899/21.34445s(draw11). Draw11 remains failed. Source helper
`repair_p3_parallel_null_shifts_v1.py`, SHA
bd4e8bcbcc6d2b0ab34999c2220c51a96ee5badc0d44bd2cf75b72946eabb51c;
checker SHA5aa7fa59ba3e1320d5bd404509e9b80d3301d691d5b2bb0cf87672bb4a33e2e2.
Root: `parallel_null_shifts_root_review_20261007T104754Z/review.json`.
Only scheduling at two saved null starts is checked, not complete-fit speed.

New separate numerical candidate `repair_p3_centered_hermite_pathwise_v1.py`
uses defensive importance quadrature: positive prior weight, sampled/refined
posterior peaks and their local/wider normal proposals. Integrand weights are
alpha_j*phi(z)/q_mix(z); the original pathwise derivatives are integrated with
proposal parameters frozen. No model, count, optimizer, threshold or variance
bound changes. Whole-line order agreement remains empirical and cannot certify
unseen peaks. Four individual fixed-point jobs22085153–22085156 compare full8,
Decimal gene2/gene5 and the new refusal point. Resources are8CPU/8GiB for each
eight-gene case and1CPU/4GiB for each Decimal case, cpu/nslab/90h, no array.
Launch `centered_hermite_pathwise_launch_20261007T110325Z/`;
checker SHA4c0b2bc868a4dd2d261078f5dc33d221927664499b244873bf8acba8c9d072c1;
pins SHA17eb8d6753156ad82390ea776a373e20f5b7a7ad91ce39be77b80cca768af6a9.
The refusal case has no independent large-variance reference yet. No full fit
or biological validation is claimed for this candidate.

All four centered-importance comparisons are now terminal failures. Full8 has
LL/g/H reference differences2.17142e-11/2.38530e-10/9.44653e-7, but order
allowances fail(maximum score4.92052e-6,Hessian0.00240232). Independent
Decimal gene2 variance-curvature error5.70330e-5 exceeds1e-5. Gene5 reference
curvature error5.84338e-6 is within1e-5, but its order allowance fails. All eight
new refusal-point genes fail order agreement. Thus none of the four qualifies;
this revision will not enter full fitting. Runtime19.97226/30.58006/29.51125/
45.08149s; batch MaxRSS0(unusable)/56904/56888/187852KiB. Complete arrays,
failed allowances and source copies are retained. Root rederived global LL/g/H,
every order difference, proposal weight sums and all source/input pins:
`centered_hermite_root_review_20261007T110620Z/review.json`.

Root also verified all three original pilot failures, counts and source copies,
and the separate refusal's per-gene order differences:
`hermite_full_pilot_failure_root_review_20261007T110514Z/review.json`.
The two current runner edits are AST-identical to their executed predecessors
after removing only the new null-snapshot wrapper and refusal-output field.
Current source hashes91e2c9f283981e612a0b39341e2bc6c7f1bc0a76414920cf392f356881d6ac52
and78c11147ee94a5ace8e7c73b749b098f0cab8973cab35a4a104b65fc8709fdad.
Executable checks of these output-only additions were not run.

New reference job22085179,8CPU/8GiB/90h,cpu/nslab,no array, evaluates all eight
genes at the newly saved large-variance refusal point with adaptive intervals.
Detected peak centers and local widths supply breakpoints only. Each integral
uses the original normal measure and native selected-RNA pathwise moments;
both z domains24/26 must meet combined difference-plus-allowance limits
LL1e-9,g1e-8,H1e-5, including derivative-tail bounds. This supplies a numerical
quadrature comparison, not an independent conditional Hessian or scientific
calibration. No fitting is performed. Source
`check_refusal_adaptive_intervals_v1.py`,SHA
5a317770f7a055f11bd3e627053eb64cca3fb649be4e4e652a5133a7558c1763;
pins SHA1cd5e84cb291775ddbfc3467d4eabb33cc1717823e8a4320f51ef2ad86ba0be2.
Launch `refusal_adaptive_intervals_launch_20261007T110805Z/`.

Original fifty collection `repair_p3_gamma_order_family_collection_20261007T110400Z/`
remains47 terminal,36 converged,11 prescribed failures,3 running,0 pending,
one known draw0 disagreement. At11:06UTC active draws18/25/30 have batch
MaxRSS313564/313452/304060KiB. No result is reclassified. B3, C5, human
identity/design and independent model validation remain incomplete. AST,
source/input inspection and bash syntax checks passed for new launchers;
Resource scope validation not run because no paper/figure/release authority changed.

Adaptive reference22085179 subsequently passed all eight genes,38.31812s
compute/49s scheduler wall,MaxRSS121956KiB. Both domains and all propagated
allowances pass. Root rederived each local raw moment, embedded gradient/Hessian,
mass lower bound, ratio/covariance error propagation, domain difference and
global sum: `refusal_adaptive_intervals_root_review_20261007T110956Z/review.json`.
Summed per-gene domain difference-plus-allowance maxima are
8.38950e-13/1.07324e-10/4.76838e-8 for LL/g/H;10248 adaptive node evaluations
over both domains. These estimated allowances do not certify floating arithmetic
or constitute independent conditional derivatives.
At this same point standard GH256 has global LL/gradient errors1.06574/46.81722.
Centered GH256 has LL/g/H errors0.00058875/0.74431/25.79578. Both failures are
therefore substantive numerical disagreements, not merely conservative order
allowances. The new reference is still a fixed point only. Next required checks
are the saved full8 GL point and both independent Decimal variance points,
followed by a separately labeled full-fit comparison if they pass. Do not
restart a fifty-draw candidate from this single-point success.

### Adaptive reference family and full-fit pilots, 2026-10-07 11:20 UTC

The unchanged adaptive worker was next checked against full8 GL and the two
independent Decimal points, jobs22085195/96/97. All three fail the combined
Hessian domain allowance: full8 gene2=1.13624e-5, Decimal gene2=1.78968e-4,
Decimal gene5=2.60316e-5, versus1e-5. Their numerical reference differences
pass: full8 LL/g/H1.92131e-11/1.62274e-10/2.60189e-8; Decimal variance
curvature4.66286e-8/2.01851e-8. The failure classifications are retained.
Source `check_adaptive_interval_references_v1.py`,SHA
c8286f45abd8a8f2c503346f84e3bb0f96346f9980a5384a8332a86d71d83d96;
launch `adaptive_interval_references_launch_20261007T111240Z/`.

Separate worker `repair_p3_adaptive_intervals_weighted_v1.py`,SHA
4d55a6bad630db3e92a6824de91cc88ed4ef816f20c1f257f8bb811396f32f5a,
changes only sampled positive conditioning units from abs(moment) to
exp(log_integrand-log_shift)*abs(moment), with minimum1. These units are
arbitrary numerical conditioning, not certified bounds. Explicit derivative
tail envelopes, original normal measure, pathwise derivatives, domains24/26,
quadrature requests and final allowances are unchanged. Static AST comparison
confirms the worker differs only inside the conditioning-scale sampling loop.

All four weighted comparisons22085200/01/02/03 pass. Full8 reference LL/g/H
differences1.92131e-11/1.62374e-10/2.64554e-8; independent Decimal variance
curvature errors5.10524e-8/1.91447e-8. Large-variance LL/g/H differences from
the prior adaptive reference are0/1.70531e-13/1.30968e-10. Summed per-gene
domain difference-plus-allowance maxima(LL,g,H) are:

| Point | LL | Gradient | Hessian | Compute seconds |
|---|---:|---:|---:|---:|
| full8 | 1.50283e-12 | 2.36093e-10 | 1.91402e-7 | 25.20053 |
| Decimal gene2 | 1.52253e-13 | 4.94272e-10 | 5.93604e-6 | 50.71143 |
| Decimal gene5 | 1.93748e-13 | 1.47875e-10 | 2.21323e-6 | 18.51913 |
| large variance, eight genes | 6.05123e-12 | 8.24832e-10 | 6.10833e-7 | 42.64940 |

Root rederived raw moments, covariance/ratio allowances, global sums and all
reference differences in `adaptive_weighted_root_review_20261007T111822Z/`.
Supplemental `adaptive_weighted_domain_review_20261007T112033Z/` rederived
every per-gene domain comparison. Pins and source copies match. Scheduler walls
26/53/19/44s; batch MaxRSS0(unusable)/54864/0(unusable)/125400KiB.
Checker SHA c111a2d92532e46e52264ae2eae11267d6858844a2ac1cd72766653cbca3b4c3;
pins SHA24cec951acb0f9be7c913f33d49fcebb22de9df64419dc17ab7f59402199e348.
Launch `adaptive_weighted_references_launch_20261007T111456Z/`.

Full-fit pilots22085227(draw0),22085228(draw11),22085229(draw18) are submitted
individually8CPU/8GiB/90h,cpu/nslab,under
`adaptive_interval_pilot_launch_20261007T112000Z/`.
They retain original P3 counts/generating null/draw AST, warm/cold starts and
selection, strict1e-6 convergence, NaN LR/exceedance failure accounting and
the previously checked original-first null refinement. Independent original
null shifts are scheduled in parallel. The separate positive backend recomputes
peak breakpoints at each point and integrates full LL/g/H on both domains,
even for score-only calls. Every gene and the summed allowances must pass;
there is no backend fallback. The original zero boundary is inherited.
Preflight requires the saved full8 endpoint, identical score/Hessian-mode LL/g,
and byte-exact original zero boundary before draw fitting. The new runner saves
null snapshots and refusal points. These are separate pilots, not a fifty-draw
candidate, calibration or human validation.

First launcher packet `adaptive_interval_pilot_launch_20261007T111929Z/` failed
metadata inspection because the draw18 count-reference path was misspelled.
No job was submitted from that packet. The root runner path was corrected and
all three metadata inspections and bash syntax checks passed before submission
from the new packet. Archived sources and existing results were preserved.
Executed runner SHA dff4af2ce976aef57671f3be34d6303287ea169991f85d4ab84cfe1110a71ef7;
pins SHA801d698dc69bd2703db05997be4b3994a10175b3f8ea3056ea876c71cc4d3108.

Original50 collection `repair_p3_gamma_order_family_collection_20261007T111540Z/`
remains47 terminal,36 converged,11 prescribed failures,3 live,one known draw0
disagreement. Older Gamma collection `repair_p3_gamma_draw_collection_20261007T111540Z/`
now has49 terminal,1 incomplete,0 collection issues; it remains a different
numerical family. The older draw18=22081587 completed0:0 in7:20:24,
MaxRSS57944KiB. B3/C5 and qualified human inputs remain unresolved.

All three new full pilots passed their preflight before fitting. Dynamically
recomputed peak breakpoints give full8 differences0/1.42109e-14/1.13687e-13
in LL/g/H. Score-only and Hessian-mode LL/g are identical; the original zero
boundary is byte-exact. Counts are independently regenerated from the supplied
uniforms and match each prior saved draw before entering the unchanged draw AST.
At11:21UTC batch MaxRSS244336/247156/242600KiB. Each is RUNNING with its
fit-start event saved; no full fit has returned yet. New helper sources are
`repair_p3_adaptive_peak_scan_v1.py`(sampled/refined peak scan only) and
`repair_p3_adaptive_interval_model_v1.py`(bounded positive integration wrapper).
The peak scan retains the previously executed scan; the full-integral worker
is the checked probability-weighted source. AST and metadata inspection,
launcher syntax and technical-document whitespace checks pass. Resource scope
validation not run because no paper/figure/release authority changed.

### Saved positive-LR endpoint, 2026-10-07 11:32 UTC

The older Gamma collection49 terminal records contains28 converged fits,
21 prescribed failures and two positive LRs: draw18=2.176064620750708,
draw25=0.7101477087253443. Draw18=22081587 has converged selected null and
marginal fits, variance0.007118859016508785 and native20/40 LL/score agreement.
This supplies a positive endpoint for the independently integrated numerical
comparison; it remains a synthetic endpoint, not a calibration result.

New fixed-endpoint job22085258 completed0:0 in18s,15.77432s compute.
Sampled MaxRSS0 is unusable. It uses the exact original173 counts and returned
parameters, recalculating adaptive full8 LL/g/H and the original zero-variance
null endpoint without any fitting. LL/g differences from the saved native fit
are1.35288e-11/4.79191e-11; every gradient component passes the specified
1e-8+1e-6*abs(reference) criterion. Projected gradient1.09980e-7 is below
1e-6; recomputed null nuisance projected gradient2.66669e-9 also passes.
The same-point LR is2.1760646207777654, differing by2.70575e-11 from the
saved fit. Summed domain difference-plus-allowance maxima are
LL1.35015e-12,g8.99578e-11,H2.31539e-8. Root rederived all global sums,
projections, LR, comparisons and count/source hashes in
`adaptive_positive18_endpoint_root_review_20261007T113129Z/review.json`.
No full-refit agreement, parameter estimation replication, independent full
conditional Hessian, calibration or human validation is established.
Runner `check_adaptive_positive18_endpoint_v1.py`,SHA
86be3c61bc8012a40ccbecc3ae741cf55fc8cbc818e8eb28b1fd8b0f5f082564;
pins SHA a581e3eedca99cd926d2ca640e7686e0575c8bfee2784410e7aada487ea868dd.
Launch `adaptive_positive18_endpoint_launch_20261007T112923Z/`.

All three new adaptive pilots have saved converged warm-selected null fits.
Inspected progress gives6 positive calls for draw0 and3 each for draws11/18,
with full-moment/domain allowances passing. Draw11/18 include multiple detected
peaks at some evaluated points. These runs have moved beyond the positive
integration refusals of the fixed-normal Hermite pilots. All remain live and
full-fit convergence/performance is still unverified.

New stdlib collector `collect_adaptive_interval_pilots_v1.py` reads only this
fixed three-job roster and pinned synthetic outputs. It verifies sources/counts,
preflight, positive allowances, saved null selection and, once terminal, original
failure/exceedance/NaN LR accounting, final null identity, marginal start selection,
recomputed LL/score and projected convergence. IntervalRefusal must include the
actual parameter/variance point and all eight gene records. The collector has
been executed on live outputs with null snapshots and positive calls; terminal
branches have not yet executed. Latest collection
`adaptive_interval_pilot_collection_20261007T113100Z/` is0/3 complete,
qualification null,incomplete,with no collection issues. No paper, figure,
release authority or protected outcome changed. B3/C5 and qualified human
inputs remain outstanding.

### Full draw11 and requested derivative modes, 2026-10-07 11:48 UTC

Adaptive draw0=22085227 is a terminal prescribed numerical failure after13
fitting-positive calls(15 progress calls including preflight). Runtime897.88596s,
scheduler wall14:59,MaxRSS247332KiB. At saved variance
3.344593761911824e-5, summed LL/g/H difference-plus-allowance maxima are
1.64365e-12/1.22644e-9/1.61278e-5. Only the aggregate Hessian allowance
fails. Each gene separately passes. This occurred during an optimizer request
for likelihood/gradient, because the first wrapper computed and required a
Hessian on every positive call. Its NaN LR and true exceedance are retained.
The saved failure includes exact parameters and all eight gene records.

A separate requested-mode worker now retains full-Hessian arithmetic when
requested and omits second moments/curvature bounds when only LL/g is requested.
Scientific LL/g/H limits, explicit tails, conditioning policy, domains and
quadrature requests remain unchanged. Source
`repair_p3_adaptive_intervals_modes_v1.py`; checker
`check_adaptive_modes_v1.py`. Five individual checks22085277–22085281 pass:
the prior full8, two Decimal and large-variance references, plus the actual
draw0 failure point. Full-Hessian outputs are identical to their previous
checked records. Score-mode LL matches exactly at all five points; maximum
mode-gradient difference2.84218e-14. At the failed point score allowances are
LL1.13699e-12,g9.11990e-10, while its original aggregate Hessian failure
remains. The fifth check qualifies only requested LL/g at that point.

Root rederived score raw moments, ratio allowances, every domain comparison,
global sums and mode differences, and compared full-Hessian arrays against
the prior outputs: `adaptive_modes_root_review_20261007T114414Z/review.json`.
All source/input pins and archived sources match. Checker SHA
a9f074ad0eb0d31c06eaf472d48d32e5002c72ce9bc74f8ee7aeb23ed6481f62;
pins SHA88e2e03aff8193b80c6b9683f9a54de06d6049c9d51adfc324420a027d8bc8e2.
Launch `adaptive_modes_launch_20261007T113906Z/`. Scheduler walls
1:47/0:33/0:29/1:13/0:32; MaxRSS125076/55248/0(unusable)/126884/114624KiB.
No speed comparison or uniform optimizer-path accuracy is established.

Separate mode-aware draw0 pilot22085287 is running,8CPU/8GiB/90h,
cpu/nslab,no array. New wrapper `repair_p3_adaptive_interval_model_v2.py`
honors the original with_hessian argument and keeps the strict full-Hessian
criteria when requested. Original fitters, null-shift scheduling, starts,
selection, convergence and failure rules remain. Runner
`run_adaptive_modes_pilot_v1.py`,SHA
64688f3bb5ba3ed8a2b4e6320f5551f483a197f25b1d2a350e45203efd9b0c95;
pins SHA9bb83f99f4a09324a790079a5b552a2d02f88292a6cf31c2efa2106c3b662359.
Launch `adaptive_modes_pilot_launch_20261007T114530Z/`.
The former failed draw0 remains intact; this is a separately labeled candidate.

Adaptive draw11=22085228 completed0:0 in24:03,MaxRSS248364KiB,
compute1432.72908s. Selected null and marginal fits converge; both retain warm
under the original tie. Variance0,LR0,loglik-512.189905853584. Fresh LL
difference0 and both projected-gradient maxima6.60281e-7<1e-6. Root verified
the null snapshot, source copies, selection arithmetic and projected scores:
`adaptive_full_draw11_root_review_20261007T114718Z/review.json`.
This does not change the earlier original draw11 disposition or qualify fifty
draws, calibration or independent human evidence.

Collector `collect_adaptive_interval_pilots_v1.py` has now exercised the actual
terminal failure and success paths. Collection
`adaptive_interval_pilot_collection_20261007T114600Z/` is2/3 complete,
1 converged,1 prescribed failure,1 live,with no collection issues. Exit1
preserves the failure. The collector separately reports preflight-inclusive
positive calls and fitting-only calls. Draw18 remains live. B3/C5 and
qualified human inputs are still missing; no paper or protected outcome changed.

### Claude source review and uniform integration comparison, 2026-10-07 12:24 UTC

Source-only Claude review completed exit0 in210.671s with tools disabled and
source hashes unchanged. Response/receipt:
`claude_adaptive_math_review_20261007T120531Z/`. It reports no demonstrated
native selected-RNA likelihood or pathwise derivative error. The adaptive
domain comparison shares an interior partition, so an additional integration
rule is needed to address a missed narrow mode at a checked point. This review
does not certify model calibration, cross-gene independence, full conditional
curvature or human validation.

Two bounded synthetic follow-ups completed and pass. Uniform GL22085533
integrates every grid interval, without the adaptive peak scan, and uses fixed-u
Gaussian density derivatives for variance. At the saved variance0.9234818606206378,
all8 genes return and global LL/g/H differences from the adaptive reference are
0/1.7053025658242404e-13/1.0186340659856796e-10. Both domain/grid refinements
cover all intervals. Runtime238.657814s, scheduler4:02,MaxRSS143372KiB.
It shares native conditional derivatives; neither floating error bounds nor
uniform optimizer-path agreement follows from this fixed-point comparison.

Boundary11 check22085532 independently sums .5*(h_uu+h_u^2) at u=0 using
native long-double conditional derivatives. Score-219.18993211473608 agrees
with the original variance score within1.219007117470028e-10; likelihood
difference3.637978807091713e-11. The raw score is negative. Runtime0.636701s,
scheduler9s; MaxRSS0 is unusable. Source inspection confirms the original v=0
score is d/dv, both implementations use eta-orientation*omega, and the fixture
loader casts selection thresholds to integers. Both orientations and173 nonzero
fixed offsets occur in the checked rows. No boundary convention change is needed.

Root independently rederived heat terms, LL/g/H sums and differences, checked
uniform interval geometry and verified all source/input pins and archives:
`claude_math_followups_root_review_20261007T122422Z/review.json`.
Checker SHA336af2a487a8ae81990d60db548476ef8193ed460be9c9be2144345a32672a5e;
pins SHAc478fc8f79f2719425a9d32b428dd5810ce007d78e9654bc854ed370006e6d4a;
launch `claude_math_followups_launch_20261007T121801Z/`.
The participant dependence limitation already stated above remains. These checks
neither define missing B3 amplitudes/allocations or C5 designs nor open protected
outcomes. The full draw0/draw18 pilots remain live at12:23UTC; no new family
qualification, calibration, biological inference or scientific adoption occurs.

### Dense five-gene comparison: tenth point complete, 2026-10-07 16:06 UTC

Older dense GL job22082313 completed0:0 in10:45:37, MaxRSS77660KiB.
The synthetic200-row gene3 point at tau0.0125 differs from the independent
scalar logmass by-7.389644451905042e-13; native20/40 comparisons also pass.
Returned17 scores and17x17 Hessian are finite and their empirical allowances
pass. An intermediate multiplication overflow warning remains recorded; this
does not provide independent gradient/Hessian accuracy or certified floating error.
Collection `repair_p3_gamma_GL_five_gene_collection_20261007T160800Z/collection.json`
is10/15 complete,0issues,5running. Root rehashed collection inputs and rederived
raw arrays, scalar differences, empirical allowances, inactive derivatives and
conditional-call bounds in
`dense_G5_gene3_quarter_root_review_20261007T160600Z/review.json`.
No job submitted/cancelled, protected data opened or paper authority changed.
Calibration, B3/C5 specifications and qualified human validation remain incomplete.

### Requested-mode draw0 full-fit result, 2026-10-07 12:41 UTC

Separate pilot22085287 completed0:0 in54:34,MaxRSS246016KiB,
compute3265.314420s. Both selected fits converge with original warm/cold
selection preserved. Null cold-minus-warm LL is approximately-4.9e-12;
marginal cold-minus-warm is approximately1.14e-11. Both retain warm under
the unchanged1e-6 tie. Variance0,LR0,LL-511.5784700492042. Root rederived
null and marginal projected-gradient maxima1.060368161331926e-9; native
variance score-541.6573915156132 points out of the nonnegative domain and
is correctly projected to zero. Active nuisance-bound gradients are also
projected under the original rule. Fresh LL difference0 and gradient equality
are recorded. All26 fitting-positive calls(28 including preflight) pass the
requested-mode allowances. No caught exception and no failed-draw flag occur.

Root evidence `adaptive_modes_full_draw0_root_review_20261007T124124Z/review.json`
also checks original173-count identity, source/input hashes and archives, null
snapshot, start-selection arithmetic, original bounds, derivative modes and
LR/exceedance accounting. This separately labeled revision does not replace the
earlier full-moment draw0 failure. It supplies neither fifty-draw agreement nor
calibration or human validation. Draw18=22085229 remains live. Quantitative B3,
C5 reduced designs and qualified participant/design inputs remain outstanding.
No scientific release, paper claim or protected outcome changed.

### Complete native-Gamma fifty-draw family, 2026-10-07 12:43 UTC

Last draw30=22081599 completed0:0 in8:59:53. The complete family has29
successful fits and21 prescribed failures. Every successful fit passes the
saved strict projected-score and native20/40 diagnostics; no collection issue
or unknown caught exception is reported. Positive LRs are draw18=2.176064620750708,
draw25=0.7101477087253443 and draw30=1.2893636931730725. All50 exceed the
synthetic observed LR0, giving original fixture-recipe p=1. Failures retain
NaN LR and true exceedance under the original rule. This is not calibration,
an independently qualified fifty-draw comparison or biological validation.

Collection `repair_p3_gamma_draw_collection_20261007T124300Z/collection.json`;
root `gamma_complete50_root_review_20261007T124255Z/review.json` independently
checks all collection input hashes and rederives successful-fit convergence
flags, warm/cold selection, LRs and the complete original recipe. The independent
comparison remains outstanding. Separate adaptive-GL family remains47/50 terminal
with3 live; neither new pilot success replaces a member of either family.
No missing quantitative B3/C5 or qualified human input is supplied by completion.

### Independent continuation review: saved draw11 precision, 2026-10-07

The separate synthetic checks22087829/22087880/22087909 are terminal and
show a floating-arithmetic sign reversal at two fixed null parameter vectors.
The saved double likelihood reports a full Newton-step decrease8.42419e-11;
independent finite-product Decimal likelihoods at80 and110 digits instead
give an increase9.572939710743758e-12. Independent projected scores are
2.127101053885565e-6 at the original failed point and2.048385519147154e-9
at the proposed point, against the unchanged1e-6 cutoff. Original failure
accounting remains unchanged; no fit or proposal was accepted in this review.

Root independently rederived all173 row likelihood/28-coordinate score sums
at both precisions, the bound projections, and the point difference. Double
and Decimal active sets agree; both score comparisons pass the specified
atol1e-8/rtol1e-6. The step infinity norm is4.979632821555668e-5. Its
endpoint-score trapezoid gives9.572717185695982e-12, differing from the
independent likelihood change by2.22525e-16. Evidence is in REC
`continuation_decimal_stationarity_review_20261007T163701642220Z/`.
This establishes arithmetic sign reversal at this point. The production
line-search rejection of this exact direction was not replayed; rowwise
double likelihood comparison and an independent Hessian were not run.

Source-only Claude review completed0:0 with tools disabled in135.288s:
`claude_continuation_spec_stationarity_review_20261007T163556586579Z/`.
It found no formula error in the probability/score checkers and recommends
reviewing precision-aware Newton polishing under v1.7 item8. Its proposed
acceptance tolerance is advisory, not adopted or certified by these checks.
Any implementation requires separately recorded complete-fit comparisons
with the original starts, selection, cutoff and failure accounting.

The same review and an independent agent confirm that the exact B3
warp-speed construction, allocation, remaining quantitative plants and P2
coverage target are absent from the supplied controlling texts. The C5
helper already exists and passed its synthetic checks; the owning agent
still needs to supply identity-qualified reduced design packets. No new
calibration recipe, private computation, human replication or paper claim
is established by this review.

### Extended-precision null-point comparison, 2026-10-07

Separate job22087992 completed0:0 in14s (13.328s internal runtime), with
self-reported maximum RSS83028KiB. Scheduler MaxRSS0 is unavailable. Both
fixed173-row null vectors satisfy the existing helper's32-row/100000-count-
entry per-gene caps. The unchanged extended-precision conditional helper
returns float64 row quantities; native predictors and chain-rule products
remain float64. Summing these contributions in longdouble gives an improving
step9.58810e-12 versus the independent Decimal9.57294e-12. Float64 summation
of the same improved row contributions also recovers the positive sign,
9.54969e-12. Thus improved internal row arithmetic resolves the sign here;
extended aggregation alone is not established as the explanation.

Both complete gradients agree with Decimal80/110 within3e-14. The reference
projection tolerance1e-10 and original null-fitter tolerance1e-8 each select
only coordinate27 as active and retain the original failed/trial-passing
score distinction. Root rederived row sums, reference differences and both
projections, with all before/after input pins matching. No independent
Hessian, optimizer-path accuracy or complete fit is established. Evidence:
REC `check_longdouble_null_points_v1_22087992/` and
`continuation_longdouble_null_points_review_20261007T165134621831Z/`.
Reviewed checker SHA8dedd618015ec67a0d34f8941678e70472b0062d01dd9f7ad799b89fe0a56615;
launch `continuation_longdouble_null_points_launch_20261007T165002191731Z/`.

### B3 method definition: primary-source check, 2026-10-07

The [Giacomini–Politis–White author manuscript](https://mathweb.ucsd.edu/~politis/PAPER/Warpspeedfinal.pdf),
sections2–3, pools one bootstrap root per independently simulated dataset.
At fixed sample size, its pooled distribution estimates the expectation of
the conditional bootstrap distribution. Its coverage justification uses
sample-size/simulation limits and conditions on the relevant quantiles; the
authors warn that pooling can miss failures from variation in conditional
bootstrap distributions. This is distinct from the one-second-level-resample
inference procedure studied by [Chang–Hall](https://arxiv.org/pdf/1408.6327).

[Inference] Setting the retained Monte Carlo p-value formula's inner B to1
gives only0.5 or1, so it cannot itself assess the specified0.05 rejection
rule. Pooled critical values define a different simulation diagnostic, and
one ratio draw does not form each dataset's prescribed90% percentile
interval. The owner must define the exact root, resampling/conditioning law,
coverage endpoint and relationship to the200 full-CRT datasets. Applicability
to the selected-count model, variance boundary, invalid lambda draws and
Holm family is not established by the method name. No calibration rule or
workload is changed by this source check.

### Independent supplied-P3 count reconstruction, 2026-10-07

Job22087964 completed0:0 in8s. Its finite-product Decimal80/110 calculation
uses the fixed generating null, supplied uniforms, inclusion-reweighted
three-genotype mixture, original integer support and smallest-CDF-index
rule. Both precisions reproduce all50 saved vectors:8650 draw-row counts
agree. Genotype classes agree between the two Decimal precisions; no saved
native genotype-class comparison is claimed. The nearest uniform-to-selected-
bin boundary is3.35593e-7. Scheduler MaxRSS0 is unavailable.

Root reviewed the complete checker and independently rechecked all input
hashes, both precision outputs, all50 vectors and their little-endian int64
hashes. Evidence: REC `p3_decimal_counts50_22087964/` and
`continuation_independent_P3_counts_review_20261007T165743256186Z/`.
Checker SHA401afa8c7d5a4ed1bafc40291b68762e326c85dab0ba65013f93ced85a6e08be;
summary SHA73e329840274236f24cb6f03ae6d60cf5432ba682ccaf26792316f8cdf062ddb.
This qualifies count reconstruction for the conditional agreement fixture.
It supplies no fit-convergence agreement, new independent simulated datasets,
statistical size/power/coverage, or human replication.

### Original-optimizer extended-precision null pilot, 2026-10-07 18:18 UTC

Separate synthetic job22088824 is RUNNING. It uses one individual cpu/nslab
allocation,1CPU/2GiB/90h, with all four numerical thread variables fixed to1,
an88h internal cap and32MiB saved-output cap. CPU partition maximum14d was
checked before submission. No array or new dependency is used.

Runner `continuation_longdouble_null_fit_v1.py`, SHA256
`af4bd5fb5f1ae19df6abcbfb26cb5c2d61113baa62ddba80fa80c139faf7a74b`,
changes only evaluation arithmetic: unchanged longdouble conditional helper,
native float64 predictor/row/chain-rule interfaces, longdouble accumulation,
float64 optimizer return. Original optimization, cold-start construction,
constrained polishing, acceptance, bounds, two starts,1e-6 likelihood tie and
strict1e-6 score cutoff are inherited. The instance fit_single wrapper only
records and delegates to the original bound method. It rejects row subsets
and nonzero random effects; no positive-variance model is fitted.

Root and a separate source reviewer found no blocking fit-path defect.
AST,47-pin metadata inspection and launcher syntax passed. Frozen launch:
REC `continuation_longdouble_null_fit_launch_20261007T181723798578Z/`.
The actual evaluator passed both saved-point comparisons at80/110 digits;
full gradient errors remain below3e-14. Output:
REC `continuation_longdouble_null_fit_v1_22088824/`.

Independent endpoint scores are recorded rather than enforced as an overall
completion condition. Root must therefore rederive endpoint agreement,
original-bound projections and start selection before claiming a successful
numerical repair. Execution completion alone will not prove these conditions.
The completed Decimal-acceptance pilot22087946 is a distinct numerical route;
neither route replaces an original failed draw or establishes complete-family
fit agreement, B3 calibration, qualified participant identity or human replication.

Root independently rechecked all48 saved source/input hashes, both parameter
vectors, both precisions' complete28-coordinate gradients, likelihood differences
and the existing comparison tolerances from saved numbers only:
REC `continuation_null_prefit_root_review_20261007T182041627170Z/review.json`.
All prefit comparisons pass. No model was called on login, and this review
does not establish terminal convergence or independent Hessian accuracy.

### Scientific-estimate stability of the Decimal null repair, 2026-10-07

Root compared the completed separate pilot22087946 with the saved failed
draw11 null, without rerunning a model. Retained warm kappa changes from
-0.10495811480799731 to-0.10495807342495057, a difference4.138304674e-8;
cold differs from the original by4.138053281e-8. Both are below the existing
1e-4 kappa agreement tolerance. Maximum parameter changes are approximately
4.98e-5 at coordinate12, a gene-dispersion coordinate. The independent
110-digit likelihood gain over the original failed warm point is positive
9.57294e-12 for each repaired start, despite negative native likelihood
differences. Previously checked independent projected scores remain below
the unchanged1e-6 cutoff.

Evidence: REC
`continuation_precision_null_estimand_review_20261007T183411701129Z/review.json`.
All34 original source/input pins and the exact baseline-reference parameter
identity were rechecked. This is saved-point estimate stability. The failed
original is not a qualified comparison fit, and this result cannot be counted
as independent complete-fit/fifty-draw agreement or calibration. Original
failure accounting and all scientific thresholds remain unchanged.

### Original-optimizer warm start returned, 2026-10-07 18:36 UTC

Job22088824 saved its warm terminal after3144 evaluations, including13 full
Hessian calls. Root rechecked all48 input/source hashes, exact observed start,
the original28-coordinate bounds and original1e-8 boundary projection.
The revised evaluator's projected maximum is6.319933598e-8, below the
unchanged1e-6 cutoff, agreeing with its saved convergence flag. Coordinate27
is active. Kappa is-0.10495807651546842, differing from the separately repaired
Decimal-polish warm estimate by-3.090517842e-9. The native likelihood difference
between these numerical routes is2.705746738e-11.

Evidence: REC
`continuation_longdouble_warm_root_review_20261007T183658575590Z/review.json`.
These checks use saved numbers only. Cold start, final warm/cold selection
and independent Decimal endpoint LL/gradient checks remain pending. No
independent endpoint Hessian, complete-fit agreement, calibration or human
validation is established by this partial result. Original failures remain.
