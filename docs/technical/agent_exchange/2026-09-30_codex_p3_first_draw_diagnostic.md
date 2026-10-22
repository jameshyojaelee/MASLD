# Codex: first supplied P3 draw diagnostic, 2026-09-30

## Observed problem

`bfix_p3_21998519/p3_000.json` contains `fit=null`, `LR=NaN`,
`exceed=true`. The role driver catches `ValueError`, `FloatingPointError`,
and `numpy.linalg.LinAlgError` without saving their cause. This record does
not establish optimizer nonconvergence or identify an exception source.

The observed-fit comparison passed with zero differences, as reported by
root. This diagnostic does not repeat or alter that comparison.

## Prepared reproduction

New files, confined to the Codex code directory:

- `diagnose_p3_first_draw_v17.py`
- `run_diagnose_p3_first_draw_v17.sbatch`

The diagnostic reproduces only supplied synthetic draw 0. It uses the exact
archived numerical modules in `bfix_p3_21998519/executed_sources/`, the
saved observed null parameters, original globals and original P3 uniforms.
All archived source hashes, role source and fixture hashes must match the
P3 role receipt. Imported numerical module paths and hashes must match
that archive. The saved observed null hash is recorded and rechecked.

It logs entry/completion for every conditional row draw, verifies the
original row ordering/uniform shape, then calls exactly:

```python
Model(rows).fit_marginal(settings['tau_max_fixture'], null_p,
                         section4_from_null=True)
```

No optimizer, starts, likelihood, bounds, convergence, draw, or statistical
rules are changed. Observational Python profiling logs entry into the
unchanged section-3/section-4 fit methods and their starting vectors. It
counts marginal calls and retains the last evaluated variance/parameters.
On an exception it saves the complete traceback, execution phase and
scalar archived-frame diagnostics such as gene, variance, mode and
curvature when present, then re-raises the actual exception. A scheduler
failure caused by that re-raise is therefore an expected diagnostic outcome.

Success instead saves the unchanged returned fit, convergence flags and LR.
Either outcome goes in a new `p3_first_draw_diagnostic_<jobid>/` directory
under the REC execution root, with receipt, executed diagnostic source,
progress and summary. Existing directories are refused. Original P3
outputs remain unchanged; the reproduction is not a replacement draw or
a scientific calibration result.

## Source inspection and checks

The archived conditional generator keeps every original row and uses
left-sided CDF search as prescribed. Its fitted marginal calculation
contains an explicit nonconcave-mode refusal, but **the saved P3 failure
does not prove that this refusal caused draw 0's exception**. No numerical
explanation is claimed before the diagnostic runs.

Executed on the login node: `diagnose_p3_first_draw_v17.py --inspect`
passed receipt/source/fixture hashes and Python syntax without importing
numerical modules, opening synthetic arrays numerically, fitting, or
creating result directories. `bash -n` passed for the launcher.
The saved observed null hash is
`a617d331fd3fb93576e5bd8588581aaab8f77946936536401905490438660181`.

## Actual submission

One individual job, name 1591; no array. Header:

```text
--account=nslab --qos=nslab --partition=cpu
--cpus-per-task=1 --mem=2G --time=90:00:00
```

```bash
sbatch /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_diagnose_p3_first_draw_v17.sbatch
```

Root submitted **21998696 (1591)** on 2026-10-01 UTC, with the header
above and logs `1591_21998696.out` / `1591_21998696.err`. Root changed
only the launcher's log names to include the actual job ID. Shell syntax
passed before submission; the diagnostic source remains unchanged.
No existing role runner, numerical source, or original draw was changed.
The numerical reproduction subsequently finished with the captured exception;
see the completed result below. The original P3 failure remains unchanged.

## Completed exact draw0 diagnostic

Job **21998696** finished after42:20 with exit1, following the intentional
re-raise of `ValueError: nonconcave mode in adaptive quadrature` at archived
`likelihood_v1.py:423`. Peak RSS was52,208KiB. The exception occurred in the
warm marginal fit's third marginal evaluation, gene2, v=0.9234818604998617,
with recorded inner point u=0 and curvature763.1283114221253. The complete
28-parameter outer vector and archived traceback are in
`p3_first_draw_diagnostic_21998696/summary.json` under the REC execution.
This establishes the refusal that terminated this exact synthetic replay.
It does not establish stationarity at u=0, whether another maximum exists,
or the accuracy of the integral. An independent bounded check is now being
prepared at those recorded parameters; no numerical source was changed.

## Source-only numerical review; actual cause pending

This section concerns the archived source, not the diagnostic job's eventual
exception. Root reports diagnostic **21998696/1591** RUNNING, with all 173
conditional rows generated and the unchanged marginal fit underway. No fit,
array calculation, job action or numerical-source edit was performed during
this review.

Line references below refer to
`bfix_p3_21998519/executed_sources/likelihood_v1.py` (SHA256
`41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6`).

### What curvature means in this model

For gene g at fixed outer parameters theta and v>0, the inner log-integrand is

```text
H_g(u) = -u²/(2v) + sum_i ell_ig(eta_ig + S_i u, r_ig).
H'_g(u) = -u/v + sum_i S_i d ell_ig/d eta.
H''_g(u) = -1/v + sum_i S_i² d² ell_ig/d eta².
```

These are the implemented quantities at lines 405–411. Each `ell` is the
log three-class mixture probability of the observed count, conditional on
that row passing the inclusion predicate (lines 175–186). The mixture
Hessian contains the positive semidefinite term
`w_H(1-w_H) score_H score_H^T` (line 138). Its inclusion-normalizer Hessian
is then subtracted (line 186). Consequently neither the row likelihood
nor H_g is guaranteed globally concave by the source model.

For a stationary univariate point, strictly positive H'' indicates a local
minimum, not an inner maximum. At a nonstationary point it only describes
local curvature. Neither observation establishes a violated biological or
statistical assumption. At exactly zero curvature the adaptive Gaussian
scale is undefined; a degenerate maximum is mathematically possible and
would require a separately justified integration method.

The proper Gaussian prior gives useful tail behavior: each conditional
count probability is at most one, so `ell_ig <= 0` and
`H_g(u) <= -u²/(2v)`. For fixed finite parameters with defined conditional
probabilities the integral is finite and the continuous log-integrand has
a finite maximum. Local nonconcavity does not invalidate that integral.
This analytic argument is not a numerical proof about the pending draw.

### What the current search verifies

The inner search starts only at u=0 (line 412). It uses
`step=score/max(-curvature,1e-8)` and at most 30 halvings (lines 415–420).
It stops on no accepted step, an accepted step smaller than 1e-10, or the
100-iteration limit (lines 413–421). There is no final residual-score test,
no other-mode comparison, and no explicit success flag for this inner search.

When curvature is nonnegative, replacing its denominator by 1e-8 can
propose a very large ascent step. Finite backtracking may fail even when
a smaller valid ascent step exists. A tiny step or exhaustion of the
iteration limit also does not prove stationarity. Negative final curvature
alone therefore does not certify either a converged mode or the dominant mode.

The outer projected-gradient check at lines 520–525 is a separate parameter
fit criterion. It does not prove inner-mode convergence or quadrature
accuracy. The exception at line 423 prevents computing an imaginary or
infinite scale. At present any such exception propagates from an outer
candidate evaluation (lines 483–490) and aborts the returned fit; the
original P3 driver subsequently counts that failure as an exceedance.

Existing `test_adaptive_integral_and_gradient` compares a separate small
synthetic example to direct scalar integration and GH40. Those useful checks
do not cover every optimizer proposal or establish the pending draw's
quadrature accuracy; no new numerical check was executed here.

### Repair recommendation, conditional on the actual diagnostic

The source establishes missing inner convergence certification, but it does
not establish which evaluation failed draw 0. No exact numerical repair is
claimed yet. If the traceback identifies line 423, first evaluate the recorded
outer parameters and gene using the same synthetic rows, inspecting H', H'',
the inner search path and competing stationary points. This distinguishes
failure to find a maximum from a degenerate or materially multimodal integrand.

For a demonstrated search failure, the smallest justified repair is confined
to locating and certifying the inner mode: retain the exact conditional
likelihood, finite outer parameters, prescribed outer starts and final fit
criteria; use a safeguarded scalar search and verify residual score and
negative curvature before applying the existing adaptive quadrature. Multiple
modes would additionally require demonstrating integration accuracy, rather
than assuming that one mode determines the whole integral. Any repair must
be checked at the actual failing parameters against higher-node and independent
scalar integration, including gradients used by the outer fit.

Do not clip curvature to negative, use its absolute value, silently return
v=0, cap the statistical variance, or replace an exception by a guessed
likelihood. No alternate quadrature fallback or new failure policy is
authorized by this source-only review. v1.7 leaves numerical polish to each
implementation but retains prescribed starts, likelihood agreement and
failure/exceedance rules.

## Fixed-call evidence from concurrent diagnostic21998880

The concurrent worker's bounded integrand check is running with the exact
recorded synthetic parameters and reconstructed original draw. Its completed
small `damping_and_ascent_probes.json` now records at u=0:
H=−93.29258296470611, score=−103.39494676363267 and
curvature=763.1283114221253. Stable log-shape arithmetic agrees locally:
H=−93.29258296462017, score=−103.39494676242259,
curvature=763.1283114064892. The scalar likelihood and model evaluator agree
at u=0 and the small probes to displayed numerical tolerance.

The original curvature-clamped proposal is−10,339,494,676.363268; even its
smallest of30 damping trials is−19.258809604427345. That final trial reduces
H to−421.1639644585133 and is rejected. By contrast u=−0.001 increases the
stable H by0.10375439716746371, and u=−0.0001 increases it by
0.010343287915219435. Thus the recorded u=0 point is nonstationary and the
fixed-call search omits demonstrably valid smaller ascent steps. Local
old/stable agreement makes extreme-tail arithmetic unnecessary to explain
this refusal; it does not certify either implementation everywhere.

The accompanying Gaussian envelope gives radius13.126614802530403 for a
global maximum whose H is at least H(0). All original damping displacements
are larger than that radius. The diagnostic covers−16..16 with finer grids;
its extrema, complete quadrature comparison and tail/error summary are still
pending. Do not infer a unique mode, integration accuracy or calibrated test
from the small-probe evidence. No production source or failed original draw
has been replaced.

The concurrent check subsequently **completed21998880,9:39,exit0**,
peakRSS57,976KiB. Both0.01 and0.005 grids have five identical sign-change
brackets. Root inspected its saved complete summary: three verified stationary
maxima were located at u=−4.3420208251883885,−0.1597680479632258 and
4.2315863428009655. The middle maximum has H=−84.00014325137812;
the other two have H=−206.924067868171 and−176.94233713770194.
The grid can miss tangential/subgrid stationary points; these are the located
maxima, not a proof that no other stationary point exists.

Direct bounded integration gives logI=−86.38066008342989. The20-node adaptive
quadrature centered at the dominant located maximum gives−86.38066010426598;
root independently subtracted these saved values, obtaining−2.083609729197633e−8.
The quadrature's reported relative error estimate is5.306317872893949e−13,
with no warning. The analytic Gaussian tail upper bound relative to the
bounded integral is2.3970013727831635e−24, rederived from its saved log bound.
Centering on either lower maximum instead gives a materially smaller integral.
The full checked grid has6,401 finite original evaluations, with old/stable
maximum differences1.4500756151392125e−10 for H,3.0807711937086424e−9 for
score and4.134244591114111e−8 for curvature.

This supports locating a valid dominant mode before using the existing formula
at this synthetic call. It does not test the marginal parameter gradients
used by the outer optimizer or other proposals. Root assigned a separate
bounded normalized-prior posterior-score/finite-difference check for those
gradients. No alternative production integration, test law, start, failure
policy or numerical-source change has been made.
