# Supplied toy draws: P1 24/25 and P2 3

This is a source and small-JSON review of the synthetic agreement fixture,
not scientific calibration. Draw numbers are zero-based. No fits, gradients,
Hessians, arrays, or numerical tests were recomputed; no source or job changed.

## Recorded conditions

The fixture has eight genes and seven design columns. Parameter indices below
are zero-based: alpha 0–7, r 8–15, kappa 16, delta 17–23, omega_S 24,
rho_S 25, rho_b 26, rho_d 27. The saved P2 source defines bounds and
predictors at `executed_sources/likelihood_v1.py:203–218`.

|Saved draw|Recorded warm/cold convergence|Kept start|Largest reported interior score relevant to failure|Recorded statistic|
|---|---|---|---|---|
|Original P1 24|false / true|warm|omega_S, index 24: +5.749339058e-6; kappa, index 16: −5.139877108e-6|z=NaN; exceed=true|
|Original P1 25|true / true|warm|All reported interior scores below 1e-6|z=−0.46427268786007; exceed=false|
|Isolated P2 3|false / true|warm|alpha(G8), index 7: −3.214674488e-6; r(G8), index 15: −2.643592057e-6|kappa=NaN; supplied lambda pair retained|

P1 24 places r(G3) and r(G5), indices 10/12, at −15 with negative
scores −0.001068998423/−0.000052075416. These point outside the permitted
range and project to zero. P1 25 has the same valid lower-bound pattern
and converges; a large raw bound score does not itself imply a failed fit.

P2 3 places r(G2)/r(G5), indices 9/12, at −15 with negative scores.
It also places rho_b/rho_d at +3 with positive scores
1.728057470/0.187397017. All four are outward scores at active bounds
and correctly project to zero. Its remaining interior scores exceed the
unchanged strict `max(abs(projected_score)) < 1e-6` criterion. The large
raw score 1.728 is therefore not the convergence defect.

These outputs do not record an optimizer termination message or a polish
history. [Unverified] They cannot establish whether iteration limits,
line-search rejection, conditioning, or arithmetic caused the residual
interior scores.

## Starts, source versions, and failure policy

Warm/cold likelihoods are respectively:

- P1 24: −533.01604200377 / −533.0160420037774.
- P1 25: −530.720741732697 / −530.7207417326919.
- P2 3: −500.1941068365454 / −500.19410683674164.

None supplies the >=1e-6 cold improvement required by the prescribed
two-start selection rule. Selecting cold because it converged would change
that rule. P2's saved source independently confirms both starts and selection
at lines 259–270, projected convergence at 305–315, and the existing constrained
polish call at 308–310. P2 3 already used that source; it is not evidence of
an omitted polish call.

Original P1 runs in `bfix_tests_21985999/`, located through `1538.out`,
not a `bfix_p1*` directory. The log confirms draws 24 and 25 were saved.
That process loaded its numerical source before the later constrained polish
and arithmetic revisions, as documented in `2026-09-30_codex_bfix_results.md`.
It has no executed numerical-source snapshot in its output directory.
Current source bytes or the P2 snapshot hash must not be attributed to P1.

Failure handling is retained in the inspected records: P1 24 remains a NaN
statistic counted as an exceedance; P1 25 remains its genuine finite statistic.
P2 3 retains its slot and lambda pair with NaN kappa. Its executed runner
lines 63–70 append every draw, and saved `tests_v1.py:70–84` includes invalid
kappa in both directional exceedance calculations without reducing the
draw denominator. A final P2 result was not present at inspection, so final
aggregation is not yet verified.

## Actionable check within existing work

P2 3's interior alpha(G8)/r(G8) coordinates and active rho_b/rho_d bounds
are useful distinct cases for a subsequent review of null-fit termination
if it is warranted. The separately submitted1598 checker remains confined
to P3's frozen gene2 marginal call and is not broadened to these different
rows/fits. A raw maximum over all coordinates would give the wrong diagnosis
here. This review does not create a duplicate gradient checker or arithmetic
replay.

Any subsequent correction must preserve both prescribed starts, their
likelihood comparison, original bounds and convergence threshold, every
supplied draw/paired reliability value, and the original files. The existing
P1 draw-9 two-start correction procedure is the relevant precedent; no new
replacement draw or convergence-based start selection is justified.

## Inspected receipts

Paths below are relative to
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/`.

|File|SHA256|
|---|---|
|`bfix_tests_21985999/p1_024.json`|`7a6ea53579890bdcc2a9d139ada3fa7b7ff0fc7d10c612d4da0db97ec5f9f79f`|
|`bfix_tests_21985999/p1_025.json`|`01c43d8dea586e144bcdae017c7bff4d06eae5b2cab06759f6d48b833826ca30`|
|`bfix_p2_21998518/p2_003.json`|`bea38141306a6ddc45ed20cc537501767fb0461f5803fb6d80a47c6d72aa660c`|
|`bfix_p2_21998518/executed_sources/likelihood_v1.py`|`41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6`|

Checks performed: bounded source/log/JSON inspection, receipt SHA256 checks,
and search for an existing draw-24/draw-25/P2-draw-3 diagnosis. Numerical
verification and scientific calibration: not run.

## P2 draw3 reconstruction and focused check prepared

The actual P2 executed source can reconstruct this conditional draw without
fitting. Its `bfix_v1.py:19–45` loads the supplied 60-person/173-row toy design;
`design_v1.py:128–160` expands the prescribed cohort-cluster resample and
recreates complete-participant centering with the original schema. The
guarded `draws/p2_resample.tsv` contains draw3's 58 cluster slots: 19 in C1,
19 in C2 and 20 in C3. These slots are not the 60 original toy individuals;
the expanded participant and row counts will be recorded on compute. No
original expanded draw-row hash was saved, so a new reconstruction hash
cannot independently prove agreement with an archived realized-row hash.

`p2_003.json` stores the complete selected warm parameter point, score,
Hessian and both start flags/likelihoods. The same role's `observed_null.json`
stores its warm initializer, SHA256
`a617d331fd3fb93576e5bd8588581aaab8f77946936536401905490438660181`.
The role receipt, SHA256
`63b377f65b4f37c8b6d60ffca2745f639ffa44519a7bdd2a0a55b682402f8cd8`,
binds all required supplied design/bootstrap inputs and source snapshots.
There is no missing input preventing reconstruction. No no-information
coordinate is assumed: the archived `informative()` mask and initializer
determine the effective fixed bounds on compute.

Prepared, **not run**:
`scripts/analysis/in_biopsy_allelic/rec_codex/check_p2_draw3_fixed_point.py` and
`run_check_p2_draw3_fixed_point.sbatch`. The helper compares saved likelihood,
score and Hessian against the archived original source and an isolated copy
using the already guarded log-shape BB probe
`codex_p3_gene2_integrand_21998880/executed_sources/codex_stable_bb_probe.py`,
SHA256 `723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488`.
No live numerical module or source is changed.

The only finite-difference coordinates are alpha(G8)7, r(G8)15, kappa16 and
the lower-bound r(G2)9/r(G5)12 controls, at absolute steps 1e-3, 3e-4, 1e-4,
3e-5. Five-point likelihood derivatives check scores; derivatives of scores
check full Hessian columns; likelihood second derivatives additionally check
their diagonals. Bound coordinates use inward one-sided stencils. This is a
focused conditional derivative check, not proof of whole-fit stationarity or
scientific calibration. Step instability remains a possible inconclusive
outcome. No fit, cold selection, draw replacement, or lambda-slot change is
performed. The original failed point and supplied lambda pair remain intact.

Launcher 1580 requests cpu/nslab, 1 CPU, 2 GB, 1 hour, actual-job-ID log and
output names, and a 32 MiB own output allowance checked after writing the
summary and source inventory. [Estimate] Memory is well below 2 GB for this
28-coordinate conditional point; the generous hour bounds the focused
derivative evaluations. AST parsing and launcher `bash -n` passed. Executable
numerical verification and submission: not run by this worker.

## Completed fixed-point check: retained P2 draw3 still fails

Job **21999265 / 1580** completed the focused check in 42 seconds, exit 0,
empty stderr. No fit, original draw replacement, cold selection, lambda-slot
change, or scientific calibration was performed.

|Job|Progress|State|Errors / memory|Downstream interpretation|
|---|---|---|---|---|
|21999265 fixed conditional point|1/1 ██████████|Completed; 0 running, 0 pending, 1 done|Exit 0; stderr 0 bytes; batch MaxRSS 54,436 KiB|Keep original failed P2 draw3 slot|

Reconstruction expands the 58 source bootstrap cluster slots to **59 sampled
toy individuals and 172 rows**, from the original 60-person/173-row fixture.
There are no no-information coordinates in this resample. Saved versus
recomputed original likelihood, gradient and Hessian differences are all
**exactly zero**. The maximum projected score remains
**3.2146744879391065e-6**, above the unchanged strict 1e-6 criterion.
Warm/cold flags and likelihoods remain as recorded above; cold does not earn
the required >=1e-6 likelihood improvement.

At this same parameter point, replacing only the BB arithmetic with the
isolated guarded log-shape implementation changes likelihood by
6.821210263296962e-11, any score by at most 1.8170149699626471e-9, and any
Hessian entry by at most 1.1142731182189891e-9. The two failing interior
coordinates remain above tolerance under that representation:

|Coordinate|Stable analytic score|Likelihood FD, h=1e-3|Likelihood FD, h=3e-4|
|---|---:|---:|---:|
|alpha(G8), index 7|−3.214491723e-6|−3.211335790e-6|−3.212884773e-6|
|r(G8), index 15|−2.643491297e-6|−2.638680977e-6|−2.641182088e-6|
|kappa, index 16|+6.099480725e-7|+6.236007266e-7|+5.272700794e-7|

Thus BB derivative precision does not explain the retained warm point's
failure at alpha(G8)/r(G8). This is a focused nonstationarity finding, not
proof about the entire fit or the reason its optimizer stopped. The kappa
score itself is below tolerance; the draw remains failed because the selected
full conditional fit did not meet its prescribed criterion. Its supplied
paired reliability values are preserved in the original failed slot.

Finite differences are not uniformly precise. Stable score-derived Hessian
columns for indices 7/15/16 at h=1e-3 differ by at most roughly 1.43e-9,
3.66e-10 and 1.12e-9 respectively, but direct likelihood second derivatives
are noisier and vary by coordinate/step. For example, stable kappa's diagonal
discrepancy is +0.000104 at h=1e-3 and +0.257841 at h=3e-5; the lower-bound
r(G2) control reaches −0.425945 at the smallest step. Small-step objective
score estimates also become noisy. [Inference] This pattern is consistent
with subtraction cancellation amplified by division by h or h². No arbitrary-
precision reference was run, so no blanket finite-difference Hessian pass is
claimed. The stable-versus-original and larger-step score evidence above
remains separate from that limitation.

All 13 recorded input/source hashes still match; copied original/stable
sources and the executed helper/launcher match their recorded hashes. Whole
output size including summary and seven source copies is **166,082 bytes**,
below the **32 MiB** allowance. The reconstruction row hash is
`c90d901a90c8655e747049e3b7fd4fc2e7291d34a0b9e9af3ffc6b1ab742c83b`;
the individual/gene/depth/count tuple hash is
`c24bb9dc0e25bd0d8d72059ecc091338aab20b115bb6d42f65cd54b14537155c`.
An original expanded-row reference hash remains unavailable.

Executed receipt:
[summary.json](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/p2_draw3_fixed_point_21999265/summary.json),
SHA256 `40234628601f4a9a3e2435de059ff35f1462c2841b56fa1178ef4e952e182fb4`.
Only this owned report was amended after the run; no new tests or jobs were
started. Whole-fit stationarity, corrective refitting and calibration are not
claimed.
