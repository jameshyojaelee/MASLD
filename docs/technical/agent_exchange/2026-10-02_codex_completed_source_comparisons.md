# Completed source sequence comparison and original synthetic P2/P3

Checked 2026-10-02 22:02 UTC. This technical report records completed fixed
analyses, without changing manuscript authorities, model recipes, or results.
`REC` below means
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/`.

## Result

The full associated-peak sequence readout has no demonstrated improvement over
the centered readout. Their joint model reduces development MSE by 0.355% relative to the
centered model, with a descriptive chromosome interval spanning zero. These
are previously inspected source development data, not independent validation.
The completed original synthetic P3 calculation retains 24 failed draws out
of 50. Its reported p=1.0 does not establish calibrated inference.

## Scheduler evidence

`sacct -j 21997998,21998518,21998519,22000146
--format=JobID,State,ExitCode,Elapsed,MaxRSS -P` returned four COMPLETED jobs,
each with exit code 0:0. A separate original serial job, 21985999, remains
running. Running 1; pending 0; done 4 in this bounded set.

| Job | Completion | Elapsed | Batch maximum RSS | Numerical errors or limits | Downstream |
|---|---|---|---|---|---|
| 21997998 sequence comparison | 1/1 ██████████ | 00:00:15 | 0 reported; measurement unavailable | No stderr | Fixed comparison complete after sequence scoring |
| 21998518 original P2 | 50/50 ██████████ | 1-02:25:45 | 53,308 KiB | 3 nonconverged fits; SLSQP bounds warning | Original summaries retained |
| 21998519 original P3 | 50/50 ██████████ | 1-07:01:57 | 52,532 KiB | 19 missing fits and 5 nonconverged fits; arithmetic warnings | Exact draw diagnostic complete |
| 22000146 exact P3 draws 12/13 | 2/2 ██████████ | 01:24:41 | 53,780 KiB | Both original exceptions reproduced; arithmetic warnings | Diagnostic only; no replacement fits |
| 21985999 original serial P1→P2→P3 | P1 50/50 ██████████; P2 50/50 ██████████; P3 4/50 █░░░░░░░░░ | Running, 2-07:50:24 at check | 55,220 KiB sampled by sstat | P1 4 failures; P2 7 failures | P3 continues in the original allocation |

Process completion is distinct from fit convergence. P2 stderr is 226 bytes;
P3 and its diagnostic each have 4,800 bytes of divide-by-zero/invalid-value
warnings. An empty error log is not claimed for those jobs.

## Sequence comparison: fixed source development

Results are in `REC/target_readout_development_comparison_21997998/`.
Of 25,477 fixed geometry rows, 24,430 passed sequence admission and 1,047
were excluded. Folds 2–4 supplied 17,989 training rows. Previously inspected
fold 1 supplied 6,441 evaluated rows across chr14, chr19, chr2, chr22 and chr9.
Rows are variant–peak associations, not independent participants.

The target is the unchanged FastQTL signed ALT-dosage accessibility beta,
a marginal association rather than a causal sequence-edit effect. Both
features average the same three native `DIFF_LOG2_SUM` track readouts: one
over the variant-centered 501-bp mask and one over the full associated
BED0 half-open peak. Intercepts and slopes were fitted by fixed unpenalized
OLS on folds 2–4, without tuning or significance filtering.

| Fixed model | Held MSE | Held MAE | Pearson correlation |
|---|---|---|---|
| Zero | 0.42540542451772034 | 0.5866834968793665 | Undefined constant vector |
| Training mean | 0.4247823861218843 | 0.5851056800676669 | Undefined constant vector |
| Center | 0.28714353316732805 | 0.45649953063984394 | 0.5697008478844123 |
| Peak | 0.2877104646663631 | 0.45877482735538494 | 0.5686350219336986 |
| Center + peak | 0.2861236591312102 | 0.45648317012340356 | 0.5718585558006511 |

Positive reference-minus-candidate MSE means lower candidate loss.

| Reference → candidate | MSE reduction | Descriptive 95% interval | Relative MSE reduction | Descriptive 95% interval |
|---|---|---|---|---|
| Center → peak | −0.0005669314990350283 | [−0.0025247369418798385, 0.0003845316836509847] | −0.1974383656777823% | [−0.9144189063596597%, 0.13145972479443602%] |
| Center → center + peak | 0.0010198740361178626 | [−0.0002544662269166231, 0.001639183705994946] | 0.35517917637502006% | [−0.09216355377972274%, 0.5603872139524513%] |
| Peak → center + peak | 0.001586805535152891 | [0.0012533851757185418, 0.0022301374261332724] | 0.5515286129731178% | [0.4237527355585247%, 0.793881631221649%] |

The 2,000 bootstrap draws use seed 20260930 and resample five whole held
chromosomes. Each draw sums chromosome SSE and row counts before taking
row-weighted MSE. Intervals are pointwise across the three contrasts,
conditional on the fixed fitted coefficients and deposited beta; training
is not refitted. Five chromosomes do not establish participant-level or
LD-block replication. No p-values, simultaneous correction, independent
validation, or practical pass criterion is supplied. Fold 0 remained closed.
The joint model beats peak alone descriptively but supplies no clear gain
over the better centered model under this uncertainty summary.

## Original synthetic P2: 47 converged, 3 nonconverged

`REC/bfix_p2_21998518/p2_result.json` reports p=0.7058823529411765.
All 50 original draw receipts exist; none has a missing fit. Nonconverged
zero-based indices are 3, 14 and 16. They remain in the failure accounting.

The two source-specific κ summaries differ in interval validity:

| Source label in fixture | Invalid indices | 90% κ interval | p upper / lower |
|---|---|---|---|
| GSE193066 | 3, 7, 14, 16 | [−0.6075636773431523, 0.2773088051165903] | 0.23529411764705882 / 0.7058823529411765 |
| PXD051911_A1 | 3, 14, 16 | [−0.5255425521180697, 0.22037582261302788] | 0.17647058823529413 / 0.6666666666666666 |

Draw 7 has a converged saved fit and finite saved pooled κ, yet is invalid
for the first source summary. The reason was not diagnosed here. Source
labels are labels in the synthetic fixture, not new evidence from those
human cohorts. Overall pκ is 0.7058823529411765; pτ is
0.2884906619232298. No original failure was replaced by a corrected fit.

## Original synthetic P3: 26 converged, 24 failed

`REC/bfix_p3_21998519/p3_result.json` reports observed LR=0, p=1.0 and
Monte Carlo SE=0.0. All 50 draws are counted as exceedances: the fixed rule
counts unsuccessful fits as exceedances, and finite nonnegative LR values
also exceed the observed zero. This single agreement fixture supplies no
calibration or power estimate.

Missing-fit indices are 0, 2, 3, 6, 12, 13, 15, 16, 17, 22, 23, 24,
27, 29, 31, 35, 37, 40 and 48. Saved but nonconverged fits are 1, 20, 21,
38 and 44. Their LR entries are NaN, and the original missingness is preserved.
Only draws 12 and 13 received the exact exception diagnostic described here;
the causes of other missing fits are unverified.

Both diagnosed draws raise `ValueError: nonconcave mode in adaptive quadrature`
inside the original archived marginal likelihood.

| Zero-based draw | Fit start from traceback | Marginal call | Gene index | Mode | Curvature | Variance v |
|---|---|---|---|---|---|---|
| 12 | Cold | 308 | 4 | 0.07658763958862473 | 595.2263730106787 | 0.0023674063950871395 |
| 13 | Warm | 16 | 2 | 0.0 | 325.17730497248385 | 0.8829298120179231 |

Both failure frames used 20 quadrature nodes and reached the explicit
`curvature >= 0` exception. The diagnostic replay used the exact supplied
uniforms and saved null start, without new draws or model changes. Positive
curvature at these attempted points does not prove that a different optimizer
or likelihood would repair every draw or deliver calibrated inference.

## Separate serial run: P1 and P2 complete, P3 still running

The original job 21985999 writes to `REC/bfix_tests_21985999/` and remains
distinct from the isolated P2/P3 jobs above. Its P1 result has p=
0.23529411764705882, Monte Carlo SE=0.05998846486579747, observed z=
−1.5398358068945364 and observed SE=0.08552062064696732. Original failures
are 2, 9, 24 and 32. Separately corrected diagnostics do not replace them.

Serial P2 reports p=0.7647058823529411 and seven failed fits: 3, 14, 16,
31, 37, 41 and 42. Its first source additionally marks draw 7 invalid.
The respective 90% κ intervals are [−0.6120049491358269,
0.2919391217782978] and [−0.5296521890177605, 0.2260778010688234].
pτ remains 0.2884906619232298. At the scheduler check there were four
saved original serial P3 draw receipts, not a completed P3 summary.

The P2 discrepancy between runs is unresolved. The declared original code
reads supplied `p2_resample.tsv` and `lambda_draws.tsv`, without generating
new P2 draws. The isolated runner preserves the original P2 statement block
by an AST comparison. Saved λ vectors match between the runs for indices
3, 31, 37, 41 and 42. Draw 31 has a serial nonconverged log likelihood of
−409.0627232668295 versus an isolated converged value of
−408.07411268500505. These observations establish different fitted outcomes,
not different realized resamples or a seed explanation. A complete executed
source/input receipt for the original serial process was not located inside
its output directory. No inter-run agreement or calibration claim is made.

## Recorded source identity and checks

The sequence protocol SHA256 is
`f1e1c8451ba4bca0d0b7dc0228c3cb383541b6e900439a6b90f5e528cda0505d`;
comparison executed source is
`4f24c10a7aba2b91377b52aa2773ec1217033f2a871b718607e6804e95a01c17`.
Native scoring source is
`586ecee076917c52098263e3313c170da587fc835ab370141c552e57a9831ca0`.
Original P2/P3 both record likelihood SHA256
`41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6`.
The exact diagnostic source is
`fdfdf55e63225135f06eea45dc8e6dc1189033aab2189cbba03fb8e206e5c63a`.
Recorded synthetic environment: Python 3.10.18, NumPy 2.2.6, SciPy 1.15.2.

| Reviewed result | SHA256 |
|---|---|
| Sequence summary | `66707f6475b08be54c81f9878f4d5bed4f32d26fe0e939c581501c9c9497d28d` |
| Original P2 result | `6d257fd560e0f7301230829a66f94ae6eb810b98ac1a067531c98222ef13c682` |
| Original P3 result | `6231d56149eee71d88f60e28674675a9e10b581a100c05be2e478c3476bafbfe` |
| Exact diagnostic summary | `4b7f5abce51e1e3280bfef371090be0296a9278be083fd1e14fd1e2dad13a2b5` |
| Draw 12 diagnostic | `b2fae226258d93c5cc1c5ac93657f915c9348d61416dad9bf962df81fb0c9ae5` |
| Draw 13 diagnostic | `de45de09a8bcf4ce40aa835c85d0c1140d2c3e32b8063bcc176bc499f7530526` |

Executed checks: scheduler accounting, small summary/contrast/coverage table
inspection, original draw receipt classification, traceback inspection and
SHA256 calculation for the listed small files. Numerical reconstruction of
the sequence comparison and cross-agent all-draw fit agreement were not run
in this read-only review. No jobs or fits were launched. Protected outcomes,
fold 0 and private allelic data were not accessed. Only this technical report
was written.

### Further static check of the serial/isolated P2 difference

Root compared the complete saved observed-null and observed-marginal
parameter lists across21985999 and21998518. Both lists match exactly at every
coordinate, not only at kappa/log likelihood; both retained warm starts.
Current `check_bfix_tests_v13.py`, `design_v1.py`, `bfix_v1.py`, `tests_v1.py`
and `likelihood_v1.py` match the isolated run's archived source bytes.
The likelihood SHA remains41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6.

Static source inspection shows participants, cohort/cluster enumeration and
resampled row construction are sorted. The warm fit copies its parameter
vector before optimization; the driver does not assign a new `null_p` during
P1. These checks do not establish what module bytes the original long-lived
serial process imported. They rule out a different saved observed parameter
list as the explanation; an executed-source/environment or numerical replay
is still needed to identify the cause of the P2 fit disagreement. No fitting,
statistical recipe change, original-fit replacement or new numerical test
was performed for this check.
