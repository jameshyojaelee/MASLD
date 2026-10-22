# Fixed synthetic gene2 marginal-gradient check

Prepared 2026-09-30; numerical job21998909 **completed, exit0**. The measured
20-node score discrepancy and numerical limits are reported below; completion
does not mean all agreement checks passed. This standalone diagnostic
uses the same frozen draw0/gene2/parameter point as completed mode check
21998880. It does not change the original mode search, outer starts, likelihood
law, test statistic, production source or treatment of failed draws. No receiving
or private data are inputs. The check is a separate numerical investigation,
not an integration fallback for a fitter.

Files:

- `scripts/analysis/in_biopsy_allelic/rec_codex/codex_p3_gene2_gradient_diagnostic.py`
- `scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_p3_gene2_gradient_diagnostic.sbatch`

The launcher requests one cpu/nslab job, one CPU, 2GB and one hour; label 1598,
job-specific `1598_%j.out/.err`, no array. Root reviews/submits separately. The
output prefix is `codex_p3_gene2_gradient_<jobid>`, distinct from the previous
mode diagnostic and RNA index despite reuse of the numeric launcher label.

## Exact target and source evidence

Let L_g(p,u) be the product of the original properly selection-conditioned row
probabilities for this synthetic gene, and let u have normalized N(0,v) density.
The quantity checked is

\[
I_g(p,v)=\int L_g(p,u)(2\pi v)^{-1/2}e^{-u^2/(2v)}\,du.
\]

For v>0, differentiation under this finite-score Gaussian integral gives

\[
\partial_p\log I_g=E_{u\mid g,p,v}[\partial_p\log L_g],\qquad
\partial_v\log I_g=E_{u\mid g,p,v}[(u^2-v)/(2v^2)].
\]

The variance score includes the derivative of the prior normalizer. Omitting it
would add 1/(2v) to this derivative. The archived
`bfix_p3_21998519/executed_sources/likelihood_v1.py:425–444` explicitly uses this
normalized prior, its variance score and posterior-weighted parameter scores.
Its predictor Jacobians and conditional gradient are reused from the immutable
source, not independently invented (`likelihood_v1.py:213–242`). The comparison
is conditional on that row likelihood. Finite differences provide a separate
check of its resulting integrated gradient; they do not independently verify
the biological law.

The completed mode-check summary contains three sign-bracketed maxima. The
central maximum u=-0.1597680479632258 has unnormalized log integrand
-84.00014325137812. Its bounded log integral is -86.38066008342989; the central
20-node approximation differs by about 2.08e-8. Agreement of masses alone does
not establish gradient agreement. These numbers were read from the completed
synthetic summary; no new numerical calculation of them was performed here.

The new helper guards the original failed-call receipt/summary, original role
receipt, original generation parameters, all executed-source hashes, the
synthetic fixture and uniform draw hashes. It additionally guards completed
21998880 summary/protocol/generation hashes and that job's archived helper and
stable beta-binomial arithmetic source. It reconstructs the 173 synthetic rows
with the original generator before using stable arithmetic, then checks the
completed reconstruction hash
`5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691`.
Only the 24 synthetic rows for gene2 enter these integrals. The original
Python/NumPy/SciPy versions must match; all imported model sources are private
copies of immutable executed files.

## Numerical comparisons

Direct vector integration evaluates normalized mass and every p score plus the
variance score over [-16,16], with breakpoints at the recorded extrema and
several distances around the dominant peak. It repeats at tighter requested
precision and over [-18,18]. Scalar mass is checked against the completed
mode-check result after adding the correct prior normalizer. Comparisons save
the requested tolerances, quadrature error estimates, warnings, evaluation
counts, phase timing and changes after precision/domain refinement.

Dominant-center 20- and 40-node Hermite quadrature then provides posterior-score
averages on every coordinate. The checker saves their mass and gradient
differences from the direct integral, without inventing a statistical margin or
declaring predictive success. Finite differences use four absolute step sizes
(1e-3,3e-4,1e-4,3e-5) for alpha_gene2, r_gene2, kappa, omega_S, rho_S, v and the
remaining p coordinate with the largest direct absolute score. That final
selection is a numerical sensitivity diagnostic, not selection of biological
effects. Both direct integrals and 20/40-node approximations are differenced.
Each perturbed approximation recomputes its central stationary maximum and
curvature in a fixed bracket; a lost bracket causes explicit failure. This
also measures whether the supplied posterior-score average equals the
derivative of the finite-node approximation with moving center/scale.

## Tail and numerical uncertainty

Because each conditional row probability is at most one, the omitted absolute
mass outside |u|<=B is bounded by T=2*Phi(-B/sqrt(v)). This remains valid even
with secondary peaks. A finite domain alone does not prove a score integral is
accurate, so the helper separately bounds score tails.

For integer beta-binomial counts, each log-shape derivative is a sum of terms
shape/(shape+j), in [0,n]. The eta and linear-dispersion BB scores are therefore
bounded in absolute value by n: the dispersion transformation has
|d log((1-rho)/rho)/d linear_r|<=1. Mixture weighting preserves these bounds;
subtracting the selection-normalizer score doubles them. Hence, for each p
coordinate j,

\[
|\partial_{p_j}\log L_g|\le C_j
=\sum_{r\in g}2n_r(|J_{\eta,rj}|+|J_{r,rj}|).
\]

These Jacobians do not depend on u for this source's predictors. The omitted
absolute p-score numerator is at most C_j*T. For the variance score, with
a=B/sqrt(v), the omitted absolute numerator is at most

\[
\int_{|u|>B}\phi_v(u)(u^2+v)/(2v^2)du
=[a\phi(a)+2\Phi(-a)]/v.
\]

The helper saves these analytic absolute tail bounds and resulting score
uncertainty after normalization, using a mass reduced by the quadrature error
estimate. These normalized bounds are conditional on that numerical error
estimate: SciPy's error estimate is not a certified bound. Componentwise ratio
errors and finite-difference mass-error amplification by 1/(2h) are saved;
step refinement measures truncation/roundoff stability rather than presuming
an exact derivative. Arithmetic error outside those integration estimates
remains an explicit limitation.

## Completion state

Root reviewed the helper and submitted **21998909/1598**, cpu/nslab,
one CPU,2GB,one hour, no array, with1598_%j logs. It subsequently completed
13:24,exit0,empty stderr,peakRSS55,408KiB. Complete results are below.
Submitted helper SHA256 is
`59d2574a804d7d3d3bf9ba6699b815c5684b8b045b81402949d1d23142b95418`.
No production numerical source, outer start, convergence threshold or original
failed draw was changed. This check remains confined to the frozen gene2 call;
it does not additionally test P1draw24 or P2draw3.

AST parsing and `bash -n` passed. Installed SciPy `_quad_vec.py` documents the
used max-norm error estimate and full-output fields. Source inspection found no
other P3 gradient helper in the bounded Codex directory at preparation. Output
is limited to 16MiB; executed-source/input bytes and actual saved-output bytes
are recorded. At preparation no numerical check had run; the subsequently
submitted check's completed result is detailed below. No install, calibration or outer fit
has been performed. Any future result applies only to this fixed
gene/draw/point. Passing it would support a subsequent bounded review of a
numerical repair; it would not establish all-gene integration, convergence or
test calibration.

## Completed fixed-call result

Job21998909 completed all seven selected directions and24,589 likelihood
evaluations. At this same synthetic call, direct normalized logI is
−87.25979655580508. Root rederived the maximum weighted-score discrepancy:
20nodes1.279247156843466e−6 at omega_S(index24);40nodes4.2518877307884395e−11.
Tighter direct integration and domain expansion±16→±18 give exactly the same
saved scores. The existing fixed-vector numerical comparison rule, used here
only as a diagnostic, is exceeded by kappa and omega_S; this does not change
the already completed registered theta1–theta4 comparison.

The completed omega_S finite differences also distinguish the weighted-score
approximation from the derivative of the finite-node objective. At h=0.001,
direct logI differentiation differs from the direct posterior score by
−8.895200437342154e−9. Recentered20-node logI differentiation differs from that
direct score by+1.9121545735689693e−9, but from the supplied20-node weighted
score by+1.281159311417035e−6. At h=0.0003 that latter difference remains
+1.2568255912115411e−6. These measured differences support a jacobian/objective
precision issue at the fixed call, while40-node score differentiation is much
closer at these steps. Smaller finite-difference steps can show more rounding
noise; the smallest step is not treated as automatically decisive.
All directions and the final summary are now saved in
`codex_p3_gene2_gradient_21998909/summary.json` under the REC execution.
The direct v score is−0.5256462990490359; finite differences approach it with
errors−1.931650576159072e−7,−1.737441635007997e−8,
−1.9082694313610205e−9 and−5.345535125655942e−10 across the four declared
steps. Alpha2 and delta0 also show truncation-error refinement. The r2 and
rho_S differences become less stable at small steps; their finite differences
do not provide a uniformly precise check at every step. No result is selected
as a new scientific success margin.

The largest direct component quadrature error estimate is1.0016737160419318e−7;
the largest analytic score-tail bound, conditional on the estimated bounded
mass, is4.3716775209054414e−20. Numerical estimates do not include a rigorous
bound on likelihood arithmetic rounding. Forty-node weighted scores reproduce
the direct integral scores closely at this fixed call; this is not authorization
to replace the prescribed20-node implementation. Exit0 establishes completed
diagnostic computation, not that every gradient comparison passed or that the
model is calibrated. No final optimum, all-gene calculation or production
repair was evaluated. Submitted source hash remained unchanged at inspection.

## Independent review of private damping replay 21998950

This bounded review inspected source and small saved synthetic diagnostic
receipts only. It did not rerun or duplicate the replay, read restricted or
receiving values, change any owner's source, submit work or perform a fit.
The inspected helper
`codex_p3_private_damping64_replay.py` has SHA256
`4c0d22e8e05efae7d618123f419947043172e5a8e0e58f5b322327c4338c4315`;
its launcher has SHA256
`9ab124bed574039d49360f148cb2afb83ed421620fdaba2ae3c3594ec614cb69`.
The launcher requests one cpu/nslab job, one CPU, 2GB and 90h, without an array.
This review makes no new live scheduler-status claim.

**The text change is exactly one inner damping loop.** The helper requires one
occurrence of `np.arange(30)`, replaces it with `np.arange(64)` and verifies
reversibility (`codex_p3_private_damping64_replay.py:96–106`). The saved
`codex_p3_private_damping64_21998950/private_source.diff` changes only archived
likelihood line417. The private likelihood hash is
`8880a8059ca18259a43b638ce9b0c60991ab7b555ec3db736d31f85d548f8ff8`.
Relative to the complete historical computation there are two numerical
changes: this private damping-loop extension and the separately checked stable
BB representation injected at runtime (helper lines148–155). Neither changes
the probability law. Describing the whole computation as solely a damping
change would omit that arithmetic substitution.

**The replay retains the original outer recipe.** The immutable replay driver
is loaded with only its source location/hash and output parent redirected to
the private copy (helper lines203–215). The private likelihood still defaults
to 20 quadrature nodes (`likelihood_v1.py:376`). Its warm start is the null p with
v=0; its cold start uses v=tau_max², with null p for this
`section4_from_null=True` draw. Both starts, shifts 0/−0.1/+0.1, bounds, boundary
restart and the retained-start rule are unchanged: cold is retained only when
its loglik exceeds warm by at least 1e-6, otherwise warm is retained (private
likelihood lines452–525). The driver still computes zero LR at v=0 and otherwise
max(0,2*(marginal−null)); declared fitting failures retain NaN LR and count as
exceedances, with the original 1e-6 comparison convention (driver lines191–218;
original `check_bfix_tests_v13.py:84–106`). A successful numerical replay does
not remove the historically recorded failed draw or establish calibration.

**The precheck proves fixed-call agreement, with limited scope.** The receipt
`fixed_gene_precheck.json` has SHA256
`85bb770af3bf814ff64b7035cff9bc106b47f165fcbbbe440310ea161159a1e8`.
Its 173-row reconstruction hash matches the completed mode and gradient checks.
The extracted private mode-search AST, with 64 proposals, reaches
u=-0.15976804796323024, score=-3.537614645665599e-12 and negative curvature.
The 20-node mass differs from the bounded reference by 2.0834448832829366e-8
in log units. It passes the helper's stated 1e-6 mode/mass checks. This is one
local-search result agreeing with the independently dominant recorded maximum
at the failed call. It is neither a global mode guarantee for new parameters
nor a gradient check. The full replay generates again from a refitted null;
its realized-count hash is checked only after the nested driver returns
(helper lines216–220). At this review's inspection neither that nested
`generation.json` nor the outer completion summary existed. [Unverified] The
live replay's actual regenerated draw identity remains pending even though
its separate precheck reconstruction agrees.

**Completed gradient evidence limits the mass result.** This review separately
read 21998909's final summary, SHA256
`f501cc10c9b65d8db20d56f2077a063aa0fe5c2b4bb8b83daf73878e0786b833`,
and its finite-difference entries. The same frozen call's 20-node weighted
scores differ from direct normalized integration by +4.827439008536771e-7
for kappa and −1.279247156843466e-6 for omega_S. At the latter coordinate the
direct component quadrature error estimate is 4.1365644209711126e-10 and its
analytic score-tail allowance is 1.614517800388724e-20, conditional on that
error estimate. Domain and requested precision refinements leave the saved
scores unchanged. The h=0.001/0.0003 differences between the recentered
20-node objective derivative and its weighted score remain approximately
1.28e-6/1.26e-6. Forty-node weights agree closely with direct scores at this
point; this is a diagnostic comparison, not authorization to change 20 nodes.
[Inference] A discrepancy comparable to the existing 1e-6 outer convergence
threshold warrants checking a completed fit before calling its gradient
stationary. It does not prove that this nonoptimal proposal caused a final
convergence failure or that a completed fit has inaccurate parameters.

**Concrete next check after a fit completes:** hold its saved parameters,
realized draw and retained start fixed. For an interior v>0, evaluate each
gene's normalized direct integral and posterior scores with a domain/tail
bound recalculated at that fitted point; sum over all genes and compare with
the unchanged 20-node objective and gradient. Include finite differences in
kappa, omega_S and v, with step refinement and recentering for the 20-node
comparison. Reapply the original bounds, fixed uninformative coordinates and
1e-10 active-bound convention when computing the projected gradient. Report
whether the direct-gradient discrepancy could change the original 1e-6
convergence conclusion, rather than fitting again or selecting a new margin.
Retain the saved start-specific objectives and original tie rule. If the
retained fit has v=0, check the original exact null gradient and right variance
score, sum_g 0.5*(s_u,g²+h_u,g), with its nonpositive-boundary criterion instead
of a positive-v quadrature calculation. This checks a completed point without
changing nodes, starts, failure handling or the fitting law. No such
completed-fit calculation was performed by this review.

Log each gene's selected mode/curvature at the completed point and every
finite-difference perturbation. Use one-sided differences where required by
the original bounds. If mode selection switches branches, that difference is
a secant across the selection change and does not establish a local derivative.
Retain the projected-bound and variance-boundary distinctions. Numerical error
estimates cannot resolve an agreement smaller than their own precision, and
abnormal termination alone does not prove the saved parameters are wrong.
Counting failed draws as exceedances preserves the declared disposition;
it does not demonstrate that the resulting test is calibrated.

## Completed retained-point check 21999267

**The actual retained point has v=0.** The completed private replay 21998950
retained the warm start with LR=0. The independent, nonfitting checker 21999267
held that point, realized synthetic draw, original bounds and 20-node source
fixed. It checked all 29 informative coordinates in 233 logged objective calls.
Its baseline loglik is −511.5784700493092; logged and unlogged likelihood values
and gradients are bitwise equal, and both reproduce the saved fit exactly.
The logging-only statement strips back to the original marginal AST. At the
baseline the mode log is empty: v=0 uses the conditional likelihood without
inner quadrature. The earlier positive-v gene2 discrepancy therefore is not
evidence of a quadrature-gradient discrepancy at this retained point.

The analytic projected maximum is 9.223438679484275e-11, below the original
1e-6 criterion. In the checker's zero-based labels, r_2 and r_4 are at their
−15 lower bounds; their negative scores point outside the feasible region
and project to zero. The exact right variance score is
−541.6573914760489. Independently summing the conditional expression
`sum_g 0.5*(s_u,g**2+h_u,g)` reproduces that score exactly. This supports the
local analytic boundary condition at this point. It does not establish a
global maximum, accuracy of other draws, or test calibration. The variance
bound remains the archived source's `[0, None]`; tau_max² is a starting value,
not an upper bound (`likelihood_v1.py:461,468,481`).

**Finite differences do not independently certify the entire 1e-6 condition.**
All four declared absolute steps were evaluated, using central differences
for interior coordinates and forward three-point differences at the bounds.
The maximum projected finite-difference magnitudes over the full checked
vector were:

| Absolute step | Maximum magnitude | Coordinate |
| --- | ---: | --- |
| 0.001 | 1.4505954482046945e-4 | alpha_2 |
| 0.0003 | 1.2996205593177972e-5 | alpha_2 |
| 0.0001 | 1.751345735101495e-6 | alpha_2 |
| 0.00003 | 2.537111261820731e-6 | rho_S |

For alpha_2 the successive central estimates approach zero, with a final
estimate of 5.712763595511205e-7. The larger-step behavior is consistent with
second-order truncation. For rho_S the four estimates are
5.113918177812593e-7, −2.00752007610087e-7, 1.5802470443304628e-7 and
−2.537111261820731e-6, while its analytic score is −1.2670073323839404e-13.
[Inference] The sign changes and amplified finest-step discrepancy are
consistent with arithmetic sensitivity during subtraction; they do not
identify a nonzero true derivative. The saved last-adjacent refinement plus
subtraction-scale estimate for rho_S is 2.700058339955641e-6. This is a
heuristic estimate, not a certified likelihood arithmetic error bound.

Bounded postprocessing of the saved central differences, without further
likelihood evaluations, used every adjacent step pair and all 26 central
coordinates. Under the explicitly unverified expansion
`D(h)=D0+C*h²+O(h⁴)`, unequal-step Richardson gives
`R(a,b)=(a²*D(b)-b²*D(a))/(a²-b²)` for coarse a and fine b. Maximum absolute
residuals against the supplied analytic scores are:

| Adjacent steps | Maximum absolute residual | Coordinate |
| --- | ---: | --- |
| 0.001 → 0.0003 | 2.9649191910021485e-7 | rho_d |
| 0.0003 → 0.0001 | 6.440009512956313e-7 | rho_d |
| 0.0001 → 0.00003 | 2.8036630438703715e-6 | rho_S |

These use all pairs, rather than selecting a pair that satisfies 1e-6.
[Inference] Removing the leading coarse-step truncation term improves those
estimates, but extrapolation amplifies the small-step sensitivity. Neither
Richardson nor the saved heuristic uncertainties provide a certified
whole-vector error bound. The r_2, r_4 and v forward differences are excluded
from this central-difference extrapolation.

The four secondary positive-v probes give right derivative estimates
−516.3867414498782, −537.5004552263363, −541.0826815136716 and
−541.6012084045482. They approach the exact negative boundary score and have
the same sign. They do not replace the exact v=0 calculation. Their mode and
curvature logs are saved per evaluation; negative curvature alone is not
proof of exact mode stationarity or globally correct integration.

**Execution and immutable evidence.** Root recorded COMPLETED, exit0, 2:46,
RSS55,904KiB for 21999267; the checker's measured elapsed time is
165.4884454011917 seconds. Logs are `REC/1590_21999267.out` and
`REC/1590_21999267.err`; stderr is empty. The results are under
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_p3_retained_fit_precision_21999267/`.
All protocol-recorded frozen input hashes match the inspected bytes. The
realized synthetic count-vector hash is
`5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691`.
The final regular-file byte sum is 750,507, below the 33,554,432-byte cap;
the saved 750,311-byte progress figure precedes its own 196-byte final line.

| File | SHA256 |
| --- | --- |
| `summary.json` | `011561653fa3f3b7361987cb6cb46bbf36fede38426cdd048f34945d65eb7601` |
| `protocol.json` | `0e985d571db5841de3aa3ffb8bb94efd2645f68aa8ae4029adf54739d318615b` |
| `baseline.json` | `1076ee85652a8997ef27185907f1de3822c53800e413a49e711c49f31f097567` |
| `zero_variance_boundary.json` | `89c78dbe66c4e0d367995bc6ec830e194e230b555096246794653a219fafc70a` |
| Executed `codex_p3_retained_fit_precision.py` | `4c81136cfc93a0b8544ab16f2d51acde4b4c3963483cdad0bf4d35fd6e5d2548` |
| `run_codex_p3_retained_fit_precision.sbatch` | `d8f8561afe5c57f5e28a285161597a4f64ec13f7da823aded0de65e4c5de0e6a` |

This review only reads the completed synthetic results and postprocesses saved
differences. No objective re-evaluation, refit, true-integral calculation,
node change or statistical disposition change was performed. The analytic
boundary-point reproduction is supported; independent whole-vector finite-
difference certification is not established at 1e-6 with these results.
