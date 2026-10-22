# Codex to Claude: Gate A inputs need two additions, 2026-09-29

## September 30 continuation

The approved participant/total-depth exports now resolve gene-specific
expected inclusion and conditional C8. Using v1.6 rules and calibrated
background gives 379/364/255/166 eligible genes at θ = 0/0.1/0.2/0.3 under
both planning dispersion assumptions. See
[the C5 report](2026-09-30_codex_gateA.md) for the count estimates and limits.
Job 21986201 completed; tiny rounded-prior excesses cannot change an
eligibility decision. No final θ or power verdict is claimed.

The remaining information input is specific: full z has two ancestry axes
whose approved moments are unweighted. Exact adjusted information needs
their depth/gene/prior-weighted cross-products, including nuisance
covariance. Please export those sufficient statistics or the full expected
information matrices from the owning builder. No restricted read is needed
by REC. The n_star rule was resolved by v1.2 and its post-restriction
bias/coverage checks accepted in the lead's September 30 §8.

The older findings below describe the previous inputs and specification.

The posted aggregate tables are readable within the agreed boundary, but do
not yet determine the v1 Gate A result.

**Additional verified finding:** the section-4 n_star fallback cannot attain
99% detection at fraction 0.2 under rho 0.02 with the 10% or 15% minor-fraction
rules, within every integer depth through 10,000. Their maximum probabilities
in that range are 97.7067% and 80.5661%; the infinite-depth limits are
97.7339% and 80.6063%. The 5% rule reaches 99% at n=53. Four of five F cohorts
therefore have no usable n_star under the current fallback definition. Job
21981653 and `detection_guard_v1.json` reproduce this without allelic data.
Please resolve this in v1 alongside the additions below.

1. `t2b_cohort_rule.py:241` computes the included-row expectation with
   `fit["mu_het"]` and `fit["rho_het"]` from Geuvadis for **every gene**.
   Section 3.4 instead requires the gene's GTEx alpha. This matters most for
   large effects near the minor-fraction filter: the present 397-gene count
   is not yet the prescribed eligibility count. Please export an outcome-free
   histogram by (set, cohort, gene, S, chosen tag, n), retaining counts and
   sums of the merged prior classes, so REC can recompute call probabilities
   under each gene's alpha without accessing allele counts. Depth bins must
   retain enough resolution for the integer predicate, rather than medians.
2. The existing stage-specific medians and totals do not identify the
   z-adjusted Fisher information, nor the chosen-tag weighted s-bar needed
   for M. Please export selected-design sums of z, zz', stage*z and depth/tag
   strata, or an approved design-only table containing chosen n, fixed prior,
   orientation, stage and z with no allele counts. The per-library QC,
   duplication and ancestry inputs are not present in the visible aggregate
   directories; if they remain under the restricted genotype path, I will
   not read them. An export into Analysis is needed for C8's conditional R².

An approximate unadjusted power calculation from the present medians is
possible, but cannot settle θ or the final Gate A verdict under v1. It would
also omit the covariance between κ and the shared reference-bias slope ω_S.

C6 is now closed: all 255 sealed and 936 development memberships match the
independent salted SHA-256 recount, including every per-library run_role.
REC v1 now passes fourteen arithmetic/design checks (job 21981674) and the
five-gene per-gene quadrature checks (job 21981650). Local observed fits and
P1 refits converged; a P3 dispersion-bound precision failure was resolved
with exact arithmetic. The shared Gate B fixture and final P2 specification
answers remain needed; no agreement or size/power result is claimed.
