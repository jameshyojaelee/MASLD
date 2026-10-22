# Codex review of Model A likelihood v0, 2026-09-29

Reviewed only `2026-09-29_model_A_likelihood_spec_v0.md` and the Model A
coordination note. No allelic outcome, sealed set, or builder implementation was
read. This is a pre-hash critique for Claude's v1; it changes no analysis.

## Decisions needed in v1 before two implementations can agree

1. **The heterogeneity likelihood must integrate a gene-level effect once per
   gene.** In section 4, one `u_g` is shared by all observations of gene `g`.
   Write the gene contribution explicitly as
   `log ∫ Normal(u;0,τ²) × ∏_{i∈g} L_ig(u) du`. Integrating each
   `L_ig` separately would be a different model and would discard the
   within-gene correlation. The truncation denominator in `L_ig(u)` belongs
   *inside* that integral. State whether the 20 Hermite nodes are standard
   physicists' or probabilists' nodes and how they are scaled; otherwise τ
   estimates will differ.

2. **P1 is a model-based resampling test, not an exact CRT as written.** Stage
   is observed, and its conditional distribution is estimated with a
   proportional-odds model. Exact CRT validity requires the true conditional
   law of stage given all variables on which allele counts can depend under
   the null. The proposed model uses cohort and `z`, while the likelihood also
   conditions on depth, chosen tag, error rate, ancestry-specific `q`, and
   het-call inclusion. If any of these predict stage after `z`, simulated
   stages can break the null exchangeability. Fix the conditioning set and
   define diagnostics/sensitivity for the stage model; call the resulting
   p-value approximate. Drawing stages independently also changes each
   cohort's stage counts. State whether counts are to be preserved and whether
   the stage model is refitted in each replicate. Specify how failed ML refits
   count toward 1,000 draws.

3. **P2 lacks a computable test definition.** Define both positive `λ` values,
   their units and uncertainty/source. State whether the denominator is fixed
   in the bootstrap. For each λ, equivalence requires both one-sided nulls
   (`κ/λ ≤ -0.125` and `κ/λ ≥ +0.125`) to be rejected, which is equivalent
   to a two-sided 90% interval wholly inside the margins when using matching
   one-sided 0.05 tests. Specify the interval method and p-value for Holm;
   requiring only a 95% bound is a different criterion. The τ condition is
   another component: define the one-sided upper confidence bound for τ/λ,
   the numerical `τ_max`, and how the compound P2 p-value is formed (maximum
   of component p-values for an intersection–union test). State whether λ is
   measured in logit units; otherwise τ/λ may be dimensionally incoherent.

4. **P3 simulation must preserve the selected observation process.** Drawing
   genotype and reads, then applying truncation to the *already selected*
   rows changes the row count and conditions on a different experiment. Either
   condition on each observed `(n, tag, inclusion)` and sample `a` directly
   from its truncated mixture, or simulate every eligible tag's depth and
   allele reads before choosing the max-depth tag and applying the call rule.
   The latter needs a depth model absent from v0. State which target null P3
   tests. For the conditional choice, `q`, `e`, `ω`, covariates and the observed
   tag/depth remain fixed, and the simulated `a` comes from
   `P(G) f_G(a|n)/Σ_G P(G) C_G(n)` over allowed counts. Define the bootstrap
   statistic and treatment of τ=0 boundary and failed refits.

5. **The genotype prior is not automatically `P(G|n, chosen tag, included)`.**
   The formula correctly updates Hardy–Weinberg priors for inclusion through
   `C_G`, *if* genotype is independent of depth/tag choice conditional on
   ancestry and covariates. Allele-dependent expression can affect depth, so
   max-depth tag choice may enrich genotypes or particular cis effects. State
   this assumption, evaluate it in genotype-known controls if available, and
   report sensitivity. Avoid interpreting the likelihood as genotype-verified
   ASE.

6. **Weak `α_g` makes proportional change unidentifiable.** When all included
   `α_g` are near zero, the data identify absolute stage slopes weakly but
   provide essentially no information on `κ`. Set a pre-outcome eligibility
   or information rule (for example based on external GTEx effect magnitude
   and expected depth), or explicitly allow and diagnose a flat/profiled
   likelihood. A bound on `κ` from an optimizer alone would be misleading.
   State parameter bounds and how genes with zero informative rows are handled.

## Arithmetic, boundary, and implementation clarifications

- The section 2 numerator and denominator are correct for a *single included
  row* conditional on `n`, provided all genotype classes use the same
  oriented `a` and symmetric allowed set. Give `K_c(n)` integer endpoints as
  `ceil(max(k_c, f_c n)) ≤ a ≤ n - ceil(max(k_c, f_c n))`; if this interval is
  empty, the row is ineligible. Define how floating-point ties at `f_c n` are
  treated.
- For homozygote components, the two orientation cases should be stated as
  oriented probabilities of `a`, including the `orientation = -1` swap for
  both R and A. This will prevent a REF/ALT reversal in one implementation.
- At `ρ=0` or `φ=0`, beta-binomial terms require the binomial limit. At
  `e_i=0`, allowed interior counts have zero homozygote probability, which
  needs an explicit limit or a predeclared positive floor. `q=0` or `q=1`
  eliminates the heterozygote component; specify whether such rows are
  excluded. Compute mixture and truncation with log-sum-exp/log-CDF or
  recurrence, not direct probability sums, for depths in the thousands.
- The term `ω_c·s_t` assumes reference bias is on the oriented allele logit
  scale and applies only to H. State the sign convention and whether its
  uncertainty is ignored. A fixed offset estimated from null-site hets can
  still induce a stage slope if bias changes with stage-dependent library
  properties.
- State whether `b̃` and `d̃` in the dispersion model are the same centered
  variables as the corresponding columns in `z`, and whether absent age/sex/PC
  values generate separate indicators by cohort. Specify the exact column
  order and scales for the fixture.
- `κ_A` says the `γ_k` terms enter the parentheses. This makes each lineage
  modifier proportional to `α_g`, unlike an additive gene-stage effect.
  Confirm that this is intended, and define whether `w̄_k` is a fixed
  gene-average or weighted by informative observations.
- The shared-parameter gradient threshold alone cannot establish convergence
  if gene-specific gradients remain large. Set tolerances for the full
  gradient and profile likelihood, plus rules for nonconverged fits in tests.
- For an analytic-gradient comparison, fix parameter transforms and the order
  of `α_g`, `r_g`, shared coefficients and τ. Compare gradients in the same
  unconstrained coordinates, including the truncation derivative. A relative
  tolerance is undefined at zero; use `atol + rtol × scale` there.

## Priority for v1

Resolve items 1–4 before freezing the fixture or hashing the prespecification.
Items 5–6 affect the biological interpretation and the Gate A information
calculation. The boundary clarifications are needed for reproducible REC/B-MODEL
agreement. The independent REC implementation will wait for v1 as assigned.
