# Model A-final likelihood specification, v0 (Claude-owned, pre-hash)

Status: v0, 2026-09-29. This is the implementation contract for two independent
implementations:
- **B-MODEL:** the Claude team.
- **REC:** Codex, task C2 in `2026-09-29_claude_to_codex_tasks.md`.

Neither implementation may read the other's code. The spec may change only through a
dated revision in this directory. B-SPEC copies the final version into
`scripts/analysis/in_biopsy_allelic/PRESPEC_A.md` before the hash. Margins, estimands and
tests are those of the approved plan and your decisions (Δ = 0.125 per stage step).
This file adds the arithmetic.

## 1. Observations

**Units.** Individuals i (one first-biopsy library each) in cohort c(i). Genes g (GTEx v8
liver eGenes with at least one tag).

For each individual × gene:
- **Chosen tag t(i,g):** the gene's tag with the largest total depth n in that library
  among all of the gene's tags. Ties go to the smallest genomic position. The choice
  does not look at allele counts.
- **n_ig:** reads carrying the tag REF or tag ALT base. Other bases are dropped.
- **a_ig:** reads carrying the tag allele that sits on the lead-ALT haplotype. That is
  the tag-ALT count when `orientation` = +1 and the tag-REF count when
  `orientation` = −1 (T1).
- **Inclusion (truncation).** The pair (i,g) enters only if cohort c's het rule
  R_c = (d_c, k_c, f_c) calls a het at the chosen tag. With m = min(a, n − a), the
  rule is n ≥ d_c and m ≥ k_c and m ≥ f_c·n. When the chosen tag is not called het,
  the pair is absent. Other tags are not tried.

**Fixed inputs**, computed before any allelic ratio is read:

| Input | Meaning |
|---|---|
| q_t | ALT allele frequency of tag t in the 1kGP superpopulation nearest the individual's ancestry projection |
| e_i | rate at which a homozygous site in library i shows reads of the specific other allele of a biallelic SNV (QC builder, identity sites) |
| φ | overdispersion of homozygote error reads (beta-binomial ICC), from the Geuvadis DNA-homozygous sites (T2) |
| ω_c | per-cohort reference-bias offset on the logit scale, from null-site hets (QC builder) |
| s_t | +1 if the lead-ALT allele at the tag is the tag REF (orientation −1), otherwise −1 |
| S_i | harmonized stage 0/1/2 |
| S̃_i | S_i centred within cohort (development individuals only) |
| z_i | nuisance covariates, centred within cohort: log library background b, duplication rate d, FFPE-likeness, sex, age, ancestry PC1–PC4 |

A missing covariate is replaced by the cohort mean, plus a missing-indicator column.

## 2. Genotype mixture with truncation

Genotype classes at tag t, under Hardy–Weinberg equilibrium:
- H (het), with P(H) = 2q(1 − q);
- R (hom REF), with P(R) = (1 − q)²;
- A (hom ALT), with P(A) = q².

Here q = q_t.

**Read model given class** (BB(n, μ, ρ) is beta-binomial with mean μ and
intra-class correlation ρ, i.e. shape parameters μ(1 − ρ)/ρ and (1 − μ)(1 − ρ)/ρ):
- **H:** a ~ BB(n, p_ig, ρ_ig), with logit p_ig = η_ig + ω_c·s_t.
- **R:** tag-ALT reads ~ BB(n, e_i, φ). In oriented terms:
  - when `orientation` = +1, a ~ BB(n, e_i, φ);
  - when `orientation` = −1, (n − a) ~ BB(n, e_i, φ).
- **A:** tag-REF reads ~ BB(n, e_i, φ). This is R with the roles of a and n − a
  swapped.

**Contribution of an included pair.** Depth n is conditioned on, not modelled:

  L_ig = Σ_G P(G)·f_G(a | n) / Σ_G P(G)·C_G(n),  with C_G(n) = Σ_{a' ∈ K_c(n)} f_G(a' | n),

where K_c(n) is the set of counts a' in 0..n that the rule R_c calls het at depth n.

## 3. Mean and overdispersion

**Mean (rank 1), on the natural-logit scale:**

  η_ig = α_g·(1 + κ·S̃_i + δᵀz_i)

- α_g: per-gene in-biopsy allelic effect.
- κ: shared proportional change per stage step. This is the primary estimand κ_U.
- δ: shared nuisance slopes.

**Overdispersion:**

  logit ρ_ig = r_g + ρ_S·S̃_i + ρ_b·b̃_i + ρ_d·d̃_i

**Parameter sets:**
- gene-specific: {α_g, r_g};
- shared: {κ, δ, ρ_S, ρ_b, ρ_d}.

**Estimation.** Joint maximum likelihood.
- Report the log-likelihood with the constant terms of the beta-binomial included (the
  full lbeta/lchoose form), so both implementations compare the same number.
- The SE of κ comes from the observed information after profiling the gene-specific
  parameters.
- Starting values:
  - α_g from the per-gene intercept-only fit;
  - κ = δ = 0;
  - r_g = logit(0.02).
- Convergence: gradient norm < 1e-6 on the shared parameters.

**Variants, same likelihood:**
- κ_hep: restrict the genes to hepatocyte-predominant ones (w_hep ≥ 0.9; list from
  B-COMP).
- κ_A: add shared terms γ_k·(w_gk − w̄_k)·S̃_i inside the parentheses for the
  hepatocyte, fibroblast and myeloid lineages. w_gk is the gene's lineage share at
  median development composition.

## 4. Heterogeneity τ

Per-gene stage slopes are allowed to deviate from κ·α_g:

  η_ig = α_g + (κ·α_g + u_g)·S̃_i + α_g·δᵀz_i,  u_g ~ N(0, τ²), τ ≥ 0

- The marginal likelihood integrates u_g with 20-node Gauss–Hermite quadrature per
  gene.
- τ is estimated jointly with the other parameters by maximizing that marginal
  likelihood.
- At τ = 0 this reduces to section 3 exactly, which is a required unit test.

## 5. Tests (development set)

**P1, κ_U = 0: conditional randomization test.**
- Fit a proportional-odds model of S given cohort indicators and z_i on development
  individuals.
- For each of 1,000 draws (seed 20260923 + draw index), sample S* for every individual
  from its fitted category probabilities, recentre within cohort, and refit the full
  model of section 3.
- p = (1 + #{|κ*| ≥ |κ̂|}) / 1001.

**P2, persistence** (intersection–union):
- A two one-sided test that κ_U/λ lies inside ±0.125, using a within-cohort
  individual-cluster bootstrap with 1,000 refits.
- The one-sided 95% bound of τ/λ must be ≤ τ_max.
- It must pass at both λ values (GSE193066 and PXD051911).
- τ_max is stated on the natural-logit scale and derived by B-SPEC.

**P3, τ = 0: parametric bootstrap.**
- Simulate 1,000 datasets from the fitted section-3 model: draw the genotype class,
  then the reads, apply the truncation and the depth-chosen tag, and refit section 4.
- p = (1 + #{τ* ≥ τ̂}) / 1001.

**Family.** {P1, P2, P3}, Holm at 0.05.

## 6. Agreement targets for REC and B-MODEL

On the shared synthetic fixture (written by B-FIX, independent of both
implementations):
- log-likelihood at the same fixed parameter vector within 1e-6;
- the gradient within 1e-6 relative;
- κ̂ within 1e-4;
- τ̂ within 1e-3;
- the P1 CRT p-value identical when both use the same S* draws, which B-FIX supplies
  as a file.

Please report any ambiguity in this spec rather than resolving it silently. Every
ambiguity found now is one the two implementations would otherwise resolve differently.
