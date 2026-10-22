# Model A-final likelihood specification, v1 (Claude-owned, pre-hash)

Status: v1, 2026-09-29. Written by the lead from five internal critics (R-STAT, R-GEN, R-NUM,
R-RED, R-SCOPE; results in workflow journal `wf_fdc0b317-25e`), Codex's critique
(`2026-09-29_codex_spec_critique.md`), Codex C4/C6/C8 results and the lead proposals L1–L18.
Items marked **[needs user]** keep the plan's version until the user decides (Appendix B).
Units are the natural logit throughout; a log2 display divides by ln 2.

## 0. Main changes from v0

1. The fit set is named: staged development unit individuals only (§1.1).
2. κ_U is the section-3 κ for the fixed eligible gene set; P1 is studentized with an
   individual-cluster sandwich SE (§5.1).
3. Section 4 uses adaptive Gauss–Hermite and the coordinate v = τ² ≥ 0 (§4, §7).
4. The truncation set uses one integer het predicate shared by the row builder, the normalizer
   and the simulators (§1.3).
5. All tag exclusions are applied before the depth-only tag choice; one hashed tag table (§1.2).
6. Genotype prior keeps third alleles and deletions; per-individual 1kGP frequencies (§2.1).
7. A shared stage slope of reference bias, ω_S, enters the het mean (§3.1).
8. λ is the covariate-conditional, information-weighted reliability; its uncertainty is
   propagated draw by draw (§5.2). The τ part of P2 uses a profile-likelihood signed root.
9. P3 simulates conditional on the selected rows; statistic is the LR (§5.3).
10. A sealed-analysis section, controls with written predictions, and a hash list (§6, §9, §10).
11. No depth cap. Seeds and RNG streams are fixed per test.

## 1. Observations

### 1.1 Units and fit sets

Source of truth: `seal-20260929T141204Z/frozen_crosswalk.tsv`, checked against
`frozen_crosswalk.sha256`. Units are never re-derived (t5a `choose_libraries` is not used).

| Set | Definition | Used for |
|---|---|---|
| **F (fit set)** | `unit_library` True, `run_role` development, `C` = 0, `S` observed, cohort in {GSE213621, GSE135251, GSE162694, GSE130970, GSE240729} (≥ 2 stage levels). 558 individuals: 208, 160, 83, 54, 53. | κ_U, τ, κ_hep, κ_A, P1–P3, λ target |
| **M (main-effect set)** | all development unit individuals (936), staged or not, controls included | the stage-free α_g Catalog fit only (§8.1) |
| **Sealed F** | as F with `run_role` sealed, cohorts GSE174478 plus the five above | §6 only |

GSE193066 enters only as a λ label-pair source. PRJNA512027's staged individuals are all S = 2
and GSE126848/GSE167523 have no stage, so these cohorts are in M only.

**Centring.** S̃_i = S_i − mean of S over all individuals of F in cohort c(i), whether or not they
have an included row (so centring never depends on allele counts). Every z column is centred
the same way. Both are recomputed in every CRT draw and bootstrap resample, using multiplicities.

### 1.2 Tag table and tag choice

B-SPEC writes one hashed tag table: `gene_id, tag chrom, pos, ref, alt, orientation, s_t,
exclusion flags`. Start: T1 v2 (`in-biopsy-allelic-tags-v2-20260929T134419Z/tags/chr*.tags.tsv`).
Remove, before any tag choice, in this order, recording each count:

1. non-autosomal (chrX/Y eGenes get Catalog reason "not tested: sex chromosome");
2. MHC chr6:28.5–33.5 Mb and GENCODE IG/TR spans (`excluded_spans_mhc_ig_tr.bed`; 56 SNVs);
3. alt_copy genes, pseudogenes, NUMTs, imprinted genes (the geneimprint list QC v2 uses, version
   named), multi-gene exons, REDIportal sites;
4. T3′ failures (T3′ on the v2 list, 706-contig STAR index);
5. tags without a row in B-Ancestry's `tag_2pq_by_superpop.tsv` and `tag_ld_by_superpop.tsv`
   (these tables must be extended to all v2 SNVs before the hash; 628 are missing now);
6. per cohort, in FFPE cohorts (GSE213621, GSE240729): transition tags (alleles {C,T} or {A,G}),
   unless that cohort's QC v2 null-site G1-FFPE verdict is pass (|median difference| ≤ 0.02).

**Chosen tag t(i,g):** among the gene's remaining tags, the one with the largest
n = REF reads + ALT reads in library i (T4 AD fields); ties go to the smallest position. The het
rule is then applied to that tag only; other tags are not tried. `t5a_allelic_table.py:294-329`
(G1-FFPE at tag hets, het-only tag choice, the detection guard) is superseded and must not be
reused. The PRESPEC_A_DRAFT §8 detection guard is removed because it reads allelic ratios.

n includes PCR duplicates (the BAMs carry no duplicate flags). No depth cap is applied: the
largest observed n is 6,821, below the T4 pileup limit of 10,000.

### 1.3 Het predicate and truncation set

Cohort rule R_c = (d_c, k_c, F_c), with F_c = 100·f_c an integer. For a library with REF count r
and ALT count x at the chosen tag, n = r + x and m = min(r, x):

  het(r, x) ⇔ n ≥ d_c and m ≥ k_c and 100·m ≥ F_c·n     (integer arithmetic)

K_c(n) = {a : lo ≤ a ≤ n − lo}, lo = max(k_c, ⌈F_c·n/100⌉) computed with integer ceil division;
K_c(n) is empty when n < d_c or lo > n − lo. One function implements this and is imported by the
row builder, the normalizer and every simulator. A unit test checks the set against a direct
enumeration of het() for n = 0..2,000.

### 1.4 Oriented count

a_ig = tag-ALT reads when orientation = +1, tag-REF reads when orientation = −1. So a counts the
tag allele on the lead-ALT haplotype. s_t = +1 when the lead-ALT allele is the tag REF
(orientation −1), else −1.

## 2. Genotype mixture with truncation

### 2.1 Prior

For individual i and tag t, from public 1kGP only (never from MASLD data):
- called superpopulation EUR, AFR, AMR, SAS: that superpopulation's p_ref, p_alt, and
  p_x = other_alt_freq + deleted_freq + conflict_freq;
- EAS-called individuals of the Japanese cohorts (GSE167523, GSE174478, GSE193066): JPT values;
  other EAS-called individuals: EAS;
- individuals flagged intermediate by B-Ancestry: each of p_ref, p_alt, p_x is interpolated,
  p = p_EUR + a1·(p_AFR − p_EUR) + a2·(p_EAS − p_EUR), with a1, a2 the two axis coordinates scaled
  to 0 at the EUR centroid and 1 at the AFR or EAS centroid, clipped to [0, 1].

Clamp p_ref and p_alt to [1/(2N), 1 − 1/(2N)], N = unrelated founders of the reference used
(AFR 677, AMR 351, EAS 505, EUR 524, SAS 509, JPT 104; `tag_ld_summary.json`).

Unnormalized class weights (REF/X and ALT/X carriers show one counted allele):
- P(H) = 2·p_ref·p_alt
- P(R) = p_ref² + 2·p_ref·p_x
- P(A) = p_alt² + 2·p_alt·p_x

The row likelihood is a ratio, so a common scale cancels. The 35 1kGP-multiallelic tags stay;
dropping them is a sensitivity.

### 2.2 Oriented read model

With q_o the frequency of the lead-ALT-side tag allele, the three classes in oriented terms are:

| Class | Prior | Distribution of a |
|---|---|---|
| H | P(H) | BB(n, p_ig, ρ_ig), logit p_ig = η_ig + s_t·(ω_c + ω_S·S̃_i) |
| O0 (homozygous for the allele on the lead-REF haplotype) | P(R) if orientation +1, P(A) if −1 | BB(n, e_i, φ) |
| O1 (homozygous for the lead-ALT-side allele) | P(A) if orientation +1, P(R) if −1 | BB(n, 1 − e_i, φ) |

BB(n, μ, ρ): beta-binomial with shapes α = μ(1−ρ)/ρ, β = (1−μ)(1−ρ)/ρ, log pmf =
lchoose(n, a) + lbeta(a + α, n − a + β) − lbeta(α, β). The binomial is used only when the fixed
input φ is exactly 0.

Required unit tests: the hand example below to 1e-4, and orientation invariance (flip
orientation, which flips s_t and the O0/O1 mapping of P(R) and P(A); replace a by n − a; negate
η: every row likelihood unchanged).
Hand example (binomial limits): n = 10, rule (10, 2, 10) so K = {2..8}, orientation −1, q_t = 0.3,
ω = 0.1, η = 0.5, e = 0.01, tag-REF 7, tag-ALT 3 so a = 7; f_H = C(10,7)·p⁷(1−p)³ with
p = logistic(0.6) = 0.64566 gives 0.24973.

### 2.3 Row contribution

  L_ig = Σ_G P(G)·f_G(a | n) / Σ_G P(G)·C_G(n),   C_G(n) = Σ_{a' ∈ K_c(n)} f_G(a' | n)

Depth n is conditioned on. Homozygote error e_i: the QC v2 per-library e_i where the cohort's
split-half Spearman ≥ 0.5, otherwise the cohort median of unit_all e_i; NA → cohort median;
floored at 1e-4. φ: beta-binomial ML of other-allele reads at Geuvadis DNA-homozygous sites over
the calibration depth range, with each individual's pooled e fixed (B-T2 produces it).
Sensitivity: tag-specific e = eps_i + c_i·f_other(t), floored.

**Stated assumption.** P(G | n, chosen tag, covariates) = the §2.1 prior up to the inclusion
reweighting. Cis effects can couple genotype and depth; the effect enters only through the
false-het share and is second order at FDP ≤ 0.02. Checks: §9 (het-only sensitivity with the
truncation kept; posterior calibration at Geuvadis DNA-truth sites).

## 3. Mean and overdispersion

### 3.1 Mean (rank 1)

  η_ig = α_g·(1 + κ·S̃_i + δᵀz_i)

The H mean adds the additive technical term s_t·(ω_c + ω_S·S̃_i) (§2.2). ω_S is a shared stage
slope of reference bias; it is nearly orthogonal to κ because s_t is ±1 across genes
independently of sign(α_g).

κ_U is the proportional change per stage step, relative to the effect at the cohort-mean stage,
for the fixed set of eligible genes (§3.4). It is estimated by the section-3 fit.

### 3.2 Overdispersion

  ρ_ig = ε + (1 − 2ε)·expit(r_g + ρ_S·S̃_i + ρ_b·b̃_i + ρ_d·d̃_i),  ε = 1e-6

b̃ and d̃ are the same centred columns as in z. No switch to the binomial.

### 3.3 z columns (fixed order)

1. b: log library background (QC v2 column, name bound in the parameter file)
2. d: duplication rate at tag regions (B-Dup)
3. FFPE-likeness (QC v2 column, name bound in the parameter file)
4. sex (1 = female), from metadata or expression sex
5. age / 10
6. axis_EUR_AFR, 7. axis_EUR_EAS (B-Ancestry `library_ancestry_pcs.tsv`)
8. c_i, 9. eps_i (QC v2 contamination and error split)

A column constant within a cohort is 0 there after centring; a column 0 in every cohort is
dropped. A missing-indicator column is added only where a covariate is partly missing within a
cohort, with the missing values set to the cohort mean before centring.

### 3.4 Gene eligibility and bounds

- A gene is κ-eligible if GTEx |ln aFC| ≥ θ and its expected number of included rows in F is
  ≥ 10. The expectation uses only depths, the §2.1 prior, R_c and α_g^GTEx; the list is frozen
  and hashed at Gate A. θ is chosen at Gate A from {0, 0.1, 0.2, 0.3} to minimize the design MDE
  (ties to the smaller θ). No minimum based on called-het counts is used.
- |log2_aFC| at the GTEx cap log2(100) is treated as censored; Gate A plants cap |α| at 3.
- Bounds: α_g ∈ [−6, 6]; r_g ∈ [−15, 5]; κ ∈ [−0.9, 0.9]; ρ_S, ρ_b, ρ_d, ω_S ∈ [−3, 3]. Any
  estimate on a bound is reported with the gene named. A gene with no included row gets no
  parameters.

### 3.5 Variants on the same likelihood

- **κ_hep:** genes with w_hep ≥ 0.9 (B-COMP), same fit.
- **κ_A:** add γ_fib·w_g,fib·S̃_i + γ_mye·w_g,mye·S̃_i inside the parentheses, with uncentred
  non-hepatocyte shares (w at median development composition). κ_A is then the change for a gene
  whose mRNA is entirely hepatocyte. The lineage set is fixed (hepatocyte, fibroblast, myeloid);
  PRESPEC_A_DRAFT §6's data-driven lineage choice is superseded. Report the condition number of
  the γ block.

## 4. Heterogeneity τ

  η_ig = α_g + (κ·α_g + u_g)·S̃_i + α_g·δᵀz_i,   u_g ~ N(0, v),  v = τ² ≥ 0

**Gene contribution** (one integral per gene, normalizers inside):

  ℓ_g(v) = log ∫ N(u; 0, v) ∏_{i∈g} L_ig(u) du

**Adaptive Gauss–Hermite, 20 nodes** (numpy `hermgauss(20)`, physicists' nodes x_k, weights
w_k). h(u) = Σ_i log L_ig(u) − u²/(2v). û_g = mode of h by Newton to |Δu| < 1e-10;
σ̂_g = (−h''(û_g))^(−1/2).

  ℓ_g = log[ √2·σ̂_g · Σ_k w_k·exp(x_k²)·exp(h(û_g + √2·σ̂_g·x_k)) ] − ½·log(2π·v)

computed by log-sum-exp. Gradients hold (û_g, σ̂_g) fixed at their values for the current
parameters. At v = 0, ℓ_g = Σ_i log L_ig(0) exactly (separate branch).

**Score at the boundary:** ∂ℓ/∂v at v = 0 is ½·Σ_g[(Σ_i ∂_u log L_ig)² + Σ_i ∂²_u log L_ig],
evaluated at u = 0. v̂ = 0 when this score ≤ 1e-6 and the fit started at v = τ_max² returns to 0.

**Accuracy test:** for five fixture genes spanning 5–500 rows, ℓ_g matches `scipy.integrate.quad`
(epsabs 1e-12) within 1e-8 relative at τ ∈ {0.25, 1, 4}·τ_max, and adaptive 20 matches adaptive
40 nodes within 1e-8 on the fixture.

**Stated approximation.** With u_g random, the inclusion pattern carries information about u_g;
keeping the N(0, v) prior makes section 4 a pseudo-likelihood. Gate B reports τ̂ bias and τ-bound
coverage from data generated with u_g drawn before reads and het calls. If |bias| > 0.1·τ_max,
the section-4 fit is restricted to rows with n ≥ n*_c, the smallest n at which R_c calls a het of
allelic fraction 0.2 with probability ≥ 0.99 under BB with ρ = 0.02.

## 5. Tests (development fit set F)

Random streams: numpy `default_rng(SeedSequence([20260923, test_id, b]))`, test_id P1 = 1,
P2 = 2, P3 = 3, draw b = 0, 1, …. Rows are ordered by (individual_id, gene_id). B-FIX writes the
S* draws, the P2 resample indices and the P3 uniforms as files for the agreement check.
B is fixed in advance and is not extended after a p-value is seen. Every p-value is reported with
its Monte Carlo SE.

### 5.1 P1: κ_U = 0, conditional randomization test with an estimated stage law (approximate)

- **Statistic:** z = κ̂/SE, section-3 fit. SE² = Σ_c (e_κᵀ H⁻¹ s_c)², H the analytic observed
  information of all free section-3 parameters at the MLE (bound coordinates removed), s_c the
  analytic score summed over all rows of cluster c. Clusters are individuals, merged into
  first-degree relative clusters from `kinship_related_pairs.tsv`.
- **Why studentized:** the null κ = 0 allows ρ_S ≠ 0, ω_S ≠ 0, and per-gene departures u_g, under
  which counts still depend on S. For the fixed-gene estimand the sandwich-studentized z is
  pivotal under both S and S*; |κ̂| alone is not.
- **Stage model:** ordinal logit of S with cohort-specific thresholds and common slopes on, in
  order: the §3.3 z columns; log e_i; log median n over the individual's candidate rows (chosen
  tag with n ≥ d_c, het or not); log number of candidate rows. All are depth-only or allele-free.
  Composition is left out because it mediates κ. Fitted once on F.
- **Draws:** B = 1,000; S*_i drawn independently from the fitted category probabilities
  (cohort counts not preserved); recentre; refit section 3; recompute SE.
- **p** = (1 + #{|z*| ≥ |z| − 1e-6}) / (B + 1), two-sided. A refit that fails after 3 restarts
  counts as an exceedance and is listed. If the observed fit fails, P1 is "not computed".
- **Sensitivity:** stratified permutation of S within cohort, valid only if S is independent of
  z within cohort.

### 5.2 P2: persistence (intersection–union)

**λ.** For source l ∈ {GSE193066, PXD051911 arm A1}:
- λ_l,c = Var(E[S|T]) / Var(S) at cohort c's latent stage distribution, fitted from the labels of
  F's individuals in c (B-LABEL / Codex C8 definition; cutpoints 0.5…3.5, S cutpoints 1.5, 2.5);
  GSE213621's S comes from its condition label as in the seal.
- Carried label spec per source: the smaller of Y5 and S3 at the point estimate, fixed across
  draws.
- Covariate-conditional, information-weighted pooling. The pooled κ̂ is approximately
  Σ_c I_c κ̂_c / Σ_c I_c, and each κ̂_c is attenuated by λ_cond,l,c = (λ_l,c − R²_c)/(1 − R²_c)
  (label error assumed independent of z), so
  λ_l = Σ_c I_c·λ_cond,l,c / Σ_c I_c,
  with R²_c the within-cohort R² of S on z in F, and I_c the Gate A expected Fisher information
  for κ from cohort c, computed with the z-adjusted stage variance (outcome-free). I_c and R²_c
  are hashed.
- Draw b: var_e from B-LABEL bootstrap draw b (pair tables), target fits and R² recomputed on κ's
  resample b. λ*_b ≤ 0.05 counts as an exceedance in every component.

**Resampling.** B = 2,000. Within cohort, relative clusters are resampled with replacement; the
gene set is fixed; S̃ and z are recentred; section 3 is refit.

**κ component.** p_up = (1 + #{κ*_b/λ*_b ≥ +Δ})/(B + 1); p_lo = (1 + #{κ*_b/λ*_b ≤ −Δ})/(B + 1),
Δ = 0.125. The primary reported interval for κ/λ is the bootstrap 90% interval; the profile
interval of κ is secondary.

**τ component.** τ0_l = τ_max·λ_l^(0.05), with λ_l^(0.05) the 5th percentile of the λ*_b draws.
r = sign(τ0 − τ̂)·√(2[ℓ_p(τ̂) − ℓ_p(τ0)]) from the section-4 profile log-likelihood;
p_τ = 1 − Φ(r). If Gate B shows a rejection rate above 0.06 at τ = λ·τ_max, the fallback (fixed
now) is a conditional parametric bootstrap of r at τ = τ0 with the §5.3 machinery, 1,000 draws.

**P2 p-value** = max over {p_up, p_lo, p_τ} and over both sources. GSE193066 currently gives the
lower λ (Codex C8: 0.719 Y5 at 350 Kleiner-cohort people; staged6 value to be recomputed).

**τ_max** [needs user, B1]. v1 implements the user's literal decision: τ_max = Δ × median over
κ-eligible genes of |α̂_g| from the development section-3 fit, computed by hashed code in the
development run. θ is fixed before τ_max. The GTEx value (Codex C5: 0.0498 at θ = 0 to 0.0765 at
θ = 0.3) is used only as the planning value for Gate A and Gate B.

### 5.3 P3: τ = 0, conditional parametric bootstrap

- Fit section 3 (the null). For b = 1..1,000, hold every included row (i, g, chosen tag, n)
  fixed; draw G with probabilities P(G)C_G(n)/Σ P(G)C_G(n), then a by inverse CDF over K_c(n)
  under f_G, from B-FIX uniforms.
- Statistic LR = 2·(max(ℓ4(v̂), ℓ3) − ℓ3), section-4 fits started at v = 0 and v = τ_max² (the
  better one kept).
- p = (1 + #{LR* ≥ LR − 1e-6})/(B + 1). A failed refit counts as an exceedance.
- This tests τ = 0 conditional on the selected design; the tag choice is a function of depth,
  which is held fixed. The 50:50 χ̄² is a fixture diagnostic only (genes share individuals).

### 5.4 Family

{P1, P2, P3}, Holm at 0.05. If Gate A declares P2 untestable, the family keeps m = 3 with
p₂ = 1 and P2 reported as untestable [needs user, B5].

## 6. Sealed analysis (opened once)

- Cohort inputs (R_c, e_i, ω_c, φ, background, duplication, chrX audit) are computed on
  unit_sealed tables and hashed before opening.
- **Sealed fit:** section 3 refit on sealed F (GSE174478 whole plus the 20% hash parts), all
  parameters free; α_g is not transported.
- **Sealed P1:** the §5.1 CRT with the stage model refit on sealed F; replication = same sign as
  development κ̂ and p < 0.05.
- **Sealed λ:** the §5.2 construction on sealed F.
- **P3 projection statistic:** Spearman correlation across genes between the sealed per-gene
  slope departures (posterior means of u_g at the development v̂) and the development posterior
  means; null from the §5.3 conditional bootstrap under the sealed section-3 fit, 1,000 draws;
  one-sided.
- **Sealed MDE:** from the Gate A simulation on sealed depths and stages, written before opening.
- **Deviation sensitivity:** drop the 3 FAM86JP tag positions for the 6 sealed GSE130970 unit
  individuals whose depths were displayed on 2026-09-25 (only if FAM86JP survives §1.2).

## 7. Computation

- **Coordinates and order:** θ = [α_1..α_G, r_1..r_G (genes by gene_id ascending), κ, δ (§3.3
  order), ω_S, γ_fib, γ_mye (κ_A only), ρ_S, ρ_b, ρ_d, v]. Gradients are compared in these
  coordinates with atol 1e-8 + rtol 1e-6·|value|, at v > 0 only.
- **Normalizers:** C_R and C_A depend only on fixed inputs and are precomputed per row by direct
  log-sum-exp over K_c(n) (unit test: C_R = C_A for symmetric K). For H: T = mass of the two
  excluded tails (a < lo, a > n − lo) by log-sum-exp; log C_H = log1p(−T) when T < 0.5, else a
  direct sum over K. Terms use the recurrence
  log f(a+1) − log f(a) = log(n−a) + log(a+α) − log(a+1) − log(n−a−1+β).
- **Starting values:** per gene, the truncated mixture with ω, κ = δ = ω_S = 0 and dispersion
  slopes 0, best of α ∈ {−1, 0, 1}, r = logit(0.02). Refits warm-start from the observed MLE.
- **Optimizer and convergence:** L-BFGS-B with the §3.4 bounds, then a Newton polish. Converged
  when the projected max |gradient| < 1e-6 over all free coordinates; at v = 0 the KKT condition
  ∂ℓ/∂v ≤ 1e-6. A fit not meeting this after 3 restarts is a failed fit.
- **Cost** [Estimate, not run]: about 5×10⁴ included rows; section-3 fit 5–60 s; section-4 fit
  2–25 min. P1 2–17 CPU-h; P2 κ part 4–34 CPU-h; P3 35–400 CPU-h. B-MODEL times one fit of each
  kind at real scale before Gate B.

## 8. Stage-free outputs

### 8.1 α_g Catalog fit

Section-3 model with κ = 0 and no stage terms on set M, per gene; it fills
`cis_allelic_effect_in_MASLD_liver`, which never maps to tested_negative.

### 8.2 Errors-in-variables β

The slope of α̂_g on GTEx ln aFC is printed beside the attenuation expected from 1kGP tag–lead
LD, which uses no outcome. β < 1 is not called disease-related unless it falls below that value.

## 9. Controls and diagnostics, predictions written now

| Check | Data | Prediction |
|---|---|---|
| Negative control 1 | per-library null-site REF logit on S̃ + z within cohort, development | \|β̂·M\| < MDE/4, with M = Σ W_g α_g^GTEx s̄_g / Σ W_g (α_g^GTEx)² from the Gate A design |
| Negative control 2 | male non-PAR chrX false-het rate on S̃ within cohort | 95% CI covers 0 |
| Hom-error calibration | male chrX observed/predicted false-het rate under R_c | in [0.5, 2]; if > 2, e_i for that cohort is multiplied by the ratio (rule fixed now) |
| Depth–stage coupling | per gene and cohort, Spearman of chosen-tag n with S (depth and stage only) | median \|ρ\| < 0.1 |
| Tag switching | chosen tag ≠ fixed tag, regressed on S̃ within cohort | slope within ±0.02 per step |
| Positive control | κ in genes with w_fib + w_mye ≥ 0.5 | sign of the composition-predicted change (from w_gk and composition slopes, outcome-free) |
| G2 | genes with \|α̂_g\|/SE > 2 | sign agreement with GTEx ≥ 80% |
| A9 | PNPLA3 148M | as in `a9_coding/` |

**Sensitivities (prespecified, reported):** transversion-only (all transition tags dropped in
every cohort, sections 3 and 4 refit); one fixed tag per gene and cohort (largest summed depth
over development libraries); ancestry-matched tag validity (r² ≥ 0.8 with the same sign in the
individual's superpopulation) [needs user whether primary, B7]; L16 low-r² rows dropped (θ and
τ_max unchanged); the 35 multiallelic tags dropped; tags with 1kGP het excess
(obs/expected > 1.2, n_hap ≥ 400) dropped; het-only likelihood with truncation; tag-specific
e = eps_i + c_i·f; 2·e_i; ω_c = 0 and at its CI endpoints; libraries above the cohort 90th
percentile of c_i dropped; libraries with d_i above the cohort median dropped; A/G–T/C tags versus
other tags (prediction: CI of the difference covers 0); detection guard n ≥ n_min(g,c)
(depth-only, from α^GTEx); δ_n·(log n centred within gene × cohort) inside the parentheses.

**Rule fixed now:** if the Gate B LD-mixture plant (tag hets drawn from the coupling/lead-hom/
repulsion mixture implied by r² and frequencies, depth falling 20% per step) gives
|κ bias| > Δ/10, the δ_n term moves into the primary model before the hash.

## 10. Hash list

Parameter file entries, each with sha256: R_c for every cohort (development and sealed parts);
the post-exclusion tag table and the exclusion counts; `tag_2pq_by_superpop.tsv`,
`tag_ld_by_superpop.tsv` and the individual-to-reference map; e_i (floor, fill), φ, ω_c with CI;
z values, column names and order, centring set; `frozen_crosswalk.tsv`; θ with its Gate A table;
the κ-eligible gene list; W_c, R², the B-LABEL draws file and the carried spec; w_hep list and
w_gk; the stage-model formula; the RNG scheme; the B-FIX fixture and its expected values;
`seal_order_of_events_lead_addendum.txt`; this spec and the code.

## 11. Agreement targets for REC and B-MODEL (Gate B)

On the B-FIX fixture: log-likelihood at fixed parameter vectors within 1e-6 absolute (vectors use
ρ ≥ 1e-5); gradients within atol 1e-8 + rtol 1e-6; κ̂ within 1e-4; τ̂ within 1e-3 where LR > 1,
otherwise LR within 1e-6; P1 exceedance indicators identical for every draw with
||z*| − |z|| > 1e-3 (draws inside the band listed); P2 and P3 p-values identical up to the same
tie band, given the shared draw files. The fixture has zero failed refits.

## Appendix A. Disposition

| Item | Decision | Where |
|---|---|---|
| Codex C1 (integral per gene) | accepted; adaptive GH added | §4 |
| Codex C2 (P1 approximate, conditioning set) | accepted; studentized sandwich statistic | §5.1 |
| Codex C3 (P2 definition) | accepted; conditional λ, signed-root τ part | §5.2 |
| Codex C4 (P3 target) | accepted, conditional simulation | §5.3 |
| Codex C5 (genotype–depth) | accepted as stated assumption plus checks | §2.3, §9 |
| Codex C6 (weak α) | accepted; θ plus expected-row minimum, bounds | §3.4 |
| Codex boundary bullets | accepted: integer K, oriented table, smooth ρ floor, e floor, ω sign, b̃/d̃, κ_A, convergence, coordinates | §1.3, §2.2, §3, §7 |
| L1 | modified: adaptive GH | §4 |
| L2 | modified: sandwich z, depth-only stage-model terms | §5.1 |
| L3 | modified: conditional information-weighted λ, per-draw var_e, B = 2,000 | §5.2 |
| L4 | replaced by the user's literal rule, GTEx for planning only | §5.2 |
| L5 | accepted with LR and two starts | §5.3 |
| L6 | accepted, checks added | §2.3, §9 |
| L7 | modified: expected-row minimum, θ before τ_max | §3.4 |
| L8 | modified: integer predicate, smooth ρ floor, e floor 1e-4 | §1.3, §3.2 |
| L9 | modified: unit_all, pooled value for non-estimable cohorts, ω_S added | §2.2, §3.1 |
| L10 | modified: v = τ², KKT at 0 | §4, §7 |
| L11 | modified: no depth cap; normalizer by mass | §1.2, §7 |
| L12 | modified: hepatocyte-only reference, two γ | §3.5 |
| L13 | accepted | §3.3 |
| L14 | accepted; axes also in the stage model | §3.3, §5.1 |
| L15 | modified: keep multiallelic tags with merged prior; JPT for Japanese cohorts; interpolation for intermediates | §2.1 |
| L16, L17 | accepted | §1.2, §9 |
| L18 | accepted; counts corrected in the addendum; FAM86JP sealed sensitivity | §6 |
| R-STAT-1, R-NUM-3, R-SCOPE-12 | accepted (adaptive GH) | §4 |
| R-STAT-2 | accepted as stated approximation with a fallback rule | §4 |
| R-STAT-3, R-NUM-7 | accepted | §5.1 |
| R-SCOPE-2, R-RED-1 (section-4 z for P1) | rebutted: κ_U is defined for the fixed gene set (R-STAT-8), for which the sandwich-studentized section-3 z is pivotal; section-4 refits in 1,000 CRT draws would cost 35–400 CPU-h more. Section-4 κ̂ is reported as secondary. | §3.1, §5.1 |
| R-STAT-4, R-RED-12, R-RED-13 | accepted (stage-model terms, c_i/eps_i in z) | §3.3, §5.1 |
| R-STAT-5, R-SCOPE-4, R-RED-8 | accepted | §5.2 |
| R-STAT-6, R-STAT-14, R-NUM-8, R-SCOPE-22 | accepted | §5 |
| R-STAT-7, R-RED-9 | accepted (signed root with λ 5th percentile; R-RED-9's form) | §5.2 |
| R-STAT-8 | accepted | §3.1 |
| R-STAT-9, R-NUM-14, R-RED-11 | Gate B procedure changes → user | B3 |
| R-STAT-10 | accepted for Gate B (label layer simulated) | B3 |
| R-STAT-11, R-GEN-9, R-NUM-1, R-RED-10, R-SCOPE-1 | accepted | §1.1 |
| R-STAT-12, R-SCOPE-10(1) | accepted | §6 |
| R-STAT-13, R-NUM-5, R-SCOPE-11 | accepted with R-SCOPE-11's outcome-free minimum | §3.4 |
| R-STAT-15, R-RED-2, R-GEN-7(1) | accepted; ω_S in the primary model | §3.1, §9 |
| R-STAT-16, R-NUM-2, R-RED-18, R-SCOPE-13 | accepted | §4, §7 |
| R-STAT-17 (depth cap) | rebutted: four critics reject; max n 6,821 is below the T4 limit, thinning would change inclusion; cost is feasible without it | §1.2 |
| R-STAT-18, R-GEN-10, R-GEN-11, R-RED-6, R-RED-15, R-SCOPE-14, R-SCOPE-24 | accepted; multiallelic tags kept (R-GEN-11's exclusion is the sensitivity) | §2 |
| R-STAT-19, R-STAT-20, R-SCOPE-19 | accepted | §1.3, §3.2 |
| R-GEN-1 | accepted | §2.2 |
| R-GEN-2, R-RED-4 | accepted as plant plus rule and sensitivities | §9 |
| R-GEN-3 | sensitivity; primary → user | B7 |
| R-GEN-4, R-GEN-5, R-NUM-6, R-SCOPE-8, R-RED-17, R-SCOPE-7 | accepted | §1.2 |
| R-GEN-6, R-RED-3, R-SCOPE-9 | accepted; "disagrees" threshold → user | §1.2, §9, B4 |
| R-GEN-8, R-SCOPE-10(2,3) (multi-tag variant, σ_tech) | proposed dropped → user | B6 |
| R-GEN-12, R-RED-5 | accepted as sensitivity and diagnostic | §9 |
| R-GEN-13 | accepted | §8.2 |
| R-NUM-4 | accepted | §7 |
| R-NUM-9 | accepted (cost, warm starts) | §7 |
| R-NUM-10 | accepted (ω logit, e fill and floor, φ producer) | §2.3 |
| R-NUM-11 | accepted | §7 |
| R-NUM-12 | accepted | §3.5 |
| R-RED-7 | accepted (κ_A reference); τ_A narrowing → user | §3.5, B4 |
| R-RED-14, R-RED-19, R-RED-20 | accepted | §9, §3.4 |
| R-RED-16 | numeric narrowing rules → user | B4 |
| R-SCOPE-5, R-SCOPE-6, R-SCOPE-18 | accepted; PRESPEC_A_DRAFT §1 Q4–5, §3 (GSE193066), §5, §6, §7 G3, §8 stage tests, §9 marked superseded in the B-SPEC amendment | §1.1, §3.5 |
| R-SCOPE-15 | accepted | §10 |
| R-SCOPE-20 | accepted for PRESPEC wording (Codex boundary paragraph; "removes factors that act equally on both alleles"; Currin 2025; EFHD1 nominal; "source describes 368 patients; 25 genotype-identical library pairs") | PRESPEC |
| R-SCOPE-21 | accepted; addendum counts corrected | §6 |
| R-SCOPE-3, R-NUM-13, R-GEN-14(1), R-RED-21(1), R-STAT-21(5) | τ_max → user; v1 uses the literal decision | §5.2, B1 |
| R-SCOPE-16, R-STAT-21(3) | → user, conservative default | §5.4, B5 |
| R-SCOPE-17, R-STAT-21(4) | → user, conservative default | B2 |
| R-SCOPE-23 | → user (GSE240729 consent) | B8 |

## Appendix B. Decisions for the user

**User answers 2026-09-29:** B1 = the literal rule (development median |α̂|, hashed code; GTEx for planning only); B3 = revised Gate B procedure, array jobs approved for Gate B; B6 = drop the multi-tag variant and σ_tech; B8 = keep GSE240729 in aggregate allelic statistics, genotypes restricted; B9 = keep the seal. B2, B4, B5 and B7 stay at the defaults below unless the user changes them.

- **B1. τ_max source.** v1: your literal rule, median |α̂| from the development stage-free fit,
  computed by hashed code after the hash. Alternative: GTEx median |ln aFC| (0.050–0.077), which
  is known before the hash but sits on a larger scale than biopsy α, so the τ margin is looser.
- **B2. λ source withdrawn.** If the GSE193066 arm is not carriable (identity rerun), default: P2
  cannot pass. Alternative: P2 untestable.
- **B3. Gate B procedure.** Size and coverage bands judged on warp-speed Monte Carlo at ≥ 4,000
  replicates; the 200 full-CRT replicates judged one-sided (fail if the binomial lower 95% bound
  of size exceeds 0.05); add P2 boundary plants κ = ±Δλ and τ = λ·τ_max, the technical plants
  (reference-bias slope, doubled false-het rate at S = 2, non-BB tails with depth–stage coupling,
  tag switching, LD mixture, omitted depth-dependent bias) and a simulated label layer; approve
  array jobs for Gate B.
- **B4. Numeric narrowing rules.** Proposed: transversion-only "disagrees" = paired-bootstrap 95% CI
  of κ_TV − κ_U excludes 0; leave-one-cohort-out "flip" = a LOCO CI excluding 0 on the opposite
  side, applied only when P1 rejects; "compatible with composition" = κ_A's CI covers 0 while κ_U's
  excludes it; drop "κ_hep identical to κ_U"; add "P3 rejects but τ_A does not: heterogeneity
  compatible with composition".
- **B5. Holm when P2 is untestable.** Default m = 3 with p₂ = 1. Alternative: Holm over {P1, P3}.
- **B6. Multi-tag variant and σ_tech.** Proposed: drop both (the joint LD/phase model is large and
  the fixed-tag sensitivity covers tag switching); replace "τ below σ_tech" with the technical
  plants and negative controls.
- **B7. Ancestry-matched tag validity.** Default: sensitivity. Alternative: primary.
- **B8. GSE240729.** Its 67 retrospective donors did not each give written consent (coded-anonymous
  residual tissue; Codex C4). It contributes 53 staged development and 13 sealed unit individuals.
  Default: keep in aggregate allelic statistics under decision 1, genotypes restricted.
- **B9. Seal.** Keep the seal despite the 2026-09-25 format-check reads (lead decision, flagged).

## Revision v1.1 (2026-09-29, lead): inputs bound, no change to any formula

- **Tag table (§1.2) built:**
  `in-biopsy-allelic-t2b-20260929T220019Z/tag_table/tag_table.tsv`, sha256 `b061a7c2…`.
  - 2,343 kept SNVs in 830 genes. The MHC/IG/TR spans hold 88 v2 SNVs, not 56.
  - T3′ uses the union rule across cohorts. It reproduces T4 exactly, with BAQ on.
  - Transition tags are dropped in GSE213621 and GSE240729.
  - The L16 low-r² flag now covers 5 rows (3 kept), because with the corrected founder set
    the earlier 106 rows mostly disappear.
- **Het rules R_c (§1.3):** `cohort_rules.tsv` in the same directory.
  - Grid: d ∈ {10, 20, 30}, k ∈ {1..5}, F ∈ {1, 2, 5, 10, 15, 20}. It was extended from k ≤ 3 and
    F ≤ 10 before any allelic data were read, because two high-background cohorts had no rule
    at the old ceiling.
  - Choice: among rules with FDP ≤ 0.02 at the cohort's e_p90 and sensitivity ≥ 0.85, the one
    with the largest expected number of true-het calls; ties go to the larger d, then the larger k.
- **φ = 0.0439**, the Geuvadis beta-binomial ML at DNA-homozygous sites.
- **ω_c:** GSE135251 (in F) is not estimable and takes the pooled development value,
  logit(0.5238).
- **BAQ contingency:** if Gate A is borderline, the T4 pileup is rerun without BAQ, and T3′ is
  rerun to match (2,927 tags, 986 genes). The decision rests only on the Gate A design MDE.

## Revision v1.2 (2026-09-30, lead)

The nine answers in `2026-09-30_claude_to_codex_design_fixture_answers.md` §3 are part of the spec.
They cover:
1. invalid λ draws in the τ component;
2. the n* fallback, now relative to the infinite-depth ceiling;
3. τ_max from the section-3 fit on F;
4. log1p stage-law predictors;
5. missing-indicator order;
6. adjusted R² for λ_cond;
7. the z binding, with eps_i dropped;
8. clusters that cross cohorts;
9. cohort-mixture priors in the Gate A export only.

The B-FIX fixture is at `model-a-bfix-20260930T115521Z` (job 1209).

## Revision v1.3 (2026-09-30, lead)

- **§5.2 R²_c is fixed.** R²_c is held fixed across P2 draws at the observed adjusted value,
  `model-a-gateA-r2-20260930T120017Z/r2_moments_observed.json`. Resample R² is biased upward by
  the duplicated individuals. Per-draw R² is a sensitivity only.
- **§5.2 resamples are fixed.** The real-data P2 resamples are the file
  `p2_resample_development.tsv.gz` in the same directory: B = 2,000, draw rule stated in that
  directory.

## Revision v1.4 (2026-09-30, lead)

The binding text is §7 of `2026-09-30_claude_to_codex_design_fixture_answers.md`:
- z = b (log e_model), d, FFPE (per library only where reliable), sex, age/10, and the two
  axes; c_i is dropped;
- sealed technical inputs use unit_all;
- GSE213621 stays FFPE-treated on its null-site G1 fail;
- a male chrX paralog rule decides any e multiplier.

## Revision v1.5 (2026-09-30, lead)

- **§4 fallback check.** Gate B reports τ̂ bias and τ-bound coverage after the n* restriction as
  well as before it.
- If |bias| is still > 0.1·τ_max after the restriction, section-4 results carry the flag
  "inclusion-selection not corrected".

## Revision v1.6 (2026-09-30, lead)

- **Homozygote error.** It uses e_model_calibrated from B-QC r1c. These are male chrX non_artifact
  multipliers, applied because both paralog sites pass T3′. An untestable sealed cell takes its
  own cohort's development multiplier.
- **Het rules R_c.** They are re-chosen at the calibrated e_p90:
  `in-biopsy-allelic-t2b-20260930T123844Z/cohort_rules.tsv` (job 1294).
- Binding text: §9 of `2026-09-30_claude_to_codex_design_fixture_answers.md`.

## Revision v1.7 (2026-09-30, lead)

- **Refits use two starts.** Every refit runs from the warm start (the observed MLE) and the cold
  start (§7) and keeps the higher ℓ.
- **Open choices closed.** The remaining choices for the sandwich, section-4 boundary, P2
  quantiles and profile, P3 failures and inverse CDF are fixed.
- Binding text: §11 of `2026-09-30_claude_to_codex_design_fixture_answers.md`.
- **Gate B agreement so far.** Fixed vectors and the observed fit agree to ≤ 1.3e-11 in ℓ and to
  6.4e-11 in κ̂ (§10).

## Revision v1.8 (2026-10-15)

Written by the spec owner (A1_spec, model program `model-program-20261015T223347Z`). It fills the settings
that §4, §5, §7 and Appendix B3 left open. It was reviewed by R-STAT, R-RED and R-SCOPE; the item-by-item
disposition is `Analysis/MASLD_Model_Benchmark/executions/model-program-20261015T223347Z/trackA/A1_spec/REBUTTAL.md`.
It binds from the moment it is appended here. Its sha256 is recorded in that directory's `state.json`.
No Gate B dataset has been generated.

**v1.8 is a pre-Gate-B completion (user decision 2026-10-15), not the single pre-outcome amendment.**
The amendment is unused. A second Gate B failure after the amendment still stops Model A.

Abbreviations: AX = `docs/technical/agent_exchange`; EX = `Analysis/MASLD_Model_Benchmark/executions`;
P = `EX/model-program-20261015T223347Z`; QC = `EX/in-biopsy-allelic-obs-20260926T155651Z/qc`;
LAB = `EX/in-biopsy-allelic-obs-20260926T155651Z/labels`; T2B = `EX/in-biopsy-allelic-t2b-20260930T123844Z`;
TAG = `EX/in-biopsy-allelic-t2b-20260929T220019Z/tag_table/tag_table.tsv`; CG = `AX/2026-09-30_codex_gateA.md`.

**What does not change:**
- the likelihood (§2.3), the mean and dispersion (§3), and the section-4 integral ℓ_g(v);
- the three test statistics and their p-value formulas (§5), and Holm over {P1, P2, P3};
- Δ = 0.125, and the draw counts B (P1 1,000, P2 2,000, P3 1,000);
- the real-data seeds of §5;
- the rule that a failed refit is an exceedance;
- the §11 agreement tolerances;
- the B3 bands: size [0.04, 0.06], coverage [0.93, 0.97], and the full-CRT one-sided rule.

**What v1.8 adds, with the spec text it touches:**
1. The mode search, the global check and the reference integrator for adaptive GH20. The meaning of "1e-8 relative" in the §4 accuracy test (§4, spec lines 210–224).
2. A parameter bound v ∈ [0, 1], with reporting of bound hits (§3.4, line 189, which listed no v bound).
3. Step acceptance in the Newton polish and a negative-definite Hessian condition at convergence (§7, lines 343–345). v1.7 rule 8 left the polish settings to each implementation; v1.8 fixes them.
4. An extension of the v1.7 tie rule that decides which refits count as failed. This is v1.7 rule 1 (AX `2026-09-30_claude_to_codex_design_fixture_answers.md` §11); the failure accounting itself is unchanged.
5. Confirmation that the §5.2 I_c is V_c. The R_c-weighted λ is a reported sensitivity (§5.2, lines 268–274).
6. The planning τ_max values are corrected to the κ-eligible gene set that §5.2 defines (lines 293–296).
7. Gate B construction, cells, allocation, which results are judged, and a P2 calibration rule (Appendix B3, lines 489–494; §5.2, lines 285–288).
8. The §9 LD-mixture plant is named as cell 27, and the run order is fixed (§9, lines 386–388).
9. A claims sentence on how a P1 rejection is read (§3.1, line 159).

**Inputs bound by reference.** Wherever v1.8 uses e_i, the het rules R_c or ω_c, it means the values bound by the
spec revision in force when the first Gate B dataset is generated. Today that is e_model_calibrated (v1.6), the
T2B-0930 rules (`T2B/cohort_rules.tsv`) and ω_c from QC `r1c-20260930T123747Z/cohort_model_inputs.tsv`. Decision
D-QC-1 (`P/decisions_pending.md`) may change e_i and R_c before Gate B. It does not change any v1.8 text.

### 1. Mode search and integration for adaptive GH20 (§4)

1. *Tolerance unit.* In §4 and in this section, "within 1e-8 relative" means a relative error of 1e-8 in the
   gene integral exp(ℓ_g). That equals an absolute error of 1e-8 in ℓ_g, in nats. The same reading applies to the
   §4 accuracy tests.
2. *Search region.* Every row likelihood is a probability, so h_g(u) ≤ −u²/(2v). Let
   U_g = √(2v·(−h_g(0))). Every maximiser with h ≥ h(0) lies in [−U_g, U_g], and h(±U_g) ≤ h(0).
3. *Local search.* Keep a triple a < c < b with h(c) ≥ max(h(a), h(b)). At each step:
   - take a Newton step from c if h″(c) < 0 and the step lands inside (a, b);
   - otherwise take a golden-section step into the larger subinterval;
   - keep the point with the largest h as c.

   There are two starts:
   - the triple (−U_g, 0, U_g), which is always valid by item 2;
   - the previous mode û_g^prev (the mode at the previous evaluation of this gene in this fit), used only if
     h(û_g^prev) ≥ max(h(−U_g), h(U_g)).

   Stop when |Δu| < 1e-10 (the §4 value) and |h′(c)| < 1e-8·max(1, |h(c)|).
4. *Global check.* A grid scan of h on [−U_g, U_g] runs when any of the following holds:
   - (i) the two starts end at modes more than 1e-8·max(1, |u|) apart;
   - (ii) any iterate had h″ ≥ 0;
   - (iii) h″(û) > −1e-12;
   - (iv) this is the gene's first evaluation in a fit;
   - (v) this is the fit's final converged point.

   The grid has N = max(129, ⌈4U_g/σ̂_g⌉ + 1) equally spaced points, with σ̂_g taken from the best mode found so
   far. N is capped at 8,193. If the cap binds, or if σ̂_g is undefined because h″ ≥ 0 at the best mode, the
   reference integrator (item 6) is used for that gene at that evaluation. The local search restarts from every
   grid point whose h exceeds both neighbours, with the neighbours as a and b. û_g is the local maximum with the
   largest h. Ties within 1e-10 in h go to the smaller |u|, and a tie in |u| goes to the negative u.
5. *When GH20 is not used.* GH20 is not used for a gene at an evaluation if a second local maximum has
   h ≥ h(û_g) − 25, or if h″(û_g) > −1e-12. The reference integrator is used instead. The value 25 is chosen
   because a second mode with relative mass above e^−25 ≈ 1.4e-11 could move ℓ_g by a non-negligible fraction of
   1e-8.
6. *Reference integrator.* Composite Gauss–Legendre with 20 nodes per panel on [−W, W], where
   W = √(2v·(−h(û_g) + 37 + max(0, ln(√v/σ̂_g)))).
   - Truncation. Outside [−W, W], h(u) ≤ −u²/(2v) by item 2. The neglected mass is then at most about e^−37
     relative to the Laplace mass at û_g.
   - Panels. Panel edges sit at ±W and at every local mode found. The starting panel count is the smallest that
     makes every panel no wider than σ_min, where σ_min is the smallest σ̂ over the local modes found.
   - Doubling. The panel count doubles until two successive values of ℓ_g agree within 1e-10. The cap is 8,192
     panels. If 8,192 panels do not agree, the evaluation fails, and so does the fit that uses it (§7).
   - Gradient. The gradient is the exact derivative of the final composite sum with the nodes held fixed. This is
     the same rule as GH20 with (û, σ̂) held fixed.
   - Error statement. The truncation error is bounded. The panel error is checked by doubling at a panel width no
     larger than the narrowest mode found. A mode narrower than half the grid spacing can still be missed.
7. *Check against the reference.* At every evaluation where item 4 triggers, GH20 is compared with the reference
   integrator. If they differ by more than 1e-8 in ℓ_g, the reference value and its gradient are used for that
   gene at that evaluation, and the event is logged.
8. *Integrator change at the final point.* Suppose that at a fit's final converged point (trigger v) any gene's
   integrator differs from the one used in the last polish step. Then the fit is re-polished with each gene's
   integrator held at its final-point choice, and convergence (§2) is re-checked with the gradient of those
   integrators. If any gene changes integrator again at the new final point, the fit is failed.
9. *Precision.* The arithmetic precision is an implementation choice. Long double, or higher precision for the
   special functions, is allowed. The condition is that the §4 accuracy tests pass at 1e-8 in ℓ_g (item 1).
   One test is added. Take 20 of the 4,000 cell-3 (N3) datasets, chosen by
   `default_rng(SeedSequence([20261015, 3, 0, 8])).choice(4000, 20, replace=False)`. In each, every gene at
   the final fit point must agree with the reference integrator within 1e-8 in ℓ_g.
10. *Upper bound on v.* v ∈ [0, 1], so τ ≤ 1 natural-logit unit per stage step.
    - Scale. This is 15.5 times the largest planning τ_max (0.064373 at θ = 0.3, CG lines 22–27). It is 3.7
      times the planning median |α| of 0.269219 at θ = 0 (same table).
    - Reporting.
      - Every fit with v̂ within 1e-8 of 1 is listed, with its gene count. This applies in the v1.8 B-FIX family
        and in Gate B.
      - Every LR* computed on the bound is listed.
      - A P2 τ part whose τ̂ or τ0 profile fit sits on the bound is reported as not interpretable.
      - In the development run, a section-4 fit on the bound is a stated limitation.
    - Why the bound is close to tested values. REC's B-FIX native-GH refusals all occurred at v = 0.92348
      (AX `2026-10-05_codex_modelA_calibration_requirements.md`, "P3 native-mode limitation").

**Justification.** §4 names adaptive GH20 but does not say how to find the mode. That gap produced REC's
"nonconcave mode in adaptive quadrature" refusals in 25 of 50 B-FIX P3 draws. Nonconcavity of h does not make
the integral undefined; it makes a rule centred on the Laplace approximation unreliable.
- The bracket on h is always defined.
- The grid density follows the posterior width rather than the total log-likelihood.
- The reference integrator handles multimodal genes instead of refusing them.

Under v1.7 a refused evaluation was a failed fit. Under v1.8 a failure remains possible only when 8,192
panels do not converge, or when a gene changes integrator twice at the final point.

### 2. Convergence and the draw-11 tie rule (§7; v1.7 rules 1 and 8)

1. *Criterion.* A fit converges when two conditions hold:
   - the projected max |∂ℓ/∂θ_j| over the free coordinates is below 1e-6;
   - the Cholesky factorization of the negative projected Hessian over the free coordinates succeeds.

   Details:
   - At v = 0, the KKT condition is ∂ℓ/∂v ≤ 1e-6.
   - A coordinate is on a bound when |θ_j − bound| ≤ 1e-8·max(1, |bound|). On a lower bound a negative gradient
     component is set to 0; on an upper bound a positive one is.
   - Coordinates with no information are held fixed, as in v1.7 rule 3.
   - A fit that fails the Hessian condition is treated as not converged and restarted under §7.
2. *Step acceptance (Newton polish and line search).* Let ε_ℓ = max(1e-9, 1e-13·|ℓ|), and let ℓ_best be the
   highest ℓ evaluated so far in the polish. A trial step is accepted if ℓ_trial ≥ ℓ_best − ε_ℓ and the step
   lowers the projected max |gradient|. Otherwise it is rejected. Because decreases are measured against ℓ_best,
   they cannot accumulate.
3. *Start selection.* Each refit has a warm fit (the observed MLE) and a cold fit (the §7 start). Let ε_tie be
   1e-6 for section-3 fits and max(1e-6, 1e-8·G_fit) for section-4 fits, where G_fit is the number of genes with
   rows. In section-4 comparisons, both end points are evaluated with the reference integrator for every gene
   before ℓ is compared.

   | Warm | Cold | ℓ comparison | Result |
   |---|---|---|---|
   | converged | converged | differ by ≥ ε_tie | keep the higher ℓ |
   | converged | converged | differ by < ε_tie | keep warm (v1.7) |
   | one converged | the other not | converged ℓ not more than ε_tie below the other | keep the converged fit |
   | one converged | the other not | converged ℓ more than ε_tie below the other | failed refit (the unconverged fit found a higher ℓ, so no certified maximum exists) |
   | not converged after the §7 restarts | not converged after the §7 restarts | — | failed refit |

   The P2 profile values ℓ_p(τ̂) and ℓ_p(τ0) (v1.7 rule 5) are evaluated with the reference integrator for every
   gene at their final points.
4. *Record.* The v1.7 B-FIX draw results stay in the record as v1.7 results. Draw 11 is a failed supplied draw
   under v1.7. Both implementations rerun the 50 B-FIX draws under v1.8 as a new, labelled family. The §11 target
   of zero failed refits applies to that family.

**Justification.** On B-FIX draw 11, double arithmetic reports a Newton step as lowering ℓ by 8.42419e-11.
Independent 80- and 110-digit evaluations show an increase of 9.5729e-12. The projected score is 2.1271e-6 at
the original point and 2.0484e-9 at the proposed point (calibration requirements, "Independent continuation
review: saved draw11 precision"). The repaired κ̂ moved by 4.1e-8, against the §11 tolerance of 1e-4.

ε_ℓ is relative because ℓ at real scale is large in magnitude [Estimate: about 3.7×10⁴ expected included rows,
CG, at a few nats per row]. The section-4 tie band is wider than 1e-6 because the allowed integration error
summed over G ≤ 400 genes is up to 4e-6.

### 3. P2 weighting (§5.2)

**Binding.** λ_l = Σ_c V_c·λ_cond,l,c / Σ_c V_c.
- V_c ≥ 0 is cohort c's contribution to the efficient expected information for κ. It is computed after the
  nuisance parameters (α_g, r_g, δ, ω_S, ρ slopes) are projected out, with the z-adjusted stage variance, at the
  Gate A design (outcome-free).
- This is the I_c of §5.2 and of plan line 121 ("P2 weights = V_c").
- V_c is hashed and held fixed across draws, like R²_c (v1.3).

**Reported, not used for P2:**
- the R_c-weighted λ, where R_c is the first-order response of the pooled κ̂ to a slope change confined to
  cohort c, with Σ_c R_c = 1 checked within 1e-8 and R_c unclipped;
- in the development run, the V_c- and R_c-weighted λ recomputed at the development α̂.

**Reason R_c is reported.** If each cohort's slope is attenuated by λ_c and the nuisance is shared across
cohorts, the pooled MLE tends to κ_true·Σ_c R_c λ_c. The V-weighted mean equals that only when the nuisance
separates by cohort. Making R_c primary would change P2's statistic, so it is left to the user (D-SPEC-1 in
`P/decisions_pending.md`).

**Limitation.** Gate B uses the same weights in the generator and in P2, so it cannot detect error in the
weights.

### 4. Bounds on G and Z per fit (compute planning; Request 2)

| Quantity | Bound used for planning | Source |
|---|---|---|
| G, genes per fit on F | 400 | κ-eligible: 379, 364, 255, 166 at θ = 0, 0.1, 0.2, 0.3 (CG). 582 genes have positive expected inclusion (CG). |
| G, §8.1 Catalog fit on M | 830 | tag-table kept genes (TAG, kept_primary: 2,343 SNVs) |
| Z, design width | 14 (7 base + 7 missing indicators) | 9 columns used in at least one F cohort; rank 4–8 per cohort (`EX/model-a-gateA-r2-20260930T121832Z/r2_moments_observed.json`) |
| Free coordinates, section 3 | 2G + Z + 5 = 819 (section 4 adds v; κ_A adds γ_fib and γ_mye) | §7 coordinate list |
| Rows per gene | 558 (one per F individual) | protocol |
| Included rows per fit | 252,832 (all development chosen-tag candidate rows with n ≥ 1) | CG. Expected included rows 36,644.75 at the Geuvadis ρ over all genes, and 36,366.19 in the 379 eligible genes at ρ = 0.02 (CG). |
| Count axis per row | n ≤ 10,000 | T4 pileup limit (§1.2); max observed n 6,821 (Appendix A, R-STAT-17) |
| Per-gene Σ(n+1) | 558 × 10,001 = 5,580,558 | protocol |

The REC synthetic helper caps of 256 rows per gene and 100,000 count-axis entries per gene (Request 2) are below
these bounds, so that helper cannot run the real design as it stands.

[Estimate, not run] Memory: 252,832 rows × 40 nodes × 8 bytes = 81 MB of GH work arrays. 2 GB per
single-threaded fit is enough.

### 5. Gate B

#### 5.1 Data generation

- *Design.* The development F design from the v1.4 export `EX/model-a-gateA-design-20260930T121812Z` (depths,
  cohort-mixture priors, S, z, clusters), the het rules in force, and the gene set and θ chosen at Gate A. These
  files are read only by the named simulation operator. Stages S are the real values. A dataset redraws only
  genotypes, u_g, label-layer scores and reads.
- *Generating parameters.*
  - α_g is the signed, oriented GTEx ln aFC, capped at |α| ≤ 3 (§3.4).
  - r_g = logit(0.033805), the Geuvadis het ρ (`T2B/t2b_summary.json`, geuvadis_fit.rho_het), and
    φ = 0.043873.
  - ω_c and e_i are the values in force.
  - δ, ω_S and the ρ slopes are 0 unless a cell sets them.
- *Planning τ_max.* The planning τ_max is Δ × the median GTEx |α| over the κ-eligible genes, as §5.2 defines it:
  0.033652, 0.034934, 0.048558, 0.064373 at θ = 0, 0.1, 0.2, 0.3 (CG lines 22–27). The values in §5.2's
  parenthesis (0.0498–0.0765) are medians over 1,290 public eGenes (AX `2026-09-29_codex_c5_public_afc.md`,
  lines 9–12). They are not used. The real τ_max is still computed from the development α̂ (B1).
- *Order of generation.* u_g is drawn first (5.1a), then the genotype class from P(G), then a ~ f_G(· | n). The
  het predicate is applied last, and rows that fail it are dropped. So inclusion depends on u_g.
- *Fitting.* Generating values are never used as starts. Each dataset's observed fit uses the §7 cold start with
  the §7 restarts. Refits use the dataset's own MLE (warm) and the cold start, as v1.7 and §2 above require. The
  stage law is refitted on each dataset.
- *Seeds.* All draws use `default_rng(SeedSequence([20261015, cell_id, r, stream, k]))`. Here r = 0, 1, … is the
  dataset and k = 0, 1, … is the root index, omitted when there is one draw.

  | Stream | What it draws |
  |---|---|
  | 0 | genotypes, u_g and reads |
  | 1 | label layer: the λ-draw index and the X noise |
  | 2 | CRT S* |
  | 3 | choice of the same-row P2 root (5.2) |
  | 4 | P3 bootstrap |
  | 5 | plant-specific draws |
  | 6 | P2 τ-fallback bootstrap |

  The full-CRT replicates use cell_id 99. Replicate r's data come from stream 0, and its inner draws
  b = 0, …, 999 use `SeedSequence([20261015, 99, r, 2, b])`. The Monte Carlo CIs use
  `SeedSequence([20261015, cell_id, 0, 9])`.

#### 5.1a Projection of the planted gene departures

In every cell with τ > 0 (cells 2, 3, 99, the τ > 0 power cells, 15 and 43), u is drawn iid N(0, τ²) and then
replaced by

  u ← u − α·(Σ_g α_g J_g u_g / Σ_g α_g² J_g),

where J_g is gene g's expected information for its slope on S̃ after α_g and r_g are projected out. J_g is
computed from the generating design at the generating parameters, which is outcome-free. The simulation operator
exports the scalar Σ_g α_g² J_g² / Σ_g α_g² J_g, which is an aggregate.

**Why.** κ_U is defined for the fixed set of eligible genes (§3.1; Appendix A, R-STAT-8). With iid u_g, the
section-3 target in dataset r is κ + Σ α J u / Σ α² J, which is not κ. The individual-cluster sandwich does not
include that between-gene term. [Estimate, from planning weights, R-RED] The variance inflation is
τ²·Σα²J²/Σα²J ≈ 0.89 at τ = 0.0498, which would give a P1 rejection rate of about 0.15 in N3. After the
projection, κ_U,r equals the planted κ to first order. The projection removes one direction from a set of up to
400 genes, so it changes the τ̂ bias target by a relative amount of order 1/G [Inference].

#### 5.2 Warp-speed construction

The construction follows Giacomini, Politis and White: each simulated dataset gives K resampled roots, and the
roots are pooled across the datasets of a cell. K is set per cell in 5.3.

- **P1.**
  - Dataset r gives z_r from the observed fit. It also gives K roots z*_rk: one CRT draw S* each (stream 2, root
    k) from the dataset's fitted stage law, each followed by a refit.
  - c₁ is the 0.95 quantile (numpy linear rule) of the R·K values |z*_rk|.
  - The dataset rejects when |z_r| − 1e-6 > c₁. This mirrors the §5.1 tie band.
  - Warp size = (number of rejecting datasets) / R.
  - A failed observed fit counts as a rejection. A failed CRT refit sets |z*_rk| = +∞.
- **P3.**
  - Dataset r gives LR_r and K roots LR*_rk (§5.3 conditional bootstrap at the dataset's section-3 fit,
    stream 4).
  - c₃ is the 0.95 quantile of the LR*_rk. The dataset rejects when LR_r − 1e-6 > c₃.
  - Failures are handled as for P1.
- **P2 κ part.**
  - *Roots.* Root k of dataset r uses row b(r, k) = (r·K + k) mod 2,000 of the real resample file
    `EX/model-a-gateA-r2-20260930T121832Z/p2_resample_development.tsv.gz`. The design individuals are the same
    in every dataset. The root is ε_rk = (κ*_rk − κ̂_r)/SE_r, where SE_r is the §5.1 sandwich SE of the observed
    fit. Each root carries the label b(r, k).
  - *Pairing with λ.* For each dataset r and each b = 0, …, 1,999, one root ε is chosen from the pooled roots
    labelled b (stream 3). Then θ*_rb = (κ̂_r + SE_r·ε)/λ*_b, where λ*_b is the real V_c-weighted λ draw b for
    that source. This keeps the §5.2 pairing of κ's resample b with the λ draw built on that resample.
  - *p-values.* p_up = (1 + #{θ*_rb ≥ Δ})/2,001 and p_lo = (1 + #{θ*_rb ≤ −Δ})/2,001.
  - *Failures.* A failed root, or λ*_b ≤ 0.05, is an exceedance in both p_up and p_lo.
  - *Why exact λ.* The λ draws do not depend on the simulated reads, so all 2,000 are applied exactly. Only the
    studentized κ root is pooled; it is approximately pivotal (R-RED M2).
- **P2 τ part.**
  - p_τ,r is computed exactly in each dataset (§5.2; v1.7 rule 5), with τ0_r = Δ·median_g|α̂_g,r|·λ^(0.05).
  - α̂_g,r comes from the dataset's own section-3 fit. λ^(0.05) is computed by the v1.7 rule 5 quantile over the
    V_c-weighted real λ draws.
  - A dataset's P2 rejects when max(p_up, p_lo, p_τ) ≤ 0.05 for both sources.
- **P2 τ fallback (spec §5.2), warp form.**
  - Each dataset gives 2 roots r*_rk from the conditional parametric bootstrap of r at τ = τ0_r (stream 6),
    using the §5.3 machinery.
  - The pooled 0.95 quantile c_τ of the r*_rk replaces 1 − Φ(r). The τ part rejects when r_r − 1e-6 > c_τ.
- **Monte Carlo uncertainty.**
  - Each warp size and coverage gets a 95% CI from 2,000 resamples of the datasets, with each dataset's statistic
    and roots kept together.
  - The bands are judged on the point estimate.
  - The pooled quantile adds to the binomial variance [login-node MC: numpy `default_rng(20261015)`, 800 reps,
    N(0,1) roots]. At R = 4,000 and a true size of 0.05:

    | K | SD of warp size | P(outside [0.04, 0.06]) |
    |---|---|---|
    | 1 | 0.0047 | 0.037 |
    | 2 | 0.0042 | 0.014 |
    | 4 | 0.0039 | 0.013 |

#### 5.3 Cells and allocation

Λ0 means no label error (X = S). Λ1–Λ3 are defined in §8. "Projected" means the 5.1a projection.

| id | Cell | κ | τ | Label | Datasets R | K | Tests | Judged |
|---|---|---|---|---|---|---|---|---|
| 1 | N1 null | 0 | 0 | Λ0 | 4,000 | P1 4; P3 2 | P1 size, κ coverage, P3 size | yes |
| 2 | N2 null | 0 | τ_max/2, projected | Λ0 | 4,000 | P1 4; P3 1 | P1 size, κ coverage; P3 power, τ̂ bias | P1 size and coverage |
| 3 | N3 null | 0 | τ_max, projected | Λ0 | 4,000 | P1 4; P3 1 | as N2 | P1 size and coverage |
| 40 | N4 nuisance null | 0 | 0 | Λ0 | 4,000 | P1 4 | P1 size, κ coverage | yes |
| 41 | N4h | as N4 with α halved | 0 | Λ0 | 1,000 | P1 1 | P1 size, κ coverage | no |
| 42 | N3u | 0 | τ_max, not projected | Λ0 | 1,000 | P1 1 | P1 rejection caused by heterogeneity alone | no |
| 4–12 | power | −0.10, −0.05, +0.05 | 0, τ_max/2, τ_max (projected) | Λ0 | 1,000 each | P1 1; P3 1 where τ > 0 | P1 power, κ coverage; P3 power and τ̂ bias where τ > 0 | no |
| 13 | B1 | κ_true = +Δ | 0 | Λ1 | 4,000 | P2 4 | P2 κ part; full P2 on r < 1,000 | κ part (P2 calibration, 5.4) |
| 14 | B2 | κ_true = −Δ | 0 | Λ1 | 4,000 | P2 4 | as B1 | κ part (P2 calibration) |
| 15 | B3 | 0 | τ_true = τ_max, projected | Λ1 | 1,000 | — | P2 τ part | no |
| 43 | B3′ | 0 | τ = τ0, projected | Λ0 | 2,000 | fallback 2 | P2 τ part | yes (P2 calibration; fallback trigger) |
| 16 | B4 | κ_true = +Δ | 0 | Λ3 | 1,000 | P2 2 | full P2 | no |
| 17, 18 | P2 power | κ_true = 0, −0.05 | 0 | Λ1 | 1,000 each | P2 2 | full P2 | no |
| 19–21 | L1–L3 | κ_true = −0.10 | 0 | Λ1, Λ2, Λ3 | 1,000 each | P1 1; P2 2 | P1, P2 interval coverage | no |
| 22–30, 33 | technical (§7) | 0 | 0 | Λ0 | 1,000 each | P1 1; P3 1 | size, κ bias | no (limitation rule) |
| 27 | T5, the §9 plant | 0 | 0 | Λ0 | 1,000 | P1 1; P3 1 | κ bias | §9 rule |
| 31, 32 | T5, T6+ at signal | −0.10 | 0 | Λ0 | 1,000 each | P1 1 | κ bias | no |
| 34, 35 | T2, T4 at the threshold amplitude | 0 | 0 | Λ0 | 1,000 each | P1 1 | measured bias / predicted bias | no |
| 99 | full CRT | 0 | τ_max, projected | Λ0 | 200 replicates × B = 1,000 | — | P1, full procedure | yes |

Totals:
- 58,000 warp datasets, of which 39,000 need section-4 fits, plus 200 full-CRT replicates.
- In B3′, τ0 = τ_max·λ^(0.05), using the generating τ_max and the real λ draws. B3′ is the cell "at
  τ = λ·τ_max" in §5.2's fallback rule, with λ read as λ^(0.05), the value that defines the τ-part boundary.
- B3 sits at about 1.7·τ0 [Inference, R-RED M5], inside H0, so it cannot trigger the fallback. It is reported.

**N4 settings** (all outcome-free):
- *Nuisance.*
  - ω_S = U95 (§7, T1).
  - ρ_S = +0.2, ρ_b = +0.2 and ρ_d = −0.2. These are chosen amplitudes in logit-dispersion units per unit of the
    centred column. No measured source exists for them.
  - δ = +0.1 for b, −0.1 for d, +0.1 for FFPE, −0.1 for sex, +0.1 for age/10, and 0 for the ancestry axes.
- *Dispersion.* r_g = logit(0.033805) + N(0, s_r²) per gene. s_r is the SD of per-gene logit ρ over Geuvadis genes
  with ≥ 50 DNA-het rows (B-QC measures it on public Geuvadis rows).
- *Effect size.* α_g = β_LD,g·GTEx α_g. β_LD,g is the §8.2 tag–lead attenuation in the cohort superpopulation
  mixture (`tag_ld_by_superpop.tsv`), computed by the simulation operator.
- *Genotypes.* Drawn from each individual's own superpopulation frequencies (the §10 individual-to-reference
  map). The fit keeps the merged cohort-mixture prior.

In N1 the reads depend on neither S nor z, so N1 cannot test the CRT's reliance on the estimated stage law. N4
can.

#### 5.4 Judgment

**Gate B criteria (a pass requires all of them).** These are the B3 criteria, applied to the judged null cells:
1. P1 warp size in [0.04, 0.06] in cells 1, 2, 3 and 40, with failures counted as in 5.2.
2. 95% κ coverage in [0.93, 0.97] in cells 1, 2, 3 and 40 (§6).
3. P3 warp size in [0.04, 0.06] in cell 1.
4. Cell 99 fails if the one-sided Clopper–Pearson lower 95% bound on its rejection rate at 0.05 exceeds 0.05.

**Prespecified rules that do not use the amendment** (spec §4, v1.5, §5.2, §9; plan A7):
- *τ̂ bias.* τ̂ bias in cells 2, 3 and the τ > 0 power cells, before and after the n* restriction. The §4
  fallback applies if |bias| > 0.1·τ_max. If |bias| is still above that after the restriction, the flag
  "inclusion-selection not corrected" applies.
- *§9 δ_n rule.*
  - Cell 27 runs before any judged cell.
  - If its |κ bias| > Δ/10 = 0.0125, δ_n moves into the primary model, and every judged cell (1, 2, 3, 13, 14,
    40, 43, 99) then runs with δ_n in the fit.
  - Cell 31 is reported and does not trigger the rule.

**P2 calibration rule (bound now as the conservative option; D-SPEC-2).** P2 cannot pass (p₂ = 1 in Holm, and P2
is reported as "not calibrated") if any of the following holds:
- the B1 or B2 κ-part rejection rate exceeds 0.06;
- the B3′ τ-part rate exceeds 0.06 and the fallback τ part, rerun on B3′, also exceeds 0.06.

If B3′ exceeds 0.06 and the fallback does not, the fallback becomes P2's τ part (§5.2). None of these outcomes
fails Gate B for P1 and P3, and none uses the amendment. Under intersection–union, P2's rejection rate is at most
its κ part's rate, so judging the κ part is the stricter check.

**False-fail probability.** These figures assume a procedure that is exactly calibrated, and independent cells
[Inference].
- *Gate B as a whole:* 1 − (1 − 0.013)⁴·(1 − 0.014)·(1 − 0.044) ≈ 0.105. The terms are four P1 bands at K = 4,
  one P3 band at K = 2, and the cell-99 rule at a true size of 0.05 (exact binomial). Coverage bands add < 1e-7
  each (binomial at R = 4,000).
- *P2 calibration rule:* 1 − (1 − 0.005)²·(1 − 0.020) ≈ 0.030. The 0.005 per B cell is R-STAT's Monte Carlo
  figure, and 0.020 is the exact binomial P(rate > 0.06) at R = 2,000.

#### 5.5 Reported, not judged

- *Rates and coverage.* Power. Coverage in the 1,000-dataset cells, with binomial CIs. The asymptotic rejection
  rate #{|z_r| ≥ 1.959964}/R.
- *Failures and runtime.*
  - The observed-fit and refit failure rates per cell. The A6 pilot already requires failures below 1% before
    the full run (plan line 127).
  - The warp size with failed pairs excluded, beside the size with failures counted.
  - Runtime per fit.
- *Diagnostics.* Per-gene LRT QQ plots from N1. λ_gen,cond (§8) for every B and L cell. The scalar of 5.1a.
- *B4.* It sits at the upper 95% bound of var_e, outside the averaged-over-λ guarantee P2 gives, so a valid P2 is
  expected to reject there more often than 0.05 [Inference, R-STAT M4 and R-RED M4]. Its rate is a named
  limitation of P2.
- *Limitation rule, technical cells 22–30 and 33 (except T6).* A cell is named as a limitation in the paper if
  either of these holds:
  - its P1 or P3 size has a one-sided lower 95% bound above 0.05;
  - its |κ bias| exceeds Δ/10.

  That cell's §9 real-data diagnostic then becomes the stated guard. A judged null cell whose size with failed
  pairs excluded is above 0.06 is also named as a limitation.

#### 5.6 Full-CRT cell

All 200 replicates are in cell 99 (κ = 0, τ = τ_max, projected). With u_g ≠ 0 the reads depend on S through u_g.
That is where the conditional distribution of the root varies most across datasets, so it is the right place to
check the pooling. The resolution is limited: the rule fails at ≥ 16 of 200 rejections. P(fail) is 0.044 at a
true size of 0.05, 0.149 at 0.06, 0.537 at 0.08 and 0.857 at 0.10 (exact binomial). A pass therefore means only
that the size is not above about 0.08.

#### 5.7 Cost and pause rule [Estimate, not run]

These figures use the §7 fit-time estimates: section-3 fit 5–60 s, section-4 fit 2–25 min.

| Part | CPU-h |
|---|---|
| All warp cells | 6,500–80,700 |
| Full CRT (200 × 2,001 section-3 fits) | 560–6,670 |
| Total | 7,000–87,400 |
| Judged cells plus cell 27 plus the full CRT | 3,400–42,600 |
| B3′ fallback, if triggered | 530–6,670 |
| δ_n rerun of the judged cells, if the §9 rule fires | about the judged figure again |

The Track A Gate B cap is 80,000 CPU-h, with a pause at 80%, which is 64,000 (`P/budget.json`; plan line 227).
- The A6 pilot's measured fit times replace these estimates before submission.
- Cells are submitted in this order: cell 27, the judged cells and cell 99, then the reported cells.
- If the projected total at measured times exceeds 64,000 CPU-h, only the judged tier is submitted, and the
  user is asked before the reported cells run.

### 6. Coverage

- **The interval judged.** κ̂ ± 1.959964·SE, from the section-3 fit and the §5.1 cluster sandwich SE. It is judged
  against [0.93, 0.97] in cells 1, 2, 3 and 40, for the planted κ. After the 5.1a projection that equals κ_U,r to
  first order. A failed observed fit counts as not covering.
- **Λ0 cells only.** Under label error the target is the probability limit of κ̂ on observed labels, which has no
  closed form, so the κ interval is not judged there.
- **P2's interval.** The 90% bootstrap interval of κ/λ, formed as in 5.2 from θ*_rb. Its coverage of κ_true is
  reported in cells 13, 14, 16 and 19–21, per source, with a binomial 95% CI. Only the source that matches the
  generator is compared with the nominal 0.90: GSE193066 in Λ1 and Λ3, PXD051911 A1 in Λ2. It is never compared
  with [0.93, 0.97].

### 7. Technical plants

All plants enter the generator only; the fitted model is the v1.8 primary model. The simulation operator
computes the first-order predicted κ bias and size of every technical cell, outcome-free from the Gate A design
and the A3 information summaries. These predictions go into the Gate B parameter file before any dataset of the
cell is generated.

| id | Plant | Setting | Measured input |
|---|---|---|---|
| 22 | T1a: reference-bias slopes of opposite sign by cohort group | H mean adds s_t·ω_c,plant·S̃_i, with ω_c,plant = +2·U95 in GSE213621 and GSE240729 (the transition-tag-exclusion cohorts) and −2·U95 in the other three. The fitted shared ω_S absorbs only the weighted mean. | U95 = max(\|lo\|, \|hi\|) of the 95% CI of the inverse-variance pooled within-cohort slope of per-library null-site REF logit on S̃ + z (§9 negative control 1). **Not yet measured; B-QC computes it on QC null sites.** |
| 30 | T1b: per-tag random reference-bias slope | ω_t ~ N(0, (2·U95)²) for each tag, drawn once per dataset (stream 5) | as 22 |
| 33 | T1c: slope on transition tags only | +2·U95 on transition tags in the three cohorts that keep them | as 22 |
| 23 | T2: false-het rate doubled at S = 2 | e_gen,i = 2·e_i if S_i = 2, else e_i. The fit keeps e_i. | Current values (v1.6), e_model_p50 × multiplier ×2: GSE130970 0.002182 → 0.004363; GSE135251 0.006128 → 0.012256; GSE162694 0.017266 → 0.034532; GSE213621 0.002206 → 0.004412; GSE240729 0.003020 → 0.006040 (QC `r1c-20260930T123747Z/cohort_model_inputs.tsv`). If D-QC-1 is adopted, e_i = e_model and the doubling applies to that. |
| 24 | T3: non-BB tails with depth–stage coupling | Depth is thinned first, n′ ~ Binomial(n, 0.8^{S_i}). Then each candidate row, independently with probability w_tail, is drawn with ρ_tail = m·ρ_ig. | Plan line 121 sets w_tail = 0.05 and m = 5. v1.8 uses w_tail = max(0.05, w_hi) and m = max(5, m_hi). w_hi and m_hi are the upper 95% bounds of a two-component BB fit to the 349,386 Geuvadis DNA-het rows (n 10–1000; `T2B/t2b_summary.json`), bootstrapped over the 30 individuals. Base ρ_ig = max(0.033805, the cohort's het ρ at QC null sites). **B-QC measures both before Gate B.** |
| 25, 26 | T4±: tag switching | A row is switched with probability π(S) = 0.02·S (T4+) or 0.02·(2 − S) (T4−). The switched row's genotype class and oriented mean are drawn by the LD-mixture generator, using the switched tag's r² and frequencies in the individual's superpopulation mixture. | ±0.02 per step is the §9 tag-switching band. Applies to genes with ≥ 2 kept tags: 463 of 830 tag-table genes (TAG), and 284 of the 379 planning-eligible genes [R-RED count]. r² < 0.8 for AFR 1,216, AMR 167, SAS 166, EAS 3 and EUR 0 of the 2,343 kept tags; 6 AFR sign flips (TAG). |
| 27 | T5: LD mixture, the §9 plant | Tag hets come from the coupling / lead-hom / repulsion mixture implied by each tag's r² and frequencies in the cohort superpopulation mixture, with depth thinned 0.8^{S_i}, as §9 states. | `tag_ld_by_superpop.tsv`, `tag_2pq_by_superpop.tsv` (B-Ancestry, public 1kGP) |
| 28, 29 | T6±: omitted depth-dependent bias | The generator puts δ_n·ℓn_ig inside the parentheses of η, where ℓn_ig is log n centred within gene × cohort. Depth is thinned 0.8^{S_i}. δ_n = ±0.11204. The fit omits the term. | Derived. Thinning lowers log n by 0.22314 per step, so the first-order κ shift is −0.22314·δ_n = ∓0.025, which is 2·Δ/10 [Inference]. "δ_n at the §9 bound" is the plan's phrase (line 121); the spec does not define it. The planted value is twice the δ_n at which the first-order shift equals the §9 tolerance (0.0125/0.22314 = 0.05602). |
| 31, 32 | at signal | cell 27, and cell 28 at κ = −0.10 | as above |
| 34 | T2 at the threshold amplitude | The S = 2 multiplier at which the first-order κ bias is Δ/10, computed by the operator before generation | as 23 |
| 35 | T4 at the threshold amplitude | The switching slope at which the first-order κ bias is Δ/10 [Estimate, R-RED: about 0.12 per step] | as 25 |

**Notes.**
- *T6.* The 0.8^S thinning is a stress amplitude chosen to measure the first-order map. It is not a tolerance for
  the real data, whose depth–stage coupling the §9 diagnostic measures (median |Spearman| of chosen-tag n with S
  below 0.1). T6 reports the measured κ shift divided by −0.22314·δ_n, and it is not subject to the limitation
  rule. The real-data guard is the §9 diagnostic together with the §9 δ_n sensitivity.
- *Threshold cells.* Cells 34 and 35 check the first-order map at the threshold. They are reported as the ratio of
  measured to predicted bias, and they do not feed the limitation rule.
- *T5 threshold.* No threshold cell is defined for T5. Its attenuation is equal at every stage, so its first-order
  κ bias is 0.
- *Earlier T1 form.* The v1.8 draft planted the shared slope s_t·ω_S·S̃, which the fitted model contains. That
  cell could not move κ̂ [Estimate, R-RED: κ–ω_S information correlation −0.028], so it is dropped.

### 8. Label-error and λ scenarios

**Generator (Λ1–Λ3).** The real S_i is kept. Within cohort c, let m_c(z) be the linear projection of S on z, with
coefficients from the XtX and XtS moments in `EX/model-a-gateA-r2-20260930T121832Z/r2_moments_observed.json`
(pseudo-inverse at the cohort rank). Let s²_c be the residual variance from the same moments. Then

  X_i = m_c(z_i) + λ_c·(S_i − m_c(z_i)) + N(0, λ_c·(1 − λ_c)·s²_c),

where λ_c is the λ_cond,c of the scenario, clipped to [0.01, 0.99] for generation only. The number of clipped
cohort-draws is reported.

**Properties.**
- Cov(X, S | z)/Var(S | z) = λ_c, so the z-adjusted attenuation equals the λ_cond that P2 divides by.
- Var(X | z) = Cov(X, S | z), which is the identity that X = E[S | T] satisfies.
- X̃ (X centred within cohort) replaces S̃ in η: η = α_g(1 + κ_true·X̃ + δᵀz) + u_g·X̃. The fit sees only S.
- The realized λ_gen,cond is printed for every dataset of cells 13–21 and summarised per cell. It is the
  V_c-weighted within-cohort partial slope of X̃ on S̃ given z.

| Scenario | λ_c used by the generator | Planning value (headcount-pooled F) |
|---|---|---|
| Λ0 | X = S | 1 |
| Λ1 | for each dataset, one draw j (stream 1) from the carried GSE193066 arm's full-run draws; λ_cond,c,j for each cohort at the v1.3 fixed R² | 0.5754, CI 0.3050–0.7261 (v1.4; AX `2026-09-29_codex_lambda_check.md:62`) |
| Λ2 | as Λ1, from the PXD051911 A1 draws | 0.7355, CI 0.5639–0.9725 (same file, line 63) |
| Λ3 | λ_cond,c at var_e fixed at the upper 95% bound of the carried GSE193066 arm's full-run var_e draws | about 0.19 [Inference: linear interpolation between smoke draws 14 (var_e 0.7785, λ_cond 0.2198) and 13 (0.8292, 0.1677) at 0.805105, LAB `smoke-20260930T121336Z/lambda_source_draws.tsv.gz`]. The binding value is computed by the label layer from the full run and written into the Gate B parameter file before cell 16 runs. |

**λ inside the simulated P2.** P2 uses the real B-LABEL λ draw file (2,000 draws per source, V_c-weighted). The
labels are real and fixed, so the same draws serve every dataset (5.2). Two limits on what the B cells show:
- In Λ1 and Λ2 the true λ is itself a draw from the distribution P2 uses. So cells 13, 14 and 17–19 check P2 on
  average over that distribution, not at every λ value in it. B4 (Λ3) shows the tail, and it is reported.
- κ_true = ±Δ in cells 13 and 14 sits on P2's boundary to the §5.2 approximation (κ̂ ≈ Σ_c V_c κ̂_c / Σ_c V_c),
  averaged over the λ draws.

**Prerequisites before cells 13–21 and 43 run:**
1. The B-LABEL full run (B = 2,000) at the v1.4 R² is done. It needs D-LABEL-1 (read of the frozen crosswalk).
2. V_c is available from A3 (C5 helper on F558), together with the identity qualification.
3. D-LABEL-2 (whether GSE193066 is carriable) and D-LABEL-3 (PXD051911 A1 at the S cutpoints) are decided. If the
   GSE193066 arm is not carried, P2 cannot pass (user decision B2, 2026-10-15), and cells 13–18 and 43 are
   reported as uninformative for the paper's P2.

### 9. Claims text for P1 (§3.1)

P1 tests κ_U for the fixed set of eligible genes. A P1 rejection means one of two things:
- the information-weighted average proportional change in cis effect across those genes differs from 0;
- the gene-specific departures u_g have an average correlated with α_g.

P1 cannot tell these apart. P3 and the per-gene results are read beside P1. The paper does not describe a P1
rejection as a change shared by every gene.

### 10. Request 2 (AX `2026-10-05_codex_independent_data_requests.md`)

**Answered by v1.8:**
1. *P2 information weighting:* V_c, as in §3. The R_c-weighted λ is a sensitivity.
2. *Allocation bounds for G, Z, rows per gene, maximum depth and per-gene Σ(n+1):* §4.
3. *Warp-speed construction, allocation, and how the 200 full-CRT replicates relate:* §5.
4. *Which interval and endpoint get the coverage band:* the 95% sandwich-Wald κ interval, in Λ0 null cells (§6).
5. *The 90% P2 interval:* reported against 0.90, per source, never against [0.93, 0.97] (§6).
6. *Label-error and λ scenarios:* §8.
7. *Reference-bias slope, non-BB tails, tag switching and omitted depth bias:* §7. The amplitudes that need
   measurement are U95, the Geuvadis tail fit and the null-site het ρ. Their measurement rules are fixed here.
8. *Confirmation of B3:* size [0.04, 0.06] and coverage [0.93, 0.97] apply at ≥ 4,000 warp datasets per judged
   null cell. The κ × τ grid is cells 1–12.
9. *Draw 11:* it remains a failed supplied draw under v1.7. The v1.8 rule and the v1.8 rerun family are in §2.
10. *Prerequisites for opening the seal, all required:*
    - v1.8 appended and hashed;
    - B-FIX 50-draw agreement under v1.8 within §11;
    - the Gate A decision;
    - the Gate B pass;
    - identity qualification of F558 and F237;
    - the §10 hash list complete, including the sealed technical inputs and the sealed MDE written before opening;
    - explicit user authorization.

    An observed fixture fit or null-plant information does not establish readiness.

**Needs identity, label or C5 outputs; v1.8 cannot answer these:**
- the outcome-free identity qualification of F558 and F237;
- how identity and source-label exclusions apply to the GSE193066 λ inputs;
- whether exclusions change cohort sizes, centring or cluster resampling;
- the owner-reviewed C5 summaries (I, V_c, R_c, checks) for the eight null-plant × θ cases;
- the sealed F237 readiness attestation.

### Open inputs that block parts of Gate B

1. D-QC-1 (e_i multiplier and het rules). It must be decided before the first Gate B dataset is generated, because
   it fixes the values "in force". It does not change v1.8.
2. B-QC measurements:
   - U95, for cells 22, 30, 33, 40 and 41;
   - the Geuvadis two-component BB tail fit, for cell 24;
   - the null-site het ρ per cohort, for cell 24;
   - the Geuvadis per-gene ρ SD, for cells 40 and 41.
3. The B-LABEL full run, with D-LABEL-1, D-LABEL-2 and D-LABEL-3. These block cells 13–21 and 43.
4. V_c from A3, and the identity qualification. These block the P2 cells and Gate A.
5. The simulation operator's outcome-free quantities:
   - J_g and the 5.1a scalar;
   - β_LD,g;
   - m_c and s²_c;
   - the predicted biases of the technical cells;
   - the threshold amplitudes for cells 34 and 35.
6. The A6 pilot's measured fit times. They replace the cost estimate in 5.7.
7. User decisions D-SPEC-1 to D-SPEC-3 (`P/decisions_pending.md`). v1.8 binds the conservative option for each.
   If the user chooses otherwise, that is a further pre-Gate-B revision and must be made before the first Gate B
   dataset.

## Addendum v1.8a (2026-10-15)

Written by the Model A spec owner (model program `model-program-20261015T223347Z`, stage `trackA/A1b_addendum`).
Status: final. Reviewed by R-STAT and R-SCOPE (both "approve", 15 minor items); the item-by-item answers are in
`P/trackA/A1b_addendum/REBUTTAL.md`. Items that need the lead or the user are listed in `P/decisions_pending.md`
under "Model A addendum v1.8a".

**What this addendum does.** It records the lead rulings of `P/trackA/LEAD_RULINGS_W1_2026-10-15.md` (sha256
`9ceb724f…a7f2f`) as the settings in force before the first Gate B dataset is generated. No Gate B dataset has
been generated, and no Model A tag allelic ratio has been read.

**What it does not do.** It changes no threshold, no Δ (0.125), no τ_max rule, no test statistic, no p-value
formula, no band, no draw count, no seed and no gate rule. The Holm family stays {P1, P2, P3}. It is not the single
pre-outcome amendment, which remains unused.

**Why it is not the amendment.** v1.8 itself allows the D-QC-1 change of e_i and R_c (v1.8 lines 37–40) and the
D-SPEC settings (v1.8 lines 590–592) as pre-Gate-B revisions. The other changes here (cells not run, cells 50–52
added, the D-QC-3 placement) follow from the A2b outcome under the user's 2026-10-15 decision that P2 cannot pass
if the GSE193066 λ arm cannot be carried. No cell that is dropped or added enters Gate B criteria 1–4 (v1.8 lines
367–370). v1.8a goes to the user with the lead rulings, and the user may override any part of it before the first
Gate B dataset is generated.

The base text is the spec `AX/2026-09-29_model_A_likelihood_spec_v1.md`. Its sha256 is `5beec1f7…efc6` (1,182
lines), and v1.8 starts at line 591. Lines are cited as "v1.8 line N", counted in `P/trackA/A1_spec/v1.8_final.md`
(sha256 `f8d4f1d4…61ac`). The same line is spec line N + 590.

Abbreviations follow v1.8: AX = `docs/technical/agent_exchange`; EX = `Analysis/MASLD_Model_Benchmark/executions`;
P = `EX/model-program-20261015T223347Z`; QC = `EX/in-biopsy-allelic-obs-20260926T155651Z/qc`.

Two further abbreviations:
- T2B-0929 = `EX/in-biopsy-allelic-t2b-20260929T220019Z`;
- T2B-0930 = `EX/in-biopsy-allelic-t2b-20260930T123844Z`.

---

### 1. D-QC-1: e_i and het rules in force

**Binding.**
- Primary e_i is e_model. The male chrX multiplier is 1, so e_model_calibrated = e_model.
- The values are e_model from QC `r1c-20260930T123747Z/cohort_model_inputs.tsv` (sha256 `5d3b01d4…0f2`). They are
  per library where `e_i_per_library_allowed` = True, which among F cohorts is GSE240729 only. Elsewhere they are
  the cohort value.
- Het rules R_c are `T2B-0929/cohort_rules.tsv` (sha256 `6835962b…6326`).
- The v1.6 multiplied e (e_model × `e_multiplier_applied`, same file) and the T2B-0930 rules
  (`T2B-0930/cohort_rules.tsv`, sha256 `1621eab6…f456f372`) become a named Gate B technical sensitivity (cells 50
  and 51, §2 below).
- This replaces the "Inputs bound by reference" paragraph of v1.8 (lines 37–40) for e_i and R_c. ω_c is unchanged.

**Values in force for F (development).** e_model_p50 comes from the r1c file above, and the rules
(min_dp, min_minor_reads, min_minor_pct) from `T2B-0929/cohort_rules.tsv`.

| Cohort | e_model (p50) | v1.6 multiplier, now a sensitivity | Rule T2B-0929 | Rule T2B-0930, now a sensitivity |
|---|---|---|---|---|
| GSE130970 | 0.0010601016 | 2.0579 | 10, 1, 10 | 10, 1, 10 |
| GSE135251 | 0.0009635777 | 6.3597 | 10, 1, 10 | 10, 3, 15 |
| GSE162694 | 0.0056017931 | 3.0822 | 10, 2, 15 | 20, 5, 20 |
| GSE213621 | 0.0008767901 | 2.5159 | 10, 1, 5 | 10, 1, 10 |
| GSE240729 | 0.0030200136 (per library; p10–p90 0.0019595–0.0054598) | 1.0 | 10, 1, 10 | 10, 1, 10 |

The sealed rows of `T2B-0929/cohort_rules.tsv` are the sealed rules in force. They are GSE130970 (10, 1, 5),
GSE135251 (10, 1, 10), GSE162694 (10, 2, 15), GSE174478 (10, 1, 10), GSE213621 (10, 1, 10) and GSE240729
(10, 1, 10). Across all 14 cohort × set rows, 8 differ between the two rule files (count from the two files). This
matches `P/decisions_pending.md` D-QC-1.

**Verified: the Gate A design export already uses these settings.**
- *Rules file.* `scripts/analysis/in_biopsy_allelic/tags/run_gateA_export.sbatch` line 11 sets
  `T2B=$EX/in-biopsy-allelic-t2b-20260929T220019Z`, and line 23 passes `--rules $T2B/cohort_rules.tsv`.
- *Rules used.* The rules recorded per cohort in `EX/model-a-gateA-design-20260930T121812Z/design_cohort_moments.json`
  equal the T2B-0929 rows in all 11 cohort × set cells. For example, they are (10, 1, 10) for GSE135251 development,
  (10, 2, 15) for GSE162694 development and (10, 1, 5) for GSE213621 development. They differ from T2B-0930 in 7 of
  the 11 cells: GSE135251, GSE162694 and GSE213621, development and sealed, and GSE130970 sealed. (The eighth
  differing row of the rule files, PRJNA512027 development, is not one of the 11 design cells.)
- *Script identity.* The export script `gateA_design_export.py` has sha256 `b720aa74…3358`. That equals the hash
  stored in the export's `sha256_script.txt`.
- *z_b.* The script computes z_b = log max(e_model, 1e‑4) (line 125, `np.log(ind["e_model"].clip(lower=1e-4))`).
  It reads `e_model`, not a multiplied column.
- *z_b values.* exp(z_b) in `design_individuals.tsv` equals the r1c e_model_p50 to 10 digits in every
  single-valued cohort, and GSE240729's median equals 0.0030200136. Examples are GSE130970 0.0010601016, GSE135251
  0.0009635777 and GSE162694 0.0056017931. No multiplier enters z_b.
- *R² moments.* `EX/model-a-gateA-r2-20260930T121832Z` (`run_r2_moments.sbatch`) reads only
  `design_individuals.tsv`, the crosswalk and the ancestry axes. It reads no rules and no e other than z_b.

So the Gate A design export and the R² moments with the fixed P2 resamples need no rerun.

**C5 and A3 need a rerun.** They used the v1.6 rules and the multiplied e.
- *Expected counts.* The C5 expected counts that A3 consumed are
  `EX/codex-rec-20260929T142434Z/expected_counts_v16_21986201/expected_gene_rows.tsv` (sha256 `043ded25…`). Its
  `summary.json` records `"spec": "v1.6"` and the rules {GSE135251: (10, 3, 15), GSE162694: (20, 5, 20),
  GSE213621: (10, 1, 10)}. These are the T2B-0930 rules.
- *Information weights.* `information_weights_v16_21986233/weights_for_full_z.tsv.gz` (sha256 `027aa5d7…`) lists
  `T2B-0930/cohort_rules.tsv` and the r1c `cohort_model_inputs.tsv` among its input hashes.
- *Which e.* AX `2026-09-30_codex_gateA.md` line 5 says these use "v1.6 rules and calibrated background", and line
  261 says "calibrated e bounds".
- *A3.* The A3 run (job 22261696, `P/trackA/A3_info/state.json`) hashed both files as inputs. Its export
  (I, V_c, R_c) therefore sits at the v1.6 rules and the multiplied e.
- *What must be rerun* before the Gate A decision, at the T2B-0929 rules and e_model:
  - the C5 expected counts;
  - the information weights;
  - the optimistic-information bound (`optimistic_information_v16_*`, which gave 0.028–0.033);
  - the A3 owner run.

  The D-QC-1 ruling already lists C5 (A3) and Gate A among the consumers to rerun.
- *Only the rerun counts at Gate A.* The Gate A decision, the κ-eligible gene list and its hash use only the MDE
  and eligibility from this rerun. Any Gate A figure computed on the v1.6 export is superseded and is not used.
- *Timing of D-QC-1 relative to the v1.6 export* (file times, America/New_York, 2026-10-15):
  - the A3 v1.6 export `P/trackA/A3_info/export/c5_information_summary.json` was written at 18:54:51;
  - the D-QC-1 recommendation (`P/trackA/A0_reviews/qc/ROUND2_ANSWERS.md`) was written at 18:56:40 from QC job
    22261693, which ran in parallel with A3. Its text says the Gate A inputs were "not checked here"
    (`P/decisions_pending.md` lines 30–31);
  - the lead ruling (`P/trackA/LEAD_RULINGS_W1_2026-10-15.md`) was written at 19:34:39;
  - the first Gate A MDE files under `P/trackA/A4_gateA/` were written from 19:41 on, after the ruling.

  So no Gate A MDE file existed when D-QC-1 was ruled. Whether anyone opened the A3 v1.6 export between 18:54 and
  19:34 is not recorded in any file [Unverified]. The ruling rests on the QC evidence at the RBMX site
  chrX:136874059, which uses no information quantity. I is a design quantity that uses no allelic ratio. The
  order of magnitude of the v1.6 information was already public on 2026-09-30 (optimistic κ detection bound
  0.028–0.033, `P/trackA/CONTEXT_modelA.md`), far from both Gate A thresholds.
- *Which v1.8 numbers this replaces.* The following v1.8 figures come from the v1.6 C5 run (CG lines 22–27 equal
  `expected_counts_v16_21986201/summary.json`; for example, planning τ_max at θ = 0 is 0.03365238):
  - the planning τ_max values 0.033652, 0.034934, 0.048558 and 0.064373 (v1.8 lines 221–222);
  - the κ-eligible counts 379, 364, 255 and 166, and the 582 genes with positive expectation (§4, line 192);
  - the expected included rows 36,644.75 and 36,366.19 (line 197);
  - "284 of the 379 planning-eligible genes" (cells 25 and 26, line 473);
  - the 15.5× and 3.7× scale statements (lines 101–102).

  Each is replaced by the value from the rerun when it is written. The rules that define them do not change:
  τ_max = Δ × median |α| over the κ-eligible genes (§5.2), and the §4 G bound.
  [Inference] The T2B-0929 rules admit more true hets (development expected true het calls rise from 54,406 to
  57,811, D-QC-1), so the eligible counts may rise. If G at the chosen θ exceeds 400, the §4 planning bound is
  restated. The ε_tie rule already scales with G_fit.

**Amplitudes in v1.8 that used the multiplied e, restated at e_model.**
- *Cell 23 (T2, v1.8 line 471).* e_gen,i = 2·e_i if S_i = 2, else e_i, with e_i = e_model. The fit keeps e_model.
  The values are e_model_p50 × 2:

  | Cohort | e_model (p50) | ×2 |
  |---|---|---|
  | GSE130970 | 0.0010601016 | 0.0021202033 |
  | GSE135251 | 0.0009635777 | 0.0019271554 |
  | GSE162694 | 0.0056017931 | 0.0112035861 |
  | GSE213621 | 0.0008767901 | 0.0017535802 |
  | GSE240729 | 0.0030200136 (p50) | 0.0060400272 (p50), applied per library |

  Source: QC `r1c-20260930T123747Z/cohort_model_inputs.tsv`, column `e_model_p50`.
  The v1.8 values (0.004363, 0.012256, 0.034532, 0.004412) were e_model × multiplier × 2 and are withdrawn. The
  GSE240729 value is unchanged, because its multiplier was 1.0.
- *Cell 34 (T2 at the threshold amplitude, line 477).* The operator computes the S = 2 multiplier on e_model. The
  definition is unchanged.
- *§5.1 (line 219).* "e_i in force" means e_model as above.

No other v1.8 amplitude uses e. Cells 22, 30 and 33 use U95; cell 24 uses ρ; cells 25–29 and 31–32 use LD, tag
switching or δ_n; and N4 uses ω_S, ρ slopes and δ.

**Measured inputs that depend on het calls (v1.8 "Open inputs" item 2, lines 577–580).** D-QC-1 re-binds R_c, so
the B-QC measurements that use RNA het calls at null sites are made at the T2B-0929 rules and e_model:
- U95 (cells 22, 30, 33, 40 and 41);
- the null-site het ρ per cohort (base ρ_ig of cell 24).

The Geuvadis per-gene ρ SD s_r (cells 40, 41) and the two-component BB tail fit (w_hi, m_hi, cell 24) use Geuvadis
DNA-het rows, so no RNA het rule enters them. They are unchanged by D-QC-1.

The simulation operator's predicted technical-cell biases and the cell 34 and 35 threshold amplitudes (v1.8
lines 461–464) use the A3 information summaries. They are computed from the A3 rerun above, not from the v1.6
export.

### 2. D-QC-2 to D-QC-5

**D-QC-2: excess-het tag screen, a prespecified A9 sensitivity.**
- *Before Gate B.* No tag-level het call is read.
- *Who runs it and when.* In the development run (plan A9), after the primary P1, P2 and P3 results and κ̂ are
  written, a named restricted-tree operator runs the screen inside `rna_genotypes_restricted/`.
- *Test (the A9 rule, defined here).* It keeps the one-sided binomial form, the minimum of 10 called libraries
  and the Bonferroni level of `scripts/analysis/in_biopsy_allelic/qc/qc_technical_v2.py` lines 165–177 (sha256
  `ee3490cc…0e9a`). Its null differs from that code, which tests every site against one rate pooled over all
  eligible sites (lines 169–170). Here each tag has its own model-based null:
  - for each kept tag t, a one-sided binomial test of its het-call count over F libraries with depth ≥ min_dp
    ("called") against p0_t;
  - p0_t = (expected het calls at t) / (called at t), with the expectation computed outcome-free from the Gate A
    design (genotype prior, depths, e_model, rules in force);
  - only tags with ≥ 10 called libraries are tested (`X_FLAG_MIN_CALLED` = 10);
  - Bonferroni over the tested tags at α = 0.05 (`X_FLAG_ALPHA`).
- *Approximation.* Libraries differ in depth, so the exact null count is a Poisson-binomial sum with mean
  called·p0_t. The binomial with the same mean has at least as large a variance. For counts at least one above the
  mean, which covers every Bonferroni rejection in practice, its upper tail is at least as large, so the test is
  conservative [Inference; Hoeffding's 1956 binomial bound for sums of unequal Bernoulli trials, recalled, not
  re-checked].
- *Sensitivity fit.* Flagged tags are dropped, the chosen tag is re-chosen by the §1.2 rule, and P1, P3 and κ̂ are
  refitted. The result is reported beside the primary. The primary does not change.
- *Stage.* The screen pools all F libraries. No stage-stratified view of the screen or of its flags is computed
  before the primary development result is written.
- *Export.* Aggregates only: the number tested, the number flagged, and flagged tag IDs with ≥ 10 het calls.

**D-QC-3: PRJNA512027 own-e_i sensitivity.**
- *As ruled.* Fit PRJNA512027 with its own per-library e_model instead of the cohort fill. The split-half Spearman
  is 0.355 and Spearman–Brown 0.524 (`P/decisions_pending.md` D-QC-3).
- *Problem with the placement.* No Gate B dataset contains PRJNA512027. Gate B redraws the F design, and F is five
  cohorts (`design_individuals.tsv`: GSE130970 54, GSE135251 160, GSE162694 83, GSE213621 208, GSE240729 53).
  Spec §1.1 lines 39–40 put PRJNA512027 in M only, because its staged individuals are all S = 2.
- *Placement bound here.* The sensitivity is bound where PRJNA512027 enters: the §8.1 α_g Catalog fit on M.
  - It runs after the primary §8.1 fit is written.
  - The refit replaces only the PRJNA512027 e_i, with its own per-library e_model; all else is as the primary.
  - Per-gene α̂ changes are reported beside the primary.
  - It is reported only. It changes no primary α̂, no Gate C input and no Catalog field.
- The lead is asked to record this placement as a correction of the ruling text "Gate B sensitivity"
  (`P/decisions_pending.md`, "Model A addendum v1.8a").

**D-QC-4: ω_c = 0 in GSE213621, a Gate B sensitivity.**
- Primary ω_c stays 0.0800427 (r1c `omega_c`, REF 0.520).
- The sensitivity is cell 52 (below).

**D-QC-5: wording.**
- *Replacement sentence.* The v1.4 sentence "GSE213621 stays FFPE-treated on its null-site G1 fail" (spec line 562)
  is read as: "Transition tags are excluded in GSE213621 and GSE240729 by prespecification (transition-tag
  exclusion). G1 is undecidable in the primary set, with a C/G-allele loss of 0.027–0.043 within the range of
  non-FFPE cohorts. GSE213621 tissue is fresh-frozen."
- *Source for 0.027–0.043.* QC `r1c-20260930T123747Z/cohort_model_inputs.tsv`, column
  `g1_ffpe_unit_development_stat`, development rows: GSE213621 −0.0272777, GSE240729 −0.0430882. The non-FFPE
  development rows in the same column run from −0.0582 (PRJNA512027) to +0.0425 (GSE135251). The per-site input is
  QC `v2-20260929T143503Z/g1_ffpe_null_sites.tsv` (`P/trackA/A0_reviews/qc/ROUND2_ANSWERS.md` line 44).
- *Source for fresh-frozen.* PMC9957752 (PMID 36535507), Methods, "Human Liver Samples": "We utilized frozen liver
  biopsies prospectively collected …" (`P/trackA/A0_ffpe/GSE213621_preservation_quote.md`).
- *Tag table.* Unchanged (sha256 `b061a7c2…`).
- *Tables and paper.* Tables that reach the manuscript call the flag "transition-tag exclusion (prespecified)", not
  FFPE.

**Sensitivity cells added (reported, not judged).** All three are κ = 0, τ = 0, Λ0 and R = 1,000, with seeds as in
v1.8 §5.1 and these cell ids.
- Cells 50 and 52 follow the §5.5 limitation rule for technical cells, because their fit is the primary fit.
- Cell 51 uses the v1.6 configuration in both generator and fit, so its size describes the v1.6 analysis, not the
  primary one. It is reported without the limitation rule. A size or bias outside the §5.5 limits there is stated
  as a property of the v1.6 configuration only.
- *Multiplier precision.* The cell 50 generator uses `e_multiplier_applied` at full precision from QC
  `r1c-20260930T123747Z/cohort_model_inputs.tsv`, development rows: GSE130970 2.0578876050101367, GSE135251
  6.359740702208149, GSE162694 3.082208914359038, GSE213621 2.5158819388863902, GSE240729 1.0. The 4-decimal
  values in the §1 table are for reading only. Cell 51 uses the same values.

| id | Sensitivity | Generator | Fit | Tests |
|---|---|---|---|---|
| 50 | D-QC-1, e | e_gen,i = v1.6 e_model_calibrated (e_model × `e_multiplier_applied`, full precision as above); rules T2B-0929 | e_model; rules T2B-0929 | P1 1; P3 1 (size, κ bias). Extra false hets can show up as gene heterogeneity, so P3 size is read. |
| 51 | D-QC-1, full v1.6 configuration | e_model_calibrated and T2B-0930 rules | e_model_calibrated and T2B-0930 rules | P1 1 (size, κ bias) |
| 52 | D-QC-4 | ω_c = 0 (REF 0.500) in GSE213621; other cohorts as primary | ω_c primary (0.0800427 in GSE213621) | P1 1 (size, κ bias) |

### 3. D-SPEC-1: P2 weights, bound rule

- *Primary.* V_c is the primary P2 weight (v1.8 §3). R_c is reported.
- *Trigger.* The trigger is evaluated per λ source on the outcome-free A3 export from the rerun of §1. It is never
  evaluated on the development-run recomputation. It fires if:
  - any signed R_c is negative; or
  - |λ_V − λ_R| > 0.05, where λ_V = Σ_c V_c λ_cond,c / Σ_c V_c and λ_R = Σ_c R_c λ_cond,c, at the point
    λ_cond,c of that source.
- *Effect.* If the trigger fires, that source's P2 p-value is max(p_Vc, p_Rc), each computed by the §5.2
  construction. P2 still rejects only if both sources reject.
- This is bound for the record. It cannot change the decision, because P2 cannot pass (§4).

### 4. A2b outcome and its consequence for P2

**Outcome.** From `P/trackA/A2b_gse193066/run-full-20261015T225326Z/summary.json` (sha256 `6596940e…6640`):
- `lambda_arm_decision` = UNDECIDABLE.
- Primary ρ = 0.08, which is the top of the grid. Of the 58 recorded pairs, 34 are SAME, 5 DIFFERENT and
  19 UNDECIDED.
- Sensitivity ρ = 0.01 gives 26 SAME, 7 DIFFERENT and 25 UNDECIDED, which is NOT_CARRIED.
- Reason recorded: "SAME 34 < 40 <= SAME + UNDECIDED 53; sensitivity rho gives NOT_CARRIED (7 recorded pairs
  called DIFFERENT (> 5))".
- Synthetic controls: C1 sensitivity 1.00 (99/99), C2 0 false SAME of 1,000, C3 0 false DIFFERENT of 99.

**Consequence.** These follow from spec Appendix B2, the user decision of 2026-10-15 and D-LABEL-2:
- P2 cannot pass. p₂ = 1 enters Holm, and the family stays {P1, P2, P3} with m = 3 (spec §5.4, lines 311–312;
  B5 default).
- P2 is reported as "not testable: GSE193066 lambda arm undecidable". This label takes precedence over "not
  calibrated" (D-SPEC-2), because the D-SPEC-2 cells are not run (§5).
- No threshold, Δ or family changes.
- Holm with p₂ = 1 compares the smaller of p₁ and p₃ with 0.05/3 and the larger with 0.05/2. That is the existing
  procedure.
- **Bound: no persistence claim from any λ source.** No κ/λ interval or ratio, no τ upper bound, and no
  "persistent", "stage-invariant" or "unchanged across stage" statement about cis effects reaches a table, figure
  or the manuscript, from GSE193066, PXD051911 A1 or any other λ source. These are the statements P2 exists to
  test, and P2 is not testable. κ̂ with its P1 interval (the slope on observed labels) and λ (label reliability)
  may each be reported, separately and without dividing one by the other. A P1 non-rejection is reported as "no
  stage-associated change detected at this power", with the MDE.
- **Bound: the A9 P2 resample refits are not run.** With p₂ = 1 fixed and no κ/λ interval reported, the 2,000
  P2 resample refits have no consumer.

### 5. Gate B cells not run

**The criterion.** A cell is not run if its only purpose is P2 calibration or P2 power. A cell that feeds a P1,
P3, κ-coverage or τ̂-bias judgment stays, even if it also reports P2.

**The judgments in v1.8 §5.4 (lines 364–397) and §6, and the cells each reads:**
- (a) P1 warp size: cells 1, 2, 3, 40.
- (b) κ coverage: cells 1, 2, 3, 40, all Λ0 (§6, line 449).
- (c) P3 warp size: cell 1.
- (d) Full CRT: cell 99.
- (e) τ̂ bias, before and after n*: cells 2, 3 and the six τ > 0 cells among 4–12.
- (f) §9 δ_n trigger: cell 27.
- (g) P2 calibration rule: cells 13, 14, 43.
- (h) §5.5 limitation rule: technical cells 22–30 and 33, and judged null cells.

| v1.8 cell | Purpose in v1.8 | Read by (a)–(f) or (h)? | Status |
|---|---|---|---|
| 13 B1 (κ_true = +Δ under Λ1, so the attenuated slope is +Δλ: the P2 boundary plant κ = +Δλ) | P2 κ-part size; full P2 on r < 1,000 | No. Read only by (g). It appears in the §9 rerun list (line 378) as a rerun target, not as a trigger input. | not run |
| 14 B2 (κ_true = −Δ, boundary plant κ = −Δλ) | as 13 | No, as 13 | not run |
| 43 B3′ (τ = τ0 = τ_max·λ^(0.05): the P2 boundary plant τ = λ·τ_max) and its τ-fallback rerun | P2 τ-part size; fallback trigger | No. Read only by (g). §9 rerun list as for 13. | not run, including the fallback (v1.8 line 438) |
| 15 B3 (τ_max under Λ1) | P2 τ part, reported | No. It is not a power cell, so (e) does not read it. The 5.1a projection applies to it only at generation. | not run |
| 16 B4 (Λ3 tail) | full P2, the named P2 limitation (§5.5) | No | not run |
| 17, 18 P2 power | full P2 | No | not run |

**Cells 19–21 (L1–L3).** These are not P2-only, because each also reports P1 power and coverage at κ = −0.10
under label error. None is read by (a)–(h), since all are "not judged" (line 335).
- **Cell 20 (L2, Λ2 = PXD051911 A1 draws) stays.** It runs its P1 part (K = 1). Its P2 part (K = 2, P2 interval
  coverage, §6 line 454) is not computed, because no P2 quantity enters any report once P2 is "not testable".
- **Cells 19 (L1, Λ1) and 21 (L3, Λ3) are not run.** Their generators are defined on "the carried GSE193066 arm"
  (§8, lines 513 and 515). No arm is carried, so the generator is undefined. They are not redefined on A1, because
  that would add a new scenario.

The P1-under-label-error report is carried by cell 20. [Inference] Label error cannot change P1 size at κ = 0,
because with κ = 0 and τ = 0 the reads do not depend on X̃ (v1.8 §8 generator: η = α_g(1 + κ_true·X̃ + δᵀz) +
u_g·X̃). So no P1 validity claim depends on 19 or 21. No L cell had τ > 0 or a P3 test, so no P3 robustness plant
is lost.

**Bound: P1 power under label error is scoped to A1 Y5.** Dropping cells 19 and 21 removes the P1 power reports at
the GSE193066 planning λ (0.5754) and the Λ3 tail (about 0.19; v1.8 lines 513–515). Any P1 power under label error
that is reported names its source, PXD051911 A1 Y5 (planning λ 0.7355), and is not extended to lower λ. The lead
is asked to confirm the drop (`P/decisions_pending.md`, "Model A addendum v1.8a").

**Unaffected.** All other cells run as in v1.8: 1–12, 22–35, 40–42 and 99. This includes:
- the λ_gen,cond report, now for cell 20 only;
- the §9 δ_n rule, whose rerun list becomes cells 1, 2, 3, 40 and 99.

**τ-bound coverage is not computed.** Spec §4 (line 227) and v1.5 (line 567) ask Gate B to report "τ-bound
coverage" beside τ̂ bias. The only τ bound the spec defines is P2's τ component (§5.2, lines 285–288: the
profile root r at τ0). P2 is not testable, and section 4 of this addendum bars any τ upper bound from the
manuscript, so this report has no object. v1.8 never assigned or costed it, so the cost figures below do not
change. τ̂ bias, before and after the n* restriction, is still computed, and the spec §4 fallback rule still
applies to it, because it bears on the section-4 fit that P3 uses.

The P2 seed streams 3 and 6 are unused.

**Revised counts.** [Estimate, not run]
- *Datasets.* v1.8 has 58,000 warp datasets. Removed: 16,000 (13 and 14: 8,000; 15: 1,000; 43: 2,000; 16: 1,000;
  17 and 18: 2,000; 19 and 21: 2,000). Added: 3,000 (cells 50–52). The result is **45,000 warp datasets**, plus the
  200 full-CRT replicates.
- *Section-4 datasets.* v1.8 has 39,000. Removed: 11,000 (13 and 14 at r < 1,000: 2,000; 15: 1,000; 43: 2,000;
  16: 1,000; 17 and 18: 2,000; 19–21 P2 τ parts: 3,000). Added: 1,000 (cell 50). The result is **29,000**.

**Revised CPU-h.** [Estimate, not run; v1.8 §5.7 method: section-3 fit 5–60 s, section-4 fit 2–25 min]

*Fit accounting per dataset.* Section-3 fits = 1 + 2·K_P1 + 2·K_P3 + 2·K_P2κ. Section-4 fits = 2 (section-4 MLE,
if any section-4 work) + 2·K_P3 + 2 (P2 τ0 profile, if a P2 τ part).

*Check of the accounting.* This accounting reproduces v1.8's judged tier exactly: 3,443–42,587 CPU-h against the
printed 3,400–42,600. For the full v1.8 warp list it gives 6,014–74,900, against the printed 6,500–80,700, so v1.8's
reported-cell figure carries about 490–5,800 CPU-h that I could not trace to a cell. Both bases are given below.

| Part | Explicit count | v1.8 total minus removed cells, plus new cells |
|---|---|---|
| Removed cells: 13–18, 43, 19, 21, and cell 20's P2 part (108,000 section-3 and 44,000 section-4 fits) | −1,617 to −20,133 | same |
| Added cells 50–52 (11,000 section-3 and 4,000 section-4 fits) | +148 to +1,850 | same |
| Warp cells after the change (297,000 section-3 and 124,000 section-4 fits) | 4,546–56,617 | not computed on this basis |
| Full CRT (unchanged) | 556–6,670 | 560–6,670 (v1.8 printed, included in its total) |
| **Total** | **5,100–63,300** | **5,600–69,100** |
| Judged tier: cells 1, 2, 3, 40, 27 and the full CRT | 2,807–34,687 | not computed on this basis |

The upper end on the conservative basis (69,100) is above the 64,000 pause (`P/budget.json`, 0.8 × 80,000). The
v1.8 §5.7 rule applies unchanged: the A6 measured times replace these figures, and if the projection exceeds 64,000,
only the judged tier is submitted before the user is asked.

### 6. Gate A branches now that P2 cannot pass

The rule is unchanged. Its three texts:
- plan A4 (line 125): "MDE ≤0.125 → full design; ≤0.25 → P2 and τ untestable; otherwise per-gene fields only and
  Model D becomes the lead";
- `P/trackA/CONTEXT_modelA.md` line 66: "MDE ≤ 0.25 → P2 and τ are untestable, and P1, P3 and per-gene results
  remain";
- the source request, AX `2026-09-29_claude_to_codex_gateA_inputs.md` lines 42–53: item 2 defines the τ quantity
  as P2's "τ-component bound at τ = 0", and the verdict reads "MDE ≤ 2Δ: P2 and τ untestable".

So "τ untestable" in the plan means P2's τ-component bound, not P3. P3 stays in the middle branch, as CONTEXT
line 66 says. With P2 already removed the rule reads:

| κ MDE (at θ chosen by Gate A) | Result |
|---|---|
| ≤ 0.125 | P1 + P3 + per-gene fields. The "P2 bound fitting" condition no longer separates any outcome, because P2 is not testable whatever its half-width. |
| > 0.125 and ≤ 0.25 | P1 + P3 + per-gene fields (P2 was already untestable in this branch) |
| > 0.25 | per-gene fields only, and Model D becomes the lead (Track B item B8) |

This is the existing rule read with P2 removed. The thresholds 0.125 and 0.25 are not changed and no threshold is
added. The MDE is computed on the C5 and A3 rerun at the T2B-0929 rules and e_model (§1), not on the v1.6
export.

### 7. λ for the retained label-draw plant and the judgments

- **Binding (D-LABEL-3, option (a)).**
  - *Which λ.* The B-LABEL label-draw plants keep A1 Y5 (PXD051911). The planning point is λ = 0.7355
    (0.7354848564; 300-draw 95% interval 0.5639–0.9725; AX `2026-09-29_codex_lambda_check.md` line 63).
  - *Generator.* The one retained label-draw plant is cell 20 (Λ2). Its generator uses the A1 Y5 draws from the
    B-LABEL full run (B = 2,000, at the v1.3 fixed R²), run on the D-LABEL-1 labels-only crosswalk copy.
  - *Caveat.* The error at the S cutpoints is extrapolated from the F0/F1 boundary. The level-dependent fit and
    the b profile are reported as a sensitivity.
- **κ coverage and τ̂ bias read no λ.** In v1.8 both judgments are made only in Λ0 cells:
  - κ coverage in cells 1, 2, 3 and 40 (§6, lines 449–453);
  - τ̂ bias in cells 2, 3 and the τ > 0 power cells, against 0.1·τ_max, where τ_max = Δ × median |α| has no λ.

  So "keep A1 Y5" holds for these two judgments only in the sense that no other λ enters them. The one judged
  quantity that used λ (B3′, τ0 = τ_max·λ^(0.05)) is not run (§5).
- **What still needs λ.**
  - the B-LABEL full run, which feeds cell 20 and is unchanged by the lead ruling;
  - V_c, for cell 20's λ_gen,cond report only.

  V_c no longer blocks any judged cell. v1.8's "Open inputs" item 4 (line 582) blocks Gate A only through I (the
  κ information), which comes from the A3 rerun.

### 8. What this addendum leaves as it was

- §1–§4 of v1.8;
- the Gate B criteria (§5.4 items 1–4);
- the false-fail figure for Gate B as a whole (≈ 0.105);
- D-SPEC-2 and D-SPEC-3 as bound in v1.8;
- the §11 agreement targets and the B-FIX v1.8 rerun family;
- the real-data seeds;
- the seal prerequisites (§10 item 10).

## Revision v1.8.2 and v1.8.3 (2026-10-16): Newton polish acceptance, the observed section-3 fit, and observed-fit failure

Written by the Model A spec owner (model program `model-program-20261015T223347Z`, stage `trackA/A6_gateB`).

**Status.** A pre-Gate-B completion, adopted by the user on 2026-10-16. The user adopted the polish rule as
v1.8.2 (`P/trackA/LEAD_RULINGS_W4_2026-10-16.md`, H-10). H-10 left open how the reported fit is chosen among
attempts. The lead closed that the same day in H-11 and H-12, and those rules are v1.8.3. H-9 and H-13 supply the
implementation notes in item 7. This is not the single pre-outcome amendment, which stays unused. It binds from the
moment it is appended. The spec's sha256 before and after the append is recorded in `P/trackA/A6_gateB/state.json`
under `spec_v182_v183`.

**No Gate B outcome was read.** The inputs were convergence flags, attempt counts, CPU time, bitwise equality of end
points, and ℓ differences between end points of the same dataset. They came from Gate B-type null datasets drawn
with measurement seeds (base [20261016, 995], not the Gate B seeds) and from the B-FIX fixture. No κ̂, SE, z, LR,
v̂, rejection indicator or coverage indicator of any Gate B dataset was read or printed.

**Self-contained.** This section is written so that an implementation can be built from the spec text alone,
without B-MODEL code. Where it says "B-MODEL uses", that describes the reference implementation's choice in a place
the spec leaves to each implementation (v1.7 rule 8). That description does not bind other implementations.

Abbreviations follow v1.8: AX = `docs/technical/agent_exchange`; EX = `Analysis/MASLD_Model_Benchmark/executions`;
P = `EX/model-program-20261015T223347Z`.

**What does not change.**
- Everything that v1.8 and v1.8a list as unchanged. That includes the likelihood (§2.3), the mean and dispersion
  (§3), the section-4 integral and its integrators (v1.8 §1), the three test statistics and their p-value formulas,
  Holm over {P1, P2, P3}, Δ = 0.125, the draw counts, the real-data seeds of §5, the §11 agreement tolerances and
  the B3 bands.
- The convergence criterion (v1.8 §2 item 1) and the final-point integrator check (v1.8 §1 item 8).
- The §7 start, the optimizer (L-BFGS-B followed by a projected Newton polish) and the limit of 3 restarts per fit.
- Every refit (P1 CRT refits, P2 κ resample refits, P2 τ-profile fits, P3 bootstrap fits). Each still has a warm fit
  and a cold fit, kept by the v1.8 §2 item 3 table with the ε_tie given there. Only the polish rule of item 1 applies
  to them.
- Section-4 fitting, including the observed section-4 fit (two v starts, kept by the v1.8 §2 item 3 table), apart
  from the polish rule of item 1.
- The §5.1 sandwich SE formula and v1.7 rule 3.
- The consequences of failure. A failed refit is an exceedance. In Gate B a failed observed fit counts as a rejection
  (v1.8 §5.2) and as not covering (v1.8 §6). On real data a failed observed fit makes P1 "not computed" (§5.1).
  Item 4 below defines what "failed observed fit" means. It does not change these consequences.

**Terms.**
- *Attempt.* One L-BFGS-B run from one start point, then one Newton polish. In section 4 this also includes the
  v1.8 §1 item 8 final-point check and any re-polish it triggers. The *end point* of an attempt is the point where
  its last polish stops.
- *Converged attempt.* Its end point meets v1.8 §2 item 1: the projected max |gradient| is below 1e-6 over the free
  coordinates, and the Cholesky factorization of the negative projected Hessian succeeds. In addition, none of the
  following occurred: the v1.8 §1 item 6 failure (8,192 panels without agreement), the v1.8 §1 item 8 failure (a
  gene's integrator changes again at the new final point), or an evaluation error at the end point.
- *Fit.* The §7 sequence of at most 4 attempts from one start point (item 2).
- *Start.* The point at which a fit begins.
- ℓ is the log-likelihood in nats. In section-4 comparisons it is evaluated with the reference integrator for every
  gene, as v1.8 §2 item 3 says.

### 1. Step acceptance in every Newton polish (v1.8.2; replaces v1.8 §2 item 2)

The rule applies to every Newton polish. That means section-3 fits and section-4 fits, both observed fits and
refits, and the v1.8 §1 item 8 re-polish.

1. *ℓ_best.* ℓ_best is the highest ℓ among the accepted iterates of the current polish. The polish start counts as
   accepted, so ℓ_best starts at ℓ there. A trial that is rejected does not change ℓ_best. A re-polish under v1.8
   §1 item 8 is a new polish: ℓ_best starts again at ℓ of the re-polish start, evaluated with the integrators held.
2. *ε_ℓ.* At the start of each Newton iteration, set ε_ℓ = max(1e-9, 1e-13·|ℓ_best|). This value is used for every
   trial in that iteration.
3. *Step.* Let pg be the projected max |gradient| at the current iterate, using the v1.8 §2 item 1 projection. The
   step d is computed and applied on the free coordinates. A coordinate that sits on a bound with its gradient
   pointing out of the box is left out. Every trial point is clipped to the §3.4 bounds.
4. *Trials.* θ_t = clip(θ + t·d) for t = 1, 1/2, 1/4, …, for as long as t > 1e-12. A trial whose evaluation fails,
   or whose ℓ is not finite, is rejected. A trial is accepted if either of these holds:
   - (a) ℓ_trial ≥ ℓ_best − ε_ℓ, and the projected max |gradient| at θ_t is strictly below pg;
   - (b) ℓ_trial > ℓ_best + ε_ℓ.

   The first accepted trial becomes the current iterate, and ℓ_best is set to max(ℓ_best, ℓ_trial).
5. *Stopping.* The polish stops at the current iterate in any of three cases:
   - no trial with t > 1e-12 is accepted;
   - pg < 1e-6 and the Newton decrement gᵀd (gradient g on the step's coordinates) is below 1e-12;
   - the implementation's iteration cap is reached.

   Convergence is then judged at the stop point under v1.8 §2 item 1.
6. *Left to each implementation* (v1.7 rule 8): the direction d, as long as it is an ascent direction built from the
   projected Newton system; the iteration cap; and the L-BFGS-B settings.
   - B-MODEL's direction: it solves (−H)d = g on the step's coordinates after symmetric Jacobi scaling, with each
     eigenvalue replaced by its absolute value and floored at 1e-10 times the largest.
   - B-MODEL's iteration cap is 60.
   - B-MODEL's L-BFGS-B ftol is 1e-15 in sections 3 and 4 (H-9 rejected 1e-10 for section 3).

**What changed from v1.8.** Under v1.8, ℓ_best was "the highest ℓ evaluated so far", so a rejected trial could raise
it. Every accepted step also had to lower the projected gradient. A polish could therefore stop with no step accepted
next to a point of higher ℓ, and the fit then failed the convergence test.
- Clause (a) is the v1.8 test, measured against accepted iterates only. An accepted iterate is never more than ε_ℓ
  below ℓ_best, and ℓ_best never falls, so losses cannot add up.
- Clause (b) accepts a step that raises ℓ by more than ε_ℓ even when the projected gradient rises.

### 2. Attempts within one fit (§7; stated in full, unchanged)

1. Attempt 0 starts at the fit's start. If it does not converge, up to 3 restarts follow (attempts 1, 2, 3). The
   sequence stops at the first converged attempt.
2. If an attempt converges, the fit is converged, and its end point is that attempt's end point.
3. If no attempt converges, the fit is not converged. It carries the end point of its highest-ℓ attempt.
4. Each restart starts from a jitter of the highest-ℓ attempt so far. The jitter is left to each implementation
   (v1.7 rule 8). B-MODEL uses the form in item 3.2: s_j computed at the fit's own start, and entropy
   [20260923, 90, a] for restart a. That entropy is the same in every dataset.
5. For a refit, this is the whole rule inside each of its warm and cold fits. An earlier unconverged attempt of
   higher ℓ does not change a converged warm or cold fit. The observed section-3 fit adds the check in item 3.4,
   step 3.

### 3. The observed section-3 fit from two starts (v1.8.3; new v1.8 §2 item 3a)

**Scope.** Item 3 applies to every observed section-3 fit:
- each Gate B warp dataset, pilot datasets included;
- each full-CRT replicate (cell 99);
- every B-FIX workflow;
- every real-data observed section-3 fit, including A9 and the sealed run.

Refits do not use it. Their warm start is the observed MLE that item 3 reports.

**3.1 Starts.** n_obs = 2.
- Start 0 is the §7 start: the §7 per-gene starting values, with v = 0 and the no-information coordinates held as
  in v1.7 rule 3.
- Start 1 is a jitter of the §7 start (3.2).
- Each start is a full fit with its own §7 restarts (item 2).

**3.2 Start 1.**
1. x0 is the §7 start, clipped to the §3.4 bounds.
2. s_j = 1/√|H_jj|, where H is the section-3 Hessian of ℓ at x0. This applies to every coordinate j that is free in
   section 3 and has a finite H_jj with |H_jj| > 1e-12. For every other coordinate, s_j = 1. If the Hessian cannot
   be evaluated at x0, s_j = 1 for all j.
3. e is a vector with one entry for every coordinate of θ, in the §7 coordinate order, free or not. It is drawn as
   numpy `default_rng(SeedSequence(entropy)).normal(0.0, 0.1, size)`, with the entropy given in 3.3.
4. Start 1 equals x0 + s_j·e_j on the coordinates that are free in section 3, clipped to the bounds. Every other
   coordinate keeps its x0 value.

**3.3 Entropy of start 1.** These values are fixed before any outcome is read. Stream 7 is added to the v1.8 §5.1
table.

| Run | Entropy |
|---|---|
| Gate B warp dataset r of cell c | [20261015, c, r, 7, 1] |
| Full-CRT replicate r | [20261015, 99, r, 7, 1]. The replicate's 1,000 inner CRT draws are refits and have no start 1. |
| Real data (every real-data observed section-3 fit) | [20260923, 7, 0]. No test uses test_id 7 (P1 = 1, P2 = 2, P3 = 3). |
| B-FIX | [20260923, 7, 0] |

New row of the v1.8 §5.1 stream table:

| Stream | What it draws |
|---|---|
| 7 | start 1 of the observed section-3 fit (k = 1) |

The first v1.8.3 build used a provisional stream label. These entries replace it (H-12 m2, m3).

**3.4 Selection.** Let C be the starts that converged with a finite ℓ, and let ℓ_0 and ℓ_1 be the ℓ of starts 0
and 1.
1. *No converged start.* If C is empty, the observed fit fails.
2. *Candidate.* If both starts converged, start 1 is the candidate when ℓ_1 > ℓ_0 + 1e-9, and start 0 otherwise.
   So a tie within 1e-9 goes to the §7 start. If only one start converged, that start is the candidate.
3. *Unconverged higher.* Let U be the end points of every attempt in either start that did not converge. U includes
   the attempts before the converged attempt inside a converged start, and every attempt of a start that did not
   converge. If any end point in U has ℓ > ℓ_candidate + 1e-6, the observed fit fails. The margin 1e-6 is the
   section-3 ε_tie of v1.8 §2 item 3. This makes it the same rule the warm/cold table applies to refits: an
   unconverged search found a higher ℓ, so no certified maximum exists. Non-finite ℓ in U is ignored.
4. *Reported fit.* Otherwise the candidate is the reported observed fit:
   - its θ is the observed MLE;
   - κ̂ and the §5.1 sandwich SE at that point give z and the κ interval;
   - refits warm-start from that θ.

The two margins do different jobs:
- 1e-9 decides which of two certified maxima is reported. In the measurements, end points of one maximum reached from
  different starts or hosts differed by at most 1.2e-10.
- 1e-6 decides failure, in step with the refit table.

**3.5 Record per observed fit.** For each start: converged or not, attempts, ℓ minus ℓ_0, and whether its θ is
bitwise equal to start 0's θ. Also the reported start, the selection rule that applied, and the CPU time over both
starts. These fields are recorded beside the fit. They are not outcomes.

### 4. Failed observed fit (definition)

Wherever this spec says that an observed fit failed, it means one of the following:
1. the observed section-3 fit fails under item 3.4 (no converged start, or an unconverged end point more than 1e-6
   above the candidate);
2. the §5.1 sandwich SE at the reported point is not finite and positive, so z is undefined;
3. for statistics that use the observed section-4 fit (P3, and any other statistic built on that fit), that fit is a
   failed fit under the v1.8 §2 item 3 table. This affects only those statistics.

Case 2 is the rule the A6 pilot already used to count observed-fit failures (`P/trackA/A6_gateB/pilot-20261016T025010Z/GO_RULE.md`,
"Failure counting"). Case 3 restates v1.8 §5.2 ("failures are handled as for P1"). The consequences are the ones
listed under "What does not change".

### 5. Refits (unchanged)

- Every refit has a warm fit, started from the reported observed θ, and a cold fit, started from the §7 start.
  Each is a fit under item 2. The v1.8 §2 item 3 table keeps one of them or declares the refit failed, with
  ε_tie = 1e-6 in section 3 and max(1e-6, 1e-8·G_fit) in section 4.
- Refits get no start 1. The check in item 3.4, step 3 does not apply inside a refit's warm or cold fit.
- The polish rule of item 1 applies to refits.

The observed fit is stricter than a refit in one respect: an unconverged attempt inside a converged start can fail
it. That can only add observed-fit failures. In Gate B these count as rejections and as not covering, so the
asymmetry cannot lower a measured size or raise a measured coverage [Inference].

### 6. Why two starts and not one

- Under one start, nothing flags an observed fit that converges at a lower maximum. Each refit already has two
  starts and the "unconverged higher means failed" rule. With one observed start, the null would not regenerate the
  observed statistic with the same search.
- Two observed starts match the two refit starts.
- n_obs = 2 was fixed on cost grounds. It was not tuned to the datasets still running (H-11 M1, v182-REVIEW m2).

### 7. Implementation notes (not likelihood changes)

These notes change no formula, threshold or selection rule. They fix how the judged runs are executed and recorded.

- *Threads.* Every Gate B, B-FIX, A9 and sealed-run job runs with one BLAS/OpenMP thread:
  OMP_NUM_THREADS = MKL_NUM_THREADS = OPENBLAS_NUM_THREADS = 1 (H-9).
- *Environment.* All Gate B, A9 and sealed runs use the conda environment `vcc26`, for CPU and GPU work alike. Each
  job writes its pip freeze (H-13). The reason: the same CPU code in two environments moved one P1 draw's z by
  2.5e-5.
- *Headers.* Every job records host, CPU model, glibc version, node type and the three thread variables. A
  per-dataset record carries these, or points to the job header that does.
- *Node type.* Gate B CPU jobs run on AMD EPYC nodes (SLURM feature v6) for cost. The node type does not define the
  reported fit; item 3 does.
  - The v6 hosts fall into two glibc groups, 2.35 and 2.39. Same-start end points were bitwise equal within a group
    (11 of 11) and in no pair across groups (0 of 19) (v182-REVIEW M2).
  - At a dataset with two maxima, one start can reach a different maximum on each group. Item 3 reports the higher
    maximum whenever either start reaches it.
  - Bitwise reproduction is expected only within one glibc group. Between implementations, agreement is judged by
    the §11 tolerances, not by bitwise equality.
- *Optimizer.* L-BFGS-B followed by the Newton polish, as §7 says. Projected Newton and trust-constr were measured
  and not adopted (H-13).
- *Evaluation failures and run errors.*
  - An evaluation failure comes from the model itself. Either a v1.8 §1 rule cannot complete for a gene (for example
    item 6 at 8,192 panels), or ℓ, its gradient or its Hessian is undefined or not finite at the point.
  - An evaluation failure rejects a polish trial (item 1). At an end point it makes the attempt not converged. It can
    lead to a failed fit only through the rules above.
  - Errors of the computing system are run errors: device out-of-memory, a driver fault, an unsupported
    configuration. A run error stops the task and is reported. It is never entered as ℓ = −∞, as a rejected trial or
    as a failed fit.
  - A GPU backend must classify outcomes the same way as the CPU code (H-13 M1).
- *GPU before judged GPU fits.* Two checks are required (H-13 m3, m4), judged against the pass criteria that were
  hashed before any GPU result (`P/trackA/A6_gateB/gpu-backend-20261016T170127Z/PASS_CRITERIA.md`):
  - a cross-node repeat of one full observed fit;
  - an order-2 comparison on a problem with both reference-integrator and GH20 genes.

### 8. Measured basis (no Gate B outcome read)

**Polish rule.** Observed section-3 fits on 80 measurement datasets per node type (cells 1, 2, 3 and 40).
- Under the v1.8 rule, 5 of 80 failed the convergence test on v6 and 3 of 80 on v4 (H-10).
- Run from the same L-BFGS-B end points, the item 1 rule converged in attempt 1 in 80 of 80 on each node type.
- Where v1.8 had converged, the item 1 rule gave bitwise-identical end points in 62 of 62 fits (v6) and 65 of 66
  (v4). The one exception moved ℓ by 5.8e-11.
- B-FIX P1 workflow: 123 of 123 fits bitwise identical between v1.8.1 and v1.8.2, with P1 p = 0.17647 under both.
- B-FIX P3 workflow: 323 of 323 fits bitwise identical, with P3 p = 1.0 under both.
- Source: `P/trackA/A6_gateB/v182-REVIEW.md`.

**Section-4 polish** (H-11 M3; run directory `P/trackA/A6_gateB/v183-20261016T163618Z`, check (b)).
- Ten measurement datasets, each fitted under v1.8.1 (section-3 ftol 1e-15) and v1.8.3 in the same task: 8 tasks on
  glibc 2.35 and 2 on 2.39.
- Section-4 fits converged in 10 of 10 under both versions, with equal attempt counts in 10 of 10.
- End points were bitwise equal in 9 of 10. The tenth dataset was the one whose two-start observed fit reported
  start 1, so its section-4 start differed. Started from the same section-3 point, its section-4 end points were
  bitwise equal.
- Clause (b) did not change any section-4 end point in these data.

**Two maxima** (v182-REVIEW, 49 measurement datasets with 4 starts each).
- The §7 start missed the highest converged point once: cell 2 r16, by 4.83 in ℓ, with both points passing the
  convergence test. The miss rate for one start is 1/49 (Clopper–Pearson 95% upper bound 0.109).
- The best of the §7 start and the first jittered start missed 0 of 49 (upper bound 0.073).

**v1.8.3 build.**
- B-FIX P1 under two starts: equal to v1.8.2 in every draw, P1 p = 0.17647. Start 1 converged 9.8e-12 above
  start 0, which is a tie, so start 0 was reported (v183-REVIEW).
- Tests: 34 passed (check (c), job 22272008).
- These checks used the provisional stream label of 3.3 and did not include the item 3.4 step 3 check.

**Cost** [Estimate, not run; H-11].
- n_obs = 2 adds about 812 CPU-h in the judged tier, which is 0.7% of the 118,046 CPU-h projection, at 170 s per
  extra start.
- At the 49-dataset median of 188 s per jittered start, the figure is about 1.1 times larger (v182-REVIEW m3).
- The v1.8 §5.7 pause rule applies unchanged.

### 9. Open before the first judged Gate B dataset

1. *Code.* The code must implement the 3.3 entropies and the item 3.4 step 3 check (H-12: v1.8.3a, in both the CPU
   and GPU copies). B-FIX and any Gate B pilot runs that set the reported fit are then rerun under that code.
2. *GPU checks.* H-13 m3 and m4 (item 7) are required before any judged GPU fit.
3. *For the lead: roots of a failed observed fit.* v1.8 §5.2 pools R·K roots, and a failed CRT refit sets
   |z*_rk| = +∞. It does not say what a dataset whose observed fit failed contributes to the pooled quantile, since
   its refits have no certified warm start. The A6 pilot counted those roots as failed refits. This revision leaves
   the point as v1.8 has it.
