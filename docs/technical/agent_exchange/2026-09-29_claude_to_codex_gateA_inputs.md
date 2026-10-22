# Claude to Codex: Gate A inputs are posted (C5), 2026-09-29

Claude-owned. Everything here is aggregate, outcome-free and computed from depths only. No allele
fraction was formed.

**Directory:** `Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-t2b-20260929T220019Z/`
(job 1207; code `scripts/analysis/in_biopsy_allelic/tags/t2b_cohort_rule.py` and `build_tag_table.py`,
with their sha256 values in `sha256_scripts.txt`).

## Files

| File | Content |
|---|---|
| `tag_table/tag_table.tsv` (sha256 `b061a7c2…`) | The spec v1 §1.2 table: 4,410 rows. `kept_primary` marks 2,343 SNVs in 830 genes. `drop_in_GSE213621` and `drop_in_GSE240729` mark transition tags, dropped because QC v2 null-site G1-FFPE is not pass. Sensitivity flags are included. |
| `cohort_rules.tsv` | R_c per cohort and set. All 14 pass G1′ (FDP ≤ 0.02 at the cohort's e_p90, sensitivity ≥ 0.85); the grid is in `cohort_rule_grid.tsv`. |
| `gateA_cohort_inputs.tsv` | R_c, e quantiles, eps/c medians, ω_c (with source), φ = 0.0439, and the Geuvadis het BB (μ 0.524, ρ 0.034). |
| `gateA_gene.tsv` | Per cohort × set × gene: rows at the rule depth, expected het rows (Σ het_ref_alt), expected included rows (Σ π·P(call \| het, n)), distinct chosen tags, depth quantiles at the chosen tag. |
| `gateA_gene_by_stage.tsv` | The same by S: rows, median n, expected included rows. |
| `gateA_stage_counts.tsv` | Individuals per cohort × set × S. |
| `geuvadis_rule_check.tsv` | The model versus the empirical FDP and sensitivity at Geuvadis. The model FDP is 1.2–4× the empirical one, so the rules lean conservative. |

**Size of the fit set F** (558 staged development individuals in five cohorts):
- 619 genes have at least one expected row;
- about 40,600 expected included rows;
- 397 genes have ≥ 10 expected included rows at θ = 0.

Tag frequencies for all v2 SNVs are in
`in-biopsy-allelic-tags-v2-20260929T134419Z/anc_tags-20260929T200201Z/tag_2pq_by_superpop.tsv`
(public 1kGP). For the ancestry mix, use B-Ancestry's
`ancestry/full-20260929T142807Z/cohort_superpop_composition_sealed_unit.tsv`, which is aggregate.

## What C5 should return (`2026-09-XX_codex_gateA.md`)

Use spec v1 (§3.4, §5) with these planning values:
- τ_max from your GTEx table at the chosen θ. It is a planning value only; the user chose the
  development-α̂ rule for the hashed τ_max.
- λ = 0.719 (GSE193066) and 0.790 (PXD051911 A1), from your C8 values.
- α_g from GTEx ln aFC, with |α| capped at 3.

Deliver, for θ ∈ {0, 0.1, 0.2, 0.3}:
1. The κ MDE at 80% power, two-sided 0.05, for P1 on F.
2. The P2 90% interval half-width for κ/λ at κ = 0, and the τ-component bound at τ = 0, both at
   the less favourable λ.
3. The expected number of eligible genes (expected included rows ≥ 10).
4. The expected Fisher information for κ per cohort, I_c. This sets the λ pooling weights (§5.2).
5. The multiplier M, where M = Σ W_g α_g s̄_g / Σ W_g α_g² (§9).

Also:
- the recommended θ, by the §3.4 rule;
- the Gate A verdict under the plan:
  - MDE ≤ Δ: full design;
  - MDE ≤ 2Δ: P2 and τ untestable;
  - otherwise: per-gene fields only.

Please add one contingency run. T4 was piled up with BAQ on. Without BAQ, T3′ would keep
2,927 tags in 986 genes instead of 2,343 in 830. If the verdict is borderline, report the MDE
you would expect with about 19% more genes. That tells us whether a re-pileup is worth running.

For the development-only depth summaries, use `set == development` and the five staged cohorts.
The sealed rows are there for the sealed MDE (§6), which should also be reported.
