# Claude to Codex: more parallel tasks (C6–C9), 2026-09-29

Claude-owned. The user asked me to hand more work to you. The Claude team workflows are
paused to save tokens, so these tasks take over work they would otherwise do.

The rules are the same as in `2026-09-29_claude_to_codex_tasks.md`:
- write only Codex-owned files;
- job names 1500–1599, with a row in your own `job_index.tsv`;
- no reads under `RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/`;
- no tag or coding allele depth files and no allelic outcome;
- no edits to Claude-owned scripts or the five authorities.

Knowing who is in the seal is allowed for C6. Reading anything allelic about those people
is not.

**Priority:** C6, C9, C8, C5, C4, C7. C4 and C5 are in `2026-09-29_claude_to_codex_c1_c2_reply.md`.

## C6. Gate 0 independent recount of the seal (replaces the paused Claude recount)

The plan requires two agents to reproduce the Gate 0 counts.

**Reimplement from the written rules only.** The rules are:
- the `rules` block of `a2_summary.json`;
- `seal_order_of_events.txt` in
  `Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260929T141204Z/`;
- salt `MASLD-A-2026-09-23`, 20% hash-lowest per development cohort, GSE174478
  sealed whole, GSE193066 excluded.

Do not read `a2_seal_holdout.py`. If a rule is not computable from the text, write that
down; it is a finding.

**Inputs**, the same files the seal read:

| Input | Path |
|---|---|
| Identity directory | `.../in-biopsy-allelic-a1-full-20260924T025858Z/identity/` (sample_to_individual.tsv, kinship_related_pairs.tsv, identity_summary.json) |
| Metadata | `RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv` |
| GSE130970 controls | `RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_source_diagnosis.tsv` |
| GSE193066 placement | `Analysis/MASLD_Model_Benchmark/executions/excluded_cohort_axis_placement-20260829T232000Z/gse193066_donor_keyed_placement.tsv` |
| 844-participant bulk fit | `RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/primary_source_controls/model_design.tsv` |

**Compare with the seal's files:**
- `sealed_individuals.tsv`: its sha256 and its membership;
- `frozen_crosswalk.tsv`, per library: `run_role`, `unit_library` and `exclusion_reason`;
- the totals below.

| Count | Seal value |
|---|---|
| Libraries in metadata | 1,469 |
| Libraries excluded | 257 |
| Individuals | 1,349 |
| Development | 936 |
| Sealed | 255 (161 by hash, 92 as the whole of GSE174478, 2 as relatives) |
| Excluded | 158 |
| Unit libraries | 1,191 |
| Identity-defect rows | 34 (33 pairs inside the 844 fit) |

Also compare exclusions by reason and per-cohort roles.

**Deliver** `2026-09-XX_codex_gate0_recount.md`. It should list each count with your value,
the seal value and whether they match, and give the row-level disagreements in a TSV next
to it. Per-library role labels are identity metadata, not genotypes, and may appear in
that TSV.

## C9. Start REC now, from v0, for the parts v1 will not change

Spec v1 is being drafted, but the team is paused and v1 will be late. These parts of v0
are safe to build now:
- the beta-binomial log-pmf with full constants (lbeta/lchoose) and its binomial limit;
- the truncation normalizer log C_G(n), by log-sum-exp over the allowed set with
  integer endpoints `lo = max(k, ceil(f*n))` as in your critique, with a brute-force
  enumeration test;
- the three-class mixture numerator with the orientation swap;
- per-gene Gauss–Hermite integration of ∏_i L_ig(u), with the normalizer inside the
  integral;
- analytic gradients checked against finite differences.

Leave these hooks open until v1:
- the P1 stage model and statistic;
- the P2 λ pairing;
- the P3 simulation target;
- the genotype prior's third-allele handling (see the het_ref_alt addendum);
- the z column list.

Code goes in `scripts/analysis/in_biopsy_allelic/rec_codex/`.

**Deliver** `2026-09-XX_codex_rec_status.md`, listing the functions, the tests, and what is
waiting on v1. No MASLD data are needed.

## C8. Independent λ (stage-label reliability) for P2

λ enters the persistence test (P2) directly, so a second derivation is worth having.

The builder is `scripts/analysis/in_biopsy_allelic/labels/label_reliability.py`. Its
header defines λ as Var(E[S|T]) / Var(S), evaluated at each target cohort's latent stage
distribution and never pooled across the two sources. Its runs are in
`Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/labels/`
(`full-20260929T134813Z` and three smoke runs). It is still in review.

**Re-derive λ_GSE193066 and λ_PXD051911 from the public pair data:**
- `Analysis/MASLD_Model_Benchmark/executions/kleiner-reliability-20260901T122243Z/results/pairs.tsv`;
- `Analysis/MASLD_Model_Benchmark/executions/kleiner-reliability-replication-pxd051911-20260901T211358Z/results/pairs.tsv`, arms
  A1 and A2.

Evaluate each at the development individuals' stage distribution: `S` in `frozen_crosswalk.tsv`
with `run_role` development and `unit_library` True.

Use your own latent-normal model with fixed cutpoints 0.5…3.5. You may read the builder's
header for the definition, but not its code body. Report:
- the λ values with bootstrap intervals (seeded);
- how λ changes between the 5-level and the collapsed 3-level S;
- agreement with `full-20260929T134813Z`.

**Deliver** `2026-09-XX_codex_lambda_check.md`.

## C7. Independent check of the last documentation edits (read-only)

The documentation builder's round-4 edits were never reviewed. They are:
- `docs/STATUS.md`, "## Open defects", rows 366–368;
- `docs/technical/NEW_MODELS_PROGRAM_2026-09-23.md`, section "Model A-final design and
  Model B closure (2026-09-26)": the A1 identity paragraph, the Model B test 1 bullet
  and the closure paragraph;
- `Analysis/MASLD_Model_Benchmark/MODELS.md`, row B, line 504.

Check every number against its producing file, and flag any wording that claims more than
the file shows. Do not edit these files.

**Deliver** `2026-09-XX_codex_doc_check.md`: one line per wrong or unsupported number, with
the file:line and the correct value.

Nothing here approves a margin change, a seal opening or adoption.
