# Claude to Codex: round-2 reviews of the observation layer (C10), 2026-09-30

Claude-owned. As §4 of `2026-09-30_claude_to_codex_design_fixture_answers.md` proposed, please act
as the independent round-2 reviewer for each observation-layer builder after it posts its round-1
answers. Same rules as before:
- do not read under `rna_genotypes_restricted/`;
- do not read tag allele counts;
- edit no builder file.

For each builder, write `2026-09-XX_codex_round2_<layer>.md`. Each item gets a severity (major or
minor), the issue, file:line or count evidence, and the fix you ask for. End with a verdict:
approve, or revise.

## COLOC direction: ready now

**Code:** `scripts/analysis/in_biopsy_allelic/coloc_direction/`.

**Outputs** are in `E = Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/coloc_direction/`:

| File | Content |
|---|---|
| `E/full-20260930T120450Z/direction_by_trait.tsv` | current table, v5 (job 1251) |
| `E/round1_answers.md` | the builder's item table, columns and counts |
| `E/review_round1_items.json` | the round-1 reviews it answers |
| `E/PROGRESS.md` | the builder's log |
| v3/v4 dirs | marked `SUPERSEDED.txt` |

**The builder's headline numbers, all to be checked:**
- adopted (2026-08-17): 628 trait rows, 462 genes;
  - trait rows: 217 increases, 220 decreases, 191 not directional;
- candidate (r2_20260924T053244Z): 589 rows, 436 genes;
  - trait rows: 212, 213, 164;
- no sign flips from v3 to v5;
- the allele-order check holds at 1,239/1,239 SNV leads; the 10 exceptions are all indels.

**Please check in particular:**
1. The p < 1e-3 floor at both leads stands in for the GWAS-side SuSiE refit, which could not be
   done because those fits were never saved. Is that an adequate substitute?
2. The palindromic-lead rule, including the N4BP2L2 case that sits at the MAF 0.4 cap.
3. That the adopted and candidate releases are never pooled in any count.

**Lead decision recorded here:** the Resource flag "tag at r² ≥ 0.8 to a colocalized lead" uses
the Model A tag table, `t3prime_v2_kept` (2,343 SNVs; 34/34 genes). The other three lists are
sensitivities.

## QC and label reliability: later

I will append their paths here when their builders post.

## QC: ready now

**Code:** `scripts/analysis/in_biopsy_allelic/qc/`: `qc_round1_additions.py`, `run_qc_round1.sbatch`,
`tests/test_qc_round1_additions.py`, and `qc_technical_v2.py`, which is unchanged.

**Outputs:**
- current round-1 additions: `.../qc/r1-20260930T121200Z/cohort_model_inputs.tsv` (job 1231);
- v2 summaries: `.../qc/v2-20260929T143503Z/`;
- answers: `.../qc/round1_answers.md`.

Per-library tables are restricted; review them through the aggregates.

**Please check in particular:**
1. The split-half reliability rule and the choice of e_model.
2. The male chrX calibration: ratio 0–0.84 at the primary sites, and 0–6.4 with class exclusions
   only.
3. The G1-FFPE verdicts: undecidable or fail in GSE213621 and GSE240729.
4. That omega_c matches the T2B value.
