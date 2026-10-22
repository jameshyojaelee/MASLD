# Remediation record, 2026-09-07

Response to `../REVIEW.md`. Section numbers are the plan's
(`~/.claude/plans/noble-percolating-curry.md`). Every number written into prose
here was recomputed from its producing file in the same session; none was
carried across from the review.

## Retired numbers removed, with what replaced them

| number removed | was at | why it had to go | replaced by |
|---|---|---|---|
| adjusted OR **2.1 against 9.1** (genetic-only vs disease-state-only at the protein layer) | `draft/fig5_discussion.md:102` | no producing artifact on any release; traces only to the July ledger | deleted. The class comparison it attempted is respecified as a failable test (plan §3) |
| **corrected gsMap** class analysis paragraph | `draft/fig5_discussion.md:104` | the July ledger records the gsMap class gate as FAIL and lists the numbers under "Retired from the active manuscript" | deleted |
| **strict drug benchmark** vs modality count paragraph | `draft/fig5_discussion.md:106` | same; the development-gradient test is recorded as invalid for want of a curated preclinical stage | deleted |
| credible-set linkage covers **75 of 473** primary genes | `draft/fig5_discussion.md:127` | 473 is the superseded SuSiE count (adopted: 462), and the sentence disclaims a coding-architecture conclusion that the adopted Figure 2D result supports at BH q = 6.0e-4 | deleted |
| "Across **Figures 1–6**" | `draft/fig5_discussion.md:70`, `05_figure_legends.md:8` | seventh figure added 2026-08-30; nothing downstream was updated | "Figures 1–7" |
| "**six-figure** story" | `STATUS.md:54,76`, `README.md:70`, `05_figure_legends.md:27` | same | "seven-figure" |
| enzyme instances collapse to **259** operational clusters | `HIGH_IMPACT_PATTERNS.md` row 11 | no producing artifact | **266** (coding test) and **263** (the five noncoding-context tests), from `n_clusters_enzyme` in `GWAS/finemapping/results/fig2_multiplicity/20260817T161128Z/fig2_family_multiplicity.tsv` |
| non-overlap **92.5%** under the TREAT comparator | `RESULTS.md` | rounding of 92.4485% | **92.4%**, and the "producing run is unrecorded" warning removed — see below |

## Numbers added because they were missing

| number added | where | source |
|---|---|---|
| **1,066** canonical DEGs inside the 14,920-gene jointly testable background | `RESULTS.md` adopted table | `overlap_by_arm.tsv`, column `n_deg_joint`. This count entered the overlap test but appeared in no authority; only the superseded 1,079 did |
| **92.9%** non-overlap expected by chance, beside the observed 95.4% | `RESULTS.md`, `draft/01_intro.md`, `draft/fig5_discussion.md`, Figure 1B candidate | `(437 − 31.2227)/437`, recomputed. This is the review's finding 1: two sets of these sizes are 92.9% disjoint by arithmetic, so the observed 95.4% is 2.6 points of signal, not the organizing observation |

## Defects the review reported that turned out not to exist

| reported defect | finding |
|---|---|
| GNMT protein **−1.136 / q 0.00642** over 6,071 proteins has no producing table | **It exists.** `Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv` has exactly 6,071 data rows and gives GNMT logFC −1.13581364725247, padj 0.00641955641526544. Producer `Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R:145-151`, design `~ group + acquisition_batch + age_z + bmi_z + sex`, which is what `RESULTS.md` describes. The 6,089-protein / −1.128 value is the separate BMI+age-only arm in `Analysis/Proteomics/`. **No number changed.** |
| the TREAT-comparator overlap row has no producing run | **It exists.** `RNA-seq/results/genetics_state_interface/20260817T162108Z_promoted_coloc_raw_treat/matched_treat_by_arm.tsv` gives `n_case_treat = 33` at `treat_lfc = 0.25` on the `susie_pp4_0.5`/`official_27638` arm, and its `input_manifest.json` carries the same three input sha256 values as the canonical run. The warning is removed and the row now cites the file. |

## Figure 1B

The rendered main panel was read with `pdftotext` and prints **34** shared,
**413** colocalized-only, **92.4%** non-overlap and **1,227** disease-state-only
under a DEG arm labelled "Canonical DEGs (interval-null FDR gate)". Those are
July values under the retired TREAT gate. Cause:
`scripts/figures/fig1_fig4_evidence_classes.R:16-17` defaults
`MANUSCRIPT_RELEASE_ID` to `2026-07-15-r2`, and no `manuscript_release/`
directory exists for the adopted release, so re-running the script cannot fix it.

Rebuilt as a candidate at
`figures/candidates/fig1b-adopted-release-rebuild-20260907T105813Z/`, reading
the adopted interface output directly: 20 shared, 417 colocalized-only, 95.4%,
1,046 disease-state-only, with a dashed line marking the 7.1% and 2.9% shared
fractions that independence predicts. Values read back off the rendered PDF and
inspected as PNG. `figures/main/` is unchanged.

**Figure 1C has the same defect and is not yet fixable.** It prints 8 + 14 + 451
= **473**, the superseded count. Rebuilding it needs the direct-disease versus
liver-enzyme stratum-scope classification, which `HIGH_IMPACT_PATTERNS.md` row
12 records as not derivable from the adopted registry alone. Neither panel
should be adopted until 1C is regenerated.

## Figure 7 surfaces created

`draft/fig7.md` and a Figure 7 section in `05_figure_legends.md` did not exist;
no draft prose mentioned the seventh figure at all. Both were written from
verified producing files:

| statement | verified against |
|---|---|
| 585 external participants; ρ 0.638–0.674 against ceilings 0.848–0.942; 0.68–0.79 of ceiling | `executions/external-validation-expansion-20260901T190000Z/results.json`; 109 + 299 + 177 = 585 recomputed |
| paired deltas −0.040 and −0.011, reference cohort +0.054; the registered primary is the weakest arm | same file, `vs_native_paired_delta` |
| a nonlinear model with no biological pretraining differs from ridge by −0.046 to +0.001 over six comparisons against MDE 0.067 | `executions/tabicl-nonlinear-comparator-20260905T000000Z/results/results.json`. **The figure "0.000 ± 0.067" quoted in the review and strategy exists nowhere in the repository**; these are the real deltas, and every one carries `exceeds_mde: false` |
| within-participant Kleiner σ_e 0.4655–0.6433 stages | `RESULTS.md:769-800`; two deposits, reported as a range because it does not pool |
| two pretrained single-cell models report no delta because zero held-out studies are clean | `release/encoder-benchmark-v1/MODEL_CARD.md:32-33` |

## Still open in this section

- **7F is an unweighted mean over studies.** `scripts/score_encoder_field_common_head_50000.py:384,648-653` averages study-level scores arithmetically, so GSE189600 (2 donors, 149 cells) carries the same weight as a 46-donor study and contributes −0.455 where every other study contributes −0.065 to −0.002. The legend now says deltas are reported unweighted and donor-weighted; the donor-weighted arm has not been computed yet.
- **The TabICL result has no write-up** in `MODEL_RESULTS_SUMMARY.md`, `MODELS.md`, or `docs/`. The Figure 7 legend and draft are currently its only prose home.
- `release/masld-severity-v1` is superseded by `v1.1` (bitwise identical scores, plus a measured input contract). Citations should move to v1.1.
- Hallmark "median 148 genes" in `BULK-PROGRAM-MAP-v9/HANDOFF.md:78` is not reproducible from `vocabulary_coverage.tsv` (median observed 163, median expected 180).
