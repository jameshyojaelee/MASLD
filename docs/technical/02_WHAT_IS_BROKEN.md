# What is broken, ordered by how much it costs

Generated 2026-10-07 from an audit of every tracked file in this repository, the
2026-08-13 manuscript draft, the June 2026 Drive framework document, and the lab
notebook. Items that cite a repository `path:line` were re-checked against the code
on 2026-10-07 (34 citations checked, corrections applied in this revision). Items
without a repository citation are marked **[Uncited]**: their numbers come from
audit-agent summaries, the manuscript or the Drive document, and have not been
re-checked against a producing file. `manuscript.txt:N` refers to a line-wrapped
plain-text export of the 2026-08-13 draft; those line numbers will not match the
Google Doc and should be located by quoted text.

**Revision 2026-10-07 (post-review).** Corrections: 1.1 (cohort count; the
sample-disjoint analysis already exists), 1.2 (palindrome handling is fixed in
current code; the guard is a regression test, severity reworded), 1.3 (the Tier-1
gate change was a missing-data fix), 2.3 (the full concordance line is now quoted),
Tier 3 (five plans, not four), and the closing paragraph (which plan items
re-derive the substrate).

This document exists because the hypothesis inventory
(`01_HYPOTHESIS_INVENTORY.md`) is not the project's problem. There are about 330
reconstructed hypotheses and no shortage of ideas. The problem is that a
substantial number of the **load-bearing** ones rest on substrate that cannot
currently support them, and the repository already knows it. Most of these were
found by this project's own audits and then not propagated.

Severity is judged on one axis: **if a reviewer found this independently, what
would it cost?**

---

## Tier 1 — a reviewer finding this independently would be fatal

### 1.1 A validation set that is inside the discovery set

The convergence score's headline calibration is `AUROC 0.939` against the
"Govaere 2020 NASH signature, a gene set not used in its construction"
(`manuscript.txt:54`). The Govaere 2020 signature is derived from **GSE135251**.
GSE135251 carries `include_in_mega: true` (`config/human_datasets.yaml:59-71`)
and is the **largest** of the pooled cohorts that produce the bulk differential
expression feeding the score's S1 channel. (`config/human_datasets.yaml` marks six
cohorts `include_in_mega: true`; `lomo_validation.R:345` describes GSE135251 as one
of ten cohorts in the dream mega-analysis. The pooled-cohort count is itself
inconsistent across the repository.)
(`RNA-seq/46d_convergence_evidence.R:357-385`).

The repository's own audit file already calls this leakage
(`RNA-seq/46d_baselines_and_lomo.R:18-19,277-280`). The leave-one-modality-out
analysis written to replace the headline is itself buggy: the S2 column name is
wrong, so `LOMO_S2` equals the full score, and the score being tested is not the
canonical one (`46d_baselines_and_lomo.R:121-128`). The "nearest gene" baseline
is the top-PP4 score (`:111`).

A sample-disjoint version already exists: `RNA-seq/scripts/wave1_d_lomo/lomo_validation.R:344-358`
re-scores Govaere with S1 rebuilt from a leave-GSE135251-out dream fit
(`dream_loo_GSE135251.csv`, overall contrast only, other four S1 sub-contrasts
zeroed). Its number is on the HPC, not in the repository, and the draft does not
use it. The Figure 5 calibration panel still scores the in-domain Govaere set
(`scripts/figures/fig5_panel_calibration.R:231-233`).

**The fix is therefore adoption, not new analysis:** report the sample-disjoint
AUROC (after fixing the `LOMO_S2` column name) and drop "a gene set not used in its
construction".

**Cost if found externally:** the central calibration claim, and with it the
argument that convergence is informative, is withdrawn under review rather than
on the authors' terms.

### 1.2 The figure's claim doubles as the pipeline's pass condition

`GWAS/finemapping/src/07_combine_susie_coloc.R:513-535` is a control-recovery
regression guard (C8). It stops the build if TM6SF2, GCKR, THRB or EXOC3L4 fall
below PP.H4 0.5, **or** if PNPLA3, HSD17B13 or MBOAT7 rise above 0.5. It does not
alter any value. The problem is that the "controls" are the same genes Figure 2
reports as its coding-versus-expression partition. A run in which the partition
fails never produces a combined table, so every shipped table agrees with the
figure by construction. This is recoverable: disclose the guard and report the
anchor values from a run with the guard disabled. Left undisclosed, a reviewer
reads it as an unfalsifiable figure.

PNPLA3's causal variant rs738409 is a C/G palindrome. Palindromes **used to be**
dropped unconditionally. Current code keeps them by default after a strand
certificate (`06_susie_coloc.R:429-470`; the old drop is opt-in via
`COLOC_DROP_PALINDROMIC=1`), a change committed 2026-08-10 (`f88ed33`), three days
before the 2026-08-13 draft. **[Unverified]** whether the draft's
`PNPLA3 PP.H4 = 0.012` (`manuscript.txt:23`) came from a run before or after that
change. Item 11 of the plan re-computes it both ways.

A separate contradiction sits in the same arm: the manuscript states SuSiE-coloc
is the primary test with ABF as fallback (`manuscript.txt:19`), while the code
states ABF is the "PRIMARY/headline arm" and SuSiE-coloc covers 22,072 of 927,225
rows, 2.4% (`06_susie_coloc.R:256-259`).

**Cost if found externally:** high if undisclosed, because it reads as the figure
being enforced. Low if disclosed with unguarded values reported alongside.

### 1.3 Clinical outcome is an input to a score claimed to use none

The manuscript states the ranking "separates the approved target THRB from the
failed target NR1H4 using no clinical information" (`manuscript.txt:11,60`). The
tier assignment includes a C4 trial-positive term (`RNA-seq/40_drug_target_validation.R:244`),
obeticholic acid is coded inconsistently (`:78`), and those tiers then serve as
machine-learning labels and validation sets
(`45b_ml_target_prediction.py:173-178`; `46a_cumulative_enrichment.R:41-46`).

The C4 term is counted directly in the tier tally:
`n_pass <- sum(c(c1_deg, c2_coloc, c3_tf_clean, c4_trial_pos))` (`:252`).

Separately, the Tier-1 gate was changed in a way that returned THRB and RORA to
Tier 1 (`46d_convergence_evidence.R:1266-1269`). The code comment gives the reason:
their INTACT score was missing (NA), which sank them to `4_Weak`, so the gate was
moved to COLOC PP.H4. That is a data-coverage fix, not tuning to a target, but it
was made with the expected answer in view, and the score was debugged against
expected genes (`:1342`). It belongs in a disclosed list of analyst decisions.

**Cost if found externally:** the C4 term alone makes "using no clinical
information" false as written. It is the most damaging single finding in this
audit. The execution plan does not re-derive a clinical-free tiering; item 20
withdraws the tier labels as inputs. Withdrawing the claim is the cheaper option.
Re-deriving the tiers without C4 and re-running the downstream ML would keep it.

### 1.4 Reported donor counts are library counts

The manuscript reports 269 QC-passing single-cell donors, 211 in stage-resolved
analyses with steatohepatitis n=101 (`manuscript.txt:38`). The true counts are
**102 analyzed donors and 73 in the stage fit**
(`Analysis/SingleCell/scripts/lib_donor_collapse.R:5-22`;
`511_donorlevel_disease_stage.R:9-16`). Every single-cell p-value computed before
2026-07-12 is pseudoreplicated.

The consequence is already measured and is severe. The manuscript's claim of 18
of 30 stage-associated hepatocyte modules with 12 replicating in bulk
(`manuscript.txt:39`) becomes **2 robust and 2 HC3-fragile of 117** in the frozen
donor-level registry (`516:27-38`).

### 1.5 The two surviving programs may be contamination **[Uncited]**

Both programs that survive donor-level testing — hepatocyte module 8
"Stromal ECM/IGFBP7" and module 20 "Ductular injury/BICC1" — are textbook
signatures of ambient RNA, doublets or mislabeled lineage inside hepatocyte
barcodes, and Figure 4D shows them scoring **higher in non-hepatocyte lineages**.
Their hepatocyte-intrinsic interpretation is unresolved and under decontX
sensitivity analysis.

### 1.6 The stage axis may be a dataset axis

99.94% of steatohepatitis cells come from GSE244832, and the
steatosis-to-steatohepatitis change "CANNOT be tested cross-cohort"
(`Analysis/SingleCell/scripts/361:5-8,238-245`). Dataset was also the scVI batch
key, so batch correction and the stage signal are entangled.

This directly undermines the manuscript's single-cell claim that reorganization
is sharpest at the steatosis-to-steatohepatitis transition
(`manuscript.txt:41,62`).

### 1.7 Control labels were wrong in the largest adopted analysis

GSE130970 controls were defined as `fibrosis_stage == 0`, which placed **17
source NAFLD biopsies (16 steatotic) in the control arm**
(`scripts/manuscript/resource_f_five/refit_corrected_controls.R:5-8`). Every
adopted bulk headline — 1,347 canonical DEGs, the Figure 3 triad cohort node,
Figure 3G, the evidence classes, the Gene Catalog, the continuum, and all
F0-referenced stage contrasts — was computed on those labels. A corrected refit
to 140 control / 704 disease is wired in the 2026-10-01 commit, and **no
corrected result is adopted**.

---

## Tier 2 — the claim survives but at a much smaller scope than stated

### 2.1 The project's own benchmark contradicts the progression narrative

`scripts/analysis/continuous_axis_benchmark/11_assemble_report.R` concluded that
no data-driven ordering adds information beyond recorded histology, that a random
expression direction explains more held-out genes than the inferred axis
(`:127`), and that remodeling across stages is dominated by rank-1 amplification
with change-point entropy below the amplification null (`:170-174,253-257`). The
benchmark is declared non-Resource (`:10-12`), while the histology-anchored
continuum was promoted into the working Figure 4 supplements.

The scope validator now bans the words "progression" and "cascade"
(`scripts/manuscript/validate_resource_scope.py:276`), which is the correct
response, but the draft still advances the cascade framing throughout.

The one longitudinal test available failed: median cosine 0.164, 95% CI
-0.033 to 0.244, P=0.218 on the locked GSE193066 paired-biopsy holdout **[Uncited]**.
Note that 2.5 below finds 11 of those 58 pairs are different people, so this
test ran on a partly invalid substrate in either direction.

### 2.2 "Multi-modal convergence" counts one modality three times

`layers_active` sums L1 (differential expression), L3 (cross-species
concordance, which requires human DE), L4 (genetics), L5, L6 (GSEA leading edges
from the same DE ranking) and L7 (`RNA-seq/27a_assemble_evidence_atlas.R:1531-1596`).
L1, L3 and L6 all derive from the same bulk DE. **L7 records only whether DepMap
has the gene**, adding about +1 to most protein-coding genes. L4 ORs about 20
correlated colocalization sources and takes an uncorrected minimum TWAS p<0.05
across many GWAS (`27a:228-231,1539-1579`), so roughly 5% of 20,000 genes acquire
a false "genetic" call.

The independence argument — max pairwise |rho| = 0.18, therefore agreement
reflects independent biology (`manuscript.txt:52,58`) — does not follow. Low rank
correlation does not establish independent failure modes, and every
expression-based column shares expression-level and literature-attention
confounders. The druggability column directly encodes prior study attention, and
DGIdb was queried **only for upregulated DEGs** (`20_drug_repurposing_v3.R:186-200,397-400`),
so the conserved-core druggability enrichment (OR 4.92) and the network-proximity
result are partly true by construction.

### 2.3 Fixed-effect pooling across strongly heterogeneous cohorts

Median I² is 73.1%
(`RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/concordance_meta_vs_canonical_C2.R:9`).
The same file records the earlier meta-versus-dream comparison as rho 0.887, 100%
direction agreement and Jaccard 0.168 for the significant sets (`:8`, labelled
"OLD"). So the gene **rankings** agree closely; what differs is which genes cross
the significance cutoff. Under that heterogeneity, fixed-effect p-values are
anti-conservative: 13,043 of 27,638 genes reach padj<0.05 **[Uncited]**. Counts
of significant genes are the affected quantity, not effect directions.
Every downstream count inherits this. The DE method itself was chosen by ranking
on colocalization, mouse and Open Targets enrichment
(`build_ranked_table.R:28-35`) — the same sources later treated as independent
convergence evidence.

### 2.4 NMF programs are seed-unstable and mislabeled

The code's own audit finds k=6 seed-unstable (cross-session ARI about 0.30),
favours k=4 (cophenetic 0.875 versus 0.757), and names the Fibrogenic program the
least stable (`RNA-seq/94c_compare_refit_vs_canonical.R:222-229`). The
manuscript's "fibrogenic stromal program (PLN, MYOCD, CNN1)" is described with
smooth-muscle markers that the code's own lineage rule labels
vascular-smooth-muscle/pericyte (`95_nmf_clean_ksweep.R:98-100,129-131`). k and
the program labels are additionally chosen using their association with fibrosis
stage (`95:387-397,140-155,461-472`), after which the manuscript reports that the
programs rise with severity.

### 2.5 Cohorts that are not independent, and donors that are not distinct

- GSE106737 and GSE83452 are **one cohort**, confirmed only by expression
  fingerprinting (≥78 reciprocal r≥0.9999), collapsing a pre-registered
  three-cohort gate to two after the fact
  (`Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/17_validate_plan43.py:99-103`).
- GSE193066: only 29 of 58 recorded "paired" biopsies genotype as one person; 11
  are **different people**; 18 are uncheckable
  (`scripts/analysis/in_biopsy_allelic/labels/label_reliability.py:953-956`).
  Paired-continuum tests were run on this substrate.
- GSE213621 holds 25 genotype-confirmed same-person library pairs, 11 with
  discordant stage labels (`label_reliability.py:977-981`), inside the 844-donor
  headline substrate, which has **no donor-uniqueness gate**: the unit-of-inference
  claim is a prose comment, never a test (`refit_corrected_controls.R:27,333,709`).
- Ghodsian includes deCODE and UKBB; the Sveinbjornsson arms overlap them; PDFF
  and UKBB enzymes share the UKBB cohort
  (`GWAS/finemapping/src/00c_compare_overlapping_gwas.R:4-9`). Counting
  colocalization support across these as independent replication overstates
  convergence.

### 2.6 The genetics arm's resolution claims rest on mismatched LD

`recommended_pip` is a **max** over SuSiE, CARMA, SuSiEx and MESuSiE with a
**union** credible set (`04_aggregate_results.R:172-183,337,340`), which is
anti-conservative and contradicts the manuscript's "retained variants supported
by both" (`manuscript.txt:17`). Those PIPs were computed on indefinite LD
matrices (minimum eigenvalue -0.19 to -0.36) that `susie_rss` does not check
(`04_audit.R:38-43`); a 19% ridge-failure rate once shipped as passing, and one
run reported complete after modelling 17 of 1,950 loci (`:86-89`). About half of
163 adjudicable palindrome re-orientations were wrong — 86 correct, 77 wrong
(`seqfunc_v2/finemap/common.R:298-340`).

UKBB (n≈344k) and Biobank Japan (n=160k) summary statistics are paired by default
with 1000 Genomes panels of a few hundred samples (`10_run_susiex.py:9-12,166`).
The AFR and SAS arms of 6.6k and 8.9k cannot support the claimed resolution gains.
A single European microarray eQTL panel serves all ancestries, so
ancestry-restricted effector claims (CYP2A6) are not separable from panel
mismatch — a point a co-author already raised and which remains open.

### 2.7 Small-n claims carrying large narrative weight

Bullets without a `path:line` are **[Uncited]**.

- GSE192741 is **4 individuals, tested 2 versus 2**, not the "n=5 donors" stated.
  Scripts keyed on `sample_id` count donor H35 twice.
- The CYP3A4 "Moran's I" 0.381 → 0.465 values are **Hotspot's C statistic** from
  donor-pooled slides (`fig4f_cyp3a4_zonation.R:53-60`).
- The pericentral-versus-periportal count (8 versus 4) includes the marker genes
  that **define** the zonation axis (`spatial_params.yaml:36-51`).
- Vu's 33 biopsies from 32 patients across 10 arrays are described as "5 donors".
- The liver proteome is n=58 with 3 F3 and 0 F4 samples.
- The approved-versus-failed drug calibration is n=2 (THRB, NR1H4).
- The plasma proteome PXD052937 was **permanently excluded on 2026-08-03** for
  wrong condition labels (`recover_proteomics_metadata.R:46-57`), yet the draft
  still reports plasma rho 0.31/0.51/0.52 and "35 of 94" (`manuscript.txt:44`).

### 2.8 The four-level regulatory cascade does not hold **[Uncited]**

The draft's THRB/HNF4A/RORA/MLXIPL/MAX cascade (`manuscript.txt:56`) conflicts
with the corrected caQTL analysis: **46.3% concordance, below chance**, and **zero**
disease-regulon TFs reach four-way support. The original "12 of 24 disrupted TFs"
came from a cell-level test treating cells as independent replicates; the
donor-level re-test gives zero regulons at n=18.

---

## Tier 3 — the project's own escalation attempts all terminated

Five pre-registered plans (40–44), each with genuine sealing and hash-pinning discipline,
terminated `complete_no_promotion` or verdict FALSE:

| plan | hypothesis | outcome |
|---|---|---|
| 40 (Myojin HLF) | evidence classes predict palmitate-lipotoxicity CRISPR survival | nonconfirmatory |
| 41 | the two frozen hepatocyte programs are induced by disease cues in PCLS/organoids | failed all six gates |
| 42 (double dissociation) | genetic-only genes are functional, state-only genes are spatially organized | refuted |
| 43 (chronic-state bridge) | the 117-program axis reverses with histologic improvement | refuted |
| 44 | lineage axes encode a fat-versus-profibrotic assembly topology | refuted |

Table outcomes are **[Uncited]**: taken from audit summaries, not re-checked here.

The project's own goal-completion audit records
`n_achieved_biological_requirements = 0` of eleven
(`scripts/source_independent_risk_state_relay/49_build_goal_completion_audit.py:121-125`).
That 0 is a hard-coded literal in the script, the same pattern item 4.3 below
criticizes; it records the authors' judgement rather than computing it.

The convergence claim itself fails directly: no convergent-versus-single-map odds
ratio reaches q<0.05, and the design cannot resolve below OR 4.4-7.8 **[Uncited]**.

**This is not the same thing as the hypotheses being wrong.** Most of these
designs were powered only for very large effects (6-versus-6 and 5-versus-9 exact
permutation ATAC designs; convergent classes of n=11-34). The honest statement is
*underpowered*, not *refuted*, and the distinction is currently collapsed —
which matters because an inflated meaningful-effect anchor, derived from the
modality's own surviving positives, mechanically converts underpowered nulls into
"tested_negative" (`catalog-mde-layer.../02_unify_and_sensitivity.py:4-7`).

---

## Tier 4 — process problems that generate the above

Items 1, 2, 4 and 5 are **[Uncited]** (drawn from the drafts and audit summaries).

1. **Numbers are not frozen across drafts of the same analysis.** Lead loci
   225→241, converged 214→233, cross-ancestry 95/116→70/108, triple-ancestry
   genes 4→3, Hotspot modules 17/13/9→18/12. RORA is "concordantly downregulated
   across all five mouse diet models, three significantly" in June and "two of
   the four diet models" in August. CYP3A4's spatial autocorrelation changes
   between drafts. At most one version of each can be right.

2. **Inconsistent denominators within one draft.** The genetics arm reports 618
   (ABF), 368 (SuSiE), 751 (coding/noncoding split), 466 (method-specific split)
   and 564 (orthogonality analysis) colocalizing genes. Percentages computed on
   different bases are not comparable.

3. **Hard-coded verdicts inside scripts described as mechanical.** Claim-ledger
   statuses are authored constants independent of the computed gates
   (`build_evidence_class_release.R:491,532`); gate_pass values and expected
   q-values are literals (`risk_state_double_dissociation/07_build_verdict.py:113-157`;
   `06_validate_outputs.R:493-503`). These validators confirm that a stored file
   matches a stored expectation.

4. **Gates that cannot fail.** The orthogonality gate "≥85% genetic-not-DE" is
   nearly guaranteed by DEG base rates. The convergent-versus-genetic-only
   enrichment is circular because convergent genes are DEGs by definition and the
   endpoints are disease-DE calls.

5. **Post hoc rescue after null results.** `281` FDR=0 → `281b` effect-size
   filter → `281c` lower floor; `284` under 10 edges → `288` consensus; `283` →
   `287` PP4 0.3 with distance as an LD proxy; empty cross-species → `297` module
   cliques. No global error control spans these.

6. **Selection-conditioned worked examples.** GNMT/MAT1A/CYP2C19 were selected in
   PXD051911 by protein significance and mRNA-protein concordance, so their
   protein and concordance nodes cannot count as validation — which the code
   itself admits. The hero programs were picked from 117 after inspection. THRB,
   RORA, NR1H4, HKDC1, CYP3A4, SERPINE1, CYP2A6 and FADS2 are each presented with
   knowledge of their clinical status, with no stated denominator.

7. **Code defects that silently corrupt inputs.** The cirrhosis regex `F4`
   matches `Fibrosis_F3F4` (`05c:106-107`); the ancestry proxy is a no-op
   (`14.5c:66-70`); the drug-class regex puts statins in JAK-STAT (`247:73`); the
   F3 sub-state "dip test" is a normality test; `282`'s LOCO compares a signed
   statistic with an unsigned one; `281c` tests a 3-term statistic against a
   cached 2-term null. The BG-013 length-offset no-op is enforced with `stop()` in
   M02/M03 but **still executes and logs as applied** in M02b and M08.

8. **Almost no result is in the repository.** Every number reachable here is a
   code comment, a hard-coded assertion, or a provenance note. Two load-bearing
   manuscript numbers (the 1,355-gene conserved core, OR 4.92) are asserted in
   `_legacy` figure scripts whose figures are not in the current main-figure
   authority index. The model benchmark dashboard states plainly that "a summary
   file is a report, not a source" and documents three numbers with no producing
   file.

---

## What this adds up to

The engineering discipline in this project is genuinely unusual: outcome
firewalls, hash pinning, sealed holdouts, condition-blindness audits, independent
re-derivation, explicit untestable states, and a long list of honestly recorded
negatives. Roughly half the scripts in some subsystems are validators.

That discipline is also where the risk now sits. Validator volume has become a
proxy for rigor while the validators check provenance and arithmetic rather than
design validity, and several of the most important checks are gates that cannot
fail or verdicts that are authored constants.

The single change that would most raise the credibility of the whole Resource is
narrow: **re-derive the headline substrate once, correctly, and recompute every
downstream claim on it exactly once.** In the execution plan that is item 1
(the donor-unique substrate and the shared fail-closed validator) and item 20
(recompute the Gene Catalog once under all corrections). Items 2 and 3 run early
because they fix the target-selection rule and the single-cell program set that
later items consume, not because they re-derive the substrate. Everything else in the top 20 is a specific claim being made
honest or being withdrawn on the authors' terms rather than a reviewer's.

The three independent judges and the completeness critic all raised the same
objection to that plan, and it is recorded in `03_TOP20_EXECUTION_PLAN.md`
rather than buried: a portfolio that is eighteen-twentieths audit protects this
paper and starves the next one. The parallel track in that document is the
answer — the long-lead data applications and the label-quality investment cost
almost nothing to start now and take months to mature.
