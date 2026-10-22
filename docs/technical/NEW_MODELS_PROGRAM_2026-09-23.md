# Three new MASLD models and the bulk ALT-contig count check (candidate work, 2026-09-23/24)

Status: candidate. Nothing here has been applied to an adopted Resource number,
class, figure or claim. Plan: `~/.claude/plans/plan-thoroughly-to-build-imperative-wren.md`.
Jobs are named by number; each execution directory has a `job_index.tsv`
(index, job ID, purpose). Executions live under
`Analysis/MASLD_Model_Benchmark/executions/`.

## Models and prespecifications

| Model | Question | Prespecification | State |
|---|---|---|---|
| A. In-biopsy allelic | How do inherited cis-regulatory variants act in MASLD liver, in which lineage, and does the effect change with fibrosis stage? | `scripts/analysis/in_biopsy_allelic/PRESPEC_A_DRAFT.md` (rewritten 2026-09-25 for observed tag hets; hashed after the A1 crosswalk, holdout seal and T2 rule) | imputation failed its gate; tag-SNV design built; allele depths at tags collected; A1 identity done; Model A-final design chosen 2026-09-26 (see that section) |
| B. Noise-aware aspects | How many separable molecular axes do the histology components carry once each score is a noisy ordinal readout? | `scripts/analysis/aspect_latent/PRESPEC_B.md` (hash in `PRESPEC.sha256`; Amendments 1–3) | **closed 2026-09-26** at the SVI result (user decision; see the 2026-09-26 section) |
| C. H&E to molecular | Can program scores be read from H&E beyond fat content and postmortem damage? | `scripts/analysis/he_molecular/PRESPEC_C.md` (hashed with the donor seal and negative-control lists; Amendments 1–5) | **sealed test done** 2026-09-25 (job 0950, user go-ahead) |

## Findings so far

**Bulk count undercount (affects adopted numbers).** The bulk STAR BAMs were
aligned to a 706-contig GRCh38.p14 (681 ALT/patch/scaffold contigs, 203 Mb), and
featureCounts drops multimappers. 4,115 primary genes (1,483 protein-coding) that
share a name with an ALT/patch gene keep about a quarter of their reads (median
multimapper-inclusive/default ratio 4.2–4.5; 2.0 in GSE193066; 1.00 for other
genes; `alt-contig-count-check-20260923T234550Z`, job 0908). In the official DEG
table APOB has average log-CPM −1.18 (APOA1 12.1) and COL3A1 −1.92 (COL1A1 5.40);
PKLR, TDO2 and GLUD2 are not tested. The five-cohort update shows the same.
Rebuild of the adopted five-cohort pooled model with ALT copies folded back
(`alt-contig-recount-f-five-rerun-20260925T203044Z`, jobs 0912/0966; 0913 failed
because `03_integrate_counts.R` writes `provenance/` at line 142 before creating it
at line 171). The unchanged counts reproduce the adopted table exactly (1,347
canonical DEGs, max logFC difference 1e-14). With ALT copies folded back:
- 1,419 canonical DEGs (89 gained, 17 lost); 967 more genes pass the expression
  filter (958 of them ALT-copy genes), 24 fewer;
- APOB average log-CPM −1.04 → 13.6, COL3A1 −1.93 → 7.49 (logFC 0.74 → 0.77);
  PKLR, TDO2 and GLUD2 are now tested (logFC 0.16, −0.27, −0.04; none a DEG);
- of the gained DEGs, 78 were not tested before (immunoglobulin, TRBC2, LSP1,
  LAIR1, DEFB1, CFHR5, MUC6 and others); IGHA2 changes sign (−0.68 → +0.47);
- non-ALT genes move by at most 0.016 in logFC; the 12 non-ALT genes that change
  call (4 gained, 8 lost, e.g. FASN 0.502 → 0.499) all sit within 0.002 of the
  |logFC| = 0.5 cut.
The same correction on the corrected-controls fit (the paper candidate, 1,725
DEGs; job 0981, `r5_refit_corrected_controls_alt.R`, a copy of the refit with only
the count input and the guard's count expectations changed; its guard reproduces
the corrected-count pipeline fit gene for gene): **1,833 canonical DEGs** (1,353 up,
480 down), 24,650 genes tested (23,702 before). 128 gained (113 not tested before;
119 folded genes), 20 lost (8 folded; 2 no longer tested). The 23,324 genes without
an ALT/patch copy move by at most 0.019 in logFC; 18 of their 19 call changes sit
within 0.01 of the |logFC| = 0.5 cut. Strict NAS-0 controls: 1,888 (was 1,775);
sex-discordant removed: 1,879 (was 1,770). Gene list:
`alt-contig-recount-f-five-rerun-20260925T203044Z/paper_candidate_primary_deg_changes.tsv`.
Recorded in `docs/STATUS.md` open defects.

**Model B design limit (B1, job 0910).** The sealed cohort GSE281797 (n = 91) cannot
test component-specific signal: even a large own-component signal (0.5 on the
latent scale) has at most 31–39% power with an oracle predictor, because the
partial-correlation null sits at 0.35–0.46. All four components are declared
untestable there in advance; component separation is assessed in development only.
The K = 2 recovery test passed (cosines 0.90–0.93).

**Model C data (C1).** GTEx liver slides are open-access DICOM in the NCI Imaging
Data Commons. Liver slides and liver RNA aliquots carry different tissue codes, so
matching is by donor: 224 donors (99.4 GB, 20×, 0.494 µm/px), 45 sealed. Phikon-v2
(primary embedding; H-optimus-0 is gated) was pretrained on 13,302 GTEx slides, so
these images were likely seen during pretraining; no RNA target was.

**Model C baseline check.** The handcrafted fat-vacuole fraction (no learning)
agrees with the GTEx pathologists' notes on the 179 development donors: median
1.6% where the notes mention steatosis versus 0.3% where they do not
(Mann–Whitney p = 2e-16), and Spearman 0.82 against the noted percentage (n = 27).
All 224 slides are embedded (2,000 tiles each; job 0929).

**Model C development result (179 donors; `he-molecular-c2-*/c4_run3`, `c4_run4_residualized`).**
Held-out R² gain of ridge on Phikon-v2 slide embeddings over the baseline (fat
fraction, ischemic time, Hardy class, autolysis, RIN, age, sex); 1,000-shuffle null
within Hardy × ischemic strata; BH over 79 programs:
- primary: 74/79 programs meet the rule (q < 0.05, gain ≥ 0.05), median gain 0.15;
- but the negative controls gain more (hypoxia 0.34, immediate-early 0.20), and
  embeddings alone predict ischemic time (R² 0.62) and Hardy class (R² 0.35) beyond
  the recorded covariates, so the images carry postmortem state;
- after removing each program's fit on the two RNA negative-control scores
  (Amendment 4), 39/79 remain (median gain 0.048); 36 meet the rule in both. The
  survivors are hepatocyte metabolic and plasma-protein programs (CYP3A5, CFH,
  SAA1, MLXIPL, PLIN5, APOA1/APOE, CYP2C). Fibroblast, ECM and hypoxia programs
  mostly do not survive.
- The first development run (job 0930) had a defect (standardized trailing PCs) and
  is provenance only (Amendment 3).

**Model C sealed test (45 donors, run once; `c5_sealed_test/`).** Primary: 61/79
programs meet the rule (median gain 0.17); residualized on the RNA negative-control
scores: 30/79 (median 0.10); 29 in both. **16 programs pass all four
(development and sealed × primary and residualized)**
(`c5_sealed_test/readability_four_criteria.tsv`): hepatocyte complement (CFH),
amino-acid (GOT1), lipid droplet (PLIN5), stress (SLC38A1); resident and
lipid-associated macrophage (LIPA, PLIN5); ECM adhesion (TNC); myeloid/immune
(CD163, DOCK2); hepatocyte-like metabolic (CYP3A5, MLXIPL, CYP2C9); acute phase
(SAA1); cell motility. An immediate-early (FOS) program also passes, so linear
residualization does not remove all postmortem signal; claims stay at "readable
from GTEx postmortem H&E beyond recorded covariates", not biopsy readouts.
Exposure: Phikon-v2 likely saw these slides during pretraining (images only).

**Model C external test (TCGA-LIHC non-tumor liver, 50 donors; Amendment 6; jobs 1600-1605, 2026-10-03).**
Frozen GTEx image model, no refit, run once (`he-molecular-c6-tcga-20261003T133755Z/c6_external_test/`).
Primary endpoint not met: median partial Spearman over 14 testable programs 0.204
(within-site shuffle p = 0.103); residualized on the RNA negative controls 0.231
(p = 0.001). One of 16 programs meets the rule in both versions (Myeloid / immune,
DOCK2). The hypoxia negative control (0.313) exceeds the family median. Two programs
are untestable (coverage below 0.80). Numbers and limits are in `docs/RESULTS.md`.
The Visium check in PRESPEC_C section 7 was not run: local GSE192741 images are
~4 µm/px previews.

**Model A gate G1 (30 Geuvadis libraries vs 1kGP DNA; `in-biopsy-allelic-a3-pilot-*`).**
RNA genotype calls are accurate at panel sites (non-reference concordance 0.969;
0.975–0.980 at depth ≥ 20; het sensitivity 0.975, but 0% where the RNA allelic
fraction is beyond 0.1/0.9). Imputing untyped variants from RNA calls fails: mean
per-site r² at GTEx liver lead SNVs 0.26–0.39 by MAF bin (threshold 0.80). The user
chose (2026-09-25) to drop imputation and full WASP realignment and use hets
observed directly at exonic tag SNVs of known liver eQTL leads.

**Model A tag SNVs (T1/T3, jobs 0953/0956).** GTEx v8 liver eGenes with SNV leads;
tags are exonic biallelic SNVs with r² ≥ 0.8 to the lead in both unrelated EUR and
EAS 1kGP founders (same LD sign, MAF ≥ 0.05 in both): 4,390 tag SNVs. A simulated
100-bp read test on the primary-chromosome index keeps 3,894 whose two alleles
both map back uniquely at equal rates; they cover 1,191 eGenes.

**Model A composition (job 0958).** InstaPrism initial posterior on kallisto
counts for all 1,281 libraries of the 9 cohorts in the kallisto matrix; the 521
libraries already in `deconv_estimand_harmonization-20260828T141918Z` reproduce to
5e-16. PRJNA512027 has no kallisto counts.

**Sample identity (A1 smoke).** Expression sex disagrees with metadata for the
three known GSE130970 libraries and for 13 GSE193066 libraries; seven GSE193066
libraries express both XIST and Y genes, consistent with RNA from two people.
The full identity panel adds a genotype-based contamination estimate.

## Model A-final design and Model B closure (2026-09-26)

**Design selection.** Five revised designs were scored before the choice. The
design-selection scores (judge totals) below are copied from the approved plan
`~/.claude/plans/plan-thoroughly-to-build-imperative-wren.md` (lines 16-24); they
rank designs and say nothing about any model's result:

| Design | Judge total |
|---|---:|
| Model A, κ/τ (skeptic version) | 100 |
| Model A, λ_hep version | 99 |
| Model D, within-stage fibrosis decomposition | 98 |
| A-T7, GTEx transport | 88 |
| PRESPEC_O, cross-layer coupling | 80 |

Model A-final combines the two Model A versions: κ_U (proportional difference in
the allelic effect per stage step, cross-sectional) and τ (between-gene SD of stage slopes) are
co-primary, and κ_hep (hepatocyte-predominant eGenes) is printed beside κ_U. It
won because it is the only design that links inherited cis variation to diseased
liver at the participant level: the two alleles in one biopsy share stage, cohort,
batch, library complexity and trans effects. Design detail, gates and exclusions
are in the plan; nothing from them is adopted.

**User decisions (2026-09-26; plan lines 49-61).**
1. Consent: proceed under the earlier approval. Genotypes stay restricted; only
   per-variant summaries with at least 10 heterozygotes are published; PRESPEC_A
   states the GEO-consent limit.
2. Margin: Δ = 0.125 per stage step, a 25% difference in the allelic effect
   between F0-1 and F3-4 biopsies (two steps; cross-sectional). τ_max is stated on the natural-logit scale from the
   development genes' median |α|.
3. Model D: fallback only. Its label-reliability layer is built now; the full
   decomposition runs only if Gate A fails.
4. Model B: closed at the SVI result below. No loadings, no sealed test.

**A1 identity, full panel (job 0989; `Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/identity/identity_summary.json`).**
1,466 libraries, 1,417 with enough genotype sites, 1,349 genetic individuals
(including 105 GSE193066 individuals, provisional until the contamination-aware
rerun; `individuals_per_cohort`). Same-person pairs within cohorts: GSE213621 25,
GSE193066 36 (provisional; see the possibly-mixed rule below), GSE135251 5,
GSE130970 3; none between cohorts. First-degree pairs 10, second-degree 4.
Expression sex disagrees with metadata in 23 libraries and is ambiguous in 57.
Three consequences for Resource numbers are recorded as `docs/STATUS.md` open
defects. Pairs come from `identity/kinship_related_pairs.tsv`. People counts map
each library to its `individual_id` in `identity/sample_to_individual.tsv`, count
distinct individuals, and add each library without a call as a separate person, so
they are upper bounds. The 844-fit pairs were counted 2026-09-26; the Figure 3D,
stage-contrast and Model B counts were counted 2026-09-29.
- the 844-participant bulk fit holds 32 same-person pairs and 1 first-degree pair;
  its 844 libraries are at most 813 people. Its Figure 3D roster (416) is 404
  people. In the stage contrasts on the same libraries, the adjacent fibrosis
  rosters (209, 195, 165, 100) are 206, 192, 161 and 97 people, the
  fibrosis-versus-F0 rosters (209, 194, 179, 129) are 206, 193, 176 and 129, and the
  adjacent NAS rosters (111, 186, 296) are 111, 186 and 289; `docs/RESULTS.md:646-658`
  prints all of these as participants. Four people have libraries on both sides of
  NAS 3–4→5–8 and one person on both sides of F0→F1 and F1 vs F0 (accessions and
  manifest lines in the STATUS row)
  (`figures/candidates/bulk-release-corrected-controls-20260925T175215Z/`
  `source_tables/fig3d_endpoint_pca_coordinates.tsv` and
  `analysis/stage_extensions/stage_extension_sample_manifest.tsv`);
- 11 of the 25 GSE213621 same-person pairs carry different stage labels;
- 132 of 164 GSE193066 libraries meet the A2 possibly-mixed rule: expression sex
  "both" (XIST log2 CPM > 5 and mean Y-gene log2 CPM > 3;
  `scripts/analysis/in_biopsy_allelic/a1_identity/03_identity.py:61-63`), or an
  other-allele fraction above 0.02 at homozygous sites. 130 libraries exceed the
  0.02 fraction and 24 express both XIST and Y genes. The rule and the 132 are
  recorded as run in `rules.possibly_mixed` of
  `in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260926T155512Z/a2_summary.json`
  (job 1000) and of `seal-20260929T042934Z/a2_summary.json` (job 1002, same rule and
  count; job-1000 script copy `seal-20260926T155512Z/a2_seal_holdout.job21926005.py`,
  sha256 `1d347064…`).

**Model B final numbers.** Executions under `Analysis/MASLD_Model_Benchmark/executions/`.
- Test 1, K = 2 against K = 1, corrected inputs
  (`aspect-latent-b3-tests-corr-20260924T043421Z/summary.json` and `observed.json`;
  jobs 0941-0983). Held-out log predictive density gain 0.0280 nats per participant
  record over 1,007 participant records in 7 cohorts. By genotype these are at most
  975 people: GSE213621 has 25, GSE135251 4 and GSE130970 3 repeat libraries; the
  provisional GSE193066 calls add 6 more (at most 969); the 99 GSE267145 records are
  not in the identity panel and 10 GSE193066 records have no call
  (`aspect-latent-b2-inputs-20260924T003835Z/b_inputs_corr.tsv.gz` joined to
  `in-biopsy-allelic-a1-full-20260924T025858Z/identity/sample_to_individual.tsv` and
  `kinship_related_pairs.tsv`). The 299 GSE213621 records are 274 people. The
  records hold 39 same-person pairs. In the Model B fibrosis labels 11 GSE213621,
  2 GSE135251 and 3 GSE193066 (provisional) pairs differ. Counting any Model B label
  (steatosis, ballooning, inflammation, NAS, fibrosis), 22 of the 39 pairs differ:
  GSE213621 11, GSE135251 5 (NAS in 4), GSE130970 2 (steatosis 2, inflammation 2,
  ballooning 1, fibrosis 0) and GSE193066 4 (provisional; NAS 4, fibrosis 3). Among
  the 898 records with an identity call, no same-person pair crosses cohorts, so for
  those records no person is in both the training and the held-out cohort of a fold;
  the 99 GSE267145 records and the 10 uncalled GSE193066 records could not be
  checked. PRESPEC_B section 2 sets the participant (first biopsy) as the unit, but
  the inputs keep these repeat libraries, and the per-record gain counts each as a
  separate participant. Parametric-bootstrap null (200 draws): mean
  −0.0052, 95th percentile 0.0124. p = 0.005, which is 1/201, the smallest value 200
  draws can give. K = 2 kept. The gain is not present in every cohort: GSE135251
  −0.0075 and GSE162694 −0.0038 per participant; the other five cohorts 0.031 to 0.067.
  The largest of the five is GSE193066 (0.067 over 106 records), where 132 of 164
  libraries meet the possibly-mixed rule. Without GSE193066 the gain is positive in
  4 of 6 cohorts, and the record-weighted gain over the other 901 records is 0.023
  nats (computed from the per-cohort `n` and `gain_per_participant` in
  `observed.json`; no null was drawn for this subset).
- Test 2, same file: K = 3 against K = 2 −0.110 (p 1.0) and K = 4 against K = 3
  −0.014 (p 0.398), in GSE130970 (n 78) and GSE267145 (n 99). K = 3 and K = 4
  rejected.
- Original-count sensitivity, Amendment 1
  (`aspect-latent-b3-tests-orig-20260926T014659Z/summary.json`; jobs 0984-0986):
  K = 2 gain 0.0217 (p 0.010, 2/201); K = 3 −0.079 (p 0.995); K = 4 −0.025
  (p 0.498). Same decisions.
- Test 3, component-specific signal
  (`aspect-latent-b5-test3-corr-20260926T140056Z/summary.json`; jobs 0990-0992).
  Statistic: Spearman partial correlation between each label's K = 2 prediction and
  the observed label given the other three labels, leave-one-cohort-out between
  GSE130970 and GSE267145. Null: 1,000 datasets from the K = 1 fit. Observed against
  null mean: fibrosis 0.111 against 0.251 (p 0.972), steatosis 0.322 against 0.314
  (p 0.465), inflammation 0.239 against 0.328 (p 0.889), ballooning 0.202 against
  0.351 (p 0.981); BH q 0.981 for all four. Test 3 did not detect
  component-specific signal. For three labels (fibrosis, inflammation, ballooning)
  the observed value is below the K = 1 null mean, so this is not evidence that the
  components lack signal of their own.
- Final NUTS fit, K = 2, 4 chains of 1,000 warmup and 1,000 draws
  (`aspect-latent-b4-final-20260926T135812Z/nuts_K2_summary.json` and
  `nuts_diagnostics.tsv`; job 0988). 1,007 participant records (same records as
  test 1), 0 divergences, but max
  R-hat 4.55 (loading `lam_2_1`; program weights `W` max R-hat 3.19, minimum
  effective sample size 2.2). The chains did not converge, so no per-program
  loadings, label-coefficient cosines or programs-excluding-zero counts from this
  fit are reported.

**Model B closure (user decision 4).** Model B stops here. What it supports: in
held-out likelihood, the histology labels are described better by two molecular
axes than by one (gain positive in 5 of 7 cohorts, one of them GSE193066; 4 of 6
without it), and not better by three or four; the two axes were not
characterized, because the NUTS fit did not converge.
Test 3 did not detect component-specific signal, and that is not evidence of its
absence. Under PRESPEC_B section 7 no component can be `supported` (test 3 failed
in development) or `tested_negative` (K = 1 was not kept, and GSE281797 was
declared untestable for all four). The proposed section 7 state for
component-specific attribution is `untestable` (component minimum detectable effect
above the declared limit in GSE281797 for all four components). In paper text the
resolved aspect is recorded as indeterminate (`docs/PAPER.md:530`); that is the
paper's wording, not the separate section 7 state `indeterminate`. It enters only as a candidate qualifier: steatosis,
inflammation, ballooning and fibrosis scores are not read as four separate
molecular readouts, and which component a molecular difference belongs to is
indeterminate. Registry row:
`Analysis/MASLD_Model_Benchmark/MODELS.md`, "Candidates 2026-09 (not adopted)".
Nothing is adopted without the user's decision.

## Infrastructure notes

- New primary-chromosome STAR index: `~/reference_genome/grch38_primary_chr_gencode_v49/`.
- New environments: `allelic_qtl` (tensorQTL 1.0.10, NumPyro 0.22, pandas-plink),
  `he_molecular` (OpenSlide 4.0.1 with DICOM, timm, transformers, idc-index).
- Downloads under `data/external/allelic_refs/` (1kGP high-coverage panel,
  REDIportal v3, GTEx v8 liver eQTL, Beagle 5.4 and GRCh38 map), 30 Geuvadis
  RNA-seq libraries under `data/external/geuvadis/`, GTEx liver slides under
  `data/external/gtex_liver_slides/`.
- Per-person genotypes only under `RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/` (chmod 750, git-ignored).
- Torch 2.6 cu124 has no kernels for some GPUs on this cluster; GPU jobs use `--constraint=l40s`.
- `ls | head` inside `set -o pipefail` scripts exits 141 (SIGPIPE); use `awk 'NR<=n'`.
- bcftools rejects `bam##idx##bai` in a `-b` list; reach BAMs through symlinks with the index beside them.
- Calling on the bulk BAMs: `-d 100000` and identity SNVs in the most abundant liver
  genes made one job seek-bound (≈75 MB/s, no output after 2 h). The identity panel
  now uses genes with median CPM 20–500, `-d 1000`, 9,976 sites thinned to 20 kb,
  and a 100-task array (cohort × genome chunk).
- REDIportal `TABLE1_hg38_v3`: chromosome and position are columns 2 and 3 (column 1 is an accession).
- `sbatch` keeps the script it was given; editing a submitted batch script does not change a pending job.
- `/scratch` is node-local: a job log written there is not visible from the login node.
- In `unified_metadata.csv`, GSE213621 `fibrosis_stage` is an ordinal code with
  Control = 0 (0/1/2/3 = Control/F0F1/F2/F3F4), not a Kleiner stage.
