# Kumak, Darde & Konu 2026: MASLD-META NETWORK

**Analyzed:** 2026-09-07  
**Status:** competitive assessment. Not a result authority.

## Identity

Kumak E, Darde T, Konu O. *Network-based meta-analysis maps stage-dependent
molecular programs in MASLD through MASLD-META NETWORK application.*
bioRxiv 10.64898/2026.08.26.747338, v1, posted 2026-08-31, CC-BY-NC,
category bioinformatics. Corresponding author Ozlen Konu, Bilkent University.
Thomas Darde is Scilicium, which hosts the portal at
https://masld.scilicium.com. Academic plus commercial-portal collaboration.

## What they did

Nine GEO bulk liver-biopsy RNA-seq cohorts with NAS annotation: GSE135251,
GSE193066, GSE225740, GSE130970, GSE192959, GSE207310, GSE162694, GSE174478,
GSE185051. NAS n = 897, fibrosis n = 856.

Four parallel per-study analyses, then a frequency vote across them:

1. `lm()` of VST counts on NAS group and on fibrosis stage, per study, BH within study.
2. GAMMA rank correlation plus Self-Organizing Map gene-disease patterns (Aibar 2016).
3. Per-study DESeq2 pairwise contrasts, combined with metaRNASeq inverse-normal
   and Fisher meta-p, with a Wasserstein-distance study-inclusion filter and a
   requirement of consistent logFC direction across all studies.
4. Global arm: ComBat-seq batch correction with `study` as batch, then DESeq2 on
   the merged matrix.

NAS bins were set by PERMANOVA majority vote across studies: group 0 = NAS 0-2,
1 = NAS 3, 2 = NAS 4, 3 = NAS 5, 4 = NAS 6-8. Fibrosis baseline is F0+F1,
contrasts F2, F3, F4.

Total score = number of occurrences across the four analyses. A 0.95 quantile
cut gives 951 NAS genes from 14,713 and 1,028 fibrosis genes from 20,119; a
0.99 cut gives 163 and 240.

Downstream: MSigDB over-representation and DoRothEA TF enrichment; STRING v12
PPI at confidence >= 0.4; a combined gene score equal-weighting the standardized
meta-analysis score and degree centrality; pathway-pathway Jaccard >= 0.2
similarity networks partitioned by Louvain (1,029 NAS and 881 fibrosis enriched
pathways, 34 modules each); GTEx v10 median TPM as a four-category tissue
lookup; and reassessment of 156 genes compiled from six published MASLD panels.

Headline genes: COL1A1, COL3A1, THBS2, FBLN5, PDGFA as shared central genes;
SULF2, MMP14, IL32, GPNMB, COL3A1 as differentially expressed in the earliest
NAS contrast; LAMA2 and LAMA3 as "previously unrecognized" central candidates.

The R Shiny portal adds user gene-list upload with logFC filtering, DIOPT
mouse and zebrafish ortholog mapping, adjustable centrality metrics and pruning
thresholds, sex-stratified networks, per-module centrality, and a transcriptome
browser across all genes.

## What is genuinely good

- Stage-resolved module emergence and reorganization is a legible narrative that
  a hepatologist can read directly off the figures.
- The panel reassessment is honest and useful: of 156 genes from six published
  panels, only 5 were reported by more than one study and none by more than two.
- Per-module centrality, rather than centrality in one global network, is a real
  increment over prior portal papers.
- Sex-stratified networks and cross-species ortholog mapping make the portal
  usable for bench groups working in mouse or zebrafish.

## Where the analysis is weak

Ranked by what the defect does to a claim.

1. **Repeat biopsies counted as independent participants.** GSE193066 is 164
   samples from roughly 105 donors, with second biopsies a median two years
   apart. Nothing in their Methods removes repeats: the `lm()` arm, the
   per-study DESeq2 arm, and the ComBat-seq merged arm all treat samples as
   rows. Their n = 897 therefore overcounts donors, and every p-value in the
   contributing study is anticonservative. We made this exact error in June 2026
   and retired it; the corrected design uses one pass-QC first biopsy per
   GSE193066 donor, 105 eligible (`docs/reference/progression_atlas.md`).

2. **The ranking statistic is confounded with literature attention.** The
   combined score equal-weights the meta score and STRING degree centrality.
   Degree in a curated PPI database is a function of how much a gene has been
   studied. Their own PubMed audit measures that confound and reports it as a
   result: all 20 top NAS hubs and all 20 top fibrosis hubs were already
   published in MASLD. The genes they call novel (LAMA3, LAMA2, DHRS9, TTYH3,
   BBOF1) are precisely the ones the centrality term demotes, several of which
   they note do not connect to the main component. There is no degree-preserving
   or annotation-matched null, so "most central" is not separable from
   "most annotated."

3. **The composite score has no calibrated error rate.** The total score sums
   occurrences across four analyses run on the same samples carrying largely the
   same signal, then cuts at an arbitrary 5% and 1% quantile. There is no
   permutation of the composite, no FDR on it, and no held-out cohort. The 951
   and 1,028 gene sets are ranked lists, not error-controlled sets. This is the
   convergent-score pattern already flagged as a crowded claim in
   `docs/literature/REVIEW.md`.

4. **Outcome bins were chosen from the same expression data later tested.**
   PERMANOVA on VST counts decided which NAS values to merge, and the resulting
   groups were then tested for differential expression on those same matrices,
   with no held-out split and no permutation of the selection step. The NAS
   contrasts are optimistically biased by an unquantified amount.

5. **The study-inclusion filter selects against signal.** Only studies whose
   p-value distribution sat within the Wasserstein distance of the most extreme
   contrast were retained, that is, studies whose p-histograms look more uniform
   and therefore more null-like. In a meta-analysis this preferentially drops
   the most informative cohorts, and the retained study set changes per
   contrast. No sensitivity analysis is reported.

6. **Nominal p at the meta step.** As written, genes are retained at "a p-value
   below 0.05 in both of the meta-p-values" across roughly 15,000 to 20,000
   genes, with direction consistency as the only additional constraint. The text
   says p-value, not adjusted p-value. [Inference: supplementary code not
   inspected. If metaRNASeq's adjusted p was used, the Methods do not say so.]

7. **ComBat-seq then DESeq2.** Feeding batch-corrected counts into a model that
   does not carry the batch term removes batch twice and understates residual
   variance. The standard alternative is a `~ study + condition` design. Both
   this arm and the per-study meta arm feed the frequency score, so a gene can
   accumulate score from the anticonservative arm.

8. **No confounder adjustment in the primary regression.** Expression on NAS
   with no covariates, not even sex. They state the cause: metadata across the
   nine cohorts did not permit controlling for age, diabetes, or medication. Sex
   enters only as a stratified supplement.

9. **The baseline group is composed differently in each cohort.** GSE174478
   (94 samples, Japanese NAFLD) and GSE193066 (164 samples, HCC-naive NAFLD)
   have no healthy controls; both are excluded from our mega-analysis for that
   reason (`docs/reference/dataset_notes.md`). Their NAS 0-2 and F0+F1 baselines
   therefore mix control-derived and disease-derived low-grade samples in
   cohort-dependent proportions.

10. **Figure callouts are swapped.** The Results text cites Fig. 2A-G for the
    stage-resolved pathway networks and Fig. 3A-F for the panel refinement,
    while the captions assign those contents to Figure 3 and Figure 2
    respectively. A reader cannot follow the Results as printed.

## Overlap with our paper

**Data.** Three of our five canonical bulk cohorts are in their nine: GSE135251,
GSE130970, GSE162694. They do not have GSE213621 (n = 361, our largest) or
GSE126848. We hold GSE174478 and GSE193066 but exclude both from the
mega-analysis. We do not use GSE225740, GSE192959, GSE207310, or GSE185051.
Their n = 897 against our 844 pooled participants is a tie, not a gap. Do not
compete on sample count.

**Analytic overlap, by figure.**

| Our figure | Their competing content | Exposure |
|---|---|---|
| Figure 3 (multicohort DE, per-cohort membership, adjacent-stage contrasts) | Per-stage DE meta-analysis across nine cohorts with NAS and fibrosis axes | Highest. Directly competing on the same axis with comparable n |
| Figure 4G/4H (43 size-ranked pathway and program systems) | 34 Louvain pathway-similarity modules per axis, stage-resolved | High. Different construction, same reader takeaway |
| Continuum layer (within-stage remodeling) | "Modules emerge and reorganize across increasing severity," citing Kamzolas 2026 | Moderate. Both are downstream of Kamzolas; ours is donor-level and stage-adjusted, theirs is contrast-level |
| Portal | masld.scilicium.com, joining SteatoSITE, MegaMASLD, LiRNA | A MASLD portal is now a commodity. It cannot carry novelty |
| Stromal ECM program biology | COL1A1, COL3A1, COL1A2, THBS1, THBS2, SPP1, TIMP1, MMP2, GPNMB, LGALS3, TREM2, MMP14 as central | Moderate. Same biology, no lineage or spatial resolution |
| GNMT-MAT1A-CYP2C19 example | F3/F4 "one-carbon and folate metabolism" and "amino acid and nitrogen transfer" modules down; HNF4A and CEBPA enriched among downregulated genes | Adjacent but they never name GNMT, MAT1A, or CYP2C19. They reached the pathway and stopped |

**What they have none of.** Verified by string search over the full text:
GWAS, variant, eQTL, colocalization, germline, and heritability appear zero
times, as do PNPLA3, TM6SF2, HSD17B13, and MBOAT7. Single-cell, spatial and
proteomics appear once each, in the limitations paragraph. GNMT, MAT1A,
CYP2C19, BICC1 and IGFBP7 appear zero times.

## What we still hold

1. **The germline-to-remodeling link.** They have no genetics at all. Kamzolas
   has no genetics. The Gaulton-group regulatory-risk paper has genetics but no
   multicohort remodeling map. That intersection is still empty.
2. **A verifiable donor unit of inference.** They pseudoreplicate GSE193066; we
   diagnosed and corrected the same defect with recorded numbers.
3. **Outcome-blind module definition.** Our 117 Hotspot programs are frozen with
   fixed weights before any outcome is seen. Their 34 modules are built from
   over-representation of genes that were themselves selected on the outcome, so
   the module structure is outcome-selected end to end.
4. **Assay-calibrated tissue claims.** Stromal ECM exceeds matched spatial
   background in both Visium sources; Ductular injury does not. They have no
   tissue context and cannot separate disease association from tissue
   organization or from lineage origin.
5. **Machine-readable evidence states.** Their limitations are a prose paragraph.
   Ours are per-object states with observability accounting: CosMx adequately
   observes 12 of 117 programs; 38,066 of 39,914 jointly testable ATAC peaks are
   indeterminate.
6. **Molecular-object separation.** Coding variant, noncoding DNA, transcription
   through a lncRNA locus, and the lncRNA molecule are distinct objects here and
   are one undifferentiated gene symbol there.
7. **Figure 7 model-exposure auditing.** Nothing comparable exists in this space.

## What they do better, or first

- They posted 2026-08-31. We are unreleased.
- Their portal offers user gene-list upload with direction filtering and DIOPT
  ortholog mapping for mouse and zebrafish. A bench group with a mouse screen
  can use their site today.
- Per-module centrality gives readers a within-module priority our systems layer
  does not currently offer.

## Strategy

**One recommendation: do not add analyses. Cite them, explicitly disclaim the
overlapping layers, and move the paper's weight onto the germline link and the
donor/assay calibration they structurally cannot do.**

The real risk is not being scooped. It is a reviewer writing "the multicohort
differential expression and stage-resolved pathway module layer is already
published." The mitigation is a citation plus an explicit non-claim, not new
compute.

Concrete moves, in order of value:

1. **One defensive number, cheap to get.** Intersect their Supplementary Data 1
   (951 NAS and 1,028 fibrosis top-scoring genes) against our 1,347 canonical
   DEGs and against our 462 SuSiE genetic set. Two results are needed before
   submission: how much of their list we recover on the shared axis, which shows
   we are not weaker there; and how little of their list touches the genetic map,
   which shows the inherited axis is genuinely orthogonal and unreachable from
   bulk expression alone. This is the highest-value new number in this analysis.

2. **Make the donor unit visible in Figure 1.** Add donor count beside sample
   count for every cohort, and name the GSE193066 decision (164 samples,
   105 first biopsies) in the Methods. This converts our rigor into a checkable
   difference from the published record at zero analytical cost.

3. **Promote the GNMT bridge to the paper's lead example.** It is the one place
   on their own axis that they approached and could not reach. They found
   one-carbon and amino acid catabolism modules falling at F3/F4 and HNF4A and
   CEBPA enriched among downregulated genes. We name GNMT, MAT1A and CYP2C19,
   show agreement across pooled RNA, adjacent fibrosis stage and liver protein,
   and attach an inherited signal at the GNMT locus. That contrast is the
   sentence a reviewer remembers.

4. **Add one Discussion sentence on centrality, non-adversarially.** Network
   centrality in a curated interaction database is a function of annotation
   depth, which is why our programs are frozen with fixed weights before outcome
   and their per-gene contributions are stated (IGFBP7 is 2.47% of the Stromal
   ECM L1 weight; BICC1 is 5.44% of Ductular injury). No need to name them.

5. **Cite Kumak in the Figure 3 and Figure 4G text and say what they
   establish.** Never present a multicohort DE gene list or a stage-resolved
   pathway module as our discovery.

6. **Portal scope.** Confirm our gene query returns evidence state, source
   dependence, and next experiment. That is the differentiator. Do not chase
   their centrality sliders or ortholog mapping.

7. **Do not build on their LAMA2/LAMA3 novelty claim.** It rests on an
   easyPubMed keyword query, which is a literature-search negative, not a
   biological one. If either gene carries independent evidence in our Catalog,
   that is a useful contrast; otherwise leave it.

**Journal implication.** This is a portal and meta-analysis paper. It is not a
Nature Metabolism competitor and does not change our target. It does raise the
bar on what the transcriptomic core of Figure 3 can claim on its own.
