# Methods — working Resource-paper draft

> `docs/PAPER.md` governs interpretation. The named
> `2026-07-15-r2` package is a historical analysis release, not the final
> standalone Resource release. BG-001 adjudication and Plan 60 adoption must
> precede submission; no method below authorizes Cas13-screen evidence or a
> mechanism-escalation branch.

## Analysis release and provenance

The latest named historical analysis release is `2026-07-15-r2`. Source paths, file hashes, modification times, thresholds, and roles for that release are recorded in `docs/manuscript/release/analysis_release_manifest.tsv`; its claim ledger governs only that read-only release. The final manuscript will cite a new Plan 60 read-only Resource release containing the selected bulk substrate, complete program registry, source graph, MASLD Gene Catalog, figure sources, and consistency products.

The r2 cut regenerated only part of the former manuscript package. Drug-target calibration and credible-set coding-architecture products retained from r1 are supplementary historical stress tests and cannot define a main claim. The two releases also differ in proteomic and single-cell evidence-class results; r1 values for those endpoints are superseded. Plan 60 must snapshot one internally consistent substrate rather than combine release versions at figure-render time.

The release distinguishes primary definitions from sensitivity analyses. Primary transcriptomic evidence is the canonical pooled limma-voom quality-weighted result under a conventional adjusted-p-value and effect-size threshold, with the effect-size-aware interval-null threshold retained as the calibrated comparator. The genetic analysis has two explicitly different coverage layers across 35 prespecified Tier-1/2 GWAS strata: approximate-Bayes-factor (ABF) colocalization is evaluated for every testable gene--GWAS pair, whereas multi-signal SuSiE colocalization is available only when the precomputed eQTL fit exists and converges and the GWAS-side LD and convergence checks pass. The fixed evidence-class interface uses the high-resolution SuSiE call; the complete ABF map and the SuSiE-or-ABF union are reported alongside it with their own tested denominators. The expanded 50-GWAS portfolio, broad prioritized gene universe, and alternate transcriptomic contrasts remain sensitivities and are not substituted into main claims.

## Human bulk RNA-sequencing cohorts

Nine human MASLD cohorts were collected and uniformly processed. A total of 1,260 samples passed technical quality control. Five cohorts containing healthy controls contributed 846 samples to the canonical pooled disease-versus-control analysis. PRJNA512027 was retained for provenance but excluded from the pooled paper analysis because library preparation was confounded with case-control status.

Technical inclusion used the `pass_technical` flag, defined by PCA and library-size quality control. Samples failing only the expression-based sex check remained eligible for non-sex analyses. Sex-stratified analyses applied the additional sex-concordance filter separately.

Reads were quantified against GENCODE v49 on GRCh38.p14 using the repository's canonical STAR and featureCounts workflow with reverse-stranded counting. Cohort metadata and exclusions were read from `config/human_datasets.yaml`.

## Canonical disease-versus-control differential expression

The canonical model was fitted with limma-voom quality weights using the design

`~ dataset + inferred_sex + group_binary`

where dataset was a fixed effect. The model included 846 samples from five control-bearing cohorts in the current named release, and 844 in the selected fragment-count substrate. The primary differential-expression threshold was a Benjamini-Hochberg adjusted p-value below 0.05 together with an absolute log2 fold change above 0.50. An effect-size-aware interval-null threshold, implemented with `limma::treat` at a log2-fold-change threshold of 0.25 and BH FDR below 0.05, was computed on the same fit and is reported throughout as the calibrated comparator.

The selected fragment-count substrate contains 1,347 records (1,003 higher and
344 lower in disease) under the primary threshold and 1,616 records (1,144 higher and
472 lower) under the interval-null comparator, all within the same 23,370-gene
universe from 844 participants. The candidate is validated but not adopted;
all final genetics–expression overlap values remain withheld until the corrected
COLOC release and the fragment-count figure package are synchronized.

Both thresholds are reported because applying a fold-change filter after a
point-null test does not preserve error-rate control, whereas the interval null
incorporates the effect threshold into the test (McCarthy & Smyth, 2009).
Duplicate gene symbols are resolved outcome-blindly within the adopted
substrate. Raw point-null, ashr, and superseded model results remain
sensitivities only.

## lncRNA identity, measurement, and disease association

Versioned Ensembl identifiers from the exact GENCODE v49 annotation used for
counting were the identity key for both protein-coding genes and lncRNAs;
symbols were display aliases only. Small RNAs, circRNAs, pseudogenes, and
transcript-isoform claims were excluded from the primary lncRNA analysis. The
five pooled cohorts were required to share the same reverse-stranded featureCounts
setting. Strand orientation was audited from deposited STAR summaries and,
for GSE162694, from a deterministic ten-BAM recount using featureCounts modes
0, 1, and 2. The approved cohort-level criterion required a median reverse-to-same
assignment ratio of at least four, at least 90% of audited libraries at four or
above, and every library above one. GSE213621 additionally required no clinical
stage group to have more than 10% of libraries below four and a maximum
between-group failure-fraction difference of five percentage points. No sample
was excluded on the basis of this audit.

lncRNAs were tested in the same limma-voom/TREAT model and the same 23,370-gene
Benjamini-Hochberg family as every other gene. A separate lncRNA-only FDR was
not calculated. Cohort-specific, leave-one-cohort-out, and equal-library-weight
fits were prespecified robustness analyses. Named lncRNA examples additionally
required a canonical chromosome, an unambiguous genomic mapping, concordant
direction in at least four cohorts, no materially opposite fifth-cohort
confidence interval, concordant direction in all five leave-one-cohort-out
fits, retained TREAT significance in at least four such fits, and concordant
equal-weight results. The secondary association between biotype and DEG status
used a gene-level logistic model adjusted for average expression, expression
variability, log gene length, and the number of cohorts in which the gene was
detected. Genes were not described as biological replicates in this descriptive
comparison.

Each lncRNA was assigned one deterministic genomic-context class using exon
overlap, gene-body overlap, strand, and opposite-strand TSS proximity to
protein-coding genes, in that priority order. Same-strand exonic-overlap and
complex mappings remained in coverage tables but could not become named
examples. Genomic context was not interpreted as evidence that the mature RNA
molecule was functional.

Prespecified Hotspot programs with at least 5% lncRNA L1 weight were rescored after
removing all uniquely mapped lncRNA members, without changing module membership
for other genes or rediscovering modules. Because exact original cell scores but
not the original neighbor graph were serialized, the source graph was
reconstructed from the original cell-type matrix, source filters, Hotspot
model, neighbor count, and gene filter. The exact stored score was retained as
the baseline and only the reconstructed gene-removal delta was added. The
technical reconstruction required correlation of at least 0.995 with the
stored all-gene score; this threshold was finalized during source audit and was
kept separate from biological acceptance. Scores were collapsed to biological
donor. A sensitivity passed when donor-level Pearson and Spearman correlations
were both at least 0.90 and the coefficient direction in the prespecified model
`score ~ stage_ordinal + dataset` was preserved.

## Fibrosis-stage and disease-severity analyses

Adjacent fibrosis stages and adjacent NAS groups were compared
cross-sectionally in the selected five-cohort fragment-count substrate. Each
analysis retained one record per participant and only cohorts carrying both
compared groups. The design was `~ dataset + inferred_sex + group` when the
multi-cohort terms were estimable. Counts were normalized with edgeR and
analyzed with `limma::voomWithQualityWeights`; empirical-Bayes P values were
adjusted by Benjamini-Hochberg separately across the complete 23,370-gene
family for every contrast. No stage-specific fold-change threshold was
applied. Direct F1, F2, F3, and F4 versus F0 models used the same substrate and
covariates. Sample count, cohort count, design rank, coefficient name, and the
complete family are serialized for every model. The five-cohort candidate
passed checksum, design-rank, family-cardinality, and FDR recomputation checks.
Arrow notation in compact panel labels denotes ordered cross-sectional
contrasts; participants were not followed longitudinally, and raw DEG counts are not
treated as comparable effect magnitudes when sample size or cohort availability
differs among contrasts.

Bulk NMF programs were interpreted as continuous axes. Cross-seed matched-factor cosine similarity was used to assess factor-subspace stability. Hard sample partitions were not treated as reproducible subtypes because sample-level adjusted Rand index was modest and no internal metric uniquely selected six factors.

## Histology-anchored molecular continuum

The molecular continuum was adapted from Kamzolas et al. (2026;
doi:10.1038/s42255-026-01543-7) as a complementary Resource evidence layer,
not as a replacement for recorded histology or the canonical disease-versus-
control differential-expression estimand. Source-overlap discovery cohorts were
GSE126848, GSE130970, and GSE135251. GSE162694 and GSE213621 were evaluated
separately as competitor-independent cohorts, and all 54 quality-controlled
paired GSE193066 participants were retained for within-person sensitivity.
Continuum and molecular scores were calculated in all five Resource cohorts.
Descriptive fixed-window figures show all five cohorts; stage-stratified
figures show GSE130970, GSE135251, GSE162694, and GSE213621 because GSE126848
lacks recorded fibrosis. Source-overlap cohorts were not added to primary
stage/sex-adjusted inference or meta-analysis.
Versioned Ensembl identifiers were reduced to GENCODE v49 base IDs. The fixed
145-gene signature contained 139 observed genes; the six unavailable genes
remained explicit in the mapping ledger.

The co-primary scores were a 145-gene cohort-specific PC1 and a fixed external
projection using prespecified discovery loadings. PC1 signs were oriented with the
prespecified discovery direction rather than target histology. A mean percentile rank
of the two scores was used for visualization only. Within GSE162694 and
GSE213621, each molecular outcome was modeled separately against each score as
`outcome_z ~ continuum_z + factor(fibrosis_stage) + inferred_sex`. Cohort
coefficients were combined by fixed-effect inverse-variance meta-analysis;
random-effects estimates and NAS substitution in GSE162694 were sensitivities.
Natural splines with three degrees of freedom were compared with the linear
model under complete-family correction. Nine prespecified 20-percentile-wide
windows with 50% overlap were calculated within each displayed cohort after
within-cohort feature standardization. They were descriptive and never used to
select features or define significance.

For the pooled molecular-systems window display, we first calculated the
collection-balanced system score separately in every cohort and window: the
median feature score was taken within each represented collection, followed by
an equal-weight mean across collections. We then averaged the five
cohort-specific system scores with equal cohort weight. Thus GSE213621 did not
receive greater weight because it contained more participants. Continuum
percentiles and feature standardization remained cohort-specific, so this was
an aligned descriptive summary rather than a jointly normalized participant
pool.

For display in Figures S4T–S4U, systems were sorted by the descriptive
equal-cohort pooled late-minus-early contrast: the mean score in windows
centered at percentiles 80 and 90 minus the mean score in windows centered at
10 and 20. The same row order was used for the five cohort-specific heatmaps,
and the contrast was shown as the final `Delta` column in the pooled heatmap.
This result-adaptive ordering was not used for feature selection, support
calls, nonlinear timing, or any inferential test.

Transcript models used count-native `limma::voomWithQualityWeights` and retained
the complete 23,370-gene Benjamini–Hochberg family. The 139 observed signature
genes were excluded from ordinary self-association inference and tested only
with target-specific leave-one-gene-out versions of both axes in a separate
139-gene family. Results were termed continuum-associated transcripts, not
continuum DEGs. Before pathway, TF-regulon, Hotspot, or fixed-basis NMF scoring,
all signature genes were removed. A pathway or program was testable only with
at least ten retained genes and at least 80% of original membership or signed
L1 weight. Multiple-testing correction used the full prespecified family: 50
Hallmark, 658 KEGG, 1,787 Reactome, 7,583 GO biological-process, 1,855 GO
molecular-function, 1,042 GO cellular-component, 268 TF, 117 Hotspot, and ten
NMF tests per declared score and model.

The ten NMF loadings retained their original k4/k6 definitions and continuous
interpretation. Inferential NMF results used signature-excluded fixed-W scores
only after full-basis reconstruction correlated at least 0.99 with released H,
at least 80% factor weight remained, and excluded versus native scores
correlated at least 0.90 in both evaluation cohorts. The 117 Hotspot programs
retained their fixed memberships and weights. Paired GSE193066 models used
`delta_outcome_z ~ delta_continuum_z + delta_fibrosis`; paired support upgraded
but did not veto adequately supported cross-sectional results. Every continuum
analysis used biological participants as the inferential unit and was described
as cross-sectional molecular remodeling, not within-person temporal ordering.

For the molecular-systems display, we first preserved a complete 13,092-row
testability ledger spanning Hallmark, KEGG, Reactome, three GO collections, and
the 117 fixed Hotspot programs. Geometry eligibility required a feature to pass
its prespecified gene-count and retained-membership criteria in both evaluation
cohorts. The 139 observed genes from the published 145-gene continuum signature
were removed before mapping member genes to GENCODE v49 identifiers. Weighted
feature-by-gene membership was inverse-document-frequency transformed, reduced
to 100 singular-vector coordinates, and used to construct a mutual 30-nearest-
neighbor cosine graph. Leiden clustering was restricted to the largest graph
component; three disconnected two-node components remained explicit as
microcomponents but were excluded from system summaries. Systems were assigned
size-ranked identifiers rather than biological names. A seeded UMAP of the
membership coordinates was used only for display. Coordinates, component
assignments, and system memberships were checksummed before any stage,
continuum, or donor-score outcome was read.

Fibrosis-stage and stage-adjusted fixed-projection coefficients were overlaid
on the same coordinate table. System summaries were calculated both as the
equal mean of within-collection medians and as the median across all member
features. Fixed-window summaries used the nine prespecified windows in all five
cohorts under both aggregation rules. Cross-system correlations, direction
counts, and window-direction concordance were treated as descriptive because
gene sets are nested, participants recur across molecular features, and
adjacent windows overlap. No independent-sign null, system-level p value,
outcome-adaptive community, or single-medoid biological system name was used.

Single-cell F-stage-specific analyses were restricted to Andrews GSE202379 samples with a documented SAF-derived F stage. Documented F stage was recoverable for 58 of the 67 sequencing runs, and those 58 runs map to exactly 46 biological donors with concordant stage labels across every multi-run donor. The accompanying scVI F-stage inference was validated at run rather than donor level and does not meet the donor-level standard applied to every other per-donor statistic reported here. scVI-predicted and disease-stage-augmented F-stage fields were therefore excluded from primary cross-cohort inference.

## Continuum extension: intervention hypotheses and measurable disease states

This implemented extension is a candidate analysis, not an adopted change to
genetic memberships, Gene Catalog evidence classes, or clinical practice. Its
code is in `scripts/analysis/histology_anchored_continuum/translation/`. Complete
tables, assay-specific contracts, statistical checks and the integrated figure
are under
`RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z/`;
the completed integration is `integrated-20260906T231500Z/`. The extension asks
whether inherited expression associations and established tissue remodeling
suggest different experiments, and which measurements could characterize
molecular variation within an existing histologic eligibility group. It does
not estimate treatment benefit or establish a clinical cutoff. The continuum
and its sliding-window adaptation remain attributed to
[Kamzolas et al.](https://www.nature.com/articles/s42255-026-01543-7).

**Inherited direction.** We checked that the adopted 50-study aggregate matched
the aggregate used for the existing signal replay. We retained all 5,362
exported SuSiE signal pairs and all 22,356,697 variant posterior records from
816 main-trait gene–study pairs. Signal support retained the existing strict
`PP.H4 > 0.5` threshold. Direct MASLD, liver fat, liver enzymes, and other traits
were separate fields; ancestry and GWAS/eQTL signal identifiers were retained.
Supported gene–trait records elsewhere in the portfolio, including ABF-only
sensitivity support, remained visible when signal posteriors were unavailable.
They did not acquire a direction from an unrelated lead variant. Rows without
a gene identifier were exported separately and could not be joined to a named
Gene Catalog object. Main-trait and full-portfolio testing denominators were
recorded separately, distinguishing missing SuSiE tests from tests without
support. These posterior support thresholds are not frequentist P values.

At matched GRCh37 chromosome and position, we matched the eQTL effect/other
allele pair to the recorded GWAS effect/other allele pair. An exact swap negated
the eQTL beta. Unmatched alleles, non-SNV alleles, ambiguous duplicate effects,
missing effects, and palindromic variants remained direction-unresolved. We did
not infer strand from an effect sign. This conservative palindrome rule affects
direction resolution only; it does not remove variants from the adopted
colocalization analysis. For each shared-variant posterior vector, the sign of
`GWAS beta × aligned eQTL beta` describes expression associated with a higher
value of the named trait. Only for direct disease phenotypes may this be
described as disease-risk-associated expression. We summed the original shared
posterior mass carrying each sign and required at least 0.95 of the original
total in one resolved direction. Missing or ambiguous mass was not
renormalized away. This is a declared descriptive direction criterion, not a
calibrated posterior probability of a causal expression effect.

The available signed effects are marginal summary-statistic effects. They do
not identify the conditional effects of individual SuSiE signals. Therefore,
directional perturbation records were resolved only for loci with one exported
GWAS signal and one exported eQTL signal; records with multiple exported
signals retain marginal-direction evidence but remain signal-direction
unresolved. Even the single-exported-pair designation does not exclude
undetected additional signals. No instrumental-variable causal effect or drug
effect was estimated. Per-study GWAS sample sizes and per-variant eQTL sample
sizes accompany the signed evidence.

We additionally tabulated consistency among resolved records for the same
gene, named trait and ancestry across studies. Conflicting signs remain in
their individual records and suspend a common directional perturbation
hypothesis. Consistent signs across studies are not called independent
replication without resolving participant overlap. Different signal locations
within a gene remain distinct possible mechanisms, not a pooled estimate.

**Expression relationship and evidence records.** We joined inherited evidence
to the complete 23,370-gene bulk testing family and retained genetics-only genes
outside that measured universe. Matching versus opposing directions required
the existing two-axis, two-cohort continuum-support rule. Genes with inherited
support but no supported remodeling were retained, and bulk-unobservable genes
were distinguished from tested genes. The 139 observed axis components retain
their separate leave-one-out results and 139-gene BH family; they do not acquire
ordinary independent continuum support from the axis they define. Opposing
directions motivate bidirectional experiments and alternative explanations;
they do not prove compensation. Lifelong genetic associations do not establish
the sign or timing of a pharmacological effect. Gene Catalog sidecars separate
the candidate transcript, perturbation hypothesis, cell-context links,
disease-associated profile, proposed readout and unresolved experiment. Shared
program membership alone does not establish regulation between a target and a
different proposed marker. No composite evidence score or target ranking was
introduced.

**Within-stage programs and composition.** We retained all 117 fixed program
memberships and weights and reused the signature-excluded participant scores.
The primary model was `program_score_z ~ continuum_z + factor(fibrosis_stage)
+ factor(inferred_sex)`, fitted separately in GSE162694 and GSE213621 for each
of the two co-primary axes. There were 109 and 361 complete participants,
respectively, for testable programs; every program has an explicit denominator
and failure state. These previously examined cohorts provide replication
evidence, not untouched prospective validation. The three source-overlap
cohorts contribute descriptive profiles only.

Composition sensitivity added 15 orthonormal log-ratio coordinates from the
16 estimated BayesPrism lineage fractions. Fractions below `1e-6` were floored
before taking logarithms; a normalized Helmert basis supplies the coordinates.
We report the frequency of floored fractions by lineage and cohort. Both models
used the same complete participants and fixed outcome units. Axis scores were
standardized within the complete cohort, and singular designs were marked
untestable. Ordinary least-squares coefficients, standard errors and two-sided
t tests were retained; HC3 standard errors and t-reference sensitivity tests
were exported separately. Fixed-effect inverse-variance summaries required
both evaluation cohorts. BH correction covered all 117 programs separately in
each cohort/axis/model family and each meta-analysis axis/model family,
including untestable members in the family size. Support requires both
co-primary meta-analysis BH q values below 0.05 and consistent signs in all
four cohort/axis fits. A separate, stronger field requires BH q below 0.05 in
each cohort and each axis. Neither criterion implies 117 independent biological
processes. Composition-adjusted associations can still reflect imperfect
deconvolution or shared expression measurements; attenuation does not estimate
a causal cell-composition contribution. We added no nonlinear test. Existing
three-degree-of-freedom spline tests and their complete-family corrections
remain unchanged. Coherence claims are excluded pending resolution of the N1
reconstruction and N2 calibration issues.

**Cell and tissue context.** Fixed single-cell source-lineage labels are
retained as discovery context, not proof of transcript origin after ambient
RNA correction or independent confirmation of a bulk RNA result. Signal-level
ATAC overlap and program-level ATAC association retain their own lineage,
coverage, donor denominators and evidence states. Spatial assay coverage is
reported for the complete program family, while formal matched-background
tests remain confined to the originally tested programs and their stated
multiplicity family. Untested observable programs do not receive a negative
spatial result. A program-level spatial result is never transferred to each
member gene. Existing paired-biopsy program and gene analyses retain their
within-person unit, complete testing families and measured changes; repeat
biopsies are not treated as technical replicates. The current integration does
not newly establish within-stage, cell-specific genetic regulation.

**Protein readouts and verified participant linkage.** Existing adjusted liver
DIA-MS disease-association results retain their log2 abundance contrast and
full-proteome BH correction, with 58 source samples (12 controls and 46 with
MASLD), adjusted for acquisition batch, age, BMI and sex. Legacy program-protein
tests retain their original subset and identifiers and are supplied separately;
they are not silently mapped to all 117 program-level tests. Complete member
protein observability is exported without treating individual member-protein
associations as a test of an entire RNA program.

An exact join through the PXD051911 source metadata identified 41 distinct
participants with initial liver and plasma measurements: 11 controls and 30
with MASLD. We verified unique participant and assay-file identifiers and
exported all included and unavailable links. In each assay separately, positive
intensities from unambiguous single-gene protein groups were log2 transformed;
groups observed in at least half of these matched participants were retained,
quantile normalized without disease labels, and duplicate gene groups collapsed
by the within-sample median. Missing abundances were not imputed. The complete
shared family contained 345 genes. For each gene, we fitted standardized plasma
abundance against standardized liver abundance, categorical fibrosis, age,
BMI, sex, liver acquisition batch and the recorded plasma batch. A separate
sensitivity omitted fibrosis. Each model required at least 20 complete
participants, full design rank and at least 10 residual degrees of freedom;
constant categorical covariates were omitted. Two-sided coefficient t tests
were BH corrected across all 345 eligible genes separately for each model,
including untestable genes in the family size. The six untestable genes remain
visible. These estimates describe participant-linked cross-compartment
association in an existing public cohort. They do not establish liver origin,
secretion, a target-to-marker regulatory mechanism, or a new independent
clinical validation cohort. The linked protein participants remain unpaired
with the bulk RNA continuum cohorts. No bulk continuum percentile was assigned
to protein or spatial samples.

**Display and reproducibility.** The integrated figure retains the original
nine overlapping descriptive windows. Its complete program map averages the
five cohort window means equally; per-cohort means and participant counts
remain available. Window membership did not select genes, transitions,
thresholds or significance. Native assay effects and their distinct error
models remain separate; repeated RNA-derived features do not count as
independent modalities. All scripts set seed `20260906`, record R session
information and write into a new candidate directory. Allele-swap, ambiguity,
posterior-mass, complete-family, participant-denominator and source-effect
reproduction checks accompany the output.

## GWAS portfolio and phenotype tiers

The primary genetic portfolio contained 35 GWAS strata listed as `placement = main` and `tier` 1 or 2 in `GWAS/finemapping/config/gwas_trait_tier.tsv`.

Tier 1 contained direct NAFLD, NASH, and PDFF traits. Tier 2 contained ALT, AST, and GGT traits. Tier-3 mixed-etiology sequelae and Tier-4 distal nonspecific traits were supplementary. Direct-disease and enzyme results were tabulated separately throughout.

The 35 primary strata covered European, African, East Asian, Admixed American, and South Asian GWAS populations. The ancestry label refers to the GWAS stratum. The Broadaway liver-eQTL panel is predominantly European, so analyses pairing non-European GWAS with this eQTL panel were described as cross-population rather than fully ancestry-matched molecular colocalization.

## Liver-eQTL colocalization

GWAS strata were integrated with the Broadaway liver-eQTL meta-analysis of 1,183 donors. Single-variant ABF colocalization was evaluated for every testable gene--GWAS pair without using LD. Multi-signal SuSiE colocalization was attempted only when a converged precomputed eQTL SuSiE fit, sufficient overlapping variants, a usable ancestry-specific GWAS LD matrix, and a converged GWAS SuSiE fit were all available. Failure or structural absence of any SuSiE prerequisite produced an explicit `abf_fallback` or `abf_only` state; it was not interpreted as a negative SuSiE result. The prespecified high-resolution support threshold was SuSiE `PP.H4 > 0.5`; the broader definition was SuSiE or ABF `PP.H4 > 0.5`.

For every final portfolio and trait tier, the release reports the complete ABF denominator, the number and fraction of pairs with a finite SuSiE result, the number of genes with at least one finite SuSiE result, method-specific positive counts, and fallback reasons. In the historical July 50-GWAS master, ABF was finite for 927,225/927,225 rows whereas SuSiE was finite for 22,072 rows (2.38%) spanning 4,945 genes. These values diagnose the historical release only; the final denominators and counts will be regenerated after the corrected 1,100-file rerun and must not be copied into the adopted manuscript from this working draft.

European analyses used the full PolyFun UK Biobank European LD reference. Non-European analyses used the configured ancestry-matched GWAS LD panel where available, but the eQTL side remained European. Colocalization was described as nominating candidate effectors and did not by itself establish a causal gene.

Rows with empty gene symbols are excluded from named-gene counts. The historical July release yielded 473 Tier-1/2 SuSiE genes and 1,030 named genes in the SuSiE-or-ABF union; the older union count of 1,031 included one blank-symbol row. **These July values are superseded. The corrected 1,100-file rerun is complete and was adopted as canonical on 2026-08-17, yielding 462 Tier-1/2 SuSiE genes, 822 Tier-1/2 ABF genes, and 1,013 named genes in the SuSiE-or-ABF union (one blank-symbol row excluded, as before); the full 50-study portfolio yields 745 SuSiE genes at PP.H4 > 0.5.** Method-availability denominators must still be regenerated from the adopted registry rather than copied from this draft.

The rerun differs from the July release in three respects, all of which affect the SuSiE arm. First, strand-ambiguous (palindromic) variants — **14.1% of all variant instances** — were previously dropped at every locus and are now retained under an ancestry-matched allele-frequency guard; across 414,720,322 testable palindromic instances the guard rejected 97,341 (0.0235%), consistent with every study in the portfolio being forward-strand. This restores PNPLA3 rs738409, itself palindromic (C/G), which had been excluded from every locus and is now the lead variant in 6 of 15 MASLD strata. Second, the eQTL SuSiE fits were regenerated with the LD panel pinned to PolyFun and provenance recorded, replacing fits whose LD reference was unrecorded. Third, the non-European 1000 Genomes LD panels were rebuilt as mean-imputed Gram matrices; this had no measurable effect on colocalization (PIP correlation ≥ 0.999998).

Because fine-mapping credible sets depend on the variant set, retaining these variants changes the SuSiE arm materially: **147 Tier-1/2 genes change their PP.H4 > 0.5 call (79 lost, 68 gained)**, against a Wakefield-ABF arm that is essentially unchanged (r = 0.990) as expected for a statistic that is invariant to allele orientation. Pinning the LD panel additionally left 204 eGenes non-convergent at the original 500-iteration cap; refitting at 5,000 iterations recovered 145, of which 107 had also failed silently in the July release. **Sixty eGenes remain unfittable against the pinned panel** — 14 exhaust 5,000 iterations and 45 trigger susieR's prior-variance check, which indicates inconsistency between the summary statistics and the LD reference and was deliberately left enabled rather than overridden.

## Joint testing universe and evidence classes

The joint testing universe comprised nonempty gene symbols with a canonical bulk differential-expression test under this threshold and at least one Tier-1/2 liver-eQTL colocalization test. Genes absent from either testing universe were labeled `not_jointly_testable`; they were not assigned to the neither class.

Primary evidence classes were defined as follows:

- `genetic_only`: SuSiE `PP.H4 > 0.5` and not a canonical DEG.
- `disease_state_only`: canonical DEG and no SuSiE `PP.H4 > 0.5`.
- `convergent`: both SuSiE `PP.H4 > 0.5` and canonical DEG.
- `neither`: jointly testable and meeting neither threshold.
- `not_jointly_testable`: missing from at least one testing universe.

Sensitivity classes substituted the SuSiE-or-ABF union for SuSiE. `genetic_trait_scope` was assigned as direct disease, enzyme, both, or none based on the primary SuSiE result. `genetic_confidence` was assigned as SuSiE, ABF-only, or none.

## Orthogonality audit

Overlap was evaluated for SuSiE and SuSiE-or-ABF definitions at PP.H4 thresholds of 0.5, 0.7, and 0.9, separately for all primary traits, direct disease/PDFF, and liver enzymes. For every analysis we reported the total genetic count, jointly tested genetic count, unique-symbol DEG count, overlap, fraction of jointly tested genetic genes without differential expression, reciprocal fraction of DEGs without genetic support, Jaccard index, Fisher odds ratio, and Fisher p-value.

Fisher tests used the joint testing universe rather than the annotated genome. The descriptive non-overlap fraction was treated as primary. Continuous sensitivity analyses used Spearman correlation between per-gene maximum PP.H4 and either the absolute canonical bulk t statistic or the negative log10 interval-null p-value.

The central non-overlap criterion required at least 85% of jointly tested genetic genes to remain non-DE under primary SuSiE, the SuSiE-or-ABF union, and the direct-disease/PDFF union at PP.H4 greater than 0.5.

## Assay-native evidence-class context tests

Evidence classes were defined only from the canonical bulk differential-expression and primary SuSiE results. Proteomic, spatial, and single-cell endpoints therefore did not contribute to class assignment.

Proteomic positivity was defined as best available protein-level adjusted p-value below 0.05 among genes measured in the integrated proteomic datasets. Spatial positivity was defined by the atlas's measured spatially variable-gene flag among genes with a spatial test. Single-cell positivity was defined as best pseudobulk adjusted p-value below 0.05 among genes tested in the integrated single-cell analysis.

Each modality used its measured-gene background. Raw positive rates were reported by evidence class. Pairwise Fisher tests compared convergent genes with genetic-only and disease-state-only genes. Adjusted logistic models compared every class with the neither reference while including standardized bulk average expression and protein-coding status. P-values were Benjamini-Hochberg corrected across the prespecified class comparisons.

A prespecified historical stress test asked whether the convergent class exceeded both single-map classes in at least two assay contexts. That criterion did not pass and is retained as a supplementary boundary; class superiority is neither a Resource completion criterion nor a reason to privilege convergence.

## Single-cell analyses

Seven public human liver single-cell or single-nucleus datasets were harmonized
at donor and cell-type levels, comprising 104 descriptive biological donors
across 275 integrated libraries and 1,232,318 cells. The analyzed atlas contains
102 donors, 273 libraries, and 1,232,285 cells after two prespecified low-cell
GSE189600 libraries are excluded. Dataset identity was retained during
integration.

The atlas used for Hotspot discovery and primary program scoring did not
receive ambient-RNA correction. A QC harmonization field that recorded decontX
as applied was inconsistent with the execution audit and producing script and
is not used as methods provenance. In the aligned candidate sensitivity,
decontX was run against raw counts without changing the 117 prespecified
program memberships, weights, or labels. It succeeded in GSE136103, GSE174748,
GSE185477, GSE202379, and GSE244832. GSE189600 and Liver Atlas could not be
corrected because an ambient batch contained only one contamination cluster;
both are excluded from the 62-donor, five-dataset primary common-universe
comparison and retained only as explicit uncorrected passthrough datasets in a
full-universe sensitivity. Corrected program scores were collapsed to donor and
tested with model-based and HC3 standard errors. Benjamini-Hochberg correction
was applied separately over all 117 programs for raw, corrected, and paired
correction-induced effects. Integer-count transport was valid for 106 programs;
the 11 prespecified T-cell programs were labeled untestable because their native
non-integer scoring substrate could not be reproduced. This sensitivity can
qualify, but cannot establish, the source of ambient-sensitive transcripts.

Sequencing runs and CD45 sort-fraction libraries were collapsed to biological donor before any per-donor statistic, because several datasets deposit multiple libraries per patient that are not independent biological replicates (GSE244832, 117 sequencing runs from 18 donors; GSE202379, 67 from 46; GSE185477, 21 from 3; GSE136103, 20 CD45-sorted libraries from 10 patients; the remaining datasets were already one library per donor). Per-cell-type counts were summed and cells pooled to the donor before testing. Cell-type differential expression was tested by pseudobulk methods at the donor level rather than treating cells, runs, or sort fractions as independent replicates. Where a disease-stage stratum fell below the pre-registered minimum donor count, or was fully confounded with a single dataset, the contrast was reported as non-estimable rather than significant. The atlas records the best cell-type effect, adjusted p-value, number of significant cell types, hepatocyte effect, and cell-type-specificity annotations.

Cell-cell communication was inferred with LIANA-based consensus analyses and modality-appropriate follow-up methods, aggregated to biological donor. Communication results were interpreted as inferred ligand-receptor coordination, not direct signaling measurements. At true-donor resolution the stage ligand-receptor communication landscape was reproducible in aggregate, but individual ligand-receptor pairs were not claimed to have stage-ordered statistical support.

Same-atlas localization was evaluated for the two primary-selected programs
using donor-paired, lineage-minus-hepatocyte score contrasts. Five lineages were
prespecified for each program, producing a ten-contrast Benjamini-Hochberg
family. The focal comparison is the fixed lineage score paired within biological
donor to that donor's hepatocyte score; model-based and HC3 standard errors are
reported. This analysis describes cellular context within the same atlas and
does not test transcript origin, secretion, or cell autonomy.

For bulk tissue-state transport, the prespecified L1-weighted program members were
projected into F1–F4-versus-F0 bulk fibrosis-stage models. Each displayed value
is a fixed-member weighted gene-level log2 fold change. Point area reports the
fraction of available member weight passing the complete 23,370-gene BH family
within that stage contrast. These summaries are cross-sectional and descriptive;
they are not a program-level significance test.

For the supplementary program-linked communication screen, donor-level scores
for the two primary-selected Hotspot programs were crossed with the fixed
13-pair communication roster. For each donor and pair, the communication score
was −log10 of the LIANA magnitude rank, lower-bounded at 10−4. Each of the 26
planned models used the communication score as the response and the standardized
program score as the focal covariate, with disease stage, dataset, log10 sender
cell count, and log10 receiver cell count included as covariates. Models with
fewer than 30 matched biological donors were labeled untestable. Model-based
and HC3 p-values were corrected separately by Benjamini–Hochberg over all 26
planned tests; support required q<0.05 in both families. Partial residuals were
used only to display the adjusted supported coefficient. Leave-one-dataset-out
directions were descriptive stability checks and did not define another
significance family. Because program and communication scores came from the
same atlas, their association was interpreted as contextual covariation rather
than signaling mechanism.

Bulk transcription-factor activity was inferred separately for each F1–F4
versus F0 contrast using `decoupleR::run_wmean`, DoRothEA confidence A/B/C
regulons, seed 42, and 1,000 permutations. Benjamini-Hochberg correction was
applied across every tested TF within each contrast before HNF4A, RELA, SMAD3,
and TCF7L2 were selected for display. Filled points require q<0.05 in that
complete contrast-specific family. The score represents inferred regulon
activity rather than TF expression or direct binding.

## Spatial transcriptomic analyses

Spatial organization, zonation, and cross-cohort annotations were integrated with each source's measured-gene background, biological unit, and provenance retained. None of the accepted spatial branches supplies a valid MASLD-versus-healthy disease-effect estimand: GSE192741 and Vu test residual spatial autocorrelation against matched genes, Yakubovsky is a non-MASLD tissue-context analysis, and GSE287826 was skipped because no authoritative donor key was public. Spatial results therefore localize or physically contextualize fixed genes and prespecified programs; they are not described as independent disease-associated change.

## Corrected gsMap analysis

gsMap was run on two human MASLD spatial cohorts using 18 liver-disease, PDFF, and liver-enzyme GWAS traits. The primary gene-level endpoint was gsMap's native PCC statistic, the correlation across spots between gene spatial specificity and the trait's spatial heritability signal. PCC was averaged across samples within each cohort and trait.

The previous risk-gene analysis selected genes with high marker enrichment in significant spots and forced a top-decile spot fallback when fewer than ten spots passed p below 0.05. That derived binary endpoint was excluded from the primary manuscript analysis.

The corrected analysis attached the fixed SuSiE and differential-expression classes to the PCC table. Convergent genes were matched separately to genetic-only and disease-state-only genes on protein-coding status, bulk-expression quintile, and median spatial gene-specificity quintile. Up to five controls were selected per convergent gene using nearest expression and gene-specificity distance. The primary effect was the mean PCC difference. One-sided empirical p-values were calculated with 10,000 label permutations and corrected across trait, cohort, and comparison tests.

Direct NAFLD, NASH, and PDFF traits were analyzed separately from ALT, AST, and GGT traits. A replicated result required a positive adjusted q below 0.05 in both cohorts. No direct trait met this criterion. Unrelated high-powered specificity-control GWAS were not available in the completed gsMap run, providing an additional reason not to elevate a spatial-exclusivity headline.

## Proteomic analyses

The active proteomic endpoint uses the PXD051911 liver DIA-MS cohort with assay-appropriate normalized abundance and differential testing. PXD052937 plasma was excluded permanently because its per-sample disease labels could not be authenticated. Protein significance, effect direction, measured-gene coverage, cohort overlap, and selection conditioning were retained explicitly. Protein results are described as assay-native tissue context or corroboration, not functional validation.

For the candidate GNMT–MAT1A–CYP2C19 state analysis, protein abundances in 46
MASLD participants were rank transformed and residualized for acquisition
batch, age, BMI, and sex. Pairwise Spearman correlations were computed between
the residuals, with Benjamini-Hochberg correction over the three prespecified
gene pairs. Sensitivities added fibrosis stage, NAS, or the four component
histology scores to the residualization model and reported both correction
within each three-pair family and across all 12 pair-by-adjustment tests.

## External cross-reference with independent liver studies

External liver protein-QTL and single-cell multiomic studies were used only to cross-reference our fixed evidence classes against source-reported prioritized genes, loci, and transcription factors; neither contributed to class definition or the historical multi-evidence score. The studies were treated as cited comparisons, not reanalyzed data sources. We compared their reported calls with our fixed per-gene values, reporting direct-MASLD and enzyme-trait anchors separately. Loci where layers colocalized but effect directions differed were reported as shared loci rather than sign-concordant validations.

For Elison et al., the genome-wide supplementary tables of cell-type quantitative-trait, gene-regulatory-network, and chromatin-loop statistics were embargoed at the time of analysis; we therefore cross-referenced the per-locus and per-transcription-factor calls stated in the preprint main text and figure legends—the featured colocalized loci, the cell-type-annotated risk loci, and the stage-specific hepatocyte gene-regulatory-network transcription factors—rather than reanalyzing primary data. No primary data from either study were reanalyzed or redistributed, and both are cited as non-peer-reviewed preprints. Each comparison is recorded as external annotation of the fixed evidence-class table and is structured to extend to the full supplementary statistics on their release.

## MASLD Gene Catalog, source dependence, and follow-up experiments

Gene Catalog entries were built only from checksum-bound output files accepted by their
upstream workstreams. Each gene or program retains a stable identifier, evidence
class, testing universe, assay-specific call, biological unit, source dataset,
source-publication relationship, and parent-file graph. Membership in a
program is not expanded into gene-level protein, spatial, chromatin, or
functional evidence.

Evidence states distinguish `supported`, `discordant`, `tested_negative`,
`indeterminate`, `untestable`, `not_applicable`, and `source_dependent` outcomes.
`tested_negative` is permitted only when a prespecified adequate-negative or
equivalence rule passed; non-significance, inadequate coverage, a failed source
check, or an assay that cannot observe the proposed mechanism remains
indeterminate or untestable.

Next-experiment routes are deterministic functions of the fixed evidence,
molecular object, and testability fields. A noncoding regulatory element routes
to allele-aware reporter or locus editing in the lineage where the element is
accessible. A disease-associated lncRNA without a genetic anchor routes to RNA
depletion with neighboring-gene and phenotype readouts. A colocalized lncRNA or
ambiguous lncRNA locus routes to paired RNA depletion and locus-level CRISPRi or
allele perturbation so that RNA-product, transcriptional, and DNA-element effects
can be distinguished. Protein-coding nominations retain the existing
activation, inhibition, protein-level, lineage-restricted, and
context-specific routes. The routing layer neither selects new scientific
results nor uses downstream outcomes. The signed Catalog candidate remains
noncanonical until synchronized Resource adoption.

## Historical multi-evidence score

The convergence output summarizes six source-level evidence channels and exposes the active-source vector, signed concordance fields, evidence-weighted heuristic score, rank, and tier. The score is not a Bayesian posterior and has no gene-level significant hits at permutation-based BH FDR below 0.05 in the fixed build. Evidence classes and source-level vectors are primary; the aggregate score is a navigational field.

Re-running the historical score after collapsing sequencing runs and CD45 sort-fraction libraries to biological donor left its Tier-1 set byte-identical and its rank correlation near one. This robustness does not validate the score as a target hierarchy; it remains a supplementary navigational field, and the Gene Catalog defaults to gene class and source-level calls rather than rank.

## Supplementary strict drug-target benchmark

MASLD-specific clinical labels were rebuilt without using broad Open Targets preclinical annotations as development stages.

- `approved`: manually curated MASLD/MASH approval.
- `clinical`: manually curated active program or MASLD-specific ClinicalTrials.gov record with phase at least 1.
- `discontinued`: manually curated discontinued MASLD program.
- `association_only`: Open Targets animal-model or literature association without qualifying clinical evidence.
- `none`: no qualifying MASLD-specific evidence.

Open Targets global drug phase and approval fields were not treated as MASLD-specific status. Strict preclinical intervention labels were unavailable and no ordered discovery-to-approval gradient was tested.

Clinical/approved targets were matched to druggable controls on protein-coding status, liver-expression quintile, and grouped Pharos target-development level. Up to five unique nearest controls were selected per positive target. Score performance was summarized by AUROC. Paired stratified bootstrap resampling with 5,000 iterations compared convergence with active-modality count, genetics-only PP.H4, and bulk evidence on the identical matched universe.

The prespecified headline criterion required at least two approved targets in the ranking universe and convergence superiority of at least 0.03 AUROC over modality count and the best constituent channel, with a paired 95% confidence interval excluding zero. That criterion did not pass.

## Fine-mapped DNA architecture and target-transcript biotype

Fine-mapped DNA and colocalized target transcripts were analyzed as separate
objects. The fine-mapping unit was one reliable, converged SuSiE credible set,
identified by study, trait, locus, ancestry, and credible-set ID. Sets with raw
PIP sums outside 0.90 to 1.10 failed the primary criterion. Within eligible sets,
SuSiE PIP was normalized and divided into mutually exclusive
protein-altering, canonical-splice, synonymous-or-UTR, other-noncoding, and
unresolved masses. Promoter proximity, liver ABC-enhancer overlap, and human
liver lineage-accessible chromatin were nonexclusive context annotations and
were not summed as probabilities.

The target-transcript analysis used corrected gene-by-GWAS COLOC rows and
classified positive targets by versioned GENCODE biotype. ABF and SuSiE
denominators, and direct MASLD/PDFF and enzyme-trait results, remained separate.
A credible set was not assigned to a colocalized transcript through a lead- or
top-SNP shortcut. Unless the corrected output explicitly identified a shared
signal pair, the credible-set-to-gene link was reported as unresolved. Variant
location within or near a lncRNA locus was not interpreted as evidence that the
lncRNA molecule mediated the association.

## Multiple testing and statistical reporting

All tests report the analysis universe, numerator, denominator, effect size, and uncertainty where available. Benjamini-Hochberg adjustment was applied within prespecified analysis families. Permutation p-values used the plus-one correction. Null or failed acceptance criteria were retained in the release and used to constrain manuscript wording.

## Data availability

This study integrates publicly available datasets, cited with their accessions in the relevant Methods subsections, source registry, and reference list. The final Resource release will include source checksums, access dates, testing-universe manifests, a data dictionary, and the source-dependency graph used by the Gene Catalog.

## Software and code availability

The historical analysis components are generated by the following files. Their
`passport` filenames are retained only because they are hash-bound to the
recorded v1 candidate; the manuscript and portal call the product the MASLD Gene
Catalog.

- `scripts/manuscript/build_evidence_class_release.R`
- `scripts/manuscript/build_strict_drug_benchmark.py`
- `Analysis/Spatial/scripts/15l_gsmap_pcc_narrowing.py`
- `scripts/manuscript/build_gsmap_class_validation.py`
- `scripts/manuscript/build_coding_architecture_release.R`
- `scripts/portal/generate_evidence_passports.py`
- `scripts/portal/build_evidence_passport_bundle.py`
- `scripts/portal/validate_evidence_passports.py`

Plan 60 alone snapshots the selected inputs, figure sources, Gene Catalog, and
manuscript ledgers into the final read-only Resource release. R analyses use the
repository's `rnaseq` environment. Python spatial, benchmark, and Gene Catalog
analyses use their recorded existing environments. Randomized analyses record
their seed and complete permutation family in the release manifest.
