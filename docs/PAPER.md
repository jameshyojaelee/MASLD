# Paper identity, novelty, and claims

**Abstract numbers and figure logic, 2026-11-12 (user-approved rewrite).**
These results are used in the manuscript text by user approval, 2026-11-12;
formal adoption under the ROADMAP final release checklist is pending. Each
number below is from the novelty section of
[RESULTS](RESULTS.md#novelty-analyses-2026-1110-11-not-adopted) or the named
`RESULT_*.md` file. Genetic shares are covariance partitions, not mechanisms.

Abstract numbers:

- ALT heritability shared with neither liver fat nor cirrhosis: 0.70 (95%
  interval 0.60 to 0.77) (ENZFAT1).
- 260 of 286 direction calls have liver-enzyme support only; 243 of the 260 are
  underpowered for liver fat, so this is mostly absence of evidence (DIR6).
- Hepatocyte-protein panel in plasma: 0.434 protein SD per ALT SD across 169
  ALT loci (PROT5).
- TNFSF10: ALT colocalizes with plasma TRAIL (PP.H4 1.00) (TRAIL4). The trial
  table has no TRAIL-pathway trial; the trial result is the JOIN7 null (no
  matched pair informative, Holm p 1.0).
- Plasma lipid species colocalize at 4 of 6 liver-fat ALT loci against 10 of
  143 others (OR 25.2, p 6.8e-4) (MQTL2). TM6SF2, GPAM and PPP1R3B are
  established lipid loci (source-owned).
- One AlphaGenome 2 number: liver accessibility benchmark, Spearman 0.719
  against caQTL (version 1 0.736; hepatocyte single-nucleus ATAC 0.747)
  (AG2X Stage A).
- Supplement only, one sentence each: the expression-direction null (sign
  agreement 0.526, CI 0.492 to 0.561) and the GSET8 cirrhosis gene-set result.

Figure logic, one sentence each:

1. Figure 1 is the user's plate and is unchanged.
2. Figure 2: fine-mapping and liver-eQTL colocalization narrow liver-trait GWAS
   signals to candidate genes for testing in other measurements (28,233 Catalog
   genes, 18,720 testable, 429 colocalized, 419 in bulk RNA); panels qualify each
   candidate (PP.H4, lead consequence, ancestry, trait source) and show which
   measurements can test it next (bulk RNA, single-cell RNA, liver and plasma
   proteome, CosMx). Colocalization marks a gene for follow-up, not a disease
   gene. Enzyme-locus MASLD liability, lipid and hepatocyte-protein results,
   ALT heritability partition and TNFSF10 are supplementary (Figure S2Z to
   S2AB, Supplementary Note 22), decided by the user 2026-11-12.
   Sequence-supported candidate variants appear only as an
   accessibility annotation beside measured evidence: 126 of 4,639
   credible-set members, of which about 101 are expected from the rule alone.
   No predicted gene, no predicted expression direction, never "causal
   variant".
3. Figure 3: of 7,214 replicated stage genes, 1,155 change inside hepatocytes.
4. Figure 4: most CCL20 RNA is in cells with hepatocyte expression profiles
   (median 0.65 to 0.83 of CCL20 UMIs after removing cholangiocyte-like cells,
   three cohorts; CCL1). This replaces the earlier label-based statement that
   hepatocytes are the main source. The prespecified primary identity
   criterion was not met (one primary cohort's control failed; the other gave
   Holm 0.078).
5. Figure 5: Complex I protein is lower without an RNA fall, and PLIN2 protein
   exceeds its RNA. Planned sentence only; its values are in
   `Analysis/Hypotheses/DIG1-20261024T113951Z/RESULT_DIG1.md` and not yet in
   the manuscript text.
6. Figure 6: Catalog genes by disease evidence, direction, hepatocyte
   expression, compounds and trial results. Figure 6G is a per-mechanism table
   (28 mechanisms, alphabetical): target-gene tiles by candidate Catalog
   genetic state, highest randomized trial phase (not approval) and four
   endpoint columns (liver fat, ALT, MASH resolution, fibrosis improvement)
   (`figures/candidates/fig6g-trial-mechanisms-20261112T140124Z/`). The JOIN7
   gene-level support is not plotted; no matched genetics-by-endpoint pair was
   informative (Holm p 1.0). No universal score or target rank.

**Organizing claim and figure plan, 2026-11-09 (user-approved review plan).**
The paper keeps its two arms but is organized around two positive partitions
instead of the gene-set overlap:

1. Genetic arm: of 167 colocalized ALT/AST loci (after the eQTL refit; 169 in r2), 18 carry MASLD liability at
   the same lead SNP and 46 raise ALT without raising MASLD risk (enzyme only;
   102 intermediate, 1 opposite), against the PNPLA3/TM6SF2 benchmark ratio
   0.891. Out of sample (FinnGen MASLD, PDFF) enzyme-only and intermediate loci
   are not separated; both carry about a quarter or less of the benchmark
   signal. This is exploratory (DIG1 job 73, post hoc).
2. State arm: of the 7,214 discovery F3–F4 genes replicated in SteatoSITE
   (493 donors), 1,155 change inside hepatocytes after ambient correction in
   four snRNA cohorts (689 up, 466 down); most of the rest cannot be separated
   from cell mix in bulk data (CONF1 C3/C3b, post hoc combination).
   GNMT, MAT1A and CYP2C19 are lower per hepatocyte in the held-out GSE289173
   cohort (CONF1 C5; controls are a separate series).

The genetics–state overlap (26 vs 21.4 expected) becomes one sentence and a
supplementary panel. Figure 7 moves to the supplement; the main paper has six
figures. The genetic class rule admits a fine-mapped coding variant (PIP ≥ 0.5)
as a second mechanism type beside colocalization and adds a disease-liability
field; this moves 20 genes, including PNPLA3 and TM6SF2, from unresolved to
anchored (CONF1 C4). The rule is approved; the rebuilt Catalog table is a
candidate until the final checklist. Exact values are in
[RESULTS](RESULTS.md#review-round-confirmatory-tests-conf1-2026-11-09-not-adopted).
Follow-up screen outcomes are not used.

**Variant identity and transportability, 2026-11-07.** The 436-gene map is a
provisional historical fit set: 90 genes use eQTL fits with duplicate variant
identifiers, and MPRIP has both an ambiguous repeat event and invalid
regularized reference LD. The remaining 346 genes have not thereby acquired
verified GWAS variant identity or source-matched LD. Genetic counts, Catalog
classes and downstream overlap/accessibility-QTL claims retain this qualification
until corrected full fits are evaluated together.

The 1,854 pooled tissue calls remain discovery associations. Their directions
are broadly concordant (1,850 same-sign random-effects means; 1,208 same-sign
effects in every available cohort), but no gene passes full-family BH in the
24,115-gene REML/modified Hartung–Knapp sensitivity (minimum q=0.171430).
This does not establish absence of tissue effects. Historical consistency
classes and within-discovery agreement are not independent replication.
The original prespecified SteatoSITE stage comparison retains its separate
contrast and source credit. Genetic cross-study overlap is gene-level recurrence
with a shared eQTL panel, not evidence that the same GWAS signal replicated.

**Physical-capture correction, 2026-11-07.** Repeated sequencing runs from one
capture do not contribute independent cells. Deeper sequencing may detect
additional molecules; joint counting must deduplicate shared molecules and
retain unique capture/barcode identities. Antibody-tag reads are excluded from
RNA counting even when deposited under the same experiment accession.
The 117 program memberships and weights stay fixed. A fixed-weight expression
projection and the historical neighborhood-smoothed Hotspot PCA score are
different measurements and must not be relabeled as each other. Historical
embeddings and discovery stability remain conditional on their original input.
Histology-error calculations require explicit assumptions and cannot establish
empirical reliability or a performance ceiling from serial-biopsy variability.

**Interpretation corrections, 2026-11-06.** The checked lineage-specific gene-DE
fit uses at least 50 contributing cells per donor and tested lineage; the
117-program stage result remains 0 selected. Intronic/exonic concordance supports
a transcription-associated tissue component but cannot exclude RNA processing,
stability or composition effects. Protein abundance does not measure metabolic
function, between-batch agreement does not establish an assay ceiling, and a
minimum detectable effect is not a bound on the true effect. Xenium disease
resampling repeats carrier selection and retains its age/source qualifications.
Near-one colocalization posteriors with LD diagnostics remain qualified rather
than treated as certain mechanisms. Current correction values and the explicit
historical-versus-candidate Catalog distinction are in
[RESULTS](RESULTS.md#review-corrections-2026-11-06). These changes do not redefine
the 117 programs, the 43 systems or the two-paper boundary.

**Single-cell reference definition, 2026-10-26.** Disease comparisons use
source-reported diagnoses with histology and etiology retained separately.
Conflicting diagnoses are excluded rather than assigned an invented disease
stage. A source-labelled nonsteatotic control is not necessarily numerically
confirmed S0A0F0 tissue; P30 requires its activity qualification and exclusion
sensitivity. Numeric fibrosis analyses require deposited numeric SAF and must
not translate free-text healthy-control/end-stage categories into measured
F0/F4. Liver Atlas's histology-only records and GSE136103 Alcohol/PBC cirrhosis
cannot become MASLD cases by harmonization. The corrected RNA comparisons
retain donor-level inference and the fixed 117-program coordinate system;
they do not verify RNA–ATAC pair identities or physical sample swaps.
The completed correction is described in
`Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/REPORT.txt`;
the user requested incorporation into main Figure 4 on 2026-11-03. Its
individual panels, compact assembly, legend and Results text now use these
corrections. Source-dependent supplementary artwork and other outputs still
await incorporation; named release r9 is unchanged.

The November 2 exploratory checks also qualify Figure 4. The focal gene sets
were not recovered by Hotspot without GSE244832; their donor-level
co-expression outside that dataset does not establish hepatocyte origin.
Their continuum associations do not exceed matched random-axis references
(Holm over four program-by-axis tests). This comparison did not test the 43
systems. Neither program specificity nor independent molecular mechanisms
follow from the displayed shared expression gradient.

**Figure 3F interpretation, updated 2026-10-14:** the installed 54-gene matrix illustrates
cross-sectional expression and exact-pair genetic support as separate quantities.
It is conditioned on stage significance in the corrected 774-unit fits;
no enrichment, mediation, or independent-donor replication follows. The panel
shows expression and PP.H4; the PIP column was removed at the user's request.
The 216 effects use four complete 24,115-gene BH families and 205/190/157/94
selected genetic-identity units. These are not verified donors. The historical
844-record matrix had 61 genes and 244 effects and is superseded. The current
source is `figures/candidates/fig3f-no-pip-774-20261012T011128Z/`, with verified
COLOC release `susie_coloc_r2_20260924T053244Z`, selected October 14.

**Historical Figures 3–4 identity qualification, 2026-09-29:** the genotype ledger identifies
repeated bulk libraries and conflicting source labels; GSE193066 identity calls
remain unresolved. Historical bulk roster counts and newly reconstructed fits
are conditional, not verified independent-donor inference. The correction
candidate includes an explicit roster join and a separate provisional identity
sensitivity. That restriction applied to the September candidates; the user
adopted the corrected-count bulk refit on genetic-identity clusters on October
11. Neither that decision nor the historical network analysis verifies
independent donors or qualifies GSE193066; see
[the disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).

**Historical Figures 3–4 correction boundary (2026-09-29).** The defensible claim is
reproducible cross-sectional molecular remodeling and tissue covariance.
Use source histology categories, not invented exact stages; the two continuum
scores are closely related methods. Separate discovery lineage, content lineage
and the lineage tested for disease effects. Ordinary 117-program selection and
HC3 sensitivity remain distinct. Metabolic abundance is not measured function;
unsigned PP.H4–signed logFC comparisons do not establish biological opposition
or complementarity. Figure 4H plots 43 system coefficients estimated in the 367
disease units and carries no outline. On the 435-unit roster that includes the
68 GSE213621 controls, one system exceeds the uniform size-matched null and
none the connected-region null (0 connected, 1 uniform), and 18 pass strict
maxT; none of these is marked in Figure 4H. No claim of
exceptional vocabulary specificity or cell-to-cell regulation follows.
[Supporting evidence and dispositions](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md)
include the corrected candidate analyses; final adoption remains separate.

**Historical correction implementation, 2026-10-05; N4 updated October 14.** Primary candidates exclude unresolved
donor/source-category conflicts and require complete counting-annotation
correspondence. Existing bulk numbers and the 18 maxT outlines remain conditional
on their historical library rosters and counts; neither is a target for the
corrected analysis. Preserve the 117 programs, 43 system memberships, external
continuum definitions and k4/k6 bases. The rechecked identity screen identifies
conservative genetic clusters, not verified independent people; GSE193066
remains excluded from primary replication, ICC and derived ceiling claims.
The N4 eligibility correction completed on 2026-10-14. Its 22 different-visit
pairs remain a retrospective description on historical counts, not qualified
primary within-person evidence. Equal N1 eligibility and exact-stage matching
leave insufficient independent null units; the historical 87 random-null calls
are withdrawn. The historical 69 stage-BH calls also came from uncalibrated
bootstrap tail fractions and are withdrawn. The untyped 28-pair comparison
has zero calls against its eligible exact-stage null (BH over all 117 fixed
programs). See [the corrected N4 results](RESULTS.md#n4-within-person-agreement-eligibility-correction-2026-10-14).
Both six-panel Figure 7 sets are prepared separately for review. See the
[implementation record](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#resource-corrections-implementation-2026-10-05).

**Status:** scientific scope authority  
**Updated:** 2026-10-14  
**Paper type:** standalone human MASLD Resource and analytical framework

## Title

**A multimodal human MASLD Resource records liver-trait genetic evidence and
disease-state liver expression side by side for each gene**

## Ultimate goal

The goal is to publish a standalone human MASLD Resource that maps five linked
dimensions of human disease biology: inherited genetic variation relevant to
MASLD and liver traits; multicellular remodeling in established disease; how
genes and gene sets relate to the three disease aspects of metabolism,
inflammation, and fibrosis across cross-sectional disease stages; and how
inherited variants, combinations of variants, and chromatin state relate to
those aspects. The Resource connects those maps to donor-aware cell programs
and protein, chromatin, and spatial tissue context, then organizes the results
in a transparent MASLD Gene Catalog with appropriate follow-up experiments. It
distinguishes protein-altering variation, supported regulatory relationships,
unresolved noncoding signals, lncRNA transcription, and the mature RNA
molecule rather than assuming that all inherited variation is regulatory. It
records where an aspect-specific or variant-combination question is supported,
indeterminate, or untestable with the available assays rather than leaving
that distinction implicit. Finally, it measures what current predictive
models, including pretrained sequence and single-cell models, can and cannot
learn from those maps, reporting every score beside an audit of what each
model already saw during training so that apparent skill is separated from
prior exposure. The paper must be scientifically useful without
follow-up screen results; the future experimental screen prospectively tests Resource
nominations but cannot redefine its classes, figures, or claims.

## The central scientific contribution

The user-approved continuum extension distinguishes an intervention hypothesis
from a measurable disease-state readout. Its candidate implementation preserves
signal-specific genetic uncertainty, continuous variation within fibrosis
stage, and assay-specific cell/tissue context. Target and readout need not be
the same molecule, but shared program membership alone cannot establish their
regulatory connection. Participant-linked liver/plasma associations in an
existing public cohort support a measurement relationship; they do not imply
bulk-RNA coupling, liver origin, treatment benefit, or a new cohort discovery.
The intended use is mechanistic experiment and trial-design hypotheses within
existing histologic eligibility groups. Neither a novelty-first claim nor
"largest atlas" is established by this analysis. Its implementation and limits
are recorded in the [technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md).

The paper records **liver-trait genetic evidence** and **disease-state
remodeling** in human MASLD side by side for each gene. The two maps describe different measurements: inherited variation
identifies trait-associated genetic signals and candidate molecular links,
whereas multicohort transcriptomics identifies the active multicellular state
of established disease. Donor-aware single-cell programs and assay-native
protein, chromatin, and spatial measurements show where those signals and
remodeling programs are expressed, accessible, or physically organized. The
MASLD Gene Catalog converts the combined map into a reusable guide to gene
interpretation and experimental design.

The paper establishes eight linked ideas:

1. Trait-associated inherited genetic variation and the established MASLD
   transcriptome describe inherited trait associations and cross-sectional
   disease remodeling. Their overlap (26 genes against 21.4
   expected within expression deciles, among 3,589 eligible genes evaluated by
   multi-signal COLOC; expression not fully controlled, so no direction is
   stated) does not establish biological complementarity or antagonism. Most
   state-only genes (935 of 1,129 in the 14,704-gene interface) have no
   multi-signal COLOC result, so their genetic call is indeterminate, not
   negative.
2. Multicohort bulk analyses describe established disease as cross-sectional
   remodeling, and donor-aware single-cell analyses place fixed programs in
   cell lineages; no single-cell program passes the 117-program family test for
   a stage association (H59). Recorded histology
   remains the clinical severity anchor. Following Kamzolas et al. (2026;
   doi:10.1038/s42255-026-01543-7), an externally anchored molecular continuum
   provides a complementary layer for describing
   transcriptional variation among participants assigned the same fibrosis
   stage; it is not a universal or longitudinal disease axis.
3. Protein, chromatin, and spatial assays place genetic candidates and disease
   programs in assay-native molecular and physical tissue contexts.

   The October 25 ATAC completion retains the current voom reference and treats
   proximity/coaccessibility links as descriptive context. Internal prediction
   robustness, cross-cohort accessibility agreement, and independent MASLD
   validation are separate claims. Separate-library RNA may corroborate disease
   effects after unique-cell reconstruction using source-native RNA diagnoses;
   the completed RNA arm is testable only in hepatocytes and remains sensitive
   to recorded batch. Donor-matched associations are untestable without a
   documented cross-assay specimen key. The Cicero gain over distance also
   persists after within-donor RNA shuffling, so correct nucleus-level pairing
   has not demonstrated added predictive value for Cicero beyond proximity.
   Manuscript inclusion was approved on October 26. The main-text finding is
   modest cross-cohort accessibility agreement with batch sensitivity, not a
   replicated MASH regulatory program. Four stellate peakwise candidates,
   descriptive links, internal prediction controls and CHST9/TMEM150B
   single-cohort RNA associations belong in the supplementary results.
   These analyses cannot redefine the 117 programs or assign causal genes.
4. Coding variants, noncoding DNA, supported regulatory relationships,
   transcription through a lncRNA locus, and the lncRNA molecule are distinct
   candidate objects whose interpretation cannot be inferred from proximity.
5. The MASLD Gene Catalog connects each nominated object to its biological
   context, evidentiary boundary, unresolved alternative, and next experiment.
6. MASLD is not one axis. Metabolic, inflammatory, and fibrotic change are
   three distinct disease aspects, and a gene or gene set may relate to one
   without relating to the others. The Resource states which aspect an
   observation speaks to, and records `indeterminate` where the available
   design resolves only a summary severity score.
7. Inherited risk is not always one nucleotide acting alone. Variant
   combinations at a locus, and the chromatin state a variant sits in, are
   distinct candidate objects from a single colocalized variant. The Resource
   keeps them separate and records which are supported, which remain
   unresolved, and which the available assays cannot test.
8. What a predictive model can learn from these maps depends on the task. Of six
   pretrained sequence models scored in seven rows, one (HyenaDNA) beats a
   sequence-free allele-identity ridge on a HepG2 reporter assay after
   correction, but not a prespecified flank and dinucleotide ridge (gain 0.049,
   Bonferroni lower bound -0.015); adding either pretrained
   embedding to distance lowers the score on endogenous eQTL effects, and no
   sequence model was scored on tissue-level severity. No score is
   interpretable without knowing
   what each model already saw during its own training, so every score is
   reported beside that record.

The practical consequence is an integrated field resource rather than a
universal target rank. Readers can move from a germline signal or a
disease-remodeling program to its cell and tissue context, determine what the
available assays support, and identify the experiment needed to test the
remaining biological explanation.

## The three disease aspects and the variant-combination dimension

Two questions the Resource is built to serve are stated here explicitly because
they are easy to lose behind a single severity score or a single lead variant.
Both are recorded with an evidence state from the vocabulary above, and neither
carries a claim today.

### Aspect resolution: metabolism, inflammation, fibrosis

**Question.** Which disease aspect does a gene, program, or gene set speak to,
and does that differ across cross-sectional disease stages?

**What exists.** Recorded fibrosis stage is the severity anchor throughout. NAS
enters as a graded sum in adjacent contrasts and as an adjustment set: the
GNMT-MAT1A-CYP2C19 proteomics audit holds after adjusting for fibrosis stage but
loses support after adding NAS or all four histology components, which is a
recorded result about the components rather than about the three aspects.

**What is missing.** The three aspects have never been resolved as separate
outcome axes. A NAS sum cannot attribute an observation to steatosis,
ballooning, or lobular inflammation, and no aspect-specific claim is available.

**Evidence state: separability `indeterminate` in one cohort; attribution
`indeterminate`.** The separability of the three NASH-CRN components was
examined in GSE267145, in 99 participants carrying steatosis, ballooning and
lobular inflammation as separately recorded columns, against an RNA
prespecification recorded before the analysis. Two of its three criteria are
properties of the histology grades alone. The grades correlate at Spearman
0.68, 0.73 and 0.82 (the rule needs every pair below 0.80; ballooning and
lobular inflammation give the 0.82) and span a participation ratio of 1.42 of a possible 3.0,
below the 2.0 floor, so the rule could not return a positive result whatever
the RNA showed. The RNA criterion also failed: the per-gene association vectors
of the three grades correlate at 0.83, 0.87 and 0.87 against a 0.80 threshold.
In simulations on these 99 grades with independent genes, in which 2% to 40%
of genes carry signal at an effect of 0.5 or 1 noise standard deviation, one
shared program
and three programs that each follow one recorded grade give the same values
(0.82 to 0.83 when a fifth of genes carry signal). Programs that follow only the
part of each grade not shared with the other two give values below 0.66, and an
equal mixture of one shared program and grade-specific programs gives values
below 0.75; the observed values
match neither. The separate H3K27ac analysis agrees: its per-region association
vectors give a participation ratio of 1.15, against a null centred at 1.46,
inside the simulated range for the first two forms (1.09 to 1.43) and below
that for grade-specific programs (1.48 to 2.08). The test therefore cannot tell
one shared program from separate programs that follow the correlated grades.
Attribution remains `indeterminate`: nothing can be assigned to steatosis,
ballooning, or lobular inflammation individually.

This concerns one cohort. It does not address the SAF partition,
which divides the same tissue differently into steatosis, activity, and
fibrosis. SAF grades in GSE202379 and PXD051911 remain unused substrate, and
GSE135251 carries only the NAS sum. An unresolved result on the NASH-CRN
partition in one bulk cohort says nothing about SAF.

### Variant combinations and chromatin state

**Question.** Do combinations of variants at a locus, or the chromatin state a
variant sits in, relate to MASLD differently from a single variant considered
alone?

**What exists.** Fine-mapping partitions credible-set mass by consequence class
rather than assuming inherited variation is regulatory. One two-variant local
haplotype is characterised at GNMT and reported as unresolved. Chromatin is
evaluated across jointly testable peaks, the large majority of which resolve to
`indeterminate` and a minority to `source_dependent`.

**What is missing.** No systematic treatment of variant combinations beyond that
single haplotype, and no epigenetic effect claim of any kind. Accessibility
supplies assay-native context, not mechanism.

**Evidence state: `indeterminate` for combinations, `indeterminate` for
chromatin effect.** The single-variant colocalization design does not test
combinations, and an accessibility context does not test causation.

### Why both are stated at goal level

An agent or reader who takes only the two original dimensions will build a map of
variants and remodeling and never ask which aspect a result belongs to, or
whether an effect requires more than one variant. Both questions determine how
later functional evidence can be interpreted.

## Disease-state vocabulary

Cross-cohort displays use the harmonized labels **Healthy**, **Steatosis**,
**Steatohepatitis**, and **Cirrhosis**. Source-native `MASL` is mapped to
Steatosis only where that source defines a steatosis group; source-native
`MASH` or `NASH` is mapped to Steatohepatitis when the source diagnosis supports
that interpretation. Generic `MASLD` is an umbrella diagnosis, not a disease
stage, and is never converted to Steatosis. Source-native terms remain in
provenance tables and source-specific contrasts. Fibrosis F0–F4 and NAS are
separate histologic axes and are not relabeled as these disease states.

## Narrative spine and recurring biological examples

The Results should follow biological questions rather than modality order:
what is inherited, what state is remodeled in established disease, where that
state is organized in tissue, and which experiment resolves the remaining
mechanism. Two outcome-blind examples carry that sequence without being
combined into a modality score.

**Working addition, 2026-11-13 (panels installed, not adopted).** Four example
sets are labeled in the same way in Figures 2 to 5 so a reader can follow them
across modalities: BCAA catabolism (PPM1K, BCKDHA, BCKDHB, DBT) and ADH1B as
known biology; CFHR4/CFHR5 and HKDC1 as less-studied candidates. They were
chosen after all modalities had been examined (XMOD1, exploratory), so they are
illustrations, not outcome-blind examples and not tests. Panels: 2B labels, new
2H (genetic evidence matrix), 3F and 3G labels, 4D program labels, 4E
(hepatocyte forest for the eight genes; the GNMT/MAT1A/CYP2C19 forest is now
Figure S4AG) and 5A labels. HKDC1 is a colocalized gene whose lead allele has
the opposite sign for cirrhosis and liver fat; it is not a risk gene. CFHR4 is
not supported in the one hepatocyte cohort that was not used to choose it.
Numbers are in `docs/RESULTS.md` (XMOD2). The body text of Figures 2 to 5 and
the abstract still name the earlier examples.

1. **A GNMT locus-to-state bridge.** In the current official genetics release,
   GNMT has a single-signal shared-signal observation for direct NAFLD, with an
   unresolved local haplotype. At nonpalindromic proxy rs2296805, lower
   expression aligns with higher disease risk; the current anchor rs2296804
   (6:42931261, SNP.PP.H4 0.258990) is palindromic and its reported labels do
   not independently establish strand. Native European r2 with the proxy is
   0.998001. Related variants have source-reported human liver
   protein and reporter evidence. The adopted corrected COLOC release retains the exact direct-NAFLD ABF
   pair (PP.H4 0.705832); the molecular mechanism remains provisional. In the
   disease map, GNMT joins MAT1A and CYP2C19 in lower multicohort RNA,
   lower RNA abundance at F2 to F4 than at F0 in three staged cohorts (the GNMT
   and MAT1A differences at F1 depend on GSE130970), and lower liver protein. The three
   proteins also covary within 46 MASLD participants after adjustment for
   acquisition batch, age, BMI, and sex; attenuation after NAS or component
   histology adjustment supports a coherent severity-linked metabolic state,
   not an autonomous three-gene mechanism. GNMT is the bridge between the two
   maps. MAT1A and CYP2C19 are state partners, not an inherited regulatory unit.
2. **Tissue organization does not determine lineage origin: Stromal ECM and
   Ductular injury programs.** The programs were prespecified in
   hepatocyte-labelled profiles, but the original atlas was not ambient-RNA
   corrected. The aligned 117-program candidate uses current source-diagnosis
   metadata and a common-dataset primary universe. Original-source checking
   identifies two GSE174748 disease explants as cirrhotic, excluding them under
   the unchanged three-stage rule. With source-diagnosis corrections adopted 2026-11-03 (GSE202379 P98
   excluded, P30 retained as a qualified reference, GSE244832 D11/D14 set to
   their deposited RNA diagnoses), a separate residual variance per dataset and
   a within-dataset permutation null on at most 59 ambient-correctable donors
   (54 hepatocyte donors fitted, all slope information from GSE202379 and
   GSE244832), no program passes the 117-program BH/stability rule (smallest
   q=0.18135): Ductular injury β=0.2065 (95% CI 0.0837 to 0.3293), q=0.18135;
   Stromal ECM β=0.1624 (0.0113 to 0.3134), q=0.242929. Ambient-corrected GLS
   selects none of 106 computable programs (smallest q=0.11115). The
   equal-variance fit on at most 61 source records is a supplementary
   sensitivity. The 60/62/64-donor estimates are historical.
   GSE189600 and Liver Atlas
   remain explicit uncorrected passthrough datasets only in the full-universe
   sensitivity. Donor-paired same-atlas contrasts descriptively localize
   Stromal ECM to fibroblasts and Ductular injury to cholangiocytes, without establishing
   transcript origin or cell autonomy. After all-cell adjustment, Stromal ECM
   exceeds its matched-gene spatial null in both Visium sources and Ductular
   injury does not (indeterminate), but 68 and 74 of 114 fixed programs also
   exceed that null, so this is not program-specific organization; no
   prespecified test compared the two programs. The null holds donors fixed, and Stromal ECM excess is positive in
   3 of 4 donors. Disease remodeling, tissue architecture, and lineage
   attribution are reported as separate claims. It does
   not establish hepatocyte origin or an IGFBP7–BICC1 pathway. IGFBP7 is one of
   49 fixed Stromal ECM genes and contributes 2.47% of its L1 weight; BICC1 is
   one of 22 Ductular injury genes and contributes 5.44%. All spatial inference
   here applies to the weighted program scores, not to either gene individually.

The corrected COLOC r2 release was independently checked and adopted on
2026-10-14. It contains 436 primary tier-1/2 SuSiE genes, 560 single-signal-only
genes, a 996-gene union and 708 full-portfolio SuSiE genes. All 50 studies are
represented; the shared-posterior correction reran 1,307 previously positive
gene-study rows across 32 studies and carried the other results unchanged.
In it GNMT has SuSiE PP.H4 0.128
and ABF PP.H4 0.706 (2021_34841290_NAFLD_EUR, where SuSiE was attempted for
GNMT but returned no credible-set result), so GNMT remains an ABF-supported provisional bridge; COLOC alone does not identify a causal gene. Its genetic
membership is adopted; its unresolved haplotype and ABF-only support remain
explicit. The ambient and
same-atlas localization candidates passed independent validation but are not
yet adopted, pending synchronized adjudication. These conditions affect the strength of
the worked examples, not the paper's two-map identity.

## Evidence discipline

The two maps are credible only when phenotype provenance, ancestry,
joint testability, assay observability, biological unit, source independence,
and applicability are explicit. Apparent negatives therefore remain
`untestable`, `source_dependent`, `indeterminate`, or `not_applicable` when the
relevant test was not adequate. A strict `tested_negative` call is reserved for
an adequately observed, applicable, properly joined assay that fails its
prespecified test. These states are the Resource's guardrails, not its central
biological headline.

Observability limits are an important Resource result, and their causes differ.
Panel content limits CosMx, which adequately observes 12 of 117 programs. GeoMx
and GSE287826 each test 0 of 117 programs because no gene axis is loaded for
either: the GeoMx release was dropped after its direction and normalization
could not be resolved, and GSE287826 has no AOI-to-donor key, so its
observability is not yet determined. The single-cell ATAC hepatocyte and
stellate contrasts have 5 normal and 9 MASH donors in GSE244832 and 6 and 6 in
GSE281367; 38,066 of 39,914 jointly testable peaks are `indeterminate`, and
1,844 are `source_dependent` (q < 0.05 in one cohort only; 1,505 of them have
the same sign in the other cohort). Against these, the four Visium/ST datasets
each observe 114 of 117 programs. The current Catalog defines no negative-call
rule, so `tested_negative` is not assessed. For bulk RNA, a
spike-in calibration (H81) measures what such a rule would need: a uniform 1.0
log2FC change is called in at least 0.98 of genes (lower bound) in expression
deciles 3 to 10 but in 0.507 and 0.796 in deciles 1 and 2, and the candidate
equivalence rule meets its error limits only in deciles 7 to 9. These are not failed validations.
They show why modality-vote counting converts absent coverage, unresolved donor
identity, and assay mismatch into false confidence.

The adopted eligibility-corrected Gene Catalog (`gene-catalog-eligibility-20261106`)
retains 28,233 gene identities: 400 genetically anchored, 26 supported in both
domains (class `concordant`; effect directions not compared), 1,118
established-state associated, 12,958 unresolved and 13,731 untested. Its bulk input uses 774 conservative genetic-identity
clusters, not verified independent donors; 855 held annotations remain
untestable. The November 6 user-approved correction applies the existing
priority eligibility rule: 188 genetic calls become untestable, and 180 classes
change (169 unresolved and 11 state-associated become untested). All RNA
effects/calls, 436 native genetic positives and 117 program definitions are
preserved. An untested interface class can
coexist with positive evidence in one measured domain. The 14,720 pre-adjudication
symbol-grouped rows, 14,704 native stable-ID interface rows and 14,502 jointly
eligible Catalog identities are different universes and must be named beside
their estimates. The molecular-object companion infers no credible-set target
link or causal mechanism from a gene-level support call; see
[the current Catalog results](RESULTS.md#catalog-eligibility-adoption-2026-11-06).

For spatial omics, the finding concerns two programs fixed in advance from the
superseded August stage fit; neither passes the corrected 117-program rule.
**Reproducibly positive spatial autocorrelation is not sufficient evidence that
a program has program-specific tissue architecture.** After residualizing for
all 16 cell2location abundances, both programs had positive aggregate residual
Moran's I in both Visium sources, and residual Moran's I was positive in 3/4 and
4/4 GSE192741 donors and 10/10 and 9/10 Vu arrays for Stromal ECM and Ductular
injury. Stromal ECM exceeded its matched-gene null in both sources (excess 0.040, q=0.001; 0.050, q=0.0006) and Ductular injury did not (0.015, q=0.079; 0.003, q=0.32). That null resamples genes with donors and arrays held fixed, so it is not a donor-level test: 3 of 4 GSE192741 donors have positive Stromal ECM excess (exact one-sided sign-flip P=0.19), and the 10 positive Vu arrays have no public donor key; no prespecified test compared the two programs. The null
resamples genes matched on expression, detection, biotype and correlation with
hepatocyte abundance, with donors held fixed. At the donor level, Stromal ECM
excess was positive in 3 of 4 donors (exact one-sided sign-flip P=0.19); four
donors cannot reach P<0.05 (smallest attainable P 0.0625). The Vu result remains
source-dependent because arrays cannot be joined to donors. Calibration against
an assay's intrinsic spatial propensity, not the sign or apparent replication of
Moran's I, is therefore required before spatial organization is called
program-specific. Matched-gene calibration alone is also not sufficient. In a
descriptive analysis that was not prespecified, 68 of 114 (GSE192741) and 74
of 114 (Vu) fixed programs exceed a matched-gene null built the same way, with
all-16 adjustment and controls matched on each program's own discovery
lineage (BH q < 0.05 over 114 programs per source). Stromal ECM ranks 42nd and
29th by excess, and programs of hepatocyte plasma-protein, acute-phase and
metabolic genes rank highest. Stromal ECM still exceeds the null with controls
matched on correlation with fibroblast abundance (q=0.011; q=0.0008), and
Ductular injury still does not with controls matched on correlation with
cholangiocyte abundance (q=0.10; q=0.50).

The Resource therefore treats spatial organization as an assay-calibrated
transportability and observability layer, not a disease-direction score or an
additional vote in a universal ranking. The Ductular injury program is indeterminate
rather than negative: the available assays do not distinguish it from spatially
matched background.

The accepted Stromal ECM program example is explicitly composition-linked: adjustment for
all 16 inferred cell-abundance factors attenuates, but does not eliminate, its
matched-null spatial organization. Associations between program scores and
cell2location factors are descriptive context because both are derived from
the same Visium RNA matrix; they do not establish lineage origin, mediation,
cell autonomy, or mechanism.

## Why the paper stands on its own

| Resource contribution | What readers receive | Why it matters |
|---|---|---|
| Harmonized human disease-state map | Donor- and cohort-aware bulk and single-cell context across heterogeneous cohorts | Makes cross-study MASLD signals comparable without pretending all cohorts measure the same phenotype |
| Inherited genetic-variation map | Phenotype-, ancestry-, method-, and source-aware GWAS, fine-mapping, coding context, and GWAS/eQTL shared-signal evidence | Connects trait-associated germline signals to candidate genes without assuming that every variant is regulatory |
| Prespecified program registry | All 117 programs, weights, testability, stability, and complete statistical families | Makes program validation reproducible and protects against outcome-driven module selection |
| Noncoding interpretation | GENCODE-versioned lncRNAs, fine-mapped DNA context, target-gene biotype, and object-specific assay coverage | Prevents a noncoding variant, a transcribed lncRNA locus, and an RNA-product mechanism from being treated as the same claim |
| Assay-native validation | Protein, chromatin, single-cell, and spatial evidence retained in their native units | Preserves what each modality can and cannot establish |
| MASLD Gene Catalog | How each gene was nominated, its biological context, source dependence, uncertainty, and suitable follow-up | Makes the Resource directly useful for gene interpretation and study design |
| Reproducible release | Read-only inputs, manifests, source tables, figure provenance, and portal requirements | Lets other groups reuse and challenge the framework |

The scale is valuable, but “largest atlas” and “most modalities” are not the
headline. The contribution is a per-gene record that keeps inherited
genetic evidence and disease-state remodeling side by side, states whether each
layer was untested or tested without support, and makes both usable through
cell and tissue context, a shared set of evidence rules, and the Gene Catalog.

The noncoding addition is a coverage and interpretation advance, not an
enrichment claim. In the validated five-cohort candidate, lncRNAs account for a
similar share of tested genes and interval-null calls, and abundance-adjusted
descriptive analysis does not support preferential lncRNA enrichment. The
contribution is that lncRNAs are tested in the same complete statistical family
as coding genes, checked across cohorts, and routed differently from regulatory
DNA elements when the mechanism remains unresolved.

## Seven-figure story

The one-sentence figure logic of 2026-11-12 at the top of this page takes
precedence over the take-home column below where they differ.

| Figure | Question | Main content | Take-home message |
|---|---|---|---|
| **1. A multimodal human MASLD Resource** | What do inherited variation, disease remodeling, and tissue context each show for a gene? | Cohorts, modalities, biological units, phenotype provenance, the two primary maps, assay coverage, and the shared evidence rules | The Resource records germline signals and multicellular disease biology side by side without treating modalities as interchangeable |
| **2. Inherited genetic variation** | Which trait-associated germline signals and candidate molecular links are supported? | GWAS portfolio, fine-mapped coding and noncoding variation, GWAS/eQTL shared signals, ancestry, method eligibility, DNA context, target-gene biotype, and a GNMT worked locus if retained after corrected COLOC | Inherited variation spans coding, supported regulatory, and unresolved noncoding mechanisms; shared-signal evidence narrows a locus but does not resolve its causal variant |
| **3. Multicohort disease remodeling** | How is established MASLD remodeled across cohorts and cross-sectional disease stages? | Five-cohort differential expression on 774 genetic-identity clusters (not verified donors) and per-cohort DEG membership estimated for all five cohorts; the restored GSE126848 model omits constant sex on 23 control and 7 disease units. GSE126848 enters per-cohort membership through a refit on the 15,592 genes that pass its own expression filter (BH over those genes; the other four cohorts use the fixed 24,115-gene family). The per-cohort union contains 12,180 genes, including 8,341 called in one cohort and 3,839 in two or more, of which 2,628 have one sign in every calling cohort (40 in all five) and 1,211 (31.5%) have opposite signs in different cohorts; it is distinct from the 1,854 pooled calls. Prespecified endpoint PCA and DEG sensitivity, common-reference fibrosis-stage contrasts, continuous NMF axes, supplementary direction/leave-one-cohort-out stability and adjacent fibrosis/NAS views, the GNMT–MAT1A–CYP2C19 RNA–stage–protein state, and the genetics–state interface | The historical fixed-family comparison agrees in direction for 78.3% to 100% of the 1,854 calls within each discovery cohort; native GSE126848 filtering instead gives 937/1,056 same-sign testable calls with 798 untestable. The current full-family REML/modified Hartung–Knapp sensitivity gives 0 BH calls (minimum q=0.171430); direction agreement is descriptive and holds under a stricter histologic endpoint; shared per-cohort membership is not replication, because 1,211 of 3,839 multi-cohort genes change sign between cohorts. Removing one cohort and refitting keeps 88.6% to 96.8% of the 1,854 calls as same-sign calls for four cohorts and 66.7% (1,236) for GSE213621, which alone calls 1,654 of the 1,854 with the same sign and supplies 326 of the 774 units, 68 of the 138 controls and 45.1% of the quality weight; of the lost calls still tested, 505 of 510 keep their sign and effects shrink (slope 0.87). Figure S3C plots the fixed 24,115-gene version (70.7% to 96.6% kept; Spearman 0.870 to 0.995); a severity-linked reduction in hepatocyte metabolic gene and protein abundance connects the GNMT locus to established disease only provisionally (one single-signal ABF result, PP.H4 0.706, in a study where SuSiE was attempted for GNMT but returned no credible-set result) and without turning the three-gene state into one inherited mechanism |
| **4. Multicellular programs and molecular continuum** | Which donor-level programs organize disease remodeling, how secure is their lineage assignment, and what variation remains within recorded fibrosis stage? | Seven-dataset atlas, all-117 source-corrected stage skyline, complete discovery/content lineage comparison, 112-program bulk fibrosis coefficients, descriptive focal trajectories, coefficient strips across all 43 outcome-blind systems and their donor-first stage-versus-continuum comparison; complete transcript/pathway/program comparisons remain supplementary | Histology anchors severity; the continuum describes a shared expression gradient without establishing signature specificity, transcript origin or temporal order |
| **5. Molecular and physical tissue context** | Where are candidate genes and disease programs observable, and in which cell type of intact tissue? | Liver protein (PXD051911) and Xenium in situ cell-type location (Tzouanas et al. 2026, 9 control livers); snATAC and Visium autocorrelation in the supplement | 188 of 211 protein-significant genes match the RNA sign; 22 of 95 programs have protein scores robustly associated with histologic burden (14 lower, 8 higher); in intact control liver the top cell type matches the snRNA lineage for 12 of 14 protein-supported and 8 of 11 colocalized genes. Chromatin enrichment sits in shared peaks and the Visium matched-gene null is passed by most programs, so neither supports a program- or lineage-specific claim. No disease change is claimed from Xenium. |
| **6. MASLD Gene Catalog** | How can the two maps guide interpretation and the next experiment? | Gene and molecular-object entries, evidence states, unresolved alternatives, and deterministic follow-up routes | The Resource turns inherited variation and multicellular remodeling into an experiment-ready field resource |
| **S7 (moved to supplement 2026-11-09). Model framework and transfer limits** | What can current predictive models learn from these maps, and how would we know it is not prior exposure? | A model evaluation framework run across tasks and model families, every score reported beside a per-cohort audit of that model's training exposure; gain over each task's own baseline; the measured limits on tissue-severity transfer and what they follow from | Of six pretrained sequence models scored in seven rows, one (HyenaDNA) beats a sequence-free control on a HepG2 reporter assay after correction; adding a pretrained embedding to distance lowers the score on endogenous eQTL effects; no sequence model was scored on severity; no score is readable without exposure control |

The single-cell model comparison maps five broad cell classes rather than
disease stages. The exact evaluated TranscriptFormer checkpoint has five
source-qualified clean studies and two training-seen studies; its pooled
score and displayed all-study gains include both. Geneformer V2 has no
qualified clean study among those seven. Source qualification changes the
interpretation of the scores, not their values. The separate five-study
TranscriptFormer comparison uses 80 source donor identifiers and the original
donor bootstrap; it remains conditional on the evaluated models and studies.
[Current source correction](../Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/README.md).

GSE126848 was not male-only at source: it included 10 women, all with disease.
The identity screen excluded 24 of 31 disease records and 2 of 25 controls for
insufficient filtered identity sites, including all 10 women. Its retained
23 controls and 7 disease genetic-identity clusters are therefore all male.
The restored per-cohort model omits constant sex and gives 6,147 calls in the
fixed 24,115-gene BH family at FDR <0.05 and |log2FC| >0.5. Native fragment
counts and an independent group-means refit reproduce the result. Unequal
eligibility may limit representativeness; these clusters are not verified donors.
Retained disease libraries are about half as deep as controls; 437 of the 6,147
calls are genes with no reads in any library, and adding log2 library size keeps
2,025 to 3,056 of the 6,147 calls (library size before or after normalization)
and none of the 437, so the GSE126848 calls are not read as disease replication.
Figure 3C instead uses a refit on the 15,592 genes that pass the cohort's own
expression filter (BH over those genes; 3,440 calls, 1,821 up and 1,619 down).

A separate, unadopted Figure 7 candidate focuses this question on two measured
molecular endpoints: held-participant H3K27ac profiles from liver RNA and
signed liver caQTL effects from sequence. Its proposed six panels show measured
variation, matched model increments, measurement-acquisition value, and one
source-qualified locus. The existing severity instrument remains a supporting
reference. Neither the source-exposed chromatin development set nor the spent
Currin fold-0 confirmation supplies independent external validation of the new
comparison. A predicted accessibility effect alone supplies no target gene or
MASLD effect. The candidate's exact task boundaries and current evidence are
in [RESULTS](RESULTS.md#two-model-molecular-value-development-2026-09-22-not-adopted).
No result from this unadopted candidate is reported in the paper.

**Figure 7 redesign plan (user, 2026-11-03; not adopted, not built).** If the
two-model candidate is adopted, Figure 7 asks what models can infer about a
liver sample that was not directly measured. It treats two tasks separately
and draws no link between them: sequence to liver caQTL effect (regulatory
potential) and liver RNA to H3K27ac profile (measured tissue state). Model
names label panels; they do not organize the figure. Every score keeps its
training-exposure label. Planned content:

- A small two-branch schematic of inputs, measured targets and held-out units.
- Main chromatin result: the 99-participant increment of global-plus-local RNA
  over the histology, sex and library comparator (+0.12162 [0.09542, 0.14682]).
  The 22 histology-matched pairs (+0.02211 local over global) are shown as the
  matched check and labelled as fibrosis 0–1 only (20 of 22 pairs are women
  with fibrosis 0; no fibrosis 2–4 pairs).
- Transfer: the released severity instrument in its external cohorts (Figure 7D
  values), with its two limits stated beside it: no advantage over a ridge
  fitted inside each cohort, and its expression transform was chosen on
  GSE268273, so only GSE213621 and GSE276114 test transfer.
- Supplement or text, not main panels, because each states a threshold or a
  confounded test rather than liver biology: the adapter versus matched 2 kb
  and native 1 Mb (3.15% against the 5% margin), sequence-by-RNA versus
  additive (below the +0.01 margin), and the B4 cross-assay H3K27ac transfer
  (both fixed recipes below the training mean in 17/17 GSE128072 normal
  transplant-donor ChIP livers, which earlier jobs had already opened).
- Supplement: single-cell foundation models, H&E and spatial models, other
  severity models, and Model A.

Building waits for adoption of the candidate and for completion of all models,
as stated under figure rules below.

This order starts with the complete Resource, then maps inherited genetic
variation and established multicellular disease remodeling before resolving
their cell and physical tissue context. It preserves the recognizable
transcriptomic core of differential gene
counts, continuous bulk axes, and donor-level gene programs without describing separate human
biopsies as longitudinal change.

The working compact artwork selection uses native slot dimensions and 6 pt
Helvetica with no footnotes. Figure 4C displays all 117 source-corrected stage
coefficients (0 selected; minimum permutation BH q=0.18135). Figure 4D compares discovery
and content-associated lineage across all 117 fixed programs. Figure 4E shows
the complete 112-program weighted linear fibrosis/NAS associations, with BH
across 224 tests; linear NAS results do not summarize nonlinear remodeling.
Figure 4E uses GSE130970, GSE135251 and GSE162694 (364 selected
genetic-identity units), the three 774-unit cohorts that record fibrosis and
NAS. No fibrosis or NAS association meets BH q<0.05 (smallest q 0.079 and
0.159), with or without composition adjustment under either fraction set. With
three cohorts the Knapp-Hartung test (t on 2 df) is low-powered, so Figure 4E
leaves these associations unresolved; it does not show that they are absent. A
sensitivity arm that adds 93 GSE174478 source records (outside the 774-unit
fit, not identity-resolved) gives 16 fibrosis and one NAS call; after
composition adjustment three fibrosis calls remain with the rebuilt fractions
and one with the earlier fractions. All depend on the untruncated
Knapp-Hartung variance (Q below its expectation; none under the truncated
variant) and are not reported as findings. The historical 469-record, 113-program panel is superseded.
For the
composition-adjusted continuum coefficients in the 367 disease units, Ductular injury is
reported as holding after composition adjustment (0.60 to 0.73 across codings
and fraction sets; unadjusted 0.707), Stromal ECM as reduced by composition
adjustment by an amount that depends on the deconvolution input (0.294 to 0.479
in the standard 16-type composition candidate; 0.203 in the separate 13-part
rare-type pooled sensitivity; unadjusted 0.760), and no adjusted system is named. The earlier two-program lineage and
member-gene transport examples remain supporting source results. Working
main Figure 4F shows five-cohort descriptive Stromal ECM and Ductular injury
trajectories over the fixed external projection. Figure 4G shows the
stage-adjusted continuum coefficients, estimated in the 367 disease units,
as two columns of system coefficient strips. The fixed systems were defined
from signature-excluded, outcome-blind gene membership spanning pathways and
Hotspot programs; the membership map remains in Figure S4Q. Its systems
are size-ranked visualization units, not an ontology, and its cross-system and
sliding-window summaries are descriptive because the features and windows are
correlated. The purpose of the systems layer is to give the continuum a unit of
interpretation between single genes and whole-cohort axes: it compresses 13,092
correlated features into a few dozen auditable groupings so that within-stage
remodeling can be described as coordinated rather than gene by gene. That
purpose is organizational, and the stress tests bound it. The groupings are a
navigable vocabulary and a demonstration that within-stage variation is shared
across pathway families and single-cell-derived programs; they are not a claim
that any one system is more disease-relevant than another, and they do not
support nomination or subtyping. Figure 4H shows all 43 systems directly against
fibrosis-stage effects. Both coefficients are estimated in the 367 stage-complete disease units;
the 68 GSE213621 controls have no recorded histology. On the donor-first estimates of the system inference tables, Spearman rho is 0.67 (bootstrap 95% interval 0.38 to 0.84) and 34/43 effects share direction (feature-first summaries: 0.63 and 32/43; 0.7265 and 37/43 when the controls entered as stage 0); both remain descriptive. The size-matched competitive null was
run on the 435-unit roster with the controls, so the panel carries no outline
until it is rerun on disease units. The historical 844-record values were 0.658 and 36/43. The four-test maxT
count no longer carries ink, because it is set by the smaller cohort's
within-stratum permutation resolution and the jointly powered cross-cohort test
supports 33 of 43 on the corrected counts. Supplementary S4V and S4W carry the null distributions and
the complete 43-system competitive check. The corrected 774-unit stress analysis
(`hac-system-stress-corrected774-20261009T171017Z`) gives 33/43 joint maxT calls
on both axes, 0/43 exceeding the connected-region competitive null on both axes
and 1/43 exceeding the unconstrained size-matched null on both axes. These are
different tests; their stage-adjusted roster has 435 units, including the 68
GSE213621 controls, whereas Figures 4G and 4H use the 367 disease units. The supplementary pooled system heatmap gives each cohort equal
weight after within-cohort score standardization; it is not raw participant
pooling or an additional inferential test. Its shared late-minus-early row
order is a display choice used in both system heatmaps, not a support or timing
statistic. A separate row-level donor analysis tests the 43 prespecified
membership-defined systems in
GSE162694 and GSE213621 after fibrosis-stage and sex adjustment: 33 meet the
two-axis fixed-effect meta-analysis BH/direction criteria, and 18 also pass maxT
   family-wise P<0.05 in all four axis-by-cohort permutation tests. Those models
   use 109 GSE162694 and 326 GSE213621 units, 435 in total, where the 326
   include 68 controls with no recorded histology that the stage code places at
   stage 0, not the 466-unit descriptive roster (140
   GSE162694 and 326 GSE213621) on the adopted 774-unit counts; on the 844
   records these were 109, 361, 470 and 501. The 33 and 18 counts are the same
   on both. All 18 strict calls are stage-and-continuum-supported systems; none of
   the nine continuum-only systems clears the four-test conjunction, which is
   limited by the smaller cohort's within-stratum permutation resolution. A
   jointly powered cross-cohort permutation that combines both cohorts before
   maxT correction gives 33 of 43 and includes seven of the nine continuum-only
   systems. The historical 844-record values were 35/43 and all nine.
   The four-test result is therefore partly a resolution limit rather than
   evidence against a continuum association.
   Both criteria are reported. The binding limit is a
   different one: against 5,000 size-matched draws from the same membership
   geometry, no system passes BH on both axes under the conservative connected
   null and one does under the permissive unconstrained null. System counts of
   this kind therefore describe a diffuse shared axis and must not be presented
   as system-specific association. S4U retains
   all 43 systems and marks only those 18 strict calls with a neutral asterisk
   (in the 367 disease units, 4 of the 18 are called in all 11 estimable
   composition-adjusted fits, 14 lose the call in at least one, and 9 are
   called in 5 or fewer);
   the 33 primary BH calls remain in the source registry. The marks are shown
   beside, not derived from, the descriptive windows. Figure 4G and S4Q–S4U
   contain no embedded explanatory footnotes; all method and interpretation text
   remains in the external caption authority. Figure 3 remains unchanged and the full
synchronized Resource release remains open. The ambient
audit is retained as source-table provenance rather than a figure panel. Figure
5 (2026-11-03) is 5A and 5B protein, 5C Xenium cell-type location, and 5D
program protein-histology forest; the former 5A, 5D, 5E and 5F (input
diagram, snATAC context, three-track program summary and Visium Moran
calibration) are supplementary. Figure
6C–D use evidence and experiment-routing matrices rather than prose-heavy
nodes. This is an assembly candidate, not an adopted release.

Spatial artwork follows the same compact grammar as the other main figures:
one decision per panel, direct evidence-state labels, and no methods paragraph
inside the plotting area. Exact estimates, unit qualifications, and claim
boundaries remain in the panel source table and legend. The full 117-program
coverage family remains visible in the supplement, while focused diagnostic
plots may show only the comparison needed to explain the main call.

### Figure priority: biology before statistics

New panels are ranked by the biological finding they show, not by the
analysis that produced them.

1. Highest: one gene or program measured in two or more assays on the same
   biology, where the assays agree or disagree in a stated way. Examples:
   lineage in snRNA, transport to bulk cohorts, location in tissue, protein
   or chromatin at the same gene.
2. Next: a single-assay biological result. It names a gene or program, the
   cell type or tissue compartment, and the direction and size of change
   with disease or stage, on donors.
3. Supplement, source table or one sentence of text: plots about the data or
   method rather than the liver. These include coverage and observability
   counts, power, null distributions, calibration, permutation histograms,
   QC and method comparisons.

Before proposing a main panel, write the biological sentence it supports. If
that sentence is about an assay, a test or a threshold, the panel is not a main
panel. An untestable or negative result goes in the text unless it changes a
biological claim. The user chose the existing main-panel exceptions, which
include the bulk-cohort overlap and endpoint PCA in Figure 3. Do not add others
without asking. Figure 1 is the user's own Illustrator figure; do not generate
panels for it. Figure 7 will be redesigned after all models are built; do not
build Figure 7 panels before then.

## Claims we can make

The authorized initial regulatory-model campaign evaluates signed molecular
effects, donor differences and allele-by-condition interactions as distinct
endpoints. Its development results and local Catalog specimen remain candidates;
they cannot change adopted genetic memberships, fixed117, main figures or
claims without review and the separate adoption decision. Training-calibrated
native-scale error and independent biological units govern new comparisons.

- The Resource records liver-trait genetic evidence and cross-sectional
  disease-state remodeling in established MASLD side by side for each gene, and
  places both in cell and tissue context where an assay can observe them.
- Pooled across ancestries, direct MASLD and liver-fat credible sets have a
  larger descriptive protein-altering plus canonical-splice posterior fraction
  than liver-enzyme proxy credible sets. The pooled contrast is confounded with
  GWAS ancestry and fine-mapping coverage and is not a trait-class claim (H58,
  2026-10-15, recomputed 2026-10-21 on the 4,107 gated credible sets whose lead has
  p < 5e-8, one threshold in every stratum: 106 of 164 direct but 82 of 3,943 enzyme credible sets are
  European; EUR 29.88% vs 6.92%, AMR 8.11% vs 6.91%; no call in any stratum;
  the pooled excess is carried mainly by three PDFF GWAS). Pooled values:
  29.65% against 6.67%, a 4.45-fold ratio (physical-cluster bootstrap
  2.05–9.88), after restoring 715 truncated positions among 6,415 gated credible sets. The
  pooled single-linkage clustering of lead positions within 1 Mb yields 120 clusters: 16 direct, 115 enzyme and 11
  shared. Each restored trailing-zero encoding is uniquely supported within
  1 Mb of every native member coordinate. At genome-wide leads, the
  robustness menu has 51 rows: one reference, 12 direct-study omissions,
  16 enzyme-study omissions, 16 direct-cluster omissions, five clustering
  windows and one LD-block comparison. Over the 50 non-reference comparisons,
  ratios range from 3.44 to 5.18, and every difference interval has a positive
  lower bound (minimum 0.0211). None of these comparisons stratifies by ancestry or trait
  or has matched fine-mapping coverage, so they are correlated sensitivities of
  the pooled value, not replications of a trait-class difference.
  Comparing the trait classes within one ancestry requires SuSiE fine-mapping
  of all 767 UKBB ALT, AST and GGT lead windows and of the MVP European ALT,
  AST and NAFLD GWAS (H58).
  Noncoding-context, trait-directness and ascertainment reruns also completed
  and were adopted for Figure 2 only; H58 (2026-10-15) withdrew the
  noncoding-context unresolved fractions, which reflected failed lookup
  (repaired unresolved 0% in both classes). The historical
  September intervals remain invalid. No permutation p-value is reported;
  the regenerated Figure 2D and unassembled noncoding-context companion use
  the adopted tables, with source-derived captions and separate interval tables;
  the October 13 decision concerned Figure 2's DNA architecture; COLOC r2 was
  separately adopted on October 14. See the
  [adopted coordinate correction and panels](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#figure-2-adoption-and-regenerated-panels-2026-10-13).
  Separately, the colocalized map is predominantly a
  liver-enzyme map: the strongest evidence for 423 of 436 genes comes from a
  liver-enzyme trait and for 13 from a direct disease or liver-fat trait, and
  22 genes have any direct-trait colocalization. Far more enzyme than
  direct-trait gene-by-GWAS pairs were tested by multi-signal colocalization
  (12,121 against 771; 6.1% and 3.6% colocalized, descriptive, pairs not
  independent). The 171 direct and 3,844 enzyme SuSiE signal pairs of the
  shared-posterior rerun were selected on first-pass colocalization, so their
  hit rates do not measure per-signal yield. Some protein-altering variants may
  act through protein function without changing steady-state RNA, while others
  change RNA through nonsense-mediated decay or splicing (Rivas et al. 2015,
  PMID 25954003). Any role of the pooled coding difference in the
  genetic/disease-state overlap is a hypothesis this comparison did not test.
- Among the 3,589 genes eligible under the adopted priority rule and evaluated
  by multi-signal COLOC (SuSiE state positive or tested at PP.H4 ≤ 0.5), the 426
  colocalized and 209 canonical disease-state genes share 26, against 21.396
  expected within 10 AveExpr quantile bins (Mantel-Haenszel OR 1.2815, 95% CI
  0.8262 to 1.9878; exact two-sided P 0.27699; one overlap comparison). This
  replaces the historical 3,769-gene comparison after the shared-posterior
  eligibility check was applied to carried negatives. In the sensitivity that
  applies the shared-posterior check to the carried negatives and counts a gene
  as tested if any tier-1/2 study stays eligible (3,776 genes before the check,
  because the seven SuSiE-untestable genes count as tested), 40 genes become
  untestable, the universe is 3,736, and 26 genes are shared against 21.15
  expected (OR 1.300, 95% CI 0.838 to 2.017; exact P 0.2746). Under the
  priority rule adopted on November 6 it is 3,589, as above. These values use the
  adopted 774-unit bulk refit and COLOC release
  `susie_coloc_r2_20260924T053244Z`, adopted on October 14. The historical
  3,769-gene expression-matched control carried 13 DEGs against 21.3 expected;
  its power analysis found less than 80% power for ORs about 0.46–1.77. Neither
  diagnostic was rerun on 3,589 genes and neither quantifies the updated
  comparison. Expression is not established as fully controlled, and the
  current interval includes 1. Do not claim equivalence or biological
  opposition from this comparison. The
  earlier values counted genes without a multi-signal COLOC result as genetic
  negatives and are withdrawn: 26 against 33.46 expected, OR 0.757, in the
  14,704-gene interface (10,928 such genes, 935 of them canonical DEGs), and
  26 against 33.4 expected, OR 0.76, in the 14,720 symbol-grouped rows (10,941
  such genes, 935 of them canonical DEGs). The 14,704-gene stable-Ensembl interface, the 14,720
  symbol-grouped rows and the 14,502-gene jointly eligible Catalog are distinct
  universes and must not be interchanged. This is reported as a
  description, never as a finding in its own right; any link to the coding
  difference above is a hypothesis. The class comparison was rerun on the
  774-unit fit on 2026-10-14 using the stable-ID joint universe of 14,704
  genes. Matched for expression, the prespecified prediction of greater
  hepatocyte specificity among colocalized genes fails; the post hoc mean
  difference is opposite and attenuated (281 pairs, −0.265515 log2 units,
  nominal paired sign-flip p=0.01482985; 126/155 positive/negative, sign test
  p=0.094677; locus-clustered sign-flip p=0.01606 over 216 loci, post hoc;
  `Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Claims/R14-10/r14_10_recompute.log`).
  Gene dependence remains a limit on that conditional comparison; the P is not
  printed in the abstract or introduction. No difference is detected
  in RNA-effect-adjusted protein response (73 pairs, nominal paired sign-flip
  p=0.487465); this is not equivalence. The preserved relabel/rematch p/q
  values are diagnostics from a mismatched null, not calibrated significance.
  The earlier 272-pair, −0.44, p=2×10⁻⁵ statement is historical.
  Most state-only genes (935 of 1,129; 931 of the 1,124 single-cell-observable
  genes in the comparison pool) have no multi-signal COLOC result.
  Restricted post hoc to the 193 single-cell-observable state-only genes that
  COLOC tested, the difference keeps its direction (129 expression-matched
  pairs, −0.341885 log2 units, nominal sign-flip p=0.01543, seed 20260923;
  sign test 57/72, p=0.218; MDE d 0.247). Within the 281 published pairs, the
  56 with a tested state partner carry most of the difference (−0.748,
  p=0.00148) and the 225 with an untested partner do not (−0.145, p=0.233).
  These P values are nominal and do not account for genes sharing loci; none is
  printed in the abstract or introduction (conflict-7 rerun, job 22257837,
  `Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave3/c7/20261015T185127Z/c7_results.tsv`).
- The established-state transcriptome includes cross-sectional stage-associated
  remodeling in bulk tissue and stable donor-aware single-cell gene programs;
  no single-cell program passes the 117-program family test for a stage
  association (H59).
- An externally anchored molecular continuum complements recorded fibrosis
  stage by describing cross-sectional transcriptional variation
  among participants with the same recorded histologic category; grouped
  categories must not be represented as exact stages. It does not replace histology or
  establish within-person temporal ordering.
- GNMT, MAT1A, and CYP2C19 provide an RNA–stage–protein example of a coherent,
  histologic-severity-linked reduction in hepatocyte metabolic gene and protein abundance in
  established disease; this does not establish a shared upstream mechanism or
  inherited regulatory program.
- The current-official GNMT direct-NAFLD shared-signal observation, aligned
  liver-eQTL direction, and related source-reported functional variants define
  a provisional worked locus, not a resolved causal variant. The adopted r2
  release retains the exact direct-NAFLD ABF pair.
- A preliminary ambient-RNA audit raises a lineage-assignment concern for the
  Stromal ECM and Ductular injury programs. An aligned rerun is independently
  validated but not adopted. Their
  matched-null spatial results are reported per program: Stromal ECM exceeds its gene-level matched null in both Visium sources and Ductular injury does not, and 68 and 74 of 114 fixed programs also exceed that null; at donor level, 3 of 4 GSE192741 donors are positive (sign-flip P=0.19) and Vu arrays have no donor key, so neither source gives a donor-level test; no
  prespecified test compared the two programs. Disease association, tissue
  organization, and lineage origin must be reported as separate claims.
- Protein, chromatin, single-cell, and spatial assays provide distinct forms of
  context rather than interchangeable validation votes.
- Prespecified programs differ in assay observability and physical organization;
  absence of an adequate score or biological-unit join is not a negative
  spatial result.
- Positive spatial autocorrelation, even with consistent direction across
  physical units, requires matched-gene calibration and comparison with the other
  fixed programs before it can be described as program-specific
  organization; most fixed programs exceed the matched-gene null.
- The MASLD Gene Catalog identifies what is supported, what is untested
  because one domain has no call (RNA not detected for 4,232 genes; no Tier-1/2
  COLOC record for 9,312; SuSiE-untestable in one Tier-1/2 study for 7), what
  was tested without a call, and which experiment can discriminate the leading
  explanations.
- The Resource keeps fine-mapped noncoding DNA architecture separate from the
  biotype of a colocalized target transcript and records when no authoritative
  link exists between them. Two candidate fields add to this without changing
  any class: a coding field lists 31 Catalog genes with a fine-mapped
  protein-altering or canonical-splice credible-set member at normalized SuSiE
  PIP ≥ 0.5 (H87), and a regulatory-region field lists 191 of 581 European
  tier-1/2 liver-trait GWAS study-lead signals that colocalize with a liver eQTL
  and share a signal with a bulk-liver accessible region from 138 donors (at
  least 45 already reported by Currin et al. 2025), with 108 Catalog genes (100 when the caQTL and eQTL
  lead variants must have r² ≥ 0.8) where
  that region also colocalizes with the gene's eQTL (H53). Neither field
  identifies a causal gene.
- Disease-associated lncRNAs are included under the same pooled model and
  all-gene multiple-testing family as coding genes; their value here is robust
  inclusion and correct follow-up, not preferential enrichment or proven
  RNA-product function.
- NAS enters as a graded sum and as an adjustment set. The three disease
  aspects of metabolism, inflammation, and fibrosis are not resolved as separate
  outcome axes in the current release, so no aspect-specific attribution is
  available and the state is recorded as indeterminate rather than negative.
- One two-variant local haplotype is characterised at GNMT and reported as
  unresolved. No general claim about combinations of variants is available,
  because the single-variant colocalization design does not test them.
- Chromatin supplies assay-native context for placing candidates. With the large
  majority of jointly testable peaks resolving to indeterminate, no epigenetic
  effect or causal claim is available from accessibility.
- Reanalysis of a published 138-donor liver caQTL map (Currin et al. 2025,
  PMID 40393811) places a colocalized accessible region at 191 of 581 European
  tier-1/2 liver-trait GWAS study-lead signals that colocalize with a liver
  eQTL; at least 45 of the 191 were already reported by the source, so this is
  a credited reanalysis, not a discovery. In 108 Catalog genes (100 with caQTL-eQTL lead r² ≥ 0.8) the same region
  also colocalizes with the gene's eQTL. Among loci selected for liver-eQTL
  colocalization, liver-only loci share a liver caQTL more often than
  platelet-only loci at equal GWAS and caQTL strength (adjusted OR 3.350, 95%
  CI 1.622 to 7.397; 172 against 78 loci), which agrees with the source's
  trait-level pattern and is not presented as new (H53).
- In 222 GTEx postmortem livers, across the 1,002 canonical disease-state genes
  measurable in GTEx liver, the slope with image-measured histologic fat
  correlates with the biopsy disease-versus-control logFC (Spearman rho 0.638,
  Freedman-Lane p 0.0002 at the permutation floor), also in the 138 donors
  outside the ventilator death class (rho 0.632), and more strongly for fat
  than for ischemic time (H55). The claim is limited to histologic steatosis of
  unrecorded cause in postmortem liver; the tissue-area control came close to
  its stop rule (rho 0.354, p 0.0564); adding tissue area to the model leaves
  rho at 0.631 (post hoc; difference -0.0065, paired bootstrap 95% interval
  -0.033 to 0.004), and PLIN2 and the 197-gene steatosis composite met their
  positive-control predictions.
- The number of bulk disease-versus-control calls is sensitive to genetic
  ancestry, while their log2FC estimates change little. Adding four ancestry
  PCs changes 233 of the 1,854 calls against a permuted-PC median of 120, and
  restricting to units of European genetic similarity changes 570 against a
  matched-subset median of 409.5 (Holm p 0.0199 for each, family of two). Over
  the calls, log2FC with the PCs correlates with the adopted values at Pearson
  0.9972 (slope 0.984, descriptive). 64 genes, GSTM1 among them, carry a
  descriptive ancestry-sensitive label against a median of 6 expected false
  labels. The adopted fit remains the release, and no wording may say that
  ancestry causes disease-state expression (H79).
- A missing bulk call is not evidence of no change in low-expression genes. In
  expression deciles 1 and 2 (4,323 of the 22,261 tested genes without a call)
  the pooled fit called a uniform 1.0 log2FC change in 0.507 and 0.796 of
  spiked genes, so those genes are below the detection limit rather than
  unchanged; in deciles 3 to 10 the lower bound is at least 0.98. The pooled fit
  called 0.0004 of unchanged spiked genes (H81).
- Leave-one-cohort-out loss mainly reports GSE213621. In the spike-in, an effect
  confined to one cohort was called in the pooled fit in 87.45% of replicates
  when that cohort was GSE213621 and in at most 0.25% for any other cohort
  (1.0 log2FC), and 99.8% of supported one-cohort effects (2.0 log2FC) lost
  their call in at least one LOCO refit. This does not show that LOCO retention
  proves cohort independence for real genes; none of the 828 calls not kept in
  all five refits is relabelled (H81).
- The Gene Catalog keeps, for each gene, whether a layer was untested or tested
  without support. For the Kumak et al. (2026) lists, 661 of 902 NAS genes and
  714 of 943 fibrosis genes present in both the bulk and colocalization tables
  have no multi-signal COLOC result. Among evaluated genes, 36 of 241 NAS genes
  and 28 of 228 fibrosis genes colocalized, against 28.9 and 24.9 expected
  within expression deciles (MH OR 1.32, 95% CI 0.91 to 1.92; OR 1.16, 0.76 to
  1.76; Holm-adjusted P 0.29 and 0.51), so their COLOC rate did not differ
  detectably from that of genes of the same expression level (H71). Planted
  shifts of OR about 1.74 or 0.48 (NAS) and 1.75 or 0.48 (fibrosis) were
  detected in 82% and 65%, and 77% and 62%, of simulations; ORs between about
  0.48 and 1.75 have below 80% power unstratified. The expression-matched
  random-list negative control missed its predicted range but stayed inside
  the stop limits.
- The Resource can prospectively define hypotheses for a separate functional
  screen.
- On a synthetic reporter assay in two HepG2 contexts, one of seven pretrained
  sequence-model rows (six models), HyenaDNA, exceeds a ridge on allele identity
  alone after Bonferroni correction across the seven rows. Against sequence
  baselines (CpG ridge, flank and dinucleotide ridge) the gain does not survive
  the same correction. On the endogenous variant-to-gene effect, adding
  either pretrained embedding to distance lowers the score, and permuting the
  correspondence between a variant and its embedding does not degrade the score.
  No pretrained sequence model was scored on tissue-level severity. These are statements about
  measured transfer on these cohorts, not about model quality in general and not
  about clinical use. The endogenous arm carries its own scope limit: its labels
  are marginal, un-fine-mapped associations, so the design cannot separate a
  model failure from a target that is not sequence-determined at the variants
  tested.
- Model scores in this Resource are reported beside a per-cohort record of what
  each model saw during its own training. Where that record is unresolved the
  score is reported as uninterpretable rather than as evidence of skill.
- The limits on tissue-level severity transfer are stated as measured bounds
  rather than as failed attempts: where an upper bound over a family of
  strategies falls below the improvement declared worthwhile in advance, the
  family is reported as excluded on these data.
- The evaluation framework reports what a model can and cannot do on a given
  input, including the sample sizes below which it declines to answer. This
  characterizes the framework's own limits and is not a performance claim for
  any clinical or diagnostic use.

## Claims we cannot make

| Avoid | Use instead |
|---|---|
| “causal gene” from COLOC | “colocalized” or “genetically anchored candidate” |
| longitudinal change from unrelated stage groups | cross-sectional stage-associated remodeling |
| universal convergence score or definitive target rank | gene class plus Catalog entry and assay-specific context |
| a molecular system is continuum-associated because its coefficient differs from zero | continuum-associated relative to size-matched regions of the same membership geometry, or else reported as part of a diffuse shared axis |
| a count of supported systems as evidence of that many findings | count plus the effective number of independent systems in the same family |
| genetic null when the eQTL or tissue is weakly observed | untested, weakly observed, or indeterminate |
| discrete, reproducible patient subtypes from k=6 NMF | continuous k4/k6 interpretive axes |
| prospective biomarker performance | exploratory biomarker application or worked vignette |
| hepatocyte-autonomous mechanism from mixed tissue | lineage-associated or tissue-context evidence |
| lineage origin from program–cell2location correlation in one RNA matrix | descriptive spatial context covariation requiring orthogonal confirmation |
| source-owned public finding as our discovery | source finding plus our prespecified cross-resource interpretation |
| variant inside or near a lncRNA proves lncRNA mediation | noncoding locus requiring RNA-depletion and locus-level disambiguation |
| disease-associated lncRNA proves RNA-product function | robust lncRNA association that requires RNA depletion, neighboring-gene readouts, and rescue or locus controls |
| “largest” or “most comprehensive” without a dated audit | broad, multimodal, reusable Resource |
| sparse genetic/disease-state overlap presented as a finding, or as “separable”, “complementary” or “biologically opposed” layers | overlap of 26 genes against 21.2 expected within expression deciles, with expression not fully controlled and no direction stated; protein-altering architecture and expression-based colocalization observability are separate descriptions, and their role in the overlap is an untested hypothesis |
| “PIP mass” or “coding mass” | proportion of fine-mapped posterior probability carried by protein-altering variants |
| aspect-specific attribution from a NAS sum or fibrosis stage | severity-associated, with the resolved aspect recorded as indeterminate |
| combination or interaction effect from single-variant colocalization | single colocalized variant; combinations untested by this design |
| epigenetic causation or mechanism from accessibility | assay-native chromatin context, with peak-level evidence state |
| “coordination of multicellular remodeling” or coordination between cell types | bulk program-pair covariance, reported in Results only with the two further H54 bulk cohorts and their limits; bulk covariance mixes cell abundance and biopsy content |
| a trait-class difference in coding or splice architecture | pooled descriptive fraction, confounded with GWAS ancestry and fine-mapping coverage; no trait-class call (H58) |
| a stage-associated single-cell program | 0 of 117 programs pass the family rule (H59) |
| direction agreement between ALT TWAS and disease-state expression | no direction stated (H60 gives no direction) |
| liver caQTL colocalization at GWAS signals as this paper's finding | credited reanalysis of Currin et al. 2025; at least 45 of 191 signals already reported |
| the 1,854 bulk calls are robust to genetic ancestry, or ancestry causes disease-state expression | the call count changes with ancestry PCs beyond permuted PCs (233 against a median of 120, Holm p 0.0199) while log2FC estimates change little; 64 genes carry a descriptive ancestry-sensitive label (H79) |
| GSE135251 per-cohort membership is ancestry matched | not matchable: 2 of its 10 controls are of European genetic similarity (H79) |
| a gene without a bulk call is unchanged in disease | below the detection limit in expression deciles 1 and 2; elsewhere not called at the 0.5 floor; no bulk `tested_negative` call exists (H81) |
| retention in all five LOCO refits proves cohort independence, or a LOCO loss marks a call source-dependent | one-cohort effects are called mainly through GSE213621 and lose a LOCO call in 99.8% of simulated cases; uniform 1.0 log2FC effects keep all five in 90.6%; no call is relabelled (H81) |

## Two-paper firewall

| Standalone Resource paper | Separate experimental follow-up screen paper |
|---|---|
| Human molecular maps and observability | Final coding and noncoding target roster |
| Prespecified programs and evidence classes | Guide pool and screen QC |
| Assay-native validation and MASLD Gene Catalog | Mouse pooled-screen phenotypes and hits |
| Experiment-routing rules | Mechanistic follow-up and validation |
| May prospectively motivate screen hypotheses | May not retroactively choose Resource classes, thresholds, hero genes, or figures |

The Resource must remain publishable if the future screen is delayed, negative,
or biologically surprising. Likewise, the screen paper should be able to say it
tested hypotheses specified prospectively by an independent Resource.

## Competitor-aware positioning

Recent MASLD atlases overlap with individual components—large RNA-seq
collections, cross-sectional stage modeling, biomarker exploration, regulatory
variant interpretation, or spatial cell-state analysis. The defensible
difference is the combination of:

- explicit separation of inherited genetic variation from established tissue
  state;
- complete observability and joint-testability accounting;
- donor-aware, assay-native integration rather than modality vote counting;
- explicit separation of noncoding DNA, lncRNA transcription, and mature-RNA hypotheses;
- prespecified program families and source-dependence labels; and
- a gene-level Catalog that records context, limitations, and suitable follow-up.

Each component has precedent: staged bulk transcriptomes and their
meta-analysis with gene-level GWAS (Piras 2024), liver caQTL colocalization
with GWAS signals (Currin 2025), a liver eQTL meta-analysis for cardiometabolic
traits (Broadaway 2024), coding plus expression gene nomination
(Sveinbjornsson 2022), and assignment of MASH risk loci to cell types (Elison
2025). None of these is claimed as new here.

Detailed comparisons and source papers live in [literature/REVIEW.md](literature/REVIEW.md).

## Journal posture

The paper should be submitted and reviewed as a **Resource/framework**, not as
an incomplete mechanism paper. Cell Genomics is the clearest working fit;
Genome Biology and Nature Communications are credible alternatives depending
on the final coherence, portal, and biological examples. Cell Metabolism or
Nature Metabolism remain aspirational Resource inquiries, but the manuscript
should not overstate causality to imitate a mechanistic Article.

Journal names and article categories must be checked against current author
instructions immediately before submission.

## Analysis triage and stop rules

An analysis belongs in the main paper only if it improves at least one of:
coverage, observability, evidence interpretation, cross-assay transportability,
or experiment routing—and if it can be completed without weakening biological
unit control or reopening outcome-driven selection.

Stop or move to the supplement when:

- the public metadata do not identify donors or authoritative assay joins;
- a result depends on cells, spots, guides, runs, or sections as biological
  replicates;
- the conclusion disappears after expression, cohort, cell-composition, or
  source-dependence adjustment;
- the only novelty is rediscovering a source paper’s headline;
- the analysis requires a new permission-dependent cohort;
- the result would only matter if a future follow-up screen outcome is known; or
- a valid null has already answered the prespecified question.

Plans 40–44 are complete supplementary stress tests. Plans 45–46C describe
possible future experiments and collaborations; they are not dependencies of
this paper.

## Release acceptance

The Resource is ready only when one synchronized release contains:

- the chosen five-cohort bulk substrate and its stage/composition derivatives;
- the corrected COLOC portfolio and rebuilt evidence classes;
- the prespecified 117-program registry and assay-native context tables;
- the MASLD Gene Catalog, data dictionary, and portal export;
- read-only figure source tables and legends;
- one reconciled number and claim ledger; and
- explicit confirmation that no follow-up screen outcome influenced the Resource.

The current state of each item is in [STATUS.md](STATUS.md); exact numbers are in
[RESULTS.md](RESULTS.md); execution order is in [ROADMAP.md](ROADMAP.md).
