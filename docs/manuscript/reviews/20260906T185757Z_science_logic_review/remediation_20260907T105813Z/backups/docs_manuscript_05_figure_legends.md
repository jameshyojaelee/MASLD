# Figure legend requirements — standalone MASLD Resource paper

> **Status: working legend skeleton, not an adopted release.** Plan 60 must
> replace bracketed panel details and provisional counts from read-only figure
> source tables, then verify exact agreement among PDFs, sidecars, Methods,
> Results, and `docs/RESULTS.md`.

Across Figures 1–6, every claim links to assay applicability, biological unit,
coverage and join gate, machine-readable evidence state, state reason,
unresolved alternative, and next discriminating experiment. Panels do not
encode modality count as confidence.

## Figure 1. A multimodal human MASLD Resource

**(A)** Cohort-and-assay map for the human MASLD Resource, including biological
units, phenotype provenance, source reuse, and assay-specific applicability.
Raw coding-versus-lncRNA coverage is descriptive; no abundance-adjusted
cross-assay comparison is made without complete versioned eligibility tables.
**(B)** Interface between
jointly testable liver-trait regulatory evidence and established disease-state
transcription; exact class counts are inserted from the adopted BG-001 arm.
**(C)** Assay-observability and evidence-state guardrails distinguishing
supported, discordant, indeterminate, untestable, not applicable,
source-dependent, and prespecified adequate-negative evidence. The panel includes CosMx 12/117
complete-program observability, GeoMx and GSE287826 0/117 testability, ATAC
state counts, and zero strict Catalog `tested_negative` calls. **(D)**
Six-figure analysis flow from read-only inputs through one shared set of claim,
coverage, state, unresolved-alternative, and next-experiment requirements to the MASLD Gene
Catalog. Percentages
use jointly testable denominators; technical units are never counted as donors.

## Figure 2. Phenotype- and ancestry-aware inherited genetic variation

**(A)** Prespecified 35-stratum Tier-1/2 GWAS portfolio across five GWAS
ancestries and six liver traits, followed by 360 trait-specific fine-mapped
loci (265 globally unique loci) and 1,013 COLOC genes. Brackets bind the
endpoint bands to their separate multi-signal COLOC (462 genes) and
single-signal COLOC-only (551 genes) denominators. The molecular-QTL panel is
European, so European GWAS are ancestry-matched and non-European GWAS are
cross-ancestry eQTL-limited.
**(B)** Fine-mapping PIP versus multi-signal COLOC posterior for genes with
finite high-resolution results; absence from this panel is not a negative
colocalization call.
**(C)** The 1,013-gene multi-signal-or-single-signal COLOC union partitioned by
the strongest consequence annotation of one lead COLOC variant per gene: 36
coding, 73 UTR, 107 promoter, 455 intronic, and 342 intergenic. This is a
gene-level lead-variant summary, not the posterior-mass composition of all
credible sets. **(D)**
Locus-balanced allocation of fine-mapping posterior probability across four
mutually exclusive consequence classes for direct MASLD/liver-fat and
liver-enzyme traits. The classes sum to 100% within each trait group. Sentinels
within 1 Mb are grouped to limit repeated-study dominance; this does not assert
LD or statistical independence. The panel describes PIP allocation, not
heritability, effect size, variant count, or gene-causality probability.
**(E)** Global ancestry partition of the 462-gene multi-signal COLOC evidence
set: 352 EUR-only, 70 shared with EUR, 26 non-European-only but sub-threshold at
the lead variant, and 14 non-European-only genes passing lead-variant
genome-wide significance (13 EAS and 1 AFR). **(F)**
Genetics-only gene-by-ancestry dot matrix showing all 13 genes with PP.H4 > 0.5
in four of the five GWAS ancestry panels; no gene is supported in all five.
Rows are ordered by multi-signal breadth, mean multi-signal PP.H4, and gene
symbol. Circle area denotes the displayed
PP.H4. Dark blue uses multi-signal SuSiE PP.H4; light blue uses single-signal
ABF PP.H4 when it alone exceeds 0.5; gray uses the larger available posterior
when neither model exceeds 0.5. Exact posterior, trait, method availability,
lead variant, LD source, GWAS ancestry, and eQTL ancestry
remain in the source table. No displayed cell is not evaluable, so that unused
category is omitted from the artwork. COLOC identifies shared-signal candidates; it does
not by itself establish causality or biological differences between
ancestry groups. **(G-H)** Matched locus views showing GWAS association,
fine-mapping, COLOC, expression, and local gene context. The canonical examples
are RORA and FABP1. Both use EUR GWAS and do not provide evidence for
population-dependent biology. **(I)** Shared legend for the two locus panels.
**(Companion panel, not assembled in canonical A-I)** Phenotype provenance of the 35 prespecified Tier-1/2 GWAS
strata, separated by direct MASLD/MASH diagnosis, MRI-PDFF or histologic
steatosis, and liver-enzyme traits, with European ancestry-matched and
cross-ancestry eQTL-limited pairings distinguished.

## Figure 3. Bulk transcriptomic remodeling and genetic overlap

**(A–B)** Cohort coverage and documented fibrosis-by-NAS support from the
synchronized five-cohort fragment-count substrate. **(C)** Per-cohort direction
agreement among the 1,347 canonical disease-state DEGs and all-gene effect
correlation after leaving out each cohort. The same-direction fraction ranges
from 79.8% to 100%; leave-one-cohort-out Spearman rho ranges from 0.867 to
0.992. Participant is the inferential unit. **(D)** Group-blind PCA of the
prespecified endpoint roster: 115 strict controls and 303 advanced-disease
participants across four cohorts. Advanced disease is NAS at least 5 or
fibrosis at least F3; histology takes precedence for one source-label conflict.
Strict controls are source-defined controls with observed NAS no greater than
1 and observed fibrosis F0; missing histology is retained only when no observed
value contradicts control status. Cohorts require at least five participants
in each arm. Expression is adjusted for cohort and inferred sex without
protecting endpoint status. Both density and scatter use the same 418
participants; controls are open gray circles and advanced disease follows the
pale-red-to-magenta fibrosis ramp. The canonical disease-state model remains
the five-cohort comparison of all 844 participants. Its all-sample PCA and the
endpoint DEG sensitivity are shown in Figure S3. **(E)** Fibrosis-stage
differential-expression counts at F1–F4 versus the common F0 reference and
continuous F0-centered NMF activity. GNMT, MAT1A, and CYP2C19 provide the recurring
RNA–stage–protein example of reduced hepatocyte metabolic function. BH correction is applied across each complete
all-gene family; NMF axes are not patient classes. The three proteins also
covary within 46 MASLD participants after adjustment for acquisition batch,
age, BMI, and sex; histology-adjusted sensitivity interprets this as a
severity-linked state, not an autonomous molecular circuit. The full two-track panel plus
donor-aware hepatocyte and non-parenchymal composition, adjacent fibrosis
contrasts, and NAS views are shown in Figure S3.
**(F)** Genes that are BH-significant in at least one adjacent fibrosis
contrast and belong to the adopted primary colocalized set. Signed stage
effects, exact-pair PP.H4, and exact-linked lead GWAS SuSiE PIP are displayed
in separate lanes. Unresolved PIP linkage is missing, not zero; PP.H4 and PIP
are never combined. Exact counts and labels await the adopted corrected-COLOC
manifest.
**(G)** Threshold-free relationship between raw RNA-seq log2 fold change
and SuSiE COLOC PP.H4 for 3,617 genes with finite estimates from the synchronized
fragment-count and adopted corrected-COLOC releases. No fitted line is shown
(Spearman rho=0.014; P=0.40). Point colors use four matched evidence labels:
Neither, Transcriptomics sig., Genetics sig., and Both. Transcriptomics sig. is
defined by the calibrated TREAT family (`lfc=0.25`, BH q<0.05), without a
post-FDR raw fold-change cutoff. The Figure S3G matched-pair matrix compares 428
Genetics sig. genes with one-to-one expression-matched controls: 30 pairs are
Transcriptomics sig. only in the Genetics sig. member, 21 only in the control,
377 in neither, and none in both (paired OR=1.43; exact McNemar P=0.26). PP.H4
is unsigned, so this is not a risk-allele direction test. Inherited regulatory
susceptibility and established disease state are complementary without evidence
of biological antagonism.

## Figure 4. Donor-level single-cell programs and multicellular context

**(A)** Seven single-cell datasets, 104 descriptive biological donors, 1.23
million cells, one annotated atlas, 102 analyzed donors, and 117 prespecified
programs across five tested lineages. **(B)** Integrated UMAP with direct
cell-type labels and donor-averaged disease density; sequencing runs and sort
fractions are not treated as donors. **(C)** The three positive and three
negative Hotspot effects with the smallest complete-family BH q per lineage,
shown across Healthy, Steatosis, and Steatohepatitis and paired with signed
disease-stage β estimates and 95% confidence intervals. The complete 117-program
registry view is shown in Figure S4. In the original uncorrected atlas, only the
Stromal ECM and Ductular injury programs pass the primary 117-program
BH/stability rule; neither passes HC3 family correction. **(D)** Donor-paired
same-atlas lineage-minus-hepatocyte score contrasts for those two programs.
Error bars are HC3 95% confidence intervals and BH correction spans the ten
prespecified contrasts. The hepatocyte reference is the zero line and is not
replotted; both programs use the same effect axis. Stromal ECM/*IGFBP7* is
highest in fibroblasts (53
paired donors across four datasets), while Ductular injury/*BICC1* is highest
in cholangiocytes (51 paired donors across five datasets). This is descriptive
localization, not transcript origin or cell autonomy. **(E)** Fixed
L1-weighted member-gene bulk log2 fold changes for F1–F4 versus F0. Uniform
points mark the eight program-by-contrast summaries; no program-significance
encoding is shown. This is cross-sectional
tissue-state transport, not longitudinal within-person change or a program-level
significance test. No inferred cross-cohort fibrosis stage enters a primary
single-cell panel. **(F)** Cohort-standardized Stromal ECM and Ductular injury
scores across nine prespecified fixed-projection windows in all five Resource
cohorts. GSE126848, GSE130970, and GSE135251 overlap continuum construction and
are shown as dashed gray descriptive trajectories; GSE162694 and GSE213621 are
the competitor-independent evaluation cohorts. The right side reports
stage/sex-adjusted fixed-effect meta-analysis estimates from only those two
evaluation cohorts for both co-primary continuum scores. Window means do not
determine statistical support or imply temporal ordering. GSE126848 lacks
recorded fibrosis and is omitted only from stage-stratified supplements.

## Figures S3–S4. Histology-versus-continuum comparison

Matched evaluation-cohort inferential panels compare standardized linear fibrosis-stage
trends with stage-adjusted fixed-projection coefficients for the complete
transcript, pathway, 117-program Hotspot, and ten-axis NMF families. Separate
descriptive panels show the nine fixed continuum-window means in all five
Resource cohorts and raw fibrosis-stage means in the four cohorts with recorded
stage for the prespecified focal programs, six fixed Hallmarks, and all ten NMF
axes. Source-overlap cohorts remain descriptive and do not enter the primary
meta-analysis. Transcript results are continuum-associated transcripts rather than
continuum DEGs. All comparisons use the same expression substrate and are
alignment analyses, not independent validation or evidence of within-person
temporal ordering.

The molecular-systems supplement additionally contains an equal-cohort pooled
heatmap of all 43 membership-defined systems across the same nine aligned
continuum windows. Collections are balanced within cohort and the five cohort
summaries then receive equal weight. This is a descriptive cross-cohort summary,
not raw participant pooling or a new inferential family; within-cohort score
standardization and percentile definitions are retained. The pooled panel is
ordered by the descriptive pooled late-minus-early contrast, defined as the
mean of windows centered at 80 and 90 minus the mean of windows centered at 10
and 20, and shows that contrast in a final `Delta` column. The cohort-specific
panel uses the same row order. This ordering is presentation only.

## Figure S4. Single-cell program and communication detail

**(C)** Complete 117-program Hotspot stage heatmap. **(D)** Fixed 13-pair
donor-collapsed communication roster across Healthy, Steatosis, and
Steatohepatitis; zero of 13 pass donor-level family FDR. **(E)** Donor-matched
screen of the two Figure 4C starred programs against the fixed communication
roster. Twelve estimable pairs are displayed; the two untestable
*LAMC3*→*DAG1* tests remain in the complete 26-test family and source table.
One association passes both model-based and HC3 BH correction: the Ductular
injury program with hepatocyte→endothelial *CDH1*→*PTPRM*. The
right subpanel shows partial residuals for the 60 biological donors underlying
that adjusted coefficient. This same-atlas association does not establish
ligand secretion, receptor activation, or mechanism.

## Figure 5. Assay-native molecular and physical context

**(A)** Input separation between gene-level context, the prespecified
program display set from the complete 117-program registry, and the
fixed descriptive protein display. **(B)** Adjusted liver DIA-MS context for
all measured prioritized genes. The scatter retains the nonsignificant
background; the adjacent bar includes protein-significant genes only and shows
concordant and discordant counts and within-significant percentages. The full
measured denominator and nonsignificant count remain in the source table and
caption. GNMT, MAT1A, and CYP2C19 are lower at both RNA and protein levels and
are interpreted as a severity-linked hepatocyte metabolic-function state, not as an inherited
regulatory program. **(C)** Selection-conditioned 25-protein descriptive process view.
**(D)** Normalized SNP.PP.H4 posterior mass for the primary signal pair of each
adopted colocalized gene-study pair intersected with static open chromatin in
two liver snATAC sources across five shared lineages. Bars separate direct
MASLD/liver-fat pairs (n=29) from liver-enzyme pairs (n=787) and retain
both-source, one-source, partial, and no-overlap states; all pairs met the 0.95
mapped-posterior testability criterion, so the zero untestable state remains in the
source table. Peak overlap localizes regulatory DNA but does not assign a
target gene or establish functional regulation. **(E)** Assay-native program
evidence with complete-family correction and testability states. **(F)** At
left, a registered H&E overlay shows the all-cell-adjusted Stromal ECM score in
one Guilliams et al. steatotic section selected outcome-blind as closest to the
dataset median eligible spot count. At right, points show all-cell-adjusted
Moran excess over matched genes for the prespecified Stromal ECM and Ductular
injury programs in Guilliams et al. and Vu et al.; labels are BH-adjusted
empirical q-values from 9,999 draws. Stromal ECM is supported in both sources,
with Vu et al. source-dependent because its arrays lack a public donor key,
whereas Ductular injury is indeterminate in both. The overlay is illustrative;
no pathologist region masks are available,
so it does not establish localization to a named histologic structure. IGFBP7
is one of 49 fixed Stromal ECM
genes (2.47% L1 weight), and BICC1 is one of 22 Ductular injury genes (5.44%);
the inference applies to weighted program scores, not either gene individually.
These spatial states do not establish lineage origin. Effects
with different units are not summed or ranked; spatial panels do not constitute
a MASLD-versus-healthy spatial disease contrast.

## Figure 6. The MASLD Gene Catalog

**(A)** Catalog structure linking each nominated biological object, including
protein-coding gene, lncRNA transcript, regulatory DNA element, and ambiguous lncRNA locus, to
phenotype, assay, biological unit, tested universe, provenance, and source
dependence. **(B)**
Coverage and observability matrix showing why missing, untestable,
source-dependent, and tested-negative evidence are distinct. **(C)**
Outcome-blind Catalog entries chosen to demonstrate different gene states,
including the GNMT–MAT1A–CYP2C19 RNA–stage–protein example and the contrasting
ambient and spatial states of the Stromal ECM and Ductular injury programs, not to
create a greatest-hits list. **(D)** Deterministic routing from the
unresolved alternative to the smallest discriminating measurement or
perturbation. RNA depletion, locus-level CRISPR interference or allele editing,
neighboring-gene measurement, and rescue are kept as distinct routes. The
molecular-object panel separates regulatory DNA, a disease-associated lncRNA,
a genetically linked or ambiguous lncRNA locus, a protein-coding gene, and an
untestable nomination. **(E)** Portal and downloadable source-graph requirements. The
historical convergence heuristic is optional supplementary provenance and does
not determine default ordering or target validity. No Cas13 screen outcome,
roster, guide, hit, or performance metric appears in this figure.
