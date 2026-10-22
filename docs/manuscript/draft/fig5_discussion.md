## The MASLD Gene Catalog records evidence class, observability and a follow-up experiment for each gene

To make the Resource usable gene by gene rather than as a ranked list, we built
the MASLD Gene Catalog (Figure 6). For every eligible gene and program, a Catalog
entry records the primary gene class, phenotype provenance, assay and biological
unit, testing universe, source dependence and direction where interpretable. It
also records a machine-readable evidence state and the reason for it, the
unresolved alternative explanation and a deterministic next-experiment rule. The
Catalog covers 28,233 gene identities and keeps the eQTL-colocalization call
separate from the bulk RNA call (Figure 6B). The genetic call is supported when
SuSiE-COLOC PP.H4 exceeds 0.5 in any of 35 Tier-1/2 GWAS. The RNA call is
supported at BH *q* < 0.05 and |log2FC| > 0.5 in one limma-voom fit of 774
genetic-identity clusters, which are not verified independent donors. The state
vocabulary distinguishes supported, discordant, indeterminate, untestable, not
applicable, source-dependent and prespecified adequate-negative states. The
Catalog gene calls use three of them (supported, indeterminate and untestable).
No nonsignificant result is converted into a negative call. In this paper,
unsupported means tested without support (indeterminate), not tested negative.

A blank evidence cell can mean that the gene was not detected, the source was
underpowered, the assay cannot observe the proposed mechanism, the dataset lacks
an authoritative biological-unit join, or the experiment has not been done
(Figure 6B). These states have different scientific consequences, and the
Catalog records them with the following limits. None of the Catalog's 56,466 gene-level evidence entries is classed as not
applicable. Every entry also carries a source-dependence flag, separate from its evidence
state, so that flag does not separate one gene from another.
The 9,312 genes with no Tier-1/2 COLOC record have an untestable genetic call, and
the Catalog gives the reason as no Tier-1/2 COLOC record. No power assessment was
made for them, so the call marks a missing record, not a measured lack of power. The Catalog defines no negative-call rule, so no gene is
recorded as tested negative. The genetics call also does not separate evidence
for distinct causal variants from absence of evidence. Among the 816 Tier-1/2
gene-study pairs rerun with a shared-posterior eligibility check, 456 genes have
at least one eligible SuSiE signal pair. Of these, 20 have PP.H3 above 0.8 and
best PP.H4 at most 0.5, the pattern expected when the trait and expression
signals have distinct causal variants. All 20 had PP.H4 above 0.5 in the first
run and lost it when the check removed signal pairs, so the pattern was
produced by the check and is not an independent observation. Eighteen of these 20 are called
indeterminate, the same label given to genes with no evidence for either a shared
or a distinct causal variant. The same pattern is common when every Tier-1/2
SuSiE result is counted, including negative results carried forward from the
first colocalization run without the shared-posterior check. Of the 3,918 genes
with any finite Tier-1/2 SuSiE result, a set that does not require a bulk RNA
test and so is larger than the 3,769-gene comparison universe, 3,335 have PP.H3
above 0.8 in at least one stratum and best PP.H4 at most 0.5. Of these 3,335,
3,319 are called indeterminate. When the shared-posterior check was also applied
to the carried negatives, 3,877 genes kept a finite Tier-1/2 SuSiE result and
3,306 showed this PP.H3 pattern. This count keeps the first-run PP.H3 of
gene-study results in which only some signal pairs failed. Their PP.H3 after
the check is unknown, so the count lies between 3,101 and 3,416. Under the
Catalog's rule for untestable genes, 3,148 of the 3,306 would be called
indeterminate (range 2,948 to 3,252). The indeterminate genetic call also covers
15,009 genes with a Tier-1/2 COLOC record but no SuSiE result; for them the
multi-signal test is missing, not negative. Of these, 333 have a single-signal
(ABF) PP.H4 above 0.5: 290 unresolved, 18 established-state and 25 untested
genes, TM6SF2 (0.64) among them.

The Catalog also carries annotation fields that describe each gene without
changing its class or route. For each of the 28,233 genes, it records whether
the protein was detected in liver (58 participants) and plasma (143
participants) in the PXD051911 proteomes. It also records whether the gene is
targeted by the CosMx panel of GSE312698 (Boesch et al., 2026), and in which of
nine liver cell lineages its RNA is detected across donors. Where an assay did
not observe a gene, absence is recorded as untestable, not negative. A liver
protein row was quantified for 173 of 1,129 established-state genes, 203 of 400
genetically anchored genes and 11 of 26 genes supported on both axes. Each gene
also carries its position along the central-portal axis in five non-steatotic
live-donor livers (Visium; Yakubovsky et al., 2026) and its hepatocyte protein
zonation (Weiss et al., 2026). Across 1,337 genes measured in both, the two
agreed in direction (Spearman *rho* 0.39, 95% interval 0.34 to 0.44; positive in
each of five donors). Weiss et al. reported a comparison of RNA and protein
zonation against single-nucleus RNA (their Supplementary Fig. 9). Each
disease-state call also records whether the gene changed with image-measured
histologic fat in 222 GTEx postmortem livers. Of the 1,002 measurable calls, 326
did so at BH *q* < 0.05, and 323 of these changed in the direction seen between
disease and control biopsies. This field concerns histologic steatosis of
unrecorded cause in postmortem liver. It does not assign any biopsy call to
steatosis.

On the genetics side, 31 Catalog genes carry a fine-mapped protein-altering or
canonical-splice credible-set member at normalized SuSiE PIP of 0.5 or above,
and 20 of them are in the unresolved class. This field does not change the
eQTL-colocalization class and does not identify the gene through which a variant
acts. Of 581 European liver-trait GWAS signals that colocalize with a liver
eQTL, 191 share a signal with an accessible region in a 138-donor bulk-liver
chromatin-accessibility QTL (caQTL) map (Currin et al., 2025). At least 45 of
these 191 were already reported by the source. In 108 Catalog genes, the same region also
colocalizes with the gene's eQTL, and in 100 when the accessibility-QTL and
eQTL lead variants must have r² ≥ 0.8. These are annotation fields on
gene records. They create no regulatory-element record and no
credible-set-to-gene link. Constraint, essentiality and secretion fields are
given for follow-up design and do not enter any class.

We illustrate the record with four worked entries: GLP1R, THRB, HKDC1 and MTARC1
(Figure 6C). GLP1R is untested because its RNA was not detected. Its genetic call is
indeterminate, with no Tier-1/2 SuSiE result, and it routes first to the missing
decisive assay. THRB has genetic
support (SuSiE PP.H4 > 0.999, UK Biobank GGT) and an indeterminate RNA call,
because its change (log2FC -0.20) is below the 0.5 floor. HKDC1 is supported on
both axes (SuSiE PP.H4 0.988, UK Biobank AST; log2FC +1.07, BH *q* = 8.0 ×
10^-7). MTARC1 remains unresolved on these expression-based axes. The printed
protein-coding routes propose context and function measurements. They do not
convert genetic support into causality, and they do not rank the entries.

The recurring examples of the Resource show further follow-up experiments. GNMT
and CYP2C19 are lower in disease at RNA and protein levels and are
established-state entries, which the Catalog routes to context-matched
perturbation and rescue. MAT1A is also lower at both levels, but its pooled RNA
effect (log2FC -0.485) is below the 0.50 floor of a canonical call. Its Catalog
class is therefore unresolved, not established-state. COLOC retains GNMT's single-signal direct-NAFLD ABF result (PP.H4 0.706), which rests on a sub-threshold GWAS association (p = 2.6 × 10⁻⁴) and does not extend to MAT1A or CYP2C19. GNMT has no multi-signal SuSiE call (best PP.H4 0.128), so its
Catalog genetics call is indeterminate and it takes the same perturbation and
rescue route. Its local regulatory haplotype is unresolved. The ECM/IGFBP7 and
ductular-injury/BICC1 entries separate three routes: ambient-resistant lineage
measurement, matched-null tissue organization and functional perturbation. With
dataset-specific variances modelled, neither program passes the 117-program BH
and stability rule for a pooled stage association. This analysis cannot exclude
the ECM/IGFBP7 effect size estimated in the earlier stage analysis. The original
single-cell atlas lacks ambient-RNA correction. An ambient-RNA sensitivity
analysis corrected five datasets. Two datasets (GSE189600 and GSE192742) could
not be corrected, and 11 T-cell programs were untestable. It does not establish cell of origin.
Within the sampled sections, the ECM/IGFBP7 score is more spatially
autocorrelated than matched control genes. This is not a donor-level result,
and most fixed programs (68 and 74 of 114 per source) exceed the same null, so
it does not show program-specific organization. Ductular-injury/BICC1 is
spatially indeterminate. Ductular-injury/BICC1 therefore routes first to a better-powered
donor-resolved tissue test, and both programs need an independent cell-of-origin
assay. PNPLA3 and MTARC1 remain boundary controls for coding or protein-level
mechanisms that expression-QTL analysis need not recover. Each carries a
fine-mapped missense variant (maximum PIP 0.706 for PNPLA3 and 0.757 for
MTARC1), and this leaves its eQTL-colocalization class unchanged. These examples
show how routing works. They are not a ranking of the strongest genes.

The routing layer assigns the smallest experiment that could distinguish the
remaining explanations (Figure 6D). Of the 28,233 genes, 13,551 are untested,
13,127 unresolved, 1,129 established-state associated, 400 genetically anchored
and 26 supported on both axes. Genetically anchored regulatory candidates route
to allele-aware editing with target-expression and locus-linked phenotypes.
Established-state candidates route to context-matched perturbation and rescue.
Candidates supported on both axes need both relationships connected in the same
model. Untested candidates (13,551) route to the missing decisive assay and
unresolved candidates (13,127) to measurement of the most decisive missing layer;
for these 26,678 genes the route names the missing layer, not a specific assay. Transcript
depletion is only one possible perturbation. Increasing gene activity, testing
regulatory-element changes, measuring protein function, lineage restriction and
context-specific challenge are separate designs that no current record
receives. For noncoding nominations,
the route starts from the molecular object. A disease-associated lncRNA routes
to RNA depletion with measurement of the neighboring gene. A colocalized or
ambiguous lncRNA locus needs paired RNA and locus perturbations to distinguish
DNA-element, transcriptional and mature-RNA effects. The Catalog records genes, not
regulatory DNA elements, and does not link credible sets to genes. Genetically
anchored genes route to allele-aware regulatory editing at the gene level, but no
record names the element to edit, and the element-level route defined in Methods
is assigned to no record.

The Catalog is deposited as machine-readable tables (Data availability;
[accession to be supplied]). Each displayed entry links to its source graph, its
complete tested universe, a data dictionary and its machine-readable
next-experiment rule. A reader can therefore see, for any gene, which layers
were tested, which were untestable and which experiment is proposed.

Across Figures 1 to 7, the Resource records inherited genetic evidence and
disease-state remodeling for each gene, and places both in cell and tissue
context where an assay can observe them. The Gene Catalog is the reusable
output. For each nominated molecular object, it records the biology, the limit
of the evidence and the next experiment.

## Discussion

In this work, we present a human MASLD Resource that records, for each gene,
liver-trait genetic evidence and cross-sectional stage-associated disease-state
remodeling, and keeps the assay boundaries of the two separate. The bulk
disease-state layer comes from 774 genetic-identity units in five cohorts. These
units are conservative clusters, not verified independent donors. The genetic
layer comes from eQTL colocalization across liver-trait GWAS. Among the 3,769
genes that had stable Ensembl IDs, were tested in bulk RNA and had a
multi-signal colocalization result, 426 were genetically supported and 220 were canonical
disease-state calls, and 26 were in both, against 21.21 expected at the same
expression level. We detected neither enrichment nor depletion. Because the
matching did not fully control expression, we state no direction
(Supplementary Note 3). A further 935 disease-state calls had no multi-signal
colocalization result, and these are not genetic negatives. The size of this
overlap does not establish biological separation, equivalence or a privileged
intersection of the two maps. It is also not evidence that the underlying
processes are independent. The signed comparison is bounded too. Genome-wide
rank agreement between ALT TWAS z and the disease contrast lies inside ±0.03
without composition adjustment, but not with hepatocyte fraction in the design,
and we state no direction (Supplementary Note 5).

Most colocalized genes take their strongest evidence from liver-enzyme traits,
largely because far more gene-trait pairs were tested in that class. The
strongest native colocalization evidence comes from liver-enzyme traits for 423
of 436 genes and from direct disease or liver-fat traits for 13. Only 22 genes
have any direct-trait colocalization, and far more enzyme than direct-trait
gene-by-GWAS pairs were tested (12,121 against 771). Pooled across ancestries,
protein-altering plus canonical-splice variants carry 29.65% of fine-mapped
posterior probability for direct MASLD and liver-fat traits and 6.67% for enzyme
traits, counting credible sets at genome-wide-significant leads. This contrast is confounded with GWAS ancestry and fine-mapping
coverage, and no ancestry stratum supports a calibrated trait-class call
(Supplementary Note 4). A comparison of the trait classes within one ancestry
needs fine-mapping of all UKBB enzyme lead windows and of the MVP European ALT,
AST and NAFLD GWAS. Some protein-altering variants may act through protein
function without changing steady-state RNA. Others change RNA abundance through
nonsense-mediated decay or splicing (Rivas et al., 2015). One hypothesis for
which genes appear in both sets is therefore that expression-based
colocalization misses risk variants that act through protein function rather
than RNA abundance. The present comparison does not test it, because it does not
hold gene eligibility, phenotype, study power and expression observability
fixed. These results set one requirement for gene lists built this way in MASLD:
the trait that supplies each gene's evidence should be stated. A gene nominated
by ALT is not the same claim as a gene nominated by liver fat.

Matching on expression did not recover the prespecified prediction that
genetically supported genes are more hepatocyte-specific than disease-state
genes. In 281 expression-matched gene pairs, the observed difference was in the
opposite direction (Supplementary Note 3). Every test of that opposite direction
is post hoc, and none is calibrated. The difference comes mainly from the 56
pairs whose disease-state partner was itself tested by colocalization. We
therefore read it as a descriptive difference conditioned on selection, not as
confirmatory class separation or independent biological replication. The
protein-response comparison of 73 pairs gives no support for a class difference,
and it does not establish equivalence. Expression remains a major observability
difference between the classes. Among genes with a single-cell hepatocyte-specificity
value (399 of the 400 genetic-only and 1,124 of the 1,129 state-only genes),
mean bulk expression of genetic-only genes exceeds that of state-only genes by
3.80 log2 units, and 281 of the 399 genetic-only genes obtain an
expression-matched partner. This 3.80 log2 gap belongs to the
hepatocyte-specificity comparison. Within the expression deciles of the
3,769-gene overlap test, colocalized genes differ from other genes by 0.07 log2
units on average (0.70 before stratification). Unequal observability could affect which genes
appear in both sets, but the present comparison does not isolate that
explanation. Germline regulatory effects may act before biopsy, in restricted
cell states or through small persistent shifts. Disease-state expression can
reflect injury, adaptation or composition. Neither class is automatically causal
or irrelevant.

The recurring examples make this distinction concrete. GNMT has a single-signal
direct-NAFLD colocalization (ABF PP.H4 0.706) but no multi-signal call. With
MAT1A and CYP2C19, it is lower in pooled RNA, in cross-sectional fibrosis
contrasts and in adjusted liver protein. In an exploratory analysis of 46 MASLD
participants, the three proteins covary. The covariation remains after
adjustment for fibrosis stage but not after adjustment for NAS or for all four
histology components (BH over 12 pair-by-adjustment tests). This pattern is
compatible with a shared association with histologic severity rather than a
three-gene circuit, but adjustment does not establish mediation. The ECM/IGFBP7 and ductular-injury/BICC1
programs show a different limit. The original atlas lacks ambient-RNA
correction, the ambient-RNA sensitivity analysis cannot resolve all datasets or
programs, and of the two only ECM/IGFBP7 exceeds a matched-gene spatial null,
which most of the 114 scorable fixed programs also exceed. Cross-assay
evidence therefore strengthens specific claims, but disease association, tissue
architecture, lineage origin and inherited mechanism remain distinct.

Each component of this Resource has precedent. Staged bulk liver transcriptomes
have been published (Hoang et al., 2019; Govaere et al., 2020; Pantano et al.,
2021), meta-analyzed with gene-level GWAS signal (Piras and DiStefano, 2024),
organized into network modules (Kumak et al., 2026) and ordered along a molecular continuum
(Kamzolas et al., 2026; doi:10.1038/s42255-026-01543-7). Single-cell and
spatial atlases resolve disease cell states (Ramachandran et al., 2019;
Guilliams et al., 2022; Boesch et al., 2026). A single-cell multiome assigns MASH
risk loci to liver cell types (Elison et al., 2025), single-cell eQTL mapping has been
done in MASLD liver biopsies (Hong et al., 2025), chromatin,
reporter-assay and CRISPR interference data assign target genes to MASLD risk
variants (Zhu et al., 2026), and large GWAS have combined
expression and coding evidence to nominate genes (Sveinbjornsson et al., 2022).
Multi-ancestry genetic association data have been used to nominate 212
candidate MASLD genes and, through genetically predicted expression, to select
drug-repurposing targets among them (Seagle et al., 2025).
This Resource adds a per-gene Catalog that keeps colocalization results beside
disease-state calls and records whether each layer was untested or tested
without support. Our reanalysis of a 138-donor liver caQTL map (Currin et al.,
2025) places a colocalized accessible region at 191 of 581 European liver-trait
GWAS signals that colocalized with a liver eQTL. At least 45 of these were
already reported by the source.

A related distinction between evidence roles has been reported at the protein layer. A
liver proteogenomic study of 504 biopsies found that differential protein
abundance and protein-QTL colocalization nominate largely different gene sets
(Gobeil et al., 2026b). Of the
five genes that study nominated, our map places GCKR on the primary
expression-colocalization side, MTTP colocalizes only in the
single-signal (ABF) sensitivity analysis (best PP.H4 0.82), ERLIN1 does not
reach that set (best ABF PP.H4 0.47; Catalog class unresolved), and LGALS1 is in the established-state class. MTARC1 lies outside our
expression map but is resolved by the protein-QTL analysis of that study. GCKR
has differing expression-QTL and protein-QTL directions. These cases show the
difference between protein and expression observability. They do not establish
a binary split into causal and reactive genes. A single-cell multiomic liver
study likewise reports that chromatin-accessibility QTL and expression QTL
colocalize with largely distinct risk loci. It supplies cell-type regulatory
context that tissue-level eQTL analysis cannot provide (Elison et al., 2025).
These external studies extend the assay boundary that the Gene Catalog records.
They are credited source comparisons, not results of this paper.

The gene and module results of Kumak et al. are likewise source findings, not
results of this paper. Of their NAS-associated and fibrosis-associated genes,
944 and 1,013 map to the Catalog, and 258 and 301 of these reach the canonical
bulk call here. Of those present in both the bulk and colocalization tables, 661
of 902 and 714 of 943 had no multi-signal colocalization result. The Catalog
keeps this state separate from a tested negative. Among the 241 NAS and 228
fibrosis genes that were evaluated, 36 and 28 colocalized, against 28.9 and 24.9
expected for genes of the same expression level (Mantel-Haenszel OR 1.32, 95% CI
0.91 to 1.92, and OR 1.16, 95% CI 0.76 to 1.76; Holm-adjusted *p* 0.29 and 0.51,
Holm over the two lists). Thus, among genes evaluated by multi-signal COLOC, the
colocalization rate of these genes did not differ detectably from that of genes
of the same expression level (power in Supplementary Note 15). Kumak et al.
prioritize genes by combining meta-analysis evidence with centrality in a
curated protein interaction network. The 117 programs and 43 systems used here
were fixed before any outcome was read and are never reweighted by network
position. A sparsely annotated gene therefore keeps its assigned weight, and an
inadequately observed gene is recorded as untestable rather than ranked low.

The stage-stratified analyses add biological context without establishing
change over time. Adjacent-stage contrasts and cell-resolved programs describe
stage-ordered remodeling across different participants. The externally anchored
molecular continuum, adapted from Kamzolas et al. (2026), adds a separate
molecular view without redefining disease. It aligns strongly with the
fibrosis-based pathway, NMF and Hotspot results, which is consistent with the
major biological programs being preserved. This alignment was computed on the
same data and is not independent validation. Stage-adjusted associations show
that molecular variation within histologic bins tracks the continuum, and which
programs carry that variation depends on adjustment for cell-type composition.
Histology therefore remains the interpretable clinical anchor, and the continuum
adds molecular resolution within a stage. All 117 Hotspot program definitions
and 43 system memberships remain fixed. An outcome-blind map extends this
comparison across complete pathway and Hotspot layers. Stage and continuum
effects align broadly, but only a minority of systems keep one sliding-window
direction across all five cohorts. This supports shared but cohort-variable
remodeling, not a universal molecular sequence or ontology. Program-pair
covariance is reported in Results. It could not be tested in GSE213621, where
the program-group permutation control failed. Bulk covariance also mixes cell
abundance with program activity, so it does not show coordination between cell
types. Neither analysis identifies transitions within individuals, a dominant
fibrosis switch or therapeutic windows. Those claims need repeated sampling of
the same participants and perturbation experiments.

For users of the Resource, the main lesson of Figure 7 is that each model score
must be read beside its exposure audit and its task-native baseline.
Our per-cohort exposure audit showed that all 109 GSE268273
validation samples were in the BulkFormer pretraining data, so that model could not
be evaluated on samples it had not seen. Against task-native baselines, the paired
difference between the severity instrument and the in-cohort ridge has a 95%
interval that includes zero, and GSE268273 was also used to choose the
expression transform. A per-region model of H3K27ac that adds the expression of genes near each
region (the cis head) scored higher than a model on global expression components
alone, but its gain over a reduced-rank RNA model like that of Figure 7B was
smaller (Supplementary Note 14).

There are several limitations to this study. The liver-eQTL panel is
reused and of European ancestry, so non-European GWAS lack ancestry-matched
molecular QTL data. The regulatory-genetic layer is based on expression QTL, so
protein-altering and protein-level mechanisms that leave the transcript
unchanged, such as MTARC1, fall outside its resolution. This layer is also
tissue-level, not cell-type-resolved. Source-reported protein-QTL and
single-cell chromatin studies supply distinct assay contexts, but no single
layer is a complete causal filter. Most genetically supported genes come from
liver-enzyme rather than direct MASLD traits. Partitioned polygenic scores
already show that MASLD risk variants do not form one class (Jamialahmadi et
al., 2024). The direct-versus-enzyme split used here is a different, trait-based
partition. The class supported on both axes is small, and every context assay
has source-specific coverage. No accepted spatial analysis provides a
MASLD-versus-healthy disease-effect estimand, and this paper contains no new
experimental perturbation. The GTEx comparison measures histologic fat of
unrecorded cause in postmortem liver. Its slide tissue-area negative control
missed its prediction (Spearman *rho* 0.354, *p* = 0.0564) and came close to its
stop rule. In a post hoc check, adding tissue area to the model left the fat
correlation at 0.631 (0.638 without it). Ischemic time was also correlated with the same genes' disease
effects, with the opposite sign and less strongly than fat (*rho* −0.350 against
0.636). The caQTL reanalysis rests on one 138-donor map. Its chr1 and chr19
positive controls were seen in a pilot and are not blind, and at least 45 of the
191 colocalized signals were already reported by its authors. We did not check
overlap with published three-way eQTL, caQTL and GWAS colocalizations
(Broadaway et al., 2024). We record these as Catalog limitations and do not
infer mechanisms from source-owned data.

Further limitations concern units of analysis and models. Genotypes called from RNA-seq
reads identified repeated libraries and conflicting labels, but these calls,
sample-key checks and expression-based QC do not give complete genotype-based
identity confirmation, which was not established here. The bulk units therefore remain genetic-identity
clusters, not verified donors. Gene-level totals can hide opposing isoform
changes. Differential gene expression concerns total gene abundance,
differential transcript expression concerns individual transcript abundance, and
differential transcript usage concerns relative isoform proportions within a
gene. This Resource tests the first, not isoform switching. Short reads shared
among isoforms add assignment uncertainty, and transcript estimates depend on
the annotation. Oakes and Pachter (2026) emphasized these distinctions and the need to
separate demonstrated responses from exploratory physiological explanations. Model validation is specific to its task, cohort and
training exposure, and rank correlation is not absolute-stage calibration. The
GSE193066 repeat-biopsy intraclass correlation (ICC) and SD 0.6433 use
source-declared pairs whose identity was not qualified. They cannot establish
primary within-person reliability or label-error calibration, and ceiling values
derived from them remain sensitivity scenarios. Serial biopsies can include
shared technical effects and true change, and uncertainty conditional on fitted
predictions omits retraining variation. The comparator models were fitted
within each test cohort, so they had already seen that cohort's data.
Comparison with them therefore cannot show that the models perform equally on
unseen cohorts, or that one model class is better in general.

Two questions stay open. We do not assign genes or programs separately to steatosis, inflammation or fibrosis; fibrosis stage is the severity anchor, and NAS enters only as a sum or as an adjustment set. Colocalization here tests single variants, so the Resource does not test whether combinations of variants act differently, and the local haplotype at GNMT is unresolved. Both questions are recorded as indeterminate, not negative.

Within these limits, the Resource gives each gene an experiment matched to its evidence. Genetically anchored candidates motivate allele-to-gene and
context-appropriate perturbation without requiring disease-state expression.
Disease-state-associated candidates motivate tests of active tissue function and
reversibility without being labeled inherited drivers. Candidates supported on
both axes motivate experiments that connect the two relationships in one model.
Untested and unresolved candidates, 26,678 of the 28,233 genes, first need the
missing decisive measurement. The Gene Catalog and the source tables needed to check each call
are deposited with this paper (Data availability), so that others can look up
these distinctions for any gene instead of relying on a single aggregate rank.
