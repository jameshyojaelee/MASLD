# Figure 2 Results — inherited genetic variation and candidate molecular links

> Scope: phenotype- and ancestry-aware regulatory evidence, not a causal-gene
> catalogue. Exact COLOC classes, counts, memberships, and overlap values remain
> withheld until the corrected 50-study by 22-chromosome rerun is audited and
> adopted; no downstream Cas13-library membership belongs in this text.

## Trait-associated inherited variation defines a genetic map distinct from established disease-state transcription

We assembled a prespecified inherited-variation portfolio spanning direct disease,
MRI-derived liver-fat, and liver-enzyme GWAS across multiple ancestry groups
(Figure 2A). Direct disease and imaging traits remain distinct from ALT, AST,
and GGT proxy phenotypes because enzyme associations report informative liver
biology but are not equivalent to direct MASLD susceptibility. The final
study, stratum, and ancestry counts will be populated only from the adopted
corrected registry.

Colocalization against the Broadaway liver-eQTL panel produces two differently
observable maps (Figure 2B). This figure applies the Figure 1 decision grammar:
phenotype applicability, eQTL and LD coverage, method eligibility, evidence
state, unresolved alternative, and next discriminating experiment.
Single-variant ABF and multi-signal SuSiE results
have distinct eligibility requirements; SuSiE requires converged eQTL and GWAS
fits plus usable LD. Missing high-resolution output is therefore an
observability state, not a negative result or a lost modality vote. Figure 2 will report complete
method-specific denominators and single-signal-only states from the adopted release.

We describe the complete result as a map of inherited genetic variation
relevant to MASLD and liver traits and reserve direct MASLD language for the
explicitly tiered disease subset. Regulatory language is restricted to
supported GWAS/eQTL relationships; protein-altering and unresolved noncoding
signals retain their own labels. Because the eQTL side reuses the Broadaway
panel, the added value is the MASLD/liver-trait tiering, multi-ancestry GWAS
integration, fine-mapped DNA context, and comparison against diseased human
tissue rather than a new molecular resource. Historical source-overlap and
credible-set counts remain supplementary until they are rebuilt against the
post-rerun registry; no July count is adopted by appearing in this working
draft (Broadaway et al., 2024).

Ancestry labels describe the GWAS populations, not ancestry-matched molecular
regulation in every analysis. Non-European GWAS paired with a European liver
eQTL resource test cross-population sharing rather than ancestry-matched liver
regulation. We report explicit LD and eQTL provenance rather than treating
these as equivalent ancestry replications (Figure 2E–F).

The genetics–expression interface will be recomputed from the corrected COLOC
release and the synchronized fragment-count bulk substrate (Figure 1B). Until
then, exact memberships, overlap fractions, enrichment tests, and genome-wide
PP.H4–expression associations remain withheld rather than inherited from the
July release.

The current official release provides a useful boundary check on the recurring
disease-remodeling examples. GNMT, MAT1A, CYP2C19, IGFBP7, and BICC1 have no
primary SuSiE PP.H4 above 0.5. GNMT alone has a single-signal ABF PP.H4 above
0.5 in one direct-NAFLD stratum, anchored at GRCh37 6:42928758; its primary
multi-signal result remains below threshold. The strongest stored ABF anchors
for MAT1A, CYP2C19, IGFBP7, and BICC1 are 10:81143387, 10:96306509,
4:57949517, and 10:60583832, respectively, but their gene–trait PP.H4 values
remain below 0.5. These coordinates are COLOC shared-variant anchors, not
causal variants or GWAS credible-set PIP. Exact results will be replaced, not
silently inherited, once the corrected COLOC release is adopted.

GNMT is the provisional worked bridge rather than evidence for a genetically
regulated three-gene program. At the current-official direct-NAFLD ABF anchor
rs2296805, the lower-GNMT-expression allele aligns with higher disease risk.
rs2296805 tags an unresolved local haplotype: rs2296804 has UKBB-European
r2=1.000, while rs9471976 and rs11752813 have r2=0.819 and 0.580 with the
anchor. The latter variants also align lower liver GNMT expression with higher
direct-NAFLD risk and were associated with hepatic GNMT protein and reporter
activity in a source study of 268 human liver biopsies. These observations
motivate allele-aware locus resolution, but they do not identify a causal
variant and will enter the final figure only if corrected COLOC retains the
exact GNMT–direct-NAFLD pair (Ji et al., 2012).

The distinction between inherited genetic evidence and established tissue
state has two implications. First, differential expression cannot serve as a
required validation filter for genetically anchored regulatory candidates.
Second, the absence of genetic support does not render a disease-state program
unimportant, because common-variant colocalization has limited power and
interrogates inherited regulation rather than tissue activity. We consequently
define genetic-only, disease-state-only, and convergent classes rather than
treating either map as a noisy approximation of the other.

Protein-altering and protein-level loci mark the outer boundary of
expression-based colocalization. PNPLA3 I148M and external MTARC1 proteogenomic
evidence illustrate mechanisms that an expression-mediated model need not
recover. Their final Resource classes will be assigned only from the adopted
corrected release (Gobeil et al., 2026).

We therefore quantify DNA architecture directly from all reliable Tier-1/2 SuSiE credible-set instances rather than restricting the analysis to credible sets that can be reconstructed from a colocalized gene's top SNP. In the validated candidate, 6,415 of 6,467 study-by-locus credible-set instances pass the prespecified PIP-sum criterion. Complete GRCh37 reference-oriented VEP annotation assigns 100% of their normalized posterior inclusion probability. After grouping repeatedly recovered sentinels within 1 Mb and weighting each genomic region equally, direct MASLD/liver-fat loci allocate 22.4% of that probability to protein-altering variants, 0.3% to canonical splice variants, 4.7% to synonymous or UTR variants, and 72.6% to other noncoding variants. Liver-enzyme loci allocate 3.4%, 0.0%, 3.6%, and 93.0%, respectively (Figure 2D). The grouping limits repeated-study dominance and does not assert LD or statistical independence. Credible-set-weighted estimates remain a sensitivity analysis in the source table.

Taking protein-altering and canonical-splice variants together, the difference between the two trait classes is 22.5% against 3.4%, a 6.5-fold ratio that survives permutation testing and Benjamini-Hochberg correction across the six-test family (q = 6.0e-4) and is stable to dropping any single study, dropping any single locus group, changing the grouping window between 0.25 and 5 Mb, and replacing the distance window with ancestry-matched LD blocks (ratio 4.99 to 9.58 across all 62 refits, every permutation P below 8.0e-4). Dropping the study that contributes 95 of the 208 direct credible-set instances raises the ratio rather than removing it.

This architecture has a direct consequence for how such loci can be mapped to genes. A variant that alters protein sequence acts on protein function rather than on transcript abundance, so it need not appear as an expression QTL and cannot be recovered by expression-based colocalization. The colocalized map reflects this: the strongest colocalization evidence for 448 of 462 genes comes from a liver-enzyme trait (ALT 150, AST 164, GGT 134) and for only 14 from a direct disease or liver-fat trait (NAFLD 8, PDFF 5, NASH 1). Because the table records one row per gene, this describes which trait class supplies each gene's best evidence rather than excluding weaker evidence from the other class. Expression-based colocalization applied to liver traits therefore resolves predominantly liver-enzyme regulatory biology, and the most disease-proximal part of inherited MASLD risk is the part it is least able to see.

Within the noncoding posterior probability, promoter, liver ABC-enhancer, and liver-lineage open-chromatin annotations are retained as overlapping contexts rather than probabilities that must sum to one. The same candidate also makes missing context visible: 52.1% of direct-trait and 60.4% of liver-enzyme noncoding posterior probability lacks resolved deposited regulatory context. Consequently, the DNA result supports a broad noncoding architecture and an assay-coverage map, not a claim that every noncoding variant has a known regulatory mechanism. Colocalized target-transcript biotype remains a separate COLOC-dependent analysis.

An external massively parallel reporter assay provides variant-level context
without contributing to our map's discovery. Any locus-level agreement or
observability gap will be linked to the final corrected COLOC membership rather
than copied from the historical map (Zhu et al., 2026).

Together, these analyses define the genetic arm as a bounded map of
trait-associated inherited variation, fine-mapped molecular context, and
colocalized liver-trait regulatory candidates. It complements, but does not
substitute for, the established disease-state transcriptome.
