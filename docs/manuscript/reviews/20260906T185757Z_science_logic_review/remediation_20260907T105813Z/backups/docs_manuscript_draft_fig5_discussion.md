## The MASLD Gene Catalog connects gene classes to discriminating experiments

> Working Figure 6 Results. The signed Gene Catalog bundle is a validated candidate,
> not an adopted manuscript release. The final panel roster and exact counts
> must come from Plan 60. Historical convergence, drug-calibration, and receptor
> follow-up panels are supplementary and do not define Figure 6.

The final Resource does not end with a ranked list. For every eligible gene and
program, a Gene Catalog entry records the primary gene class, phenotype
provenance, assay and biological unit, testing universe, source dependence,
direction where interpretable, machine-readable evidence state, state reason,
the unresolved alternative, and deterministic next-experiment rule (Figure 6A).
Calls distinguish supported, discordant, indeterminate, untestable, not
applicable, source-dependent, and prespecified adequate-negative states. A
nonsignificant result is never converted into a negative call merely because it
is convenient for the narrative.

Figure 6B makes observability visible. A blank evidence cell can mean that the
gene was not detected, the source was underpowered, the assay cannot observe
the proposed mechanism, the dataset lacks an authoritative biological-unit
join, or the experiment simply has not been performed. These states have
different scientific consequences and remain distinct in the downloadable
tables and source graph. In the current signed candidate, zero rows satisfy the
strict `tested_negative` rule; that absence is a valid Resource result, not a
software failure.

Worked Catalog entries show how the same Resource produces different next steps rather
than different degrees of confidence (Figure 6C). GNMT, MAT1A, and CYP2C19
carry concordant disease-state RNA, cross-sectional stage, and tissue-protein
evidence, routing the hepatocyte metabolic-function hypothesis to
context-matched perturbation and rescue. If retained by corrected COLOC, GNMT
also routes to allele-aware resolution of its local regulatory haplotype; this
does not extend to MAT1A or CYP2C19. The ECM/IGFBP7 and
ductular-injury/BICC1 entries separate three routes: ambient-resistant lineage
measurement, matched-null tissue organization, and functional perturbation.
The current source-diagnosis refit supports both hepatocyte-labelled programs,
but the original atlas lacks ambient correction and the first ambient audit is
not synchronized to that refit. ECM/IGFBP7 is spatially supported, whereas
BICC1 is spatially indeterminate. BICC1 therefore routes first to a
better powered donor-resolved tissue test, while both require an independent
cell-of-origin assay. PNPLA3 and MTARC1 remain boundary controls for coding or protein-level
mechanisms that expression-QTL analysis need not recover. These entries are
routing demonstrations, not a curated greatest-hits ranking.

The routing layer assigns the smallest experiment capable of distinguishing the
remaining explanations (Figure 6D). Genetically anchored regulatory candidates
route to allele-aware editing with target-expression and locus-linked
phenotypes; established-state candidates route to context-matched perturbation
and rescue; candidates supported on both axes require both relationships to be
connected in the same model; and unobservable candidates first require the
missing decisive assay. Transcript depletion is only one possible perturbation.
For noncoding nominations, the plotted route begins from the molecular object:
a regulatory DNA element routes to allele-aware element testing, a
disease-associated lncRNA routes to RNA depletion with neighboring-gene
measurement, and a colocalized or ambiguous lncRNA locus requires paired RNA
and locus perturbations to distinguish DNA-element, transcriptional, and mature
RNA effects.

Figure 6E closes the Resource loop by linking each displayed Catalog entry to its
source graph, complete tested universe, data dictionary, and machine-readable
next-experiment rule. The historical evidence-weighted heuristic is retained as
an optional backward-compatible annotation, not a posterior, a default sort, or
a target-validity claim. Its failed matched clinical benchmark and the terminal
public perturbation stress tests remain visible in the supplement as
falsification boundaries. The publication-facing portal must default to the
Gene Catalog and gene-class view before Plan 60 can accept it.

Across Figures 1–6, the Resource connects inherited genetic variation to
multicellular disease remodeling and then places both in cell and tissue
context. The Gene Catalog is the reusable endpoint: it preserves the biology,
evidence boundary, and next experiment for each nominated molecular object.

## Discussion

This Resource maps two connected dimensions of MASLD biology: trait-associated
inherited genetic variation and multicellular remodeling in established
disease. In the current official release, the resulting gene maps are largely
non-overlapping among jointly testable genes. Exact updated memberships and
overlap values await adoption of the corrected COLOC portfolio and
synchronized bulk release. Their separation is biologically informative, but
it does not imply that the maps are independent, opposed, or that their
intersection is a privileged truth set.

The separation is mechanistically plausible. Germline regulatory effects can act before biopsy, in restricted cell states, under exposures not represented by the sampled tissue, or through small persistent shifts that do not meet a disease-state DEG threshold. Conversely, strong disease-state expression may reflect tissue injury, adaptation, fibrosis, immune recruitment, or altered cellular composition without detectable common-variant regulatory support. Neither class is automatically causal or irrelevant.

The recurring examples make that distinction biological rather than abstract.
GNMT is a provisional direct-disease locus that intersects a coordinated
GNMT–MAT1A–CYP2C19 loss across pooled RNA, cross-sectional fibrosis contrasts,
and adjusted liver protein. Protein covariation and its attenuation after
histology adjustment identify a severity-linked state rather than a
three-gene circuit. IGFBP7 and BICC1 then expose a different boundary: the
original atlas lacks ambient-RNA correction, and the first sensitivity is not
synchronized to the current stage metadata, while only ECM/IGFBP7 exceeds
matched spatial background. Cross-assay evidence
therefore strengthens specific claims, but disease association, tissue
architecture, lineage origin, and inherited mechanism remain distinct.

The same evidence-role separation is reported at the protein layer. A liver proteogenomic study of 504 biopsies found that differential protein abundance and protein-QTL colocalization nominate largely different sets (Gobeil et al., 2026). Cross-referencing its five source-reported effectors with our fixed map captures GCKR on the primary expression-colocalization side, ERLIN1 and MTTP in the sensitivity union, and LGALS1 in the established-state class. MTARC1 remains outside our expression map but is resolved by the source's protein-QTL analysis, while GCKR has differing expression- and protein-QTL directions. These are examples of protein-versus-expression observability, not proof of a causal-versus-reactive binary. A separate single-cell multiomic liver study likewise reports that chromatin-accessibility QTL and expression QTL colocalize with largely distinct risk loci and supplies cell-type regulatory context unavailable to tissue-level eQTL analysis (Elison et al., 2025). These external studies therefore extend the assay boundary represented by the Gene Catalog; they are credited source comparisons, not results generated by this paper.

Assay-native tissue measurements reinforce this interpretation. Disease-state-only and convergent genes are enriched for proteomic and single-cell changes. Genetic-only genes are also enriched at the protein layer, but far more weakly, with an adjusted odds ratio of 2.1 against 9.1 for disease-state-only genes, and they show no detectable enrichment in the spatial or single-cell endpoints. The separation is therefore a gradient in tissue-state expression rather than a categorical absence of it. The convergent class often has the highest raw positive rate but is small and is not detectably superior to disease-state-only genes. Thus, the Resource supports complementary coverage rather than a monotonic evidence hierarchy.

The corrected gsMap analysis further limits the inference. Native spatial concordance identifies trait-associated tissue patterns, but it does not show that direct MASLD risk localizes preferentially to convergent genes after expression and spatial-specificity matching. The earlier exclusivity claim was produced by broader, partially expression-defined sets and an invalid PCC arm label. Spatial evidence remains useful for localization, but it does not convert the intersection into a privileged truth set.

The strict drug benchmark reaches a similar conclusion. Broad MASLD-associated target annotations are enriched near the top of the heuristic score, but the weighted score does not outperform modality count or genetics in a matched clinical endpoint. Sparse approved-target coverage and the absence of a curated intervention-level preclinical stage prevent a valid development-gradient test. Clinical examples therefore illustrate observability and evidence-role boundaries; they do not calibrate the score as a predictor of success.

The stage-stratified analyses add biological context without establishing
temporal change. Adjacent-stage contrasts and cell-resolved programs describe
stage-ordered remodeling across different patients. The externally anchored
molecular continuum, adapted from Kamzolas et al. (2026;
doi:10.1038/s42255-026-01543-7), adds a complementary view rather than a
competing disease definition. Its strong alignment with fibrosis-based pathway, NMF, and Hotspot
results shows that the major biological programs are preserved, while
stage-adjusted associations show that histologic bins compress coordinated
molecular variation. Histology therefore remains the interpretable clinical
anchor and the continuum supplies within-stage molecular resolution. An
outcome-blind membership map extends this comparison across complete pathway
and Hotspot layers: stage and continuum effects align broadly, but only a
minority of systems retain one sliding-window direction across all five
cohorts. This argues for coordinated yet cohort-variable remodeling rather than
a universal molecular sequence or ontology. Neither
analysis identifies transitions within individuals, a dominant fibrosis
switch, or therapeutic windows. Longitudinal sampling and perturbational
experiments are required for those claims.

Several limitations define the scope of the Resource. The liver-eQTL panel is reused and predominantly European; some non-European GWAS therefore lack ancestry-matched molecular QTL data. The regulatory-genetic layer is expression-QTL based, so protein-altering and protein-level mechanisms that leave the transcript unchanged, such as MTARC1, fall outside its resolution. It is also tissue-level rather than cell-type-resolved; source-reported protein-QTL and single-cell chromatin studies supply complementary contexts without making any single layer a complete causal filter. Most genetically supported genes derive from liver-enzyme rather than direct MASLD traits. The convergent class is small, and every context assay has source-specific coverage. The credible-set linkage available from aggregated COLOC outputs covers only 75 of 473 primary genes, precluding a global coding-architecture conclusion. No accepted spatial branch provides a MASLD-versus-healthy disease-effect estimand, and this paper contains no new experimental perturbation. Those are explicit Catalog limitations rather than reasons to manufacture a mechanism from source-owned data.

Within these boundaries, the Resource provides a practical experiment-routing framework. Genetically anchored candidates motivate allele-to-gene and context-appropriate perturbation without requiring disease-state expression. Disease-state-associated candidates motivate tests of active tissue function and reversibility without being labeled inherited drivers. Candidates supported on both axes motivate experiments that connect the relationships in one model. Untestable or source-dependent candidates first require the missing decisive measurement. Preserving these distinctions—and releasing the provenance needed to audit them—is more useful than interpreting a single aggregate rank as target truth.
