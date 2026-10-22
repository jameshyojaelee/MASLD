# MASLD Resource paper: documentation guide

**Current corrections, 2026-10-15.** The paper's working figures now use the
adopted 774-unit bulk fit and verified COLOC r2 (436 primary, 996 broader
genetic candidates). Figure 3 includes all five per-cohort contrasts and
54-gene/216-effect stage artwork; Figure 4 includes the atlas schematic,
the three-cohort 364-unit Figure 4E arm (review round 3) and Figure 4G/H on
the 367 disease units, with no outline in 4H. Figure 5B/C use
774-unit RNA effects. All five individual Figure 6 panels use the selected
r8 Catalog, whose local review, molecular-object companion and native-size
review assembly are complete. Figure 5's vector assembly uses current panels;
all 66 fixed twelve-trait genetic-correlation pairs now have checked estimates.
Original-source diagnosis correction also refits the fixed 117-program
stage models. With per-dataset variances on the 60-donor ambient-correctable
universe (H59), zero programs are selected; the 62-donor ductular-injury/BICC1
call (parametric q 0.038, permutation q 0.123) is a sensitivity only.
The two GSE174748 disease explants are cirrhotic, not simple steatosis;
the 60-donor ambient comparison retains unchanged weights and eligibility rules.
The corrected N4 null is untestable for different-visit pairs; historical
87 random-null and 69 stage-BH calls are withdrawn. Exact values and remaining
claim limits are in [RESULTS](RESULTS.md) and [STATUS](STATUS.md).
Historical paragraphs below retain earlier inspection states; they do not
describe current missing panels. The complete publication bundle still follows
the [ROADMAP checklist](ROADMAP.md#final-release-checklist).

**Figure 2 coordinate correction adopted, 2026-10-13.** The coding/splice
posterior ratio is 6.026671 (95% cluster-bootstrap interval 2.574593–13.510165),
with 233 physical clusters and 57 robustness comparisons. This value is pooled
across ancestries and confounded with GWAS ancestry and fine-mapping coverage
(H58, 2026-10-15): the coding/splice posterior fraction is 29.88% for direct-trait
vs 6.92% for liver-enzyme traits within EUR and 6.73% vs 4.54% within AMR. These
are descriptive estimates, with no trait-class call. Current estimates and
the scope of this adoption are in
[RESULTS](RESULTS.md#adopted-corrected-coordinate-results-2026-10-13).
Figure 2D and its unassembled noncoding-context companion have been regenerated
and installed with corrected source tables and captions.
This October 13 decision concerned Figure 2 coordinates. COLOC r2 is now
selected on October 14; overall publication-bundle approval remains separate.

**Historical Figure 4 working artwork, 2026-10-05, superseded for E–H:** the reviewed compact B–H PDFs were
at the indexed individual-panel paths. D covers all 117 discovery/content
lineage assignments; E covers all 113 eligible weighted linear bulk associations
with historical donor-identity qualification. C retains 24 modules. This does
not adopt the synchronized scientific release.

**Historical Figure 3F artwork, 2026-10-05, superseded October 14:** the previously missing adjacent-stage matrix
is now built at `figures/main/fig3_bulk_transcriptomics/panels/fig3f_stage_deg_genetics_matrix.pdf`.
It shows 61 genes from the completed September genetics candidate and historical
stage fits. The PIP column was removed at the user's request. This display does not resolve
counting, donor identity, or synchronized scientific adoption.

**Historical Figures 3–4 identity qualification, 2026-09-29:** the genotype ledger identifies
repeated bulk libraries and conflicting source labels; GSE193066 identity calls
remain unresolved. Historical bulk roster counts and newly reconstructed fits
are conditional, not verified independent-donor inference. The correction
candidate includes an explicit roster join and a separate provisional identity
sensitivity. Do not adopt its bulk/network significance claims until identity
is resolved; see [the disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).

The September Figures 3–4 correction candidates separated exact from grouped histology,
cohort replication from pooled evidence, and metabolic abundance from function.
See [the scientific dispositions and supporting analyses](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).
The 25 corrected PDFs were installed at working manuscript paths on 2026-09-30
at the user's request. Bulk/network inference remains conditional on identity
resolution; this artwork update does not adopt the synchronized scientific
release. The panel sequence and fixed program definitions are unchanged.

The historical 2026-10-05 correction work excluded unresolved identity/source-label records
from primary inference and requires full ALT/scaffold annotation correspondence
before the corrected bulk fit. The identity candidates have been independently
rechecked, and both six-panel Figure 7 sets exist for review. The annotation
candidate now has 855 unresolved correspondences after official-source checks.
The authorized reference policy excludes 63 contamination/false-duplication
annotations while retaining their records; TRPM3 remains included. A revised
fixed-basis NMF projection has been prepared and tested on historical counts;
the primary bulk refit
and its dependent results remain on hold. Candidate outputs are under
`figures/candidates/resource-corrections-20261005-v1/`; completed changes,
checks and unresolved limits are logged in the
[existing correction outcomes](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#resource-corrections-implementation-2026-10-05).
No final bulk count or synchronized release is adopted by this work.

This repository supports a standalone, human-centered MASLD Resource and
analytical-framework paper. Start here; do not reconstruct the project from old
plans or dated progress logs.

The Resource analyzes existing human datasets to map inherited-risk candidates
and disease-associated liver biology. Genetic association and colocalization
nominate candidates; they do not alone establish that a gene causes MASLD. The
separate Cas13 screen is outside this Resource's analysis and evidence claims.

## Read these five documents

| Document | Answers | Read when |
|---|---|---|
| [PAPER.md](PAPER.md) | What is the paper, what is novel, and what may we claim? | Deciding analyses, figures, framing, or journal strategy |
| [STATUS.md](STATUS.md) | What is official, what has been validated, and what remains open? | Starting any work session |
| [RESULTS.md](RESULTS.md) | Which numbers may be used now, and from which release? | Writing text, legends, slides, or reviews |
| [ROADMAP.md](ROADMAP.md) | What should happen next, in what order, and what is out of scope? | Planning or coordinating work |
| [manuscript/README.md](manuscript/README.md) | Where are the current manuscript drafts? | Editing the paper |

The order of authority is: **PAPER → STATUS → RESULTS → ROADMAP → technical
records → archived history**. A dated progress note never overrides a current
summary or a read-only release.

## The paper in one paragraph

**Title:** *A multimodal human MASLD Resource records liver-trait genetic evidence and
disease-state liver expression side by side for each gene.* (Proposed 2026-10-15, pending the user's approval. The
2026-09-07 coordination clause was removed under the joint H52/H54 rule, and
"resolves the architecture" was dropped after H58.) The Resource records trait-associated
germline variation, multicohort cross-sectional disease remodeling, and
donor-aware cell and tissue context side by side. It distinguishes coding, regulatory, and
unresolved noncoding genetic signals; localizes remodeling across cell
programs, protein, chromatin, and spatial assays; and delivers the results in a
MASLD Gene Catalog with appropriate follow-up experiments. Assay observability,
biological-unit control, provenance, and machine-readable evidence states make
these maps interpretable, but they are supporting rigor rather than the paper's
headline. A provisional direct-NAFLD GNMT locus provides the bridge from
inherited variation to a severity-linked GNMT–MAT1A–CYP2C19 RNA–stage–protein
state, without treating the three genes as one inherited mechanism. The
ECM/IGFBP7 and ductular-injury/BICC1 programs then separate disease association,
matched-null tissue organization, and lineage attribution. The aligned
117-program ambient candidate decreases the ECM/IGFBP7 estimate while changing BICC1 little,
although neither paired effect change passes family correction. Same-atlas
donor-paired contrasts descriptively localize ECM/IGFBP7 to fibroblasts and
BICC1 to cholangiocytes, while only ECM/IGFBP7 exceeds matched spatial
background. Following Kamzolas et al. (2026;
doi:10.1038/s42255-026-01543-7), a validated histology-anchored molecular
continuum candidate adds stage-adjusted donor-level resolution without
replacing recorded fibrosis stage. Its sliding windows are fixed descriptive
summaries, not inferential or temporal bins. The exact GNMT direct-NAFLD ABF pair (PP.H4 0.7058, 2021_34841290_NAFLD_EUR) is in COLOC r2; its local haplotype remains unresolved.
The paper must be scientifically useful without results from the separate
follow-up study, which may prospectively test Resource nominations but cannot
redefine its classes, figures, or claims.

The compact aging-inspired Figure 4–6 rebuild is a validated, visually reviewed
candidate under `figures/candidates/aging-inspired-resource-rebuild-candidate-2026-08-15-v8/`.
It has not been adopted into `figures/main`.

## Folder map

The user-requested [integrated MASLD regulatory model and database plan](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md)
connects inherited molecular effects, cellular context and disease-associated
remodeling. It reconciles prior evidence and retains T01–T24 and their parameter
menus across six campaigns: accessibility prediction, context, assay transfer,
remodeling, indels/targets/haplotypes, and generalization. Each task defines its
question, evidence, comparison, evaluation and completion criterion.
The plan includes public-data expansion, a reproducible public core with
separately identified restricted research models, and the authorized initial-week coordinator
execution. Rigorous, useful MASLD adaptations of methods from other diseases
also qualify as potential contributions. The expanded
[literature review](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#22-literature-informed-expansion-and-experiments-worth-testing)
adds 30 annotated references and 16 bounded comparisons with concrete variables,
controls and stopping rules. New experiments are not dependencies.
The initial one-week campaign is now authorized and running under the
[execution envelope](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#194-authorized-initial-week-2026-09-15).
Its code is in `scripts/analysis/alphagenome_campaign/`; candidate outputs and
actual compute accounting are in `GWAS/finemapping/results/alphagenome_campaign/week1-20260915/`.
Protected outcomes remain closed. This adds no adopted result; earlier executed
candidate findings remain in the linked
[evidence table](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#32-integrated-evidence-and-next-action-dispositions)
and [program results](technical/ALPHAGENOME_PROGRAM_RESULTS_2026-09-14.md).

The current bounded two-model candidate asks whether liver RNA recovers measured
H3K27ac differences beyond recorded histology and whether a short liver allele
adapter adds information beyond native 1-Mb sequence context. Its fixed design,
source eligibility and development results are in the [integrated program](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#23-two-useful-models-and-a-compact-figure-7-candidate)
and [RESULTS](RESULTS.md#two-model-molecular-value-development-2026-09-22-not-adopted).
The completed development comparisons retain global-plus-local RNA and
calibrated native 1 Mb under the fixed complexity margins. Candidate scoring
interfaces and reproducible examples are described in the
[analysis README](../scripts/analysis/masld_two_models/README.md).
No candidate here changes a working panel or the synchronized release.

The implemented continuum translation extension is documented in the
[technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md)
and [Methods](manuscript/METHODS.md#continuum-extension-intervention-hypotheses-and-measurable-disease-states).
It connects inherited direction, within-stage program variation, assay context
and readouts as a candidate analysis; it does not replace the adopted evidence
classes or establish treatment-response prediction.

| Location | Contents | Authority |
|---|---|---|
| `docs/manuscript/` | Draft text, legends, and release notes | Writing surface; numbers come from RESULTS |
| `docs/technical/` | Dataset, pipeline, release, number-history, and figure indexes | Technical reference |
| `docs/literature/` | Competitor and literature map | Scientific context, not result authority |
| `docs/reference/` | Source notes and imported references | Supporting material |
| `docs/competitor_analysis/` | Detailed paper-by-paper analyses | Supporting material |
| `docs/archive/` | Superseded plans, progress, audits, and verbatim originals | Provenance only |
| `figures/main/INDEX.md` | Current seven-figure directories, panel requirements, candidate states, and what has to be true before a panel is adopted | Figure-routing authority |

The paper uses Gene Catalog r8: 28,233 genes from the adopted 774-unit bulk
fit and COLOC r2. Its unchanged scientific tables now have a completed local
review record in the [read-only reviewed copy](../figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/PASS06_VALIDATED).
The [updated local interface](../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/access_v2/index.html)
and [28,233-record molecular-object companion](../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/index.html)
show the current unit and held-annotation restrictions. These local outputs
do not publish an external portal. Exact classes and distinct overlap
denominators are in [RESULTS](RESULTS.md#catalog-review-and-molecular-object-companion-2026-10-14).

## Plain-language status labels

| Label | Meaning |
|---|---|
| **Current official release** | Read-only release currently allowed to support exact manuscript claims |
| **Validated update awaiting release** | Checked analysis intended for the paper, but downstream numbers and output files have not yet been adopted together |
| **Continuum working figures adopted** | The Kamzolas-adapted co-primary scores and complete molecular families met their candidate criteria; hash-pinned Figure 4F–4G and S3/S4 panels are in the working figure authorities, while Figure 3 and the synchronized release remain unchanged |
| **Validated spatial candidate** | Spatial registry, coverage, effects, portal tables, release-linked Figure 1/Figure 6 panels, and bounded diagnostics met their criteria; Figure 5A–F and its vector assembly are adopted working outputs, while the synchronized Resource release remains pending |
| **In progress** | Analysis or compute is still running; memberships and exact counts can change |
| **Supplementary boundary** | Valid stress test that did not support a main-paper upgrade |
| **Could not be tested** | Public data lacked an authoritative biological-unit or assay join |
| **Historical** | Preserved rationale or provenance; not an active instruction |
| **Aspect-specific** | An observation attributed to metabolism, inflammation, or fibrosis individually. Not available in the current release |
| **Aspect indeterminate** | The design resolves a summary severity score only, so which aspect an observation belongs to is unresolved rather than negative |
| **Variant combination** | Two or more variants at a locus considered jointly. Distinct from a single colocalized variant, and untested by the colocalization design |

Internal labels such as `F_five`, `BG-001`, candidate release IDs, job IDs, and
hashes remain in technical records and the archive. Human-facing documents use
“updated five-cohort bulk analysis,” “bulk-counting decision,” and “validated
update awaiting release.”

**A NAS sum is not aspect resolution.** NAS is the sum of steatosis, ballooning,
and lobular inflammation. A result that tracks the sum cannot say which of the
three it belongs to. Say "severity-associated" and record the aspect as
indeterminate; do not say "inflammation-associated" without a component-level
analysis.

## Common tasks

| If you need to… | Use |
|---|---|
| cite a sample, gene, DEG, program, or overlap count | [RESULTS.md](RESULTS.md) |
| understand why two numbers differ | [technical/NUMBERS_HISTORY.md](technical/NUMBERS_HISTORY.md) |
| check a live dependency or decide what can run now | [STATUS.md](STATUS.md) and [ROADMAP.md](ROADMAP.md) |
| edit the story or a figure section | [PAPER.md](PAPER.md), then [manuscript/README.md](manuscript/README.md) |
| locate a current main panel or generator | [`figures/main/INDEX.md`](../figures/main/INDEX.md), then the figure-specific `CANDIDATE_PANEL_INDEX.tsv` |
| inspect pipeline or dataset ownership | [technical/DATASETS_AND_PIPELINES.md](technical/DATASETS_AND_PIPELINES.md) |
| inspect the noncoding upgrade | [STATUS.md](STATUS.md) for the current checks, [RESULTS.md](RESULTS.md) for validated candidate counts, `scripts/manuscript/noncoding_resource/` for the producers, and `Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-2026-08-11/` for the human-readable Catalog rules |
| inspect or hand off the validated spatial candidate | Start with the concise [spatial omics handoff](technical/SPATIAL_OMICS_HANDOFF_2026-08-12.md), then inspect `Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/`, `spatial-diagnostic-bundle-candidate-2026-08-11/`, and the result-led `spatial-impact-figures-candidate-2026-08-12/`; canonical Figure 5F now pairs the accepted matched-gene calibration with one outcome-blind registered Guilliams et al. H&E overlay, and accepted results and layout decisions are summarized in [RESULTS.md](RESULTS.md) |
| inspect the validated single-cell ATAC candidate | Use `Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/`; `NON_GENETIC_READY`, `GENETIC_READY`, and `FULL_READY` mark the non-genetic and adopted-COLOC replay outputs as complete. Corrected program SE/p/q and primary/legacy peak SE tables are in the separate `atac-context-v3-statistics-correction-20261022/` candidate beside it; see [RESULTS.md](RESULTS.md#validated-single-cell-atac-context-v3-candidate-not-adopted). Canonical Figure 5D uses the completed primary-pair context without adopting the entire ATAC candidate package |
| inspect the histology-anchored molecular continuum | Start with `scripts/analysis/histology_anchored_continuum/README.md`. `hac-molecular-layers-20260818T173348Z` retains transcript/pathway/program inference; `hac-system-inference-20260824T232135Z` retains the participant-level 43-system models and permutations. The corrected atlas and pooled candidates retain Figure 4G/S4Q–S4U values. The active display revision uses `histology-continuum-molecular-systems-caption-free-20260825T130145Z` and `histology-continuum-molecular-systems-s4u-caption-free-20260825T130145Z`: embedded explanatory prose is removed from Figure 4G/S4Q–S4U and retained in the external caption authority. S4U keeps all 43 systems, nine descriptive windows, and 18 strict maxT asterisks; the late-minus-early contrast remains the retained display-only row-order statistic but is not displayed, and the 33 primary BH calls remain in the source registry. Mirrored source tables, hashes, and audits are in the Figure 4 authority. |
| audit a release or recover an old plan | [technical/RELEASE_PROVENANCE.md](technical/RELEASE_PROVENANCE.md) and [archive/INDEX.md](archive/INDEX.md) |
| compare recent competitors | [literature/REVIEW.md](literature/REVIEW.md) |

## Safety rules

- Never mix numbers from the current official release and the validated update
  in one claim without labeling both.
- Describe stage contrasts only as cross-sectional stage-associated remodeling.
- Treat the molecular continuum as complementary to recorded histology. Fixed
  sliding windows are descriptive only and cannot select features, define
  significance, orient an axis, or imply within-person temporal ordering.
- Cross-cohort plots use Healthy, Steatosis, Steatohepatitis, and Cirrhosis;
  retain MASL/MASH/NASH only for source-specific provenance or contrasts, and
  never treat generic MASLD as a stage.
- COLOC identifies shared-signal candidates; it does not prove a causal gene.
- Cells, spots, arrays, sections, guides, or sequencing runs are not biological
  donors.
- In spatial results, zero or missing coverage is not negative evidence;
  distinguish `untestable`, `indeterminate`, and `source_dependent` states.
- Do not use future follow-up study outcomes to select Resource genes, classes, figures,
  thresholds, or claims.
- Do not reopen broad public-data mechanism hunting as a submission dependency.

## What happened to the older documents?

The former current-state, scope, narrative, outline, number ledger, plan pack,
and progress logs were retired from the active documentation tree on
2026-08-11. Their verbatim, checksummed originals remain recoverable through
[archive/INDEX.md](archive/INDEX.md) and the migration manifest. Historical
files are provenance, not current instructions.
