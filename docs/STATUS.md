# Current paper and codebase status

Navigation updated November 9, 2026. This page summarizes the latest recorded
manuscript state; it is not a live SLURM report. Exact values and their
qualifications remain in [RESULTS](RESULTS.md), with scientific scope in
[PAPER](PAPER.md) and authorized next work in [ROADMAP](ROADMAP.md).

## Executive assessment

The working manuscript and supplement have editable LaTeX sources and compiled
PDFs. Figures 1–7 are embedded as complete plates in the main Results. The
selected genetic map remains provisional after the variant-identity and LD
checks, and full unique-cell verification is incomplete. Formatting completion
does not make the manuscript submission-ready or adopt new scientific results.

## What is official, updated, or open

| Area | Latest recorded state | Detail |
|---|---|---|
| Manuscript | Working main paper and supplement compiled November 8–9; citations and unassigned panels remain open | [Manuscript workspace](manuscript/README.md) |
| Named publication bundle | `2026-10-21-r9` remains a named historical bundle; later working changes are separate and do not overwrite it | [Release history](technical/RELEASE_PROVENANCE.md) |
| Bulk and stage | Validated update selected for the paper: adopted 774 genetic-identity units, not verified independent donors; transportability qualifications are incorporated | [Bulk results](RESULTS.md#human-bulk-rna-seq) · [Transportability](RESULTS.md#variant-identity-and-transportability-2026-11-07) |
| Genetic map | Historical selected r2 map is provisional; duplicated eQTL identities and invalid LD require qualified full fits before new dependent memberships or claims | [Current genetic qualifications](RESULTS.md#variant-identity-and-transportability-2026-11-07) |
| Single-cell counts and programs | Full unique-cell verification remains incomplete. All 117 program memberships and weights remain fixed; verified count files alone do not qualify donor effects or native Hotspot scores | [Physical-capture correction](RESULTS.md#physical-capture-correction-2026-11-07-in-progress) |
| Gene Catalog | Eligibility correction adopted November 6; RNA effects and native genetic positives retained, with later genetic-identity qualifications still applicable | [Selected Catalog](../RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/index.html) · [Eligibility results](RESULTS.md#catalog-eligibility-adoption-2026-11-06) |
| SteatoSITE | Prespecified stage comparison adopted October 30; extension B/C adopted November 3, with eQTL direction retained as exploratory | [SS results](RESULTS.md#steatosite-held-out-replication-ss-series-adopted-2026-10-30-user-approval) · [SSX results](RESULTS.md#steatosite-extension-analyses-ssx-series-b-and-c-adopted-2026-11-03-user-approval) |
| Molecular context and model evaluation | Working figure-specific decisions and their assay/source limits apply; candidate existence alone does not select a result | [Figures 1–7](../figures/main/INDEX.md) · [Model registry](../Analysis/MASLD_Model_Benchmark/MODELS.md) |
| Experimental follow-up | Separate future paper; its outcomes cannot redefine Resource evidence, figures, thresholds or claims | [Two-paper firewall](PAPER.md#two-paper-firewall) |

## What is still missing

1. Qualified genetic identities, adequate LD and full fits before revising
   genetic membership or dependent overlap, accessibility and Catalog claims.
2. Complete unique-cell checks and qualified donor effects; historical graph/PCA
   Hotspot measurements remain conditional on their original inputs.
3. Remaining references, unassigned panels, scientific reconciliation and the
   complete [final release checklist](ROADMAP.md#final-release-checklist).
4. Explicit limits where evidence remains unavailable: GSE193066 identity,
   HMSMA donor/histology/metabolite joins, independent protein confirmation and
   sufficient verified RNA–chromatin donor pairs. Future evidence is distinct
   from unfinished analyses required for this manuscript.

## Figure readiness

The [main figure index](../figures/main/INDEX.md) records the working panel
selections for all seven figures. The [manuscript workspace](manuscript/README.md)
links the current embedded plates and PDFs. The
[figure directory](../figures/README.md) separates main artwork, supplementary
material and alternatives without changing any panel or source path.

## Active compute

Query SLURM for live jobs. Dated job IDs and completed attempts are retained
with each analysis and in the [earlier status record](archive/folder_navigation_2026-11-09/STATUS.md#active-compute).
The November unique-capture work below is a recorded scientific state, not a
claim that a particular job is currently running. See the
[hypothesis directory index](../Analysis/Hypotheses/README.md) for its files.

## Work that can proceed now

Follow the [current roadmap](ROADMAP.md), starting with its
[genetic-identity follow-up](ROADMAP.md#variant-identity-and-transportability-follow-up-2026-11-07).
Use the [analysis routes](technical/DATASETS_AND_PIPELINES.md) to locate the
existing scripts. Earlier completed-work lists are retained in the snapshot;
they are not instructions to repeat fits or submit new work.

## Known interpretation boundaries

- Bulk genetic-identity clusters are not verified independent donors. Cells,
  capture runs, sections and spots do not supply independent biological replication.
- The 1,854 pooled bulk calls remain discovery associations; native-cohort
  direction agreement does not establish independent disease/control replication.
  The prespecified SteatoSITE stage comparison remains a separate contrast.
- Historical genetic counts, overlap and accessibility links retain the current
  variant-identity and LD qualifications. Subtraction does not establish a
  validated unaffected gene set.
- Fixed-weight expression projections do not replace historical Hotspot
  neighborhood/PCA measurements. The 117 program definitions and weights stay fixed.
- Cross-sectional stage and continuum associations do not establish temporal
  change, a disease mechanism or a reliable patient-level performance ceiling.
- Retain assay observability, source dependence and the distinction between
  negative tests, unavailable evidence and non-applicability in every interpretation.

## Latest recorded updates

The following November 6–9 entries are retained from the previous status page.
They describe their stated dates and scope.

**Review-round figure restructure installed in the manuscript, 2026-11-10 (candidate, not adopted).**
Per the approved plan. Figure 1 is the user's plate (`figure1_JL.pdf`); old 1B/1C are S1A/S1B and the
overlap result is one paragraph (detail in Supplementary Note 3). Figure 2 adds enzyme-only vs
disease-liability (2F), out-of-sample check (2G) and TNFSF10 (2H); old 2D is S2X/S2Y. Figure 3 adds the
SteatoSITE scatter (3F), composition adjustment 1.463 → 0.256 (3G) and the within-hepatocyte partition
(3H); old 3E strip and 3G are S3Y/S3Z. Figure 4 now shows SteatoSITE programs (98/112, 4D), the trio per
hepatocyte (4E), the CCL20+ NF-kB state (4F), hepatocyte share of CCL20/FCAMR/TRIM31 (4G) and PNPLA3 148M
(4H); old 4C and 4E–4H are S4AB–S4AF and their text is Supplementary Note 20. Figure 6 shows known genes
(6B), the 29 both-supported genes (6C/6D), clinical targets (6E), recovered therapeutic axes (6F) and
drug-target calibration (6G) from the candidate Catalog
`RNA-seq/results/evidence_passports/candidates/gene-catalog-conf1-20261110T023225Z/`; old 6B/6D/6C are
S6A–S6C. Figure 7 is S7. Unique-cell recount: 650,083 cells, 0/117 programs (q 0.232, ≤54 donors); main
text still prints 1,232,318 with a note. eQTL duplicate-ID refit (90 genes) applied to the manuscript
on user approval 2026-11-10 (`Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/RESULT.md`): 83 stay
positive, CDK5RAP3 lost, IRF5/CTSS/LRRC37A/KANSL1-AS1/ABHD12/MPRIP untestable; primary 436 → 429;
anchored 400 → 393 (CONF1 rule 420 → 413); overlap 3,583 genes, 419/209/26, 21.0 expected, MH OR 1.317
(0.848–2.044), exact p 0.223. Catalog `candidates/gene-catalog-conf1-eqtlfix-20261110T130857Z/`; panels
S1A, S3Z, S6A, S6B redrawn (`figures/candidates/eqtl-refit-panels-20261110T131833Z/`). Descriptive
breakdowns of the 436 (Fig 2A, 2F/2G, S1B, S2G, p12 sensitivity) are labeled r2 and not recomputed. Build
`CONF1-20261109T221846Z/latex_build_23151594`: 0 errors. Backups in `CONF1.../latex_backup_*`.

**Review-round confirmatory tests, 2026-11-09 (not adopted).** Six prespecified tests in
`Analysis/Hypotheses/CONF1-20261109T221846Z/`; numbers in
[RESULTS](RESULTS.md#review-round-confirmatory-tests-conf1-2026-11-09-not-adopted). Positive: the
ambient-corrected within-hepatocyte partition (1,155 of 7,214 replicated stage genes) and the
candidate Catalog fields (coding fine-mapping; disease liability). Not confirmed or not replicated:
ballooning link of the inflamed-hepatocyte score; per-hepatocyte GNMT/MAT1A/CYP2C19 in the
replication cohorts. Manuscript text updated in `fig3.tex` (trio per-hepatocyte paragraph, wrong
Figure 3E pointer) and Supplementary Note 19; `catalog.tex` now says Figures 1 to 6. GSE289173
hepatocytes opened with user approval (C5, job 23144496): the trio is lower per hepatocyte (Holm ≤
5e-5 each; controls are a separate series). The inflamed-hepatocyte state is not testable there
(C1b, CCL20 detection gate failed). Candidate Fig 2 and Fig 3 panels in
`CONF1-20261109T221846Z/figures_candidate_20261109T183739/`, not placed in `figures/main`.

**Main figure layout corrected, 2026-11-09.** Each of Figures 1–7 now keeps all
of its panels together as one complete plate on one page within the Results.
Full captions follow the artwork. Source panels, scientific text and existing
qualifications are preserved.

**LaTeX manuscript and Overleaf, 2026-11-08.** The working manuscript uses the
supplied reference paper's single-column format and Times-style typography.
Main figures and their complete captions appear within the Results;
supplementary notes, legends and artwork compile separately. The Summary and
Introduction explain the biological question, observations and use without
metric inventories. Both documents compile without missing glyphs or overfull
boxes, and all figure PDFs match their current source files. Direct Overleaf
Git access is authenticated; manuscript locations and synchronization are
recorded in [the manuscript workspace](manuscript/README.md). This formatting
change preserves the scientific qualifications and adopts no new results.

**Manuscript reconciliation, 2026-11-08.** Working Introduction, overview,
Catalog Results, Discussion, Methods and legends now carry the genetic-identity
and physical-capture qualifications already recorded below. RORA and FABP1
historical H4 values are explicitly not valid-LD controls. No new genetic
membership or unique-cell effect is adopted. Full unique-cell verification
remains incomplete; completed count files alone do not establish corrected
donor effects or native Hotspot scores.

The linked Google manuscript's main text and section outline now match the
current working draft, with original author details and reference entries
preserved. Exact text readback and the Resource scope check pass. Citation
placeholders and unassigned panels remain; this update does not make the
manuscript submission-ready. See [the manuscript workspace](manuscript/README.md).

**Genetic identity and replication interpretation, 2026-11-07.** The all-fit
screen is complete: 3,633/18,889 saved eQTL fits have duplicate variant IDs,
including 90/436 primary positive genes and 137/708 full-portfolio positive
genes. MPRIP and DIP2B have duplicate positions inside a supported eQTL credible
set; other affected fits still need refitting because duplicated regressors can
alter conditional signals. The current genetic map is explicitly provisional.
MPRIP's original MVP allele labels and frequency disagree, and both rounded
and original unrounded PolyFun LD fail the regularized matrix check. New fitting
entrypoints refuse ambiguous identities and invalid LD; the historical named
outputs are preserved. Exact sensitivities and remaining source requirements are
in [RESULTS](RESULTS.md#variant-identity-and-transportability-2026-11-07).

The completed transportability results are now incorporated into the abstract,
Results, Discussion, Methods, Supplementary Note 19 and figure legends. Two
checked supplementary PDFs show native-cohort direction agreement and the
three specified tissue examples with REML uncertainty. The 1,854 pooled calls
are unchanged; independent disease/control replication is not established by
this sensitivity. SteatoSITE's original stage comparison remains separate.

The genotype-derived MPRIP sensitivity is complete. All three GWAS fits converge
on valid LD; both eQTL fits exhaust 500 iterations without convergence, including
whole-position exclusion. No qualified MPRIP H4 is obtained. Four bounded jobs
completed successfully; no further genetic job is pending from this task.

**Unique-nucleus correction in progress, 2026-11-07 (user request).** Source
experiment identifiers and deposited filenames resolve 275 atlas run records
into 145 physical captures: 271 RNA runs and four antibody-tag runs excluded
from RNA processing. There are 24 repeated RNA captures. The 18 GSE244832 joint
recounts are reused; six additional recounts cover three GSE185477 and three
GSE202379 captures. The 117 programs and their stored membership weights remain
fixed. Corrected count, donor-lineage and fixed-weight projection outputs are
being verified under `Analysis/Hypotheses/NUCLEUS-FIX-20261107/`; historical
neighborhood-dependent Hotspot scores and embeddings are not thereby corrected.
The leftover cohort-disagreement attribution is removed. The active histology
error calculator is now `Analysis/MASLD_Model_Benchmark/kleiner_sensitivity.py`:
it requires an explicit assumption, has no empirical error default, and passed
formula, invalid-input and no-default checks. Named historical packages remain
intact and their `kleiner.py` is retired in the model registry.

**Catalog eligibility adopted, 2026-11-06 (explicit user request).** The selected
Catalog is `RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/`.
It contains 28,233 genes: 400 genetically anchored, 26 supported in both domains,
1,118 established-state associated, 12,958 unresolved and 13,731 untested.
Exactly 188 genetic states become untestable and 180 classes change; eight of
those states were already in untested classes. All RNA measurements/calls,
436 native genetic positives, 117 programs and 6,633 memberships are preserved.
The Catalog has 14,502 jointly eligible identities; the separate native-ID
overlap comparison has 3,589 eligible genes, 426 genetic and 209 RNA calls,
with 26 shared (MH OR 1.2815, 95% CI 0.8262–1.9878, exact p=0.27699).
Historical r8/r9 copies are unchanged. Detailed adoption checks and jobs are
recorded under `Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/`.

**Highest-priority research started, 2026-11-06; new results are candidates.**
The five-cohort random-effects sensitivity completed in job `22811901` after
all five source fits reproduced. All 24,115 models were estimable, but none
passed full-family BH with modified Hartung–Knapp intervals and four/five
cohorts. This does not change the fixed 1,854 pooled calls or establish absent
effects: 1,850/1,854 random-effects means retain the pooled sign and 1,208 have
the same sign in every available cohort; the minimum full-family q is 0.17143.
It limits claims about transfer across cohorts. Median I-squared among
the fixed calls is 60.35%; 908/1,854 nominal heterogeneity p values are below
0.05 (descriptive, not corrected discoveries).

The saved genetic inventory completed in `22811613` and the original-reference
three-fit pilot in `22811698`. All three posterior values reproduce, but the
MPRIP input has two opposite-allele indel records assigned the same position ID
and incompatible LD entries. Original fit reproduction therefore does not
validate MPRIP's genetic support. Its variant identity and source-matched LD
remain unresolved. Guarded pilot `22812217` completed: MPRIP is explicitly
untestable under the position-ID representation, and RORA/FABP1 reproduce
exactly. No variants were arbitrarily deleted and no additional genetic
membership change is adopted. Source-matched LD and allele resolution remain
the next requirements; all submitted jobs for this bounded task are complete.

**Confirmed review corrections completed, 2026-11-06 (user request).**
The corrected donor-lineage filter gives 71 hepatocyte donor identifiers
(13 references, 58 cases), 29,948 tested genes and 90 BH calls; the P30-excluded
fit gives 90 calls in 70 identifiers. Nineteen of 27 former coarse-stage
comparisons remain estimable; eight are explicitly untestable. All-cell effects
are reproduced and its hepatocyte comparison is updated. Main Figure 4 program
models and their 0/117 result are unchanged.

The Xenium correction repeats all label-dependent selection in 48,620
permutations and 10,000 donor bootstraps. SCGN remains the only Holm call among
the 12 colocalized genes. Adjusted exon/intron inference gives 2,711 calls and
114 highlighted genes; its recomputed global permutation p is 0.004975. The
corrected local Catalog preserves all 28,233 named identities/classes and fixes
4,218 RNA applicability states. The eligibility candidate made the 188
genetic-state and 180 class changes implied by the existing priority rule
explicit; these are now adopted above.
LD-sensitive gene qualifications and separate SuSiE/ABF values are provided.

Working Figure 3 is assembled from the current individual panels; its exon/intron
supplement and Figure 7E are corrected and visually checked. Code fixtures,
independent roster/BH calculations, saved-estimate reconstruction, independent
Xenium resampling checks and the Resource scope check pass. Detailed commands,
failed attempts and final checks are in the
[technical correction record](../Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/REPORT.md).
New research in the [analysis plan](ROADMAP.md#scientific-analysis-plan-2026-11-06)
was subsequently authorized; current run details are recorded with the Catalog
adoption above. Remaining scientific limits include LD qualification,
independent protein/lineage validation and orthogonal composition evidence.
Named release r9 is preserved; final scientific adoption remains separate.

## Earlier status and execution records

The [full preceding status page](archive/folder_navigation_2026-11-09/STATUS.md)
retains per-analysis details, dated jobs and earlier figure inventories. Its
older states remain subject to the current qualifications above and in RESULTS.
Historical labels **Corrected rerun PROMOTED to canonical 2026-08-17** and
**Five-cohort synchronized candidate** describe earlier steps, not a new
scientific decision or the current genetic qualification.

### Historical missing-work inventory, superseded 2026-10-14

See the [dated inventory](archive/folder_navigation_2026-11-09/STATUS.md#historical-missing-work-inventory-superseded-2026-10-14).
Use the current missing-work list and ROADMAP rather than this earlier queue.

## Canonical DEG threshold change (2026-08-12)

The current bulk definition is BH-adjusted `padj < 0.05` with
`|log2FC| > 0.50`. The interval-null result at `lfc = 0.25` is a distinct
sensitivity. [RESULTS](RESULTS.md#human-bulk-rna-seq) records the release-specific
counts; the [threshold history](archive/folder_navigation_2026-11-09/STATUS.md#canonical-deg-threshold-change-2026-08-12)
retains earlier comparisons and does not select a current overlap estimate.

## Historical defect inventory and scientific dispositions

See the [preserved inventory](archive/folder_navigation_2026-11-09/STATUS.md#historical-defect-inventory-and-scientific-dispositions).
Its dated findings and completed corrections do not certify later figures.

## Source of truth

| Question | Authority |
|---|---|
| Scientific identity and claims | [PAPER](PAPER.md) |
| Exact numbers and qualifications | [RESULTS](RESULTS.md) |
| Next work and final checklist | [ROADMAP](ROADMAP.md) |
| Current panels | [Main figure index](../figures/main/INDEX.md) |
| Detailed history | [Release history](technical/RELEASE_PROVENANCE.md) · [Archive index](archive/INDEX.md) |
