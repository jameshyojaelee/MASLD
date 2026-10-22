# Agent-context wording review — 2026-09-23

This is a maintenance record, not a paper scope, result, or release authority.
The five current documents in `docs/` remain authoritative.

## Outcome

Simplified wording in 29 existing files during the broader review. References
to the separate experimental study now use its role rather than repeatedly
naming its platform. The human Resource remains independent of its outcomes.
The older RNA-seq README now identifies its February 2026 summaries as historical
and links to the current scope, status, results, and separate-study documentation.
No dataset, named release, analysis result, figure, guide sequence, or parameter
value was rewritten. Existing unrelated changes were preserved.

[Unverified] No phrase in this review is an established classifier trigger.
No requests were submitted to Claude to test classification, and these edits
provide no evidence that future flags will stop.

## Search coverage and limits

The initial bounded text scan covered documentation, scripts, configuration,
agent-local text and the library-design directory, to four directory levels.
It read 13,101 files, including duplicate files in existing agent worktrees.
After excluding agent-local copies, the inventory contained 602 matching lines
in 159 main-tree files. These are lexical matches, not 602 safety findings;
some are paths, identifiers, ordinary statistical language or legitimate methods.

A supplementary `git grep` examined tracked Python, R, Markdown, shell, YAML
and TOML sources across the repository for narrower experimental phrases.
It found nine matching lines before the final two comment edits. Raw data,
generated results, environments, third-party source, historical archives,
release packages and other worktrees were not edited. The scan was bounded,
not an exhaustive inspection of every file or nested experiment directory.

## Changes made

- Rephrased incidental platform references as the separate experimental study
  or follow-up screen in the five authorities, active manuscript instructions,
  figure guidance, selected model handoffs and implementation notes.
- Retained explicit exclusions on using future screen results to choose Resource
  genes, classes, thresholds, figures, or claims.
- Kept mouse context and RNA-reduction meaning in study-specific documentation.
- Clarified the historical RNA-seq introduction and restored a working link from
  the separate-study README to the current paper scope.
- Updated the configuration's mirrored goal prose and the scope checker's exact
  wording marker together. Configuration keys and enforcement logic are intact.
- Limited source edits to docstrings/comments, apart from that wording marker.
  The Catalog recommendation text edited in the preceding task remains in place;
  this turn changed only its module/function documentation.

## Precise terms retained

| Location | Reason to retain the wording |
|---|---|
| `Analysis/MASLD_Model_Benchmark/config/datasets/gse275483.toml` | The dataset title identifies the source experiment; an in vivo study is not interchangeable with a cultured-cell experiment. |
| `Analysis/downstream_analysis/essentiality/scripts/extended_essentiality_analysis.py` and `run_essentiality_pipeline.py` | The platform names distinguish the published Liang data from DepMap data. |
| `RNA-seq/Mouse/Western_Diet_Datasets/GSE246088/scripts/run_de_analysis.R` | The note distinguishes the animal intervention from whole-liver RNA sampling. |
| `scripts/portal/test_gene_catalog_v2.py` and `gene_catalog_v2.py` | Original prohibited field names, input examples, matching rules and error codes are needed to detect screen-derived evidence. |
| `Cas13_Library_Design/scripts/guides/` | Model names, score fields, reference sequences, file paths and specificity criteria are necessary to interpret or execute the existing workflow. |
| Source analyses and reference annotations | Perturbation, knockout, knockdown, overexpression and similar terms are retained only where they specify measured conditions or model assumptions. |
| Archives, frozen outputs and other worktrees | These preserve historical or independently owned records and were not rewritten. |

## Verification

- `python3 scripts/manuscript/validate_resource_scope.py`: PASS, 37 required documents.
- `git diff --check`: PASS for tracked changes.
- Python AST inspection: syntax valid; executable structure identical for all
  three Python files changed this turn after excluding docstrings and accounting
  for the scope checker's single wording marker.
- JSON comparison: all configuration values except the mirrored goal prose unchanged.
- R/YAML inspection: changes confined to comments.
- Numerical-token comparison: all existing numbers in changed `docs/` Markdown files unchanged.
- New navigation destinations and the paper-scope heading: checked and present.
- Runtime tests and scientific reruns: not run; no workload or job was submitted.
- Provider classifier behavior: not tested.

## Files changed in this broader pass

- [Analysis/MASLD_Model_Benchmark/README.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/README.md)
- [Cas13_Library_Design/README.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/README.md)
- [Cas13_Library_Design/scripts/build_independent_tier_candidates.py](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/scripts/build_independent_tier_candidates.py)
- [Cas13_Library_Design/scripts/guides/GUIDE_DESIGN.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/scripts/guides/GUIDE_DESIGN.md)
- [Liver_Atlas/README.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Liver_Atlas/README.md)
- [README.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/README.md)
- [RNA-seq/Mouse/Unified_Integration/scripts/M02_mouse_per_diet_de.R](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/scripts/M02_mouse_per_diet_de.R)
- [RNA-seq/README.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/README.md)
- [config/pipeline_params.yaml](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/config/pipeline_params.yaml)
- [config/resource_paper_scope.config.json](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/config/resource_paper_scope.config.json)
- [docs/FIGURE_GUIDELINES.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/FIGURE_GUIDELINES.md)
- [docs/PAPER.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/PAPER.md)
- [docs/README.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/README.md)
- [docs/RESULTS.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/RESULTS.md)
- [docs/ROADMAP.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/ROADMAP.md)
- [docs/STATUS.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/STATUS.md)
- [docs/manuscript/05_figure_legends.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/05_figure_legends.md)
- [docs/manuscript/METHODS.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/METHODS.md)
- [docs/manuscript/README.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/README.md)
- [docs/manuscript/draft/01_intro.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/draft/01_intro.md)
- [docs/manuscript/draft/fig1.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/draft/fig1.md)
- [docs/manuscript/draft/fig2.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/draft/fig2.md)
- [docs/manuscript/draft/title.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/draft/title.md)
- [docs/technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md)
- [docs/technical/MASLD_MODEL_STRATEGY_HANDOFF_2026-09-07.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/MASLD_MODEL_STRATEGY_HANDOFF_2026-09-07.md)
- [docs/technical/OAKES_PACHTER_RIGOR_INTEGRATION_HANDOFF_2026-09-08.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/technical/OAKES_PACHTER_RIGOR_INTEGRATION_HANDOFF_2026-09-08.md)
- [scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md)
- [scripts/manuscript/validate_resource_scope.py](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/validate_resource_scope.py)
- [scripts/portal/gene_catalog_v2.py](/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/gene_catalog_v2.py)

Most documentation files are outside Git tracking in this checkout; the file list
above records those edits as well as the tracked ones. No commit was created.

## Subsequent deeper documentation pass

At the user's request, a further bounded scan read 173 first-party Markdown,
QMD and RST files to six directory levels across the analysis, RNA-seq, genetics,
liver-atlas, portal, documentation and script trees. It found 57 matching lines
in 31 files before editing. These are search matches, not established safety
triggers. Third-party dependencies, generated results, archived/released
material and other worktrees were excluded from edits.

Ten additional documents were edited. Incidental platform labels were replaced
with the study's role or the dataset's in-house label. The patient-cohort diagram
now identifies its separate follow-up branch and directs readers to current
Resource authorities. The integration script index now references `docs/STATUS.md`
and `docs/RESULTS.md` instead of the retired status and number documents.

| Changed document | Change |
|---|---|
| [Deconvolution README](../../Analysis/Deconvolution/README.md) | Descriptive in-house dataset label; dataset key and counts preserved. |
| [Deconvolution plan](../../Analysis/Deconvolution/plan.md) | Project name in the heading. |
| [Continual integration README](../../Analysis/SingleCell/continual_integration/README.md) | Follow-up-study wording; selection exclusions preserved. |
| [SeqFunc README](../../GWAS/finemapping/src/seqfunc_v2/README.md) | Separate-library wording; write restrictions preserved. |
| [Saturation README](../../GWAS/finemapping/src/seqfunc_v2/saturation/README.md) | Separate-library wording; write restrictions preserved. |
| [Patient-cohort README](../../RNA-seq/Human/Patient_Cohorts/README.md) | Historical workflow context and diagram label. |
| [Integration script index](../../RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/README.md) | Current authority references and separate-library description. |
| [Isoform analysis README](../../RNA-seq/scripts/isoform_diversity/README.md) | Transcript-specific RNA-reduction context; literal environment variable retained. |
| [Portal README](../../masld-atlas-v2/README.md) | Separate-study wording; no-deployment and evidence exclusions preserved. |
| [Release guide](RELEASE_PROVENANCE.md) | Plain-language study-separation requirement. |

Retained matches include the literal input-name rejection description in the
continual-integration README, the `CAS13GROUP` environment variable, named source
assays, and historical or generated records. The phrase “cloning this repo” is
ordinary software terminology and was left unchanged. The Western diet-group
analysis combines the DIAMOND, WD and WD_Fructose groups according to its script;
it was not relabeled as an in-house cohort.

Verification: the manuscript scope check passed for 37 required documents;
`git diff --check` passed; all numerical tokens in the ten edited documents were
preserved; new authority destinations, the documented environment variable and
the retained input-filter description were checked against local files.
No executable code, parameter values, data, results or figures changed in this
pass. Runtime tests and scientific reruns were not run. Provider classification
was not tested, and the effect on future flags remains unverified.

## Five independent reviewer pass

Five subagents completed another broad first-party scan, with separate edit
areas for active manuscript documentation, RNA-seq/Analysis, GWAS and analysis
scripts, library/portal text, and an independent repository-wide review. They
edited 15 on-disk files in total. Edits simplified incidental references to the
separate experimental study while keeping the human Resource boundary explicit.

| Review area | Files edited | Result |
|---|---|---|
| Root context | `README.md` | Two incidental references to the future screen became follow-up-study language. |
| Active documentation | `docs/README.md`, `docs/FIGURE_GUIDELINES.md`, `docs/manuscript/METHODS.md`, `docs/manuscript/05_figure_legends.md`, `docs/technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md`, `docs/technical/OAKES_PACHTER_RIGOR_INTEGRATION_ASSESSMENT_2026-09-08.md`, `docs/technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23.md`, `docs/technical/MASLD_MODEL_STRATEGY_HANDOFF_2026-09-07.md`, `docs/technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md` | Incidental scope references simplified; scientific methods, numbers and the contract-bound paragraph in `docs/PAPER.md` retained. |
| Analysis context | `Analysis/Cross_Species_Concordance/README.md`, `Analysis/MASLD_Model_Benchmark/REVIEW_MULTIMODAL_MODEL_STRATEGY_2026-09-07.md` | Separate-study wording clarified in three prose locations. |
| Legacy portal | `masld-atlas-v2/README.md`, `masld-atlas-v2/hf_dataset_stage/README.md`, `masld-atlas-v2/src/app/layout.tsx` | Publication boundary and legacy UI notice simplified; legacy status remains explicit. |
| GWAS and analysis scripts | None | The matches examined were actual assay labels, method descriptions, biological interpretation, or executable input matching. |

Reviewers retained guide-design terminology in the actual separate-study
workflow, source CRISPRi and DepMap assay names, functional-effect annotations, biological
replicate metadata, and identifiers such as `CAS13GROUP` and
`CAS13_SCREEN_FIREWALL`. They also retained the description of names rejected
by the continual-integration input filter. Removing these would obscure an
actual source or make the workflow and study boundaries harder to understand.

The independent scan excluded data, results, archives, releases and other
worktrees. A few nested sample/input directories were inaccessible, so coverage
was broad but not exhaustive. Most documentation and portal files are ignored by
this checkout's `.gitignore`; the edits above exist on disk but do not appear in
ordinary `git diff` or `git status` output.

Final checks: `python3 scripts/manuscript/validate_resource_scope.py` passed
for 37 required documents; `git diff --check` passed for tracked files; all 15
reported edit targets exist; and key scope identifiers and the Resource/follow-up
separation were inspected. Runtime/scientific tests and provider classification
were not run. [Unverified] These changes may or may not affect future flags.

## Second five-reviewer pass

At the user's renewed request, five fresh subagents scanned the active
documentation, RNA-seq and analysis text, portal and library text, GWAS and
figure scripts, and the root instructions. They made targeted wording edits in
five files. The independent reviewer made no edit but identified a useful
method distinction lost during an earlier cleanup; the primary agent restored
that distinction in a sixth file.

| File | Result |
|---|---|
| `docs/PAPER.md` | A broad planning sentence now refers to later functional evidence. The exact contract-bound goal paragraph and scientific claims were retained. |
| `RNA-seq/README.md` | Historical MR/TWAS and GWAS-overlay bullets no longer claim causal gene identification or validation from association alone. |
| `scripts/portal/stage_hf_dataset.py` | Generated dataset-card prose identifies mouse content with the separate experimental study; source data and code behavior are unchanged. |
| `masld-atlas-v2/src/app/downloads/page.tsx` | DepMap CHRONOS is described as supplementary liver-cell-line viability context rather than functional validation. |
| `scripts/figures/load_figure_data.R` | A comment describes the separate follow-up study; the figure directory constant is unchanged. |
| `scripts/portal/gene_catalog_v2.py` | Restored two distinct RNA-depletion options: an RNA-targeting nuclease or an antisense oligonucleotide, plus rescue when feasible. An earlier generic description had lost this useful choice. The screen-evidence rejection logic is unchanged. |

The reviewers retained actual library-design terms, named public assays,
functional-effect annotations, ML `transductive` terminology, source disease-etiology
language, paths, environment variables, and literal validation strings.
Replacing those would lose biological or executable meaning. Existing
uncommitted analytical work was preserved. Some older documents and the portal
tree are ignored by Git, so the edited text was also checked directly on disk.

Final checks after the restored Catalog wording: Resource scope validation
passed for 37 required documents; Python compilation of the two edited portal
modules passed; `git diff --check` passed; direct inspection confirmed the
method distinction, firewall identifier, and reviewed prose edits. Runtime
tests, scientific reruns and provider classification were not run.
[Unverified] The effect on future model flags remains unknown.
