# MASLD Gene Catalog Resource

Human MASLD Resource paper: liver-trait genetic evidence, cross-sectional
disease-state liver expression, and cell and tissue context, organized in the
MASLD Gene Catalog. Colocalization nominates candidates; it does not establish
causal genes. The Cas13 experimental follow-up is a separate paper.

## Start here

| Task | Open |
|---|---|
| Read the working paper | [Main PDF](docs/manuscript/latex/build/main.pdf) · [Supplement PDF](docs/manuscript/latex/build/supplement.pdf) |
| Edit the manuscript | [Manuscript workspace](docs/manuscript/README.md) · [Main LaTeX](docs/manuscript/latex/main.tex) · [Supplement LaTeX](docs/manuscript/latex/supplement.tex) |
| Check scientific scope and current limitations | [Paper identity](docs/PAPER.md) · [Current status](docs/STATUS.md) |
| Cite an exact number | [Results and release qualifications](docs/RESULTS.md) |
| Find a main panel or source table | [Figures 1–7 index](figures/main/INDEX.md) |
| Browse analyses | [Analysis directory](Analysis/README.md) · [Hypothesis folders](Analysis/Hypotheses/README.md) · [Dataset and pipeline routes](docs/technical/DATASETS_AND_PIPELINES.md) |
| Browse figure alternatives | [Figure folders](figures/README.md) · [Candidate directory index](figures/candidates/README.md) |
| Open the selected local Gene Catalog | [Catalog](RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/index.html) · [Eligibility decision and qualifications](docs/RESULTS.md#catalog-eligibility-adoption-2026-11-06) |
| Find the next manuscript task | [Roadmap](docs/ROADMAP.md) |

## Read first

The [documentation guide](docs/README.md) explains the five scientific
authorities. Their current qualifications apply to the working PDFs and Catalog:
the selected genetic map remains provisional after variant-identity and LD
checks, and full unique-cell verification remains incomplete. Formatting and
directory indexes do not establish submission readiness.

## Main figures

| Figure | Panel directory | Subject |
|---|---|---|
| 1 | [Atlas overview](figures/main/fig1_atlas_overview/README.md) | Resource design, cohorts, assays and evidence domains |
| 2 | [Genetics](figures/main/fig2_genetics/README.md) | Fine-mapping, colocalization and ancestry |
| 3 | [Bulk transcriptomics](figures/main/fig3_bulk_transcriptomics/README.md) | Multicohort disease and cross-sectional stage remodeling |
| 4 | [Single-cell programs](figures/main/fig4_singlecell_programs/README.md) | Fixed programs, lineage context and molecular continuum |
| 5 | [Molecular context](figures/main/fig5_molecular_context/README.md) | Protein and cell location; chromatin and spatial supplements |
| 6 | [Gene Catalog](figures/main/fig6_gene_catalog/README.md) | Evidence states and follow-up experiments |
| 7 | [Model framework](figures/main/fig7_model_framework/README.md) | Predictive models and training-data exposure |

## Models

<a id="models-we-built"></a>
<a id="models-that-failed"></a>
<a id="external-models-benchmarked"></a>

Use the [model registry](Analysis/MASLD_Model_Benchmark/MODELS.md) for model
states, scores and transfer limits, and the [analysis overview](Analysis/MASLD_Model_Benchmark/README.md)
for its files. Earlier model inventories in this README are retained in the
[documentation snapshot](docs/archive/folder_navigation_2026-11-09/ROOT_README.md#models).

## Repository map

| Location | Purpose |
|---|---|
| [docs/](docs/README.md) | Scientific authorities, manuscript, technical references and documentation history |
| [RNA-seq/](RNA-seq/) | Bulk cohort processing and downstream integration |
| [GWAS/](GWAS/) | Genetic data, fine-mapping and colocalization |
| [Analysis/](Analysis/README.md) | Single-cell, chromatin, spatial, protein and other analyses |
| [data/](data/REGISTRY.md), [Liver_Atlas/](Liver_Atlas/README.md) | Dataset access and existing physical storage |
| [figures/](figures/README.md) | Manuscript panels, supplementary artwork and alternatives |
| [scripts/](scripts/) | Figure, manuscript, Catalog and shared analysis scripts |
| [config/](config/), [envs/](envs/README.md), [tools/](tools/) | Parameters, environment specifications and external tools |
| [masld-atlas-v2/](masld-atlas-v2/README.md) | Portal source; scientific values follow the selected Catalog |
| [logs/](logs/), [results/](results/), [outputs/](outputs/), [archive/](archive/) | Existing logs, other outputs and historical material |

Dataset and named-release paths are kept in place. The separate experimental
follow-up is in [Cas13_Library_Design/](Cas13_Library_Design/README.md).
