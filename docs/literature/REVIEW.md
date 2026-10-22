# Competitor landscape and Resource positioning

**Purpose:** compact framing guide. Detailed dataset sizes, panels, and source
quotes remain in the paper-specific reviews and imported PDFs.

## Closest recent papers

| Paper/resource | Main overlap with this paper | What it does especially well | What remains distinctive here |
|---|---|---|---|
| *A data-driven framework reconstructs the molecular continuum of human MASLD progression* | Cross-sectional stage ordering, sliding-window/continuum analyses, and biomarker evaluation | Focused continuum narrative and clinically legible stage structure | Broader regulatory-genetic, single-cell, chromatin, protein, spatial, observability, and experiment-routing framework; our stage language remains cross-sectional |
| *LiRNA: An interactive atlas of human liver RNA-seq databases* | Large-scale liver RNA-seq harmonization and public exploration | Accessible liver transcriptomic reference and portal-like reuse | Disease-specific analysis rules, donor-aware multimodal context, regulatory genetics, frozen programs, source dependence, and the MASLD Gene Catalog |
| *Network-based meta-analysis maps stage-dependent molecular programs in MASLD* (Kumak, Darde & Konu 2026; doi:10.64898/2026.08.26.747338) | Nine-cohort bulk DE meta-analysis on NAS and fibrosis, stage-resolved pathway modules, and a public network portal | Legible stage-resolved module emergence, an honest reassessment of six published panels, and a usable Shiny portal with ortholog mapping | No genetics of any kind, no single-cell, spatial, protein or chromatin context, donor unit not enforced (GSE193066 repeat biopsies counted as independent), and modules built from outcome-selected genes; see [the detailed analysis](../competitor_analysis/kumak_konu_2026_masld_meta_network.md) |
| *Integrative analyses elucidate transcriptional regulatory functions of risk alleles for metabolic liver disease* | GWAS/eQTL/regulatory-variant interpretation and experimental regulatory follow-up | Variant-to-regulatory-function depth and direct perturbational evidence | Connection of inherited genetic variation to multicohort disease remodeling, cell/tissue context, and explicit experimental routing |

The first paper overlaps strongly with individual Figure 3 techniques; LiRNA
overlaps with the transcriptomic resource role; and the regulatory-risk paper
overlaps with Figure 2. None of those facts invalidates this manuscript. They
make it important to avoid claiming novelty for a sliding window, an atlas, or
COLOC alone. The defensible novelty is the connection of **inherited genetic
variation, multicellular disease remodeling, and cell-resolved tissue
context** in a reusable Gene Catalog. The evidence rules and experiment
routing make that biological map trustworthy and actionable.

## Broader competitor classes

| Competitor class | Crowded claim | Better positioning for this paper |
|---|---|---|
| Large MASLD bulk meta-analyses | “We found robust DEGs” | Show how disease-state signals are made cohort-aware, cross-stage, multimodal, and experimentally routable |
| Single-cell/spatial liver atlases | “We mapped disease cell states” | Use donor-aware lineage and physical context to qualify gene evidence, not claim ownership of every cell state |
| GWAS fine-mapping and regulatory studies | “We found causal risk genes” | Provide phenotype/ancestry/source-aware shared-signal evidence and clearly separate nomination from causality |
| Multi-omics target rankings | “Our convergent score finds the best targets” | Reject universal vote counting; expose testability, discordance, source dependence, and next experiment |
| MASLD biomarker studies | “Our signature predicts disease” | Keep biomarker analysis as a worked application without prospective-performance claims |
| Perturbation screens | “We functionally validate targets” | Freeze a Resource independently, then let the separate Cas13 paper test it prospectively |

## The novelty sentence

This work provides a multimodal human MASLD Resource that maps inherited
genetic variation and multicellular disease remodeling, resolves their cell
and tissue contexts, and converts the connected biology into transparent Gene
Catalog entries specifying the next discriminating experiment.

## Reviewer question and answer

**“Is this just another multi-omics atlas?”**

No. The Resource connects trait-associated inherited variation to multicohort
disease remodeling and resolves both across donor-aware cell programs and
tissue contexts. It also records whether each signal was testable, how much it
depends on a source or proxy phenotype, and which assay should come next. The
paper's reusable product is the MASLD Gene Catalog and its reproducible
interpretation rules.

## What would weaken the positioning

- leading with total modalities, dataset count, or “largest atlas”;
- presenting a source paper’s signature, variant, or cell state as ours;
- treating sliding-window, NMF, COLOC, or a portal as individually novel;
- hiding negative or untestable public-data branches;
- collapsing heterogeneous evidence into a single target score; or
- implying that the separate Cas13 screen was used to optimize this Resource.

## Journal implications

The strongest honest identity is a high-value Resource/framework paper. Cell
Genomics is the clearest current fit; Genome Biology and Nature Communications
are credible alternatives. A metabolism-journal inquiry should emphasize the
biological utility and prospective experiment-routing value, not retrofit a
causal mechanism that the public data do not establish.

Recheck literature, journal formats, and author instructions immediately before
submission; this page is a positioning summary, not a permanently current
bibliography.
