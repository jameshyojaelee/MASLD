# Introduction — working Resource-paper draft

> This draft is governed by `docs/PAPER.md`; it does not make
> Plans 41--46C or the separate Cas13 screen part of the current manuscript.

Metabolic dysfunction-associated steatotic liver disease (MASLD) is a common
chronic liver disorder spanning steatosis, metabolic dysfunction-associated
steatohepatitis (MASH), fibrosis, and cirrhosis. Fibrosis stage is strongly
associated with liver-related and overall outcomes, yet the relationship
between inherited genetic variation and the multicellular remodeling of
established human disease remains incompletely mapped. Liver samples with
similar histological diagnoses can contain different mixtures of hepatocyte
stress, immune activation, stromal remodeling, and vascular dysfunction,
making tissue and cell context essential for biological interpretation.
Histology is therefore an essential clinical anchor but not a complete
description of coordinated molecular state. Following the data-driven
framework of Kamzolas et al. (2026; doi:10.1038/s42255-026-01543-7), a
molecular continuum can add resolution within recorded stages without
replacing pathological assessment.

Human genetics and disease-tissue transcriptomics approach this biology from
different directions. GWAS and fine-mapping identify inherited
trait-associated variation, which can include protein-altering variants,
noncoding variants with regulatory evidence, and signals whose molecular
mechanism remains unresolved. Molecular QTL colocalization can connect some
trait associations to gene regulation, but it does not by itself prove a
causal gene or require persistent differential expression after disease is
established. Conversely, differential expression captures active tissue-state
changes but cannot distinguish inherited drivers from adaptive responses,
cell-composition shifts, or consequences of chronic injury. The two maps
therefore answer different but connected biological questions.

Recent studies have produced large liver-eQTL resources, multicohort
transcriptomic analyses, single-cell and spatial MASLD atlases, and multilayer
gene-prioritization frameworks. These resources have identified candidate
loci, disease-associated cell states, fibrosis programs, and potential
therapeutic mechanisms. What remains missing is a human MASLD Resource that
connects inherited variation to multicohort disease remodeling and then places
both in donor-aware cell and tissue context. This connection must preserve the
distinct biological meaning, coverage, provenance, and inferential unit of each
assay rather than collapsing them into one modality count.

We assembled a multimodal human MASLD Resource to make this connection. The
transcriptomic arm contains 1,257 quality-controlled bulk RNA-sequencing
samples across nine cohorts, with 844 samples from five control-bearing cohorts
entering a pooled, cohort-adjusted disease-versus-control analysis. The genetic
arm integrates prespecified GWAS strata across ancestry groups with liver eQTL
colocalization and fine-mapped coding and noncoding context. Direct disease and
imaging phenotypes remain distinct from liver-enzyme proxy traits. Donor-aware
single-cell analyses and prespecified gene programs resolve multicellular
remodeling, while protein, chromatin-accessibility, and spatial assays provide
assay-native molecular and physical context. We retain fibrosis stage in the
primary severity analyses and add an externally anchored transcriptomic
continuum as a complementary layer for testing coordinated molecular variation
among participants assigned the same stage.

The organizing observation is that trait-associated inherited genetic evidence
and the established disease-state transcriptome provide complementary, largely
non-overlapping maps in the current official release. Exact memberships and
updated overlap values will be stated only after the corrected COLOC portfolio
and synchronized bulk release are adopted. Sparse overlap does not imply that
the two maps are biologically opposed or that their intersection is a
privileged truth set. It shows why inherited variation and active disease
remodeling should be interpreted as distinct dimensions of MASLD biology.

We use two recurring examples to make those dimensions concrete. GNMT begins
as a provisional direct-disease shared-signal locus with an unresolved local
regulatory haplotype, then joins MAT1A and CYP2C19 in coordinated lower
disease-state RNA, cross-sectional fibrosis-stage abundance, and liver protein.
This locus-to-state bridge connects the maps without claiming that the
three-gene metabolic state is one inherited mechanism. Separately,
ECM/IGFBP7 and ductular-injury/BICC1 test whether a disease-associated program,
its cell-of-origin assignment, and its intact-tissue organization support the
same conclusion. They do not: a preliminary ambient-RNA audit raises an
unresolved lineage-assignment concern, while matched spatial nulls support
program-specific organization only for ECM/IGFBP7. The ambient effects require
a synchronized source-diagnosis rerun. The examples were selected
from prespecified analyses, not from a count of supporting modalities.

We therefore mapped these dimensions through donor-resolved cellular programs,
liver protein, accessible chromatin, and intact-tissue organization. Each
assay retains its native effect unit and biological unit. Coverage,
applicability, provenance, and donor identity determine whether an apparent
absence is interpretable, providing guardrails around the biological maps
rather than serving as their headline. The resulting MASLD Gene Catalog links
protein-coding genes, lncRNA transcripts, regulatory DNA elements, and
unresolved loci to inherited variation, multicellular remodeling, cell and
tissue context, evidentiary boundaries, and the next discriminating
experiment. This connected map is the standalone Resource contribution.
