# Figure 1 Results — a multimodal human MASLD Resource connects inherited variation and disease remodeling

> Working prose under `docs/PAPER.md`. Values are the adopted 2026-08-17
> release. Evidence-state definitions, coverage rules and biological-unit rules
> are in Methods; complete observability inventories are in the Supplement.

## The Resource connects two biological maps to cell and tissue context

We assembled a multimodal human MASLD Resource that connects trait-associated
inherited genetic variation, multicohort disease remodeling, and donor-aware
cell and tissue context (Figure 1A). Every scientific claim is
linked to its source cohort, assay, phenotype, biological unit, testing
universe, reuse relationship, evidence state, unresolved alternative, and next
discriminating experiment (Figure 1A). The transcriptomic foundation
contains 1,257 technically passing samples from nine retained cohorts; 844
samples from five control-bearing cohorts contribute to the pooled
cohort-adjusted disease-versus-control model. The regulatory-genetic arm spans
prespecified liver GWAS strata across ancestry groups, with
direct MASLD or imaging traits kept distinct from ALT, AST, and GGT proxy
phenotypes. Donor-collapsed single-cell data, prespecified gene programs,
liver protein, chromatin-accessibility, and spatial measurements supply
assay-native context rather than additional votes in one target score.

The organizing interface compares jointly testable inherited genetic and
established disease-state evidence (Figure 1B). Among 14,920 genes assayed in
both, 437 are colocalized and 1,066 are disease-state genes, and 20 are both
against 31.2 expected under independent sampling. The two maps therefore
intersect at close to the rate chance produces for sets of these sizes, and
Figure 2 shows why. Throughout, replication is defined by participants or
biological donors; cells, spots, sequencing runs, guides and arrays never
inflate the biological sample size.

The output is a source-dependency graph and an analysis-to-decision workflow
(Figure 1D). Figure 2 maps phenotype- and ancestry-aware inherited variation
and its candidate molecular links. Figure 3 defines multicohort,
cross-sectional disease remodeling and its interface with genetics. Figure 4
resolves donor-level multicellular programs and regulatory activity. Figure 5
places candidate genes and programs in molecular and physical tissue context.
Figure 6 organizes the connected maps in the MASLD Gene Catalog and identifies
the smallest experiment capable of resolving each remaining alternative. No
step requires agreement across all modalities, and no aggregate rank is used as
a calibrated measure of biological truth.

Two examples recur through that sequence: GNMT–MAT1A–CYP2C19 connects the
RNA, cross-sectional stage, and protein views of hepatocyte metabolic
remodeling, while ECM/IGFBP7 and ductular-injury/BICC1 connect prespecified
donor-level programs to matched-null spatial interpretation. Their current
official COLOC states remain visible but do not determine their
disease-remodeling roles.

This architecture makes the Resource independently useful: readers can trace
why a candidate was nominated, determine what has genuinely been tested, and
select a discriminating experiment without access to any future Cas13-screen
outcome. The separate prospective perturbation study may consume a read-only
Resource release, but its roster, guides, phenotypes, hits, and mechanisms do
not enter this paper.
