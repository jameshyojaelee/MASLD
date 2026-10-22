## Donor-level single-cell programs resolve cell-state and regulatory context

> Working Figure 4 prose. The additive Figure 4F–4G continuum panels are
> adopted as working figures; synchronized Resource release remains open. Communication
> and TF activity remain supplementary.

Seven datasets comprising 1,232,318 cells from 104 descriptive biological
donors were integrated into one annotated liver atlas (Figure 4A); 102 donors
remain in the analyzed atlas after two prespecified low-cell library
exclusions. The analysis retains 117
prespecified cell-type-specific gene programs across five tested lineages and
uses donor-level summaries for program, communication, and regulator analyses.
Sequencing runs, sort fractions, and cells are not biological replicates.

The same shared Resource evidence rules govern the atlas. Donor identity, retained gene
weight, cell-state coverage, and complete multiple-testing families are all
required before a program receives an evidence state. Failure of any one is
`untestable`, not evidence that the biological program is absent. Each selected
program retains its unresolved alternative and next discriminating experiment.

The integrated UMAP directly labels cell types and omits the redundant cell-type
legend (Figure 4B). Disease-density surfaces are estimated per biological donor
after sequencing runs are mapped to donor identity, then averaged within each
disease group. The density scale is
the sole legend. The UMAP is descriptive and is not used as the inferential
unit.

A program's name records the cell-type subset it was discovered in, not what
its genes are. Scored against donor-level per-lineage expression, **80 of 117
programs match their declared lineage and 37 do not**, and the mismatch is
systematic: 17 are programs discovered outside hepatocytes whose content is
hepatocyte, several named for it outright, including a "Hepatocyte xenobiotic
(CYP2C)" program labelled fibroblast and a "Bile-acid transport (ABCB11)"
program labelled cholangiocyte. Ambient RNA accounts for it: those 17 carry a
median ambient burden of 0.432 against 0.061 for correctly labelled
non-hepatocyte programs, and hepatocyte transcripts are the dominant ambient
species in liver droplet data, so contamination is reading as cell-type
identity. Panels are therefore labelled by content, with the discovery subset
kept in the legend; the program definitions are unchanged.

Figure 4C shows the three positive and three negative program effects with the
smallest complete-family BH q per lineage, with all 117 programs retained in
the supplementary view. ECM/IGFBP7 and ductular-injury/BICC1, both discovered in
hepatocyte profiles, pass the primary BH/stability rule but not HC3 family
correction; both are among the 37 whose content contradicts that label. Their
ordinal-stage models use 64 biological donors across five datasets with dataset
as a covariate. The aligned 117-program ambient audit attenuates ECM/IGFBP7 from
beta 0.401 to 0.130 and BICC1 from 0.396 to 0.389, and neither
correction-induced change passes its BH family, so both remain
`attenuation_indeterminate`. Taken with the content audit, this prevents any
hepatocyte-intrinsic claim from the original atlas.

Figure 4D asks the separate question of cellular context. In donor-paired
same-atlas contrasts, Stromal ECM/IGFBP7 is highest in fibroblasts and Ductular
injury/BICC1 is highest in cholangiocytes. The ten prespecified contrasts pass
both ordinary and HC3 family criteria for those focal lineages. This is
descriptive localization, not transcript origin or cell autonomy.

Figure 4E transports the fixed members of those two programs into independent
bulk F1–F4-versus-F0 tissue contrasts. Both weighted member-gene effects rise
with fibrosis stage. This is cross-sectional tissue-state transport, not a
program-level significance test. Communication and TF activity remain in Figure
S4 because they do not strengthen this lineage-localization result.

Figure 4F then asks whether the same focal programs vary along a
continuous molecular axis after recorded fibrosis stage is retained in the
model. In GSE162694 and GSE213621, both Stromal ECM and Ductular injury increase
along the fixed external projection in stage- and sex-adjusted donor-level
models, with concordant effects for the independently fitted cohort-PC1 score.
The left side shows descriptive fixed-window trajectories in all five Resource
cohorts. GSE126848, GSE130970, and GSE135251 overlap continuum construction and
are visually separated from the two competitor-independent cohorts. The right
side shows the two-score stage/sex-adjusted meta-analysis restricted to
GSE162694 and GSE213621. GSE126848 is omitted from stage-stratified panels
because recorded fibrosis is unavailable.
The complete 117-program family, transcript comparisons, pathway collections,
and ten fixed NMF axes remain in Figures S3–S4. Their alignment with
fibrosis-based results supports a complementary interpretation: histology
anchors severity, whereas the continuum resolves coordinated molecular
heterogeneity among participants sharing a stage. It does not establish a
temporal trajectory or supersede the fixed program definitions.

Figure 4G gives the continuum a unit of interpretation between single genes and
whole-cohort axes. It compresses 13,092 correlated features into a few dozen
auditable groupings so that remodeling within a recorded stage can be described
as coordinated rather than gene by gene. The panel therefore broadens the focal
view to the complete testable pathway and Hotspot
layers. Before any stage or continuum outcome was read, signature-excluded gene
membership placed 5,439 nodes in one graph component and 43 size-ranked systems;
six nodes in three disconnected two-node microcomponents remain explicit but
outside system summaries. Stage-adjusted continuum coefficients are then
overlaid on the fixed coordinates. Across nodes, fibrosis-stage and continuum
effects have Spearman rho=0.683 and 80.5% share direction. At the system level,
equal-collection summaries give rho=0.658 and 36 of 43 shared directions;
feature-balanced medians give rho=0.701 and 37 of 43. Fixed-window direction
agrees between both evaluation cohorts for 39 of 43 systems under either
aggregation, but only 13 of 43 or 12 of 43 retain one direction across all five
cohorts. Thus, broad biological systems tend to remodel in aligned directions
along stage and continuum, while cohort-specific window shapes remain common.
These are descriptive summaries over nested gene sets and overlapping windows,
not system-level significance tests, temporal ordering, or a universal ontology.

Pathway similarity networks partitioned into modules are likewise an established
way to summarize MASLD transcriptomes (Kumak et al., 2026;
doi:10.64898/2026.08.26.747338), and Figure 4G claims no novelty for that
construction. The distinction here is the order of operations. The 43 systems
are built from signature-excluded gene membership before any stage or continuum
outcome is read, and the 117 Hotspot programs keep memberships and weights fixed
in advance, so the geometry cannot have been selected by the outcome it is used
to describe. Centrality in a curated protein interaction database is not used to
rank genes or programs at any point, because degree in such a database reflects
annotation depth as well as biology.

The complete supplement therefore shows both the five cohort-specific window
maps and an equal-cohort pooled view of all 43 systems. The pooled heatmap
balances collections within each cohort and then gives each cohort one-fifth
weight. It summarizes the shared cross-cohort shape without allowing the
largest cohort to dominate, while retaining the discordant cohort-specific
trajectories in the adjacent panel. Both panels use the same descriptive pooled
late-minus-early row order, and the pooled panel displays that contrast as a
final column. This makes increasing and decreasing shared profiles legible but
does not turn the order into a support or timing statistic.

Together, these panels provide donor-aware cell-labelled and continuous
molecular context for cross-sectional stage-associated remodeling without
redefining the program registry or treating either histologic groups or the
continuum as longitudinal states. Current
official COLOC provides no positive genetic call for either IGFBP7 or BICC1,
and the unresolved ambient sensitivity prevents assigning either program to hepatocytes from
these data alone. Figure 5 then asks the independent question of whether their
tissue organization exceeds assay-matched spatial background.
