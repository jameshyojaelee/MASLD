## Liver protein, chromatin and tissue data place genetic and disease-state evidence in each assay's own units

We next asked how far the genetic and disease-state evidence extends beyond bulk
RNA. We placed both in liver protein, accessible chromatin and intact tissue, and
kept each assay in its own units (Figure 5). Each panel uses its own input set.
Model specifications, correction families and evidence-state definitions are in
Methods, and panel construction and denominators are in the legend. For every
panel, the source table records applicability, evidence state, the unresolved
alternative explanation and the next discriminating experiment.

We first compared RNA and protein changes in the adjusted PXD051911 liver
proteome [ref: PXD051911 liver proteome source study], which includes 58
participants (46 with MASLD and 12 controls) (Figure 5B). The display universe
of 10,044 genes was fixed before the bulk fit of 774 genetic-identity
units, and 3,054 of these genes were measured as proteins. Across these 3,054
genes, RNA and protein changes had Spearman *rho* = 0.458 and 66% direction
agreement (exact values in Supplementary Note 13). Of the 211
protein-significant genes, 188 had the same sign as the RNA change and 23 the
opposite sign. Sign agreement depended on the RNA call. It was 170 of 177 when
the RNA change was itself significant, against 18 of 34, a chance level, when it
was not. We did not read a nonsignificant protein change as falsifying an RNA
call. Across all 5,836 measured proteins, RNA-protein Spearman *rho* was 0.346.
The display universe is partly selected on RNA. *rho* was 0.487 in the 2,390
genes selected from RNA-based analyses and 0.305 in the 664 genes selected on
genetic evidence alone. The RNA-based selection included differential-expression
contrasts from the first bulk fit of the same five cohorts (846 samples), which
preceded the 844-record fit used for the transfer checks and the 774-unit fit.

The protein data extend the GNMT–MAT1A–CYP2C19 example from Figure 3. All three
proteins were lower in MASLD than in control liver (BH *q* < 0.05), in the same
direction as pooled RNA (Figure 5B). The protein model compared 46 MASLD with 12
control participants and had no severity or hepatocyte-fraction term, so it
measures lower bulk protein in disease. In a post hoc analysis that was not
prespecified, all three proteins were lower at higher NAS across the 46 MASLD
participants (partial Spearman *rho* -0.43 to -0.60; BH *q* ≤ 0.007 over 3
proteins and 2 histology scores). Lower bulk abundance can reflect a smaller
hepatocyte share of the biopsy as well as lower abundance per cell. These data
are not evidence that the three genes share one inherited regulatory mechanism.

We next summarized a fixed 25-protein process view, which carries nine canonical
RNA calls and 23 protein calls (Figure 5C). The same protein cohort contributed
to row selection, so this view is descriptive, not independent validation. CD74
is among the nine RNA calls. HLA-DQA1 has a significant RNA change (BH
*p*_adj 0.025), but its log2FC of 0.38 is below the 0.5 threshold for a canonical
call. Its protein change is not significant (log2FC 1.39, BH *p*_adj 0.080). It
is therefore neither an RNA call nor a protein call in Figure 5C.

To ask whether colocalized variants lie in accessible chromatin of liver cell
lineages, we intersected the colocalization posterior of 436 genes with open
chromatin from two liver snATAC sources [ref: two liver snATAC source studies]
across five shared lineages (Figure 5D). These genes form 767 gene-study pairs:
22 genes have a direct MASLD or liver-fat pair and 424 a liver-enzyme pair, and
ten have both. For 12 of 22 direct-trait genes and 167 of 424 liver-enzyme genes,
less than 5% of the posterior lay in peaks in every lineage. We compared the
posterior mass in peaks with the expectation from the tested variants of the same
window. Averaged over liver-enzyme genes, mass in peaks was 2.5 to 3.1 times this
expectation in each lineage. In no lineage, however, did more than half of these
genes exceed their own expectation (94 to 212 of 424). For the 22 direct-trait
genes, the 95% interval of this ratio included 1. Peak overlap does not assign a
variant to a neighboring gene. Absent overlap is unresolved context, not a
negative test.

We then placed all 117 fixed Hotspot programs in chromatin and protein context
without selecting them by disease association (Figure 5E). In total, 31 programs
contained at least one gene from a colocalized pair that placed at least half of
its posterior in open chromatin of the program's own lineage in at least one
snATAC source (genes counted once per primary-assembly gene). In 3 of these
programs (hepatocyte Ion/membrane potential, cholangiocyte FXR/bile-acid and
fibroblast Elastic-fiber ECM), that gene was RORA, which is open through its
direct MASLD pair. In the other 28, genes were open only through liver-enzyme
pairs. RORA is also open in the T/NK Memory/MAIT program, but only through
liver-enzyme pairs. Its direct pair places 0.37 of its posterior in T/NK peaks of
one source, below the 0.5 threshold. In liver protein from 46 MASLD
participants, 22 of 95 scorable programs were robustly associated with a
histologic-burden axis across participants (BH over 95 programs; members counted
once per primary-assembly gene). Of these, 14 were negative, including
oxidative-phosphorylation, xenobiotic and amino-acid degradation programs, and 8
were positive, including fibroblast ECM proteoglycan and TGFB/collagen programs.
This protein cohort did not contribute to program discovery. Neither spatial
program in Figure 5F contains a colocalized gene or can be scored in protein.

Finally, we examined tissue organization for two programs, ECM/IGFBP7 and
ductular-injury/BICC1, labeled Stromal ECM and Ductular injury in the panel
(Figure 5F). Both were chosen in advance from an earlier version of the single-cell
stage model (Methods). Neither passes the 117-program BH and
stability rule for a pooled stage association in the single-cell stage analysis (Figure 4C). In two
spatial transcriptomics sources (Guilliams et al., 2022; Vu et al., 2026), the
prespecified design adjusted spot scores for cell2location hepatocyte abundance,
the discovery lineage of both programs [ref: cell2location].
We then compared residual spatial autocorrelation (Moran's I [ref: Moran's I])
with a matched null built from genes of equivalent expression, detection,
biotype and correlation with hepatocyte abundance. ECM/IGFBP7 exceeded this null
in both sources (excess 0.139 and 0.147). Ductular-injury/BICC1 did not (0.026
and 0.019) and is indeterminate, not negative. Spot scores use only the measured member genes: 26 and 25 of the 49 ECM/IGFBP7 genes (55% and 53% of its weight) and 9 of the 22 ductular-injury/BICC1 genes (41%). A model added afterwards as a
supplementary test of composition dependence, which adjusts for all 16
cell2location abundances and is plotted in Figure 5F, gave the same evidence
states with smaller excesses (ECM/IGFBP7 0.040 and 0.050; ductular-injury/BICC1
0.015 and 0.003). It removed 71% and 66% of the ECM/IGFBP7 excess, so most of
the prespecified excess is accounted for by the estimated abundances of the
other cell types. The two sources were deconvolved against different atlas
reference builds (Methods). Under that
model, residual Moran's I was positive
for ductular-injury/BICC1 in 4 of 4 donors and 9 of 10 arrays. The raw statistic
therefore does not separate the two programs, and no prespecified test compared
them. The matched-gene null holds donors fixed, so it is not a donor-level test. It also
draws each control gene independently, so control sets lack the co-expression
of a fixed program, and it matches on correlation with hepatocyte abundance, not
with fibroblast or cholangiocyte abundance, the lineages in which these two
programs score highest (Figure 4). The all-16 model removes only the linear part
of that mismatch. We therefore applied the matched-gene null, under the
all-16 model and with each program's controls matched on its own discovery
lineage, to all 114 fixed programs with enough measured genes, as a
descriptive analysis that was not prespecified. In Guilliams et al., 68 of 114
programs exceeded it, and in Vu et al., 74 of 114 did (BH q < 0.05 over 114
programs within each source). ECM/IGFBP7 ranked 42nd and 29th of 114 by
excess, and programs of hepatocyte plasma-protein, acute-phase and metabolic
genes had the largest excess. With control genes matched on correlation with
fibroblast abundance instead, ECM/IGFBP7 still exceeded the null (excess
0.035, q = 0.011; 0.047, q = 0.0008), although its member genes remained more
correlated with fibroblast abundance than the matched controls. With controls
matched on correlation with cholangiocyte abundance, ductular-injury/BICC1
still did not (0.015, q = 0.10; -0.001, q = 0.50). Exceeding this null
therefore does not show program-specific tissue organization.
At the donor level, excess over the matched null was positive in 3 of 4 donors
for ECM/IGFBP7 (exact one-sided sign-flip *p* = 0.19) and for
ductular-injury/BICC1 (*p* = 0.125); ductular-injury/BICC1 excess was also
positive in 7 of 10 Vu et al. arrays. Figure 5F shows these excess counts, not
the raw residual counts above. With four donors, no sign-flip test can reach *p* < 0.05. The four Guilliams et al. donors were sampled during surgery for other conditions: two gastric bypasses for obesity, a cholecystectomy, and a resection of colorectal-cancer liver metastases in a donor who reported six standard alcohol units daily (Guilliams et al., 2022, Table S6). Recorded steatosis ranged from none to 70%, and fibrosis from none to clear pericellular fibrosis. These donors are not a MASLD cohort, and two of them (H37 and H38) also contribute single-nucleus profiles to the atlas, so the spatial and single-cell sources share participants. Methods describe the donor joins, array-level source dependence and
matched-null construction.

Taken together, liver protein shows lower bulk GNMT, MAT1A and CYP2C19 abundance
in MASLD but does not resolve whether this reflects hepatocyte content or change
per cell. Lineage chromatin places colocalized posterior in accessible regions
for a subset of genes without assigning target genes. Of the two spatial
programs, only ECM/IGFBP7 exceeds the matched-gene null within sections, but
most fixed programs (68 and 74 of 114) also exceed it, so the spatial result
does not show program-specific organization. Disease association,
cell-of-origin assignment and tissue architecture are separate claims, and each
needs its own calibration. We next organized these results into Gene Catalog
entries, each with a follow-up experiment that could distinguish the remaining
explanations (Figure 6).
