## Liver protein tests the disease-state RNA changes, and intact tissue places the genes in cell types

We next asked how far the genetic and disease-state evidence extends beyond bulk
and single-nucleus RNA, using liver protein and intact tissue in their own units
(Figure 5). Each panel uses its own input set. Model specifications, correction
families and evidence-state definitions are in Methods, and panel construction
and denominators are in the legend. For every panel, the source table records
applicability, evidence state, the unresolved alternative explanation and the
next discriminating experiment.

We first compared RNA and protein changes in the adjusted PXD051911 liver
proteome [ref: PXD051911 liver proteome source study], which includes 58
participants (46 with MASLD and 12 controls) (Figure 5A). The display universe
of 10,044 genes was fixed before the bulk fit of 774 genetic-identity
units, and 3,054 of these genes were measured as proteins. Across these 3,054
genes, RNA and protein changes had Spearman *rho* = 0.458 and 66% direction
agreement (exact values in Supplementary Note 13). Protein-effect agreement between acquisition batches was Spearman *rho* 0.34.
The batches contain different participants, so this combines sampling error,
case mix and technical variation and does not estimate an assay reliability ceiling. Of the 211
protein-significant genes, 188 had the same sign as the RNA change and 23 the
opposite sign. Sign agreement depended on the RNA call. It was 170 of 177 when
the RNA change was itself significant, against 18 of 34 when it
was not. We did not read a nonsignificant protein change as falsifying an RNA
call. Across all 5,836 measured proteins, RNA-protein Spearman *rho* was 0.346.
The display universe is partly selected on RNA. *rho* was 0.487 in the 2,390
genes selected from RNA-based analyses and 0.305 in the 664 genes selected on
genetic evidence alone. The RNA-based selection included differential-expression
contrasts from the first bulk fit of the same five cohorts (846 samples), which
preceded the 844-record fit used for the transfer checks and the 774-unit fit.

The protein data extend the GNMT–MAT1A–CYP2C19 example from Figure 3. All three
proteins were lower in MASLD than in control liver (BH *q* < 0.05), in the same
direction as pooled RNA (Figure 5A). The protein model compared 46 MASLD with 12
control participants and had no severity or hepatocyte-fraction term, so it
measures lower bulk protein in disease. In a post hoc analysis that was not
prespecified, all three proteins were lower at higher NAS across the 46 MASLD
participants (partial Spearman *rho* -0.43 to -0.60; BH *q* ≤ 0.007 over 3
proteins and 2 histology scores). Lower bulk abundance can reflect a smaller
hepatocyte share of the biopsy as well as lower abundance per cell. These data
are not evidence that the three genes share one inherited regulatory mechanism.

We next summarized a fixed 25-protein process view, which carries nine canonical
RNA calls and 23 protein calls (Figure 5B). The same protein cohort contributed
to row selection, so this view is descriptive, not independent validation. CD74
is among the nine RNA calls. HLA-DQA1 has a significant RNA change (BH
*p*_adj 0.025), but its log2FC of 0.38 is below the 0.5 threshold for a canonical
call. Its protein change is not significant (log2FC 1.39, BH *p*_adj 0.080). It
is therefore neither an RNA call nor a protein call in Figure 5B.

We next asked which cell type has the highest measured expression of these genes in intact liver. In the
Xenium in situ data of Tzouanas et al. (2026) [ref: Tzouanas et al. 2026], 14
genes with a protein change in the same direction as RNA (Figure 5A) and 12
colocalized genes are on the 477-gene panel. We used the 9 control livers and
asked which cell type had the highest expression per cell, after removing signal
that leaks from neighboring hepatocytes (Figure 5C; 11 of the 12 colocalized
genes were detected). For 12 of 14 protein-supported genes and 8 of 11
colocalized genes, this cell type matched the gene's top lineage in the
single-nucleus atlas. TAT, GLS2 and CYP3A4 and the colocalized THRB, ADH4 and
CFHR1 had their highest measured expression in hepatocytes, KRT8 and KRT18
in cholangiocytes, MZB1 in B and plasma cells, ACTA2 in stellate and
smooth-muscle cells and KLRC1 in NK and innate T cells. Five genes differed. LYVE1 was highest in sinusoidal
endothelium rather than macrophages, ANXA5 in cholangiocytes rather than
monocytes, and TNFSF10, ANPEP and MEST in hepatocytes rather than macrophages,
cholangiocytes and fibroblasts. Counts per cell favor large cells, so the three
hepatocyte calls may partly reflect cell size. The panel locates genes in
tissue; it does not test disease change. In the same source, the change within
each gene's cell type between MASH-cirrhosis explants and controls passed the
correction family for 0 of 14 protein-supported genes (Holm) and 1 of 12
colocalized genes (SCGN in cholangiocytes, selection-aware exact p = 0.002509,
Holm 0.03011; log2 ratio -2.164, pointwise donor-bootstrap 95% interval
-3.218 to -0.946). The corrected resampling repeats the data-dependent carrier
selection. It does not remove the age and tissue-source differences between
controls and end-stage cirrhosis explants (Methods).

We then scored all 117 fixed Hotspot programs in the same liver proteome
(Figure 5D). This protein cohort did not contribute to program discovery. In 46
MASLD participants, 22 of 95 scorable programs were robustly associated with a
histologic-burden axis (BH over 95 programs; members counted once per
primary-assembly gene). Fourteen were lower with greater burden, including
hepatocyte oxidative-phosphorylation, xenobiotic and amino-acid degradation
programs, and 8 were higher, including hepatocyte stress and NRF2, cholangiocyte
junction and TGFB, and fibroblast hypoxia, ECM proteoglycan and TGFB/collagen
programs. These directions describe lower metabolic-protein abundance and higher
matrix-associated protein abundance with greater histologic burden, in a
cohort independent of program discovery. They do not distinguish tissue
composition from abundance per cell or measure metabolic function. Panel D shows only the robust programs, so it is selected on the
plotted statistic; all 117 are in the source table.

Two further layers did not add a biological claim and are reported in the
supplement. In liver snATAC, colocalized posterior lay in open chromatin 1.8 to
2.4 times more than expected, but mostly in peaks shared by all lineages:
only hepatocyte-specific peaks were enriched (1.77), and a platelet trait showed
the same shared-peak enrichment. MASH-associated accessibility did not pass the
prespecified three-lineage replication family in two donor-level cohorts
(Supplementary Figure S5P–AA). Among the 22 protein-robust programs, the share
with an open-chromatin colocalized gene (8 of 22) was close to that among the
other 73 scorable programs (21 of 73). In Visium tissue (Guilliams et al., 2022;
Vu et al., 2026), the Stromal ECM program exceeded a matched-gene spatial null
after adjustment for all 16 cell2location abundances, but so did 68 and 74 of
114 fixed programs, and 55 of 117 programs were supported in both sources, so
this null does not show program-specific tissue organization (Supplementary
Figure S5, former 5F; Methods).

Taken together, liver protein shows lower bulk GNMT, MAT1A and CYP2C19 abundance
in MASLD and fixed programs whose protein scores track histologic burden, but
does not resolve whether these reflect cell composition or change per cell.
Healthy intact tissue and single-nucleus RNA give the same top lineage for
most Figure 5 genes and different assignments for five. Cell size, normalization
and source composition differ, so these disagreements do not establish which
assignment is correct. Chromatin and Visium autocorrelation
gave context without a program- or lineage-specific result. We next organized
these results into Gene Catalog entries, each with a follow-up experiment that
could distinguish the remaining explanations (Figure 6).
