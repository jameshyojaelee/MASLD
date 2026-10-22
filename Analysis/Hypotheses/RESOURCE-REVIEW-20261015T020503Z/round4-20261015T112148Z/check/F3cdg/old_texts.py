E=[
("R03-01","docs/manuscript/05_figure_legends.md","119-122","""family. The union contains 13,761 genes: 9,065 called in one cohort, 3,462 in
two, 922 in three, 261 in four and 51 in all five; ribbon widths represent
20,054 cohort memberships."""),
("R03-01","docs/manuscript/draft/fig3.md","24-27","""the same canonical threshold. Their union contains 13,761 genes: 9,065 called
in one cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.
These are 20,054 cohort memberships, distinct from the 1,854 pooled calls."""),
("R03-01","docs/RESULTS.md","504-507","""13,761 distinct genes and 20,054 cohort memberships: 9,065 genes called in
one cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five."""),
("R03-01","docs/PAPER.md","392","The per-cohort union contains 13,761 genes, including 9,065 called in one cohort and 51 in all five; it is distinct from the 1,854 pooled calls."),
("R03-01","docs/PAPER.md","392","A broad established-disease signal recurs across cohorts and holds under a stricter histologic endpoint."),
("R03-01","docs/STATUS.md","301","The five-cohort union contains 13,761 genes: 9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts."),
("R03-01","docs/STATUS.md","814","Its 6,147 calls join a 13,761-gene per-cohort union (9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts), distinct from pooled membership."),
("R03-02","docs/manuscript/05_figure_legends.md","123-127","""all male, with constant sex omitted. The restored model yields 6,147 calls,
independently reproduced against native fragments. Its identity screen removed
24 of 31 source disease records, including all ten women, and 2 of 25 controls;
[Inference] unequal eligibility limits representativeness."""),
("R03-02","docs/manuscript/draft/fig3.md","32-33","and 2 of 25 controls. [Inference] Its unequal eligibility limits representativeness."),
("R03-02","docs/RESULTS.md","521-523","""The other four cohort estimates are unchanged. [Inference] Unequal identity
eligibility may affect representativeness; genetic clusters are not verified donors."""),
("R03-02","docs/PAPER.md","413-414","""counts and an independent group-means refit reproduce the result. [Inference]
Unequal eligibility limits representativeness; these clusters are not verified donors."""),
("R03-02","docs/manuscript/METHODS.md","144-145","""source women, and 2 of 25 controls. [Inference] Unequal eligibility limits
representativeness and is not resolved by an estimable design."""),
("R03-09","docs/manuscript/draft/fig3.md","38-39","""correlations of 0.870–0.995 and threshold retention 70.7%–96.6%. Neither
procedure is independent replication of a pooled result."""),
("R03-09","docs/RESULTS.md","120","| Strict-depth sensitivity, corrected counts (four-cohort disease contrast: GSE126848 keeps 5 controls and no disease unit) | 738 / 738 | 24,509 | 1,906 (1,344 / 562) | 2,135 |"),
("R14-07","docs/manuscript/05_figure_legends.md","158-159","3,776 finite native RNA–genetics pairs (Spearman rho=0.022, nominal P=0.18),"),
("R14-07","docs/manuscript/draft/fig3.md","139-141","""log2FC. Its 3,776 finite native RNA–SuSiE pairs have Spearman rho=0.02175
(nominal P=0.181), without a fitted line."""),
("R14-07","docs/manuscript/draft/fig2.md","59-61","""interface, PP.H4 and raw disease-state log2 fold change have Spearman
rho=0.02175 (nominal P=0.18146)."""),
("R14-07","docs/RESULTS.md","1045-1046","testable genes with 426 genetic. Spearman rho is 0.022, nominal P=0.18"),
("R14-07","docs/STATUS.md","830","Both 26; rho 0.0217507, nominal p 0.181458."),
("R16-02","docs/manuscript/05_figure_legends.md","147-149","""**(F)** Fifty-four genes with BH FDR <0.05 in at least one adjacent fibrosis
contrast and membership in the native 436-gene primary SuSiE set of the corrected
COLOC r2 input."""),
("R16-02","docs/manuscript/draft/fig3.md","124-126","""keeping the measurements separate. Fifty-four genes belong to the native 436-gene
primary SuSiE set and pass BH FDR <0.05 in at least one adjacent fibrosis
contrast."""),
("R16-02","docs/RESULTS.md","32-33","PP.H4 > 0.5 and BH FDR < 0.05 in at least one adjacent fibrosis contrast. Of the"),
("R14-09","docs/manuscript/draft/abstract.md","38-42","""A worked GNMT example connects a provisional direct
disease shared-signal locus to a severity-linked GNMT–MAT1A–CYP2C19 state,
defined by coordinated decreases across disease-state RNA, fibrosis-stage
contrasts, and liver protein. This bridge does not imply that all three genes
share one inherited mechanism."""),
("R14-09","docs/manuscript/draft/fig3.md","133-137","""0.128); GNMT has one direct-NAFLD single-signal ABF result above 0.5 (PP.H4
0.706, 2021_34841290_NAFLD_EUR) and is therefore the provisional bridge between
the genetic and state maps."""),
("R14-09","docs/manuscript/draft/fig2.md","65-68","""GRCh37 6:42931261 (rs2296804, SNP.PP.H4 0.258990); its primary multi-signal result remains below threshold
(SuSiE PP.H4 0.128)."""),
("R14-09","docs/PAPER.md","292-294","""and ABF PP.H4 0.706 (2021_34841290_NAFLD_EUR), so GNMT remains an ABF-supported
provisional bridge;"""),
]
import os
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
for fid,f,ln,t in E:
    s=open(R+f).read(); n=s.count(t)
    if n:
        i=s.index(t); a=s[:i].count('\n')+1; b=a+t.count('\n')
        print(fid,f,ln,'FOUND x%d at %d-%d'%(n,a,b))
    else:
        first=t.split('\n')[0]; 
        hits=[k+1 for k,l in enumerate(s.split('\n')) if first[:40] in l]
        print(fid,f,ln,'NOT FOUND; first-line hits',hits)
