import json
E = {}
E['R03-01'] = [
 dict(file='docs/manuscript/05_figure_legends.md', line='119-122', claim_level=False,
  old_text="""family. The union contains 13,761 genes: 9,065 called in one cohort, 3,462 in
two, 922 in three, 261 in four and 51 in all five; ribbon widths represent
20,054 cohort memberships.""",
  new_text="""family. The union contains 13,761 genes: 9,065 called in one cohort and 4,696
in two or more. Of the 4,696, 2,831 have the same sign in every calling cohort
(1,933 in two, 636 in three, 218 in four and 44 in all five) and 1,865 are
called up in one cohort and down in another (Mixed sign); 1,632 of the 1,865
include a GSE126848 call. Ribbon widths represent 20,054 cohort memberships.
Shared membership is not replication; direction agreement of the pooled calls
is shown in Figure S3C."""),
 dict(file='docs/manuscript/draft/fig3.md', line='24-27', claim_level=False,
  old_text="""the same canonical threshold. Their union contains 13,761 genes: 9,065 called
in one cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.
These are 20,054 cohort memberships, distinct from the 1,854 pooled calls.""",
  new_text="""the same canonical threshold. Their union contains 13,761 genes: 9,065 called
in one cohort and 4,696 in two or more. Only 2,831 of the 4,696 have the same
sign in every cohort that calls them (1,933 in two, 636 in three, 218 in four
and 44 in all five); 1,865 (39.7%) are called up in one cohort and down in
another, and 1,632 of these include a GSE126848 call. Without GSE126848, 346 of
2,565 multi-cohort genes (13.5%) have mixed signs. These are 20,054 cohort
memberships, distinct from the 1,854 pooled calls."""),
 dict(file='docs/RESULTS.md', line='504-507', claim_level=False,
  old_text="""13,761 distinct genes and 20,054 cohort memberships: 9,065 genes called in
one cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.""",
  new_text="""13,761 distinct genes and 20,054 cohort memberships: 9,065 genes called in
one cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.
Counted by direction (review R03-01, 2026-10-15), 1,865 of the 4,696
multi-cohort genes (39.7%) are called up in one cohort and down in another
(1,529 / 286 / 43 / 7 in the two- to five-cohort tiers); 2,831 have one sign in
every calling cohort (1,933 / 636 / 218 / 44). 1,632 of the 1,865 include a
GSE126848 call; without GSE126848, 346 of 2,565 (13.5%) are mixed. Same-sign
tiers counted as the larger of up and down calls are 2,223 / 675 / 225 / 44.
Source: `figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-01_fig3c/`."""),
 dict(file='docs/PAPER.md', line='392', claim_level=False,
  old_text="The per-cohort union contains 13,761 genes, including 9,065 called in one cohort and 51 in all five; it is distinct from the 1,854 pooled calls.",
  new_text="The per-cohort union contains 13,761 genes, including 9,065 called in one cohort and 4,696 in two or more, of which 2,831 have one sign in every calling cohort (44 in all five) and 1,865 have opposite signs in different cohorts; it is distinct from the 1,854 pooled calls."),
 dict(file='docs/PAPER.md', line='392', claim_level=True,
  old_text="A broad established-disease signal recurs across cohorts and holds under a stricter histologic endpoint.",
  new_text="The pooled established-disease signal keeps its direction in each cohort (78.3% to 100% of the 1,854 calls) and holds under a stricter histologic endpoint; shared per-cohort membership is not replication, because 1,865 of 4,696 multi-cohort genes change sign between cohorts."),
 dict(file='docs/STATUS.md', line='301', claim_level=False,
  old_text="The five-cohort union contains 13,761 genes: 9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts.",
  new_text="The five-cohort union contains 13,761 genes: 9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts; 1,865 of the 4,696 multi-cohort genes have opposite signs between cohorts (1,632 involve GSE126848), leaving 1,933 / 636 / 218 / 44 single-sign genes in two to five cohorts."),
 dict(file='docs/STATUS.md', line='814', claim_level=False,
  old_text="Its 6,147 calls join a 13,761-gene per-cohort union (9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts), distinct from pooled membership.",
  new_text="Its 6,147 calls join a 13,761-gene per-cohort union (9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts; 1,865 of the 4,696 multi-cohort genes have opposite signs between cohorts), distinct from pooled membership."),
]
E['R03-02'] = [
 dict(file='docs/manuscript/05_figure_legends.md', line='123-127', claim_level=False,
  old_text="""all male, with constant sex omitted. The restored model yields 6,147 calls,
independently reproduced against native fragments. Its identity screen removed
24 of 31 source disease records, including all ten women, and 2 of 25 controls;
[Inference] unequal eligibility limits representativeness.""",
  new_text="""all male, with constant sex omitted. The restored model yields 6,147 calls,
independently reproduced against native fragments. Its identity screen removed
24 of 31 source disease records, including all ten women, and 2 of 25 controls;
[Inference] unequal eligibility limits representativeness. Retained disease
libraries are about half as deep as controls (median 19.5 versus 35.3 million
fragments; disease deeper in 5 of 161 pairs). 437 of the 6,147 calls are genes
with no reads in any of the 30 libraries, and adding log2 library size to the
model keeps 2,025 of them. The GSE126848 calls therefore mix disease with
sequencing depth."""),
 dict(file='docs/manuscript/draft/fig3.md', line='32-33', claim_level=False,
  old_text="""and 2 of 25 controls. [Inference] Its unequal eligibility limits representativeness.""",
  new_text="""and 2 of 25 controls. [Inference] Its unequal eligibility limits representativeness.
The retained disease libraries are also about half as deep as the controls
(median 19.5 versus 35.3 million fragments; disease deeper in 5 of 161 pairs,
Wilcoxon P=1.9e-5). Because voom adds a fixed count before log-CPM, 437 of the
6,147 calls are genes with zero reads in all 30 libraries, each with log2FC 0.554,
the effective library-size ratio. Adding log2 library size as a covariate keeps
2,025 of the 6,147 calls (2,678 in total); a cohort-specific expression filter
gives 3,440 calls on 15,592 genes. GSE126848 per-cohort calls are therefore not
read as disease replication."""),
 dict(file='docs/RESULTS.md', line='521-523', claim_level=False,
  old_text="""The other four cohort estimates are unchanged. [Inference] Unequal identity
eligibility may affect representativeness; genetic clusters are not verified donors.""",
  new_text="""The other four cohort estimates are unchanged. [Inference] Unequal identity
eligibility may affect representativeness; genetic clusters are not verified donors.
Library depth (review R03-02, 2026-10-15; job 22252035): retained disease
libraries have median 19.48M fragments against 35.30M for controls (disease
deeper in 5 of 161 pairs; exact Wilcoxon P=1.87e-5; r(group, log2 depth)=-0.709).
437 of the 6,147 calls (all up, log2FC 0.554) are genes with zero reads in all
30 libraries. With log2 library size added, the model gives 2,678 calls and keeps
2,025 of the 6,147 (all-gene log2FC Spearman 0.799). A per-cohort filterByExpr
refit gives 3,440 calls (1,821 up, 1,619 down) on 15,592 genes, a different BH
family. Without the 437 zero-read calls the Figure 3C union is 13,447 genes
(8,843 unique to one cohort, 19,617 memberships). Per-cohort quality weights are
0.970 (control) and 0.935 (disease); the 1.14 versus 0.72 contrast is from the
pooled 774-unit fit. Source:
`figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-02_gse126848_depth/`."""),
 dict(file='docs/PAPER.md', line='413-414', claim_level=False,
  old_text="""counts and an independent group-means refit reproduce the result. [Inference]
Unequal eligibility limits representativeness; these clusters are not verified donors.""",
  new_text="""counts and an independent group-means refit reproduce the result. [Inference]
Unequal eligibility limits representativeness; these clusters are not verified donors.
Retained disease libraries are about half as deep as controls; 437 of the 6,147
calls are genes with no reads in any library, and a depth-adjusted model keeps
2,025 of them, so the GSE126848 calls are not read as disease replication."""),
 dict(file='docs/manuscript/METHODS.md', line='144-145', claim_level=False,
  old_text="""source women, and 2 of 25 controls. [Inference] Unequal eligibility limits
representativeness and is not resolved by an estimable design.""",
  new_text="""source women, and 2 of 25 controls. [Inference] Unequal eligibility limits
representativeness and is not resolved by an estimable design. Because retained
disease libraries are shallower than controls, a sensitivity refit added log2
library size to `~ group_binary`, and calls on genes with zero reads in all 30
libraries were counted separately; neither changes the displayed calls."""),
]
E['R03-09'] = [
 dict(file='docs/manuscript/draft/fig3.md', line='38-39', claim_level=False,
  old_text="""correlations of 0.870–0.995 and threshold retention 70.7%–96.6%. Neither
procedure is independent replication of a pooled result.""",
  new_text="""correlations of 0.870–0.995 and threshold retention 70.7%–96.6%. Neither
procedure is independent replication of a pooled result. With the release sample
weights fixed, adding within-cohort standardized log library size, top-100-gene
share or expression entropy keeps 1,760, 1,738 or 1,759 of the 1,854 calls;
adjusting for log library size, top-100 share and detected genes together keeps
1,563 (84.3%), and up calls fall from 1,333 to 1,085. All-gene log2FC Spearman
correlations with the baseline are 0.984 to 0.992. Whether these library
properties are technical or reflect disease biology is not resolved."""),
 dict(file='docs/RESULTS.md', line='120', claim_level=False,
  old_text="""| Strict-depth sensitivity, corrected counts (four-cohort disease contrast: GSE126848 keeps 5 controls and no disease unit) | 738 / 738 | 24,509 | 1,906 (1,344 / 562) | 2,135 |""",
  new_text="""| Strict-depth sensitivity, corrected counts (four-cohort disease contrast: GSE126848 keeps 5 controls and no disease unit) | 738 / 738 | 24,509 | 1,906 (1,344 / 562) | 2,135 |
| Library-complexity sensitivity, release sample weights fixed, + within-cohort z(log library size) (review R03-09) | 774 / 774 | 24,115 | 1,823 (1,351 / 472); 1,760 of 1,854 kept | not run |
| Same, + z(top-100-gene share) | 774 / 774 | 24,115 | 1,815 (1,280 / 535); 1,738 kept | not run |
| Same, + z(expression entropy) | 774 / 774 | 24,115 | 1,814 (1,287 / 527); 1,759 kept | not run |
| Same, + z(log library size) + z(top-100 share) + z(detected genes) | 774 / 774 | 24,115 | 1,605 (1,085 / 520); 1,563 kept (84.3%); all-gene Spearman 0.984 | not run |"""),
]
E['R14-07'] = [
 dict(file='docs/manuscript/05_figure_legends.md', line='158-159', claim_level=False,
  old_text="""3,776 finite native RNA–genetics pairs (Spearman rho=0.022, nominal P=0.18),""",
  new_text="""3,776 finite native RNA–genetics pairs. PP.H4 is unsigned, so its association
is tested against |log2FC| (Spearman rho=0.001, nominal P=0.93),"""),
 dict(file='docs/manuscript/draft/fig3.md', line='139-141', claim_level=False,
  old_text="""log2FC. Its 3,776 finite native RNA–SuSiE pairs have Spearman rho=0.02175
(nominal P=0.181), without a fitted line.""",
  new_text="""log2FC. PP.H4 is unsigned, so it is compared with the magnitude of change:
across 3,776 finite native RNA–SuSiE pairs, Spearman rho between PP.H4 and
|log2FC| is 0.0014 (nominal P=0.93). The 426 genes with PP.H4 >0.5 split 218 up
and 208 down; their median |log2FC| is 0.125 against 0.107 for the rest (nominal
Mann-Whitney P=0.042, uncorrected; genes are not independent units)."""),
 dict(file='docs/manuscript/draft/fig2.md', line='59-61', claim_level=False,
  old_text="""interface, PP.H4 and raw disease-state log2 fold change have Spearman
rho=0.02175 (nominal P=0.18146).""",
  new_text="""interface, PP.H4 and the absolute disease-state log2 fold change have Spearman
rho=0.0014 (nominal P=0.93); PP.H4 is unsigned, so it is not tested against the
signed fold change."""),
 dict(file='docs/RESULTS.md', line='1045-1046', claim_level=False,
  old_text="""testable genes with 426 genetic. Spearman rho is 0.022, nominal P=0.18""",
  new_text="""testable genes with 426 genetic. Spearman rho of PP.H4 with |log2FC| is 0.0014,
nominal P=0.93 (signed log2FC 0.022, P=0.18, is the wrong form for an unsigned
predictor; review R14-07; source
`figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R14-07_fig3g/fig3g_magnitude_statistics.tsv`)"""),
 dict(file='docs/STATUS.md', line='830', claim_level=False,
  old_text="""Both 26; rho 0.0217507, nominal p 0.181458.""",
  new_text="""Both 26; Spearman rho of PP.H4 with |log2FC| 0.00141, nominal p 0.931 (signed rho 0.0217507 withdrawn as the wrong form for an unsigned predictor)."""),
]
E['R16-02'] = [
 dict(file='docs/manuscript/05_figure_legends.md', line='147-149', claim_level=False,
  old_text="""**(F)** Fifty-four genes with BH FDR <0.05 in at least one adjacent fibrosis
contrast and membership in the native 436-gene primary SuSiE set of the corrected
COLOC r2 input.""",
  new_text="""**(F)** Fifty-four genes with BH FDR <0.05 in at least one adjacent fibrosis
contrast and membership in the native 436-gene primary SuSiE set of the corrected
COLOC r2 input. FDR is controlled within each contrast, not over the union of
four contrasts; with the budget split across contrasts (BH q <0.0125 in each),
which bounds the FDR of the union at 5%, 15 of the 54 genes remain."""),
 dict(file='docs/manuscript/draft/fig3.md', line='124-126', claim_level=False,
  old_text="""keeping the measurements separate. Fifty-four genes belong to the native 436-gene
primary SuSiE set and pass BH FDR <0.05 in at least one adjacent fibrosis
contrast.""",
  new_text="""keeping the measurements separate. Fifty-four genes belong to the native 436-gene
primary SuSiE set and pass BH FDR <0.05 in at least one adjacent fibrosis
contrast. Each contrast is its own BH family, so the union carries no joint 5%
FDR; 15 genes pass when the budget is split across the four contrasts (BH
q <0.0125 in each; 11 / 2 / 2 / 1 marks)."""),
 dict(file='docs/RESULTS.md', line='32-33', claim_level=False,
  old_text="""PP.H4 > 0.5 and BH FDR < 0.05 in at least one adjacent fibrosis contrast. Of the""",
  new_text="""PP.H4 > 0.5 and BH FDR < 0.05 in at least one adjacent fibrosis contrast, each
contrast its own 24,115-gene BH family (no joint FDR over the four; at BH
q < 0.0125 per contrast, which bounds the union FDR at 5%, 15 of the 54 genes and
16 marks remain: ACVR2B, BRI3, CDH6, CORO1C, CREB3L4, CTSH, CYP2A6, EPSTI1,
F2RL1, FAM177A1, GCAT, MTARC2, OMG, SLC17A2, TNFSF10; review R16-02). Of the"""),
]
E['R03-11'] = []
E['R14-09'] = [
 dict(file='docs/manuscript/draft/abstract.md', line='38-42', claim_level=True,
  old_text="""A worked GNMT example connects a provisional direct
disease shared-signal locus to a severity-linked GNMT–MAT1A–CYP2C19 state,
defined by coordinated decreases across disease-state RNA, fibrosis-stage
contrasts, and liver protein. This bridge does not imply that all three genes
share one inherited mechanism.""",
  new_text="""GNMT, MAT1A and CYP2C19 form a severity-linked
state, defined by coordinated decreases across disease-state RNA, fibrosis-stage
contrasts, and liver protein. GNMT's link to inherited NAFLD risk rests on one
single-signal colocalization (PP.H4 0.71) in one direct-NAFLD study where
multi-signal colocalization could not be run, so it remains provisional, and
the state does not imply that all three genes share one inherited mechanism."""),
 dict(file='docs/manuscript/draft/fig3.md', line='133-137', claim_level=False,
  old_text="""0.128); GNMT has one direct-NAFLD single-signal ABF result above 0.5 (PP.H4
0.706, 2021_34841290_NAFLD_EUR) and is therefore the provisional bridge between
the genetic and state maps.""",
  new_text="""0.128, from UKBB_ALT); GNMT has one direct-NAFLD single-signal ABF result above
0.5 (PP.H4 0.706, 2021_34841290_NAFLD_EUR). Multi-signal SuSiE was not run for
GNMT in that study (no SuSiE result for 18,939 of its 18,961 gene rows), so
the result is untested by the primary method, not contradicted by it. GNMT is
therefore a provisional bridge between the genetic and state maps."""),
 dict(file='docs/manuscript/draft/fig2.md', line='65-68', claim_level=False,
  old_text="""GRCh37 6:42931261 (rs2296804, SNP.PP.H4 0.258990); its primary multi-signal result remains below threshold
(SuSiE PP.H4 0.128).""",
  new_text="""GRCh37 6:42931261 (rs2296804, SNP.PP.H4 0.258990). SuSiE was not run for GNMT
in that stratum; its best SuSiE PP.H4 in any primary study is 0.128 (UKBB_ALT)."""),
 dict(file='docs/PAPER.md', line='292-294', claim_level=False,
  old_text="""and ABF PP.H4 0.706 (2021_34841290_NAFLD_EUR), so GNMT remains an ABF-supported
provisional bridge;""",
  new_text="""and ABF PP.H4 0.706 (2021_34841290_NAFLD_EUR, where SuSiE was not run for
GNMT), so GNMT remains an ABF-supported provisional bridge;"""),
]
json.dump(E, open('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F3cdg/text_edits.json','w'), indent=1, ensure_ascii=False)
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
for k,v in E.items():
    for e in v:
        s=open(root+e['file']).read(); n=s.count(e['old_text'])
        print(k, e['file'], e['line'], 'matches', n)
