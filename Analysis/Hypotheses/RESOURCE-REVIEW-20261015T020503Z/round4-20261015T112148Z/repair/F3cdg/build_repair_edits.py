"""F3cdg repair round: revised text edits. Writes repair_text_edits.json and checks each
old_text occurs exactly once in the current file, and that no new_text carries an em dash or
an evidence tag. Run from the repository root."""
import json, re, sys

C = "figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg"
LEG = "docs/manuscript/05_figure_legends.md"
F3 = "docs/manuscript/draft/fig3.md"
F2 = "docs/manuscript/draft/fig2.md"
AB = "docs/manuscript/draft/abstract.md"
RES = "docs/RESULTS.md"
PAP = "docs/PAPER.md"
STA = "docs/STATUS.md"
MET = "docs/manuscript/METHODS.md"

E = {}

# ---------------- R03-01 (Fig 3C direction tiers) ----------------
E["R03-01"] = [
 dict(file=LEG, line="120-122", claim_level=False,
  old_text="family. The union contains 13,761 genes: 9,065 called in one cohort, 3,462 in\ntwo, 922 in three, 261 in four and 51 in all five; ribbon widths represent\n20,054 cohort memberships.",
  new_text="family. The union contains 13,761 genes: 9,065 called in one cohort and 4,696\nin two or more. Of the 4,696, 2,831 have the same sign in every calling cohort\n(1,933 in two, 636 in three, 218 in four and 44 in all five) and 1,865 are\ncalled up in one cohort and down in another (Mixed sign); 1,632 of the 1,865\ninclude a GSE126848 call. Ribbon widths represent 20,054 cohort memberships.\nShared membership is not replication; direction agreement of the pooled calls\nis shown in Figure S3C."),
 dict(file=F3, line="25-27", claim_level=False,
  old_text="the same canonical threshold. Their union contains 13,761 genes: 9,065 called\nin one cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.\nThese are 20,054 cohort memberships, distinct from the 1,854 pooled calls.",
  new_text="the same canonical threshold. Their union contains 13,761 genes: 9,065 called\nin one cohort and 4,696 in two or more. Only 2,831 of the 4,696 have the same\nsign in every cohort that calls them (1,933 in two, 636 in three, 218 in four\nand 44 in all five); 1,865 (39.7%) are called up in one cohort and down in\nanother, and 1,632 of these include a GSE126848 call. Without GSE126848, 346 of\n2,565 multi-cohort genes (13.5%) have mixed signs. These are 20,054 cohort\nmemberships, distinct from the 1,854 pooled calls."),
 dict(file=RES, line="505-506", claim_level=False,
  old_text="13,761 distinct genes and 20,054 cohort memberships: 9,065 genes called in\none cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.",
  new_text="13,761 distinct genes and 20,054 cohort memberships: 9,065 genes called in\none cohort, 3,462 in two, 922 in three, 261 in four and 51 in all five.\nCounted by direction (review R03-01, 2026-10-15), 1,865 of the 4,696\nmulti-cohort genes (39.7%) are called up in one cohort and down in another\n(1,529 / 286 / 43 / 7 in the two- to five-cohort tiers); 2,831 have one sign in\nevery calling cohort (1,933 / 636 / 218 / 44). 1,632 of the 1,865 include a\nGSE126848 call; without GSE126848, 346 of 2,565 (13.5%) are mixed. Same-sign\ntiers counted as the larger of up and down calls are 2,223 / 675 / 225 / 44.\nSource: `figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-01_fig3c/`."),
 dict(file=PAP, line="392", claim_level=False,
  old_text="The per-cohort union contains 13,761 genes, including 9,065 called in one cohort and 51 in all five; it is distinct from the 1,854 pooled calls.",
  new_text="The per-cohort union contains 13,761 genes, including 9,065 called in one cohort and 4,696 in two or more, of which 2,831 have one sign in every calling cohort (44 in all five) and 1,865 have opposite signs in different cohorts; it is distinct from the 1,854 pooled calls."),
 dict(file=PAP, line="392", claim_level=True,
  old_text="A broad established-disease signal recurs across cohorts and holds under a stricter histologic endpoint.",
  new_text="The pooled established-disease signal agrees in direction for 78.3% to 100% of the 1,854 calls within each cohort (each cohort is also part of the pooled fit) and holds under a stricter histologic endpoint; shared per-cohort membership is not replication, because 1,865 of 4,696 multi-cohort genes change sign between cohorts."),
 dict(file=STA, line="301", claim_level=False,
  old_text="The five-cohort union contains 13,761 genes: 9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts.",
  new_text="The five-cohort union contains 13,761 genes: 9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts; 1,865 of the 4,696 multi-cohort genes have opposite signs between cohorts (1,632 involve GSE126848), leaving 1,933 / 636 / 218 / 44 single-sign genes in two to five cohorts."),
 dict(file=STA, line="814", claim_level=False,
  old_text="Its 6,147 calls join a 13,761-gene per-cohort union (9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts), distinct from pooled membership.",
  new_text="Its 6,147 calls join a 13,761-gene per-cohort union (9,065 / 3,462 / 922 / 261 / 51 called in one through five cohorts; 1,865 of the 4,696 multi-cohort genes have opposite signs between cohorts), distinct from pooled membership."),
]

# ---------------- R03-02 (GSE126848 depth) ----------------
E["R03-02"] = [
 dict(file=LEG, line="123-126", claim_level=False,
  old_text="all male, with constant sex omitted. The restored model yields 6,147 calls,\nindependently reproduced against native fragments. Its identity screen removed\n24 of 31 source disease records, including all ten women, and 2 of 25 controls;\n[Inference] unequal eligibility limits representativeness.",
  new_text="all male, with constant sex omitted. The restored model yields 6,147 calls,\nindependently reproduced against native fragments. Its identity screen removed\n24 of 31 source disease records, including all ten women, and 2 of 25 controls.\nUnequal eligibility may limit representativeness. Retained disease libraries\nhave about half the fragment count of controls on the tested genes (median 19.5\nversus 35.3 million; a disease library is deeper in 5 of 161 control-disease\npairs). 437 of the 6,147 calls are genes with no reads in any of the 30\nlibraries. Adding log2 library size to the model keeps 2,025 to 3,056 of the\n6,147 calls, depending on whether library size is taken before or after\nnormalization, and none of the 437. The GSE126848 calls therefore mix disease\nwith sequencing depth."),
 dict(file=F3, line="33", claim_level=False,
  old_text="and 2 of 25 controls. [Inference] Its unequal eligibility limits representativeness.",
  new_text="and 2 of 25 controls. Unequal eligibility may limit representativeness.\nThe retained disease libraries also have about half the fragment count of the\ncontrols on the tested genes (median 19.5 versus 35.3 million; a disease library\nis deeper in 5 of 161 control-disease pairs, exact Wilcoxon P=1.9e-5). Because\nvoom adds 0.5 to every count before computing log-CPM, a gene with no reads has\nhigher log-CPM in shallower libraries; 437 of the 6,147 calls are such genes, with\nzero reads in all 30 libraries. Each has log2FC 0.554, which equals the\nprecision-weighted control-minus-disease difference in mean log2 normalized\nlibrary size. Adding log2 library size as a covariate keeps 2,025 of the 6,147\ncalls (2,678 in total) when library size is counted before normalization and\n3,056 (3,649 in total) when it is counted after; neither keeps any of the 437.\nA cohort-specific expression filter gives 3,440 calls on 15,592 genes. GSE126848\nper-cohort calls are therefore not read as disease replication."),
 dict(file=RES, line="521-522", claim_level=False,
  old_text="The other four cohort estimates are unchanged. [Inference] Unequal identity\neligibility may affect representativeness; genetic clusters are not verified donors.",
  new_text="The other four cohort estimates are unchanged. Unequal identity eligibility\nmay limit representativeness; genetic clusters are not verified donors.\nLibrary depth (review R03-02, 2026-10-15; jobs 22252035 and 22256115): retained\ndisease libraries have median 19.48M fragments counted on the 24,115 tested\ngenes against 35.30M for controls (disease deeper in 5 of 161 pairs; exact\nWilcoxon P=1.87e-5; r(group, log2 depth)=-0.709); after RLE normalization the\nmedians are 13.85M and 20.40M. 437 of the 6,147 calls (all up) are genes with\nzero reads in all 30 libraries; each has log2FC 0.554, the precision-weighted\ncontrol-minus-disease difference in mean log2(normalized library size + 1).\nWith log2 library size added, the model gives 2,678 calls and keeps 2,025 of\nthe 6,147 (all-gene log2FC Spearman 0.799); with log2 normalized library size\nit gives 3,649 and keeps 3,056 (Spearman 0.867). Neither keeps any of the 437.\nA per-cohort filterByExpr refit gives 3,440 calls (1,821 up, 1,619 down) on\n15,592 genes, a different BH family. Without the 437 zero-read calls the\nFigure 3C union is 13,447 genes (8,843 unique to one cohort, 19,617\nmemberships). Median per-cohort quality weights are 0.970 (control) and 0.935\n(disease); the 1.14 versus 0.72 medians are from the pooled 774-unit fit.\nDepth sources: `" + C + "/R03-02_gse126848_depth/`\nand `" + C + "/repair/repair_libdef_log.txt`."),
 dict(file=PAP, line="414-415", claim_level=False,
  old_text="counts and an independent group-means refit reproduce the result. [Inference]\nUnequal eligibility limits representativeness; these clusters are not verified donors.",
  new_text="counts and an independent group-means refit reproduce the result. Unequal\neligibility may limit representativeness; these clusters are not verified donors.\nRetained disease libraries are about half as deep as controls; 437 of the 6,147\ncalls are genes with no reads in any library, and adding log2 library size keeps\n2,025 to 3,056 of the 6,147 calls (library size before or after normalization)\nand none of the 437, so the GSE126848 calls are not read as disease replication."),
 dict(file=MET, line="144-145", claim_level=False,
  old_text="source women, and 2 of 25 controls. [Inference] Unequal eligibility limits\nrepresentativeness and is not resolved by an estimable design.",
  new_text="source women, and 2 of 25 controls. Unequal eligibility may limit\nrepresentativeness and is not resolved by an estimable design. Retained disease\nlibraries are shallower than controls. Figure 3C shows the unadjusted\n`~ group_binary` calls. Two sensitivities are reported in the text and not used\nfor the panel: refits of `~ group_binary + log2(library size)` with library size\ncounted before or after RLE normalization, and the count of calls on genes with\nzero reads in all 30 libraries."),
]

# ---------------- R03-09 (pooled library-complexity sensitivity) ----------------
E["R03-09"] = [
 dict(file=F3, line="39-40", claim_level=False,
  old_text="correlations of 0.870–0.995 and threshold retention 70.7%–96.6%. Neither\nprocedure is independent replication of a pooled result.",
  new_text="correlations of 0.870–0.995 and threshold retention 70.7%–96.6%. Neither\nprocedure is independent replication of a pooled result. With the release sample\nweights fixed, adding within-cohort standardized log library size keeps 1,760 of\nthe 1,854 calls when library size is the total count on the tested genes and\n1,648 when it is taken after RLE normalization; top-100-gene share or expression\nentropy keeps 1,738 or 1,759. Adjusting for log library size, top-100 share and\ndetected genes together keeps 1,563 (84.3%) or 1,569 (84.6%) under the two\nlibrary-size definitions, and up calls fall from 1,333 to 1,085 or 1,100.\nAll-gene log2FC Spearman correlations with the baseline are 0.979 to 0.992.\nWhether these library properties are technical or reflect disease biology is\nnot resolved."),
 dict(file=RES, line="120", claim_level=False,
  old_text="| Strict-depth sensitivity, corrected counts (four-cohort disease contrast: GSE126848 keeps 5 controls and no disease unit) | 738 / 738 | 24,509 | 1,906 (1,344 / 562) | 2,135 |",
  new_text="| Strict-depth sensitivity, corrected counts (four-cohort disease contrast: GSE126848 keeps 5 controls and no disease unit) | 738 / 738 | 24,509 | 1,906 (1,344 / 562) | 2,135 |\n| Library-complexity sensitivity, release sample weights fixed, + within-cohort z(log library size; total counts on the tested genes) (review R03-09) | 774 / 774 | 24,115 | 1,823 (1,351 / 472); 1,760 of 1,854 kept | not run |\n| Same, + z(log normalized library size; library size x RLE factor) | 774 / 774 | 24,115 | 1,715 (1,289 / 426); 1,648 kept | not run |\n| Same, + z(top-100-gene share) | 774 / 774 | 24,115 | 1,815 (1,280 / 535); 1,738 kept | not run |\n| Same, + z(expression entropy) | 774 / 774 | 24,115 | 1,814 (1,287 / 527); 1,759 kept | not run |\n| Same, + z(log library size) + z(top-100 share) + z(detected genes) | 774 / 774 | 24,115 | 1,605 (1,085 / 520); 1,563 kept (84.3%); all-gene Spearman 0.984 | not run |\n| Same, + z(log normalized library size) + z(top-100 share) + z(detected genes) | 774 / 774 | 24,115 | 1,590 (1,100 / 490); 1,569 kept (84.6%); all-gene Spearman 0.988 | not run |"),
]

# ---------------- R14-07 (Fig 3G statistic) ----------------
E["R14-07"] = [
 dict(file=LEG, line="159", claim_level=False,
  old_text="3,776 finite native RNA–genetics pairs (Spearman rho=0.022, nominal P=0.18),",
  new_text="3,776 finite native RNA–genetics pairs. Because PP.H4 is unsigned, its\nassociation is tested against |log2FC| (Spearman rho=0.001, nominal P=0.93).\nPoints are plotted"),
 dict(file=F3, line="140-141", claim_level=False,
  old_text="log2FC. Its 3,776 finite native RNA–SuSiE pairs have Spearman rho=0.02175\n(nominal P=0.181), without a fitted line.",
  new_text="log2FC. Because PP.H4 is unsigned, it is compared with the magnitude of change:\nacross 3,776 finite native RNA–SuSiE pairs, Spearman rho between PP.H4 and\n|log2FC| is 0.0014 (nominal P=0.93), shown without a fitted line. The 426 genes\nwith PP.H4 >0.5 split 218 up and 208 down."),
 dict(file=F2, line="60-61", claim_level=False,
  old_text="interface, PP.H4 and raw disease-state log2 fold change have Spearman\nrho=0.02175 (nominal P=0.18146).",
  new_text="interface, PP.H4 and the absolute disease-state log2 fold change have Spearman\nrho=0.0014 (nominal P=0.93); PP.H4 is unsigned, so it is not tested against the\nsigned fold change."),
 dict(file=RES, line="1045-1046", claim_level=False,
  old_text="testable genes with 426 genetic. Spearman rho is 0.022, nominal P=0.18\n(`fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_matched_treat_summary.tsv`).",
  new_text="testable genes with 426 genetic. Because PP.H4 is unsigned, the test is its\nSpearman rho with |log2FC|: 0.0014, nominal P=0.93 (0.0010, P=0.95 among the\n3,769 genes of the H57 universe; review R14-07;\n`" + C + "/R14-07_fig3g/fig3g_magnitude_statistics.tsv`\nand `" + C + "/repair/R14-07_fig3g_universes.tsv`). The signed-log2FC rho 0.022\n(nominal P=0.18;\n`fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_matched_treat_summary.tsv`)\nis withdrawn as the wrong form for an unsigned predictor."),
 dict(file=STA, line="830", claim_level=False,
  old_text="Both 26; rho 0.0217507, nominal p 0.181458.",
  new_text="Both 26; Spearman rho of PP.H4 with \\|log2FC\\| 0.00141, nominal p 0.931 (signed rho 0.0217507 withdrawn as the wrong form for an unsigned predictor)."),
]

# ---------------- R16-02 (Fig 3F family) ----------------
E["R16-02"] = [
 dict(file=LEG, line="147-149", claim_level=False,
  old_text="**(F)** Fifty-four genes with BH FDR <0.05 in at least one adjacent fibrosis\ncontrast and membership in the native 436-gene primary SuSiE set of the corrected\nCOLOC r2 input.",
  new_text="**(F)** Fifty-four genes with BH FDR <0.05 in at least one adjacent fibrosis\ncontrast and membership in the native 436-gene primary SuSiE set of the corrected\nCOLOC r2 input. FDR is controlled within each contrast, not over the union of\nfour contrasts. With the 5% budget split across contrasts (BH q <0.0125 in each),\nwhich bounds the FDR over all discoveries in the four 24,115-gene families at 5%\nbut not specifically within the SuSiE subset, 15 of the 54 genes remain."),
 dict(file=F3, line="124-126", claim_level=False,
  old_text="keeping the measurements separate. Fifty-four genes belong to the native 436-gene\nprimary SuSiE set and pass BH FDR <0.05 in at least one adjacent fibrosis\ncontrast.",
  new_text="keeping the measurements separate. Fifty-four genes belong to the native 436-gene\nprimary SuSiE set and pass BH FDR <0.05 in at least one adjacent fibrosis\ncontrast. Each contrast is its own BH family, so the union carries no joint 5%\nFDR; 15 genes pass when the budget is split across the four contrasts (BH\nq <0.0125 in each; 11 / 2 / 2 / 1 marks)."),
 dict(file=RES, line="33", claim_level=False,
  old_text="PP.H4 > 0.5 and BH FDR < 0.05 in at least one adjacent fibrosis contrast. Of the",
  new_text="PP.H4 > 0.5 and BH FDR < 0.05 in at least one adjacent fibrosis contrast, each\ncontrast its own 24,115-gene BH family (no joint FDR over the four). At BH\nq < 0.0125 per contrast, which bounds the FDR over all discoveries in the four\nfamilies at 5% but not specifically within the SuSiE subset, 15 of the 54 genes\nand 16 marks remain: ACVR2B, BRI3, CDH6, CORO1C, CREB3L4, CTSH, CYP2A6, EPSTI1,\nF2RL1, FAM177A1, GCAT, MTARC2, OMG, SLC17A2, TNFSF10 (review R16-02). Of the"),
]

# ---------------- R14-09 (GNMT bridge) ----------------
E["R14-09"] = [
 dict(file=AB, line="38-42", claim_level=True,
  old_text="A worked GNMT example connects a provisional direct\ndisease shared-signal locus to a severity-linked GNMT–MAT1A–CYP2C19 state,\ndefined by coordinated decreases across disease-state RNA, fibrosis-stage\ncontrasts, and liver protein. This bridge does not imply that all three genes\nshare one inherited mechanism.",
  new_text="GNMT, MAT1A and CYP2C19 form a severity-linked\nstate, defined by coordinated decreases across disease-state RNA, fibrosis-stage\ncontrasts, and liver protein. GNMT's colocalization with inherited NAFLD risk\nrests on one single-signal result (PP.H4 0.71) in one direct-NAFLD study, where\nmulti-signal colocalization was not run for GNMT, so it remains provisional, and\nthe state does not imply that all three genes share one inherited mechanism."),
 dict(file=F3, line="134-136", claim_level=False,
  old_text="0.128); GNMT has one direct-NAFLD single-signal ABF result above 0.5 (PP.H4\n0.706, 2021_34841290_NAFLD_EUR) and is therefore the provisional bridge between\nthe genetic and state maps.",
  new_text="0.128, from UKBB_ALT); GNMT has one direct-NAFLD single-signal ABF result above\n0.5 (PP.H4 0.706, 2021_34841290_NAFLD_EUR). Multi-signal SuSiE was not run for\nGNMT in that study (no SuSiE result for 18,939 of its 18,961 gene rows), so\nthe result is untested by the primary method, not contradicted by it. GNMT is\ntherefore a provisional bridge between the genetic and state maps."),
 dict(file=F2, line="67-68", claim_level=False,
  old_text="GRCh37 6:42931261 (rs2296804, SNP.PP.H4 0.258990); its primary multi-signal result remains below threshold\n(SuSiE PP.H4 0.128).",
  new_text="GRCh37 6:42931261 (rs2296804, SNP.PP.H4 0.258990). SuSiE was not run for GNMT\nin that stratum; its best SuSiE PP.H4 in any primary study is 0.128 (UKBB_ALT)."),
 dict(file=PAP, line="293-294", claim_level=False,
  old_text="and ABF PP.H4 0.706 (2021_34841290_NAFLD_EUR), so GNMT remains an ABF-supported\nprovisional bridge;",
  new_text="and ABF PP.H4 0.706 (2021_34841290_NAFLD_EUR, where SuSiE was not run for\nGNMT), so GNMT remains an ABF-supported provisional bridge;"),
 dict(file=PAP, line="392", claim_level=True,
  old_text="connects the GNMT locus to established disease without turning the three-gene state into one inherited mechanism",
  new_text="connects the GNMT locus to established disease only provisionally (one single-signal ABF result, PP.H4 0.706, in a study where SuSiE was not run for GNMT) and without turning the three-gene state into one inherited mechanism"),
]

E["R03-11"] = []

ok = True
for fid, eds in E.items():
    for t in eds:
        s = open(t["file"]).read()
        n = s.count(t["old_text"])
        i = s.find(t["old_text"])
        start = s[:i].count("\n") + 1
        end = start + t["old_text"].count("\n")
        rng = str(start) if start == end else f"{start}-{end}"
        if n != 1 or rng != t["line"]:
            ok = False
        bad = [x for x in ["—", "[Inference]", "[Unverified]", "[Estimate]", "[Speculation]"] if x in t["new_text"]]
        if bad:
            ok = False
        print(fid, t["file"], "count", n, "line", t["line"], "actual", rng, "bad", bad)
json.dump(E, open(sys.argv[1], "w"), indent=1, ensure_ascii=False)
print("ALL OK" if ok else "PROBLEM")
