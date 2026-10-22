"""F1 repair text edits. Writes text_edits_repair.json and checks each old_text against the
CURRENT files (exactly one occurrence), overlap among F1 edits, and overlap with other
workstreams' located edits (excluding those F1 supersedes)."""
import json
import os
import re

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HERE = os.path.dirname(os.path.abspath(__file__))
LOC = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/reconcile/all_edits_located.json"
SUPERSEDED = {"E101", "E103", "E104", "E235"} | {f"E{i}" for i in range(109, 126)}

E = []


def add(file, line, old, new, claim=False, tag=""):
    E.append(dict(file=file, line=line, old_text=old, new_text=new, claim_level=claim, tag=tag))


LEG = "docs/manuscript/05_figure_legends.md"
add(LEG, "21",
    "transcription in the native 774-unit interface: 14,704 jointly testable stable-Ensembl genes, 426 genetic, 1,155 canonical disease-state and 26 shared. Canonical RNA calls use BH q<0.05 across 24,115 genes and absolute log2FC>0.5; expected overlap is 33.46 under independence. The native interface retains seven SuSiE-untestable genes as not positive; they are not eligible SuSiE negatives.",
    "transcription in the native 774-unit interface, restricted to the 3,769 genes whose multi-signal COLOC state is positive (SuSiE PP.H4 > 0.5) or tested at PP.H4 ≤ 0.5: 426 colocalized, 220 canonical disease-state and 26 shared. Canonical RNA calls use BH q<0.05 across 24,115 genes and absolute log2FC>0.5. Each bar is a fraction of its own gene list. Dashed lines mark the overlap expected within 10 AveExpr quantile bins (expression deciles), 21.21 genes (4.98% of 426 and 9.64% of 220); Mantel-Haenszel OR 1.295 (95% CI 0.835 to 2.007), exact two-sided stratified P=0.2755 (one prespecified test, no correction). A prespecified list of 426 non-colocalized genes matched one-to-one to the colocalized genes on AveExpr carried 13 DEGs against 21.3 expected within the same bins, so expression is labelled not fully controlled and no direction is read from this comparison. Without stratification 24.87 are expected (5.84% and 11.30%). The 95% CI includes 1; ORs between about 0.46 and 1.77 had less than 80% power to be detected. A further 10,928 genes with a bulk test and a tier-1/2 COLOC record, including 935 canonical DEGs, have only an ABF result and no multi-signal COLOC result; they are not counted as negatives. Seven genes with a SuSiE-untestable pair are also excluded.",
    tag="E109 revised")
add(LEG, "29-30",
    "Catalog. Percentages\nuse jointly testable denominators; technical units are never counted as donors.",
    "Catalog. Percentages in (B) use each gene list as the denominator; technical units are never counted as donors.",
    tag="E110 unchanged")
add(LEG, "160-161",
    "without a fitted line. Colors distinguish Neither (3,156), Transcriptomics sig.\n(194), Genetics sig. (400) and Both (26).",
    "without a fitted line. Colors distinguish Neither (3,156, including seven SuSiE-untestable genes; 3,149 among the 3,769 genes of Figure 1B), Transcriptomics sig.\n(194), Genetics sig. (400) and Both (26).",
    tag="E124 unchanged")

ABS = "docs/manuscript/draft/abstract.md"
add(ABS, "16-17",
    "colocalization. Among 14,704 jointly testable stable gene identifiers, the 426\ncolocalized and 1,155 disease-state genes shared 26 members.",
    "colocalization. Among 3,769 genes evaluated by multi-signal colocalization, the 426\ncolocalized and 220 disease-state genes shared 26, against 21.2 expected within\nexpression deciles (Mantel-Haenszel odds ratio 1.29, 95% CI 0.84 to 2.01).",
    claim=True, tag="E111 revised: OR 1.29, deciles")
add(ABS, "35-38",
    " Some protein-altering variants may\nact through protein function without changing steady-state RNA, while others\nchange RNA through nonsense-mediated decay or splicing; whether this contributes\nto the small overlap is a hypothesis these data do not test.",
    "",
    claim=True, tag="new, conflict 3: delete sentence")

INTRO = "docs/manuscript/draft/01_intro.md"
add(INTRO, "62-68",
    "Connecting the two maps exposes a limit in how inherited liver-disease risk is\nusually mapped. Trait-associated genetic evidence and the established\ndisease-state transcriptome name largely separate gene sets: among 14,704\njointly testable stable gene identifiers, the 426 colocalized and 1,155\ndisease-state genes share 26 members. That overlap is not\nitself informative, because two sets of these sizes are mostly disjoint whether\nor not the underlying biology differs.",
    "The two maps can be compared gene by gene. Among 3,769 genes evaluated by\nmulti-signal colocalization, the 426 colocalized and 220 disease-state genes\nshare 26, against 21.2 expected within expression deciles (Mantel-Haenszel odds\nratio 1.29, 95% CI 0.84 to 2.01). That overlap is not itself informative: the\ninterval includes 1, odds ratios between about 0.46 and 1.77 had less than 80%\npower to be detected, and two sets of these sizes are mostly disjoint whether\nor not the underlying biology differs.",
    claim=True, tag="E112 revised and widened to 62-68 (conflict 3 lead)")
add(INTRO, "87-89",
    "Whether this architecture contributes to\nthe sparse overlap with the disease-state transcriptome is a hypothesis; the\npresent comparison does not test it.",
    "Whether this architecture affects\nwhich genes appear in both the colocalized and disease-state sets is a\nhypothesis; the present comparison does not test it.",
    claim=True, tag="new, conflict 3")

FIG1 = "docs/manuscript/draft/fig1.md"
add(FIG1, "26-33",
    "The organizing interface compares jointly testable inherited genetic and\nestablished disease-state evidence (Figure 1B). In the stable-Ensembl interface\nuniverse of 14,704 jointly testable genes, 426 are SuSiE-colocalized, 1,155\nare canonical disease-state genes and 26 are both, against 33.5 expected under\nindependent sampling. This interface is distinct from the raw categorical\nuniverse of 14,720 and the annotation-adjudicated Catalog universe of 14,682;\ntheir denominators are not interchangeable. The two maps therefore\nintersect at close to the rate chance produces for sets of these sizes;\n",
    "The organizing interface compares inherited genetic and established\ndisease-state evidence (Figure 1B). Among the 3,769 genes evaluated by\nmulti-signal COLOC in the stable-Ensembl interface, 426 are SuSiE-colocalized,\n220 are canonical disease-state genes and 26 are both, against 21.2 expected\nwithin expression deciles (Mantel-Haenszel OR 1.295, 95% CI 0.835 to 2.007;\nexact stratified P=0.2755). A prespecified list of non-colocalized genes\nmatched to the colocalized genes on expression carried 13 disease-state genes\nagainst 21.3 expected within the same deciles, so expression is not fully\ncontrolled and no direction is read from this comparison. The interval\nincludes 1, and odds ratios between about 0.46 and 1.77 had less than 80%\npower to be detected. A further 935 canonical disease-state genes had no\nmulti-signal COLOC result; they are not counted as genetic negatives. This\ninterface is distinct from the raw categorical universe of 14,720 and the\nannotation-adjudicated Catalog universe of 14,682; their denominators are not\ninterchangeable.\n",
    tag="E113 revised")

FIG2 = "docs/manuscript/draft/fig2.md"
add(FIG2, "54-56",
    "release and the adopted 774-unit bulk fit (Figure 1B). In the stable-Ensembl\ninterface universe, 14,704 jointly testable genes include 426 SuSiE-colocalized\nand 1,155 canonical disease-state genes, with 26 shared against 33.5 expected.",
    "release and the adopted 774-unit bulk fit (Figure 1B). Among the 3,769\nstable-Ensembl genes evaluated by multi-signal COLOC, 426 are SuSiE-colocalized\nand 220 canonical disease-state genes, with 26 shared against 21.2 expected\nwithin expression deciles. A further 10,928 jointly tested genes, including 935\ncanonical disease-state genes, have no multi-signal COLOC result and are not\ncounted as negatives.",
    tag="E114 revised")
add(FIG2, "86-87",
    "The distinction between inherited genetic evidence and established tissue\nstate has two implications. First,",
    "The colocalized and disease-state gene sets share 26 genes against 21.2\nexpected within expression deciles (Mantel-Haenszel odds ratio 1.29, 95% CI\n0.84 to 2.01); this is a description, not evidence that the two layers are\nbiologically separate. Two practical points follow. First,",
    tag="new, conflict 4: replaces F2d E235")
add(FIG2, "105",
    "or to the small overlap with the disease-state gene set,",
    "or to which genes appear in both the colocalized and disease-state sets,",
    tag="new, conflict 3 extension")

add("docs/manuscript/draft/fig3.md", "142-144",
    "Genetics sig. and Both contain 3,156 / 194 / 400 / 26 genes. The native jointly\ntested background contains 14,704 genes; it is distinct from the Gene Catalog's\nannotation-adjudicated universe and is not silently replaced by it.",
    "Genetics sig. and Both contain 3,156 / 194 / 400 / 26 genes; Neither includes\nseven SuSiE-untestable genes, so among the 3,769 genes evaluated by multi-signal\nCOLOC (Figure 1B) Neither is 3,149. This universe is distinct from the Gene\nCatalog's annotation-adjudicated universe and is not silently replaced by it.",
    tag="E115 unchanged")

DISC = "docs/manuscript/draft/fig5_discussion.md"
add(DISC, "93-98",
    "The native RNA–genetics interface contains 14,704 genes: 426 genetically\nsupported and 1,155 canonical disease-state calls, with 26 in both compared\nwith 33.46 expected under independent membership (odds ratio 0.757, nominal\nFisher P=0.200). Non-overlap is 93.90%, against 92.14% expected from these set\nsizes. This is an arithmetic comparison, not evidence that the biological\nprocesses are independent.",
    "Among the 3,769 genes of the native RNA–genetics interface evaluated by\nmulti-signal colocalization, 426 are genetically supported and 220 are\ncanonical disease-state calls, with 26 in both against 21.21 expected within\nexpression deciles (Mantel-Haenszel odds ratio 1.295, 95% CI 0.835 to 2.007;\nexact stratified P=0.2755, one prespecified test). Expression was not fully\ncontrolled by the deciles, so no direction is read from this; the interval\nincludes 1, and odds ratios between about 0.46 and 1.77 had less than 80%\npower to be detected. A further 935 disease-state calls had no multi-signal\ncolocalization result and are not genetic negatives. The overlap is a\ndescription, not evidence that the biological processes are independent.",
    claim=True, tag="E116 revised")
add(DISC, "124-128",
    "One hypothesis for the sparse overlap is\nthat expression-based colocalization detects the coding component of direct-trait\nrisk poorly; the present comparison does not test it, because it does not hold\ngene eligibility, phenotype, study power and expression observability fixed.\nThe sparse overlap is not evidence that inherited risk and tissue remodeling",
    "One hypothesis about which genes appear in\nboth sets is that expression-based colocalization detects the coding component\nof direct-trait risk poorly; the present comparison does not test it, because\nit does not hold gene eligibility, phenotype, study power and expression\nobservability fixed. The overlap of 26 genes is not evidence that inherited\nrisk and tissue remodeling",
    claim=True, tag="new, conflict 3")
add(DISC, "164-166",
    "Unequal assay\nobservability could contribute to sparse overlap, but the present comparison\ndoes not isolate that explanation.",
    "Unequal assay\nobservability could affect which genes appear in both sets, but the present\ncomparison does not isolate that explanation.",
    claim=True, tag="new, conflict 3")

PAPER = "docs/PAPER.md"
add(PAPER, "120-121",
    "Their sparse overlap does not establish biological\n   complementarity or antagonism.",
    "Their overlap (26 genes against 21.2\n   expected within expression deciles, among 3,769 genes evaluated by\n   multi-signal COLOC) does not establish biological complementarity or\n   antagonism.",
    claim=True, tag="new, conflict 3")
add(PAPER, "572-574",
    "Any role of the coding difference\n  in the small genetic/disease-state overlap is a hypothesis this comparison did\n  not test.",
    "Any role of the coding difference\n  in which genes are both genetic and disease-state is a hypothesis this\n  comparison did not test.",
    claim=True, tag="new, conflict 3 extension")
add(PAPER, "575-584",
    "- The named genetic and bulk disease-state sets are largely non-overlapping, at\n  approximately the rate independent sampling produces for two sets of their\n  sizes (26 shared against 33.4 expected in the raw symbol-grouped,\n  pre-adjudication universe of 14,720 jointly testable genes,\n  426 genetic and 1,155 disease-state; 93.9% non-overlap against 92.2%\n  expected; odds ratio 0.76, nominal p 0.20). These values use the adopted\n  774-unit bulk refit and COLOC release `susie_coloc_r2_20260924T053244Z`,\n  adopted on October 14. The primary stable-Ensembl interface contains 14,704\n  genes and the annotation-adjudicated Catalog contains 14,682; these universes\n  must not be interchanged.",
    "- Among the 3,769 genes evaluated by multi-signal COLOC (SuSiE state positive\n  or tested at PP.H4 ≤ 0.5), the 426 colocalized and 220 canonical\n  disease-state genes share 26, against 21.21 expected within 10 AveExpr\n  quantile bins (Mantel-Haenszel OR 1.295, 95% CI 0.835 to 2.007; exact\n  two-sided P 0.2755; H57, one prespecified test). These values use the\n  adopted 774-unit bulk refit and COLOC release\n  `susie_coloc_r2_20260924T053244Z`, adopted on October 14. A prespecified\n  list of non-colocalized genes matched on expression carried 13 DEGs against\n  21.3 expected within the same bins, so expression is labelled not fully\n  controlled and no direction may be stated. The 95% CI includes 1; ORs\n  between about 0.46 and 1.77 had less than 80% power to be detected. \"At\n  chance\" may be written only with the interval and this power limit. The\n  earlier values counted genes without a multi-signal COLOC result as genetic\n  negatives and are withdrawn: 26 against 33.46 expected, OR 0.757, in the\n  14,704-gene interface (10,928 such genes, 935 of them canonical DEGs), and\n  26 against 33.4 expected, OR 0.76, in the 14,720 symbol-grouped rows (10,941\n  such genes, 935 of them canonical DEGs). The 14,704-gene stable-Ensembl interface, the 14,720\n  symbol-grouped rows and the 14,682-gene adjudicated Catalog are distinct\n  universes and must not be interchanged.",
    claim=True, tag="E117 revised")

RES = "docs/RESULTS.md"
add(RES, "1981",
    "| Primary SuSiE-nominated genes | 436 | COLOC r2; 7 SuSiE-untestable genes kept in denominators as not positive |",
    "| Primary SuSiE-nominated genes | 436 | COLOC r2; the 7 SuSiE-untestable genes are not positive and are excluded from the overlap background below |",
    tag="new, checker text problem")
add(RES, "1982-1987",
    "| Jointly testable background | 14,720 | Raw symbol-grouped, pre-adjudication categorical universe; distinct from the native stable-Ensembl interface (14,704) and annotation-adjudicated Catalog (14,682) |\n| Jointly testable genetic genes | 426 | SuSiE PP.H4 > 0.5 within the background |\n| Disease-state genes in that background | 1,155 | Canonical threshold, of the 1,854 calls |\n| Overlap with disease-state set | 26 | Against 33.4 expected (426 × 1,155 / 14,720); Fisher OR 0.76, nominal p = 0.20, not below 0.05 |\n| Non-overlap | 93.9% | `(426 − 26) / 426` |\n| Non-overlap expected by chance | 92.2% | `1 − 1,155 / 14,720`. Two sets of these sizes are mostly disjoint whether or not the biology differs; do not present non-overlap as a finding |",
    "| Background evaluated by multi-signal COLOC | 3,769 | Native stable-Ensembl interface genes with SuSiE state positive or tested at PP.H4 ≤ 0.5 (H57 U_eval). Excludes 10,928 genes with only an ABF result (935 canonical DEGs) and 7 SuSiE-untestable genes. Distinct from 14,704 (native interface), 14,720 (symbol-grouped) and 14,682 (Catalog) |\n| Genetic genes | 426 | SuSiE PP.H4 > 0.5 |\n| Disease-state genes in that background | 220 | Canonical threshold; 1,155 in the 14,704 interface |\n| Overlap with disease-state set | 26 | Against 21.21 expected within 10 AveExpr quantile bins; MH OR 1.295 (95% CI 0.835 to 2.007), exact two-sided P 0.2755 (H57, one test). Expression not fully controlled: a list of 426 non-colocalized genes matched on AveExpr carried 13 DEGs against 21.3 expected (O/E 0.610); no direction stated. The CI includes 1; ORs 0.46 to 1.77 had less than 80% power to be detected. Unstratified: 24.87 expected, Fisher OR 1.055, nominal P 0.83 |\n| Non-overlap | 93.9% | `(426 − 26) / 426` |\n| Non-overlap expected without stratification | 94.2% | `1 − 220 / 3,769`. Two sets of these sizes are mostly disjoint whether or not the biology differs; do not present non-overlap as a finding. The withdrawn 14,720-row values were 33.4 expected, OR 0.76, 92.2% |",
    tag="E118 revised")
add(RES, "1989-1991",
    "SuSiE-or-ABF on the same basis: 959 genetic, 57 overlap, 94.1% non-overlap, OR\n0.73, nominal p 0.022 (one of several definitions; no multiplicity correction\napplied).",
    "Rows above are from H57\n(`Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/results.json`);\nthe audit file named above holds only the withdrawn 14,720-row values.\nSuSiE-or-ABF in those 14,720 rows, without expression stratification, gave 959\ngenetic, 57 overlap, 94.1% non-overlap, OR 0.73, nominal p 0.022 (one of several\ndefinitions; no multiplicity correction applied); these values do not use the\n3,769-gene basis of the table. As a diagnostic, the single-signal ABF map (ABF\nPP.H4 ≥ 0.5, 788 genes) in the 14,704-gene interface within 10 AveExpr quantile bins gives 49 shared against\n43.6 expected, MH OR 1.15 (95% CI 0.85 to 1.57), nominal exact P 0.41 (H57 D8).",
    tag="new, conflict 17")
add(RES, "2051-2053",
    "**Do not upgrade the depletion to a finding.** On the 774 units with COLOC r2\n(not adopted) the overlap test gives OR 0.76, nominal p = 0.20 in the\nevidence-class universe.",
    "**Do not report a depletion.** On the 774 units with COLOC r2 (adopted\n2026-10-14), restricted to the 3,769 genes evaluated by multi-signal COLOC and\nstratified by 10 AveExpr quantile bins, no depletion is detected: 26 shared\nagainst 21.21 expected, MH OR 1.295 (95% CI 0.835 to 2.007), exact P 0.2755\n(H57; expression not fully controlled, so no direction is stated). The earlier\nunstratified test gave OR 0.76, nominal p = 0.20 in the 14,720-row\nevidence-class universe (OR 0.757, 33.46 expected, in the 14,704-gene\ninterface); both counted genes without a multi-signal COLOC result as\nnegatives (10,941 and 10,928).",
    tag="E119 revised; absorbs Methods E103 (conflict 25)")
add(RES, "2155-2157",
    "Its stable-ID joint universe contains 14,704 genes: 426 genetic positives,\n1,155 canonical RNA calls and 26 shared, hence 400 genetic-only and 1,129\nstate-only genes.",
    "Its stable-ID joint universe contains 14,704 genes: 426 genetic positives,\n1,155 canonical RNA calls and 26 shared, hence 400 genetic-only and 1,129\nstate-only genes. Of the 1,129 state-only genes, 935 have no multi-signal COLOC\nresult and 194 were tested at PP.H4 ≤ 0.5 (H57 D6); this comparison was not\nrerun on the 194.",
    tag="E125 unchanged")
add(RES, "2221-2222",
    "with COLOC r2, not adopted, 92.2% expected against 93.9% observed for 426, 1,155\nand 14,720).",
    "with COLOC r2 (adopted 2026-10-14), among the 3,769 genes evaluated by\nmulti-signal COLOC, 94.2% expected without stratification against 93.9%\nobserved for 426 and 220).",
    tag="E120 revised; absorbs Methods E104 (conflict 25)")

ST = "docs/STATUS.md"
add(ST, "827",
    "Native joint stable-ID comparison: 14,704 genes, 426 genetic, 1,155 RNA calls, 26 both;",
    "Native joint stable-ID comparison: 14,704 genes, 426 genetic, 1,155 RNA calls, 26 both; Figure 1B overlap uses the 3,769 genes evaluated by multi-signal COLOC (426 genetic, 220 RNA calls, 26 both; 21.21 expected within AveExpr deciles; H57);",
    tag="E122 revised wording")
add(ST, "998-1002",
    "The 2026-09-25 column is on the 844-record bulk fit. On the adopted 774-unit\nrefit with the unadopted COLOC r2 release `susie_coloc_r2_20260924T053244Z`\n(426 / 14,720; 1,155 disease-state genes; 26 overlap against 33.4 expected):\nnon-overlap 93.9% against 92.2% expected, overlap odds ratio 0.76, nominal\np = 0.20.",
    "The 2026-09-25 column is on the 844-record bulk fit. On the adopted 774-unit\nrefit with COLOC release `susie_coloc_r2_20260924T053244Z` (adopted 2026-10-14),\nrestricted to the 3,769 genes evaluated by multi-signal COLOC (H57): 426\ngenetic, 220 disease-state, 26 overlap against 21.21 expected within 10\nAveExpr quantile bins, MH OR 1.295 (95% CI 0.835 to 2.007), exact P 0.2755;\nexpression not fully controlled, no direction stated. The 14,720-row values\n(33.4 expected, OR 0.76) are withdrawn.",
    tag="E121 revised; absorbs Methods E101 (conflict 25)")

add("docs/manuscript/METHODS.md", "593",
    "This interface does not require a finite SuSiE result for every gene.",
    "This interface does not require a finite SuSiE result for every gene. The Figure 1B overlap is computed on the 3,769 genes whose multi-signal COLOC state is positive (SuSiE PP.H4 > 0.5) or tested at PP.H4 ≤ 0.5; the 10,928 genes with only an ABF result and the seven genes with a SuSiE-untestable pair are excluded rather than counted as negatives. The expected overlap is stratified by 10 AveExpr quantile bins (ties broken by Ensembl ID) and tested with the exact convolution of per-bin hypergeometric distributions (two-sided), with the Mantel-Haenszel OR and its Robins-Breslow-Greenland 95% CI (H57 prespecification). Expression control was checked with a prespecified list of 426 non-colocalized genes matched one-to-one to the colocalized genes on AveExpr by greedy nearest-neighbour matching; all 1,000 prespecified visiting orders returned 13 DEGs against 21.3 expected within the same bins, so the result is labelled expression not fully controlled. Power was computed before the outcome from set sizes for an unstratified exact equal-tailed test: ORs between about 0.46 and 1.77 had less than 80% power to be detected.",
    tag="E123 revised")

add("figures/main/fig1_atlas_overview/README.md", "8-10",
    "Figure 1B was also synchronized on 2026-10-14\nto the adopted 774-unit stable-Ensembl interface (14,704 jointly testable genes,\n426 genetic, 1,155 canonical disease-state and 26 shared).",
    "Figure 1B uses the adopted 774-unit fit and\nCOLOC r2, restricted to the 3,769 genes evaluated by multi-signal COLOC (426\ngenetic, 220 canonical disease-state, 26 shared against 21.2 expected within\nAveExpr deciles; H57); 935 further canonical DEGs with no multi-signal COLOC\nresult are given as a count.",
    tag="new, conflict 15")

# ---------------------------------------------------------------- checks
problems = []
spans = {}
for e in E:
    txt = open(os.path.join(ROOT, e["file"])).read()
    n = txt.count(e["old_text"])
    if n != 1:
        problems.append(f"{e['file']} {e['line']}: old_text found {n} times")
        continue
    a = txt.index(e["old_text"]); b = a + len(e["old_text"])
    l0 = txt.count("\n", 0, a) + 1; l1 = txt.count("\n", 0, b - 1 if b > a else b) + 1
    e["computed_lines"] = f"{l0}-{l1}" if l1 != l0 else f"{l0}"
    spans.setdefault(e["file"], []).append((a, b, e["line"]))
    for bad in ["—", "[Inference]", "[Estimate]", "[Unverified]", "[Speculation]"]:
        if bad in e["new_text"]:
            problems.append(f"{e['file']} {e['line']}: new_text contains {bad!r}")
    for w in ["same expression level", "matched expression", "not evaluated", "never-evaluated", "unevaluated",
              "arithmetic comparison", "there is no depletion", "showed that"]:
        if w in e["new_text"]:
            problems.append(f"{e['file']} {e['line']}: new_text contains {w!r}")
    if re.search(r"(?<![\d.])1\.30(?!\d)", e["new_text"]):
        problems.append(f"{e['file']} {e['line']}: new_text prints OR 1.30")
for f, ss in spans.items():
    ss.sort()
    for (a1, b1, l1), (a2, b2, l2) in zip(ss, ss[1:]):
        if a2 < b1:
            problems.append(f"F1 overlap in {f}: {l1} and {l2}")
loc = json.load(open(LOC))
for o in loc:
    if o["id"] in SUPERSEDED or o["ws"] == "F1" or not o.get("spans"):
        continue
    for (a2, b2, _, _) in o["spans"]:
        for (a1, b1, l1) in spans.get(o["file"], []):
            if a2 < b1 and a1 < b2:
                problems.append(f"overlap with {o['id']} ({o['ws']} {o['fid']}) in {o['file']} at F1 line {l1}")
for e in E:
    if e.get("computed_lines") and e["computed_lines"] != e["line"]:
        problems.append(f"line label {e['file']} given {e['line']} computed {e['computed_lines']}")
json.dump(E, open(os.path.join(HERE, "text_edits_repair.json"), "w"), indent=1, ensure_ascii=False)
print(f"{len(E)} edits; problems: {len(problems)}")
for p in problems:
    print(" -", p)
