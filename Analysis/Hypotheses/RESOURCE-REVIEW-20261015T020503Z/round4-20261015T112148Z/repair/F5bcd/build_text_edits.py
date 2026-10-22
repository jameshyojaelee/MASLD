"""F5bcd repair round: assemble text edits and check each old_text against the CURRENT files.
Writes text_edits_repair.json beside this script. Read-only on docs/ and figures/main/."""
import json, os, re, sys
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HERE = os.path.dirname(os.path.abspath(__file__))
REP = "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair"
LEG = "docs/manuscript/05_figure_legends.md"
F5 = "docs/manuscript/draft/fig5.md"
DISC = "docs/manuscript/draft/fig5_discussion.md"
RES = "docs/RESULTS.md"
STA = "docs/STATUS.md"
MET = "docs/manuscript/METHODS.md"
README = "figures/main/fig5_molecular_context/README.md"

E = {}
# ---------------- R06-06 ----------------
E["R06-06"] = [
 dict(file=LEG, line="384",
  old_text="with mRNA; 2,843 are not protein-significant. Spearman rho is 0.457835 and",
  new_text="with mRNA; 2,843 are not protein-significant. Sign agreement depends on the RNA call: 170 of 177 protein-significant genes with RNA BH padj <0.05 agree in sign (96.0%), against 18 of 34 (52.9%; binomial p=0.86 against one half) when RNA padj is 0.05 or higher. Spearman rho is 0.457835 and"),
 dict(file=LEG, line="385",
  old_text="direction agreement is 65.9136%, rounded to 66% on the panel.",
  new_text="direction agreement is 65.9136%, rounded to 66% on the panel. The display universe is partly selected on RNA: its transcriptomic arm (2,390 of 3,054 genes) was selected from RNA-based analyses, including differential-expression contrasts of the retired 846-sample bulk fit of the same five cohorts. Rho is 0.487 in that arm, 0.305 in the 664 genes selected on genetic evidence alone and 0.346 across all 5,836 measured proteins."),
 dict(file=F5, line="13-15",
  old_text="liver proteome (58 patients). Of the genes reaching protein significance, most\nchange in the same direction as their transcript; nonsignificance is not\ninterpreted as falsification.",
  new_text="liver proteome (58 patients). Protein-significant genes agree in sign with\ntheir transcript when the RNA change is itself significant (170 of 177), but\nonly at chance level when it is not (18 of 34); nonsignificance is not\ninterpreted as falsification. RNA-protein Spearman rho is 0.346 across all\n5,836 measured proteins. In the display universe it is 0.487 in the 2,390 genes\nselected from RNA-based analyses of the same cohorts and 0.305 in the 664 genes\nselected on genetic evidence alone."),
 dict(file=RES, line="3064",
  old_text="65.9136%, rounded to 66% on the panel. OLS slope=0.455025,",
  new_text="65.9136%, rounded to 66% on the panel (rho 0.346 across all 5,836 measured proteins, 0.487 in the 2,390 transcriptomic-arm genes and 0.305 in the 664 genetic-only genes; the transcriptomic arm was selected from RNA-based analyses that include the retired 846-sample bulk contrasts). Among protein-significant genes, sign agreement is 170/177 when RNA BH padj <0.05 and 18/34 otherwise. OLS slope=0.455025,"),
]
# ---------------- R06-08 (includes conflict 20/21 merges where noted) ----------------
E["R06-08"] = [
 dict(file=LEG, line="403-404",
  old_text="MASLD/liver-fat pairs (n=28) from liver-enzyme pairs (n=739) in the adopted\nOctober 14 r2 PDF. The historical panel had 29 and 787 pairs.",
  new_text="MASLD/liver-fat genes (n=22; 28 gene-study pairs) from liver-enzyme genes\n(n=424; 739 pairs) in the adopted October 14 r2 release; ten genes appear in both\nclasses. Each gene takes its best state over its study pairs, so genes with many\nliver-enzyme GWAS are counted once. The historical panel counted 29 and 787 pairs."),
 dict(file=F5, line="35-39",
  old_text="Figure 5D intersects the regulatory posterior of 767 colocalized gene-study\npairs with open chromatin from two liver snATAC sources across five shared\nlineages, keeping direct MASLD/liver-fat and liver-enzyme traits separate. Peak\noverlap localizes regulatory DNA; it does not assign a variant to a neighboring\ngene, and absent overlap is unresolved context rather than a negative test.",
  new_text="Figure 5D intersects the colocalization posterior of 436 genes (22 with a\ndirect MASLD/liver-fat pair and 424 with a liver-enzyme pair, ten in both; 767\ngene-study pairs) with open chromatin from two liver snATAC sources across five\nshared lineages. For 12 of 22 direct-trait genes and 167 of 424 liver-enzyme\ngenes, less than 5% of the posterior lies in peaks in every lineage. Averaged\nover liver-enzyme genes, posterior mass in peaks is 2.5 to 3.1 times the\nexpectation from the same window's tested variants in each lineage, but in no\nlineage do more than half of these genes exceed their own expectation (94 to\n212 of 424); for the 22 direct-trait genes the 95% interval of the ratio\nincludes 1. Overlap does not assign a variant to a neighboring gene, and absent\noverlap is unresolved context rather than a negative test.\n\n"
           "Figure 5E places all 117 frozen programs in two assay-native contexts without\nselecting them by disease association. Twenty-eight programs contain at least\none gene from an adopted colocalized pair that places at least half of its\nposterior in open chromatin of the program's own lineage. In three of them\n(hepatocyte Ion/membrane potential, cholangiocyte FXR/bile-acid and fibroblast\nElastic-fiber ECM) the gene is RORA, open through its direct MASLD pair; in the\nother 25, open genes are open only through liver-enzyme pairs. RORA is also open\nin the T/NK Memory/MAIT program, but only through liver-enzyme pairs; its direct\npair places 0.37 of its posterior in T/NK peaks of one source, below the 0.5 rule. In liver protein from 46 MASLD\nparticipants, 26 of 94 scorable programs are robustly associated with a\nhistologic-burden axis across participants (BH over 94): 16 negatively,\nincluding oxidative-phosphorylation, xenobiotic and amino-acid degradation\nprograms, and 10 positively, including fibroblast ECM proteoglycan and\nTGFB/collagen programs. The same protein cohort did not contribute to program\ndiscovery. Neither spatial program in Figure 5F contains a colocalized gene or\ncan be scored in protein."),
 dict(file=RES, line="2940",
  old_text="28 direct MASLD/liver-fat pairs from 739 liver-enzyme pairs (767 total), read",
  new_text="28 direct MASLD/liver-fat pairs from 739 liver-enzyme pairs (767 total; 22 and 424 genes, ten genes in both classes; the review-round-3 candidate counts genes), read"),
 dict(file=STA, line="829",
  old_text="5D uses 767 r2 gene-study pairs (28 direct, 739 enzyme).",
  new_text="5D uses 767 r2 gene-study pairs (28 direct, 739 enzyme; 22 and 424 genes). A gene-level candidate with a 5% posterior floor is in figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/."),
 dict(file=README, line="28-29",
  old_text="colocalized gene-study pairs (436 genes; 28 direct MASLD/liver-fat and 739\nliver-enzyme pairs), then sums normalized SNP.PP.H4",
  new_text="colocalized gene-study pairs (436 genes; 28 direct MASLD/liver-fat and 739\nliver-enzyme pairs, or 22 and 424 genes with ten in both), then sums normalized SNP.PP.H4"),
 dict(file=README, line="33-36",
  old_text="traits remain separate, and\nboth-source, one-source, partial, no-overlap, and untestable states remain\nexplicit. Peak overlap localizes regulatory DNA; it does not assign a target\ngene or establish functional regulation.",
  new_text="traits remain separate. Bars count genes, each at its best state over its study\npairs: at least 50% of the posterior in peaks of both sources, at least 50% in\none source, 5% to under 50%, and under 5% in both sources (including no\noverlap); the untestable state is empty and stays in the source table. Peak\noverlap is regulatory-DNA context; it does not assign a target gene or establish\nfunctional regulation."),
]
# ---------------- R11-02 ----------------
E["R11-02"] = [
 dict(file=LEG, line="405",
  old_text="both-source, one-source, partial, and no-overlap states; all pairs met the 0.95",
  new_text="states with at least 50% of the posterior in peaks of both sources, at least 50% in one source, 5% to under 50% in the better source, and under 5% in both sources (including no overlap); the 0.05 floor keeps trace overlap out of the partial state. Segments under 6.5% of genes are not labeled; all counts are in the source table. All pairs met the 0.95"),
 dict(file=LEG, line="407-408",
  old_text="source table. Peak overlap localizes regulatory DNA but does not assign a\ntarget gene or establish functional regulation.",
  new_text="source table. The locus-matched expectation is the fraction of a pair's mapped\ntested variants that lie in a lineage peak of either source. Averaged over the\n424 liver-enzyme genes, posterior mass in peaks is 2.5 to 3.1 times this\nexpectation in each lineage (gene-bootstrap 95% intervals above 1); for the 22\ndirect-trait genes the intervals include 1. In each lineage at most half of the\nliver-enzyme genes (94 to 212 of 424) exceed their own expectation, so the\naverage excess comes from a subset of genes. Peak overlap is regulatory-DNA\ncontext; it does not assign a target gene or establish functional regulation."),
 dict(file=RES, line="2946",
  old_text="are partial rather than both-source accessible. Peak overlap is regulatory-DNA",
  new_text="are partial rather than both-source accessible, and more than half of the 3,409 partial pair-lineage rows (1,943) place under 1% of the posterior in peaks of either source (239 of 767 pairs, 13 of 28 direct, are under 1% in every lineage). Against a locus-matched expectation (the fraction of each pair's tested variants in peaks), liver-enzyme genes carry 2.5 to 3.1 times the expected mass on average in every lineage (n=424 genes, gene-bootstrap intervals above 1), but in each lineage at most half of them (94 to 212 of 424) exceed their own expectation; for the 22 direct-trait genes the intervals include 1 (review-round-3 candidate, F5bcd). Peak overlap is regulatory-DNA"),
 dict(file=MET, line="1121-1122",
  old_text="separate. Partial or absent peak overlap is unresolved context, not a negative\nfunctional test.",
  new_text="separate, and counts use the gene as the unit (best state over a gene's study\npairs). Overlap below 5% of the posterior in both sources is grouped with absent\noverlap. The locus-matched expectation for a pair and lineage is the fraction of\nits mapped tested variants inside a peak of that lineage in either source;\nobserved and expected masses are averaged per gene and compared by 2,000 gene\nbootstrap resamples (seed 20260923). This expectation matches the window, its\nvariant density and its peak coverage, but not LD or allele frequency. Partial\nor absent peak overlap is unresolved context, not a negative\nfunctional test. For Figure 5E, a program member gene counts as open when one of\nits primary pairs places at least half of its posterior in open chromatin of the\nprogram's own lineage in at least one source (states `replicated_accessible` or\n`source_dependent`); hepatocyte, cholangiocyte, macrophage, fibroblast and T-cell\nprograms map to hepatocyte, cholangiocyte, macrophage, stellate and T/NK\nchromatin. The count uses the same 767 COLOC r2 pairs that Figure 5D summarizes\nper gene (22 direct, 424 enzyme), including the three MVP European-ancestry\nstudies."),
]
# ---------------- R06-09 ----------------
E["R06-09"] = [
 dict(file=DISC, line="36-39",
  old_text="GNMT, MAT1A, and CYP2C19\ncarry concordant disease-state RNA, cross-sectional stage, and tissue-protein\nevidence, routing the hepatocyte metabolic-function hypothesis to\ncontext-matched perturbation and rescue.",
  new_text="GNMT and CYP2C19 are lower in disease at RNA and protein levels and are\nestablished-state entries in the Gene Catalog, which routes them to\ncontext-matched perturbation and rescue. MAT1A is also lower at both levels, but\nits pooled RNA effect (log2FC -0.485) is below the 0.50 floor of a canonical\ncall, so its Catalog class is unresolved rather than established-state."),
]
# ---------------- R11-07 ----------------
E["R11-07"] = [
 dict(file=LEG, line="379-381",
  old_text="caption. GNMT, MAT1A, and CYP2C19 are lower at both RNA and protein levels and\nare interpreted as severity-linked reduction in hepatocyte metabolic gene and\nprotein abundance.",
  new_text="caption. GNMT, MAT1A, and CYP2C19 are lower in MASLD than in control liver at\nboth RNA and protein levels. The protein model compares 46 MASLD with 12 control\nparticipants, adjusted for acquisition batch, age, BMI and sex, and has no\nseverity or hepatocyte-fraction term. In a post hoc analysis that was not\nprespecified, protein abundance of all three is negatively associated with NAS\nacross the 46 MASLD participants (partial Spearman rho -0.49 to -0.66, same\ncovariates; BH q <=0.002 over 3 proteins and 2 histology scores), and with\nfibrosis stage for GNMT and MAT1A (rho -0.33 and -0.55, q 0.040 and 0.0004;\nCYP2C19 rho -0.27, q 0.084). Lower bulk abundance can reflect a smaller\nhepatocyte share of the biopsy as well as lower per-cell abundance, and\nmetabolic function is not measured."),
 dict(file=F5, line="23-25",
  old_text="as pooled RNA. This is assay-native support for a severity-linked reduction\nin hepatocyte metabolic gene and protein abundance, not evidence that the three genes share one\ninherited regulatory mechanism.",
  new_text="as pooled RNA. The protein model compares 46 MASLD with 12 control participants\nand has no severity or hepatocyte-fraction term, so it shows lower bulk protein\nin disease. In a post hoc check that was not prespecified, all three proteins\nare lower at higher NAS across the 46 MASLD participants (partial Spearman rho\n-0.49 to -0.66). Lower bulk abundance can reflect fewer hepatocytes in the\nbiopsy as well as lower per-cell abundance. This is not evidence that the three\ngenes share one inherited regulatory mechanism."),
 dict(file=RES, line="553",
  old_text="### Loss of hepatocyte metabolic function: GNMT–MAT1A–CYP2C19",
  new_text="### Lower hepatocyte metabolic gene and protein abundance: GNMT–MAT1A–CYP2C19"),
 dict(file=RES, line="576-577",
  old_text="a recurring example of reduced abundance of hepatocyte metabolic genes and proteins in\nestablished disease. It does not establish",
  new_text="a recurring example of reduced abundance of hepatocyte metabolic genes and proteins in\nestablished disease. The protein model compares 46 MASLD with 12 control\nparticipants and has no severity term. In a post hoc check that was not\nprespecified (review round 3, F5bcd), protein abundance within the 46 MASLD\nparticipants is negatively associated with NAS (rank-residual partial Spearman,\nbatch, age, BMI and sex: GNMT -0.486, MAT1A -0.660, CYP2C19 -0.514; BH q 0.0017,\n1.2e-5 and 0.0010 over 3 proteins x 2 histology scores) and with fibrosis stage\n(-0.329, -0.554 and -0.270; q 0.040, 4.3e-4 and 0.084). Lower bulk abundance\ncan reflect a smaller hepatocyte share of the biopsy as well as lower per-cell\nabundance. It does not establish"),
 dict(file=LEG, line="399",
  old_text="rows are descriptive rather than independent validation.",
  new_text="rows are descriptive rather than independent validation. The histology column\nshows partial Spearman rho over all 58 participants, including 12 controls with\nfibrosis 0-1 and NAS 0-1, so it combines disease status with severity."),
]
E["R06-12"] = []

problems = []
for fid, eds in E.items():
    for e in eds:
        txt = open(os.path.join(ROOT, e["file"]), encoding="utf-8").read()
        n = txt.count(e["old_text"])
        if n != 1:
            problems.append(f"{fid} {e['file']}:{e['line']} old_text count={n}")
        else:
            start = txt[:txt.index(e["old_text"])].count("\n") + 1
            e["checked_start_line"] = start
        if "—" in e["new_text"]:
            problems.append(f"{fid} em dash in new_text {e['file']}")
        if re.search(r"\[(Inference|Unverified|Speculation|Estimate)", e["new_text"]):
            problems.append(f"{fid} evidence tag in new_text {e['file']}")
        e["claim_level"] = False
# same-file overlap check (spans must not overlap within a file)
spans = {}
for fid, eds in E.items():
    for e in eds:
        if "checked_start_line" not in e: continue
        txt = open(os.path.join(ROOT, e["file"]), encoding="utf-8").read()
        i = txt.index(e["old_text"]); spans.setdefault(e["file"], []).append((i, i + len(e["old_text"]), fid, e["line"]))
for f, s in spans.items():
    s.sort()
    for a, b in zip(s, s[1:]):
        if b[0] < a[1]:
            problems.append(f"overlap in {f}: {a[2]}:{a[3]} and {b[2]}:{b[3]}")
json.dump(E, open(os.path.join(HERE, "text_edits_repair.json"), "w"), indent=1, ensure_ascii=False)
for fid, eds in E.items():
    for e in eds:
        print(fid, e["file"], e["line"], "->", e.get("checked_start_line"))
print("PROBLEMS:", problems if problems else "none")
