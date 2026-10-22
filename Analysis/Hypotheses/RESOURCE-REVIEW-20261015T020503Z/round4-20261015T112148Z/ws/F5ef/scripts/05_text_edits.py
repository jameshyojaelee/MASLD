#!/usr/bin/env python3
"""F5ef: assemble proposed text edits as data. Reads documents only; writes out/text_edits.json.

old_text is taken verbatim from the current files (by substring), and each must occur exactly once.
Nothing is installed; the JSON is the return payload.
"""
from __future__ import annotations

import json
from pathlib import Path

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
WS = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef"
LEG = "docs/manuscript/05_figure_legends.md"
FIG5 = "docs/manuscript/draft/fig5.md"
ABS = "docs/manuscript/draft/abstract.md"
PAPER = "docs/PAPER.md"
RES = "docs/RESULTS.md"
METH = "docs/manuscript/METHODS.md"
STAT = "docs/STATUS.md"
CAND = "figures/candidates/review-r3-fixes-20261015T112148Z/F5ef"


def span(path: str, start: str, end: str) -> str:
    text = (P / path).read_text()
    i = text.index(start)
    assert text.count(start) == 1, (path, start)
    j = text.index(end, i) + len(end)
    return text[i:j]


def line_of(path: str, old: str) -> str:
    text = (P / path).read_text()
    i = text.index(old)
    a = text[:i].count("\n") + 1
    return f"{a}-{a + old.count(chr(10))}"


E = {}

# ---------------------------------------------------------------- R01-03 (Fig 5E roster)
E["R01-03"] = [
    (LEG, span(LEG, "**(E)** Assay-native program", "testability states."),
     "**(E)** Genetic-chromatin and protein-histology context for the complete\n"
     "registry of 117 frozen programs; rows are not selected by disease association.\n"
     "Rows show the 45 programs with at least one member gene in an adopted COLOC r2\n"
     "primary pair (the 767 pairs of panel D) that places at least half of its\n"
     "colocalization posterior in open chromatin of the program's own lineage in at\n"
     "least one snATAC source (28 programs; in 3 the gene is RORA from a direct\n"
     "MASLD trait), or with a robust protein-histology association (26 programs).\n"
     "All 117 rows are in the source table; display filtering triggers no\n"
     "re-testing. Bars count genes once each; dark shading marks a gene with a direct\n"
     "MASLD/liver-fat pair, and the other genes carry liver-enzyme pairs only.\n"
     "Intervals are rank-residual partial Spearman correlations between\n"
     "the PXD051911 program protein score and a histologic-burden axis (first\n"
     "principal component of steatosis, ballooning, inflammation and fibrosis) in 46\n"
     "MASLD participants, adjusted for acquisition batch, age, BMI and sex, with\n"
     "participant-bootstrap 95% CIs. Freedman-Lane permutation P (9,999) is\n"
     "BH-corrected over the 94 programs with at least eight measured member proteins\n"
     "and 20% retained weight; robust means q<0.05 with the equal-weight score, the\n"
     "score without its top gene and both acquisition batches agreeing in sign.\n"
     "Stromal ECM and Ductular injury are absent: neither contains a colocalized gene\n"
     "and neither can be scored in protein. Spatial organization is in panel F.",
     False),
    (FIG5, "gene, and absent overlap is unresolved context rather than a negative test.",
     "gene, and absent overlap is unresolved context rather than a negative test.\n\n"
     "Figure 5E places all 117 frozen programs in two assay-native contexts without\n"
     "selecting them by disease association. Twenty-eight programs contain at least\n"
     "one gene from an adopted colocalized pair whose posterior lies mainly in open\n"
     "chromatin of the program's own lineage; in three of them the gene is RORA from\n"
     "a direct MASLD trait, and the others carry liver-enzyme pairs. In liver protein\n"
     "from 46 MASLD participants, 26 of 94 scorable programs track a histologic-burden\n"
     "axis robustly (BH over 94): 16 negatively, including oxidative-phosphorylation,\n"
     "xenobiotic and amino-acid metabolism programs, and 10 positively, including\n"
     "fibroblast ECM and TGFB/collagen programs. The same protein cohort did not\n"
     "contribute to program discovery. Neither spatial program in Figure 5F contains a\n"
     "colocalized gene or can be scored in protein.",
     False),
    (METH, "concordance in this same cohort, so it is descriptive and is not an independent\n"
           "validation of the bulk or single-cell programs.",
     "concordance in this same cohort, so it is descriptive and is not an independent\n"
     "validation of the bulk or single-cell programs.\n\n"
     "For Figure 5E, single-gene protein groups with at most 50% missing values are\n"
     "log2-transformed, quantile-normalized, collapsed to one value per symbol by the\n"
     "median and z-scored per protein across the 58 samples. Each of the 117 frozen\n"
     "programs is scored per participant as the L1-weighted mean of its measured\n"
     "member proteins; a program is scorable with at least eight measured members and\n"
     "at least 20% of its L1 weight retained (94 of 117). Histologic burden is the\n"
     "first principal component of steatosis, ballooning, inflammation and fibrosis\n"
     "grades in the 46 MASLD participants with complete histology; NAS is excluded\n"
     "because it sums three of these. The association is a rank-residual partial\n"
     "Spearman correlation adjusted for acquisition batch, age, BMI and sex, with a\n"
     "batch-stratified Freedman-Lane permutation P (9,999 permutations), a\n"
     "participant bootstrap 95% CI resampled within batch (4,999 draws) and\n"
     "Benjamini-Hochberg correction over the 94 scorable programs. A program is\n"
     "robust when q<0.05 and the equal-weight score, the score without its\n"
     "top-weighted gene and both acquisition batches agree in sign. Seeds are\n"
     "20260923 plus the program index (permutation) and 30260923 plus the index\n"
     "(bootstrap).",
     False),
    (STAT, "5A/5E/5F retain assay-native context and all-cell-adjusted spatial calibration.",
     "5A retains assay-native context. Candidate 5E (not installed; "
     f"`{CAND}/`) replaces the 13 rows chosen from 22 programs selected on the superseded July "
     "disease fit with the complete 117-program registry: 45 rows with a lineage-open COLOC r2 gene "
     "(28 programs, 3 with direct-trait RORA) or a robust protein-histology association (26 of 94 "
     "scorable, BH over 94), and no spatial column. 5F keeps the all-cell-adjusted calibration; its "
     "candidate source table moves Vu Ductular injury from `source_dependent` to `indeterminate`.",
     False),
    (STAT, span(STAT, "- Retain the adopted program-level Figure 5E/5F panels; positive Moran",
                "direction alone is not the main spatial result."),
     "- Review the 117-program Figure 5E candidate and the corrected Figure 5F source\n"
     f"  table in `{CAND}/`; positive Moran\n"
     "  direction alone is not the main spatial result.",
     False),
]

# ---------------------------------------------------------------- R11-03 (credible sets vs r2)
E["R11-03"] = [
    (METH, "separate. Partial or absent peak overlap is unresolved context, not a negative\nfunctional test.",
     "separate. Partial or absent peak overlap is unresolved context, not a negative\n"
     "functional test. For Figure 5E, a program member gene counts as open when one of\n"
     "its primary pairs places at least half of its posterior in open chromatin of the\n"
     "program's own lineage in at least one source (states `replicated_accessible` or\n"
     "`source_dependent`); hepatocyte, cholangiocyte, macrophage, fibroblast and T-cell\n"
     "programs map to hepatocyte, cholangiocyte, macrophage, stellate and T/NK\n"
     "chromatin. The count uses the same 767 COLOC r2 pairs as Figure 5D, including the\n"
     "three MVP European-ancestry studies.",
     False),
]

# ---------------------------------------------------------------- R06-02 ("equally robust")
fig5_old = span(FIG5, "Figure 5F resolves tissue organization for the two prespecified robust",
                "dependence and the matched-null construction are in Methods.")
fig5_tail = (
    "After adjustment for all 16 cell2location abundances, ECM/IGFBP7\n"
    "exceeds a matched null built from genes of equivalent expression, detection,\n"
    "biotype and lineage correlation in both sources (excess 0.040 and 0.050),\n"
    "whereas ductular-injury/BICC1 does not (0.015 and 0.003) and is indeterminate\n"
    "rather than negative. Residual Moran's I is positive for ductular-injury/BICC1\n"
    "in 4 of 4 donors and 9 of 10 arrays, so the raw statistic does not separate the\n"
    "programs. The matched-gene null holds donors fixed. At the donor level,\n"
    "ECM/IGFBP7 excess is positive in 3 of 4 donors, and the between-program\n"
    "difference is positive in 3 of 4 donors (exact sign-flip P=0.19); with four\n"
    "donors no sign-flip test can reach P<0.05. Disease association, cell-of-origin\n"
    "assignment and tissue architecture are separate claims requiring separate\n"
    "calibration. Donor joins, array-level source dependence and the matched-null\n"
    "construction are in Methods.")
E["R06-02"] = [
    (FIG5, fig5_old,
     "Figure 5F examines tissue organization for two prespecified programs; neither\n"
     "was selected by a stage-association rule. " + fig5_tail, False, "H59 null; H90 fail or untestable"),
    (FIG5, fig5_old,
     "Figure 5F examines tissue organization for two prespecified programs. Neither\n"
     "passed the corrected 117-program stage rule; Ductular injury is a prespecified\n"
     "program whose stage association replicated in one held-out cohort. " + fig5_tail,
     False, "H59 null; H90 pass (apply instead of the row above)"),
    (RES, "lncRNA-independent. Both source-diagnosis-refit primary-selected programs passed this",
     "lncRNA-independent. Both prespecified programs (Stromal ECM and Ductular injury) passed this",
     False),
]

# ---------------------------------------------------------------- G3-03 (lineage-only vs all-16 in PAPER)
paper_old = span(PAPER, "For spatial omics, the biological finding is that robust MASLD programs differ",
                 "interpretable.**")
paper_tail = (
    "**Reproducibly positive spatial autocorrelation is not sufficient evidence\n"
    "that a disease program has program-specific tissue architecture.** After\n"
    "residualizing for all 16 cell2location abundances, both programs had positive\n"
    "aggregate residual Moran's I in both Visium sources, and residual Moran's I was\n"
    "positive in 3/4 and 4/4 GSE192741 donors and 10/10 and 9/10 Vu arrays for\n"
    "Stromal ECM and Ductular injury. Matched expression, detection, biotype and\n"
    "lineage-correlated gene nulls separated them at the program level: Stromal ECM\n"
    "exceeded expected spatial background in both sources (excess 0.040, q=0.001;\n"
    "0.050, q=0.0006), whereas Ductular injury did not (0.015, q=0.079; 0.003,\n"
    "q=0.32). This null resamples genes with donors held fixed. At the donor level,\n"
    "Stromal ECM excess was positive in 3 of 4 donors (exact one-sided sign-flip\n"
    "P=0.19), as was the Stromal ECM minus Ductular injury difference (P=0.19); four\n"
    "donors cannot reach P<0.05 (smallest attainable P 0.0625). The Vu result remains\n"
    "source-dependent because arrays cannot be joined to donors. Calibration against\n"
    "an assay's intrinsic spatial propensity, not the sign or apparent replication of\n"
    "Moran's I, is therefore required before spatial organization is called\n"
    "program-specific.")
E["G3-03"] = [
    (PAPER, paper_old,
     "For spatial omics, the finding is that two prespecified MASLD programs differ in\n"
     "whether their spatial pattern exceeds matched-gene background. " + paper_tail,
     True, "H59 null; H90 fail or untestable"),
    (PAPER, paper_old,
     "For spatial omics, the finding is that two prespecified MASLD programs, one of\n"
     "which (Ductular injury) replicated its stage association in a held-out cohort,\n"
     "differ in whether their spatial pattern exceeds matched-gene background. " + paper_tail,
     True, "H59 null; H90 pass (apply instead of the row above)"),
    (RES, "sections had residual Moran's I of approximately 0.080 and 0.328, but they",
     "sections had residual Moran's I of 0.080 and 0.328 under lineage-only adjustment\n"
     "(0.003 and 0.077 under the all-cell adjustment of Figure 5F), but they",
     False),
]

# ---------------------------------------------------------------- G3-04 ("distinguished")
E["G3-04"] = [
    (ABS, span(ABS, "physical tissue contexts. Matched spatial nulls distinguished an ECM/IGFBP7",
               "and cellular origin rather than treating them as interchangeable evidence."),
     "physical tissue contexts. In two Visium sources, an ECM/IGFBP7 program\n"
     "exceeded a matched-gene spatial null after adjustment for cell composition,\n"
     "whereas a ductular-injury/BICC1 program did not; with four donors, neither the\n"
     "ECM/IGFBP7 result nor the difference between the programs reaches P<0.05 at the\n"
     "donor level. Disease association, tissue architecture and cellular origin are\n"
     "therefore reported as separate evidence.",
     True),
    (PAPER, span(PAPER, "Matched spatial nulls independently\n   support program-specific organization for Stromal ECM, whereas",
                 "disease remodeling, tissue architecture, and lineage attribution."),
     "After all-cell adjustment, Stromal ECM\n"
     "   exceeds a matched-gene spatial null in both Visium sources, whereas Ductular\n"
     "   injury remains indeterminate; the null holds donors fixed, Stromal ECM excess\n"
     "   is positive in 3 of 4 donors, and no donor-level test separates the two\n"
     "   programs. Disease remodeling, tissue architecture, and lineage attribution are\n"
     "   therefore reported as separate claims.",
     True),
    (PAPER, "| Stromal ECM exceeds matched spatial background while Ductular injury remains "
            "indeterminate, showing that disease association, tissue architecture, and cellular "
            "origin are separate claims |",
     "| Stromal ECM exceeds a matched-gene spatial null after all-cell adjustment while Ductular "
     "injury remains indeterminate; with four donors no donor-level test separates them, so disease "
     "association, tissue architecture, and cellular origin are reported as separate claims |",
     True),
]

# ---------------------------------------------------------------- R12-07 (P floor) + all-16 primary in RESULTS
res_old = span(RES, "The confirmatory family was recorded before execution and contains only the Stromal",
               "source-dependent rather than donor-level replication.")
E["R12-07"] = [
    (RES, res_old,
     "The confirmatory family was recorded before execution and contains only the Stromal\n"
     "ECM and Ductular injury programs. The values below are those plotted in Figure 5F:\n"
     "spot scores residualized for all 16 cell2location q05 abundances plus log library\n"
     "size and log detected genes, expressed as excess over the matched-gene null mean,\n"
     "with 9,999 matched-gene draws (seed 42 for GSE192741, 1042 for Vu) and BH over the\n"
     "two programs within each source. In GSE192741, Stromal ECM had donor-first excess\n"
     "0.04050 (matched-null SD 0.01098; 4 of 9,999 draws at or above it, P=0.0005,\n"
     "q=0.001) and is `supported`. Ductular injury had excess 0.01532 (SD 0.01060,\n"
     "P=q=0.0793) and is `indeterminate`, not negative.\n\n"
     "In Vu, Stromal ECM had excess 0.05034 (SD 0.01033; 2 of 9,999 draws, P=0.0003,\n"
     "q=0.0006). Its evidence state is `source_dependent` because the ten arrays cannot\n"
     "be joined to biological donors. Ductular injury had excess 0.002654 (SD 0.008227,\n"
     "P=q=0.3158) and is `indeterminate`.\n\n"
     "With the program's own lineage abundance as the only composition covariate, the\n"
     "excesses were larger: Stromal ECM 0.1390 in GSE192741 and 0.1471 in Vu, each with\n"
     "no null draw at or above it (P ≤ 1e-4, the resolution floor of 9,999 draws;\n"
     "q ≤ 0.0002), and Ductular injury 0.02587 (P=q=0.0755) and 0.01922 (P=q=0.1158).\n"
     "This sensitivity changes no evidence state.\n\n"
     "The matched-gene null resamples genes; donors and arrays are held fixed. Under the\n"
     "plotted all-cell model, residual Moran's I was positive for Stromal ECM in 3/4\n"
     "GSE192741 donors and 10/10 Vu arrays, and for Ductular injury in 4/4 donors and\n"
     "9/10 arrays, so neither positive sign nor cross-unit consistency distinguishes the\n"
     "programs. Stromal ECM excess over each donor's own matched-gene null was positive\n"
     "in 3 of 4 donors (0.0156, 0.0657, 0.1023 and -0.0216 for H35 to H38; exact\n"
     "one-sided sign-flip P=0.19; the floor with four donors is 0.0625). The per-donor\n"
     "Stromal ECM minus Ductular injury excess was positive in 3 of 4 donors (P=0.19)\n"
     "and in 10 of 10 Vu arrays (P=0.00098; arrays are not donors). The program-level\n"
     "matched-gene calibration therefore separates the programs, but no donor-level\n"
     "test does. **Generic spatial autocorrelation can mimic program transportability,\n"
     "and an assay-matched null is required to separate the two.** This is not a\n"
     "negative Ductular injury result, and the Vu contrast remains source-dependent\n"
     "rather than donor-level replication.",
     False),
]

# ---------------------------------------------------------------- R12-03 + R06-11 (5F legend)
E["R12-03"] = [
    (LEG, span(LEG, "**(F)** At\nleft, a registered H&E overlay", "whereas Ductular injury is indeterminate in both."),
     "**(F)** At\n"
     "left, a registered H&E overlay shows the all-cell-adjusted Stromal ECM score in\n"
     "one Guilliams et al. steatotic section (JBO019, donor H37) selected outcome-blind\n"
     "as closest to the dataset median eligible spot count; its own excess, 0.102, is\n"
     "the largest of the five sections (range -0.022 to 0.102). At right, points show\n"
     "Moran excess over matched genes for the two prespecified programs in Guilliams\n"
     "et al. and Vu et al., after residualizing spot scores for all 16 cell2location\n"
     "abundances, log library size and log detected genes; labels are BH-adjusted\n"
     "empirical q-values from 9,999 matched-gene draws, BH over the two programs within\n"
     "each source. The null resamples genes with donors and arrays held fixed. Stromal\n"
     "ECM exceeds its null in both sources (supported in Guilliams et al.;\n"
     "source-dependent in Vu et al. because its arrays lack a public donor key), whereas\n"
     "Ductular injury is indeterminate in both. Stromal ECM excess is positive in 3 of 4\n"
     "Guilliams et al. donors (exact one-sided sign-flip P=0.19; the smallest attainable\n"
     "P with four donors is 0.0625) and in 10 of 10 Vu et al. arrays, which are not\n"
     "donors; no donor-level test separates the two programs.",
     False),
]

# ---------------------------------------------------------------- G3-07 (overlay section value)
E["G3-07"] = [
    (RES, "library size and detected genes. Re-derivation from the spot matrix reproduced\n"
          "the recorded residual Moran's I exactly (0.0582734890).",
     "library size and detected genes. Re-derivation from the spot matrix reproduced\n"
     "the recorded dataset-level residual Moran's I exactly (0.0582734890); this is the\n"
     "four-donor collapsed value, not the overlay section's. Section JBO019 itself has\n"
     "residual Moran's I 0.1213 and excess 0.1023 over its matched-gene null, the\n"
     "largest of the five GSE192741\n"
     "sections (excess -0.0151, 0.0462, 0.0657, 0.1023 and -0.0216 for JBO014, JBO015,\n"
     "JBO018, JBO019 and JBO022): it was chosen by spot count, not by score, but it\n"
     "shows the strongest section.",
     False),
]

# ---------------------------------------------------------------- G3-08 (Methods)
E["G3-08"] = [
    (METH, span(METH, "Residual Moran's I for a program score is compared against a null set of genes",
                "donor-level inference."),
     "Each Visium spot receives a program score from its measured fixed member genes,\n"
     "weighted by their original L1 weights. The score is residualized for all 16\n"
     "cell2location q05 abundance factors, log1p library size and log1p detected genes\n"
     "(Figure 5F); a sensitivity analysis replaces the 16 abundances with the program's\n"
     "own lineage abundance. Residual Moran's I on a six-neighbour spot graph is\n"
     "compared with 9,999 random sets of genes matched on expression, detection rate,\n"
     "biotype and lineage correlation, scored and residualized identically (seed 42 for\n"
     "GSE192741, 1042 for Vu). The reported excess is the observed statistic minus the\n"
     "null mean; P=(k+1)/(9,999+1), where k counts draws at or above the observed value,\n"
     "so the smallest attainable P is 1e-4. Benjamini-Hochberg correction covers the\n"
     "two prespecified programs within each source. The null resamples genes and holds\n"
     "donors and arrays fixed. Donor-first estimates are used where a donor key exists\n"
     "(five GSE192741 sections collapse to four donors); an exact sign-flip test of\n"
     "per-donor excess has a smallest attainable P of 0.0625 with four donors and is\n"
     "reported descriptively. Array-level results without a donor key are reported as\n"
     "`source_dependent` and never as donor-level inference.",
     False),
]

# ---------------------------------------------------------------- R06-10 (V4 adjustment)
E["R06-10"] = [
    (RES, span(RES, "This result is spatial consistency under proxy abundance adjustment; it does\n"
                    "not establish a disease/control or stage contrast.",
               "not establish a disease/control or stage contrast."),
     "This result is spatial consistency under proxy abundance adjustment; it does\n"
     "not establish a disease/control or stage contrast. Its adjustment is weaker than\n"
     "Figure 5F's: HMSMA scores are residualized for a 15-marker hepatocyte proxy, and\n"
     "the GSE192741 and Vu rows for the cell2location hepatocyte abundance only, not all\n"
     "16 abundances. Its 0.0983 therefore cannot be compared with the 5F excess (0.040\n"
     "and 0.050) and does not replicate it. The Vu units are arrays holding several\n"
     "biopsies each, not donors.",
     False),
    (STAT, "GSE192741, 4 donors: 0.139, exact p 0.0625, the smallest possible with 4 donors, outside the family.",
     "GSE192741, 4 donors: 0.139, exact p 0.0625, the smallest possible with 4 donors, outside the family. "
     "These values adjust for the lineage abundance (GSE192741, Vu) or a 15-marker hepatocyte proxy "
     "(HMSMA), not the 16 abundances of Figure 5F (Stromal ECM excess 0.040 and 0.050), so V4 does not "
     "replicate the 5F statistic; Vu arrays are not donors.",
     False),
]

out = {}
for fid, edits in E.items():
    rows = []
    for e in edits:
        path, old, new, claim = e[:4]
        cond = e[4] if len(e) > 4 else ""
        text = (P / path).read_text()
        assert text.count(old) == 1, (fid, path, old[:80], text.count(old))
        assert "\u2014" not in new and "\u2013" not in new, (fid, "dash")
        ln = line_of(path, old)
        rows.append({"file": path, "line": ln + (f" [{cond}]" if cond else ""), "old_text": old,
                     "new_text": new, "claim_level": claim})
    out[fid] = rows
(WS / "out" / "text_edits.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
for fid, rows in out.items():
    for r in rows:
        print(fid, r["file"], r["line"], len(r["old_text"]), "->", len(r["new_text"]), "claim" if r["claim_level"] else "")
