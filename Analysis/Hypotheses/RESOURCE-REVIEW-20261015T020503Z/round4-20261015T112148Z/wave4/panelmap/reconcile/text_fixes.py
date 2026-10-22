"""Proposed staged-text fixes after panel install. Checks each old_text occurs exactly once.
Writes text_fixes.json next to this script. Edits nothing."""
import json, pathlib

S = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/"
                 "RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/staged-20261015T184532Z/final")
OUT = pathlib.Path(__file__).with_name("text_fixes.json")

F = []
def add(file, old, new, reason):
    F.append(dict(file=file, old_text=old, new_text=new, reason=reason))

# ---- Fig 1B (gated on H57)
add("docs/STATUS.md",
    "the installed Figure 1B still prints the earlier universe (26 of 1,155 DEGs, 2.3%) until the repaired panel is installed",
    "the installed Figure 1B prints this universe (26 of 220 DEGs, 11.8%, against 9.6% expected within expression deciles; 935 further canonical DEGs had no multi-signal COLOC result)",
    "Fig 1B after install prints 26 observed (11.8%) of 220 DEGs, expected 9.6%, and the 935-DEG note; the sentence says it prints 26 of 1,155 (2.3%). Only if H57 is approved and the panel installed.")

# ---- Fig 2B
add("docs/RESULTS.md",
    "Panel B (round-3 review candidate in\n`figures/candidates/review-r3-fixes-20261015T112148Z/F2abce/repair/`, not adopted): x is the",
    "Panel B (installed 2026-10-15 from the round-3 review repair): x is the",
    "After install this 327-gene panel is the main Fig 2B; the text calls it a non-adopted candidate.")
add("docs/RESULTS.md",
    "The installed panel\nshows 695 points, 152 and 56, with x taken as the maximum over four methods, all\nstrata and all leads.",
    "The previous panel\nshowed 695 points, 152 and 56, with x taken as the maximum over four methods, all\nstrata and all leads.",
    "Installed Fig 2B will print 327 genes (corner 18, top-right 18); 695/152/56 describe the replaced panel.")

# ---- Fig 2C
add("docs/STATUS.md",
    "Installed Figure 2C prints 35/70/104/450/337",
    "Installed Figure 2C prints 14/21/70/104/450/337 (coding split into own gene 14 and other gene 21)",
    "Installed Fig 2C will print 'coding, this gene (14)' and 'coding, other gene (21)', not one 35 segment.")

# ---- Fig 2 noncoding companion
add("docs/STATUS.md",
    "The installed companion still prints unresolved fractions of 52.11% and\n60.43%, computed by summing context posterior over all credible sets of a\nclass; H58 (2026-10-15) withdrew them (failed lookup; repaired value 0%). A\ncorrected candidate companion is\n`figures/candidates/review-r3-fixes-20261015T112148Z/F2d/repair/Fig2_companion_noncoding_context_repair.pdf`.",
    "The installed companion, summed over all credible sets of a class, prints\npromoter 13.84% / 14.67%, liver-lineage open chromatin 9.26% / 8.44%, liver ABC\nenhancer 0.43% / 3.70% and no deposited context 77.90% / 77.69% (direct / enzyme).\nThe earlier unresolved fractions of 52.11% and 60.43% reflected failed lookup\nand were withdrawn by H58 (2026-10-15; repaired value 0%).",
    "After replacement the companion prints no unresolved row and the repaired context values; the text says it still prints 52.11%/60.43%.")
add("docs/RESULTS.md",
    "class: unresolved 52.11% / 60.43%, promoter 7.87% / 4.46%, accessible\n7.32% / 4.15%, ABC 0.23% / 1.57% (direct / enzyme). H58 (2026-10-15) withdrew\nthe unresolved values (failed lookup). Repaired values, same estimator\n(candidate panel\n`figures/candidates/review-r3-fixes-20261015T112148Z/F2d/repair/Fig2_companion_noncoding_context_repair.pdf`):\nunresolved 0% / 0%,",
    "class (direct / enzyme). The previous companion printed unresolved 52.11% / 60.43%,\npromoter 7.87% / 4.46%, accessible 7.32% / 4.15% and ABC 0.23% / 1.57%; H58\n(2026-10-15) withdrew these (failed lookup). Repaired values, same estimator,\nas printed on the installed panel:\nunresolved 0% / 0%,",
    "Text says the installed companion displays 52.11%/60.43%, 7.87%/4.46%, 7.32%/4.15%, 0.23%/1.57%; after replacement it prints 13.84%/14.67%, 9.26%/8.44%, 0.43%/3.70%, 77.90%/77.69%.")
add("figures/main/fig2_genetics/README.md",
    "Its unresolved\nfractions, 52.11% and 60.43%, reflected failed context lookup and are withdrawn\n(H58 repaired value 0% in both classes).",
    "It prints promoter\n13.84% / 14.67%, liver-lineage open chromatin 9.26% / 8.44%, liver ABC enhancer\n0.43% / 3.70% and no deposited context 77.90% / 77.69% (direct / enzyme). The\nearlier unresolved fractions, 52.11% and 60.43%, reflected failed context lookup\nand were withdrawn (H58 repaired value 0% in both classes).",
    "README says the file at this path has unresolved fractions 52.11%/60.43%; the replaced file prints no unresolved row and gives the repaired values.")

# ---- Fig 2A support companion (window unit)
add("docs/manuscript/05_figure_legends.md",
    "(Figure 2A support companion, not assembled in canonical A-I)",
    "(Figure 2A support companion, not assembled in canonical A-I; its Windows column counts fine-mapped loci per study, summed over the studies in each trait-ancestry cell, 1,665 in all, so it does not equal the 360 trait-specific loci above)",
    "Companion prints a Windows total of 1,665 (e.g. ALT 516 over five cells) beside Fig 2A's 360 loci (ALT 122) with no stated unit; UKBB EUR rows 6/9/15 equal the per-study loci named in this legend.")

# ---- Fig 4H outline (gated on user reversal of 2026-10-13 choice)
add("docs/RESULTS.md",
    "The current 774-unit main panel has one uniform-null\n  outline; the stage unit is score SD per cohort stage-code SD.",
    "The current main panel uses the 367 disease units (rho 0.63, 32/43 same\n  direction) and has no outline; the stage unit is score SD per cohort stage-code SD.",
    "Installed Fig 4H will print rho = 0.63, 32/43 and no outline or 'Exceeds size-matched null' key.")

# ---- Fig 4E arm (gated on user arm choice)
add("docs/RESULTS.md",
    "The installed Figure 4E keeps the four-cohort arm until the user chooses between the two arms.",
    "The installed Figure 4E uses the three-cohort 364-unit arm: no fibrosis or NAS program reaches BH q<0.05 (smallest q 0.079 and 0.159); the four-cohort arm is a sensitivity.",
    "Installed Fig 4E will show 0 points above the q=0.05 line; the text says it shows the four-cohort arm with 16 fibrosis and 1 NAS call.")

# ---- Fig 4 draft note
add("docs/manuscript/draft/fig4.md",
    "> Figure 4A and Figure 4F were installed on 2026-10-14. The Figure 4C, 4E, 4G\n> and 4H text below describes candidate panels that are not yet installed.",
    "> Figure 4A and Figure 4F were installed on 2026-10-14; Figure 4C, 4E, 4G and\n> 4H were installed on 2026-10-15 and match the text below.",
    "After install the note says the 4C/4E/4G/4H panels described are not installed. Apply only for the panels actually installed (4C needs H59, 4E the arm choice, 4H the outline reversal).")

# ---- S4 native forest unit
add("docs/manuscript/05_figure_legends.md",
    "shows hepatocyte disease-grade effects in the corrected 62-donor",
    "shows hepatocyte native-score effects per stage step in the corrected 62-donor",
    "Replaced S4 native forest axis reads 'Native score per stage step (HC3 95% CI)'; the legend calls the unit disease grade.")

# ---- Fig 5D
add("docs/RESULTS.md",
    "The installed October 14 display separates\n28 direct MASLD/liver-fat pairs from 739 liver-enzyme pairs (767 total; 22 and 424 genes, ten genes in both classes; the review-round-3 candidate counts genes), read\nfrom the rendered panel.",
    "The installed display counts genes, each at its best state over its study pairs\nwith a 5% posterior floor: 22 direct MASLD/liver-fat genes and 424 liver-enzyme genes\n(from 28 and 739 pairs, 767 total; ten genes in both classes), read\nfrom the rendered panel.",
    "Installed Fig 5D will print n = 22 genes and n = 424 genes and 'Genes (%)', not 28 and 739 pairs.")
add("docs/RESULTS.md",
    "(review-round-3 candidate panel\n`figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/fig5d_candidate_locus_matched_background.pdf`,\nnot part of the installed Figure 5D)",
    "(Figure 5D locus-matched background companion)",
    "Legend (lines 523-529) describes this panel as part of Fig 5D; the text says it is not.")
add("docs/STATUS.md",
    "5D uses 767 r2 gene-study pairs (28 direct, 739 enzyme; 22 and 424 genes). A gene-level candidate with a 5% posterior floor is in figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/repair/.",
    "5D counts genes with a 5% posterior floor: 22 direct and 424 enzyme genes, from 767 r2 gene-study pairs (28 direct, 739 enzyme); a locus-matched background companion uses the same genes.",
    "Installed Fig 5D will print 22/424 genes; the text calls the gene-level panel a candidate.")

# ---- Fig 5E
add("docs/STATUS.md",
    "Candidate 5E (not installed; `figures/candidates/review-r3-fixes-20261015T112148Z/F5ef/`) replaces",
    "Installed 5E (2026-10-15) replaces",
    "Numbers match the panel (45, 28, 3, 26, 94); only the 'not installed' state is wrong after install.")

# ---- Fig 5 assembly
add("figures/main/fig5_molecular_context/README.md",
    "The current vector one-page assembly `fig5_molecular_context.pdf` was rebuilt\non 2026-10-14 from unchanged 5A/E/F, corrected 774-unit 5B/C and verified\nr2 5D (job 22236878).",
    "The vector one-page assembly `fig5_molecular_context.pdf` (job 22236878,\n2026-10-14) still prints the pair-level 5D (28 and 739 pairs) and the 13-row 5E;\nit must be rebuilt from the installed gene-level 5D and 45-row 5E.",
    "After install the one-page Fig 5 still prints 28/739 pair bars and 13 5E rows while the panels print 22/424 genes and 45 rows.")

# ---- Fig 6 assemblies
add("figures/main/fig6_gene_catalog/README.md",
    "is complete on\n2026-10-14 (468 × 970 pt; job 22237086), with all source panels unchanged.",
    "was built on\n2026-10-14 (468 × 970 pt; job 22237086) from the earlier B, C and D; it still\nprints the tested-negative zeros, '26 concordant', the unsplit 13,551 bar and no\nworked-entry estimates, and must be rebuilt from the installed panels.",
    "Review assembly prints a tested-negative row/column of zeros and '26 concordant'; installed B/C/D print a 3 x 3 grid, 'both supported', native estimates and the 4,232 + 9,319 split.")
add("docs/manuscript/05_figure_legends.md",
    "is a candidate assembly; it preserves these panel values and is not adopted.",
    "is a candidate assembly built from the October 14 panels; it labels the 26 genes 'concordant', prints no native estimates in (C) and no 4,232 + 9,319 split in (D), and is not adopted.",
    "Compact layout prints HKDC1 'Concordant' with no PP.H4 or log2FC and '26 concordant'; the legend says it preserves the panel values.")

# ---- Fig 7F legend
add("docs/manuscript/05_figure_legends.md",
    "Study-level deltas are reported both unweighted and\nweighted by donors.",
    "Per-study deltas are not drawn in this panel; the gains\ndrawn are all-study values.",
    "Panel 7F draws no per-study deltas or weighted means; the legend says it reports them.")
add("docs/manuscript/05_figure_legends.md",
    "it contributes -0.45 to -0.49 where\nevery other study contributes -0.065 to -0.002.",
    "it contributes -0.45 to -0.49 in three\nblocks, where the other six studies contribute -0.043 to +0.097 (43 of 72 deltas\npositive).",
    "fig7f_per_study_delta_ranges.tsv: other six studies span -0.042681 to +0.096671 with 43 of 72 positive; -0.065 to -0.002 is GSE189600's own range in the other nine blocks.")

# ---- Fig 7D units
add("docs/manuscript/draft/fig7.md",
    "Its three primary evaluations comprise 109 GSE268273,\n299 GSE213621 and 177 GSE276114 participant records, with source-declared\nidentity rather than independent genetic confirmation. Raw Spearman\ncorrelations are 0.6631, 0.6382 and 0.6739.",
    "In the inferential units drawn in Figure 7D, Spearman is 0.6631 in 109\nGSE268273 participant records (source-declared identity), 0.6680 in 263\nGSE213621 genetic identities (11 identities with conflicting stage labels\nremoved) and 0.7860 in 81 GSE276114 participants with MASLD; in 82 GSE276114\nparticipants with chronic viral hepatitis it is 0.5443. On the published 299\nGSE213621 and 177 GSE276114 participant records, raw Spearman correlations are\n0.6382 and 0.6739.",
    "Fig 7D draws GSE213621 as 263 identities (0.668) and GSE276114 as 81 MASLD (0.786) and 82 CVH (0.544); the text leads with 299/177 records (0.638/0.674).")
add("docs/RESULTS.md",
    "| Same instrument, two additional external cohorts and three contrasts, Spearman / tie ceiling | **0.677 / 0.774 / 0.794** | GSE213621 299, GSE213621 4-level 367, GSE276114 177 participant records;",
    "| Same instrument, two additional external cohorts, Spearman / tie ceiling in the Figure 7D units | **0.709 / 0.866 / 0.680** (record-level 0.677 / 0.774 / 0.794) | GSE213621 263 genetic identities, GSE276114 81 MASLD and 82 chronic viral hepatitis participants (Spearman 0.6680, 0.7860, 0.5443); record-level: GSE213621 299, GSE213621 4-level 367, GSE276114 177 participant records;",
    "Fig 7D ratios implied by the drawn rows are 0.709, 0.866 and 0.680; the table gives only the record-level 0.677/0.794.")
add("docs/manuscript/METHODS.md",
    "The subsequent GSE213621 disease-only and\ncontrol-inclusive contrasts use 299 and 367 records from the same cohort;\nGSE276114 contributes 177.",
    "Figure 7D evaluates GSE213621 in 263 genetic identities (one record\neach; 11 identities whose two libraries carry different stage labels removed)\nand GSE276114 as 81 MASLD and 82 chronic viral hepatitis participants. The\npublished record-level contrasts use 299 disease-only and 367 control-inclusive\nGSE213621 records and 177 GSE276114 records.",
    "Methods gives only 299/367/177 records; Fig 7D draws 263 identities and the 81/82 etiology split.")
add("docs/STATUS.md",
    "two additional external cohorts, three contrasts, Spearman/ceiling 0.677–0.794",
    "two additional external cohorts, Spearman/ceiling 0.709 (263 GSE213621 identities), 0.866 (81 GSE276114 MASLD) and 0.680 (82 GSE276114 chronic viral hepatitis) as drawn in Figure 7D (record-level 0.677–0.794)",
    "Fig 7D rows give 0.709, 0.866 and 0.680; STATUS gives only the record-level range.")

bad = []
for f in F:
    txt = (S / f["file"]).read_text()
    n = txt.count(f["old_text"])
    f["count"] = n
    idx = txt.find(f["old_text"])
    f["line"] = txt[:idx].count("\n") + 1 if idx >= 0 else -1
    if n != 1:
        bad.append((f["file"], n, f["old_text"][:80]))
OUT.write_text(json.dumps(F, indent=1, ensure_ascii=False))
print(len(F), "fixes;", "bad:", bad)
for f in F:
    print(f["count"], f["file"], f["line"])
