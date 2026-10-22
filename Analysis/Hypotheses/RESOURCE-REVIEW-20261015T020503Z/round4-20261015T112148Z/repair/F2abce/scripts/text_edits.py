"""F2abce repair text edits. Writes out/text_edits_repair.json and checks that
every old_text occurs exactly once in the current file and that no two edits in
this workstream overlap."""
import json

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OUT = BASE + "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/out/"
LEG = "docs/manuscript/05_figure_legends.md"
RES = "docs/RESULTS.md"
STA = "docs/STATUS.md"
CAND = "figures/candidates/review-r3-fixes-20261015T112148Z/F2abce"

E = {}
E["R02-04"] = [
    dict(file=LEG, line="43-45", claim_level=False,
         old_text="**(B)** Fine-mapping PIP versus multi-signal COLOC posterior for genes with\n"
                  "finite high-resolution results; absence from this panel is not a negative\n"
                  "colocalization call.",
         new_text="**(B)** For each gene, the largest multi-signal COLOC posterior (PP.H4) across\n"
                  "the 35 GWAS strata against the SuSiE fine-mapping PIP of that same lead\n"
                  "variant in the same GWAS; when the lead variant has more than one\n"
                  "fine-mapping record in that GWAS, the largest SuSiE PIP is used. Of 3,918\n"
                  "genes with a multi-signal posterior, 866 have a SuSiE PIP at that lead variant\n"
                  "and 327 are shown (PP.H4 > 0.05 or PIP ≥ 0.99). The three MVP European strata\n"
                  "were not fine-mapped, so none of the 1,348 genes whose largest posterior comes\n"
                  "from them is shown. Eighteen genes have PIP ≥ 0.9 with PP.H4 < 0.5 and 18 have\n"
                  "PIP ≥ 0.9 with PP.H4 > 0.5. Color gives the hg19 consequence of the plotted lead\n"
                  "variant: noncoding (311), coding sequence or splice site of the plotted gene\n"
                  "(5), or of another gene only (11); coding includes synonymous positions.\n"
                  "Absence from this panel is not a negative colocalization call."),
    dict(file=RES, line="2654-2655", claim_level=False,
         old_text="single-signal only could not be tested by multi-signal COLOC and are not\n"
                  "negatives. Panel C consequence counts:",
         new_text="single-signal only could not be tested by multi-signal COLOC and are not\n"
                  "negatives. Panel B (round-3 review candidate in\n"
                  f"`{CAND}/repair/`, not adopted): x is the\n"
                  "SuSiE PIP of the same GWAS and lead variant that gives each gene's maximum\n"
                  "PP.H4, taking the largest SuSiE PIP when that variant has several fine-mapping\n"
                  "records. Of 3,918 genes, 883 have a fine-mapping record at that variant and 866\n"
                  "a SuSiE PIP; 327 are shown, 18 with PIP ≥ 0.9 and PP.H4 < 0.5 and 18 with\n"
                  "PIP ≥ 0.9 and PP.H4 > 0.5. The three MVP European strata have no fine-mapping,\n"
                  "so the 1,348 genes whose maximum PP.H4 comes from them have no x. For 58 of the\n"
                  "327 the records disagree; taking the smallest SuSiE PIP gives 326 shown and 17\n"
                  "(POLL falls below 0.99; RORA x 0.433 instead of 0.440). The installed panel\n"
                  "shows 695 points, 152 and 56, with x taken as the maximum over four methods, all\n"
                  "strata and all leads. Panel C consequence counts:"),
]
E["R02-05"] = []  # legend sentence merged into the R02-04 legend edit (conflict 31)
E["R02-06"] = [
    dict(file=LEG, line="71-72", claim_level=False,
         old_text="the lead variant, and 9 non-European-only genes passing lead-variant\n"
                  "genome-wide significance (8 EAS and 1 AFR).",
         new_text="the lead variant, and 9 non-European-only genes passing lead-variant\n"
                  "genome-wide significance (8 EAS and 1 AFR). The partition is per gene, not per\n"
                  "locus: 19 of the 32 non-European-only genes have their lead variant within\n"
                  "1 Mb of the lead variant of a gene with European multi-signal PP.H4 > 0.5. The\n"
                  "9 genome-wide-significant genes fall in 8 loci (COQ6 and ENTPD5 share one); 5\n"
                  "of them, in 4 loci (COQ6/ENTPD5, ELAPOR2, LUZP1, RGS12), have no such European\n"
                  "neighbor."),
    dict(file=RES, line="2657-2658", claim_level=False,
         old_text="23 non-European-only below the genome-wide-significance threshold, and 9\n"
                  "non-European-only passing it (8 EAS, 1 AFR).",
         new_text="23 non-European-only below the genome-wide-significance threshold, and 9\n"
                  "non-European-only passing it (8 EAS, 1 AFR). These are gene-level categories:\n"
                  "19 of the 32 non-European-only genes have their lead variant within 1 Mb of the\n"
                  "lead variant of a European multi-signal-positive gene (the same 19 with each\n"
                  "gene's minimum-p lead or maximum-PP.H4 lead), and the 9 passing genes form 8\n"
                  "loci, 4 of which (5 genes) have no European multi-signal-positive neighbour\n"
                  "within 1 Mb."),
]
E["R02-07"] = [
    dict(file=LEG, line="47-48", claim_level=False,
         old_text="the strongest consequence annotation of one lead COLOC variant per gene: 35\n"
                  "coding, 70 UTR, 104 promoter, 450 intronic, and 337 intergenic.",
         new_text="the strongest consequence annotation of one lead COLOC variant per gene: 35\n"
                  "coding, 70 UTR, 104 promoter, 450 intronic, and 337 intergenic. Of the 35\n"
                  "coding genes (29 distinct lead variants), 14 have the lead variant in their own\n"
                  "coding sequence or splice site and 21 have it only in another gene's, for\n"
                  "example C2orf16 at the GCKR variant, SAMM50 at the PNPLA3 variant and MAU2 at\n"
                  "the TM6SF2 variant. Coding includes synonymous positions."),
    dict(file=RES, line="2656", claim_level=False,
         old_text="450, intergenic 337 (996).",
         new_text="450, intergenic 337 (996); of the 35 coding genes, 14 have the lead variant in\n"
                  "their own coding sequence or splice site and 21 have it only in another gene's\n"
                  "(29 distinct lead variants; round-3 review candidate panel\n"
                  f"`{CAND}/Fig2C_coding_noncoding_split_ownership.pdf`)."),
]
E["R16-08"] = [
    dict(file=STA, line="859-861", claim_level=False,
         old_text="⚠ **Task `19648525_176` hit TIMEOUT after 3d18h on 2026-08-14 and is not covered\n"
                  "by that audit.** Its disposition is unrecorded, so the ROADMAP box \"all 1,100\n"
                  "tasks have a terminal, audited disposition\" cannot yet be ticked.",
         new_text="Task `19648525_176` (2023_36280732_NAFLD_UKBB_EUR chr1) hit TIMEOUT after 3d18h\n"
                  "on 2026-08-14 and was re-executed as `20633389_176` (COMPLETED 2026-08-16, 3 min).\n"
                  "The canonical chr1 file is byte-identical to the rerun output (md5\n"
                  "aa913ed3b92c65d0c2724b7dd7dfc0dd), and COLOC r2 carries 2,035 chr1 rows for this\n"
                  "study (6 with a SuSiE posterior, none above 0.5)."),
]

ok = True
spans = {}
for fid, eds in E.items():
    for t in eds:
        txt = open(BASE + t["file"]).read()
        c = txt.count(t["old_text"])
        i = txt.find(t["old_text"])
        ln = txt[:i].count("\n") + 1 if i >= 0 else None
        print(fid, t["file"], "count", c, "line", ln)
        ok &= c == 1
        for bad in ["—", "[Inference]", "[Unverified]", ">=", "olour"]:
            if bad in t["new_text"]:
                print("  BAD TOKEN", repr(bad)); ok = False
        spans.setdefault(t["file"], []).append((i, i + len(t["old_text"]), fid))
for f, s in spans.items():
    s.sort()
    for a, b in zip(s, s[1:]):
        if a[1] > b[0]:
            print("OVERLAP", f, a, b); ok = False
print("ALL OK" if ok else "PROBLEMS")
json.dump(E, open(OUT + "text_edits_repair.json", "w"), indent=1, ensure_ascii=False)
