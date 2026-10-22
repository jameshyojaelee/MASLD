import json, re, sys
B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
edits = [
 ("R02-04","docs/manuscript/05_figure_legends.md","**(B)** Fine-mapping PIP versus multi-signal COLOC posterior for genes with\nfinite high-resolution results; absence from this panel is not a negative\ncolocalization call."),
 ("R02-04","docs/RESULTS.md","single-signal only could not be tested by multi-signal COLOC and are not\nnegatives. Panel C consequence counts:"),
 ("R02-06","docs/manuscript/05_figure_legends.md","the lead variant, and 9 non-European-only genes passing lead-variant\ngenome-wide significance (8 EAS and 1 AFR)."),
 ("R02-06","docs/RESULTS.md","23 non-European-only below the genome-wide-significance threshold, and 9\nnon-European-only passing it (8 EAS, 1 AFR)."),
 ("R02-07","docs/manuscript/05_figure_legends.md","the strongest consequence annotation of one lead COLOC variant per gene: 35\ncoding, 70 UTR, 104 promoter, 450 intronic, and 337 intergenic."),
 ("R02-07","docs/RESULTS.md","450, intergenic 337 (996)."),
 ("R16-08","docs/STATUS.md","⚠ **Task `19648525_176` hit TIMEOUT after 3d18h on 2026-08-14 and is not covered\nby that audit.** Its disposition is unrecorded, so the ROADMAP box \"all 1,100\ntasks have a terminal, audited disposition\" cannot yet be ticked."),
]
for fid,f,old in edits:
    t = open(B+f, encoding="utf-8").read()
    n = t.count(old)
    if n:
        i = t.index(old); line = t[:i].count("\n")+1
    else:
        line = None
    print(fid, f, "count", n, "start_line", line)
