#!/usr/bin/env python3
"""Check that each proposed old_text exists verbatim at the stated line(s); scan new_text for em dashes."""
import json, re, sys
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
E = [
 ("R03-03","docs/manuscript/05_figure_legends.md","113","All panels use the corrected 774-unit substrate. Per-cohort membership,"),
 ("R03-03","docs/manuscript/05_figure_legends.md","116-117","**(A–B)** Five-cohort coverage (326, 207, 140, 71 and 30 units) and documented\nfibrosis-by-NAS support."),
 ("R03-03","docs/manuscript/draft/fig3.md","16-17","On the installed 774-unit stage analysis, Figure 3E contains 1,338 / 4,266 / 6,190 /\n7,043 BH-significant genes at F1 / F2 / F3 / F4 versus F0."),
 ("R03-03","docs/manuscript/draft/fig3.md","58","F0-centered NMF activity. Each contrast uses eligible cohorts and fixed cohort"),
 ("R03-03","docs/manuscript/METHODS.md","242-245","Adjacent fibrosis stages and adjacent NAS groups were compared\ncross-sectionally in the selected five-cohort fragment-count substrate."),
 ("G2-01","docs/manuscript/05_figure_legends.md","138-141","**(E)** F1–F4-versus-F0 counts are 1,338 / 4,266 / 6,190 / 7,043 at BH q <0.05\nin each complete 24,115-gene family, without a fold-change filter. This rule\ndiffers from the canonical disease/control gate. The F0 reference contains\n103 units (14 source controls and 89 disease biopsies)."),
 ("G2-01","docs/manuscript/draft/fig3.md","63-64","depend on sample size, effect and uncertainty and cannot locate remodeling onset."),
 ("G2-01","docs/RESULTS.md","1112","the same totals are in `figures/candidates/fig3e-s3-774-20261011T165944Z/source_tables/fig3e_deg_counts.tsv`."),
 ("G2-01","docs/STATUS.md","539-541","4,266, 6,190 and 7,043 for F1 to F4, and versus NAS0 5,338, 9,337 and 9,841 for\nNAS 1-2, 3-4 and 5-8."),
 ("G2-01","docs/manuscript/METHODS.md","266","Non-estimable branches remain explicit. Early-stage GNMT and MAT1A claims that"),
 ("G2-03","docs/manuscript/05_figure_legends.md","149-152","COLOC r2 input. Four complete 24,115-gene contrast families give 216 displayed\neffects and 59 significant marks (30 / 13 / 11 / 5). Colors show\nhigher-minus-lower stage log2FC, clipped at ±0.688 for display; black dots mark\nBH significance."),
 ("G2-03","docs/manuscript/draft/fig3.md","126-127","contrast. Four complete 24,115-gene families give 216 displayed effects and 59\nsignificant marks (30 / 13 / 11 / 5)."),
 ("G2-03","docs/RESULTS.md","37-38","(65 before): 30, 13, 11 and 5 for F0→F1, F1→F2, F2→F3 and F3→F4 (32, 21, 7 and 5\nbefore)."),
 ("G2-03","docs/STATUS.md","110","157 / 94 units; ten genes leave and three enter (list in `docs/RESULTS.md`)."),
 ("R03-06","docs/manuscript/05_figure_legends.md","145-146","change. The full fibrosis/NAS, composition and F0-reference sensitivities are\nshown in Figure S3."),
 ("R03-06","docs/manuscript/05_figure_legends.md","174-175","endpoint DEG agreement, the full fibrosis/NAS remodeling panel, and hepatocyte\nand non-parenchymal composition."),
 ("R12-04","docs/manuscript/05_figure_legends.md","178-179","NAS counts use BH q <0.05 without a fold-change filter, and continuous\nNAS activity is centered on NAS0 (46 units with fixed scores)."),
 ("R12-04","docs/RESULTS.md","1106","| NAS versus NAS0 | 5–8 vs 0 | 46 / 167 | 3 | 4,571 / 5,270 |"),
 ("G2-04","docs/manuscript/draft/fig3.md","101-104","control-status term. GNMT decreases from F0 to F1 across these three fits\n(q=0.00913 / 0.0426 / 0.0430); its F1-to-F2 q is 0.0562 and is not a current\nstage call. MAT1A decreases from F0 to F1 (q=2.45e-4 / 9.5e-4 / 1.6e-3), and\nCYP2C19 decreases from F1 to F2 (q=0.0477)."),
 ("G2-04","docs/RESULTS.md","1129-1131","F0), 0.0426 (disease-only) and 0.0430 (control-status term); MAT1A: 2.45e-4,\n9.5e-4 and 1.6e-3 (`f0_arm_hero_genes.tsv`)."),
 ("G2-04","docs/PAPER.md","253-254","   disease map, GNMT joins MAT1A and CYP2C19 in lower multicohort RNA,\n   cross-sectional fibrosis-stage abundance, and liver protein."),
 ("G2-05","docs/RESULTS.md","1123","| F3→F4 | 281 | 281 | n/a |"),
]
rows = []
for fid, f, line, old in E:
    txt = open(ROOT + f, encoding="utf-8").read()
    n = txt.count(old)
    lines = txt.split("\n")
    a, b = (int(x) for x in (line.split("-") if "-" in line else (line, line)))
    pos = txt.find(old)
    start_line = txt[:pos].count("\n") + 1 if pos >= 0 else None
    end_line = start_line + old.count("\n") if pos >= 0 else None
    ok_line = pos >= 0 and start_line >= a and end_line <= b
    rows.append(dict(fid=fid, file=f, line=line, n_occ=n, found_at=f"{start_line}-{end_line}" if pos >= 0 else "NOT FOUND", within_stated_lines=ok_line))
for r in rows:
    print(json.dumps(r, ensure_ascii=False))
