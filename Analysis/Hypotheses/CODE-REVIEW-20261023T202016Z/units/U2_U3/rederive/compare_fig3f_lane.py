# Compare my Fig 3F gene set, selection group, PP.H4 and trait against the r9 panel text and r9 lane table.
import csv, re, sys
run = sys.argv[1]
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
U = R + "/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U2_U3/rederive"
txt = open(U + "/r9_fig3f_panel_text.txt").read().splitlines()
panel, grp = {}, None
for line in txt:
    if "Selected at F1-F2 or later" in line: grp = "later"; continue
    if "Selected only by F0-F1" in line: grp = "only01"; continue
    m = re.match(r"^\s*(\S+)\s{5,}(GGT|AST|ALT)\s*$", line)
    if m and grp: panel[m.group(1)] = (grp, m.group(2))
lanes = {r["gene"]: r for r in csv.DictReader(open(R + "/figures/main/fig3_bulk_transcriptomics/source_tables/fig3f_exact_pair_lanes.tsv"), delimiter="\t")}
mine = {r["symbol"]: r for r in csv.DictReader(open(run + "/u2_fig3f_genes.tsv"), delimiter="\t")}
print("panel genes", len(panel), "lane rows", len(lanes), "my genes", len(mine))
print("in panel not mine", sorted(set(panel) - set(mine)), "in mine not panel", sorted(set(mine) - set(panel)))
bad = []
for g, (grp, trait) in panel.items():
    if g not in mine: continue
    m = mine[g]
    mg = "only01" if m["only_F0F1"] == "TRUE" else "later"
    mtrait = m["study"].split("_")[1]
    pp = float(m["pp4"]); lp = float(lanes[g]["PP.H4"]) if g in lanes else float("nan")
    if mg != grp or mtrait != trait or abs(pp - lp) > 1e-8 or m["study"] != lanes.get(g, {}).get("gwas"):
        bad.append((g, grp, mg, trait, mtrait, lp, pp, lanes.get(g, {}).get("gwas"), m["study"]))
print("genes differing in group/trait/PP.H4/study:", len(bad))
for b in bad: print("  ", b)
