#!/usr/bin/env python3
"""Genomic spans excluded from Model A het calibration and tags (reads no sample data).

MHC: chr6:28,500,000-33,500,000, the definition already used by the A1 identity panel
(a1_identity/01_targets.py:18). IG and TR loci: for each chromosome, GENCODE v49 genes
with gene_type IG_* or TR_* are clustered (genes within 1 Mb of each other join one
cluster) and each cluster's full span is excluded. Output: BED (0-based start).
"""
import sys

gtf, out = sys.argv[1], sys.argv[2]
genes = {}
for line in open(gtf):
    if line.startswith("#"):
        continue
    f = line.split("\t", 9)
    if f[2] != "gene":
        continue
    gt = f[8].split('gene_type "', 1)[1].split('"', 1)[0]
    if gt.startswith("IG_") or gt.startswith("TR_"):
        genes.setdefault((f[0], gt[:2]), []).append((int(f[3]), int(f[4])))
rows = [("chr6", 28_500_000, 33_500_000, "MHC_A1_definition")]
for (c, kind), iv in sorted(genes.items()):
    iv.sort()
    s0, e0 = iv[0]
    for s, e in iv[1:]:
        if s - e0 <= 1_000_000:
            e0 = max(e0, e)
        else:
            rows.append((c, s0 - 1, e0, f"{kind}_locus")); s0, e0 = s, e
    rows.append((c, s0 - 1, e0, f"{kind}_locus"))
with open(out, "w") as fh:
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
for r in rows:
    print(*r, f"{(r[2] - r[1]) / 1e6:.2f} Mb")
