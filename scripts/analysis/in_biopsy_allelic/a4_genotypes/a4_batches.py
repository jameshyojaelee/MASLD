#!/usr/bin/env python3
"""Assign genes on one chromosome to 20 leave-gene-out batches (PRESPEC_A section 6).

Each gene's mask window is its body +-100 kb. Genes are taken in start order and
each goes to the first batch (from a hash-rotated starting batch) whose last window
ends before this window starts, so windows inside a batch never overlap. Genes that
fit no batch go to an overflow list (reported; imputed in extra batches).
Writes <out>/<chrom>.batch_genes.tsv (gene_id, batch, win_start, win_end) and, per
batch, the mask regions and the +-1 Mb cis regions used to trim masked output.
"""
import hashlib
import sys
from pathlib import Path

gtf, chrom, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
FLANK, CIS, NB = 100_000, 1_000_000, 20
genes = []
for line in open(gtf):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[0] != chrom or f[2] != "gene":
        continue
    gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0]
    genes.append((int(f[3]), int(f[4]), gid))
genes.sort()
last_end = [0] * NB
assign = []
overflow = []
for s, e, g in genes:
    ws, we = max(1, s - FLANK), e + FLANK
    start = int(hashlib.sha256(g.encode()).hexdigest(), 16) % NB
    for k in range(NB):
        b = (start + k) % NB
        if last_end[b] < ws:
            last_end[b] = we
            assign.append((g, b, ws, we, s, e))
            break
    else:
        overflow.append((g, ws, we, s, e))
# overflow genes go to extra batches with the same non-overlap rule
extra_end = []
for g, ws, we, s, e in overflow:
    for i, le in enumerate(extra_end):
        if le < ws:
            extra_end[i] = we
            assign.append((g, NB + i, ws, we, s, e))
            break
    else:
        extra_end.append(we)
        assign.append((g, NB + len(extra_end) - 1, ws, we, s, e))
out.mkdir(parents=True, exist_ok=True)
with open(out / f"{chrom}.batch_genes.tsv", "w") as fh:
    fh.write("gene_id\tbatch\twin_start\twin_end\tgene_start\tgene_end\n")
    for row in assign:
        fh.write("\t".join(map(str, row)) + "\n")
n_batches = NB + len(extra_end)
for b in range(n_batches):
    rows = sorted((r for r in assign if r[1] == b), key=lambda r: r[2])
    with open(out / f"{chrom}.batch{b:02d}.mask.tsv", "w") as fh:
        for r in rows:
            fh.write(f"{chrom}\t{r[2]}\t{r[3]}\n")
    with open(out / f"{chrom}.batch{b:02d}.cis.tsv", "w") as fh:
        merged = []
        for r in rows:
            cs, ce = max(1, r[4] - CIS), r[5] + CIS
            if merged and cs <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], ce)
            else:
                merged.append([cs, ce])
        for cs, ce in merged:
            fh.write(f"{chrom}\t{cs}\t{ce}\n")
print(f"{chrom}: genes {len(genes)}, batches {n_batches} (overflow batches {len(extra_end)})")
