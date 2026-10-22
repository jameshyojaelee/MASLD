"""Gate L02-02: map canonical bulk calls to GSE267119 H3K27ac regions near the TSS.

Reads only saved files. Region keys are 1-based inclusive (dataset toml blockers[0]).
Outputs go to this scratch directory.
"""
import gzip
import bisect
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-02"
H3 = ROOT + "/Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE267119_H3K27ac.txt.gz"
JOIN = ROOT + "/Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"
DEG = ROOT + "/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
GTF = "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
SUPP3 = ROOT + "/Analysis/Hypotheses/H25-20261009T020111Z/inputs/supp3_liver_NORvsNASH_regions.tsv"
WIN = 1000

# ---- H3K27ac matrix: header (quoted, space separated), then key + 99 counts
with gzip.open(H3, "rt") as fh:
    header = [h.strip('"') for h in fh.readline().split()]
    keys, rows = [], []
    for line in fh:
        parts = line.split()
        keys.append(parts[0].strip('"'))
        rows.append(np.array(parts[1:], dtype=np.int64))
counts = np.vstack(rows)
print("matrix regions", counts.shape[0], "samples", counts.shape[1], "header", len(header))
pid = [h.split("_")[0] for h in header]
join = pd.read_csv(JOIN, sep="\t", dtype=str)
print("header participants in join:", sum(p in set(join.participant_id) for p in pid), "of", len(pid))
# label in header vs join stage
lab_hdr = {h.split("_")[0]: "_".join(h.split("_")[1:-2]).upper() for h in header}
j = join.set_index("participant_id")
mism = [p for p in pid if lab_hdr[p].replace("NASH_F23", "NASH_F23") != j.loc[p, "stage"]]
print("header stage label mismatches vs join:", len(mism), mism[:5])

reg = pd.DataFrame({"key": keys})
m = reg.key.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
reg["chrom"], reg["start"], reg["end"] = m[0], m[1].astype(int), m[2].astype(int)
reg["width"] = reg.end - reg.start + 1
reg["total"] = counts.sum(axis=1)
libsize = counts.sum(axis=0)
print("region width median", int(reg.width.median()), "IQR", reg.width.quantile([.25, .75]).tolist())
print("library size (counts in regions) min/median/max", libsize.min(), int(np.median(libsize)), libsize.max())
stage = np.array([j.loc[p, "stage"] for p in pid])
sex = np.array([j.loc[p, "sex"] for p in pid])
dis = stage != "NOR"
print("disease", dis.sum(), "NOR", (~dis).sum())
for s in ["NOR", "NAFL", "NASH_F0", "NASH_F1", "NASH_F23"]:
    ls = libsize[stage == s]
    print(f"  libsize {s}: n={len(ls)} median={int(np.median(ls))}")

# ---- GENCODE v49 TSS: gene-level and transcript-level, primary chromosomes
gene_tss, tx_tss, gene_type, gene_name = {}, defaultdict(set), {}, {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] not in ("gene", "transcript"):
            continue
        if not re.match(r"^chr([0-9]+|X|Y)$", f[0]):
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        tss = int(f[3]) if f[6] == "+" else int(f[4])
        if f[2] == "gene":
            gene_tss[gid] = (f[0], tss, f[6])
            gene_type[gid] = re.search(r'gene_type "([^"]+)"', f[8]).group(1)
            gene_name[gid] = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
        else:
            tx_tss[gid].add((f[0], tss))
print("GTF primary-chromosome genes", len(gene_tss))

# ---- interval index per chromosome
idx = {}
for c, g in reg.groupby("chrom"):
    g = g.sort_values("start")
    idx[c] = (g.start.values, g.end.values, g.index.values)
maxw = int(reg.width.max())


def regions_near(chrom, pos, win=WIN):
    if chrom not in idx:
        return []
    st, en, ix = idx[chrom]
    lo = bisect.bisect_left(st, pos - win - maxw)
    hi = bisect.bisect_right(st, pos + win)
    return [ix[k] for k in range(lo, hi) if en[k] >= pos - win and st[k] <= pos + win]


deg = pd.read_csv(DEG, keep_default_na=False, na_values=[""])
deg["canonical"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
print("deg genes", len(deg), "canonical", deg.canonical.sum(),
      "up", (deg.canonical & (deg.logFC > 0)).sum(), "down", (deg.canonical & (deg.logFC < 0)).sum())
deg["in_gtf"] = deg.gene.isin(gene_tss)
print("deg genes with primary-chromosome GTF TSS", deg.in_gtf.sum(), "; canonical:", (deg.canonical & deg.in_gtf).sum())

recs = []
for r in deg.itertuples():
    if r.gene not in gene_tss:
        continue
    c, t, s = gene_tss[r.gene]
    g_regs = regions_near(c, t)
    tx_regs = set()
    for (cc, tt) in tx_tss.get(r.gene, []):
        tx_regs.update(regions_near(cc, tt))
    best = None
    if g_regs:
        best = max(g_regs, key=lambda k: reg.total.iat[k])
    recs.append(dict(gene=r.gene, symbol=r.symbol, logFC=r.logFC, padj=r.padj, AveExpr=r.AveExpr,
                     canonical=r.canonical, gene_type=gene_type[r.gene],
                     n_reg_gene_tss=len(g_regs), n_reg_any_tx_tss=len(tx_regs),
                     best_region=reg.key.iat[best] if best is not None else "",
                     best_region_mean_count=reg.total.iat[best] / counts.shape[1] if best is not None else np.nan))
gm = pd.DataFrame(recs)
gm.to_csv(OUT + "/gene_promoter_map.tsv", sep="\t", index=False)
can = gm[gm.canonical]
hit = can[can.n_reg_gene_tss > 0]
print(f"canonical with region within {WIN} bp of gene-level TSS: {len(hit)} (up {(hit.logFC>0).sum()}, down {(hit.logFC<0).sum()})")
hit2 = can[can.n_reg_any_tx_tss > 0]
print(f"canonical with region within {WIN} bp of any transcript TSS: {len(hit2)} (up {(hit2.logFC>0).sum()}, down {(hit2.logFC<0).sum()})")
print("canonical gene types:", can.gene_type.value_counts().head(6).to_dict())
print("canonical-with-promoter gene types:", hit.gene_type.value_counts().head(6).to_dict())
nonc = gm[(~gm.canonical) & (gm.n_reg_gene_tss > 0)]
print("non-call genes with promoter region (negative-control pool):", len(nonc))
# region sharing among hit genes
shared = hit.best_region.value_counts()
print("distinct best regions among hit genes:", shared.size, "; regions used by >1 gene:", (shared > 1).sum(),
      "genes on shared regions:", int(shared[shared > 1].sum()))
print("best-region mean count per sample among hit genes: quantiles",
      hit.best_region_mean_count.quantile([0.05, 0.25, 0.5, 0.75]).round(1).tolist())
print("hit genes with best-region mean count < 10:", int((hit.best_region_mean_count < 10).sum()))

# ---- positive control COL1A1
for sym in ["COL1A1", "COL1A2", "COL3A1"]:
    row = gm[gm.symbol == sym]
    print(sym, row[["gene", "logFC", "padj", "n_reg_gene_tss", "n_reg_any_tx_tss", "best_region", "best_region_mean_count"]].to_dict("records"))
    if len(row) and row.best_region.iat[0]:
        c, t, s = gene_tss[row.gene.iat[0]]
        print("   TSS", c, t, s, "regions within 1kb:", [reg.key.iat[k] for k in regions_near(c, t)])
        k = reg.index[reg.key == row.best_region.iat[0]][0]
        v = counts[k]
        cpm = v / libsize * 1e6
        print("   per-group median counts:", {s2: float(np.median(v[stage == s2])) for s2 in ["NOR", "NAFL", "NASH_F0", "NASH_F1", "NASH_F23"]})
        print("   zero counts in samples:", int((v == 0).sum()))

# source supplementary table S3 (NASH_F2/3 vs NOR, DESeq2, FDR<=0.01, |log2FC|>0.5)
s3 = pd.read_csv(SUPP3, sep="\t")
s3["key"] = s3.seqnames + ":" + s3.start.astype(str) + "-" + s3.end.astype(str)
print("supp S3 rows", len(s3), "keys found in deposited matrix:", s3.key.isin(set(reg.key)).sum(),
      "gained", (s3.nashvsnor_log2FoldChange > 0).sum(), "lost", (s3.nashvsnor_log2FoldChange < 0).sum())
for sym in ["COL1A1", "COL1A2", "COL3A1", "LUM", "THBS2"]:
    row = gm[gm.symbol == sym]
    if len(row) == 0:
        continue
    c, t, s = gene_tss[row.gene.iat[0]]
    near = [reg.key.iat[k] for k in regions_near(c, t)]
    sub = s3[s3.key.isin(near)]
    print("S3", sym, sub[["key", "nashvsnor_log2FoldChange", "nashvsnor_padj"]].to_dict("records"))
# how many of the candidate's promoter regions are already source S3 differential (F2/3 vs NOR)
hs = set(s3.key)
hit_in_s3 = hit.best_region.isin(hs)
print("hit genes whose best promoter region is in source S3:", int(hit_in_s3.sum()), "of", len(hit))
sub = hit[hit_in_s3].merge(s3[["key", "nashvsnor_log2FoldChange"]], left_on="best_region", right_on="key")
print("  of those, same sign as bulk logFC:", int((np.sign(sub.logFC) == np.sign(sub.nashvsnor_log2FoldChange)).sum()))
np.save(OUT + "/tmp/libsize.npy", libsize)
