#!/usr/bin/env python3
"""A2 step 4: annotate variants with hg38 position, TSS distance, lineage peaks, DA class,
ancestry-matched MAF and LD score.

Writes (work/):
  pos_annot.parquet   one row per hg19 position (chr, pos19): pos38, TSS distance, nearest
                      protein-coding gene, peak flags per lineage (any / both cohorts),
                      hepatocyte and stellate DA class
  anc_annot.parquet   one row per (pop, chr, pos19, allele pair): MAF, LD score
  da_peaks.tsv        jointly testable DA peaks with combined t and class
No enrichment statistic is computed here.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
A2 = ROOT / "Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage"
W = A2 / "work"
CAND = ROOT / "Analysis/Multimodal_Program_Projection/candidates"
PEAKS = CAND / "atac-context-v3-candidate-2026-08-11-r1/consensus_peak_manifest.tsv"
DA = CAND / "atac-context-v3-statistics-correction-20261022/da/da_peak_results.tsv.gz"
GENES = ROOT / ("RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/"
                "annotation/GENCODE_v49/run_19749483/gene_identity.tsv")
AF = ROOT / "GWAS/finemapping/data/ld_ref/panel_af"
LIN = ["hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"]
POPS = {"EUR": "eur", "AFR": "afr", "AMR": "amr", "EAS": "eas", "SAS": "sas"}

# ---- positions: windows + members + coloc ---------------------------------------------
win = pd.read_parquet(W / "windows.parquet", columns=["study", "chromosome", "position",
                                                      "allele1", "allele2"])
mem = pd.read_parquet(W / "cs_members.parquet")
cv = pd.read_parquet(W / "coloc_vars.parquet")
lift = pd.read_csv(W / "lift_out.tsv.gz", sep="\t")
lift["chromosome"] = lift.chrom.str.replace("chr", "", regex=False).astype(int)
lift = lift.rename(columns={"pos19": "position"})[["chromosome", "position", "pos38"]]

pos = pd.concat([win[["chromosome", "position"]].astype(np.int64),
                 mem[["chromosome", "position"]].astype(np.int64)]).drop_duplicates()
pos = pos.merge(lift, on=["chromosome", "position"], how="left", validate="1:1")
pos["src"] = "susie"
# COLOC variants carry their own recorded liftover (same chain, rtracklayer)
cvp = cv.loc[cv.mapping_status == "mapped", ["chromosome", "position", "hg38_chrom",
                                             "hg38_position_1based"]].dropna()
cvp = cvp.astype({"chromosome": np.int64, "position": np.int64})
cvp = cvp[cvp.hg38_chrom == "chr" + cvp.chromosome.astype(str)]
cvp = cvp.drop_duplicates(["chromosome", "position"]).rename(
    columns={"hg38_position_1based": "pos38"})[["chromosome", "position", "pos38"]]
both = pos.merge(cvp, on=["chromosome", "position"], how="inner", suffixes=("", "_cv"))
disagree = int((both.pos38.notna() & (both.pos38 != both.pos38_cv)).sum())
print("positions shared with COLOC liftover", len(both), "disagreeing", disagree, flush=True)
extra = cvp[~cvp.set_index(["chromosome", "position"]).index.isin(
    pos.set_index(["chromosome", "position"]).index)].assign(src="coloc")
pos = pd.concat([pos, extra], ignore_index=True)
pos = pos.dropna(subset=["pos38"]).copy()
pos["pos38"] = pos.pos38.astype(np.int64)
pos = pos.sort_values(["chromosome", "pos38"]).reset_index(drop=True)
print("annotated positions", len(pos), flush=True)


def lookup(starts, ends, p):
    """index of interval containing 1-based p (start0 < p <= end), else -1."""
    j = np.searchsorted(starts, p, side="left") - 1
    ok = j >= 0
    hit = np.full(len(p), -1, np.int64)
    hit[ok] = np.where(ends[j[ok]] >= p[ok], j[ok], -1)
    return hit


# ---- TSS -------------------------------------------------------------------------------
g = pd.read_csv(GENES, sep="\t", usecols=["gene_id_base", "gene_name", "gene_type",
                                          "chromosome", "start_1based", "end_1based",
                                          "strand", "is_canonical_chromosome"])
g = g[(g.gene_type == "protein_coding") & (g.is_canonical_chromosome.astype(str).str.lower() == "true")]
g["tss"] = np.where(g.strand == "+", g.start_1based, g.end_1based)
g["chr"] = g.chromosome.str.replace("chr", "", regex=False)
pos["tss_dist"] = np.int64(-1)
pos["nearest_gene"] = ""
for c, sub in pos.groupby("chromosome"):
    gg = g[g.chr == str(c)].sort_values("tss")
    t = gg.tss.to_numpy()
    p = sub.pos38.to_numpy()
    j = np.clip(np.searchsorted(t, p), 1, len(t) - 1)
    dl, dr = np.abs(p - t[j - 1]), np.abs(t[j] - p)
    k = np.where(dl <= dr, j - 1, j)
    pos.loc[sub.index, "tss_dist"] = np.minimum(dl, dr)
    pos.loc[sub.index, "nearest_gene"] = gg.gene_name.to_numpy()[k]

# ---- peaks -----------------------------------------------------------------------------
pk = pd.read_csv(PEAKS, sep="\t", usecols=["lineage", "peak_id", "chrom", "start0", "end",
                                          "supported_both"])
pk["both"] = pk.supported_both.astype(str).str.upper() == "TRUE"
pk["chr"] = pk.chrom.str.replace("chr", "", regex=False)
pk = pk[pk.chr.isin([str(i) for i in range(1, 23)])].reset_index(drop=True)
pk["pk_index"] = np.arange(len(pk), dtype=np.int64)
pk[["pk_index", "lineage", "peak_id", "chrom", "start0", "end", "both"]].to_parquet(
    W / "peaks_index.parquet", index=False)
for lin in LIN:
    pos[f"in_{lin}"] = False
    pos[f"inb_{lin}"] = False
    pos[f"pk_{lin}"] = np.int32(-1)
for (lin, c), sub in pk.groupby(["lineage", "chr"]):
    sub = sub.sort_values("start0")
    sel = pos.index[pos.chromosome == int(c)]
    p = pos.loc[sel, "pos38"].to_numpy()
    hit = lookup(sub.start0.to_numpy(), sub.end.to_numpy(), p)
    ok = hit >= 0
    pos.loc[sel[ok], f"in_{lin}"] = True
    pos.loc[sel[ok], f"inb_{lin}"] = sub.both.to_numpy()[hit[ok]]
    pos.loc[sel[ok], f"pk_{lin}"] = sub.pk_index.to_numpy()[hit[ok]].astype(np.int32)

# ---- DA --------------------------------------------------------------------------------
da = pd.read_csv(DA, sep="\t")
da = da[da.lineage.isin(["hepatocyte", "stellate"])].copy()
crd = da.peak_coordinate.str.extract(r"^chr(\w+):(\d+)-(\d+)$")
da["chr"], da["start0"], da["end"] = crd[0], crd[1].astype(np.int64), crd[2].astype(np.int64)
print("DA peaks dropped as non-autosomal", int((~da.chr.isin([str(i) for i in range(1, 23)])).sum()),
      flush=True)
da = da[da.chr.isin([str(i) for i in range(1, 23)])].copy()
da["t244"] = da.logFC_gse244832 / da.standard_error_gse244832
da["t281"] = da.logFC_gse281367 / da.standard_error_gse281367
da["t_c"] = (da.t244 + da.t281) / np.sqrt(2)
da["abs_tc"] = da.t_c.abs()
da["da_class"] = "middle"
for lin, sub in da.groupby("lineage"):
    hi = sub.abs_tc.quantile(0.90)
    lo = sub.abs_tc.quantile(0.50)
    da.loc[sub.index[sub.abs_tc >= hi], "da_class"] = "responsive"
    da.loc[sub.index[sub.abs_tc <= lo], "da_class"] = "stable"
da["same_sign"] = np.sign(da.logFC_gse244832) == np.sign(da.logFC_gse281367)
da["direction"] = np.where(da.t_c > 0, "up", "down")
da.drop(columns=["promoter_genes"]).to_csv(W / "da_peaks.tsv", sep="\t", index=False)
print(da.groupby(["lineage", "da_class"]).size().to_dict(), flush=True)
for lin in ["hepatocyte", "stellate"]:
    pos[f"da_{lin}"] = ""
    pos[f"dadir_{lin}"] = ""
    for c, sub in da[da.lineage == lin].groupby("chr"):
        sub = sub.sort_values("start0")
        sel = pos.index[pos.chromosome == int(c)]
        hit = lookup(sub.start0.to_numpy(), sub.end.to_numpy(), pos.loc[sel, "pos38"].to_numpy())
        ok = hit >= 0
        pos.loc[sel[ok], f"da_{lin}"] = sub.da_class.to_numpy()[hit[ok]]
        pos.loc[sel[ok], f"dadir_{lin}"] = np.where(sub.same_sign.to_numpy()[hit[ok]],
                                                    sub.direction.to_numpy()[hit[ok]], "mixed")
pos["nearest_gene"] = pos.nearest_gene.astype("category")
for lin in ["hepatocyte", "stellate"]:
    pos[f"da_{lin}"] = pos[f"da_{lin}"].astype("category")
    pos[f"dadir_{lin}"] = pos[f"dadir_{lin}"].astype("category")
pos.to_parquet(W / "pos_annot.parquet", index=False)
print("in-peak fraction by lineage",
      {l: round(float(pos[f"in_{l}"].mean()), 4) for l in LIN}, flush=True)

# ---- ancestry-matched MAF and LD score --------------------------------------------------
need = []
wv = win.merge(mem[["study", "ancestry"]].drop_duplicates(), on="study", how="left",
               validate="m:1")
need.append(wv[["ancestry", "chromosome", "position", "allele1", "allele2"]])
need.append(cv[["ancestry", "chromosome", "position", "allele1", "allele2"]].dropna())
need = pd.concat(need).drop_duplicates()
need["pop"] = need.ancestry.map(POPS)
assert need["pop"].notna().all(), need.ancestry.unique()
need = need.astype({"chromosome": np.int64, "position": np.int64})
need["allele1"] = need.allele1.astype(str).str.upper()
need["allele2"] = need.allele2.astype(str).str.upper()
outs = []
for (pop, c), sub in need.groupby(["pop", "chromosome"]):
    af = pd.read_csv(AF / pop / f"chr{c}.af.tsv.gz", sep="\t",
                     usecols=["position", "bim_a1", "bim_a2", "af_a1"])
    af = af[af.position.isin(sub.position.unique())].drop_duplicates(
        ["position", "bim_a1", "bim_a2"])
    s1 = sub.merge(af, left_on=["position", "allele1", "allele2"],
                   right_on=["position", "bim_a1", "bim_a2"], how="left")
    s2 = sub.merge(af, left_on=["position", "allele2", "allele1"],
                   right_on=["position", "bim_a1", "bim_a2"], how="left")
    a = s1.af_a1.to_numpy()
    a = np.where(np.isnan(a), s2.af_a1.to_numpy(), a)
    sub = sub.copy()
    sub["maf"] = np.minimum(a, 1 - a)
    ldf = W / "ldsc" / pop / f"chr{c}.l2.ldscore.gz"
    if not ldf.exists() and os.environ.get("A2_ALLOW_NO_LD") == "1":
        ld = pd.Series(dtype=float)  # test runs only
    else:
        ld = pd.read_csv(ldf, sep="\t", usecols=["BP", "L2"]).groupby("BP").L2.mean()
    sub["ld"] = sub.position.map(ld).to_numpy()
    outs.append(sub)
anc = pd.concat(outs, ignore_index=True)
anc.to_parquet(W / "anc_annot.parquet", index=False)
print("MAF missing by pop", anc.groupby("pop").maf.apply(lambda x: round(float(x.isna().mean()), 4)).to_dict())
print("LD missing by pop", anc.groupby("pop").ld.apply(lambda x: round(float(x.isna().mean()), 4)).to_dict())
