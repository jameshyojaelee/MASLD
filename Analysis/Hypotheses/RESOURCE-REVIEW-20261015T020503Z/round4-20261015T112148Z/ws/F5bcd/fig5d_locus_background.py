#!/usr/bin/env python3
"""Fig 5D locus-matched background for colocalization posterior overlap with snATAC peaks.

For each adopted COLOC r2 primary signal pair (767 pairs, 3,835 pair x lineage rows):
  observed  = SNP.PP.H4 posterior mass on mapped variants inside a lineage peak
              (reproduces any/244/281/shared masses in the r2 primary-pair table)
  expected  = fraction of the pair's mapped tested variants inside a lineage peak
              (the overlap expected if the posterior were spread evenly over the
              same variants in the same cis window: matched on locus, variant
              density and local chromatin, not on LD or allele frequency)
No randomness. Writes one row per pair x lineage.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates")
R1 = ROOT / "atac-context-v3-candidate-2026-08-11-r1"
PRIMARY = ROOT / "atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv"
AUDIT = R1 / "genetics/context/variant_liftover_audit.tsv"
PEAKS = R1 / "consensus_peak_manifest.tsv"
OUT = Path(sys.argv[1])
LINEAGES = ["hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"]

prim = pd.read_csv(PRIMARY, sep="\t")
assert len(prim) == 3835
keys = prim[["gwas_name", "ensembl", "signal_pair_index"]].drop_duplicates()
assert len(keys) == 767
keyset = set(map(tuple, keys.astype(str).values))

parts = []
cols = ["gwas_name", "ensembl", "signal_pair_index", "snp_pp_h4", "mapping_status",
        "hg38_chrom", "hg38_position_1based"]
for chunk in pd.read_csv(AUDIT, sep="\t", usecols=cols, chunksize=2_000_000,
                         dtype={"gwas_name": str, "ensembl": str, "signal_pair_index": str,
                                "mapping_status": str, "hg38_chrom": str}):
    k = list(zip(chunk.gwas_name, chunk.ensembl, chunk.signal_pair_index))
    m = np.fromiter((x in keyset for x in k), bool, len(k))
    if m.any():
        parts.append(chunk.loc[m])
var = pd.concat(parts, ignore_index=True)
del parts
print("variants kept", len(var), flush=True)
var["mapped"] = var.mapping_status == "mapped"
var["pos"] = pd.to_numeric(var.hg38_position_1based, errors="coerce").fillna(0).astype(np.int64)

peaks = pd.read_csv(PEAKS, sep="\t", usecols=["lineage", "chrom", "start0", "end",
                                              "supported_gse244832", "supported_gse281367"])
for c in ["supported_gse244832", "supported_gse281367"]:
    peaks[c] = peaks[c].astype(str).str.upper() == "TRUE"


def merged(df):
    out = {}
    for chrom, g in df.groupby("chrom"):
        iv = g[["start0", "end"]].to_numpy()
        iv = iv[np.argsort(iv[:, 0])]
        s, e = [iv[0, 0]], [iv[0, 1]]
        for a, b in iv[1:]:
            if a <= e[-1]:
                e[-1] = max(e[-1], b)
            else:
                s.append(a); e.append(b)
        out[chrom] = (np.array(s), np.array(e))
    return out


def inpeak(idx, chrom, pos):
    res = np.zeros(len(pos), bool)
    for c in np.unique(chrom):
        if c not in idx:
            continue
        s, e = idx[c]
        sel = chrom == c
        p = pos[sel]
        j = np.searchsorted(s, p, side="left") - 1  # last start0 < pos
        ok = j >= 0
        r = np.zeros(len(p), bool)
        r[ok] = e[j[ok]] >= p[ok]
        res[sel] = r
    return res


idx = {}
for lin in LINEAGES:
    pl = peaks[peaks.lineage == lin]
    idx[(lin, "244")] = merged(pl[pl.supported_gse244832])
    idx[(lin, "281")] = merged(pl[pl.supported_gse281367])

chrom = var.hg38_chrom.fillna("").to_numpy()
pos = var.pos.to_numpy()
mp = var.mapped.to_numpy()
pp = var.snp_pp_h4.to_numpy()
gkey = var.gwas_name + "|" + var.ensembl + "|" + var.signal_pair_index
rows = []
for lin in LINEAGES:
    a = inpeak(idx[(lin, "244")], chrom, pos) & mp
    b = inpeak(idx[(lin, "281")], chrom, pos) & mp
    df = pd.DataFrame({"k": gkey, "pp": pp, "mapped": mp, "a": a, "b": b})
    df["any"] = df.a | df.b
    df["shared"] = df.a & df.b
    g = df.groupby("k")
    res = pd.DataFrame({
        "n_variants": g.size(),
        "n_mapped": g.mapped.sum(),
        "obs_244": g.apply(lambda x: x.pp[x.a].sum()),
        "obs_281": g.apply(lambda x: x.pp[x.b].sum()),
        "obs_shared": g.apply(lambda x: x.pp[x.shared].sum()),
        "obs_any": g.apply(lambda x: x.pp[x["any"]].sum()),
        "n_any_peak": g["any"].sum(),
        "n_shared_peak": g.shared.sum(),
    }).reset_index()
    res["exp_any"] = res.n_any_peak / res.n_mapped
    res["exp_shared"] = res.n_shared_peak / res.n_mapped
    res["lineage"] = lin
    rows.append(res)
out = pd.concat(rows, ignore_index=True)
out[["gwas_name", "ensembl", "signal_pair_index"]] = out.k.str.split("|", expand=True)
out = out.drop(columns="k")
p2 = prim.assign(signal_pair_index=prim.signal_pair_index.astype(str))
out = p2[["gwas_name", "gene", "ensembl", "trait_class", "signal_pair_index", "lineage", "pp_h4",
          "any_accessible_mass", "gse244832_accessible_mass", "gse281367_accessible_mass",
          "shared_accessible_mass", "evidence_state"]].merge(
    out, on=["gwas_name", "ensembl", "signal_pair_index", "lineage"], how="left", validate="1:1")
assert out.obs_any.notna().all()
for o, t in [("obs_any", "any_accessible_mass"), ("obs_244", "gse244832_accessible_mass"),
             ("obs_281", "gse281367_accessible_mass"), ("obs_shared", "shared_accessible_mass")]:
    d = (out[o] - out[t]).abs().max()
    print("max abs diff", o, d, flush=True)
    assert d < 1e-9, o
out.to_csv(OUT, sep="\t", index=False)
print("wrote", OUT, len(out))
