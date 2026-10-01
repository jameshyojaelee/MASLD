#!/usr/bin/env python3
"""A3 diagnostics after gate G1 failed (the gate itself is unchanged).

1. RNA-call errors split by type and by whether the site is in the imputation
   panel: false het (RNA het, DNA hom-ref), missed het (RNA hom, DNA het), other.
2. Non-reference concordance at panel sites only, by RNA depth and site QUAL (no per-sample GQ is written by bcftools call).
3. Imputed dosage vs DNA at liver lead SNVs: mean per-site r^2 by DNA MAF in the 30.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def q(vcf, fmt, samples=None):
    cmd = ["bcftools", "query", "-f", fmt] + (["-s", ",".join(samples)] if samples else []) + [str(vcf)]
    return [l.split("\t") for l in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines()]


def dos(g):
    return np.nan if g.startswith(".") else sum(int(x) for x in g.replace("|", "/").split("/"))


w = Path(sys.argv[1]) / "work"
eg = pd.read_csv(sys.argv[2], sep="\t", usecols=["variant_id"])["variant_id"].str.split("_", expand=True)
eg = eg[(eg[2].str.len() == 1) & (eg[3].str.len() == 1)]
rows, site_r2 = [], []
for n in range(1, 23):
    c = f"chr{n}"
    samples = subprocess.run(["bcftools", "query", "-l", str(w / f"calls/{c}.rna.vcf.gz")], capture_output=True, text=True, check=True).stdout.split()
    dna = {(int(r[0]), r[1], r[2]): r[3:] for r in q(w / f"truth/{c}.dna.vcf.gz", "%POS\t%REF\t%ALT[\t%GT]\n", samples)}
    panel = {(int(r[0]), r[1], r[2]) for r in q(w / f"panel/{c}.ref.vcf.gz", "%POS\t%REF\t%ALT\n")}
    for r in q(w / f"calls/{c}.rna.vcf.gz", "%POS\t%REF\t%ALT\t%QUAL[\t%GT:%DP]\n"):
        key = (int(r[0]), r[1], r[2]); truth = dna.get(key); inp = key in panel
        for i, cell in enumerate(r[4:]):
            gt, dp = cell.split(":")
            gq = r[3]
            g = dos(gt)
            if np.isnan(g):
                continue
            t = dos(truth[i]) if truth else 0.0
            rows.append((inp, int(dp) if dp not in (".", "") else 0, float(gq) if gq not in (".", "") else 0.0, g, t))
    imp = {(int(r[0]), r[1], r[2]): r[3:] for r in q(w / f"imputed/{c}.standard.vcf.gz", "%POS\t%REF\t%ALT[\t%DS]\n")}
    for _, L in eg[eg[0] == c].iterrows():
        key = (int(L[1]), L[2], L[3])
        if key in imp and key in dna:
            t = np.array([dos(x) for x in dna[key]], float); d = np.array(imp[key], float)
            maf = min(np.nanmean(t) / 2, 1 - np.nanmean(t) / 2)
            if np.nanstd(t) > 0 and np.nanstd(d) > 0:
                site_r2.append((maf, np.corrcoef(d, t)[0, 1] ** 2))
x = pd.DataFrame(rows, columns=["in_panel", "dp", "gq", "rna", "dna"])
x["type"] = np.select([(x.rna == 1) & (x.dna == 0), (x.rna != 1) & (x.dna == 1), x.rna != x.dna], ["false_het", "missed_het", "other"], "concordant")
nonref = (x.rna > 0) | (x.dna > 0)
out = {"error_types_by_panel": x[x.type != "concordant"].groupby(["in_panel", "type"]).size().unstack(fill_value=0).to_dict(),
       "nonref_concordance_panel_sites": float((x[nonref & x.in_panel].rna == x[nonref & x.in_panel].dna).mean()),
       "nonref_concordance_non_panel_sites": float((x[nonref & ~x.in_panel].rna == x[nonref & ~x.in_panel].dna).mean()),
       "fraction_nonref_genotypes_at_non_panel_sites": float((nonref & ~x.in_panel).sum() / nonref.sum())}
for lo, hi in [(10, 20), (20, 50), (50, 10**9)]:
    for gq in (30, 100, 300):
        m = nonref & x.in_panel & x.dp.between(lo, hi - 1) & (x.gq >= gq)
        out[f"nonref_conc_panel_dp{lo}-{hi}_siteQUAL>={gq}"] = [float((x[m].rna == x[m].dna).mean()), int(m.sum())]
s = pd.DataFrame(site_r2, columns=["maf", "r2"])
s["maf_bin"] = pd.cut(s.maf, [0, 0.05, 0.1, 0.2, 0.5])
out["lead_snv_mean_per_site_r2_by_maf"] = {str(k): [float(v.r2.mean()), int(len(v))] for k, v in s.groupby("maf_bin", observed=True)}
(Path(sys.argv[1]) / "g1_diagnostics.json").write_text(json.dumps(out, indent=2, default=str))
print(json.dumps(out, indent=2, default=str))
