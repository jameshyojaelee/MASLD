#!/usr/bin/env python3
"""Model A pivot, step T1 (one chromosome): exonic tag SNVs for GTEx v8 liver eQTL leads.

For every GTEx v8 liver eGene (qval <= 0.05) whose lead is a biallelic SNV, find
biallelic SNVs inside the gene's GENCODE v49 exons that are in LD r^2 >= 0.8 with
the lead in BOTH the unrelated EUR and the unrelated EAS founders of the 1kGP
high-coverage phased panel, with minor allele frequency >= 0.05 in both and the
same LD sign in both. r is computed from phased haplotypes (alt = 1). The sign
orients the tag: +1 means the tag's alt allele sits on the lead's alt haplotype.
The lead itself counts as a tag when it is exonic. Reads no MASLD data.

Founders (v2, 2026-09-29): every sample of the superpopulation with no recorded
parents in the pedigree (FatherID = MotherID = 0), trio parents included; children
are dropped; samples in --exclude-samples (one of each pair with KING kinship > 0.0884,
from t1_founders_king.sh) are dropped. v1 also dropped everyone listed as a parent,
which left EUR mostly FIN and removed most CEU, IBS and CHS samples (Codex review,
docs/technical/agent_exchange/2026-09-29_codex_t1_orientation_review.md).
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

MIN_R2, MIN_MAF = 0.8, 0.05


def founders(ped, superpop, exclude=frozenset()):
    p = pd.read_csv(ped, sep=r"\s+")
    f = p[(p["Superpopulation"] == superpop) & (p["FatherID"].astype(str) == "0") & (p["MotherID"].astype(str) == "0")]
    return [s for s in f["SampleID"] if s not in exclude]


def haplotypes(vcf, regions, samples):
    """dict (pos, ref, alt) -> haplotype vector (2 x n_samples, alt = 1)."""
    cmd = ["bcftools", "query", "-R", str(regions), "-s", ",".join(samples), "-i", 'TYPE="snp" && N_ALT=1',
           "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)]
    out = {}
    for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        f = line.split("\t")
        gts = f[3:]
        if any("|" not in g for g in gts):
            continue
        h = np.array([[int(g[0]), int(g[2])] for g in gts], dtype=np.int8).T.reshape(-1)
        out[(int(f[0]), f[1], f[2])] = h
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--egenes", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--ped", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude-samples", default=None, help="file of 1kGP sample IDs to drop (related pairs)")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    eg = pd.read_csv(a.egenes, sep="\t", usecols=["gene_id", "gene_name", "variant_id", "qval"])
    eg = eg[eg["qval"] <= 0.05].copy()
    v = eg["variant_id"].str.split("_", expand=True)
    eg["chrom"], eg["lead_pos"], eg["lead_ref"], eg["lead_alt"] = v[0], v[1].astype(int), v[2], v[3]
    eg = eg[(eg["chrom"] == a.chrom) & (eg["lead_ref"].str.len() == 1) & (eg["lead_alt"].str.len() == 1)]
    eg["gene_key"] = eg["gene_id"].str.split(".").str[0]

    exons = []
    for line in open(a.gtf):
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[0] != a.chrom or f[2] != "exon":
            continue
        gk = f[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
        exons.append((gk, int(f[3]), int(f[4])))
    ex = pd.DataFrame(exons, columns=["gene_key", "start", "end"]).drop_duplicates()
    ex = ex[ex["gene_key"].isin(eg["gene_key"])]

    reg = pd.concat([ex[["start", "end"]], pd.DataFrame({"start": eg["lead_pos"], "end": eg["lead_pos"]})])
    reg = reg.sort_values("start")
    reg_path = out / f"{a.chrom}.regions.tsv"
    reg.assign(chrom=a.chrom)[["chrom", "start", "end"]].to_csv(reg_path, sep="\t", header=False, index=False)

    exclude = frozenset(open(a.exclude_samples).read().split()) if a.exclude_samples else frozenset()
    haps = {}
    for pop in ("EUR", "EAS"):
        haps[pop] = haplotypes(a.panel, reg_path, founders(a.ped, pop, exclude))
    shared = sorted(k for k in haps["EUR"] if k in haps["EAS"])
    pos = np.array([k[0] for k in shared])
    exons_by_gene = {g: d[["start", "end"]].to_numpy() for g, d in ex.groupby("gene_key")}
    rows = []
    for _, L in eg.iterrows():
        lead = (L["lead_pos"], L["lead_ref"], L["lead_alt"])
        if lead not in haps["EUR"] or lead not in haps["EAS"]:
            continue
        idx = set()
        for s0, e0 in exons_by_gene.get(L["gene_key"], np.empty((0, 2), int)):
            idx.update(range(np.searchsorted(pos, s0, "left"), np.searchsorted(pos, e0, "right")))
        for key in (shared[i] for i in sorted(idx)):
            stats = {}
            for pop in ("EUR", "EAS"):
                t, l = haps[pop][key].astype(float), haps[pop][lead].astype(float)
                maf_t, maf_l = min(t.mean(), 1 - t.mean()), min(l.mean(), 1 - l.mean())
                if maf_t < MIN_MAF or maf_l < MIN_MAF or t.std() == 0 or l.std() == 0:
                    break
                stats[pop] = (float(np.corrcoef(t, l)[0, 1]), maf_t)
            if len(stats) < 2:
                continue
            (r_eur, maf_eur), (r_eas, maf_eas) = stats["EUR"], stats["EAS"]
            if r_eur ** 2 >= MIN_R2 and r_eas ** 2 >= MIN_R2 and np.sign(r_eur) == np.sign(r_eas):
                rows.append((L["gene_id"], L["gene_name"], L["variant_id"], a.chrom, key[0], key[1], key[2],
                             r_eur, r_eas, maf_eur, maf_eas, int(np.sign(r_eur)), key == lead))
    cols = ["gene_id", "gene_name", "lead_variant_id", "chrom", "pos", "ref", "alt", "r_eur", "r_eas",
            "maf_eur", "maf_eas", "orientation", "tag_is_lead"]
    t = pd.DataFrame(rows, columns=cols)
    t.to_csv(out / f"{a.chrom}.tags.tsv", sep="\t", index=False)
    print(f"{a.chrom}: eGenes with SNV leads {len(eg)}, eGenes with >= 1 tag {t['gene_id'].nunique()}, tags {len(t)}")


if __name__ == "__main__":
    main()
