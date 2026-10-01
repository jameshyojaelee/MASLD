#!/usr/bin/env python3
"""A3 pilot gate G1 (PRESPEC_A section 7) on Geuvadis RNA-seq against 1kGP DNA.

1. RNA-called genotypes vs DNA: non-reference concordance at DP >= 10; het
   sensitivity at DNA-het sites by RNA depth and by RNA allelic fraction.
2. Imputed dosage vs DNA at GTEx v8 liver lead SNVs (MAF >= 0.05 in the 30), by
   distance to the nearest RNA-called site: aggregate r^2 = squared Pearson
   correlation of pooled imputed and true dosages.
3. Leave-gene-out imputation: r^2 at RNA-called sites inside masked gene windows.
4. STAR WASP mode with two-pass: vW tags present in the 5 WASP libraries.
"""
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

THRESH = {"nonref_concordance": 0.97, "het_sensitivity_dp20": 0.90, "r2_leads_within_100kb": 0.80}


def query(vcf, fmt, region=None, samples=None):
    cmd = ["bcftools", "query", "-f", fmt]
    if samples:
        cmd += ["-s", ",".join(samples)]
    if region:
        cmd += ["-r", region]
    out = subprocess.run(cmd + [str(vcf)], capture_output=True, text=True, check=True).stdout
    return [l.split("\t") for l in out.rstrip("\n").split("\n") if l]


def gt_dosage(g):
    if g in ("./.", ".|.", "."):
        return np.nan
    return sum(int(a) for a in g.replace("|", "/").split("/") if a != ".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exec-dir", required=True)
    ap.add_argument("--gtex-egenes", required=True)
    ap.add_argument("--geuvadis-bam", required=True)
    a = ap.parse_args()
    w = Path(a.exec_dir) / "work"
    inds = [l.strip() for l in open(w / "truth/individuals.txt") if l.strip()]
    sel = pd.read_csv(Path(a.geuvadis_bam).parent / "selected_30_eur.tsv", sep="\t")
    run_to_ind = dict(zip(sel["run_accession"], sel["individual"]))

    conc_rows, het_rows, lead_rows, mask_rows = [], [], [], []
    eg = pd.read_csv(a.gtex_egenes, sep="\t", usecols=["variant_id"])
    lead = eg["variant_id"].str.split("_", expand=True)
    lead = lead[(lead[2].str.len() == 1) & (lead[3].str.len() == 1)]
    lead.columns = ["chrom", "pos", "ref", "alt", "build"]
    lead["pos"] = lead["pos"].astype(int)

    for n in range(1, 23):
        c = f"chr{n}"
        rna_samples = subprocess.run(["bcftools", "query", "-l", str(w / f"calls/{c}.rna.vcf.gz")],
                                     capture_output=True, text=True, check=True).stdout.split()
        # RNA sample names are the BAM read-group SM = individual IDs
        rna = query(w / f"calls/{c}.rna.vcf.gz", "%POS\t%REF\t%ALT[\t%GT:%DP:%AD]\n")
        dna = {(int(r[0]), r[1], r[2]): r[3:] for r in query(w / f"truth/{c}.dna.vcf.gz", "%POS\t%REF\t%ALT[\t%GT]\n", samples=rna_samples)}
        called_pos = []
        for r in rna:
            key = (int(r[0]), r[1], r[2])
            called_pos.append(key[0])
            truth = dna.get(key)
            for i, cell in enumerate(r[3:]):
                gt, dp, ad = cell.split(":")
                t = gt_dosage(truth[i]) if truth else 0.0  # absent from the DNA SNV set -> hom-ref in DNA
                g = gt_dosage(gt)
                dpv = int(dp) if dp not in (".", "") else 0
                if not np.isnan(g) and dpv >= 10:
                    conc_rows.append((g, t))
                if t == 1 and dpv >= 10:
                    x = ad.split(",")
                    af = int(x[1]) / max(1, int(x[0]) + int(x[1])) if len(x) == 2 and "." not in x else np.nan
                    het_rows.append((dpv, af, g == 1))
        called = np.array(sorted(set(called_pos)))
        # imputed dosages at liver lead SNVs
        ll = lead[lead["chrom"] == c]
        imp = {(int(r[0]), r[1], r[2]): r[3:] for r in query(w / f"imputed/{c}.standard.vcf.gz", "%POS\t%REF\t%ALT[\t%DS]\n")}
        for _, L in ll.iterrows():
            key = (L["pos"], L["ref"], L["alt"])
            if key not in imp or key not in dna:
                continue
            truth = np.array([gt_dosage(x) for x in dna[key]], float)
            if np.nanmean(truth) / 2 < 0.05 or np.nanmean(truth) / 2 > 0.95:
                continue
            ds = np.array(imp[key], float)
            dist = int(np.min(np.abs(called - L["pos"]))) if len(called) else 10**9
            lead_rows.append((c, L["pos"], dist, ds, truth))
        # masked sites
        masked_ids = set(l.strip() for l in open(w / f"imputed/{c}.excluded_markers.txt"))
        mimp = query(w / f"imputed/{c}.masked.vcf.gz", "%ID\t%POS\t%REF\t%ALT[\t%DS]\n")
        for r in mimp:
            key = (int(r[1]), r[2], r[3])
            if f"{c}:{r[1]}:{r[2]}:{r[3]}" in masked_ids and key in dna:
                truth = np.array([gt_dosage(x) for x in dna[key]], float)
                mask_rows.append((np.array(r[4:], float), truth))

    conc = np.array(conc_rows)
    nonref = (conc[:, 0] > 0) | (conc[:, 1] > 0)
    het = pd.DataFrame(het_rows, columns=["dp", "rna_af", "called_het"])
    het["dp_bin"] = pd.cut(het["dp"], [10, 20, 50, 10**9], right=False, labels=["10-19", "20-49", ">=50"])
    het["imbalance_bin"] = pd.cut((het["rna_af"] - 0.5).abs(), [0, 0.1, 0.2, 0.3, 0.4, 0.51], include_lowest=True)

    def agg_r2(rows):
        if not rows:
            return None, 0
        ds = np.concatenate([r[-2] for r in rows]); tr = np.concatenate([r[-1] for r in rows])
        ok = np.isfinite(ds) & np.isfinite(tr)
        return float(np.corrcoef(ds[ok], tr[ok])[0, 1] ** 2), len(rows)

    bins = [(0, 0), (1, 10_000), (10_001, 50_000), (50_001, 100_000), (100_001, 10**9)]
    lead_r2 = {}
    for lo, hi in bins:
        r2, n = agg_r2([r for r in lead_rows if lo <= r[2] <= hi])
        lead_r2[f"{lo}-{hi}"] = {"r2": r2, "n_sites": n}
    r2_100, n_100 = agg_r2([r for r in lead_rows if r[2] <= 100_000])
    mask_r2, mask_n = agg_r2(mask_rows)

    wasp = {}
    for f in sorted(Path(a.geuvadis_bam).glob("*/*.wasp_tag_counts.txt")):
        wasp[f.name.split(".")[0]] = dict(l.split() for l in open(f) if l.strip())

    res = {
        "rna_calls": {"nonref_concordance_dp10": float((conc[nonref, 0] == conc[nonref, 1]).mean()),
                      "overall_concordance_dp10": float((conc[:, 0] == conc[:, 1]).mean()),
                      "n_genotypes": int(len(conc))},
        "het_sensitivity": {"dp>=20": float(het.loc[het["dp"] >= 20, "called_het"].mean()),
                            "by_depth": het.groupby("dp_bin", observed=True)["called_het"].mean().round(4).to_dict(),
                            "by_rna_imbalance": {str(k): float(v) for k, v in het[het["dp"] >= 20].groupby("imbalance_bin", observed=True)["called_het"].mean().items()}},
        "imputation_liver_leads": {"by_distance_bp": lead_r2, "r2_within_100kb": r2_100, "n_within_100kb": n_100},
        "leave_gene_out_masked_sites": {"r2": mask_r2, "n_sites": mask_n},
        "wasp_vw_tags": wasp,
    }
    res["gate"] = {
        "nonref_concordance": res["rna_calls"]["nonref_concordance_dp10"] >= THRESH["nonref_concordance"],
        "het_sensitivity_dp20": res["het_sensitivity"]["dp>=20"] >= THRESH["het_sensitivity_dp20"],
        "r2_leads_within_100kb": (r2_100 or 0) >= THRESH["r2_leads_within_100kb"],
        "wasp_two_pass": bool(wasp) and all(int(v.get("reads_with_vW", 0)) > 0 for v in wasp.values()),
    }
    res["gate"]["pass"] = all(res["gate"].values())
    out = Path(a.exec_dir) / "gate_g1.json"
    out.write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
