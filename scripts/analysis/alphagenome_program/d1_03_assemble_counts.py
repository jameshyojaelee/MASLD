#!/usr/bin/env python3
"""D1 task 3: lineage-resolved allele counts, the reliability filter, and the D1.1 feasibility verdict.

Allele parsing, indel-record rejection, the heterozygote call and the site statistic are copied
verbatim from scripts/analysis/alphagenome_atlas/73_atac_allelic_analysis.py so the lineage streams
and the deposit are read by one rule. That script is not modified.

No model is fitted here. The package delivers labels and their bounds.

Writes into <out>/tables:
  lineage_allelic_counts.tsv.gz   cohort, donor, lineage, uid, n_ref, n_alt, n_other, is_het
  lineage_site_summary.tsv        scope, lineage, uid, n_het_donors, mean_log2_alt_over_ref, se, block
  lineage_feasibility.tsv         the "which lineage supports an allele-effect label" table
  d1_feasibility.json             D1.1, D1.2, D1.3 verdicts and the diagnostics
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import pathlib
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
P3 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z"
TARGETS = P3 / "tables/atac_targets_in_peaks.tsv"
DEPOSIT_SITES = {"GSE281367": P3 / "tables/allelic_sites.tsv",
                 "GSE244832": P3 / "tables/gse244832_allelic_sites.tsv"}
STAR = "<*>"
BLOCK_BP = 1_000_000
MIN_TOTAL, MIN_ALLELE, MAX_OTHER_FRAC, FRAC_WINDOW, MIN_HET_DONORS = 10, 3, 0.10, (0.15, 0.85), 3

D1_1_CLAUSES = [("Hepatocyte", ">=", 200), ("Stellate_Cell", "<", 50),
                ("Macrophage", "<", 10), ("Cholangiocyte", "<", 10)]


# ---------------------------------------------------------------- copied from 73_atac_allelic_analysis.py
def parse_ad(ref: str, alt_field: str, ad_field: str, target_alt: str) -> tuple[int, int, int]:
    alts = alt_field.split(",")
    counts = [int(x) for x in ad_field.split(",")]
    alleles = [ref] + alts
    if len(counts) != len(alleles):
        raise ValueError(f"AD length {len(counts)} != allele count {len(alleles)} ({ref} {alt_field} {ad_field})")
    n_ref, n_alt, n_other = counts[0], 0, 0
    for allele, c in zip(alleles[1:], counts[1:]):
        if allele == target_alt:
            n_alt += c
        else:
            n_other += c
    return n_ref, n_alt, n_other


def is_indel_record(ref: str, alt_field: str) -> bool:
    if len(ref) != 1:
        return True
    return any(len(a) != 1 for a in alt_field.split(",") if a != STAR)


def call_het(n_ref: int, n_alt: int, n_other: int) -> bool:
    total = n_ref + n_alt
    if total < MIN_TOTAL or min(n_ref, n_alt) < MIN_ALLELE:
        return False
    if n_other > MAX_OTHER_FRAC * (total + n_other):
        return False
    frac = n_alt / total
    return FRAC_WINDOW[0] <= frac <= FRAC_WINDOW[1]


def site_statistic(donor_counts: list[tuple[int, int]]) -> dict:
    ratios = [math.log2((a + 0.5) / (r + 0.5)) for r, a in donor_counts]
    arr = np.asarray(ratios, dtype=float)
    n = len(arr)
    sd = float(arr.std(ddof=1)) if n > 1 else float("nan")
    return {"mean_log2": float(arr.mean()), "sd_log2": sd, "n_donors": n,
            "se_log2": (sd / math.sqrt(n)) if n > 1 else float("nan")}


# ---------------------------------------------------------------- loader
def load_targets() -> dict:
    tgt = {}
    with open(TARGETS) as h:
        for line in h:
            c, p, ref, alt = line.rstrip("\n").split("\t")
            tgt[(c, int(p))] = {"ref": ref, "alt": alt, "uid": f"{c}:{p}:{ref}:{alt}"}
    return tgt


def load_counts(counts_dir: pathlib.Path, targets: dict) -> pd.DataFrame:
    rows = []
    for f in sorted(counts_dir.glob("*.tsv.gz")):
        with gzip.open(f, "rt") as h:
            for r in csv.DictReader(h, delimiter="\t"):
                key = (r["chrom"], int(r["pos"]))
                t = targets.get(key)
                if t is None or is_indel_record(r["ref"], r["alt"]):
                    continue
                if r["ref"] != t["ref"]:
                    raise ValueError(f"pileup REF {r['ref']} != target REF {t['ref']} at {key}")
                n_ref, n_alt, n_other = parse_ad(r["ref"], r["alt"], r["ad"], t["alt"])
                rows.append({"cohort": r["cohort"], "donor": r["donor"], "lineage": r["lineage"],
                             "uid": t["uid"], "chrom": key[0], "pos": key[1],
                             "n_ref": n_ref, "n_alt": n_alt, "n_other": n_other,
                             "is_het": call_het(n_ref, n_alt, n_other)})
    return pd.DataFrame(rows)


def site_table(het: pd.DataFrame, scope_label: str) -> pd.DataFrame:
    out = []
    for (lin, uid), g in het.groupby(["lineage", "uid"], sort=True):
        if len(g) < MIN_HET_DONORS:
            continue
        st = site_statistic(list(zip(g["n_ref"].tolist(), g["n_alt"].tolist())))
        c, p, ref, alt = uid.split(":")
        out.append({"scope": scope_label, "lineage": lin, "uid": uid, "chrom": c, "pos": int(p),
                    "ref": ref, "alt": alt, "n_het_donors": st["n_donors"],
                    "mean_log2_alt_over_ref": st["mean_log2"], "sd_log2": st["sd_log2"],
                    "se_log2": st["se_log2"],
                    "total_reads": int(g["n_ref"].sum() + g["n_alt"].sum()),
                    "block": f"{c}~{int(p) // BLOCK_BP}"})
    return pd.DataFrame(out)


def main() -> None:
    out = pathlib.Path(sys.argv[1]).resolve()
    tables = out / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    targets = load_targets()
    counts = load_counts(out / "raw/counts", targets)
    if counts.empty:
        raise SystemExit("no allele counts matched the target list")
    counts.to_csv(tables / "lineage_allelic_counts.tsv.gz", sep="\t", index=False)

    het = counts[counts["is_het"]]
    frames = [site_table(het, "POOLED_30_DONORS")]
    for coh in sorted(counts["cohort"].unique()):
        frames.append(site_table(het[het["cohort"] == coh], coh))
    sites = pd.concat(frames, ignore_index=True)
    sites.to_csv(tables / "lineage_site_summary.tsv", sep="\t", index=False)

    # ---- feasibility table
    feas = []
    for (scope, lin), g in sites.groupby(["scope", "lineage"]):
        sub = het[(het["lineage"] == lin)]
        if scope != "POOLED_30_DONORS":
            sub = sub[sub["cohort"] == scope]
        cov = counts[counts["lineage"] == lin]
        if scope != "POOLED_30_DONORS":
            cov = cov[cov["cohort"] == scope]
        feas.append({
            "scope": scope, "lineage": lin,
            "n_sites_passing_filter": int(len(g)),
            "n_blocks": int(g["block"].nunique()),
            "n_donors_with_any_read": int(cov.loc[cov["n_ref"] + cov["n_alt"] > 0, "donor"].nunique()),
            "n_het_donor_site_obs": int(len(sub)),
            "median_reads_per_het_donor_site": float((sub["n_ref"] + sub["n_alt"]).median()) if len(sub) else float("nan"),
            "median_reads_per_donor_site_all": float((cov["n_ref"] + cov["n_alt"]).median()) if len(cov) else float("nan"),
            "reference_bias_mean_log2": float(g["mean_log2_alt_over_ref"].mean()),
            "supports_allele_effect_label": bool(len(g) >= 50),
        })
    # lineages with zero passing sites are missing from the groupby; add them back
    seen = {(f["scope"], f["lineage"]) for f in feas}
    for scope in sites["scope"].unique().tolist() + ["POOLED_30_DONORS"]:
        for lin in sorted(counts["lineage"].unique()):
            if (scope, lin) in seen:
                continue
            cov = counts[counts["lineage"] == lin]
            sub = het[het["lineage"] == lin]
            if scope != "POOLED_30_DONORS":
                cov, sub = cov[cov["cohort"] == scope], sub[sub["cohort"] == scope]
            feas.append({"scope": scope, "lineage": lin, "n_sites_passing_filter": 0, "n_blocks": 0,
                         "n_donors_with_any_read": int(cov.loc[cov["n_ref"] + cov["n_alt"] > 0, "donor"].nunique()),
                         "n_het_donor_site_obs": int(len(sub)),
                         "median_reads_per_het_donor_site": float((sub["n_ref"] + sub["n_alt"]).median()) if len(sub) else float("nan"),
                         "median_reads_per_donor_site_all": float((cov["n_ref"] + cov["n_alt"]).median()) if len(cov) else float("nan"),
                         "reference_bias_mean_log2": float("nan"),
                         "supports_allele_effect_label": False})
            seen.add((scope, lin))
    feasdf = pd.DataFrame(feas).sort_values(["scope", "n_sites_passing_filter"], ascending=[True, False])
    feasdf.to_csv(tables / "lineage_feasibility.tsv", sep="\t", index=False)

    # ---- D1.1
    pooled = feasdf[feasdf["scope"] == "POOLED_30_DONORS"].set_index("lineage")["n_sites_passing_filter"]
    clauses = []
    for lin, op, thr in D1_1_CLAUSES:
        n = int(pooled.get(lin, 0))
        met = (n >= thr) if op == ">=" else (n < thr)
        clauses.append({"lineage": lin, "rule": f"{op} {thr}", "observed_sites": n, "met": bool(met)})
    d11 = all(c["met"] for c in clauses)

    # ---- D1.2  ALL_READS reproduces the deposit
    d12 = []
    for coh, path in DEPOSIT_SITES.items():
        dep = pd.read_csv(path, sep="\t")
        mine = sites[(sites["scope"] == coh) & (sites["lineage"] == "ALL_READS")]
        shared = set(dep["uid"]) & set(mine["uid"])
        if shared:
            j = dep.set_index("uid").loc[sorted(shared), "mean_log2_alt_over_ref"].to_numpy(float)
            k = mine.set_index("uid").loc[sorted(shared), "mean_log2_alt_over_ref"].to_numpy(float)
            r = float(pd.Series(j).corr(pd.Series(k), method="spearman"))
            maxabs = float(np.max(np.abs(j - k)))
        else:
            r, maxabs = float("nan"), float("nan")
        d12.append({"cohort": coh, "deposit_sites": int(len(dep)), "all_reads_sites": int(len(mine)),
                    "difference": int(len(mine) - len(dep)), "shared_uids": int(len(shared)),
                    "spearman_on_shared": r, "max_abs_difference_in_mean_log2": maxabs,
                    "met": bool(abs(len(mine) - len(dep)) <= 5)})

    # ---- D1.3 reference bias sign
    d13 = []
    for _, r in feasdf[(feasdf["scope"] == "POOLED_30_DONORS") & (feasdf["n_sites_passing_filter"] >= 20)].iterrows():
        d13.append({"lineage": r["lineage"], "n_sites": int(r["n_sites_passing_filter"]),
                    "reference_bias_mean_log2": float(r["reference_bias_mean_log2"]),
                    "negative": bool(r["reference_bias_mean_log2"] < 0)})

    summary = {
        "n_targets": len(targets), "n_donors": int(counts["donor"].nunique()),
        "n_lineage_streams": int(counts["lineage"].nunique()),
        "filter": {"min_reads_per_donor_at_site": MIN_TOTAL, "min_reads_per_allele": MIN_ALLELE,
                   "max_other_allele_fraction": MAX_OTHER_FRAC,
                   "het_alt_fraction_window": list(FRAC_WINDOW), "min_het_donors": MIN_HET_DONORS},
        "D1.1": {"prediction": "Hepatocyte >= 200, Stellate_Cell < 50, Macrophage < 10, Cholangiocyte < 10 "
                               "sites clearing the reliability filter pooled over 30 donors",
                 "clauses": clauses, "met": bool(d11)},
        "D1.2": {"prediction": "the ALL_READS stream reproduces the deposit site count to within 5 sites",
                 "per_cohort": d12, "met": bool(all(x["met"] for x in d12))},
        "D1.3": {"prediction": "reference bias is negative in every lineage with >= 20 het sites",
                 "per_lineage": d13, "met": bool(all(x["negative"] for x in d13)) if d13 else None},
    }
    diagnostics = []
    hep = int(pooled.get("Hepatocyte", 0))
    if hep < 20:
        diagnostics.append(f"Hepatocyte clears only {hep} sites; the deepest lineage is under the 20-site floor")
    for row in d13:
        if abs(row["reference_bias_mean_log2"]) > 0.25:
            diagnostics.append(f"{row['lineage']} reference bias |{row['reference_bias_mean_log2']:.3f}| exceeds 0.25")
    summary["diagnostics_triggered"] = diagnostics
    json.dump(summary, (tables / "d1_feasibility.json").open("w"), indent=1, default=float)

    print(feasdf[feasdf["scope"] == "POOLED_30_DONORS"].to_string(index=False))
    print(json.dumps(summary["D1.1"], indent=1))
    print(json.dumps(summary["D1.2"], indent=1))


if __name__ == "__main__":
    main()
