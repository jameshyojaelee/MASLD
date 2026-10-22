#!/usr/bin/env python3
"""H25 stage 1 census: callable heterozygous credible-set sites per participant and the stop rule.

Het-callable: n = ref + alt >= 10, ref >= 3, alt >= 3, other <= 10 % of ref + alt + other.
Site QC: drop blacklist sites; drop sites where >= 20 participants have n >= 10 and > 75 % of
them are het-callable. Stop rule: H25 is testable only if >= 20 high-grade (NAS >= 4) and
>= 20 low-grade (NAS <= 3) participants each have >= 10 callable CS_gained sites.
No allelic fraction, imbalance statistic or grade association is computed here.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

MIN_N, MIN_EACH, MAX_OTHER = 10, 3, 0.10
PARALOG_MIN_PART, PARALOG_MAX_HET = 20, 0.75
STOP_MIN_SITES, STOP_MIN_PART = 10, 20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--stage1-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    runs = pd.read_csv(a.runs, sep="\t")
    tg = pd.read_csv(a.targets, sep="\t")
    # The same hg38 SNV can enter twice when studies list its alleles in opposite order;
    # collapse to one row per (chrom, pos, cls) so no site is counted twice.
    tg = (tg.groupby(["chrom", "pos", "cls", "ref", "alt", "blacklist"], as_index=False)
            .agg(cs_id=("cs_id", lambda s: ";".join(sorted(set(s))))))
    assert not tg.duplicated(["chrom", "pos"]).any()
    parts, missing = [], []
    for srr in runs["srr"]:
        f = Path(a.stage1_dir) / srr / "counts.tsv"
        if not (Path(a.stage1_dir) / srr / "DONE").exists():
            missing.append(srr)
            continue
        parts.append(pd.read_csv(f, sep="\t"))
    if missing:
        sys.exit(f"stage 1 incomplete for {len(missing)} runs: {missing[:5]}")
    c = pd.concat(parts, ignore_index=True).drop_duplicates(["sample", "chrom", "pos", "cls"])
    c = c.merge(tg[["chrom", "pos", "cls", "blacklist", "cs_id"]], on=["chrom", "pos", "cls"], how="left")
    c["n"] = c["ref_n"] + c["alt_n"]
    tot = c["n"] + c["other_n"]
    c["het"] = ((c["n"] >= MIN_N) & (c["ref_n"] >= MIN_EACH) & (c["alt_n"] >= MIN_EACH)
                & (c["other_n"] <= MAX_OTHER * tot)).astype(int)
    c["deep"] = (c["n"] >= MIN_N).astype(int)

    site = c.groupby(["chrom", "pos", "cls", "cs_id", "blacklist"]).agg(
        participants_n_ge10=("deep", "sum"), participants_het=("het", "sum"),
        median_n=("n", "median"), mean_n=("n", "mean")).reset_index()
    site["het_frac_among_deep"] = site["participants_het"] / site["participants_n_ge10"].where(site["participants_n_ge10"] > 0)
    site["paralog_flag"] = ((site["participants_n_ge10"] >= PARALOG_MIN_PART)
                            & (site["het_frac_among_deep"] > PARALOG_MAX_HET)).astype(int)
    site["site_pass"] = ((site["blacklist"] == 0) & (site["paralog_flag"] == 0)).astype(int)
    site.to_csv(out / "site_census.tsv", sep="\t", index=False)

    c = c.merge(site[["chrom", "pos", "site_pass"]], on=["chrom", "pos"], how="left")
    c["callable"] = c["het"] * c["site_pass"]
    per = c.pivot_table(index="sample", columns="cls", values="callable", aggfunc="sum", fill_value=0)
    per.columns = [f"callable_{x}" for x in per.columns]
    depth = c.pivot_table(index="sample", columns="cls", values="n", aggfunc="median", fill_value=0)
    depth.columns = [f"median_n_{x}" for x in depth.columns]
    deep = c.pivot_table(index="sample", columns="cls", values="deep", aggfunc="sum", fill_value=0)
    deep.columns = [f"sites_n_ge10_{x}" for x in deep.columns]
    per = per.join(depth).join(deep).reset_index().rename(columns={"sample": "srr"})
    per = runs[["srr", "participant", "sex", "nas", "high_grade", "stage", "spots"]].merge(per, on="srr", how="left")
    per.to_csv(out / "participant_census.tsv", sep="\t", index=False)

    ok = per["callable_CS_gained"] >= STOP_MIN_SITES
    n_hi = int((ok & (per["high_grade"] == 1)).sum())
    n_lo = int((ok & (per["high_grade"] == 0)).sum())
    ok_l = per["callable_CS_lost"] >= STOP_MIN_SITES
    res = dict(
        sites_gained=int((tg["cls"] == "CS_gained").sum()), sites_lost=int((tg["cls"] == "CS_lost").sum()),
        sites_gained_pass_qc=int(site.query("cls == 'CS_gained'")["site_pass"].sum()),
        sites_lost_pass_qc=int(site.query("cls == 'CS_lost'")["site_pass"].sum()),
        participants=int(len(per)), high_grade=int(per["high_grade"].sum()),
        low_grade=int((per["high_grade"] == 0).sum()),
        high_grade_with_ge10_callable_gained=n_hi, low_grade_with_ge10_callable_gained=n_lo,
        high_grade_with_ge10_callable_lost=int((ok_l & (per["high_grade"] == 1)).sum()),
        low_grade_with_ge10_callable_lost=int((ok_l & (per["high_grade"] == 0)).sum()),
        callable_gained_median=float(per["callable_CS_gained"].median()),
        callable_gained_q25_q75=[float(per["callable_CS_gained"].quantile(q)) for q in (0.25, 0.75)],
        callable_gained_max=int(per["callable_CS_gained"].max()),
        callable_lost_median=float(per["callable_CS_lost"].median()),
        median_depth_gained_sites=float(per["median_n_CS_gained"].median()),
        rule=f">= {STOP_MIN_PART} participants per grade group with >= {STOP_MIN_SITES} callable CS_gained sites",
        testable=bool(n_hi >= STOP_MIN_PART and n_lo >= STOP_MIN_PART),
    )
    json.dump(res, open(out / "stop_rule.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
