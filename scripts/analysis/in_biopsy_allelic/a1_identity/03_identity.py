#!/usr/bin/env python3
"""A1 step 3: identity crosswalk from RNA-called identity-panel genotypes.

Merges per-cohort VCFs, runs plink2 KING-robust kinship, chrX heterozygosity and
autosomal PCA, infers expression sex (XIST vs Y-linked genes), and writes the
library -> genetic-individual crosswalk. Kinship >= 0.354: same individual;
0.177-0.354: first degree; 0.0884-0.177: second degree. A pair is called only
when >= --min-shared-snps sites are called in both libraries. Reads no outcome.
"""
import argparse
import gzip
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y"]
PLINK2 = os.environ.get("PLINK2", "plink2")  # the plink/2.0a5.13 module names its binary "plink"
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
COUNTS = ROOT / "RNA-seq/Human/Patient_Cohorts/results/{cohort}/counts/featurecounts/gene_counts.txt"


def run(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), check=True)


def gene_names(gtf):
    names = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
            names[gid] = f[8].split('gene_name "', 1)[1].split('"', 1)[0]
    return names


def expression_sex(cohorts, names):
    rows = []
    for c in cohorts:
        fc = pd.read_csv(str(COUNTS).format(cohort=c), sep="\t", comment="#", index_col=0)
        mat = fc.iloc[:, 5:]
        cpm = mat / mat.sum(axis=0) * 1e6
        sym = pd.Series(fc.index.str.split(".").str[0], index=fc.index).map(names)
        lc = np.log2(cpm + 1)
        xist = lc[(sym == "XIST").to_numpy()].mean(axis=0)
        ymean = lc[sym.isin(Y_GENES).to_numpy()].mean(axis=0)
        for col in mat.columns:
            run_id = Path(col).name.replace(".Aligned.sortedByCoord.out.bam", "")
            rows.append((run_id, c, float(xist[col]), float(ymean[col])))
    d = pd.DataFrame(rows, columns=["run", "cohort", "xist_log2cpm", "y_log2cpm"])
    # Both markers are bimodal (log2 CPM ~0 vs 5-14). XIST and Y genes both high is
    # not a sex call: it flags possible RNA from two individuals (smoke test, GSE193066).
    x_hi, y_hi = d["xist_log2cpm"] > 5, d["y_log2cpm"] > 3
    x_lo, y_lo = d["xist_log2cpm"] < 3, d["y_log2cpm"] < 2
    d["expr_sex"] = np.select([y_hi & x_lo, x_hi & y_lo, x_hi & y_hi], ["M", "F", "both"], "ambiguous")
    return d


def contamination(merged, samples, min_dp=20):
    """Per library: fraction of reads carrying the other allele at autosomal sites
    called homozygous with DP >= min_dp. Near the sequencing error rate for a clean
    library; raised when RNA from a second individual is present."""
    q = subprocess.run(["bcftools", "query", "-t", ",".join(f"chr{i}" for i in range(1, 23)),
                        "-f", "[%GT:%AD\t]\n", str(merged)], capture_output=True, text=True, check=True).stdout
    alt_frac_num = np.zeros(len(samples))
    alt_frac_den = np.zeros(len(samples))
    n_sites = np.zeros(len(samples), dtype=int)
    for line in q.splitlines():
        for i, cell in enumerate(line.rstrip("\t").split("\t")):
            gt, _, ad = cell.partition(":")
            if gt not in ("0/0", "1/1") or ad in ("", "."):
                continue
            a = ad.split(",")
            if len(a) != 2 or "." in a:
                continue
            r, x = int(a[0]), int(a[1])
            if r + x < min_dp:
                continue
            alt_frac_num[i] += x if gt == "0/0" else r
            alt_frac_den[i] += r + x
            n_sites[i] += 1
    return pd.DataFrame({"run": samples, "hom_sites_dp20": n_sites,
                         "other_allele_fraction_at_hom": np.where(alt_frac_den > 0, alt_frac_num / np.maximum(alt_frac_den, 1), np.nan)})


def union_find(pairs, items):
    parent = {i: i for i in items}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    return {i: find(i) for i in items}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf-dir", required=True, help="restricted dir with <cohort>.vcf.gz")
    ap.add_argument("--cohorts", required=True, nargs="+")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--summary-dir", required=True, help="non-genotype outputs")
    ap.add_argument("--min-sites", type=int, default=2000)
    ap.add_argument("--min-shared-snps", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    vdir, sdir = Path(a.vcf_dir), Path(a.summary_dir)
    sdir.mkdir(parents=True, exist_ok=True)

    merged = vdir / "identity_panel.merged.vcf.gz"
    run(["bcftools", "merge", "-m", "none", "--threads", a.threads, "-Oz", "-o", merged]
        + [vdir / f"{c}.vcf.gz" for c in a.cohorts])
    run(["bcftools", "index", "-t", merged])
    pf = vdir / "identity_panel"
    # autosomes only: plink2 refuses chrX without sex; kinship and PCA use chr1-22 and the
    # chrX het rate below is read with bcftools from the merged VCF
    run([PLINK2, "--vcf", merged, "--vcf-half-call", "m", "--max-alleles", "2", "--chr", "1-22",
         "--set-all-var-ids", "@:#", "--make-pgen", "--threads", a.threads, "--out", pf])
    run([PLINK2, "--pfile", pf, "--chr", "1-22", "--missing", "sample-only", "--out", pf])
    miss = pd.read_csv(f"{pf}.smiss", sep=r"\s+")
    miss["called"] = miss["OBS_CT"] - miss["MISSING_CT"]
    keep = miss.loc[miss["called"] >= a.min_sites, "#IID"]
    keep.to_frame().rename(columns={"#IID": "#IID"}).to_csv(f"{pf}.keep", sep="\t", index=False)

    run([PLINK2, "--pfile", pf, "--keep", f"{pf}.keep", "--chr", "1-22",
         "--make-king-table", "counts", "--king-table-filter", "0.0884", "--threads", a.threads, "--out", pf])
    king = pd.read_csv(f"{pf}.kin0", sep=r"\s+")
    king = king[king["NSNP"] >= a.min_shared_snps]
    king["relation"] = pd.cut(king["KINSHIP"], [0.0884, 0.177, 0.354, 1.0],
                              labels=["second_degree", "first_degree", "same_individual"], right=False)

    if len(keep) >= 20:  # PCA is meaningless on the smoke-test handful of libraries
        run([PLINK2, "--pfile", pf, "--keep", f"{pf}.keep", "--chr", "1-22", "--maf", "0.05",
             "--indep-pairwise", "200kb", "0.2", "--out", pf])
        run([PLINK2, "--pfile", pf, "--keep", f"{pf}.keep", "--extract", f"{pf}.prune.in",
             "--pca", "10", "--threads", a.threads, "--out", pf])

    gtx = subprocess.run(["bcftools", "query", "-r", "chrX:2781480-155701382", "-f", "[%GT\t]\n", str(merged)],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    samples = subprocess.run(["bcftools", "query", "-l", str(merged)], capture_output=True, text=True,
                             check=True).stdout.split()
    g = np.array([r.rstrip("\t").split("\t") for r in gtx]) if gtx else np.empty((0, len(samples)))
    het = ((g == "0/1") | (g == "1/0")).sum(axis=0)
    called = (g != "./.").sum(axis=0)
    xhet = pd.DataFrame({"run": samples, "chrX_called": called,
                         "chrX_het_rate": np.where(called > 0, het / np.maximum(called, 1), np.nan)})

    meta = pd.read_csv(a.metadata, keep_default_na=False, na_values=[""])
    esex = expression_sex(a.cohorts, gene_names(a.gtf))
    contam = contamination(merged, samples)
    lib = (esex.merge(xhet, on="run", how="left")
               .merge(contam, on="run", how="left")
               .merge(miss[["#IID", "called"]].rename(columns={"#IID": "run", "called": "autosomal_called"}), on="run", how="left")
               .merge(meta[["sample_id", "sex"]].rename(columns={"sample_id": "run", "sex": "metadata_sex"}), on="run", how="left"))
    lib["passes_min_sites"] = lib["run"].isin(set(keep))

    same = king[king["relation"] == "same_individual"]
    runs = lib.loc[lib["passes_min_sites"], "run"].tolist()
    root = union_find(zip(same["#IID1"], same["IID2"]), runs)
    series_num = lib.set_index("run")["cohort"].str.extract(r"(\d+)")[0].astype(int)
    groups = pd.DataFrame({"run": runs, "group": [root[r] for r in runs]})
    groups["series_num"] = groups["run"].map(series_num)
    order = groups.groupby("group").agg(first_series=("series_num", "min"), first_run=("run", "min")).sort_values(
        ["first_series", "first_run"])
    ind_id = {grp: f"IND{i + 1:05d}" for i, grp in enumerate(order.index)}
    groups["individual_id"] = groups["group"].map(ind_id)
    xwalk = lib.merge(groups[["run", "individual_id"]], on="run", how="left")
    xwalk.to_csv(sdir / "sample_to_individual.tsv", sep="\t", index=False)
    king.to_csv(sdir / "kinship_related_pairs.tsv", sep="\t", index=False)

    coh = xwalk.set_index("run")["cohort"]
    king["cohort1"], king["cohort2"] = king["#IID1"].map(coh), king["IID2"].map(coh)
    dup = king[king["relation"] == "same_individual"]
    summary = {
        "libraries": int(len(xwalk)),
        "libraries_passing_min_sites": int(xwalk["passes_min_sites"].sum()),
        "genetic_individuals": int(xwalk["individual_id"].nunique()),
        "same_individual_pairs_within_cohort": {c: int(n) for c, n in dup[dup["cohort1"] == dup["cohort2"]].groupby("cohort1").size().items()},
        "same_individual_pairs_between_cohorts": {f"{c1}|{c2}": int(n) for (c1, c2), n in dup[dup["cohort1"] != dup["cohort2"]].groupby(["cohort1", "cohort2"]).size().items()},
        "first_degree_pairs": int((king["relation"] == "first_degree").sum()),
        "second_degree_pairs": int((king["relation"] == "second_degree").sum()),
        "expr_sex_vs_metadata_mismatch": int(((xwalk["metadata_sex"].isin(["M", "F"])) & (xwalk["expr_sex"].isin(["M", "F"])) & (xwalk["metadata_sex"] != xwalk["expr_sex"])).sum()),
        "expr_sex_ambiguous": int((xwalk["expr_sex"] == "ambiguous").sum()),
        "expr_sex_both_markers_high": {c: int(n) for c, n in xwalk[xwalk["expr_sex"] == "both"].groupby("cohort").size().items()},
        "other_allele_fraction_at_hom_quantiles": {str(q): float(v) for q, v in xwalk["other_allele_fraction_at_hom"].quantile([0.5, 0.9, 0.99]).items()} if xwalk["other_allele_fraction_at_hom"].notna().any() else {},
        "libraries_other_allele_fraction_gt_0.02": int((xwalk["other_allele_fraction_at_hom"] > 0.02).sum()),
        "individuals_per_cohort": {c: int(n) for c, n in xwalk.groupby("cohort")["individual_id"].nunique().items()},
    }
    (sdir / "identity_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
