"""SS S3f step 1 (F-SS6; fork of DIG1 jobs 62b and 65). Genotype table, site QC and HSD17B13 junction readout for the
included libraries, from s3_allele_counts.py records and the S2 STAR SJ.out.tab files. Stage is not read.
Count-based call (job 62b, unchanged): depth = reads in the SNP record (REF + all ALT); called only if depth >= 10;
  alt fraction <= 0.10 -> 0; 0.25-0.75 -> 1; >= 0.90 -> 2; otherwise not called. HSD17B13 insertion as in job 62b.
Site QC (job 62b, one cohort): call rate >= 0.90 over included libraries; exact HWE p >= 1e-4, tested once over all
  included libraries; agreement with the bcftools call -m genotype >= 0.95 where both exist. If the S2 VCFs carry no
  GT (mpileup without call), agreement is not applicable and that criterion is skipped (addendum S3f).
  --dataset-term (774-unit reproduction only): HWE per dataset with >= 20 calls, as job 62b.
Calls are also made at the TNFSF10 tag sites (S3g section 45 cis test) with the same rule.
HSD17B13 junctions (job 65, unchanged): skip6 = chr4 87305309-87313822, canon67 = 87305309-87310242, STAR unique
  reads (column 7); skip_ratio = skip6 / (skip6 + canon67).
Usage: python -I s3f_genotype_table.py --input <S2 dir> --counts <s3_allele_counts dir> --targets <out_s0d_targets dir>
       --out <new dir> [--dataset-term]
"""
import argparse
import hashlib
import importlib.metadata as im
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True)
ap.add_argument("--counts", required=True)
ap.add_argument("--targets", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--dataset-term", action="store_true")
a = ap.parse_args()
IN, CN, TG, OUT = Path(a.input), Path(a.counts), Path(a.targets), Path(a.out)
OUT.mkdir(parents=True, exist_ok=False)
for line in (TG / "SHA256SUMS").read_text().splitlines():
    h, f = line.split(maxsplit=1)
    assert hashlib.sha256(Path(f).read_bytes()).hexdigest() == h, f"sha256 mismatch: {f}"
U = pd.read_csv(IN / "units.tsv", sep="\t")
U = U[U.include.astype(str).str.upper() == "TRUE"].copy()
U["sample_alias"] = U.sample_alias.astype(str)
if a.dataset_term:
    assert "dataset" in U.columns
else:
    U["dataset"] = "one_cohort"
A = pd.read_csv(TG / "s1b_targets_annotated.tsv", sep="\t")
T = A[A.role.isin(["dig1_62", "tnfsf10"])].drop_duplicates(["chrom", "pos"]).copy()
T["rsid"] = np.where(T.role == "dig1_62", T.id, "TNFSF10_tag_" + T.chrom + ":" + T.pos.astype(str))
R = pd.read_csv(CN / "geno_records.tsv.gz", sep="\t", dtype={"sample_alias": str, "gt": str, "alts": str, "ad": str})
recs = {}
for r in R.itertuples(index=False):
    ad = [int(x) if x != "." else 0 for x in str(r.ad).split(",")]
    recs.setdefault((r.sample_alias, r.chrom, int(r.pos)), []).append(
        dict(ref=r.ref, alts=str(r.alts).split(","), ad=ad, gt=r.gt if isinstance(r.gt, str) else ".", indel=bool(r.indel)))
has_gt = bool((R["gt"].astype(str) != ".").any())


def call(frac, depth):
    if depth < 10:
        return np.nan
    if frac <= 0.10:
        return 0
    if 0.25 <= frac <= 0.75:
        return 1
    if frac >= 0.90:
        return 2
    return np.nan


def gt_dose(gt, alts, target_alt):
    if "." in gt:
        return np.nan
    al = [int(x) for x in gt.replace("|", "/").split("/")]
    idx = [i + 1 for i, x in enumerate(alts) if x == target_alt]
    return float(sum(x in idx for x in al)) if idx else 0.0


rows = []
for s in U.itertuples(index=False):
    for t in T.itertuples():
        rr = recs.get((s.sample_alias, t.chrom, int(t.pos)), [])
        snp = [r for r in rr if not r["indel"]]
        ind = [r for r in rr if r["indel"]]
        dp = alt_n = ref_n = 0
        gt = np.nan
        if len(t.alt) == len(t.ref):
            if snp:
                r = snp[0]; ref_n = r["ad"][0]; dp = sum(r["ad"])
                alt_n = sum(n for x, n in zip(r["alts"], r["ad"][1:]) if x == t.alt)
                gt = gt_dose(r["gt"], r["alts"], t.alt)
        else:
            if snp:
                dp = sum(snp[0]["ad"]); gt = 0.0 if "." not in snp[0]["gt"] else np.nan
            for r in ind:
                n = sum(x for al, x in zip(r["alts"], r["ad"][1:]) if len(al) - len(r["ref"]) == 1 and al[:len(r["ref"])] == r["ref"])
                alt_n += n
                gt = gt_dose(r["gt"], r["alts"], next((al for al in r["alts"] if len(al) - len(r["ref"]) == 1), "?"))
            ref_n = dp - alt_n
        frac = alt_n / dp if dp else np.nan
        rows.append(dict(sample_alias=s.sample_alias, dataset=s.dataset, rsid=t.rsid, gene=t.gene, depth=dp, ref_n=ref_n,
                         alt_n=alt_n, alt_frac=frac, call=call(frac, dp) if dp else np.nan, gt_bcftools=gt))
L = pd.DataFrame(rows)
L.to_csv(OUT / "genotype_long.tsv.gz", sep="\t", index=False)


def hwe_exact(n_aa, n_ab, n_bb):
    """Wigginton 2005 exact HWE test, two-sided (sum of probabilities <= observed); job 62b code."""
    n_aa, n_ab, n_bb = int(n_aa), int(n_ab), int(n_bb)
    n = n_aa + n_ab + n_bb
    if n == 0:
        return np.nan
    rare = 2 * min(n_aa, n_bb) + n_ab
    hets = np.arange(rare % 2, rare + 1, 2)
    probs = np.zeros(len(hets))
    mid = int(rare * (2 * n - rare) / (2 * n))
    if mid % 2 != rare % 2:
        mid += 1
    i_mid = int(np.where(hets == mid)[0][0])
    probs[i_mid] = 1.0
    h = mid; homr = (rare - mid) // 2; homc = n - h - homr
    for i in range(i_mid, 0, -1):
        probs[i - 1] = probs[i] * h * (h - 1) / (4.0 * (homr + 1) * (homc + 1)); h -= 2; homr += 1; homc += 1
    h = mid; homr = (rare - mid) // 2; homc = n - h - homr
    for i in range(i_mid, len(hets) - 1):
        probs[i + 1] = probs[i] * 4.0 * homr * homc / ((h + 2) * (h + 1)); h += 2; homr -= 1; homc -= 1
    probs /= probs.sum()
    obs = probs[int(np.where(hets == n_ab)[0][0])]
    return float(min(1.0, probs[probs <= obs + 1e-12].sum()))


Q = []
groups = list(L.groupby(["rsid", "dataset"])) if a.dataset_term else []
for (rs, co), d in groups + [((rs, "ALL"), d) for rs, d in L.groupby("rsid")]:
    c = d.call.dropna()
    n0, n1, n2 = (c == 0).sum(), (c == 1).sum(), (c == 2).sum()
    both = d.dropna(subset=["call", "gt_bcftools"])
    Q.append(dict(rsid=rs, cohort=co, n=len(d), n_called=len(c), call_rate=len(c) / len(d), median_depth=d.depth.median(),
                  n0=n0, n1=n1, n2=n2, alt_freq=(n1 + 2 * n2) / (2 * len(c)) if len(c) else np.nan,
                  hwe_p=hwe_exact(n0, n1, n2) if (co != "ALL" or not a.dataset_term) else np.nan,
                  agree_bcftools=(both.call == both.gt_bcftools).mean() if len(both) else np.nan,
                  n_ambiguous=int(((d.depth >= 10) & d.call.isna()).sum())))
Q = pd.DataFrame(Q)
site = []
for rs, d in Q.groupby("rsid"):
    al = d[d.cohort == "ALL"].iloc[0]
    per = d[(d.cohort != "ALL") & (d.n_called >= 20)] if a.dataset_term else d[d.n_called >= 20]
    agree_ok = (al.agree_bcftools >= 0.95) if has_gt else True
    ok = (al.call_rate >= 0.90) and bool((per.hwe_p >= 1e-4).all()) and bool(agree_ok)
    site.append(dict(rsid=rs, call_rate=al.call_rate, min_cohort_hwe_p=per.hwe_p.min(), agree_bcftools=al.agree_bcftools,
                     alt_freq=al.alt_freq, gt_present=has_gt, pass_qc=ok))
Q.to_csv(OUT / "site_qc_by_cohort.tsv", sep="\t", index=False)
pd.DataFrame(site).to_csv(OUT / "site_qc.tsv", sep="\t", index=False)
W = L.pivot_table(index=["sample_alias", "dataset"], columns="rsid", values="call").reset_index()
W = U[["sample_alias", "dataset"]].merge(W, on=["sample_alias", "dataset"], how="left")
W.to_csv(OUT / "genotype_calls_wide.tsv", sep="\t", index=False)

# HSD17B13 junctions (job 65)
P = pd.read_csv(IN / "sj_paths.tsv", sep="\t", dtype=str)
P = P[P.sample_alias.isin(U.sample_alias)]
jrows, summ = [], []
for s in P.itertuples(index=False):
    skip6 = canon67 = 0
    with open(s.path) as fh:
        for line in fh:
            if not line.startswith("chr4\t"):
                continue
            c = line.split("\t"); x, y, n = int(c[1]), int(c[2]), int(c[6])
            if 87300000 <= x <= 87335000 and 87300000 <= y <= 87335000:
                jrows.append(dict(sample_alias=s.sample_alias, start=x, end=y, unique=n, annotated=int(c[5])))
                if (x, y) == (87305309, 87313822):
                    skip6 = n
                if (x, y) == (87305309, 87310242):
                    canon67 = n
    summ.append(dict(sample_alias=s.sample_alias, skip6=skip6, canon67=canon67,
                     skip_ratio=skip6 / (skip6 + canon67) if skip6 + canon67 else float("nan")))
pd.DataFrame(jrows).to_csv(OUT / "hsd17b13_junctions_long.tsv.gz", sep="\t", index=False)
J = pd.DataFrame(summ); J.to_csv(OUT / "hsd17b13_skip6.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250)
print("GT present in VCFs:", has_gt, "; libraries:", len(U), "; SJ files:", len(J))
print(pd.DataFrame(site).round(4).to_string())
print(J.describe().to_string())
