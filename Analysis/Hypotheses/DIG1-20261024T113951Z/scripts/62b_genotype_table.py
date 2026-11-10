"""DIG1 job 62b. Genotype table and site QC from the job 62 per-sample VCFs. Written before any VCF was read.
Count-based call per sample and site (FORMAT/AD; MAPQ >= 20, base Q >= 20 in job 62):
  depth = reads for REF + target ALT + other alleles; called only if depth >= 10.
  alt fraction <= 0.10 -> 0; 0.25-0.75 -> 1; >= 0.90 -> 2; otherwise not called (allelic imbalance or error).
  For the HSD17B13 insertion: alt = reads in the A -> AA indel record; depth = reads in the SNP record at the
  same position (all reads covering it). No indel record means alt = 0.
Site QC (fixed now): call rate >= 0.90 over the 774 samples; exact HWE p >= 1e-4 in every cohort with >= 20
calls; agreement of the count-based call with the bcftools call -m genotype >= 0.95 where both exist.
Usage: python -I 62b_genotype_table.py <job62_out_dir> <new_out_dir>
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

assert len(sys.argv) == 3
G62, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=False)
T = pd.read_csv(G62 / "62_targets_annotated.tsv", sep="\t", header=None,
                names=["chrom", "start", "end", "rsid", "gene", "ref", "alt", "note"])
S = pd.read_csv(G62 / "samples.tsv", sep="\t")
key = {(r.chrom, int(r.start)): r for r in T.itertuples()}


def parse(f):
    recs = {}
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts = c[0], int(c[1]), c[3], c[4].split(",")
            fmt = dict(zip(c[8].split(":"), c[9].split(":")))
            ad = [int(x) if x != "." else 0 for x in fmt.get("AD", "0").split(",")]
            indel = "INDEL" in c[7].split(";")
            recs.setdefault((chrom, pos), []).append(dict(ref=ref, alts=alts, ad=ad, gt=fmt.get("GT", "./."), indel=indel))
    return recs


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
    a = [int(x) for x in gt.replace("|", "/").split("/")]
    idx = [i + 1 for i, x in enumerate(alts) if x == target_alt]
    return float(sum(x in idx for x in a)) if idx else 0.0


rows = []
for s in S.itertuples(index=False):
    f = G62 / "vcf" / f"{s.sample_id}.vcf.gz"
    recs = parse(f) if f.exists() else {}
    for t in T.itertuples():
        rr = recs.get((t.chrom, int(t.start)), [])
        snp = [r for r in rr if not r["indel"]]
        ind = [r for r in rr if r["indel"]]
        dp = alt_n = ref_n = 0
        gt = np.nan
        if len(t.alt) == len(t.ref):
            if snp:
                r = snp[0]; ref_n = r["ad"][0]; dp = sum(r["ad"])
                alt_n = sum(n for a, n in zip(r["alts"], r["ad"][1:]) if a == t.alt)
                gt = gt_dose(r["gt"], r["alts"], t.alt)
        else:
            if snp:
                dp = sum(snp[0]["ad"]); gt = 0.0 if "." not in snp[0]["gt"] else np.nan
            for r in ind:
                n = sum(x for a, x in zip(r["alts"], r["ad"][1:]) if len(a) - len(r["ref"]) == 1 and a[:len(r["ref"])] == r["ref"])
                alt_n += n
                gt = gt_dose(r["gt"], r["alts"], next((a for a in r["alts"] if len(a) - len(r["ref"]) == 1), "?"))
            ref_n = dp - alt_n
        frac = alt_n / dp if dp else np.nan
        rows.append(dict(sample_id=s.sample_id, dataset=s.dataset, rsid=t.rsid, gene=t.gene, depth=dp, ref_n=ref_n,
                         alt_n=alt_n, alt_frac=frac, call=call(frac, dp) if dp else np.nan, gt_bcftools=gt))
L = pd.DataFrame(rows)
L.to_csv(OUT / "genotype_long.tsv.gz", sep="\t", index=False)


def hwe_exact(n_aa, n_ab, n_bb):
    """Wigginton 2005 exact HWE test, two-sided (sum of probabilities <= observed)."""
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
for (rs, co), d in list(L.groupby(["rsid", "dataset"])) + [((rs, "ALL"), d) for rs, d in L.groupby("rsid")]:
    c = d.call.dropna()
    n0, n1, n2 = (c == 0).sum(), (c == 1).sum(), (c == 2).sum()
    both = d.dropna(subset=["call", "gt_bcftools"])
    Q.append(dict(rsid=rs, cohort=co, n=len(d), n_called=len(c), call_rate=len(c) / len(d), median_depth=d.depth.median(),
                  n0=n0, n1=n1, n2=n2, alt_freq=(n1 + 2 * n2) / (2 * len(c)) if len(c) else np.nan,
                  hwe_p=hwe_exact(n0, n1, n2) if co != "ALL" else np.nan,
                  agree_bcftools=(both.call == both.gt_bcftools).mean() if len(both) else np.nan,
                  n_ambiguous=int(((d.depth >= 10) & d.call.isna()).sum())))
Q = pd.DataFrame(Q)
site = []
for rs, d in Q.groupby("rsid"):
    a = d[d.cohort == "ALL"].iloc[0]; per = d[(d.cohort != "ALL") & (d.n_called >= 20)]
    ok = (a.call_rate >= 0.90) and bool((per.hwe_p >= 1e-4).all()) and (a.agree_bcftools >= 0.95)
    site.append(dict(rsid=rs, call_rate=a.call_rate, min_cohort_hwe_p=per.hwe_p.min(), agree_bcftools=a.agree_bcftools,
                     alt_freq=a.alt_freq, pass_qc=ok))
Q.to_csv(OUT / "site_qc_by_cohort.tsv", sep="\t", index=False)
pd.DataFrame(site).to_csv(OUT / "site_qc.tsv", sep="\t", index=False)
W = L.pivot_table(index=["sample_id", "dataset"], columns="rsid", values="call").reset_index()
W.to_csv(OUT / "genotype_calls_wide.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250)
print(pd.DataFrame(site).round(4).to_string())
print(Q.round(4).to_string())
