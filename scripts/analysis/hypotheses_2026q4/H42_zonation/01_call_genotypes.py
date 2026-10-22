#!/usr/bin/env python3
"""H42 step 1: call PNPLA3 rs738409 (I148M, ref C, alt G) genotypes per person from RNA read counts.

Calling rule is fixed in PRESPEC_H42.md. Libraries in one analysis_unit_id (conservative genetic identity
cluster) are pooled. No expression, stage or score is read here.

Usage: 01_call_genotypes.py READS_ROOT ROSTER_TSV OUT_DIR
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ERR = 0.01          # per-read error for homozygous genotypes
BB_A, BB_B = 10.0, 10.0   # beta-binomial shape for heterozygotes (mean 0.5)
MIN_DEPTH = 10
POST_MIN = 0.95
RSID = "rs738409"
SEED = 20260923


def log_binom(k, n, p):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))


def log_beta(a, b):
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def log_betabinom(k, n, a, b):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + log_beta(k + a, n - k + b) - log_beta(a, b))


def posterior(ref, alt):
    """Posterior over CC, CG, GG with a flat prior."""
    n = ref + alt
    if n == 0:
        return (1 / 3, 1 / 3, 1 / 3)
    ll = np.array([log_binom(alt, n, ERR),
                   log_betabinom(alt, n, BB_A, BB_B),
                   log_binom(ref, n, ERR)])
    ll -= ll.max()
    p = np.exp(ll)
    return tuple(p / p.sum())


def call(ref, alt, min_depth=MIN_DEPTH):
    n = ref + alt
    post = posterior(ref, alt)
    best = int(np.argmax(post))
    ok = n >= min_depth and post[best] >= POST_MIN
    return post, (best if ok else None)


def hwe_p(n0, n1, n2):
    n = n0 + n1 + n2
    if n == 0:
        return float("nan")
    q = (n1 + 2 * n2) / (2 * n)
    exp = np.array([n * (1 - q) ** 2, 2 * n * q * (1 - q), n * q ** 2])
    obs = np.array([n0, n1, n2])
    if (exp == 0).any():
        return float("nan")
    chi2 = float(((obs - exp) ** 2 / exp).sum())
    return math.erfc(math.sqrt(chi2 / 2.0))   # 1 df


def main():
    reads_root, roster_path, out_dir = map(Path, sys.argv[1:4])
    if out_dir.exists():
        raise SystemExit(f"refusing to overwrite {out_dir}")
    out_dir.mkdir(parents=True)
    roster = pd.read_csv(roster_path, sep="\t", keep_default_na=False, na_values=["", "NA"])
    roster = roster[["sample_id", "dataset", "analysis_unit_id", "primary_selected"]].copy()
    roster["primary_selected"] = roster["primary_selected"].astype(str).str.lower().eq("true")

    rows = []
    for r in roster.itertuples(index=False):
        f = reads_root / "libraries" / r.dataset / r.sample_id / "variants.tsv"
        if not f.exists():
            raise SystemExit(f"missing extraction output {f}")
        v = pd.read_csv(f, sep="\t")
        v = v[v["rsid"] == RSID]
        if len(v) != 1:
            raise SystemExit(f"{f}: expected one {RSID} row, found {len(v)}")
        v = v.iloc[0]
        ref = int(v.n_ref_plus + v.n_ref_minus)
        alt = int(v.n_alt_plus + v.n_alt_minus)
        post, g = call(ref, alt)
        rows.append(dict(sample_id=r.sample_id, dataset=r.dataset, analysis_unit_id=r.analysis_unit_id,
                         primary_selected=r.primary_selected, ref=ref, alt=alt, depth=ref + alt,
                         post_CC=post[0], post_CG=post[1], post_GG=post[2],
                         lib_dosage=np.nan if g is None else g))
    lib = pd.DataFrame(rows)
    lib.to_csv(out_dir / "h42_library_genotypes.tsv", sep="\t", index=False, na_rep="NA")

    # Pool within person
    units = []
    excluded = []
    for uid, g in lib.groupby("analysis_unit_id", sort=True):
        prim = g[g.primary_selected]
        if len(prim) == 0:
            # identity-excluded cluster (no qualified primary library); not a person in any analysis set
            excluded.append(dict(analysis_unit_id=uid, dataset=g.dataset.iloc[0], n_libraries=len(g),
                                 ref=int(g.ref.sum()), alt=int(g.alt.sum())))
            continue
        if len(prim) != 1:
            raise SystemExit(f"{uid}: {len(prim)} primary libraries")
        prim = prim.iloc[0]
        ref, alt = int(g.ref.sum()), int(g.alt.sum())
        post, d = call(ref, alt)
        post20, d20 = call(ref, alt, min_depth=20)
        called = g["lib_dosage"].dropna().astype(int).unique()
        units.append(dict(analysis_unit_id=uid, primary_library=prim.sample_id, dataset=prim.dataset,
                          n_libraries=len(g), ref=ref, alt=alt, depth=ref + alt,
                          alt_fraction=(alt / (ref + alt)) if ref + alt else np.nan,
                          post_CC=post[0], post_CG=post[1], post_GG=post[2],
                          dosage=np.nan if d is None else d,
                          dosage_depth20=np.nan if d20 is None else d20,
                          n_libraries_callable=len(g["lib_dosage"].dropna()),
                          libraries_discordant=bool(len(called) > 1)))
    u = pd.DataFrame(units)
    u["callable"] = u["dosage"].notna()
    u.to_csv(out_dir / "h42_unit_genotypes.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(excluded).to_csv(out_dir / "h42_identity_excluded_units_not_used.tsv", sep="\t", index=False, na_rep="NA")
    print(f"persons with a primary library: {len(u)}; identity-excluded clusters not used: {len(excluded)}")

    # Callable fraction and QC by cohort
    qc = []
    for ds, g in u.groupby("dataset", sort=True):
        c = g[g.callable]
        n0, n1, n2 = [(c.dosage == k).sum() for k in (0, 1, 2)]
        multi = g[g.n_libraries_callable >= 2]
        qc.append(dict(dataset=ds, persons=len(g), depth_ge10=int((g.depth >= 10).sum()),
                       callable=int(len(c)), callable_fraction=len(c) / len(g),
                       callable_depth20=int(g.dosage_depth20.notna().sum()),
                       dosage0=int(n0), dosage1=int(n1), dosage2=int(n2),
                       carriers=int(n1 + n2),
                       alt_allele_freq=((n1 + 2 * n2) / (2 * len(c))) if len(c) else np.nan,
                       hwe_p=hwe_p(n0, n1, n2), median_pooled_depth=float(g.depth.median()),
                       persons_with_2plus_callable_libraries=int(len(multi)),
                       persons_discordant=int(multi.libraries_discordant.sum())))
    pd.DataFrame(qc).to_csv(out_dir / "h42_callable_by_cohort.tsv", sep="\t", index=False, na_rep="NA")
    print(pd.DataFrame(qc).to_string(index=False))
    with open(out_dir / "run_parameters.txt", "w") as fh:
        fh.write(f"RSID={RSID}\nERR={ERR}\nBB={BB_A},{BB_B}\nMIN_DEPTH={MIN_DEPTH}\nPOST_MIN={POST_MIN}\n"
                 f"SEED={SEED}\nnumpy={np.__version__}\npandas={pd.__version__}\n")


if __name__ == "__main__":
    main()
