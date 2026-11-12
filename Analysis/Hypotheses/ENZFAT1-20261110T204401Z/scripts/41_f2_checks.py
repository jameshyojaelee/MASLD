#!/usr/bin/env python3
"""ENZFAT1 F2: orientation and build checks (Addendum 2 section 3). Exit 1 on a hard failure.
Usage: 41_f2_checks.py <out_dir>"""
import sys

import pandas as pd

out = sys.argv[1]
TR = ["bmi", "tg", "fi", "dpw", "t2d", "height"]
HARD = [("rs1558902", "A", "bmi", 1), ("rs964184", "G", "tg", 1), ("rs1260326", "T", "tg", 1),
        ("rs7903146", "T", "t2d", 1), ("rs1229984", "T", "dpw", -1)]
SOFT = [("rs2943641", "C", "fi", 1), ("rs1042725", "C", "height", 1)]
raw = pd.concat([pd.read_csv(f"{out}/checks/{t}.tsv", sep="\t") for t in TR])
rows, failed = [], False
for kind, lst in (("hard", HARD), ("soft", SOFT)):
    for snp, al, t, sign in lst:
        r = raw[(raw["SNP"] == snp) & (raw["trait"] == t)]
        if len(r) != 1:
            st = "FAIL_missing" if kind == "hard" else "missing"
            rows.append(dict(snp=snp, allele=al, trait=t, kind=kind, expected=sign, z=None, status=st))
            failed |= kind == "hard"
            continue
        r = r.iloc[0]
        s = 1.0 if r["A1"] == al else -1.0 if r["A2"] == al else 0.0
        if s == 0:
            rows.append(dict(snp=snp, allele=al, trait=t, kind=kind, expected=sign, z=None, status="FAIL_allele"))
            failed |= kind == "hard"
            continue
        z = s * r["Z"]
        if kind == "soft" and abs(z) < 2:
            st = "not_asserted_absz_lt_2"
        else:
            st = "pass" if z * sign > 0 else "FAIL"
        failed |= st == "FAIL"
        rows.append(dict(snp=snp, allele=al, trait=t, kind=kind, expected=sign, z=round(z, 2), status=st))
b = pd.concat([pd.read_csv(f"{out}/checks/{t}_build.tsv", sep="\t") for t in TR])
for _, r in b.iterrows():
    if r["trait"] in ("tg", "height", "bmi", "dpw"):
        st = "pass" if r["value"] >= 0.95 else "FAIL"
    elif r["trait"] == "fi":
        st = "pass" if r["value"] >= 900000 else "FAIL"
    else:
        st = "reported"
    failed |= st == "FAIL"
    rows.append(dict(snp="build", allele="", trait=r["trait"], kind="build", expected=None, z=r["value"],
                     status=f"{st}: {r['build_check']}"))
o = pd.DataFrame(rows)
o.to_csv(f"{out}/f2_checks.tsv", sep="\t", index=False)
print(o.to_string(index=False))
if failed:
    sys.exit("F2 check failed; traits are not munged")
