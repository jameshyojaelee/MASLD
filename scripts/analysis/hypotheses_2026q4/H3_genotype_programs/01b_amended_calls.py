#!/usr/bin/env python3
"""H3 / H31 step 1b: amended SERPINA1 calls (PRESPEC_H3.md, Amendment 1) and HSD17B13 depth checks.

Reads the step 1 output. Only the rs17580 and rs28929474 calls and pi_status change; every other
call is copied unchanged. No phenotype is read.
Usage: 01b_amended_calls.py GENO_DIR OUT_DIR   (OUT_DIR must not exist)
"""
import collections
import csv
import hashlib
import json
import math
import os
import statistics
import sys

SERP = ("rs17580", "rs28929474")
D_MIN_SERP = 100
RR_MAX, RA_MIN, RA_MAX, AA_MIN = 0.05, 0.20, 0.80, 0.95
HSD = "rs72613567"
# genotype-only counts recorded in PRESPEC_H3.md, Amendment 1
EXPECTED_PI = {"MM": 1176, "MZ": 41, "MS": 102, "SS": 2, "SZ": 1}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def call_by_fraction(alt, depth):
    if depth < D_MIN_SERP:
        return "NA"
    f = alt / depth
    if f < RR_MAX:
        return "RR"
    if RA_MIN <= f <= RA_MAX:
        return "RA"
    if f > AA_MIN:
        return "AA"
    return "NA"


def hwe(n_rr, n_ra, n_aa):
    n = n_rr + n_ra + n_aa
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    q = (n_ra + 2 * n_aa) / (2 * n)
    exp = [n * (1 - q) ** 2, n * 2 * q * (1 - q), n * q ** 2]
    obs = [n_rr, n_ra, n_aa]
    chi2 = sum((o - e) ** 2 / e for o, e in zip(obs, exp) if e > 0)
    return q, chi2, math.erfc(math.sqrt(chi2 / 2))


def main():
    geno_dir, out = sys.argv[1], sys.argv[2]
    if os.path.exists(out):
        raise SystemExit("OUT_DIR exists: " + out)
    os.makedirs(out)

    counts = {}
    with open(os.path.join(geno_dir, "person_variant_counts.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            counts[(r["person_id"], r["rsid"])] = (int(r["ref_reads"]), int(r["alt_reads"]), int(r["depth"]))

    with open(os.path.join(geno_dir, "person_genotypes.tsv")) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = list(rd.fieldnames)
        people = list(rd)

    for p in people:
        for rs in SERP:
            ref, alt, depth = counts[(p["person_id"], rs)]
            p[rs + "_call"] = call_by_fraction(alt, depth)
            p[rs + "_post"] = "NA"
        a, b = p["rs17580_call"], p["rs28929474_call"]
        if a == "NA" or b == "NA":
            detail = "NA"
        elif a == "RR" and b == "RR":
            detail = "MM"
        elif a == "RR" and b == "RA":
            detail = "MZ"
        elif a == "RA" and b == "RR":
            detail = "MS"
        elif a == "AA" and b == "RR":
            detail = "SS"
        elif a == "RA" and b == "RA":
            detail = "SZ"
        else:
            detail = "other"
        p["pi_detail"] = detail
        p["pi_status"] = detail if detail in ("MM", "MZ", "NA") else "other"

    tot = collections.Counter(p["pi_detail"] for p in people)
    print("Pi detail, all cohorts:", dict(tot))
    for k, v in EXPECTED_PI.items():
        assert tot[k] == v, "Pi class %s: %d, recorded in PRESPEC %d" % (k, tot[k], v)

    with open(os.path.join(out, "person_genotypes.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + ["pi_detail"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(people)

    # Pi counts by cohort
    cohorts = sorted({p["cohort"] for p in people})
    labels = ["MM", "MZ", "MS", "SS", "SZ", "other", "NA"]
    with open(os.path.join(out, "pi_counts_by_cohort.tsv"), "w") as fh:
        fh.write("cohort\tpersons\t" + "\t".join(labels) + "\n")
        for co in cohorts + ["ALL"]:
            sel = [p for p in people if co == "ALL" or p["cohort"] == co]
            c = collections.Counter(p["pi_detail"] for p in sel)
            fh.write("%s\t%d\t%s\n" % (co, len(sel), "\t".join(str(c[l]) for l in labels)))

    # consolidated callable-fraction table: original rows, SERPINA1 rows recomputed
    with open(os.path.join(geno_dir, "callable_fraction.tsv")) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cf_cols = list(rd.fieldnames)
        cf = list(rd)
    out_rows = []
    for r in cf:
        if r["rsid"] in SERP:
            sel = [p for p in people if p["cohort"] == r["cohort"]]
            depths = [counts[(p["person_id"], r["rsid"])][2] for p in sel]
            calls = [p[r["rsid"] + "_call"] for p in sel]
            n = collections.Counter(calls)
            q, chi2, pv = hwe(n["RR"], n["RA"], n["AA"])
            r = dict(r)
            r["persons"] = len(sel)
            r["persons_depth_ge8"] = sum(d >= 8 for d in depths)
            r["median_depth_persons"] = statistics.median(depths) if depths else "nan"
            r["called"] = n["RR"] + n["RA"] + n["AA"]
            r["n_RR"], r["n_RA"], r["n_AA"] = n["RR"], n["RA"], n["AA"]
            r["allele_freq_or_carrier_fraction"] = q
            r["hwe_chi2"], r["hwe_p"] = chi2, pv
            r["note"] = "amended rule (alt fraction), PRESPEC Amendment 1"
        out_rows.append(r)
    with open(os.path.join(out, "callable_fraction_amended.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cf_cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    # HSD17B13: depth at the site by genotype class (PRESPEC Amendment 2c)
    rows = []
    for co in cohorts:
        groups = collections.defaultdict(list)
        for p in people:
            if p["cohort"] != co:
                continue
            depth = counts[(p["person_id"], HSD)][2]
            call = p[HSD + "_call"]
            if call in ("NC", "CARRIER"):
                g = call
            elif depth >= 8:
                g = "uncalled_depth_ge8"
            else:
                g = "depth_lt8"
            groups[g].append(depth)
        for g in ("NC", "CARRIER", "uncalled_depth_ge8", "depth_lt8"):
            d = sorted(groups[g])
            if d:
                q1, q3 = d[len(d) // 4], d[(3 * len(d)) // 4]
                rows.append((co, g, len(d), statistics.median(d), q1, q3, sum(x >= 15 for x in d)))
            else:
                rows.append((co, g, 0, "NA", "NA", "NA", 0))
    with open(os.path.join(out, "hsd17b13_depth_by_class.tsv"), "w") as fh:
        fh.write("cohort\tclass\tpersons\tmedian_depth\tdepth_q25\tdepth_q75\tpersons_depth_ge15\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    info = {
        "rule": "SERPINA1: depth>=%d; RR f<%.2f; RA %.2f<=f<=%.2f; AA f>%.2f; else NA" % (D_MIN_SERP, RR_MAX, RA_MIN, RA_MAX, AA_MIN),
        "input_dir": geno_dir,
        "input_sha256": {f: sha256(os.path.join(geno_dir, f)) for f in ("person_variant_counts.tsv", "person_genotypes.tsv", "callable_fraction.tsv")},
        "pi_detail_all": dict(tot),
        "persons": len(people),
    }
    with open(os.path.join(out, "run_info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
