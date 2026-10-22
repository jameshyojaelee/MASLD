#!/usr/bin/env python3
"""H3 / H31 step 1: person-level RNA-derived genotypes from the READS extraction.

Rule is fixed in PRESPEC_H3.md. No phenotype, stage or expression value is read here.

Usage: 01_person_genotypes.py OUT_DIR   (OUT_DIR must not exist)
"""
import csv
import gzip
import json
import math
import os
import sys
from collections import defaultdict, Counter

import numpy as np

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
READS = f"{REPO}/Analysis/Hypotheses/READS-20261009T141017Z/libraries"
N1 = f"{REPO}/Analysis/Hypotheses/N1-20261007T210855Z"
SEED = 20260923
D_MIN = 8
ERR = 0.01
BB_A = BB_B = 5.0
POST_MIN = 0.95
CARRIER_GRID = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
MIN_COHORT_CALLABLE = 20

VARIANTS = ["rs72613567", "rs738409", "rs58542926", "rs2642438", "rs1260326", "rs17580", "rs28929474"]
GENE = {"rs72613567": "HSD17B13", "rs738409": "PNPLA3", "rs58542926": "TM6SF2", "rs2642438": "MTARC1",
        "rs1260326": "GCKR", "rs17580": "SERPINA1", "rs28929474": "SERPINA1"}
TWO_CLASS = {"rs72613567"}

out_dir = sys.argv[1]
if os.path.exists(out_dir):
    sys.exit(f"refusing to overwrite {out_dir}")
os.makedirs(out_dir)
np.random.seed(SEED)


def log_binom(a, d, p):
    a = np.asarray(a, float)
    d = np.asarray(d, float)
    lc = (np.vectorize(math.lgamma)(d + 1) - np.vectorize(math.lgamma)(a + 1) - np.vectorize(math.lgamma)(d - a + 1))
    return lc + a * math.log(p) + (d - a) * math.log1p(-p)


def log_betabinom(a, d, al, be):
    a = np.asarray(a, float)
    d = np.asarray(d, float)
    lg = np.vectorize(math.lgamma)
    lc = lg(d + 1) - lg(a + 1) - lg(d - a + 1)
    lb = lambda x, y: math.lgamma(x) + math.lgamma(y) - math.lgamma(x + y)
    return lc + lg(a + al) + lg(d - a + be) - lg(d + al + be) - lb(al, be)


def logsumexp(arrs):
    m = np.max(arrs, axis=0)
    return m + np.log(np.sum(np.exp(arrs - m), axis=0))


def three_class_loglik(a, d):
    return np.vstack([log_binom(a, d, ERR), log_betabinom(a, d, BB_A, BB_B), log_binom(a, d, 1 - ERR)])


def two_class_loglik(a, d):
    l0 = log_binom(a, d, ERR)
    l1 = logsumexp(np.vstack([log_binom(a, d, p) for p in CARRIER_GRID])) - math.log(len(CARRIER_GRID))
    return np.vstack([l0, l1])


def em_three(ll):
    q = 0.2
    for _ in range(100):
        pri = np.log(np.array([(1 - q) ** 2, 2 * q * (1 - q), q ** 2]))[:, None]
        lp = ll + pri
        post = np.exp(lp - logsumexp(lp))
        qn = float(np.clip((post[1].sum() + 2 * post[2].sum()) / (2 * post.shape[1]), 0.001, 0.999))
        if abs(qn - q) < 1e-9:
            q = qn
            break
        q = qn
    pri = np.log(np.array([(1 - q) ** 2, 2 * q * (1 - q), q ** 2]))[:, None]
    lp = ll + pri
    post = np.exp(lp - logsumexp(lp))
    return q, post


def em_two(ll):
    pi = 0.3
    for _ in range(100):
        pri = np.log(np.array([1 - pi, pi]))[:, None]
        lp = ll + pri
        post = np.exp(lp - logsumexp(lp))
        pin = float(np.clip(post[1].mean(), 0.001, 0.999))
        if abs(pin - pi) < 1e-9:
            pi = pin
            break
        pi = pin
    pri = np.log(np.array([1 - pi, pi]))[:, None]
    lp = ll + pri
    post = np.exp(lp - logsumexp(lp))
    return pi, post


# ---------------------------------------------------------------- library eligibility
n1 = {r["run"]: r for r in csv.DictReader(open(f"{N1}/library_table.tsv"), delimiter="\t")}
calls193066 = {r["run"]: r["n1_status"] for r in csv.DictReader(open(f"{N1}/gse193066_library_calls.tsv"), delimiter="\t")}
libs = []
for cohort in sorted(os.listdir(READS)):
    for lib in sorted(os.listdir(f"{READS}/{cohort}")):
        f = f"{READS}/{cohort}/{lib}/variants.tsv"
        if os.path.isfile(f):
            libs.append((cohort, lib, f))
elig = []
for cohort, lib, f in libs:
    ind = (n1.get(lib, {}).get("individual_id") or "").strip()
    reason = "ok"
    if not ind:
        reason = "no_individual_id"
    elif cohort == "GSE193066" and calls193066.get(lib) in ("remains_called_mixed", "untestable"):
        reason = "gse193066_" + calls193066.get(lib)
    elig.append({"library": lib, "cohort": cohort, "person_id": ind, "eligible": reason == "ok", "reason": reason,
                 "input_reads": n1.get(lib, {}).get("input_reads", "")})
with open(f"{out_dir}/library_eligibility.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(elig[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(elig)
elig_by_lib = {e["library"]: e for e in elig}

# ---------------------------------------------------------------- counts per library
lib_counts = {}   # (lib, rsid) -> (a, r)
with gzip.open(f"{out_dir}/library_variant_counts.tsv.gz", "wt") as lf:
    lf.write("library\tcohort\trsid\tgene\tref_reads\talt_reads\tother_reads\tdel_reads\n")
    for cohort, lib, f in libs:
        for row in csv.DictReader(open(f), delimiter="\t"):
            rs = row["rsid"]
            a = int(row["n_alt_plus"]) + int(row["n_alt_minus"])
            r = int(row["n_ref_plus"]) + int(row["n_ref_minus"])
            lib_counts[(lib, rs)] = (a, r)
            lf.write(f"{lib}\t{cohort}\t{rs}\t{row['gene']}\t{r}\t{a}\t{row['n_other']}\t{row['n_del']}\n")

# ---------------------------------------------------------------- person sums
person_libs = defaultdict(list)
person_cohort = {}
for e in elig:
    if e["eligible"]:
        person_libs[e["person_id"]].append(e["library"])
        person_cohort[e["person_id"]] = e["cohort"]
persons = sorted(person_libs)
with open(f"{out_dir}/person_variant_counts.tsv", "w") as fh:
    fh.write("person_id\tcohort\tn_libraries\trsid\tref_reads\talt_reads\tdepth\n")
    for p in persons:
        for rs in VARIANTS:
            a = sum(lib_counts[(l, rs)][0] for l in person_libs[p])
            r = sum(lib_counts[(l, rs)][1] for l in person_libs[p])
            fh.write(f"{p}\t{person_cohort[p]}\t{len(person_libs[p])}\t{rs}\t{r}\t{a}\t{a + r}\n")

# ---------------------------------------------------------------- calling
cohorts = sorted(set(person_cohort.values()))
calls = {p: {} for p in persons}           # person -> rs -> (class label, posterior of class)
callable_rows = []
freq_rows = []
geno_label3 = ["RR", "RA", "AA"]
geno_label2 = ["NC", "CARRIER"]
for rs in VARIANTS:
    for co in cohorts:
        pp = [p for p in persons if person_cohort[p] == co]
        libs_co = [e["library"] for e in elig if e["cohort"] == co]
        libs_ok = [e["library"] for e in elig if e["cohort"] == co and e["eligible"]]
        lib_depth = [sum(lib_counts[(l, rs)]) for l in libs_ok]
        a = np.array([sum(lib_counts[(l, rs)][0] for l in person_libs[p]) for p in pp], float)
        r = np.array([sum(lib_counts[(l, rs)][1] for l in person_libs[p]) for p in pp], float)
        d = a + r
        ok = d >= D_MIN
        row = {"rsid": rs, "gene": GENE[rs], "cohort": co, "libraries_total": len(libs_co), "libraries_eligible": len(libs_ok),
               "libraries_depth_ge8": int(sum(1 for x in lib_depth if x >= D_MIN)), "persons": len(pp),
               "persons_depth_ge8": int(ok.sum()), "median_depth_persons": float(np.median(d)) if len(pp) else float("nan"),
               "called": 0, "n_RR": 0, "n_RA": 0, "n_AA": 0, "n_NC": 0, "n_CARRIER": 0, "allele_freq_or_carrier_fraction": float("nan"),
               "hwe_chi2": float("nan"), "hwe_p": float("nan"), "note": ""}
        if ok.sum() < MIN_COHORT_CALLABLE:
            row["note"] = "fewer_than_20_persons_with_depth"
            callable_rows.append(row)
            continue
        idx = np.where(ok)[0]
        if rs in TWO_CLASS:
            ll = two_class_loglik(a[idx], d[idx])
            pi, post = em_two(ll)
            row["allele_freq_or_carrier_fraction"] = pi
            for j, i in enumerate(idx):
                if post[1, j] >= POST_MIN:
                    calls[pp[i]][rs] = ("CARRIER", float(post[1, j]))
                elif post[1, j] <= 1 - POST_MIN:
                    calls[pp[i]][rs] = ("NC", float(post[0, j]))
        else:
            ll = three_class_loglik(a[idx], d[idx])
            q, post = em_three(ll)
            row["allele_freq_or_carrier_fraction"] = q
            for j, i in enumerate(idx):
                k = int(np.argmax(post[:, j]))
                if post[k, j] >= POST_MIN:
                    calls[pp[i]][rs] = (geno_label3[k], float(post[k, j]))
        cnt = Counter(calls[p][rs][0] for p in pp if rs in calls[p])
        for k in ("RR", "RA", "AA", "NC", "CARRIER"):
            row["n_" + k] = cnt.get(k, 0)
        row["called"] = sum(cnt.values())
        if rs not in TWO_CLASS and row["called"] >= MIN_COHORT_CALLABLE:
            n = row["called"]
            qh = (cnt["RA"] + 2 * cnt["AA"]) / (2 * n)
            exp = np.array([(1 - qh) ** 2, 2 * qh * (1 - qh), qh ** 2]) * n
            obs = np.array([cnt["RR"], cnt["RA"], cnt["AA"]], float)
            m = exp > 0
            chi2 = float(np.sum((obs[m] - exp[m]) ** 2 / exp[m]))
            row["hwe_chi2"] = chi2
            row["hwe_p"] = math.erfc(math.sqrt(chi2 / 2))
        callable_rows.append(row)

fields = list(callable_rows[0].keys())
with open(f"{out_dir}/callable_fraction.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for row in callable_rows:
        row = dict(row)
        w.writerow(row)

# wide person table
with open(f"{out_dir}/person_genotypes.tsv", "w") as fh:
    cols = ["person_id", "cohort", "n_libraries", "libraries"]
    for rs in VARIANTS:
        cols += [f"{rs}_call", f"{rs}_post"]
    cols += ["pi_status"]
    fh.write("\t".join(cols) + "\n")
    for p in persons:
        z = calls[p].get("rs28929474")
        s = calls[p].get("rs17580")
        if z and s:
            if z[0] == "RR" and s[0] == "RR":
                pi = "MM"
            elif z[0] == "RA" and s[0] == "RR":
                pi = "MZ"
            else:
                pi = "other"
        else:
            pi = "NA"
        vals = [p, person_cohort[p], str(len(person_libs[p])), ",".join(person_libs[p])]
        for rs in VARIANTS:
            c = calls[p].get(rs)
            vals += [c[0] if c else "NA", f"{c[1]:.4f}" if c else "NA"]
        vals.append(pi)
        fh.write("\t".join(vals) + "\n")

# ---------------------------------------------------------------- library concordance within persons
conc_rows = []
multi = [p for p in persons if len(person_libs[p]) >= 2]
for rs in VARIANTS:
    n_pairs = n_conc = 0
    for p in multi:
        co = person_cohort[p]
        # per-library calls with the same likelihood and a flat-ish prior from the cohort estimate
        row = next(x for x in callable_rows if x["rsid"] == rs and x["cohort"] == co)
        fq = row["allele_freq_or_carrier_fraction"]
        if not (fq == fq):
            continue
        labs = []
        for l in person_libs[p]:
            a_, r_ = lib_counts[(l, rs)]
            d_ = a_ + r_
            if d_ < D_MIN:
                continue
            if rs in TWO_CLASS:
                ll = two_class_loglik([a_], [d_])
                pri = np.log(np.array([1 - fq, fq]))[:, None]
                lp = ll + pri
                po = np.exp(lp - logsumexp(lp))[:, 0]
                if po[1] >= POST_MIN:
                    labs.append("CARRIER")
                elif po[1] <= 1 - POST_MIN:
                    labs.append("NC")
            else:
                ll = three_class_loglik([a_], [d_])
                pri = np.log(np.array([(1 - fq) ** 2, 2 * fq * (1 - fq), fq ** 2]))[:, None]
                lp = ll + pri
                po = np.exp(lp - logsumexp(lp))[:, 0]
                k = int(np.argmax(po))
                if po[k] >= POST_MIN:
                    labs.append(geno_label3[k])
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                n_pairs += 1
                n_conc += int(labs[i] == labs[j])
    conc_rows.append({"rsid": rs, "gene": GENE[rs], "persons_with_2plus_libraries": len(multi), "library_pairs_both_called": n_pairs,
                      "concordant_pairs": n_conc, "concordance": (n_conc / n_pairs) if n_pairs else float("nan")})
with open(f"{out_dir}/library_concordance.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(conc_rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(conc_rows)

# summary across cohorts
summ = defaultdict(lambda: Counter())
for row in callable_rows:
    s = summ[row["rsid"]]
    s["persons"] += row["persons"]
    s["persons_depth_ge8"] += row["persons_depth_ge8"]
    s["called"] += row["called"]
    for k in ("RR", "RA", "AA", "NC", "CARRIER"):
        s["n_" + k] += row["n_" + k]
with open(f"{out_dir}/callable_summary_all_cohorts.tsv", "w") as fh:
    fh.write("rsid\tgene\tpersons\tpersons_depth_ge8\tcalled\tn_RR\tn_RA\tn_AA\tn_NC\tn_CARRIER\n")
    for rs in VARIANTS:
        s = summ[rs]
        fh.write(f"{rs}\t{GENE[rs]}\t{s['persons']}\t{s['persons_depth_ge8']}\t{s['called']}\t{s['n_RR']}\t{s['n_RA']}\t{s['n_AA']}\t{s['n_NC']}\t{s['n_CARRIER']}\n")

json.dump({"seed": SEED, "d_min": D_MIN, "error_rate": ERR, "beta_binomial": [BB_A, BB_B], "posterior_min": POST_MIN,
           "carrier_grid": CARRIER_GRID, "persons": len(persons), "libraries_eligible": sum(1 for e in elig if e["eligible"]),
           "libraries_total": len(elig), "numpy": np.__version__, "python": sys.version},
          open(f"{out_dir}/run_info.json", "w"), indent=1)
print("persons", len(persons), "eligible libraries", sum(1 for e in elig if e["eligible"]), "of", len(elig))
