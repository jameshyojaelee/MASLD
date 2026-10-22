"""N2a power gate (PRESPEC_N2a.md). Uses depth, heterozygous calls and AUTOSOMAL cross-fit deltas only.
It never computes or writes an X-linked allele ratio or an X-linked delta.
Usage: python -I 05_gate.py RUN_ROOT OUT_DIR [SETS]   (OUT_DIR must not exist; SETS e.g. primary or primary,secondary)
Writes OUT_DIR/gate_per_woman.tsv, gate_summary.json, gate_decision.txt."""
import sys, os, csv, json, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from n2a_lib import load_counts, site_record

root, outdir = sys.argv[1], sys.argv[2]
if os.path.exists(outdir):
    sys.exit("refusing to overwrite " + outdir)
os.makedirs(outdir)
SEED, REPS = 20260923, 1000
SETS = sys.argv[3].split(",") if len(sys.argv) > 3 else ["primary", "secondary"]
G1_WOMEN, G1_SITES, G2_MDC, Z = 15, 5, 0.20, 2.80

cls_of = {}
with open(os.path.join(root, "inputs", "sites.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cls_of[(r["chrom"], int(r["pos"]))] = r["cls"]


pairs = [r for r in csv.DictReader(open(os.path.join(root, "inputs", "n2a_pairs.tsv")), delimiter="\t") if r["include"] == "True" and r["set"] in SETS]
pairs.sort(key=lambda r: (r["set"], r["cohort"], r["person_id"]))
rng = np.random.default_rng(SEED)
rows = []
for pr in pairs:
    c1 = load_counts(os.path.join(root, "counts", pr["lib1"] + ".counts.tsv.gz"))
    c2 = load_counts(os.path.join(root, "counts", pr["lib2"] + ".counts.tsv.gz"))
    xs, au = [], []
    for key in sorted(c1.keys() & c2.keys()):   # sorted: set order varies with the string-hash seed and would change the seeded draws
        rec = site_record(c1[key], c2[key])
        if rec is not None:
            (xs if cls_of[key].startswith("X_") else au).append(rec)
    kx, ka = len(xs), len(au)
    row = dict(set=pr["set"], cohort=pr["cohort"], person_id=pr["person_id"], lib1=pr["lib1"], lib2=pr["lib2"],
               k_x_informative=kx, k_auto_informative=ka, sd_auto_null=float("nan"), mean_auto_null=float("nan"), mdc=float("nan"),
               auto_pool_shortfall_bins=0)
    if kx >= 1 and ka >= 1:
        xb = collections.Counter(r[0] for r in xs)
        A = np.array(au, dtype=float)
        sx1 = np.zeros(REPS); sn1 = np.zeros(REPS); sx2 = np.zeros(REPS); sn2 = np.zeros(REPS)
        for b, m in xb.items():
            pool = A[A[:, 0] == b]
            if len(pool) == 0:
                pool = A
                row["auto_pool_shortfall_bins"] += 1
            if len(pool) >= m:
                idx = np.argsort(rng.random((REPS, len(pool))), axis=1)[:, :m]
            else:
                idx = rng.integers(0, len(pool), size=(REPS, m))
                row["auto_pool_shortfall_bins"] += 1
            sx1 += pool[idx, 1].sum(1); sn1 += pool[idx, 2].sum(1); sx2 += pool[idx, 3].sum(1); sn2 += pool[idx, 4].sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            delta = sx2 / sn2 - sx1 / sn1
        delta = delta[np.isfinite(delta)]
        row["sd_auto_null"] = float(np.std(delta, ddof=1)); row["mean_auto_null"] = float(np.mean(delta))
        row["mdc"] = Z * row["sd_auto_null"]
    rows.append(row)

cols = ["set", "cohort", "person_id", "lib1", "lib2", "k_x_informative", "k_auto_informative", "sd_auto_null", "mean_auto_null", "mdc", "auto_pool_shortfall_bins"]
with open(os.path.join(outdir, "gate_per_woman.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

summ = {}
for st in SETS:
    rr = [r for r in rows if r["set"] == st]
    ok = [r for r in rr if r["k_x_informative"] >= G1_SITES]
    mdcs = [r["mdc"] for r in ok if np.isfinite(r["mdc"])]
    sds = [r["sd_auto_null"] for r in ok if np.isfinite(r["sd_auto_null"])]
    s = dict(n_women_in_set=len(rr), n_women_ge5_sites=len(ok), n_women_ge3_sites=sum(r["k_x_informative"] >= 3 for r in rr),
             n_women_ge10_sites=sum(r["k_x_informative"] >= 10 for r in rr),
             k_x_informative_by_woman=sorted(r["k_x_informative"] for r in rr),
             G1_pass=len(ok) >= G1_WOMEN,
             median_mdc=float(np.median(mdcs)) if mdcs else None, max_mdc=float(np.max(mdcs)) if mdcs else None,
             cohort_mean_mdc=float(Z * np.sqrt(np.sum(np.square(sds))) / len(sds)) if sds else None,
             mean_of_autosomal_null_means=float(np.mean([r["mean_auto_null"] for r in ok])) if ok else None)
    s["G2_pass"] = (s["median_mdc"] is not None) and s["median_mdc"] <= G2_MDC
    s["gate_pass"] = bool(s["G1_pass"] and s["G2_pass"])
    summ[st] = s
summ["parameters"] = dict(seed=SEED, reps=REPS, g1_women=G1_WOMEN, g1_sites=G1_SITES, g2_mdc=G2_MDC, z=Z)
json.dump(summ, open(os.path.join(outdir, "gate_summary.json"), "w"), indent=1)
with open(os.path.join(outdir, "gate_decision.txt"), "w") as fh:
    for st in SETS:
        s = summ[st]
        fh.write("%s: women %d | with >=5 informative X sites %d (G1 needs %d) | median MDC %s (G2 needs <= %.2f) | gate %s\n" % (
            st, s["n_women_in_set"], s["n_women_ge5_sites"], G1_WOMEN, "NA" if s["median_mdc"] is None else "%.3f" % s["median_mdc"], G2_MDC,
            "PASS" if s["gate_pass"] else "FAIL (untestable)"))
print(open(os.path.join(outdir, "gate_decision.txt")).read())
print(json.dumps({k: v for k, v in summ["primary"].items()}, indent=1))
