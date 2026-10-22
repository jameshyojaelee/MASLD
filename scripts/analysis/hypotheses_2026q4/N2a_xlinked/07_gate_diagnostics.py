"""Why the gate counts are what they are: coverage versus heterozygous-call losses, for the primary set.
Depth and the pooled heterozygous rule only (the gate's own definitions); no X-linked per-biopsy ratio and no Delta.
Usage: python -I 07_gate_diagnostics.py RUN_ROOT OUT_TSV SET      (OUT_TSV must not exist)"""
import sys, os, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from n2a_lib import load_counts, site_record, MIN_DEPTH

root, out, SET = sys.argv[1:4]
if os.path.exists(out):
    sys.exit("refusing to overwrite " + out)
cls_of = {}
with open(os.path.join(root, "inputs", "sites.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cls_of[(r["chrom"], int(r["pos"]))] = r["cls"]
n_sites = {"X": sum(v.startswith("X_") for v in cls_of.values()), "auto": sum(v == "AUTO" for v in cls_of.values())}
pairs = [r for r in csv.DictReader(open(os.path.join(root, "inputs", "n2a_pairs.tsv")), delimiter="\t") if r["include"] == "True" and r["set"] == SET]
pairs.sort(key=lambda r: (r["cohort"], r["person_id"]))
cols = ["cohort", "person_id", "X_sites_listed", "X_with_reads_b1", "X_with_reads_b2", "X_depth_ge20_b1", "X_depth_ge20_b2", "X_depth_ge20_both", "X_informative",
        "X_informative_over_both", "auto_sites_listed", "auto_depth_ge20_both", "auto_informative", "auto_informative_over_both"]
with open(out, "x") as fh:
    fh.write("\t".join(cols) + "\n")
    for pr in pairs:
        c = [load_counts(os.path.join(root, "counts", pr[k] + ".counts.tsv.gz")) for k in ("lib1", "lib2")]
        depth = lambda a: a[0] + a[1] + a[2] + a[3]
        row = dict(cohort=pr["cohort"], person_id=pr["person_id"], X_sites_listed=n_sites["X"], auto_sites_listed=n_sites["auto"])
        for tag, test in (("X", lambda k: cls_of[k].startswith("X_")), ("auto", lambda k: cls_of[k] == "AUTO")):
            keys = [k for k in cls_of if test(k)]
            both = [k for k in keys if k in c[0] and k in c[1] and depth(c[0][k]) >= MIN_DEPTH and depth(c[1][k]) >= MIN_DEPTH]
            inf = [k for k in both if site_record(c[0][k], c[1][k]) is not None]
            row[tag + "_depth_ge20_both"] = len(both)
            row[tag + ("_informative")] = len(inf)
            row[tag + "_informative_over_both"] = round(len(inf) / len(both), 3) if both else float("nan")
            if tag == "X":
                row["X_with_reads_b1"] = sum(1 for k in keys if k in c[0]); row["X_with_reads_b2"] = sum(1 for k in keys if k in c[1])
                row["X_depth_ge20_b1"] = sum(1 for k in keys if k in c[0] and depth(c[0][k]) >= MIN_DEPTH)
                row["X_depth_ge20_b2"] = sum(1 for k in keys if k in c[1] and depth(c[1][k]) >= MIN_DEPTH)
        fh.write("\t".join(str(row[k]) for k in cols) + "\n")
print(open(out).read())
