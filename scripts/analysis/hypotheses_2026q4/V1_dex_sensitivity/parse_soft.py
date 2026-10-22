"""V1: parse GSE320593 SOFT metadata into a per-sample design table.

Usage: python3 -I parse_soft.py <soft.gz> <out.tsv>

Treats the SOFT file as data only. Donor and procedure are taken from the
sample title, <donor>_<tissue>_<non-surgical|surgical>_<DEX|NO_DEX>, and
cross-checked against the `tissue` and `treatment` characteristics.
"""
import csv
import gzip
import re
import sys
from collections import defaultdict

soft_path, out_path = sys.argv[1], sys.argv[2]

samples = []
cur = None
with gzip.open(soft_path, "rt", encoding="utf-8") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            cur = {"gsm": line.split(" = ", 1)[1]}
            samples.append(cur)
        elif cur is None:
            continue
        elif line.startswith("!Sample_title = "):
            cur["title"] = line.split(" = ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = "):
            key, val = line.split(" = ", 1)[1].split(": ", 1)
            cur[key.strip()] = val.strip()
        elif line.startswith("!Sample_description = Library name: "):
            cur["library"] = line.split("Library name: ", 1)[1].strip()

title_re = re.compile(
    r"^(?P<donor>[A-Z]+\d+)_(?P<tissue>liver_tissue|subcutaneous_adipose_tissue)_"
    r"(?P<proc>non-surgical|surgical)_(?P<arm>DEX|NO_DEX)$"
)

rows = []
for s in samples:
    m = title_re.match(s["title"])
    if m is None:
        sys.exit(f"title does not match expected pattern: {s['title']!r}")
    if m["tissue"] != s["tissue"]:
        sys.exit(f"tissue mismatch for {s['gsm']}")
    if m["arm"] != s["treatment"]:
        sys.exit(f"treatment mismatch for {s['gsm']}")
    rows.append({
        "gsm": s["gsm"],
        "library": s["library"],
        "title": s["title"],
        "donor": m["donor"],
        "tissue": "liver" if m["tissue"] == "liver_tissue" else "adipose",
        "procedure": "pre" if m["proc"] == "non-surgical" else "surg",
        "arm": "DEX" if m["arm"] == "DEX" else "NODEX",
        "batch": s["batch"],
    })

# Design checks: one pre and one surg biopsy per donor per tissue, arm constant.
by_donor = defaultdict(list)
for r in rows:
    by_donor[(r["tissue"], r["donor"])].append(r)
for key, rs in by_donor.items():
    procs = sorted(r["procedure"] for r in rs)
    arms = {r["arm"] for r in rs}
    if procs != ["pre", "surg"] or len(arms) != 1:
        sys.exit(f"unexpected design for {key}: {procs} {arms}")
if len({r["library"] for r in rows}) != len(rows):
    sys.exit("library names not unique")

with open(out_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"wrote {len(rows)} samples, {len(by_donor)} donor-tissue pairs to {out_path}")
