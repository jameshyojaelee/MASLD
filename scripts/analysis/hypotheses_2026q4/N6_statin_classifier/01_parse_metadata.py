#!/usr/bin/env python3
"""N6: parse GSE130991 sample characteristics into one row per array and report gate G1.

Usage: python3 -I 01_parse_metadata.py <series_matrix.txt.gz> <out_dir>
Reads only the header block (stops at the data table).
"""
import collections
import csv
import gzip
import os
import re
import sys


def main(matrix, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    gsm, title = [], []
    char_rows = []
    with gzip.open(matrix, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            f = line.rstrip("\n").split("\t")
            if f[0] == "!Sample_geo_accession":
                gsm = [x.strip('"') for x in f[1:]]
            elif f[0] == "!Sample_title":
                title = [x.strip('"') for x in f[1:]]
            elif f[0] == "!Sample_characteristics_ch1":
                char_rows.append([x.strip('"') for x in f[1:]])
    n = len(gsm)
    assert n == len(title) and all(len(r) == n for r in char_rows), "ragged header"
    recs = [collections.OrderedDict(gsm=g, title=t) for g, t in zip(gsm, title)]
    dup_keys = collections.Counter()
    for r in char_rows:
        for j, cell in enumerate(r):
            if not cell:
                continue
            k, _, v = cell.partition(": ")
            k = re.sub(r"\s+", "_", k.strip().lower())
            if k in recs[j]:
                dup_keys[k] += 1
            recs[j][k] = v.strip()
    keys = ["gsm", "title"]
    for r in recs:
        for k in r:
            if k not in keys:
                keys.append(k)
    m = re.compile(r"^(\d+)[_-]")
    pats = [("single", r"\d+_\(HTA-2_0\)"), ("dash", r"\d+-\d+_\(HTA-2_0\)"), ("double", r"\d+_\d+_\(HTA-2_0\)")]
    for r in recs:
        mm = m.match(r["title"])
        r["array_index"] = mm.group(1) if mm else "NA"
        r["title_pattern"] = next((n for n, pt in pats if re.fullmatch(pt, r["title"])), "other")
    keys.extend(["array_index", "title_pattern"])
    with open(os.path.join(out_dir, "samples.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="NA")
        w.writeheader()
        w.writerows(recs)

    lines = []
    lines.append(f"arrays\t{n}")
    lines.append(f"unique_gsm\t{len(set(gsm))}")
    lines.append(f"unique_titles\t{len(set(title))}")
    lines.append(f"duplicate_key_collisions\t{dict(dup_keys)}")
    for k in ["tissue", "statin_treatment", "statin", "matched_cohort_statin_treatment", "gender", "diabetes_treatment"]:
        c = collections.Counter(r.get(k, "NA") for r in recs)
        lines.append(f"{k}\t{dict(c)}")
    present = {k: sum(1 for r in recs if k in r) for k in keys if k not in ("gsm", "title")}
    lines.append("fields_nonmissing\t" + str(present))
    with open(os.path.join(out_dir, "gate_G1_metadata_summary.tsv"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
