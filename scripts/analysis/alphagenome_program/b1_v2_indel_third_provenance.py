#!/usr/bin/env python3
"""Correction 9, run separately when the third dbSNP orientation deposit lands after the layer was built.

b1_response_layer_v2.py reads the rerun deposit if it is already on disk with a non-empty tables/ and no
SUPERSEDED marker. When the rerun is still in flight at layer-build time, this script adds the third
provenance afterwards WITHOUT rebuilding or overwriting the layer: it reads the deposited
indel_orientation_provenance.tsv, joins the rerun's calls, and writes two new files.

It refuses to run if the rerun deposit is missing, empty, or SUPERSEDED, and it never edits anything under
GWAS/finemapping/results/alphagenome_atlas/.

Usage: b1_v2_indel_third_provenance.py <v2 deposit directory> [<rerun deposit directory>]
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS = PROJECT / "GWAS/finemapping/results/alphagenome_atlas"
DEFAULT_RERUN = ATLAS / "p0-dbsnp-orientation-20260914T192341Z"


def read_tsv(p):
    with open(p) as h:
        return list(csv.DictReader(h, delimiter="\t"))


def variant_class(ref: str, alt: str) -> str:
    r, a = str(ref).upper(), str(alt).upper()
    if len(r) == len(a):
        return "snv" if len(r) == 1 else "mnv"
    return "insertion" if len(a) > len(r) else "deletion"


def main() -> None:
    dep = pathlib.Path(sys.argv[1]).resolve()
    rerun = pathlib.Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else DEFAULT_RERUN
    tsv = rerun / "tables/dbsnp_indel_orientation.tsv"
    if (rerun / "SUPERSEDED.txt").exists():
        raise SystemExit(f"{rerun.name} carries SUPERSEDED.txt; not adding it as a provenance")
    if not tsv.exists() or tsv.stat().st_size == 0:
        raise SystemExit(f"{tsv} is absent or empty; the rerun has not landed")

    calls = {r["source_variant_id"]: r["resolved_variant_uid"]
             for r in read_tsv(tsv) if r["resolved_variant_uid"]}
    rows = read_tsv(dep / "indel_orientation_provenance.tsv")
    out = []
    for r in rows:
        uid = r["variant_uid"] or calls.get(r["source_variant_id"], "") or r["hg38_uid_recovered"]
        cls = variant_class(*uid.split(":")[2:4]) if uid else ""
        out.append({**r, "uid_dbsnp_rerun": calls.get(r["source_variant_id"], ""), "class_rerun": cls})

    dest = dep / "indel_orientation_provenance_rerun.tsv"
    if dest.exists():
        raise SystemExit(f"refusing to overwrite {dest}")
    with dest.open("w") as h:
        w = csv.DictWriter(h, fieldnames=list(out[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    def counts(k):
        return {v: sum(1 for r in out if r[k] == v) for v in sorted({r[k] for r in out})}

    def changed(a, b):
        return sum(1 for r in out if r[a] and r[b] and r[a] != r[b])
    summary = {
        "rerun_deposit": rerun.name, "n_rerun_resolved_variants": len(calls),
        "n_indel_rows": len(out),
        "n_indel_rows_the_rerun_resolves": sum(1 for r in out if r["uid_dbsnp_rerun"]),
        "class_counts": {"union_deposited": counts("class_union_deposited"),
                         "later_alone": counts("class_later_alone"),
                         "earlier_alone": counts("class_earlier_alone"),
                         "rerun": counts("class_rerun")},
        "n_rows_changing_class": {
            "rerun_vs_union_deposited": changed("class_rerun", "class_union_deposited"),
            "rerun_vs_later_alone": changed("class_rerun", "class_later_alone"),
            "rerun_vs_earlier_alone": changed("class_rerun", "class_earlier_alone")},
        "note": "reported only; the deposited variant_class column is unchanged",
    }
    json.dump(summary, (dep / "indel_orientation_rerun_summary.json").open("w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
