#!/usr/bin/env python3
"""Which numbers moved between the wave-1 B1 deposit and this correction, and which did not.

The corrections are to labels, denominators, units and one recoverable key. The analysis is unchanged, so
almost every number must be identical; anything that moved has to be named. This flattens both deposits'
JSON summaries, compares every shared path, and writes one row per path.

A path that exists in only one deposit is reported as added or removed, because a renamed column is a moved
label even when the value behind it is the same number.

Usage: b1_v2_compare_to_wave1.py <wave1 directory> <v2 directory> <output tsv>
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

FILES = ["response_layer_summary.json", "b1_checks.json", "direction_track_panel_summary.json",
         "mpra_rebuild_summary.json"]


def flatten(obj, prefix=""):
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.update(flatten(v, f"{prefix}/{k}"))
    elif isinstance(obj, list):
        if all(not isinstance(x, (dict, list)) for x in obj):
            out[prefix] = ";".join(str(x) for x in obj)
        else:
            for i, v in enumerate(obj):
                out.update(flatten(v, f"{prefix}[{i}]"))
    else:
        out[prefix] = obj
    return out


def load(d: pathlib.Path) -> dict:
    out = {}
    for f in FILES:
        p = d / f
        if p.exists():
            out.update(flatten(json.load(p.open()), f"/{f}"))
    return out


def main() -> None:
    w1, v2, dest = (pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]))
    a, b = load(w1), load(v2)
    rows = []
    for k in sorted(set(a) | set(b)):
        if k in a and k in b:
            same = str(a[k]) == str(b[k])
            rows.append({"path": k, "state": ("identical" if same else "MOVED"),
                         "wave1": a[k], "v2": b[k]})
        elif k in b:
            rows.append({"path": k, "state": "added_in_v2", "wave1": "", "v2": b[k]})
        else:
            rows.append({"path": k, "state": "absent_from_v2", "wave1": a[k], "v2": ""})
    with dest.open("w") as h:
        w = csv.DictWriter(h, fieldnames=["path", "state", "wave1", "v2"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    counts = {s: sum(1 for r in rows if r["state"] == s) for s in
              ("identical", "MOVED", "added_in_v2", "absent_from_v2")}
    print(json.dumps(counts, indent=1))
    for r in rows:
        if r["state"] == "MOVED":
            print(f"MOVED {r['path']}\n   wave1 {r['wave1']}\n   v2    {r['v2']}")


if __name__ == "__main__":
    main()
