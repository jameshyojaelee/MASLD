#!/usr/bin/env python3
"""Check that Figure 7E qualification changes preserve every estimate and interval."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    old = read(a.project_root / "figures/main/fig7_model_framework/data/panel_deep_20261017/fig7E.tsv")
    new = read(a.output / "fig7/source_tables/fig7E.tsv")
    assert len(old) == len(new) == 20
    key = lambda r: (r["row_label"], r["mark"])
    old = {key(r): r for r in old}
    new = {key(r): r for r in new}
    assert old.keys() == new.keys()
    changed = []
    for k, before in old.items():
        after = new[k]
        for col in ("x", "lo", "hi", "unit", "n", "interval", "source", "field"):
            assert before[col] == after[col], (k, col)
        if before["block"] == "tested negative":
            assert after["block"] == "untestable"
            if after["mark"] != "bound":
                assert after["cls"] == "untestable" and after["glyph"] == "circle_open" and after["colour"] == "#9E9E9E"
            changed.append(k)
        else:
            assert before == after
    assert len(changed) == 3
    assert all(r["block"] != "tested negative" and r["cls"] != "tested negative" for r in new.values())
    result = {"status": "PASS", "rows": len(new), "estimates_intervals_units_sources_unchanged": True,
              "changed_marks": changed, "negative_calls_remaining": 0,
              "reason": "Treatment-response power and available sequence/routing bounds do not qualify a tested-negative conclusion. MDE remains a sensitivity reference, not an upper confidence limit."}
    (a.output / "fig7/checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
