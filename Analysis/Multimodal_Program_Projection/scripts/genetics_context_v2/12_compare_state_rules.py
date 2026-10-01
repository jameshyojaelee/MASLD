#!/usr/bin/env python3
"""Compare frozen evidence classes built under the canonical and TREAT state rules.

Both tables come from 01_freeze_and_rederive.py on the same inputs, one with
--state-rule canonical (primary) and one with --state-rule treat (sensitivity).
Writes, into a new directory:
  class_transition_counts.tsv   canonical class x TREAT class, with gene counts
  class_differences.tsv         one row per gene whose class differs
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", type=Path, required=True)
    parser.add_argument("--treat", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def read_classes(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_symbol = {row["gene_symbol"]: row for row in rows}
    if len(by_symbol) != len(rows):
        raise SystemExit(f"duplicate gene_symbol in {path}")
    return by_symbol


def main() -> None:
    args = parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"refusing to write into existing {args.out_dir}")
    canonical = read_classes(args.canonical)
    treat = read_classes(args.treat)
    if set(canonical) != set(treat):
        raise SystemExit("the two rules must classify the same gene universe")

    transitions = Counter((canonical[s]["static_class"], treat[s]["static_class"]) for s in canonical)
    differences = []
    for symbol in sorted(canonical):
        c, t = canonical[symbol], treat[symbol]
        if c["static_class"] == t["static_class"]:
            continue
        differences.append(
            {
                "gene_symbol": symbol,
                "class_canonical": c["static_class"],
                "class_treat": t["static_class"],
                "state_canonical": c["established_state_associated"],
                "state_treat": t["established_state_associated"],
                "primary_genetic": c["primary_genetic"],
                "bulk_logFC": c["bulk_logFC"],
                "bulk_padj": c.get("bulk_padj", ""),
                "bulk_treat_fdr": c["bulk_treat_fdr"],
            }
        )

    args.out_dir.mkdir(parents=True)
    with (args.out_dir / "class_transition_counts.tsv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["class_canonical", "class_treat", "n_genes"])
        for (c, t), n in sorted(transitions.items()):
            writer.writerow([c, t, n])
    fields = list(differences[0]) if differences else ["gene_symbol", "class_canonical", "class_treat"]
    with (args.out_dir / "class_differences.tsv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(differences)
    changed = Counter((d["class_canonical"], d["class_treat"]) for d in differences)
    print(f"PASS: genes={len(canonical)}; class_changes={len(differences)}")
    for (c, t), n in sorted(changed.items()):
        print(f"  canonical={c}\ttreat={t}\tn={n}")


if __name__ == "__main__":
    main()
