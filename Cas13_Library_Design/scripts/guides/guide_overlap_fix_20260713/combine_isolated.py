#!/usr/bin/env python3
"""Combine isolated target/control outputs with existing non-targeting controls."""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

RUN = Path(__file__).resolve().parent
G = RUN / "output"
NT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
          "Cas13_Library_Design/data/guides/cas13_nt_controls_vM38.csv")
OUT = G / "library_sequence_nonoverlap.csv"

COLS = [
    "guide_id", "control_class", "guide_seq", "gene_symbol_mouse", "gene_id_mouse",
    "gene_symbol_human", "biotype", "region", "target_seq", "tiger_score",
    "cas13_score", "combined_score", "n_isoforms_targeted", "tx_id_set", "tx_id_pos",
    "position", "rank_within_gene", "n_available_pool", "expected_direction", "tier",
    "library_arm", "has_human_evidence", "has_coloc", "mouse_hep_substrate",
    "mouse_untestable", "sc_disease_celltype", "is_positive_control",
    "pos_control_direction", "essentiality_chronos_liver", "n_liver_lines",
    "is_essential_liver",
]


def rows_from(path: Path, kind: str):
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            out = {column: row.get(column, "") for column in COLS}
            if kind == "target":
                positive = row.get("is_positive_control", "").lower() in {"true", "1"}
                out["control_class"] = "positive_control" if positive else "target"
                out["expected_direction"] = row.get("pos_control_direction", "") if positive else ""
            elif kind == "essential":
                out["control_class"] = "essential"
                out["gene_symbol_human"] = row.get("human_symbol", "")
                out["expected_direction"] = "depletion"
            else:
                out["control_class"] = "non_targeting"
                out["expected_direction"] = "neutral"
            yield out


def main() -> None:
    sources = [
        (G / "cas13_library_guides_vM38.csv", "target"),
        (G / "cas13_control_guides_vM38.csv", "essential"),
        (NT, "non_targeting"),
    ]
    rows = []
    seen_ids = set()
    seen_sequences = set()
    dropped_id = 0
    for path, kind in sources:
        for row in rows_from(path, kind):
            if row["guide_id"] in seen_ids:
                dropped_id += 1
                continue
            if row["guide_seq"] in seen_sequences:
                raise RuntimeError(f"duplicate guide sequence with a different ID: {row['guide_id']}")
            if len(row["guide_seq"]) != 23:
                raise RuntimeError(f"non-23-nt guide: {row['guide_id']}")
            seen_ids.add(row["guide_id"])
            seen_sequences.add(row["guide_seq"])
            rows.append(row)
    with OUT.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[combine] wrote {len(rows):,} unique guides -> {OUT}")
    print(f"[combine] dropped {dropped_id} control IDs already present among targets")
    print(f"[combine] classes: {dict(Counter(r['control_class'] for r in rows))}")


if __name__ == "__main__":
    main()
