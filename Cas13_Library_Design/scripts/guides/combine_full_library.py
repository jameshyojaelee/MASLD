#!/usr/bin/env python3
"""
combine_full_library.py — merge the three guide files into ONE synthesis-ready
master oligo list with a single `control_class` column and a harmonised schema.
Pure stdlib (csv) so it runs instantly with system python3 (no pandas/numpy).

Sources (Cas13_Library_Design/data/guides/):
  cas13_library_guides_vM38.csv   targets (positive controls flagged inline)
  cas13_control_guides_vM38.csv   essential-gene depletion controls
  cas13_nt_controls_vM38.csv      non-targeting controls

control_class: target | positive_control | essential | non_targeting
Output: Cas13_Library_Design/data/guides/library.csv  (one row per guide)
"""
import csv
from collections import Counter
from pathlib import Path

G = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/data/guides")
OUT = G / "library.csv"

COLS = [
    "guide_id", "control_class", "guide_seq",
    "gene_symbol_mouse", "gene_id_mouse", "gene_symbol_human", "biotype",
    "region", "target_seq", "tiger_score", "cas13_score", "combined_score",
    "n_isoforms_targeted", "tx_id_set", "tx_id_pos", "position",
    "rank_within_gene", "n_available_pool", "expected_direction",
    "tier", "library_arm", "has_human_evidence", "has_coloc",
    "mouse_hep_substrate", "mouse_untestable", "sc_disease_celltype",
    "is_positive_control", "pos_control_direction",
    "essentiality_chronos_liver", "n_liver_lines", "is_essential_liver",
]

def rows_from(path, classify):
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            out = {c: "" for c in COLS}
            for c in COLS:
                if c in r:
                    out[c] = r[c]
            classify(r, out)
            yield out

def cls_target(r, o):
    pc = str(r.get("is_positive_control", "")).strip() in ("True", "TRUE", "true", "1")
    o["control_class"] = "positive_control" if pc else "target"
    o["expected_direction"] = r.get("pos_control_direction", "") if pc else ""

def cls_essential(r, o):
    o["control_class"] = "essential"
    o["gene_symbol_human"] = r.get("human_symbol", "") or r.get("gene_symbol_human", "")
    o["expected_direction"] = r.get("expected_direction", "depletion") or "depletion"

def cls_nt(r, o):
    o["control_class"] = "non_targeting"
    o["expected_direction"] = "neutral"

# Targets first; then controls, SKIPPING any guide already present as a target.
# A few essential-QC genes (Pcna, Polr2l) are ALSO disease targets -> same guide_id;
# keep the richer target row (already flagged is_essential_liver) and don't double-list.
all_rows = []
seen = set()
dropped = 0
for src, fn in [("cas13_library_guides_vM38.csv", cls_target),
                ("cas13_control_guides_vM38.csv", cls_essential),
                ("cas13_nt_controls_vM38.csv", cls_nt)]:
    for r in rows_from(G / src, fn):
        if r["guide_id"] in seen:
            dropped += 1
            continue
        seen.add(r["guide_id"])
        all_rows.append(r)
print(f"(de-dup: dropped {dropped} control rows already present as targets)")

ids = [r["guide_id"] for r in all_rows]
assert len(ids) == len(set(ids)), "duplicate guide_id across files"
bad = [r["guide_id"] for r in all_rows if len(r["guide_seq"]) != 23]
assert not bad, f"non-23bp guides: {bad[:5]}"

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    w.writerows(all_rows)

print(f"wrote {len(all_rows):,} guides -> {OUT}")
for k, v in Counter(r["control_class"] for r in all_rows).most_common():
    print(f"  {k:16s} {v:>5,}")
print(f"unique guide_seq: {len(set(r['guide_seq'] for r in all_rows)):,} / {len(all_rows):,}")
