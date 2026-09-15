#!/usr/bin/env python3
"""B2 (spec section 4): rebuild the per-context MPRA differential-allelic-variant calls from the source
reproduction, and quantify what the p3a deposit's sheet/context swap does to the numbers it printed.

The p3a deposit is read only. Nothing under GWAS/finemapping/results/alphagenome_atlas/ is modified.

Outputs into the response-layer directory:
  mpra_dav_calls_rebuilt.tsv        one row per (context, element) official DAV, context taken from the
                                    source reproduction that produced it
  mpra_swap_consequences.tsv        per p3a-printed quantity: the label it carries, the label set it was
                                    actually scored against, and the corrected positive count

Usage: b1_mpra_rebuild.py <response layer directory>
"""

from __future__ import annotations

import csv
import gzip
import json
import pathlib
import sys

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MPRA_SRC = PROJECT / "GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction"
P3A = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z/tables"
CONTEXTS = [("HepG2", "control", "HepG2_control", "dav_HepG2_ctrl"),
            ("HepG2", "PAOA", "HepG2_PAOA", "dav_HepG2_PAOA"),
            ("LX2", "control", "LX2_control", "dav_LX2_ctrl"),
            ("LX2", "TGFb", "LX2_TGFb", "dav_LX2_TGFB")]
# Each cell line's control is the constitutive baseline for its stimulated context.
STIMULUS_PAIRS = [("HepG2_PAOA", "HepG2_control"), ("LX2_TGFb", "LX2_control")]


def read_tsv(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv_once(path, rows, cols):
    path = pathlib.Path(path)
    if path.exists():
        raise RuntimeError(f"refusing to overwrite: {path}")
    with open(path, "w") as h:
        w = csv.DictWriter(h, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    outdir = pathlib.Path(sys.argv[1]).resolve()
    lib = read_tsv(P3A / "mpra_library_intervals.tsv")
    lib_ids = {r["interval"] for r in lib}

    official, rebuilt = {}, []
    for cl, ctx, name, p3a_col in CONTEXTS:
        rows = read_tsv(MPRA_SRC / f"{cl}.{ctx}.official_comparison.tsv.gz")
        official[name] = {r["element_id"] for r in rows}
        for r in rows:
            rebuilt.append({
                "context": name, "cell_line": cl, "condition": ctx, "element_id": r["element_id"],
                "rsid": r["rsid"], "official_log2FC_alt_over_ref": r["official_log2FC"],
                "official_fdr": r["official_fdr"],
                "reproduced_log2FC_alt_over_ref": r["log2FC_alt_over_ref"], "reproduced_fdr": r["fdr"],
                "reproduced_dav": r["reproduced_dav"], "direction_match": r["direction_match"],
                "in_p3a_scored_library": str(r["element_id"] in lib_ids),
                "source": f"{MPRA_SRC.name}/{cl}.{ctx}.official_comparison.tsv.gz"})
    rebuilt.sort(key=lambda r: (r["context"], r["element_id"]))
    write_tsv_once(outdir / "mpra_dav_calls_rebuilt.tsv", rebuilt, list(rebuilt[0].keys()))

    # what the p3a columns actually hold
    p3a_true = {}
    for _, _, _, col in CONTEXTS:
        p3a_true[col] = {r["interval"] for r in lib if r[col] == "True"}
    in_lib = {name: (s & lib_ids) for name, s in official.items()}

    cons = []
    for cl, ctx, name, col in CONTEXTS:
        matched = [n for n, s in in_lib.items() if s == p3a_true[col]]
        cons.append({"p3a_quantity": f"mpra_library_intervals.{col}",
                     "p3a_label": name, "p3a_n_positive": len(p3a_true[col]),
                     "label_set_actually_used": (matched[0] if matched else "no exact match"),
                     "correct_n_positive": len(in_lib[name]),
                     "label_is_correct": str(matched == [name])})
    for stim, base in STIMULUS_PAIRS:
        stim_col = dict((n, c) for _, _, n, c in CONTEXTS)[stim]
        base_col = dict((n, c) for _, _, n, c in CONTEXTS)[base]
        printed = len(p3a_true[stim_col] - p3a_true[base_col])
        correct = len(in_lib[stim] - in_lib[base])
        cons.append({"p3a_quantity": f"external_benchmark_mpra.{stim}|stimulus_specific_vs_constitutive n_pos",
                     "p3a_label": f"{stim} specific vs {base}", "p3a_n_positive": printed,
                     "label_set_actually_used": f"{stim_col} minus {base_col} as deposited",
                     "correct_n_positive": correct,
                     "label_is_correct": str(printed == correct)})
    write_tsv_once(outdir / "mpra_swap_consequences.tsv", cons, list(cons[0].keys()))

    summary = {
        "n_rebuilt_rows": len(rebuilt),
        "official_calls_by_context": {n: len(s) for n, s in official.items()},
        "official_calls_inside_p3a_scored_library": {n: len(s) for n, s in in_lib.items()},
        "p3a_positive_counts": {c: len(s) for c, s in p3a_true.items()},
        "consequences": cons,
        "p3a_deposit_read_only": True,
    }
    json.dump(summary, (outdir / "mpra_rebuild_summary.json").open("w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
