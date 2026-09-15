#!/usr/bin/env python3
"""Row-level proof that the v2 corrections changed labels, not values.

The summary-path comparison (b1_v2_compare_to_wave1.py) shows that the deposits' JSON summaries agree. This
goes further and compares the two variant tables cell by cell on the (signal_uid, variant_key) key: every
column the two deposits share, plus each renamed column against the wave-1 name it replaced. Any non-zero
count here would mean a correction moved a value, which it must not.

Usage: b1_v2_row_level_equality.py <wave1 directory> <v2 directory>
"""

from __future__ import annotations

import json
import pathlib
import sys

RENAMES = {"atac_liver_quantile": "atac_primary_liver_quantile",
           "dnase_liver_quantile": "dnase_primary_liver_quantile",
           "h3k27ac_liver_quantile": "h3k27ac_primary_liver_quantile",
           "rna_target_gene_quantile": "rna_target_gene_primary_liver_quantile",
           "rna_strongest_abs_quantile": "rna_strongest_abs_primary_liver_quantile",
           "splice_site_usage_liver_quantile": "splice_site_usage_primary_liver_quantile"}


def main() -> None:
    import pandas as pd
    w1, v2 = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    key = ["signal_uid", "variant_key"]
    a = pd.read_csv(w1 / "response_layer_variants.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    b = pd.read_csv(v2 / "response_layer_variants.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    if len(a) != len(b):
        raise SystemExit(f"row counts differ: {len(a)} vs {len(b)}")
    a, b = a.set_index(key), b.set_index(key)
    if set(a.index) != set(b.index):
        raise SystemExit("the two deposits do not carry the same (signal_uid, variant_key) rows")
    a = a.loc[b.index]
    shared = [c for c in a.columns if c in b.columns]
    diff_shared = {c: int((a[c].values != b[c].values).sum()) for c in shared}
    diff_renamed = {f"{k} -> {v}": int((a[k].values != b[v].values).sum()) for k, v in RENAMES.items()}
    out = {
        "n_rows": int(len(b)), "key": key,
        "n_columns_wave1": int(a.shape[1] + 2), "n_columns_v2": int(b.shape[1] + 2),
        "n_shared_columns": len(shared),
        "n_shared_columns_with_any_differing_cell": sum(1 for v in diff_shared.values() if v),
        "differing_cells_per_shared_column": {c: v for c, v in diff_shared.items() if v},
        "renamed_columns_differing_cells": diff_renamed,
        "columns_added_in_v2": sorted(set(b.columns) - set(a.columns) - set(RENAMES.values())),
        "columns_renamed_in_v2": RENAMES,
        "columns_absent_from_v2_and_not_renamed": sorted(set(a.columns) - set(b.columns) - set(RENAMES)),
        "verdict": ("every shared and renamed column is identical cell for cell"
                    if not any(diff_shared.values()) and not any(diff_renamed.values())
                    else "AT LEAST ONE VALUE MOVED"),
    }
    dest = v2 / "row_level_equality_vs_wave1.json"
    if dest.exists():
        raise SystemExit(f"refusing to overwrite {dest}")
    json.dump(out, dest.open("w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
