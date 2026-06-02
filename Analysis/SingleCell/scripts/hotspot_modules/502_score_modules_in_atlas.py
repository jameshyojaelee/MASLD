"""Aggregate per-cell-type module scores into a single canonical table.

Reads cell_scores.parquet from each of the 7 CT runs and writes:
  Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet
  Analysis/SingleCell/results_gpu_v2/hotspot_modules/donor_scores_all.tsv

Module IDs are namespaced as <cell_type>__<module_int> to keep them unique.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import RESULTS, RUN_ORDER


def main() -> None:
    cell_frames, donor_frames, skipped = [], [], []
    for ct in RUN_ORDER:
        cp = RESULTS / ct / "cell_scores.parquet"
        dp = RESULTS / ct / "donor_scores.tsv"
        if not cp.exists() or not dp.exists():
            print(f"[WARN] missing outputs for {ct}; skipping "
                  f"(cell_scores.parquet={cp.exists()}, donor_scores.tsv={dp.exists()})")
            skipped.append(ct)
            continue
        cf = pd.read_parquet(cp)
        cf["cell_type"] = ct
        cf["module"] = ct + "__" + cf["module"].astype(str)
        cell_frames.append(cf)

        df = pd.read_csv(dp, sep="\t")
        df["cell_type"] = ct
        df["module"] = ct + "__" + df["module"].astype(str)
        donor_frames.append(df)

    if not cell_frames:
        raise RuntimeError("No cell_scores.parquet found in any CT directory.")

    cell_all = pd.concat(cell_frames, ignore_index=True)
    donor_all = pd.concat(donor_frames, ignore_index=True)

    cell_all.to_parquet(RESULTS / "cell_scores_all.parquet", index=False)
    donor_all.to_csv(RESULTS / "donor_scores_all.tsv", sep="\t", index=False)

    n_done = len(RUN_ORDER) - len(skipped)
    print(f"cell_scores_all: {len(cell_all):,} rows, "
          f"{cell_all['module'].nunique()} modules across {n_done}/{len(RUN_ORDER)} CTs")
    print(f"donor_scores_all: {len(donor_all):,} rows, {donor_all['sample'].nunique()} donors")
    if skipped:
        print(f"[NOTE] skipped CTs (re-run 502 after these complete): {skipped}")


if __name__ == "__main__":
    main()
