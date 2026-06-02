#!/usr/bin/env python3
"""
03b_chromvar_per_donor_export.py — Aggregate per-cell chromVAR deviations
to donor x cell_type x TF (B5 of the ATAC improvement plan).

Reads:
  - results/chromvar_v2/chromvar_deviations.h5ad   (88,814 cells x 1,019 motifs)

Writes:
  - results/chromvar_v2/chromvar_per_donor.tsv.gz
      Long-format: donor_id_atac | cell_type | TF | mean_deviation | n_cells | sd_deviation

Notes
-----
The chromVAR h5ad has duplicated motif names (1,019 columns, 964 unique).
We retain all columns but deduplicate by suffix-tagging duplicates so that
joining downstream stays unambiguous (e.g., ATF3, ATF3.1, ...). The base
gene symbol (everything before the last '.') is preserved in a separate
column TF_symbol for cross-method joins.

Environment: snapatac2 (with PYTHONNOUSERSITE=1)
"""
from __future__ import annotations

import gzip
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATAC_DIR = PROJECT_ROOT / "Analysis" / "ATAC" / "Human_Multiome"
DEV_H5AD = ATAC_DIR / "results" / "chromvar_v2" / "chromvar_deviations.h5ad"
OUT_TSV = ATAC_DIR / "results" / "chromvar_v2" / "chromvar_per_donor.tsv.gz"


def main() -> None:
    import anndata as ad

    t0 = time.time()
    log.info("=" * 72)
    log.info("chromVAR per-donor x cell-type x TF export (B5)")
    log.info("=" * 72)
    log.info(f"Input:  {DEV_H5AD}")
    log.info(f"Output: {OUT_TSV}")

    if not DEV_H5AD.exists():
        log.error(f"Missing input h5ad: {DEV_H5AD}")
        sys.exit(1)

    # Load entirely in memory (~350 MB float32 dense) — fine on 64 GB node.
    adata = ad.read_h5ad(DEV_H5AD)
    n_cells, n_motifs = adata.shape
    log.info(f"Loaded {n_cells:,} cells x {n_motifs:,} motifs")

    # Build a unique-per-column TF identifier so duplicated symbols don't
    # collide in long-form output. The base symbol is preserved separately.
    base_names = list(adata.var_names)
    seen: dict[str, int] = {}
    tf_unique: list[str] = []
    tf_symbol: list[str] = []
    for nm in base_names:
        k = seen.get(nm, 0)
        seen[nm] = k + 1
        tf_unique.append(nm if k == 0 else f"{nm}.{k}")
        tf_symbol.append(nm)
    log.info(f"Unique TF column count after suffix tagging: {len(set(tf_unique))}")

    # Required obs columns
    needed = ("donor_id", "cell_type")
    for c in needed:
        if c not in adata.obs.columns:
            log.error(f"Missing obs column: {c}")
            sys.exit(1)

    # Materialize X as dense float32 (already saved that way in 03)
    X = adata.X
    if sp.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float32)
    else:
        X = np.asarray(X, dtype=np.float32)

    donors = adata.obs["donor_id"].astype(str).values
    cts = adata.obs["cell_type"].astype(str).values

    # Build (donor, cell_type) groups
    df_idx = pd.DataFrame({"donor_id": donors, "cell_type": cts})
    groups = df_idx.groupby(["donor_id", "cell_type"], observed=True).indices
    log.info(f"Donor x cell_type groups: {len(groups):,} (max possible {df_idx['donor_id'].nunique() * df_idx['cell_type'].nunique()})")

    # Aggregate
    records: list[tuple] = []
    for (donor, ct), idx in groups.items():
        sub = X[idx, :]  # (n_cells_g, n_motifs)
        n_cells_g = int(sub.shape[0])
        # Single-cell groups produce sd = nan; report nan explicitly.
        means = sub.mean(axis=0).astype(np.float32)
        if n_cells_g > 1:
            sds = sub.std(axis=0, ddof=1).astype(np.float32)
        else:
            sds = np.full(n_motifs, np.nan, dtype=np.float32)

        for k in range(n_motifs):
            records.append(
                (donor, ct, tf_unique[k], tf_symbol[k], float(means[k]), n_cells_g, float(sds[k]))
            )

    out = pd.DataFrame.from_records(
        records,
        columns=[
            "donor_id_atac",
            "cell_type",
            "TF",
            "TF_symbol",
            "mean_deviation",
            "n_cells",
            "sd_deviation",
        ],
    )
    log.info(f"Long-form rows: {len(out):,}")

    # Round for compactness (h5ad source was float32 anyway)
    out["mean_deviation"] = out["mean_deviation"].round(6)
    out["sd_deviation"] = out["sd_deviation"].round(6)

    out.to_csv(OUT_TSV, sep="\t", index=False, compression="gzip")
    log.info(f"Wrote {OUT_TSV} ({OUT_TSV.stat().st_size / 1e6:.2f} MB, {len(out):,} rows)")

    # Sanity preview
    log.info("Preview:")
    log.info("\n" + out.head(8).to_string(index=False))

    log.info(f"Done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
