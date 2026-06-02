"""Build atlas_cnmf_tcells.h5ad by subsetting the global atlas to T-lineage cells.

Mirrors Analysis/SingleCell/scripts/400_mcp/410_prep_atlas_cnmf.py conventions.

Atlas inspection (2026-05-13) shows the integrated CellTypist atlas already uses
the single merged label "T cells" (186,335 cells) rather than CD4/CD8/Treg
sub-types.  T_LABELS therefore contains just that one string; cell_type_original
is preserved for any downstream sub-labeling.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import ATLAS_GLOBAL, INPUTS  # noqa: E402

# Exact label strings present in atlas obs["cell_type"] (confirmed 2026-05-13).
# The integrated atlas already merges all T sub-types under one label.
T_LABELS: list[str] = [
    "T cells",
]


def main(smoke: bool) -> None:
    print(f"[1/4] Loading atlas: {ATLAS_GLOBAL}")
    adata = ad.read_h5ad(ATLAS_GLOBAL)
    print(f"      {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    print(f"[2/4] Subsetting to T labels: {T_LABELS}")
    mask = adata.obs["cell_type"].isin(T_LABELS)
    n_T = int(mask.sum())
    assert n_T > 0, f"No cells matched T_LABELS: {T_LABELS}"
    print(f"      {n_T:,} T cells matched ({100 * n_T / adata.n_obs:.2f}% of atlas)")
    tcells = adata[mask].copy()

    # Preserve original label (identical to "T cells" here, but retained for
    # consistency with the other subset scripts and downstream sub-labeling).
    tcells.obs["cell_type_original"] = tcells.obs["cell_type"].astype(str)
    tcells.obs["cell_type"] = "T cells"

    if smoke:
        import numpy as np
        n = min(5000, tcells.n_obs)
        rng = np.random.default_rng(42)
        idx = rng.choice(tcells.n_obs, size=n, replace=False)
        tcells = tcells[idx].copy()
        print(f"      [smoke] keeping {n} cells")

    INPUTS.mkdir(parents=True, exist_ok=True)
    out_path = INPUTS / "atlas_cnmf_tcells.h5ad"
    print(f"[3/4] Writing: {out_path}")
    tcells.write_h5ad(out_path, compression="gzip")

    n_donors = tcells.obs["sample"].nunique()
    print(f"[4/4] Done. {tcells.n_obs:,} cells, {n_donors} donors")
    print(f"      cell_type unique: {tcells.obs['cell_type'].unique().tolist()}")
    print(f"      X_scvi present: {'X_scvi' in tcells.obsm}")

    # Sanity assertions (skip count checks in smoke mode)
    if not smoke:
        assert tcells.n_obs >= 50_000, f"Too few T cells: {tcells.n_obs}"
        assert n_donors >= 100, f"Too few donors: {n_donors}"
    assert list(tcells.obs["cell_type"].unique()) == ["T cells"], \
        f"Unexpected cell_type values: {tcells.obs['cell_type'].unique().tolist()}"
    print("[4/4] All assertions passed.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Build atlas_cnmf_tcells.h5ad")
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="Subsample to 5,000 cells for fast smoke test",
    )
    main(ap.parse_args().smoke)
