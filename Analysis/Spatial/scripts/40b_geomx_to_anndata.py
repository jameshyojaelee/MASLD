#!/usr/bin/env python
# ============================================================================
# Govaere 2026 GeoMx WTA: TSV -> AnnData wrapper
# ----------------------------------------------------------------------------
# Reads the Q3-normalized + raw segment x gene tables produced by
# 40_govaere2026_geomx_load.R and packages them as an AnnData object.
#
# .X        = Q3-normalized + log1p (matches paper's downstream model input)
# .layers['raw_counts']   = pre-normalization integer counts
# .layers['q3']           = Q3-normalized counts (no log)
# .layers['loq_bool']     = above-LOQ boolean indicator per segment x gene
# .obs                    = segment metadata
# .var                    = gene-level annotation (TargetName)
#
# OUTPUT:
#   Analysis/Spatial/results/preprocessed/geomx_govaere2026.h5ad
# ============================================================================

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse


def main() -> int:
    project_root = Path(
        os.environ.get(
            "MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
        )
    )
    in_dir = project_root / "Analysis/Spatial/results/govaere2026"
    out_dir = project_root / "Analysis/Spatial/results/preprocessed"
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_tsv = in_dir / "geomx_raw_counts.tsv"
    q3_tsv = in_dir / "geomx_expression.tsv"
    meta_tsv = in_dir / "geomx_metadata.tsv"
    loq_tsv = in_dir / "geomx_loq_matrix.tsv"

    for p in (raw_tsv, q3_tsv, meta_tsv):
        if not p.exists():
            sys.exit(f"missing input file: {p}")

    print(f"[INFO] Reading raw counts {raw_tsv}")
    raw_df = pd.read_csv(raw_tsv, sep="\t", index_col=0)
    print(f"[INFO] Reading Q3-normalized counts {q3_tsv}")
    q3_df = pd.read_csv(q3_tsv, sep="\t", index_col=0)
    print(f"[INFO] Reading metadata {meta_tsv}")
    meta_df = pd.read_csv(meta_tsv, sep="\t", index_col=0)

    # Align indices: TSV index = segment, columns = genes
    assert raw_df.shape == q3_df.shape, "raw and Q3 matrices differ in shape"

    # Reorder to match metadata index
    common_segments = meta_df.index.intersection(raw_df.index)
    raw_df = raw_df.loc[common_segments]
    q3_df = q3_df.loc[common_segments]
    meta_df = meta_df.loc[common_segments]
    common_genes = raw_df.columns.tolist()
    assert q3_df.columns.tolist() == common_genes

    print(f"[INFO] Aligned shape: {len(common_segments)} segments x {len(common_genes)} genes")

    # Build var
    var_df = pd.DataFrame({"gene": common_genes}, index=common_genes)

    # Build .X = log1p(Q3) as in standard NanoString downstream practice
    q3_mat = q3_df.values.astype(np.float32)
    raw_mat = raw_df.values.astype(np.float32)
    logq3 = np.log1p(q3_mat).astype(np.float32)

    adata = ad.AnnData(
        X=logq3,
        obs=meta_df.copy(),
        var=var_df.copy(),
        layers={"raw_counts": raw_mat, "q3": q3_mat},
    )

    # Optional LOQ matrix layer (only if the TSV has actual data; the R-side
    # writer was buggy in the first iteration and could produce a 1-column
    # file — skip cleanly in that case).
    if loq_tsv.exists():
        try:
            loq_df = pd.read_csv(loq_tsv, sep="\t", index_col=0)
            if loq_df.shape[1] < 2:
                print(f"[WARN] LOQ matrix only has {loq_df.shape[1]} data columns — skipping layer")
            else:
                loq_df = loq_df.transpose()  # rows = segments, cols = genes
                loq_df = loq_df.reindex(index=common_segments, columns=common_genes)
                # Map R-style TRUE/FALSE strings (or bools) to numpy bool.
                # Use vectorised replace + astype to avoid the deprecated
                # DataFrame.applymap path.
                def _to_bool(v):
                    if pd.isna(v):
                        return False
                    if isinstance(v, (bool, np.bool_)):
                        return bool(v)
                    s = str(v).strip().upper()
                    return s in ("TRUE", "T", "1")
                bool_mat = loq_df.map(_to_bool).astype(bool)
                adata.layers["loq_bool"] = bool_mat.values
                print(f"[INFO] Added LOQ boolean layer ({int(adata.layers['loq_bool'].sum())} TRUE cells)")
        except Exception as e:  # pragma: no cover
            print(f"[WARN] Could not load LOQ matrix: {e}")

    # Save
    out_path = out_dir / "geomx_govaere2026.h5ad"
    adata.write_h5ad(out_path, compression="gzip")
    print(f"[OK] Wrote {out_path}  size_MB={out_path.stat().st_size / 1e6:.2f}")
    print(f"[OK] AnnData: X.shape={adata.X.shape} obs={adata.obs.shape} var={adata.var.shape}")

    # Sanity: shape consistency
    assert adata.X.shape == (adata.obs.shape[0], adata.var.shape[0]), \
        f"shape mismatch: X={adata.X.shape}, obs={adata.obs.shape}, var={adata.var.shape}"

    # Print breakdowns
    print("\n[INFO] Patient breakdown:")
    print(adata.obs["Patient"].value_counts().sort_index().to_string())
    print("\n[INFO] Region type breakdown:")
    print(adata.obs["region_type"].value_counts().to_string())
    print("\n[INFO] Segment marker breakdown:")
    print(adata.obs["segment_marker"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
