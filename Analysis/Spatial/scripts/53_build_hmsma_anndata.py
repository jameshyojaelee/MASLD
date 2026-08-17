#!/usr/bin/env python3
"""Build AnnData objects from HMSMA (HRA007511) STARsolo output.

Why this exists
---------------
`02_build_anndata.py` reads SpaceRanger `outs/filtered_feature_bc_matrix.h5`.
HMSMA was processed with STARsolo (no H&E images were available at processing
time), which emits `Solo.out/Gene/raw/{matrix.mtx.gz,barcodes.tsv.gz,
features.tsv.gz}` and **no spatial coordinates**. This script bridges that gap.

Spatial coordinates are reconstructed from the fixed Visium v1 slide layout:
every spot barcode maps to a permanent (array_row, array_col) on the array, so
coordinates need no H&E image. Verified 4,992/4,992 barcodes join exactly.

⛔ PHENOTYPE LABELS ARE NOT ATTACHED. Audited 2026-08-09: no key links the GSA
sample identifiers used here (`HRA_11` = manifest `sample_name` 11) to the
clinical identifiers carrying diagnosis/Kleiner/NAS (`CTRL-5113`, `MASL-1479`).
The GSA run manifest labels ALL 40 spatial runs `case set` / `fatty liver
disease`; the paper's supplementary tables carry 90 clinical IDs and zero
accessions; the GSA submission export carries 591 accessions and zero clinical
IDs. Until the HMSMA processed bundle (whose filenames ARE clinical IDs) or an
author-supplied key is obtained, these objects support unsupervised analyses
only — NOT disease contrasts, stage-ordering, or any labelled comparison.
`adata.uns['phenotype_key_available'] = False` records this.

⚠ TISSUE MASKING IS A PROXY. Without H&E there is no authoritative on-tissue
call. Spots are retained by a UMI threshold, recorded in `uns` as an explicit
analysis choice — it is not an imaging-derived fact.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
STAR_ROOT = PROJ / "Analysis/Spatial/results/starsolo_visium/HRA007511"
SAMPLES_TXT = PROJ / "Analysis/Spatial/metadata/HRA007511_samples.txt"
COORDS = Path(
    "/nfs/sw/easybuild/software/SpaceRanger/2.0.0-GCC-11.2.0"
    "/lib/python/cellranger/barcodes/visium-v1_coordinates.txt"
)
OUT_DIR = PROJ / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo"

VISIUM_V1_N_SPOTS = 4992

# Visium v1 is a hexagonal lattice and its array indices are ANISOTROPIC. Within a
# row, adjacent spots are 2 array_col units apart at a 100 um centre-to-centre
# pitch, so one col unit is 50 um. Adjacent rows are 1 array_row unit apart and
# offset by 1 col unit, so the diagonal neighbour distance sqrt(50^2 + r^2) = 100
# gives one row unit = 100*sqrt(3)/2 um.
#
# Left as raw indices, the six true neighbours are NOT equidistant: the diagonals
# land at 1.414 and the (col +/- 2) pair at 2.0, which ties against the
# (row +/- 2) spots — a full row-pair away and not neighbours at all. A
# 6-nearest-neighbour graph then breaks that tie arbitrarily and silently swaps
# real neighbours for spots ~173 um away. Emitting micrometres makes all six true
# neighbours land at exactly 100 um and the next shell at 173 um.
VISIUM_V1_SPOT_PITCH_UM = 100.0
COL_UNIT_UM = VISIUM_V1_SPOT_PITCH_UM / 2.0
ROW_UNIT_UM = VISIUM_V1_SPOT_PITCH_UM * np.sqrt(3.0) / 2.0


def load_coordinates() -> pd.DataFrame:
    """barcode -> (array_row, array_col) for the fixed Visium v1 layout."""
    df = pd.read_csv(COORDS, sep="\t", header=None,
                     names=["barcode", "array_row", "array_col"])
    df["barcode"] = df["barcode"].str.split("-").str[0]
    assert len(df) == VISIUM_V1_N_SPOTS, f"expected {VISIUM_V1_N_SPOTS} coords, got {len(df)}"
    assert df["barcode"].is_unique, "duplicate barcodes in coordinate map"
    return df.set_index("barcode")


def read_starsolo(sample_dir: Path) -> ad.AnnData:
    """Read one STARsolo Solo.out/Gene/raw matrix into AnnData (spots x genes)."""
    raw = sample_dir / "Solo.out" / "Gene" / "raw"
    mtx = raw / "matrix.mtx.gz"
    if not mtx.exists():
        mtx = raw / "matrix.mtx"

    # STARsolo writes features x barcodes; AnnData wants obs=spots, var=genes.
    with (gzip.open(mtx, "rb") if str(mtx).endswith(".gz") else open(mtx, "rb")) as fh:
        m = scipy.io.mmread(fh).tocsr()

    def _read_tsv(stem: str) -> pd.DataFrame:
        p = raw / f"{stem}.tsv.gz"
        if not p.exists():
            p = raw / f"{stem}.tsv"
        return pd.read_csv(p, sep="\t", header=None)

    features = _read_tsv("features")
    barcodes = _read_tsv("barcodes")

    assert m.shape == (len(features), len(barcodes)), (
        f"{sample_dir.name}: matrix {m.shape} vs "
        f"features {len(features)} x barcodes {len(barcodes)}"
    )

    X = m.T.tocsr()  # -> spots x genes
    var = pd.DataFrame(index=features[0].astype(str).values)
    if features.shape[1] > 1:
        var["gene_symbol"] = features[1].astype(str).values
    obs = pd.DataFrame(index=barcodes[0].astype(str).str.split("-").str[0].values)

    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    return a


def build_sample(sample: str, coords: pd.DataFrame, min_umi: int) -> ad.AnnData | None:
    sdir = STAR_ROOT / sample
    raw = sdir / "Solo.out" / "Gene" / "raw"
    if not (list(raw.glob("matrix.mtx*"))):
        print(f"  SKIP {sample}: no matrix", flush=True)
        return None

    a = read_starsolo(sdir)

    # --- coordinate join: must be exact, this is the whole point of the script
    hit = a.obs_names.isin(coords.index)
    n_hit = int(hit.sum())
    assert n_hit == a.n_obs, (
        f"{sample}: only {n_hit}/{a.n_obs} barcodes found in the Visium v1 "
        f"coordinate map — wrong slide version?"
    )
    a.obs["array_row"] = coords.loc[a.obs_names, "array_row"].values
    a.obs["array_col"] = coords.loc[a.obs_names, "array_col"].values
    # squidpy/scanpy convention: obsm['spatial'] as (x, y), here in micrometres so
    # that Euclidean distance is physical and the hex neighbourhood is recovered.
    a.obsm["spatial"] = np.column_stack(
        [a.obs["array_col"].values * COL_UNIT_UM,
         a.obs["array_row"].values * ROW_UNIT_UM]
    ).astype(float)

    a.obs["sample_id"] = sample
    a.obs["total_counts"] = np.asarray(a.X.sum(axis=1)).ravel()
    a.obs["n_genes"] = np.asarray((a.X > 0).sum(axis=1)).ravel()

    n_before = a.n_obs
    a = a[a.obs["total_counts"] >= min_umi].copy()
    print(
        f"  {sample}: {n_before} spots -> {a.n_obs} retained "
        f"(UMI>={min_umi}); median UMI {np.median(a.obs['total_counts']):.0f}",
        flush=True,
    )
    return a


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-umi", type=int, default=500,
                    help="on-tissue proxy threshold (no H&E available)")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    # Sealed candidates (the program arm and the mask-sensitivity diagnostic) read
    # the existing build directly, so a rebuild must land somewhere new.
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise SystemExit(f"refusing to overwrite populated build directory: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    coords = load_coordinates()
    samples = [s.strip() for s in SAMPLES_TXT.read_text().split() if s.strip()]
    print(f"samples listed: {len(samples)}; UMI threshold {args.min_umi}", flush=True)

    built, missing = [], []
    for s in samples:
        a = build_sample(s, coords, args.min_umi)
        if a is None:
            missing.append(s)
            continue
        a.write_h5ad(args.out_dir / f"{s}.h5ad")
        built.append(a)

    if not built:
        raise SystemExit("no samples built")

    merged = ad.concat(built, label="sample_id", keys=[a.obs["sample_id"][0] for a in built],
                       index_unique="-", merge="same")
    merged.uns["dataset"] = "HRA007511_HMSMA"
    merged.uns["processing"] = "STARsolo (GENCODE v49) — matches the bulk RNA-seq reference"
    merged.uns["coordinate_source"] = "visium-v1_coordinates.txt (fixed slide layout; no H&E)"
    merged.uns["coordinate_units"] = (
        f"micrometres; array indices rescaled to the physical hex lattice "
        f"(col x {COL_UNIT_UM:g}, row x {ROW_UNIT_UM:.4f}) so that all six "
        f"neighbours sit at {VISIUM_V1_SPOT_PITCH_UM:g} um"
    )
    merged.uns["tissue_mask"] = f"UMI>={args.min_umi} proxy — analysis choice, NOT imaging-derived"
    merged.uns["phenotype_key_available"] = False
    merged.uns["phenotype_key_note"] = (
        "No key links GSA sample_name (HRA_xx) to clinical IDs (CTRL-/MASL-). "
        "Unsupervised analyses only; no disease contrast or stage-ordering."
    )
    merged.uns["samples_missing_matrix"] = missing
    merged.write_h5ad(args.out_dir / "HRA007511_merged.h5ad")

    summary = {
        "n_samples_built": len(built),
        "n_samples_missing": len(missing),
        "samples_missing": missing,
        "n_spots_total": int(merged.n_obs),
        "n_genes": int(merged.n_vars),
        "min_umi": args.min_umi,
        "phenotype_key_available": False,
    }
    (args.out_dir / "build_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    print(f"WROTE {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
