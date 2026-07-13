#!/usr/bin/env python3
"""Landing-page mini-viz previews for the MASLD atlas home bento.

These three tiny JSONs are DECORATIVE PREVIEWS derived from already-baked
portal parquets (not new statistics). They back the animated landing tiles:

  - volcano_preview.json  -> the "Atlas" bento tile: a bulk-DE volcano cloud
       (~400 points {symbol, x=bulk_logFC, y=-log10(bulk_padj), sig}). All
       significant genes are eligible; capped by a seeded subsample so the
       classic two-wing shape survives, then padded with a seeded subsample of
       non-significant genes to ~400 total.
  - umap_thumbnail.json   -> the "Single-cell" bento tile: a ~2,500-point UMAP
       scatter {x, y, celltype}, resampled stratified-by-cell-type (seeded, with
       a per-type floor) from the downsampled UMAP so rare types stay visible.
  - landing_ticker.json   -> the gene marquee: top ~24 significant DEGs by
       |bulk_logFC| {symbol, logFC, dir}.

Sources (columns verified against the parquet schema before coding):
  atlas_core.parquet          human_symbol, bulk_logFC, bulk_padj, is_deg (bool)
  sc_umap_downsampled.parquet  umap_x, umap_y, cell_type

Env: any env with pandas + pyarrow (rnaseq or spatial). Honors
MASLD_PROJECT_ROOT and --output-dir. Compact JSON, NaN/Inf -> None (via
_portal_io.dump_json). All sampling uses a fixed seed for reproducibility.

Run (compute node; never the login node):
  srun --partition=cpu --qos=interactive --ntasks=1 --mem=8G --cpus-per-task=2 \
    --time=1:00:00 --job-name=portal \
    micromamba run -n rnaseq python scripts/portal/generate_landing_previews.py
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from _portal_io import dump_json, project_root

SEED = 42

# Volcano preview
VOLCANO_MAX_SIG = 250       # cap on significant points kept
VOLCANO_TOTAL = 400         # target total points (sig + non-sig)
VOLCANO_Y_CLIP = 50.0       # clip -log10(padj) so padj==0 -> inf doesn't blow up

# UMAP thumbnail
UMAP_TOTAL = 2500           # target total points
UMAP_FLOOR_PER_CT = 20      # min points per cell type (so rare types survive)

# Ticker
TICKER_N = 24               # top-N DEGs by |logFC|


def build_volcano(atlas: pd.DataFrame) -> list[dict]:
    """~400-point volcano: all eligible sig genes (seeded-capped at 250) plus a
    seeded subsample of non-sig genes to reach ~400."""
    df = atlas[["human_symbol", "bulk_logFC", "bulk_padj", "is_deg"]].copy()
    df["human_symbol"] = df["human_symbol"].astype(str)
    df = df[df["human_symbol"].str.len().gt(0)
            & ~df["human_symbol"].isin(["nan", "NA", "None"])]
    df["bulk_logFC"] = pd.to_numeric(df["bulk_logFC"], errors="coerce")
    df["bulk_padj"] = pd.to_numeric(df["bulk_padj"], errors="coerce")
    # need finite logFC and a padj we can -log10 (padj > 0..1; drop <=0/NaN)
    df = df[np.isfinite(df["bulk_logFC"]) & df["bulk_padj"].notna()]
    df = df[df["bulk_padj"] > 0]
    df["sig"] = df["is_deg"].astype(bool)

    sig = df[df["sig"]]
    nonsig = df[~df["sig"]]

    # cap sig via a seeded subsample (keeps the natural two-wing density)
    n_sig = min(VOLCANO_MAX_SIG, len(sig))
    sig_keep = sig if len(sig) <= n_sig else sig.sample(n=n_sig, random_state=SEED)

    # pad with a seeded subsample of non-sig to reach the ~400 target
    n_nonsig = max(0, VOLCANO_TOTAL - len(sig_keep))
    n_nonsig = min(n_nonsig, len(nonsig))
    nonsig_keep = (nonsig if len(nonsig) <= n_nonsig
                   else nonsig.sample(n=n_nonsig, random_state=SEED))

    keep = pd.concat([sig_keep, nonsig_keep], ignore_index=True)

    out = []
    for _, row in keep.iterrows():
        y = -math.log10(row["bulk_padj"])
        if y > VOLCANO_Y_CLIP:
            y = VOLCANO_Y_CLIP
        out.append({
            "symbol": row["human_symbol"],
            "x": round(float(row["bulk_logFC"]), 3),
            "y": round(y, 2),
            "sig": bool(row["sig"]),
        })
    return out


def build_umap(umap: pd.DataFrame) -> list[dict]:
    """~2,500-point UMAP scatter, stratified by cell type with a per-type floor so
    rare types survive; seeded."""
    df = umap[["umap_x", "umap_y", "cell_type"]].copy()
    df["umap_x"] = pd.to_numeric(df["umap_x"], errors="coerce")
    df["umap_y"] = pd.to_numeric(df["umap_y"], errors="coerce")
    df = df[np.isfinite(df["umap_x"]) & np.isfinite(df["umap_y"])]
    df["cell_type"] = df["cell_type"].astype(str)
    df = df[df["cell_type"].str.len().gt(0)
            & ~df["cell_type"].isin(["nan", "NA", "None"])]
    df = df.reset_index(drop=True)

    n_total = len(df)
    counts = df["cell_type"].value_counts()

    # proportional target per type, bumped to a floor so rare types stay visible
    picks = []
    for ct, grp in df.groupby("cell_type", sort=True):
        size = len(grp)
        prop = round(UMAP_TOTAL * size / n_total)
        target = min(max(prop, min(UMAP_FLOOR_PER_CT, size)), size)
        picks.append(grp.sample(n=target, random_state=SEED))

    sampled = pd.concat(picks, ignore_index=True)
    sampled = sampled.sort_values(["cell_type"]).reset_index(drop=True)

    return [
        {
            "x": round(float(r_["umap_x"]), 2),
            "y": round(float(r_["umap_y"]), 2),
            "celltype": r_["cell_type"],
        }
        for _, r_ in sampled.iterrows()
    ]


def build_ticker(atlas: pd.DataFrame) -> list[dict]:
    """Top ~24 significant DEGs by |bulk_logFC|."""
    df = atlas[["human_symbol", "bulk_logFC", "is_deg"]].copy()
    df["human_symbol"] = df["human_symbol"].astype(str)
    df["bulk_logFC"] = pd.to_numeric(df["bulk_logFC"], errors="coerce")
    df = df[df["is_deg"].astype(bool)
            & np.isfinite(df["bulk_logFC"])
            & df["human_symbol"].str.len().gt(0)
            & ~df["human_symbol"].isin(["nan", "NA", "None"])
            # drop unnamed genes whose symbol is a raw ENSEMBL id (e.g. ENSG000...)
            & ~df["human_symbol"].str.match(r"^ENSG\d", case=False, na=False)]
    df = df.reindex(df["bulk_logFC"].abs().sort_values(ascending=False).index)
    top = df.head(TICKER_N)
    return [
        {
            "symbol": row["human_symbol"],
            "logFC": round(float(row["bulk_logFC"]), 3),
            "dir": "up" if row["bulk_logFC"] > 0 else "down",
        }
        for _, row in top.iterrows()
    ]


def main() -> None:
    root = project_root()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", default=str(root / "masld-atlas-v2/public/data"))
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = root / out_dir

    atlas_path = out_dir / "atlas_core.parquet"
    umap_path = out_dir / "sc_umap_downsampled.parquet"
    for p in (atlas_path, umap_path):
        if not p.exists():
            raise SystemExit(f"Missing source parquet: {p}")

    print(f"Reading {atlas_path.name} ...")
    atlas = pd.read_parquet(
        atlas_path, columns=["human_symbol", "bulk_logFC", "bulk_padj", "is_deg"]
    )
    print(f"  {len(atlas):,} rows")
    print(f"Reading {umap_path.name} ...")
    umap = pd.read_parquet(umap_path, columns=["umap_x", "umap_y", "cell_type"])
    print(f"  {len(umap):,} rows")

    jobs = {
        "volcano_preview.json": build_volcano(atlas),
        "umap_thumbnail.json": build_umap(umap),
        "landing_ticker.json": build_ticker(atlas),
    }

    print("\n=== outputs ===")
    for name, data in jobs.items():
        path = out_dir / name
        kb = dump_json(data, path)
        print(f"  {path}  |  {len(data):,} rows  |  {path.stat().st_size:,} B "
              f"({kb:.1f} KB)")

    # console verification: first 2 elements of each
    print("\n=== first 2 elements per file ===")
    for name, data in jobs.items():
        print(f"  {name}: {data[:2]}")

    # extra sanity on the volcano split
    n_sig = sum(1 for d in jobs["volcano_preview.json"] if d["sig"])
    print(f"\nvolcano_preview: {n_sig} sig / "
          f"{len(jobs['volcano_preview.json']) - n_sig} non-sig")
    cts = {}
    for d in jobs["umap_thumbnail.json"]:
        cts[d["celltype"]] = cts.get(d["celltype"], 0) + 1
    print(f"umap_thumbnail: {len(cts)} cell types; per-type counts:")
    for ct in sorted(cts):
        print(f"    {ct}: {cts[ct]}")


if __name__ == "__main__":
    main()
