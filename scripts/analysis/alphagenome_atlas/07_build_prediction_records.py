#!/usr/bin/env python3
"""Step 07: flatten archived Atlas responses to panel-track records and liver summaries.

Outputs (tables/):
  prediction_records/<scorer>.parquet   rows = variant × gene/junction; columns = obs + per-panel-track raw/quantile
  liver_summaries.tsv.gz                per variant × scorer × gene: median/min/max over primary-liver+hepatocyte tracks
  avi_records.tsv.gz                    AVI score, quantile, 18 feature importances, 18 model features
  atlas_availability.tsv                requested vs returned per variant per scorer
  pilot_anchor.json                     SORT1 anchor check from the pilot archive
"""

from __future__ import annotations

import json
import math
import pathlib
from collections import defaultdict

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import atlas_archive as aa
import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
PANEL = json.load((la.SCRIPT_DIR / "tissue_panel.json").open())
ANCHOR_UID, ANCHOR_GENE = "chr1:109274968:G:T", "ENSG00000134243"


def panel_group(curie: str, name: str) -> str | None:
    for group, curies in PANEL["groups"].items():
        if curie in curies:
            return group
    for group, sub in PANEL["name_fallback"].items():
        if sub in str(name).lower():
            return group
    return None


def track_columns(var: pd.DataFrame, scorer: str) -> list[tuple[int, str, str]]:
    """(column index, track label, group) for panel tracks; all tracks for always-keep scorers."""
    keep = []
    names = var["name"].astype(str).tolist() if "name" in var else [str(i) for i in var.index]
    curies = var["ontology_curie"].astype(str).tolist() if "ontology_curie" in var else [""] * len(var)
    bios = var["biosample_name"].astype(str).tolist() if "biosample_name" in var else [""] * len(var)
    for j in range(len(var)):
        group = panel_group(curies[j], bios[j])
        if scorer in PANEL["always_keep_scorers"]:
            group = group or "feature"
        if group:
            keep.append((j, names[j], group))
    return keep


OBS_COLUMNS = ["variant_uid", "gene_id", "gene_name", "strand", "junction_Start", "junction_End"]


def obs_frame(adata: anndata.AnnData) -> pd.DataFrame:
    if adata.shape[0] == 0 or "variant" not in adata.obs.columns:
        # a chunk can carry zero variants for a scorer (no junction / polyA site near any queried variant)
        return pd.DataFrame({c: pd.Series(dtype=str) for c in OBS_COLUMNS})
    obs = adata.obs.reset_index(drop=True).copy()
    obs["variant_uid"] = [aa.variant_uid_from_str(str(v)) for v in obs["variant"]]
    for col in ("gene_id", "gene_name", "strand", "junction_Start", "junction_End"):
        if col not in obs:
            obs[col] = ""
    obs["gene_id"] = obs["gene_id"].astype(str).str.split(".").str[0].replace({"nan": "", "None": ""})
    return obs[OBS_COLUMNS].astype(str)


def main() -> None:
    out_dir = TABLES / "prediction_records"
    if out_dir.exists():
        raise la.ContractError(f"refusing to overwrite: {out_dir}")
    out_dir.mkdir()
    meta = {r["scorer"]: r for r in la.read_tsv(TABLES / "scorer_metadata.tsv")}
    chunks = sorted((RAW / "atlas_direct").glob("chunk_*")) + sorted((RAW / "atlas_enzyme").glob("chunk_*"))
    if not chunks:
        raise la.ContractError("no archived chunks under raw/atlas_direct or raw/atlas_enzyme")
    writers: dict[str, pq.ParquetWriter] = {}
    avi_rows, avail = [], defaultdict(lambda: defaultdict(int))
    liver_cols = ["variant_uid", "scorer", "is_signed", "gene_id", "gene_name", "gene_strand", "junction_start", "junction_end",
                  "n_liver_tracks", "liver_median_raw", "liver_min_raw", "liver_max_raw", "liver_median_quantile",
                  "liver_min_quantile", "liver_max_quantile"]
    liver_path = TABLES / "liver_summaries.tsv.gz"
    if liver_path.exists():
        raise la.ContractError(f"refusing to overwrite: {liver_path}")
    liver_handle = la.open_text(liver_path, "wt")
    liver_writer = __import__("csv").DictWriter(liver_handle, fieldnames=liver_cols, delimiter="\t", lineterminator="\n")
    liver_writer.writeheader()
    n_liver = 0
    liver_groups = set(PANEL["liver_summary_groups"])
    seen_variants: set[str] = set()
    for ci, chunk in enumerate(chunks):
        req = json.load((chunk / "request.json").open())
        for u in req["variants"]:
            if u in seen_variants:
                raise la.ContractError(f"variant archived twice: {u} in {chunk}")
            seen_variants.add(u)
        for h5 in sorted(chunk.glob("*.h5ad")):
            scorer = h5.stem
            adata = anndata.read_h5ad(h5)
            obs = obs_frame(adata)
            cols = track_columns(adata.var, scorer)
            if not cols or len(obs) == 0:
                continue                                  # nothing scored for this scorer in this chunk
            x = np.asarray(adata.X, dtype=np.float32)
            q = np.asarray(adata.layers["quantiles"], dtype=np.float32) if "quantiles" in adata.layers else None
            idx = [c[0] for c in cols]
            extra = {}
            for j, name, group in cols:
                extra[f"raw|{group}|{name}"] = x[:, j]
                if q is not None:
                    extra[f"q|{group}|{name}"] = q[:, j]
            frame = pd.concat([pd.DataFrame({"scorer": scorer}, index=obs.index), obs, pd.DataFrame(extra, index=obs.index)], axis=1)
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if scorer not in writers:
                writers[scorer] = pq.ParquetWriter(out_dir / f"{scorer}.parquet", table.schema, compression="zstd")
            if table.schema != writers[scorer].schema:
                table = table.select(writers[scorer].schema.names).cast(writers[scorer].schema)   # same columns, fixed dtypes
            writers[scorer].write_table(table)
            for u in obs["variant_uid"].unique():
                avail[u][scorer] += 1
            # liver summary per row
            liver_idx = [c[0] for c in cols if c[2] in liver_groups]
            if liver_idx:
                lx = x[:, liver_idx]
                lq = q[:, liver_idx] if q is not None else None
                for i in range(len(frame)):
                    vals = lx[i]
                    if np.all(np.isnan(vals)):
                        continue
                    n_liver += 1
                    liver_writer.writerow({
                        "variant_uid": obs.at[i, "variant_uid"], "scorer": scorer, "is_signed": meta[scorer]["is_signed"],
                        "gene_id": obs.at[i, "gene_id"], "gene_name": obs.at[i, "gene_name"], "gene_strand": obs.at[i, "strand"],
                        "junction_start": obs.at[i, "junction_Start"], "junction_end": obs.at[i, "junction_End"],
                        "n_liver_tracks": int(np.sum(~np.isnan(vals))), "liver_median_raw": float(np.nanmedian(vals)),
                        "liver_min_raw": float(np.nanmin(vals)), "liver_max_raw": float(np.nanmax(vals)),
                        "liver_median_quantile": float(np.nanmedian(lq[i])) if lq is not None else math.nan,
                        "liver_min_quantile": float(np.nanmin(lq[i])) if lq is not None else math.nan,
                        "liver_max_quantile": float(np.nanmax(lq[i])) if lq is not None else math.nan,
                    })
            if scorer.startswith("AVI_SCORE"):
                for i in range(len(frame)):
                    row = {"variant_uid": obs.at[i, "variant_uid"], "scorer": scorer}
                    for j, name, _ in cols:
                        row[name] = float(x[i, j])
                    if q is not None and scorer == "AVI_SCORE":
                        row["AVI_SCORE_quantile"] = float(q[i, 0])
                    avi_rows.append(row)
        if ci % 50 == 0:
            la.log(f"chunk {ci + 1}/{len(chunks)}: liver rows {n_liver}")
    for w in writers.values():
        w.close()
    liver_handle.close()
    avi = pd.DataFrame(avi_rows)
    avi_wide = avi.drop(columns="scorer").groupby("variant_uid").first().reset_index()
    avi_wide.to_csv(TABLES / "avi_records.tsv.gz", sep="\t", index=False)
    scorers = sorted(meta)
    avail_rows = [{"variant_uid": u, **{s: avail[u].get(s, 0) for s in scorers}} for u in sorted(seen_variants)]
    la.write_tsv_once(TABLES / "atlas_availability.tsv.gz", avail_rows, ["variant_uid", *scorers])
    # Anchor from the pilot archive.
    anchor = {}
    for chunk in sorted((RAW / "pilot").glob("chunk_*")):
        req = json.load((chunk / "request.json").open())
        if ANCHOR_UID not in req["variants"]:
            continue
        a = anndata.read_h5ad(chunk / "RNA_SEQ.h5ad")
        sub = a[(a.obs["variant"] == ANCHOR_UID.replace(":G:T", ":G>T")) & (a.obs["gene_id"].astype(str).str.startswith(ANCHOR_GENE))]
        liver = sub.var["ontology_curie"].isin(PANEL["groups"]["primary_liver"] + PANEL["groups"]["hepatocyte"]).values
        vals = np.asarray(sub.X)[0, liver]
        anchor = {"variant": ANCHOR_UID, "gene": ANCHOR_GENE, "liver_rna_raw": vals.tolist(), "all_positive": bool(np.all(vals > 0)),
                  "liver_rna_quantile": np.asarray(sub.layers["quantiles"])[0, liver].tolist()}
    json.dump(anchor, (TABLES / "pilot_anchor.json").open("w"), indent=1)
    la.log(f"done: {len(seen_variants)} variants, {n_liver} liver summary rows, {len(avi_wide)} AVI rows; anchor positive={anchor.get('all_positive')}")


if __name__ == "__main__":
    main()
