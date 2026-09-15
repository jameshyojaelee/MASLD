"""Raw Atlas response archive: write-once per chunk, one .h5ad per scorer + request.json.

Also flattens an AnnData into long prediction records (variant × gene/junction × track).
No network access here.
"""

from __future__ import annotations

import json
import math
import pathlib
import re
from typing import Iterator, Mapping

import anndata
import numpy as np

import lib_atlas as la

VARIANT_RE = re.compile(r"^(chr[^:]+):(\d+):([ACGTN]+)>([ACGTN]+)$")


def _variant_str(value) -> str:
    """genome.Variant renders as 'chr1:5:A>G'; keep strings as they are."""
    return str(value)


def variant_uid_from_str(value: str) -> str:
    m = VARIANT_RE.match(value)
    if not m:
        raise la.ContractError(f"unexpected variant string: {value}")
    return la.variant_uid(m.group(1), int(m.group(2)), m.group(3), m.group(4))


def archive_scores(results: Mapping[str, anndata.AnnData], out_dir: pathlib.Path, request: dict) -> None:
    out_dir = pathlib.Path(out_dir)
    if out_dir.exists():
        raise la.ContractError(f"refusing to overwrite archive chunk: {out_dir}")
    out_dir.mkdir(parents=True)
    for scorer, adata in results.items():
        a = adata.copy()
        if "variant" in a.obs:
            a.obs["variant"] = [_variant_str(v) for v in a.obs["variant"]]
        for col in list(a.obs.columns):
            if a.obs[col].dtype == object:
                a.obs[col] = a.obs[col].astype(str)
        for col in list(a.var.columns):
            if a.var[col].dtype == object:
                a.var[col] = a.var[col].astype(str)
        a.write_h5ad(out_dir / f"{scorer}.h5ad")
    with (out_dir / "request.json").open("w") as handle:
        json.dump(request, handle, indent=1)


def load_archive(chunk_dir: pathlib.Path) -> dict[str, anndata.AnnData]:
    chunk_dir = pathlib.Path(chunk_dir)
    return {p.stem: anndata.read_h5ad(p) for p in sorted(chunk_dir.glob("*.h5ad"))}


def long_records(adata: anndata.AnnData, scorer: str, is_signed: bool) -> Iterator[dict]:
    x = np.asarray(adata.X, dtype=np.float32)
    q = np.asarray(adata.layers["quantiles"], dtype=np.float32) if "quantiles" in adata.layers else None
    obs = adata.obs.reset_index(drop=True)
    var = adata.var.reset_index(drop=True)
    track_names = var["name"] if "name" in var else var.index.astype(str)
    curie = var["ontology_curie"] if "ontology_curie" in var else [""] * len(var)
    bios = var["biosample_name"] if "biosample_name" in var else [""] * len(var)
    btype = var["biosample_type"] if "biosample_type" in var else [""] * len(var)
    assay = var["Assay title"] if "Assay title" in var else ([""] * len(var))
    strand_t = var["strand"] if "strand" in var else [""] * len(var)
    classes = [la.track_class(c, b, t) for c, b, t in zip(curie, bios, btype)]
    for i in range(x.shape[0]):
        o = obs.iloc[i]
        uid = variant_uid_from_str(str(o["variant"]))
        base = {
            "variant_uid": uid, "scorer": scorer, "is_signed": is_signed,
            "gene_id": str(o["gene_id"]) if "gene_id" in obs and not _isna(o["gene_id"]) else "",
            "gene_name": str(o["gene_name"]) if "gene_name" in obs and not _isna(o["gene_name"]) else "",
            "gene_strand": str(o["strand"]) if "strand" in obs and not _isna(o["strand"]) else "",
            "junction_start": o["junction_Start"] if "junction_Start" in obs and not _isna(o["junction_Start"]) else "",
            "junction_end": o["junction_End"] if "junction_End" in obs and not _isna(o["junction_End"]) else "",
        }
        for j in range(x.shape[1]):
            yield {
                **base, "track_name": str(track_names[j]), "ontology_curie": str(curie[j]), "biosample_name": str(bios[j]),
                "biosample_type": str(btype[j]), "assay": str(assay[j]), "track_strand": str(strand_t[j]),
                "track_class": classes[j], "raw_score": float(x[i, j]),
                "quantile": float(q[i, j]) if q is not None else math.nan, "quantile_available": q is not None,
            }


def _isna(v) -> bool:
    try:
        return v is None or (isinstance(v, float) and math.isnan(v)) or str(v) in ("nan", "None", "")
    except Exception:  # noqa: BLE001
        return False
