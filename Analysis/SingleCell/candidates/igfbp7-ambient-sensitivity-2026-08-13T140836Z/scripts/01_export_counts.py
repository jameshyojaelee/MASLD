#!/usr/bin/env python
"""Step 1/3 of the IGFBP7 ambient-RNA SENSITIVITY ANALYSIS.

SCOPE (binding): this is a sensitivity analysis only. No Hotspot program is
rediscovered, refit, reweighted, renamed or re-selected anywhere in this
candidate. The frozen 117-program registry
(program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv) is
read only, never written.

What this step does: export RAW integer counts from the integrated human liver
atlas, split BY SEQUENCING DATASET (never pooled across batches), covering ALL
annotated cell types so that the decontX ambient pool is the real ambient pool
of that run — in particular so fibroblast / stellate ECM transcripts are
available to be modelled as contamination of hepatocyte barcodes.

Transport format matches the already-validated in-repo precedent
(Analysis/SingleCell/scripts/glp1ra_decontx_export.py): column-chunked CSC
binaries + a per-gene checksum, because the `decontx` R env has no HDF5 reader.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SEED = 20260813
np.random.seed(SEED)

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
WORK = CAND / "work"
RES = CAND / "results"
WORK.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

GLOBAL_H5 = ROOT / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
HEP_H5 = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_hepatocytes.h5ad"
HS_ROOT = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
HEP_MODULE_GENES = HS_ROOT / "hepatocytes/module_genes.tsv"
FIB_MODULE_GENES = HS_ROOT / "fibroblasts/module_genes.tsv"

CHUNK = 250_000  # cells per chunk; keeps per-chunk nnz below the int32 limit

# Diagnostic control panels (task 4). These are NOT used to define, reweight or
# re-select any program; they only calibrate what an ambient signature looks
# like in this atlas.
FIBRO_POS_CONTROLS = ["COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "ACTA2", "PDGFRB", "TAGLN"]
HEP_NEG_CONTROLS = ["ALB", "APOA1", "APOB", "CYP2E1", "HP", "TF", "TTR", "SERPINA1"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    # ---- frozen hepatocyte scoring universe (cells) ----------------------
    print(f"[export] hepatocyte scoring atlas: {HEP_H5}", flush=True)
    hep = ad.read_h5ad(HEP_H5, backed="r")
    hep_cells = pd.Index(hep.obs_names.astype(str))
    n_hep_universe = len(hep_cells)
    del hep
    print(f"[export] hepatocyte scoring universe: {n_hep_universe:,} cells", flush=True)

    # ---- frozen hepatocyte program gene universe -------------------------
    hep_mod = pd.read_csv(HEP_MODULE_GENES, sep="\t")
    hep_program_genes = sorted(set(hep_mod["gene"].astype(str)))
    print(f"[export] hepatocyte program gene universe: {len(hep_program_genes)} genes", flush=True)

    fib_mod = pd.read_csv(FIB_MODULE_GENES, sep="\t")
    fib_program_genes = sorted(set(fib_mod["gene"].astype(str)))

    print(f"[export] reading {GLOBAL_H5} (backed)", flush=True)
    a = ad.read_h5ad(GLOBAL_H5, backed="r")
    if a.raw is not None:
        var_names = np.asarray(a.raw.var_names).astype(str)
        raw_source = "raw.X"
    else:
        var_names = np.asarray(a.var_names).astype(str)
        raw_source = "X"
    ngenes = len(var_names)
    (WORK / "genes.txt").write_text("\n".join(var_names) + "\n", encoding="utf-8")
    gene_pos = {g: i for i, g in enumerate(var_names)}

    # genes whose reconstruction is checksum-verified on the R side
    check_genes = [
        g for g in (["IGFBP7", "BICC1", "PDGFRA", "COL1A1", "DCN", "LUM", "ALB", "CYP2E1"])
        if g in gene_pos
    ]
    check_idx = {g: gene_pos[g] for g in check_genes}

    # genes for which R must report ambient diagnostics + write corrected counts
    report_genes = sorted(
        set(hep_program_genes)
        | set(FIBRO_POS_CONTROLS)
        | set(HEP_NEG_CONTROLS)
        | set(fib_program_genes)
    )
    report_genes = [g for g in report_genes if g in gene_pos]
    (WORK / "report_genes.txt").write_text("\n".join(report_genes) + "\n", encoding="utf-8")

    # genes whose CORRECTED counts are transported back for re-scoring: exactly
    # the frozen hepatocyte program member genes.
    score_genes = [g for g in hep_program_genes if g in gene_pos]
    missing_score_genes = [g for g in hep_program_genes if g not in gene_pos]
    (WORK / "score_genes.txt").write_text("\n".join(score_genes) + "\n", encoding="utf-8")
    print(
        f"[export] score genes present in atlas: {len(score_genes)}/{len(hep_program_genes)} "
        f"(missing {len(missing_score_genes)})",
        flush=True,
    )

    obs = a.obs
    ct = obs["cell_type"].astype(str).to_numpy()
    ds = obs["dataset"].astype(str).to_numpy()
    sample = obs["sample"].astype(str).to_numpy()
    cell_ids = obs.index.astype(str).to_numpy()
    # pandas Index.isin is hash-based; np.isin on object-dtype string arrays
    # degrades to a quadratic scan at this scale.
    in_hep_universe = np.asarray(pd.Index(cell_ids).isin(hep_cells))
    print(
        f"[export] global atlas {len(cell_ids):,} cells; "
        f"{int(in_hep_universe.sum()):,} are in the frozen hepatocyte scoring universe",
        flush=True,
    )

    # z for decontX: cell_type -> consecutive integers, defined once GLOBALLY so
    # the same integer means the same cell type in every dataset.
    ct_levels = sorted(set(ct))
    z_map = {c: i + 1 for i, c in enumerate(ct_levels)}
    (WORK / "z_levels.json").write_text(json.dumps(z_map, indent=2), encoding="utf-8")

    manifest = {
        "seed": SEED,
        "ngenes": int(ngenes),
        "raw_source": raw_source,
        "chunk": CHUNK,
        "check_genes": check_genes,
        "n_report_genes": len(report_genes),
        "n_score_genes": len(score_genes),
        "missing_score_genes": missing_score_genes,
        "z_map": z_map,
        "n_hep_universe_cells": int(n_hep_universe),
        "n_hep_universe_cells_found_in_global": int(in_hep_universe.sum()),
        "datasets": {},
        "skipped": {},
        "global_atlas": str(GLOBAL_H5),
        "hepatocyte_atlas": str(HEP_H5),
        "fibro_pos_controls": FIBRO_POS_CONTROLS,
        "hep_neg_controls": HEP_NEG_CONTROLS,
    }

    for dataset in sorted(set(ds)):
        sel = ds == dataset
        idx = np.where(sel)[0]
        ncells = len(idx)
        n_hep_here = int(in_hep_universe[idx].sum())
        n_ct = len(set(ct[idx]))
        if n_hep_here == 0:
            # Nothing in this dataset enters the frozen hepatocyte scoring
            # universe, so it cannot change any program score.
            print(f"[export] SKIP {dataset}: no hepatocyte-universe cells", flush=True)
            manifest["skipped"][dataset] = {"ncells": ncells, "n_hep_universe": 0}
            continue

        dsdir = WORK / dataset
        dsdir.mkdir(parents=True, exist_ok=True)
        print(
            f"[export] {dataset}: {ncells:,} cells, {n_ct} cell types, "
            f"{n_hep_here:,} in hepatocyte scoring universe",
            flush=True,
        )

        meta = pd.DataFrame(
            {
                "cell_id": cell_ids[idx],
                "cell_type": ct[idx],
                "sample": sample[idx],
                "z": [z_map[c] for c in ct[idx]],
                "in_hep_universe": in_hep_universe[idx].astype(int),
            }
        )
        meta.to_csv(dsdir / "meta.csv.gz", index=False)

        chunks = []
        check_rowsum = {g: 0.0 for g in check_idx}
        total_sum = 0.0
        for ci, start in enumerate(range(0, ncells, CHUNK)):
            cidx = idx[start : start + CHUNK]
            Xsub = a.raw.X[cidx, :] if a.raw is not None else a.X[cidx, :]
            if not sp.isspmatrix_csr(Xsub):
                Xsub = sp.csr_matrix(Xsub)
            Xgc = Xsub.T.tocsc()
            Xgc.sort_indices()
            Xgc.eliminate_zeros()
            nnz = int(Xgc.nnz)
            if nnz >= 2**31:
                raise RuntimeError(f"chunk nnz {nnz} exceeds int32; lower CHUNK")
            if Xgc.shape[0] != ngenes:
                raise RuntimeError("gene axis drift during export")
            if not np.allclose(Xgc.data, np.round(Xgc.data)):
                raise RuntimeError("raw counts are not integral; refusing to export")
            if Xgc.data.max() > np.iinfo(np.int32).max:
                raise RuntimeError("count exceeds int32")
            Xgc.data = np.round(Xgc.data)
            Xgc.data.astype("<i4").tofile(dsdir / f"chunk{ci}_data.bin")
            Xgc.indices.astype("<i4").tofile(dsdir / f"chunk{ci}_indices.bin")
            Xgc.indptr.astype("<i4").tofile(dsdir / f"chunk{ci}_indptr.bin")
            chunks.append({"chunk": ci, "ncells": int(Xgc.shape[1]), "nnz": nnz})
            total_sum += float(Xgc.data.sum())
            for g, gi in check_idx.items():
                check_rowsum[g] += float(Xsub[:, gi].sum())
            del Xsub, Xgc

        info = {
            "dataset": dataset,
            "ncells": int(ncells),
            "ngenes": int(ngenes),
            "n_cell_types": int(n_ct),
            "n_hep_universe": n_hep_here,
            "chunks": chunks,
            "total_count_sum": total_sum,
            "check_rowsum": {g: float(v) for g, v in check_rowsum.items()},
        }
        (dsdir / "dims.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
        manifest["datasets"][dataset] = {
            "ncells": int(ncells),
            "n_cell_types": int(n_ct),
            "n_hep_universe": n_hep_here,
            "n_chunks": len(chunks),
        }
        print(f"[export]   {dataset} done: total_count_sum={total_sum:.0f}", flush=True)

    (WORK / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # input checksum manifest (task: input/output checksum manifest)
    inputs = [GLOBAL_H5, HEP_H5, HEP_MODULE_GENES, FIB_MODULE_GENES]
    rows = [
        {
            "role": "input",
            "path": str(p),
            "bytes": p.stat().st_size,
            "sha256": sha256(p),
        }
        for p in inputs
    ]
    pd.DataFrame(rows).to_csv(RES / "01_input_checksums.tsv", sep="\t", index=False)

    print(
        f"[export] manifest written: datasets={sorted(manifest['datasets'])} "
        f"skipped={sorted(manifest['skipped'])}",
        flush=True,
    )
    print("[export] DONE", flush=True)


if __name__ == "__main__":
    main()
