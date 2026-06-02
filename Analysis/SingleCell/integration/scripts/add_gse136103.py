#!/usr/bin/env python3
"""Convert GSE136103 (Ramachandran 2019) processed count matrices to h5ad files.

The GEO-provided matrices (MTX format) have proper cell counts from merged lanes.
Per-SRR CellRanger outputs had 1-11 cells due to unmerged lanes — unusable.

Each GSM = one library (GEM well). Multiple GSMs per donor (CD45+/CD45- fractions).

Output:
  - Per-GSM h5ad files in integration/input_h5ad/human/
  - Updated sample_manifest.csv with GSE136103 entries
"""

import csv
import gzip
import logging
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import scipy.io as sio
import anndata as ad

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
PROCESSED_DIR = BASE / "data/GSE136103/processed"
OUTPUT_DIR = BASE / "Analysis/SingleCell/integration/input_h5ad/human"
MANIFEST_PATH = BASE / "Analysis/SingleCell/integration/input_h5ad/sample_manifest.csv"

# GSM metadata from GEO SOFT file + ENA
# Format: GSM_ID -> (sample_title, condition, donor, preparation_method)
GSM_METADATA = {
    "GSM4041150": ("Healthy1_Cd45+", "Healthy", "Healthy1", "cd45_positive"),
    "GSM4041151": ("Healthy1_Cd45-A", "Healthy", "Healthy1", "cd45_negative"),
    "GSM4041152": ("Healthy1_Cd45-B", "Healthy", "Healthy1", "cd45_negative"),
    "GSM4041153": ("Healthy2_Cd45+", "Healthy", "Healthy2", "cd45_positive"),
    "GSM4041154": ("Healthy2_Cd45-", "Healthy", "Healthy2", "cd45_negative"),
    "GSM4041155": ("Healthy3_Cd45+", "Healthy", "Healthy3", "cd45_positive"),
    "GSM4041156": ("Healthy3_Cd45-A", "Healthy", "Healthy3", "cd45_negative"),
    "GSM4041157": ("Healthy3_Cd45-B", "Healthy", "Healthy3", "cd45_negative"),
    "GSM4041158": ("Healthy4_Cd45+", "Healthy", "Healthy4", "cd45_positive"),
    "GSM4041159": ("Healthy4_Cd45-", "Healthy", "Healthy4", "cd45_negative"),
    "GSM4041160": ("Healthy5_Cd45+", "Healthy", "Healthy5", "cd45_positive"),
    "GSM4041161": ("Cirrhotic1_Cd45+", "Cirrhotic", "Cirrhotic1", "cd45_positive"),
    "GSM4041162": ("Cirrhotic1_Cd45-A", "Cirrhotic", "Cirrhotic1", "cd45_negative"),
    "GSM4041163": ("Cirrhotic1_Cd45-B", "Cirrhotic", "Cirrhotic1", "cd45_negative"),
    "GSM4041164": ("Cirrhotic2_Cd45+", "Cirrhotic", "Cirrhotic2", "cd45_positive"),
    "GSM4041165": ("Cirrhotic2_Cd45-", "Cirrhotic", "Cirrhotic2", "cd45_negative"),
    "GSM4041166": ("Cirrhotic3_Cd45+", "Cirrhotic", "Cirrhotic3", "cd45_positive"),
    "GSM4041167": ("Cirrhotic3_Cd45-", "Cirrhotic", "Cirrhotic3", "cd45_negative"),
    "GSM4041168": ("Cirrhotic4_Cd45+", "Cirrhotic", "Cirrhotic4", "cd45_positive"),
    "GSM4041169": ("Cirrhotic5_Cd45+", "Cirrhotic", "Cirrhotic5", "cd45_positive"),
    # Blood (PBMC) — excluded from liver integration
    # "GSM4041170": ("Blood1", "Blood", "Blood1", "unsorted"),
    # "GSM4041171": ("Blood2", "Blood", "Blood2", "unsorted"),
    # "GSM4041172": ("Blood3", "Blood", "Blood3", "unsorted"),
    # "GSM4041173": ("Blood4", "Blood", "Blood4", "unsorted"),
    # Mouse — excluded (different species)
    # "GSM4041174": ("Mouse_healthy", "Healthy", "Mouse_H", "unsorted"),
    # "GSM4041175": ("Mouse_fibrotic", "Mouse", "Mouse_F", "unsorted"),
}


def read_mtx_triplet(gsm_id: str) -> ad.AnnData:
    """Read MTX + barcodes + genes for a GSM sample."""
    # Find files matching this GSM
    prefix = None
    for f in PROCESSED_DIR.glob(f"{gsm_id}_*_matrix.mtx.gz"):
        prefix = str(f).replace("_matrix.mtx.gz", "")
        break

    if prefix is None:
        raise FileNotFoundError(f"No matrix file found for {gsm_id}")

    mtx_file = f"{prefix}_matrix.mtx.gz"
    bc_file = f"{prefix}_barcodes.tsv.gz"
    genes_file = f"{prefix}_genes.tsv.gz"

    # Read matrix
    mat = sio.mmread(mtx_file).T.tocsr()  # cells x genes

    # Read barcodes
    with gzip.open(bc_file, "rt") as f:
        barcodes = [line.strip() for line in f]

    # Read genes (tab-separated: gene_id \t gene_symbol)
    gene_ids = []
    gene_symbols = []
    with gzip.open(genes_file, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            gene_ids.append(parts[0])
            gene_symbols.append(parts[1] if len(parts) > 1 else parts[0])

    # Create AnnData
    import pandas as pd
    obs_df = pd.DataFrame(index=barcodes)
    var_df = pd.DataFrame({"gene_ids": gene_ids, "gene_symbols": gene_symbols}, index=gene_symbols)
    var_df.index.name = None

    adata = ad.AnnData(X=mat, obs=obs_df, var=var_df)
    adata.var_names_make_unique()

    return adata


def main():
    t_start = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    log.info(f"Processing {len(GSM_METADATA)} human liver samples from GSE136103")
    log.info(f"Output: {OUTPUT_DIR}")

    new_manifest_rows = []
    total_cells = 0

    for gsm_id, (title, condition, donor, prep) in sorted(GSM_METADATA.items()):
        out_path = OUTPUT_DIR / f"{gsm_id}.h5ad"

        if out_path.exists():
            log.info(f"  {gsm_id} ({title}): exists, skipping")
            # Still add to manifest
            new_manifest_rows.append({
                "srr": gsm_id,
                "dataset": "GSE136103",
                "species": "human",
                "condition": condition,
                "h5_path": str(out_path),
            })
            continue

        try:
            adata = read_mtx_triplet(gsm_id)

            # Add metadata
            adata.obs["sample"] = gsm_id
            adata.obs["dataset"] = "GSE136103"
            adata.obs["species"] = "human"
            adata.obs["condition"] = condition
            adata.obs["donor"] = donor
            adata.obs["preparation_method"] = prep
            adata.obs["sample_title"] = title

            adata.write_h5ad(out_path)
            total_cells += adata.n_obs
            log.info(f"  {gsm_id} ({title}): {adata.n_obs:,} cells, {adata.n_vars:,} genes → {out_path.name}")

            new_manifest_rows.append({
                "srr": gsm_id,
                "dataset": "GSE136103",
                "species": "human",
                "condition": condition,
                "h5_path": str(out_path),
            })

        except Exception as e:
            log.error(f"  {gsm_id}: FAILED — {e}")

    # Update manifest: read existing, remove old GSE136103 entries, append new
    existing_rows = []
    if MANIFEST_PATH.exists():
        with MANIFEST_PATH.open() as f:
            reader = csv.DictReader(f)
            existing_rows = [r for r in reader if r.get("dataset") != "GSE136103"]

    all_rows = existing_rows + new_manifest_rows
    with MANIFEST_PATH.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["srr", "dataset", "species", "condition", "h5_path"])
        writer.writeheader()
        writer.writerows(all_rows)

    dt = time.time() - t_start
    n_healthy = sum(1 for _, (_, c, _, _) in GSM_METADATA.items() if c == "Healthy")
    n_cirrhotic = sum(1 for _, (_, c, _, _) in GSM_METADATA.items() if c == "Cirrhotic")
    log.info(f"\nDone in {dt:.0f}s: {len(new_manifest_rows)} samples ({n_healthy} healthy, {n_cirrhotic} cirrhotic)")
    log.info(f"Total cells: {total_cells:,}")
    log.info(f"Manifest updated: {MANIFEST_PATH} ({len(all_rows)} total samples)")

    # Summary by condition x prep
    from collections import Counter
    prep_counts = Counter((c, p) for _, (_, c, _, p) in GSM_METADATA.items())
    log.info("\nCondition x Preparation:")
    for (c, p), n in sorted(prep_counts.items()):
        log.info(f"  {c} / {p}: {n}")


if __name__ == "__main__":
    main()
