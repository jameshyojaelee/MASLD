#!/usr/bin/env python
"""19_pseudobulk_gse281367.py  --  Layer A A3/A4 bridge.

Build GSE281367 donor x peak pseudobulk COUNT matrices on the SAME peak
coordinate set as GSE244832, so the two cohorts pool in 18_pooled_da_gse281367.R.

Reuses the GSE244832 tile->peak->donor logic verbatim (07b_corrected_da_hepatocytes.py
+ utils_pseudobulk.py), so GSE281367 counts have the IDENTICAL definition:
  per-donor SUM of 500bp tile counts, mapped to each GSE244832 cell-type peak via a
  peak x tile overlap indicator. Peak column names = "chr:start-end" (== GSE244832
  count columns). condition from GSE281367 donor_pairing.csv (MASH=1 / NORMAL=0).

Input : Analysis/ATAC/Human_External/snapatac2/snapatac2_processed.h5ad
        (produced by reusing 02_snapatac2_processing.py on GSE281367; obs must carry
         'cell_type' [Hepatocyte/Stellate_Cell/Macrophage/Cholangiocyte...] + 'donor_id';
         .X = 500bp tile matrix; var_names = tile coords "chr:start-end")
Output: Analysis/ATAC/Human_External/pseudobulk/{hep,stellate,macrophage,cholangiocyte}_
        pseudobulk_counts_GSE281367.tsv.gz + _coldata_GSE281367.tsv
Env   : snapatac2 (anndata + scipy)
"""
import os, sys, logging
import numpy as np
import pandas as pd
import anndata as ad

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MM_SCRIPTS = os.path.join(ROOT, "Analysis/ATAC/Human_Multiome/scripts")
sys.path.insert(0, MM_SCRIPTS)
# reuse the GSE244832 donor-collapse helper verbatim (identical count definition).
# 07b's module name starts with a digit (not importable), so its tiny pure helpers
# (load_peaks / parse_chrom_sizes_from_tiles / build_peak_to_tile_mapping) are
# re-defined below byte-for-byte from 07b_corrected_da_hepatocytes.py.
from utils_pseudobulk import aggregate_cells_to_donors                       # noqa: E402

def load_peaks(bed):
    P = []
    with open(bed) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3:
                P.append((p[0], int(p[1]), int(p[2])))
    return P

def parse_chrom_sizes_from_tiles(var_names):
    d, cur, s = {}, None, 0
    for i, nm in enumerate(var_names):
        c = nm.split(":", 1)[0]
        if c != cur:
            if cur is not None:
                d[cur] = (s, i - 1)
            cur, s = c, i
    d[cur] = (s, len(var_names) - 1)
    return d

TILE_SIZE = 500
def build_peak_to_tile_mapping(peaks, chrom_to_range, n_tiles):
    from scipy import sparse
    rows, cols = [], []
    for pi, (chrom, start, end) in enumerate(peaks):
        if chrom not in chrom_to_range:
            continue
        c0, c1 = chrom_to_range[chrom]
        gs = max(c0 + start // TILE_SIZE, c0)
        ge = min(c0 + (end - 1) // TILE_SIZE, c1)
        if ge < gs:
            continue
        for t in range(gs, ge + 1):
            rows.append(pi); cols.append(t)
    return sparse.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)),
                             shape=(len(peaks), n_tiles), dtype=np.float32)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("pb281367")

H5AD = os.environ.get("MASLD_ATAC_H5AD",
                      os.path.join(ROOT, "Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad"))
PEAK_DIR = os.path.join(ROOT, "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2")
META = os.path.join(ROOT, "data/GSE281367/metadata/donor_pairing.csv")
OUT = os.path.join(ROOT, "Analysis/ATAC/Human_External/pseudobulk")
os.makedirs(OUT, exist_ok=True)

# 02 cell_type label -> (GSE244832 peak BED, output tag matching 18_pooled_da)
CT_MAP = {
    "Hepatocyte":    ("Hepatocytes_peaks.bed",   "hep"),
    "Stellate_Cell": ("Fibroblasts_peaks.bed",   "stellate"),
    "Macrophage":    ("Macrophages_peaks.bed",   "macrophage"),
    "Cholangiocyte": ("Cholangiocytes_peaks.bed","cholangiocyte"),
}

def main():
    meta = pd.read_csv(META)
    cond_map = dict(zip(meta.donor_id.astype(str), meta.condition.astype(str)))
    log.info("Reading %s", H5AD)
    A = ad.read_h5ad(H5AD)
    assert "cell_type" in A.obs and "donor_id" in A.obs, "h5ad needs cell_type + donor_id"
    var_names = list(A.var_names)
    chrom_to_range = parse_chrom_sizes_from_tiles(var_names)
    n_tiles = len(var_names)
    log.info("Cells=%d tiles=%d donors=%s celltypes=%s",
             A.n_obs, n_tiles, sorted(A.obs.donor_id.unique()),
             A.obs.cell_type.value_counts().to_dict())

    for ct_label, (bed, tag) in CT_MAP.items():
        sub = A[A.obs.cell_type == ct_label]
        if sub.n_obs == 0:
            log.warning("[%s] 0 cells annotated -> skip", ct_label); continue
        peaks = load_peaks(os.path.join(PEAK_DIR, bed))
        S = build_peak_to_tile_mapping(peaks, chrom_to_range, n_tiles)   # (n_peak x n_tile)
        donors, M_tiles, ncells = aggregate_cells_to_donors(sub.X, sub.obs.donor_id.values, agg="sum")
        M_counts = np.asarray(S @ M_tiles.T).T                            # (n_donor x n_peak)
        peak_names = [f"{c}:{s}-{e}" for (c, s, e) in peaks]
        counts_int = np.rint(M_counts).astype(np.int64)
        df = pd.DataFrame(counts_int, index=list(donors), columns=peak_names)
        df.index.name = "donor_id"
        df.to_csv(os.path.join(OUT, f"{tag}_pseudobulk_counts_GSE281367.tsv.gz"),
                  sep="\t", compression="gzip")
        col = pd.DataFrame({
            "donor_id": list(donors),
            "condition": [cond_map.get(str(d), "NA") for d in donors],
            "libsize": counts_int.sum(axis=1),
            "n_cells": ncells,
        })
        col.to_csv(os.path.join(OUT, f"{tag}_pseudobulk_coldata_GSE281367.tsv"), sep="\t", index=False)
        log.info("[%s] wrote %s: %d donors x %d peaks (cells/donor med=%.0f)",
                 ct_label, tag, len(donors), len(peak_names), float(np.median(ncells)))

    log.info("Done. Pseudobulk in %s", OUT)

if __name__ == "__main__":
    main()
