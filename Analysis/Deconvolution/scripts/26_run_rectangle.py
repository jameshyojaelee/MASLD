#!/usr/bin/env python
"""Rectangle deconvolution runner (rectanglepy 1.5.0).

Deconvolutes a bulk RNA-seq TPM matrix against a single-cell reference AnnData
using the all-in-one ``rectanglepy.rectangle`` entry point.

API (introspected from the installed rectanglepy 1.5.0, not assumed):
  rectangle.rectangle(
      adata: AnnData,                 # single-cell reference, .X = raw counts,
                                      #   obs[cell_type_col] = annotations,
                                      #   var_names = gene symbols (must be unique)
      bulks: DataFrame,               # samples x genes, values in TPM
      cell_type_col: str = "cell_type",
      *, optimize_cutoffs=True, n_cpus=None, correct_mrna_bias=True,
      gene_expression_threshold=0.5, ...
  ) -> (estimations: DataFrame, signature: RectangleSignatureResult)

  estimations: rows = bulks.index (sample ids), columns = the reference cell
  types PLUS a trailing "Unknown" column; each row (incl. Unknown) sums to 1.0.

Coarse vs fine mode is inferred from ``cell_type_col``:
  - coarse (default, cell_type_col == "cell_type"):
      * separate the raw "Unknown" fraction  -> <DATASET>_rectangle_unknown.tsv
      * reindex the remaining fractions to the 16 CANONICAL cell types IN ORDER
        (0-fill any absent type), then RE-NORMALIZE the 16 to sum to 1 AFTER
        removing Unknown                        -> <DATASET>_rectangle_proportions.tsv
      * pickle the signature                     -> <DATASET>_rectangle_signature.pkl
  - fine (cell_type_col == "cell_type_fine"):
      * write raw estimations (all fine labels + Unknown), no reindex
                                                 -> <DATASET>_rectangle_fine_proportions.tsv
      * capture the Unknown fraction             -> <DATASET>_rectangle_fine_unknown.tsv
      * pickle the signature                     -> <DATASET>_rectangle_fine_signature.pkl

CLI:
  python 26_run_rectangle.py <ref_h5ad> <bulk_tpm_tsv> <output_dir> <dataset> \
      [cell_type_col=cell_type] [--n-cpus N] [--seed 42] [--no-optimize-cutoffs]

NEVER modifies existing files; only writes the outputs listed above.
"""

from __future__ import annotations

import argparse
import logging
import os
import pickle
import random
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import rectanglepy as rectangle

# The 16 canonical coarse cell-type strings, in the EXACT output column order
# mandated by the shared contract. Do not reorder.
CANONICAL_CELL_TYPES = [
    "Endothelial cells",
    "Hepatocytes",
    "Plasma cells",
    "T cells",
    "Cholangiocytes",
    "Fibroblasts",
    "Macrophages",
    "Circulating NK/NKT",
    "Resident NK",
    "Mono+mono derived cells",
    "Basophils",
    "B cells",
    "cDC1s",
    "cDC2s",
    "pDCs",
    "Neutrophils",
]

UNKNOWN_COL = "Unknown"  # exact estimations column name emitted by rectanglepy

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("rectangle_runner")


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Rectangle deconvolution of a bulk TPM matrix against a scRNA reference.",
    )
    p.add_argument("ref_h5ad", help="Single-cell reference AnnData (.h5ad); .X raw counts, var=symbols.")
    p.add_argument("bulk_tpm_tsv", help="Bulk TPM tsv: rows=SYMBOL (header 'SYMBOL'), cols=sample ids, TPM floats.")
    p.add_argument("output_dir", help="Directory to write proportion/unknown/signature outputs into.")
    p.add_argument("dataset", help="Dataset id used as the output filename prefix.")
    p.add_argument(
        "cell_type_col",
        nargs="?",
        default="cell_type",
        help="obs column holding cell-type labels. 'cell_type' => coarse mode; "
        "'cell_type_fine' => fine mode. Default: cell_type.",
    )
    p.add_argument("--n-cpus", type=int, default=None,
                   help="CPUs for the DE step. Default: SLURM_CPUS_PER_TASK if set, else all available.")
    p.add_argument("--seed", type=int, default=42, help="RNG seed (numpy/random) for reproducibility. Default 42.")
    p.add_argument("--no-optimize-cutoffs", dest="optimize_cutoffs", action="store_false",
                   help="Disable log-fold-change gridsearch (faster; uses fixed p/lfc). "
                        "Default: optimize_cutoffs=True (rectanglepy default).")
    p.set_defaults(optimize_cutoffs=True)
    return p.parse_args(argv)


def load_reference(ref_h5ad: str, cell_type_col: str) -> ad.AnnData:
    """Load the reference AnnData and validate the contract (raw counts, label col, unique symbols)."""
    log.info("Loading reference AnnData: %s", ref_h5ad)
    adata = ad.read_h5ad(ref_h5ad)
    log.info("Reference: %d cells x %d genes", adata.n_obs, adata.n_vars)

    if cell_type_col not in adata.obs.columns:
        raise KeyError(
            f"cell_type_col '{cell_type_col}' not in reference obs columns: {list(adata.obs.columns)}"
        )

    if not adata.var_names.is_unique:
        n_dup = int(adata.var_names.duplicated().sum())
        raise ValueError(
            f"Reference var_names are not unique ({n_dup} duplicate gene symbols). "
            "rectanglepy asserts uniqueness; fix the reference builder (script 25). "
            "This runner will not silently mangle gene identities."
        )

    # Validate .X is raw, non-negative, integer-like counts (DESeq2 requires counts).
    X = adata.X
    sub = X[: min(500, adata.n_obs)]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    if sub.min() < 0:
        raise ValueError("Reference .X contains negative values; expected raw counts.")
    if not np.allclose(sub, np.round(sub)):
        raise ValueError(
            "Reference .X does not look like integer raw counts (values are not integer-like). "
            "rectanglepy/DESeq2 needs raw counts in .X."
        )
    log.info("Reference .X validated as non-negative integer counts (dtype=%s).", X.dtype)

    n_types = adata.obs[cell_type_col].nunique()
    log.info("Cell-type column '%s': %d labels.", cell_type_col, n_types)
    return adata


def load_bulk(bulk_tpm_tsv: str) -> pd.DataFrame:
    """Load a genes x samples TPM tsv and return a samples x genes float DataFrame.

    Collapses duplicate gene symbols (rows) by summing TPM so that the transposed
    column index is unique, which rectanglepy asserts.
    """
    log.info("Loading bulk TPM: %s", bulk_tpm_tsv)
    df = pd.read_csv(bulk_tpm_tsv, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if not df.index.is_unique:
        n_dup = int(df.index.duplicated().sum())
        log.warning("Bulk has %d duplicate gene symbols; collapsing by summing TPM.", n_dup)
        df = df.groupby(level=0).sum()
    bulks = df.T.astype(float)  # samples x genes
    bulks.index = bulks.index.astype(str)
    assert bulks.columns.is_unique, "bulk gene columns not unique after collapse"
    log.info("Bulk: %d samples x %d genes.", bulks.shape[0], bulks.shape[1])
    return bulks


def report_gene_overlap(bulks: pd.DataFrame, adata: ad.AnnData) -> None:
    ref_genes = set(adata.var_names)
    bulk_genes = set(bulks.columns)
    overlap = bulk_genes & ref_genes
    pct_bulk = 100.0 * len(overlap) / max(1, len(bulk_genes))
    pct_ref = 100.0 * len(overlap) / max(1, len(ref_genes))
    log.info(
        "Gene overlap bulk<->reference: %d genes (%.1f%% of bulk genes, %.1f%% of ref genes).",
        len(overlap), pct_bulk, pct_ref,
    )
    if pct_bulk < 60.0 or pct_ref < 60.0:
        log.warning(
            "LOW gene overlap (<60%%): pct_bulk=%.1f%% pct_ref=%.1f%%. "
            "Deconvolution may be unreliable; check symbol harmonization.",
            pct_bulk, pct_ref,
        )


def write_coarse_outputs(est: pd.DataFrame, output_dir: str, dataset: str) -> None:
    """Coarse mode: separate Unknown, reindex to the 16 canonical types, re-normalize."""
    # 1) Raw Unknown fraction (kept as-is, before any re-normalization).
    unknown = est[UNKNOWN_COL] if UNKNOWN_COL in est.columns else pd.Series(0.0, index=est.index)
    unk_df = unknown.rename("Unknown").rename_axis("sample_id").reset_index()
    unk_path = os.path.join(output_dir, f"{dataset}_rectangle_unknown.tsv")
    unk_df.to_csv(unk_path, sep="\t", index=False)
    log.info("Wrote Unknown fractions -> %s", unk_path)

    # 2) Drop Unknown, warn if rectanglepy emitted any unexpected non-canonical type.
    known = est.drop(columns=[UNKNOWN_COL], errors="ignore")
    extra = [c for c in known.columns if c not in CANONICAL_CELL_TYPES]
    if extra:
        log.warning("Dropping %d non-canonical cell-type column(s) not in the 16: %s", len(extra), extra)

    # 3) Reindex to the 16 canonical types IN ORDER (0-fill any absent type).
    prop = known.reindex(columns=CANONICAL_CELL_TYPES, fill_value=0.0)

    # 4) RE-NORMALIZE the 16 to sum to 1 AFTER removing Unknown. Rows whose known
    #    mass is 0 (all Unknown) are left as zeros to avoid divide-by-zero.
    row_sums = prop.sum(axis=1)
    nonzero = row_sums > 0
    prop.loc[nonzero] = prop.loc[nonzero].div(row_sums[nonzero], axis=0)
    if (~nonzero).any():
        log.warning("%d sample(s) had all mass in Unknown; left as zeros after re-normalization.",
                    int((~nonzero).sum()))

    prop.index = prop.index.astype(str)
    prop.index.name = "sample_id"
    prop_path = os.path.join(output_dir, f"{dataset}_rectangle_proportions.tsv")
    # index_label=False writes the R write.table(row.names=TRUE) "headerless-corner"
    # format (header = N cell-type fields, data rows = N+1 fields) so this file is
    # byte-compatible with the existing *_music_prop_weighted.tsv / *_bayesprism_proportions.tsv
    # and is read IDENTICALLY by read.table(header=TRUE) in script 25 / 16_compare_methods.
    prop.to_csv(prop_path, sep="\t", index=True, index_label=False)
    log.info("Wrote 16-type re-normalized proportions -> %s  (rows=%d, cols=%d)",
             prop_path, prop.shape[0], prop.shape[1])


def write_fine_outputs(est: pd.DataFrame, output_dir: str, dataset: str) -> None:
    """Fine mode: write raw estimations (all fine labels + Unknown), plus the Unknown fraction."""
    fine = est.copy()
    fine.index = fine.index.astype(str)
    fine.index.name = "sample_id"
    fine_path = os.path.join(output_dir, f"{dataset}_rectangle_fine_proportions.tsv")
    # index_label=False -> R headerless-corner format (matches the MuSiC/BayesPrism convention).
    fine.to_csv(fine_path, sep="\t", index=True, index_label=False)
    log.info("Wrote fine proportions (all labels + Unknown) -> %s  (rows=%d, cols=%d)",
             fine_path, fine.shape[0], fine.shape[1])

    unknown = est[UNKNOWN_COL] if UNKNOWN_COL in est.columns else pd.Series(0.0, index=est.index)
    unk_df = unknown.rename("Unknown").rename_axis("sample_id").reset_index()
    unk_path = os.path.join(output_dir, f"{dataset}_rectangle_fine_unknown.tsv")
    unk_df.to_csv(unk_path, sep="\t", index=False)
    log.info("Wrote fine Unknown fractions -> %s", unk_path)


def main(argv=None) -> int:
    args = parse_args(argv)
    t0 = time.time()

    # Reproducibility: rectanglepy exposes no seed kwarg, so seed the global RNGs
    # it may draw from (bootstrap/consensus). PYTHONHASHSEED must be set in the
    # environment before interpreter start (see the sbatch wrapper).
    random.seed(args.seed)
    np.random.seed(args.seed)

    fine_mode = args.cell_type_col == "cell_type_fine"
    mode = "FINE" if fine_mode else "COARSE"

    n_cpus = args.n_cpus
    if n_cpus is None and os.environ.get("SLURM_CPUS_PER_TASK"):
        n_cpus = int(os.environ["SLURM_CPUS_PER_TASK"])

    log.info("=== Rectangle runner: dataset=%s mode=%s cell_type_col=%s ===",
             args.dataset, mode, args.cell_type_col)
    log.info("n_cpus=%s optimize_cutoffs=%s seed=%d", n_cpus, args.optimize_cutoffs, args.seed)

    os.makedirs(args.output_dir, exist_ok=True)

    adata = load_reference(args.ref_h5ad, args.cell_type_col)
    bulks = load_bulk(args.bulk_tpm_tsv)
    report_gene_overlap(bulks, adata)

    log.info("Running rectanglepy.rectangle() ...")
    t_run = time.time()
    est, sig = rectangle.rectangle(
        adata,
        bulks,
        cell_type_col=args.cell_type_col,
        optimize_cutoffs=args.optimize_cutoffs,
        n_cpus=n_cpus,
    )
    log.info("Deconvolution finished in %.1f s. estimations shape=%s, columns=%s",
             time.time() - t_run, est.shape, list(est.columns))

    # Sanity: rows should sum to ~1 (including Unknown).
    rs = est.sum(axis=1)
    log.info("estimations row-sum range: [%.4f, %.4f] (expect ~1.0 incl Unknown).",
             float(rs.min()), float(rs.max()))

    if fine_mode:
        write_fine_outputs(est, args.output_dir, args.dataset)
        sig_name = f"{args.dataset}_rectangle_fine_signature.pkl"
    else:
        write_coarse_outputs(est, args.output_dir, args.dataset)
        sig_name = f"{args.dataset}_rectangle_signature.pkl"

    sig_path = os.path.join(args.output_dir, sig_name)
    try:
        with open(sig_path, "wb") as fh:
            pickle.dump(sig, fh)
        log.info("Wrote signature object -> %s", sig_path)
    except Exception as exc:  # pragma: no cover - defensive
        log.warning("Could not pickle signature object (%s): %s", type(sig).__name__, exc)

    log.info("=== DONE dataset=%s mode=%s total=%.1f s ===",
             args.dataset, mode, time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
