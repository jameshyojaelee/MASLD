#!/usr/bin/env python
"""
magma_to_gs.py
==============
Convert MAGMA .genes.out files (per-gene z-scores) to a single scDRS .gs
file with multiple traits. Following Zhang 2022 scDRS convention: top 1000
genes per trait by ZSTAT, gene weights = ZSTAT.

Inputs:
  tools/magma/work/gene_analysis_out/{gwas_id}.genes.out

Outputs:
  Analysis/SingleCell/results_gpu_v2/disease_signatures/MAGMA_anchored_scdrs.gs
  Analysis/SingleCell/results_gpu_v2/disease_signatures/MAGMA_anchored_scdrs_summary.csv
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
GENES_OUT_DIR = BASE / "tools/magma/work/gene_analysis_out"
GENE_LOC_PATH = BASE / "tools/magma/ref/NCBI37.3.gene.loc"
N_TABLE_PATH  = BASE / "tools/magma/work/gwas_inputs/_N_table.tsv"
SIG_DIR       = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
GS_OUT        = SIG_DIR / "MAGMA_anchored_scdrs.gs"
SUMMARY_OUT   = SIG_DIR / "MAGMA_anchored_scdrs_summary.csv"

TOP_N = 1000  # Zhang 2022 default top-N gene-set size


def log(msg: str) -> None:
    print(f"[magma_to_gs] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def load_gene_id_to_symbol(loc_path: Path) -> dict:
    """NCBI37.3.gene.loc cols: entrez_id, chr, start, end, strand, symbol."""
    df = pd.read_csv(
        loc_path, sep="\t", header=None,
        names=["entrez_id", "chr", "start", "end", "strand", "symbol"],
        dtype=str,
    )
    return dict(zip(df["entrez_id"].astype(str), df["symbol"].astype(str)))


def main() -> None:
    log(f"loading gene_id → symbol map from {GENE_LOC_PATH}")
    id_to_sym = load_gene_id_to_symbol(GENE_LOC_PATH)
    log(f"  {len(id_to_sym):,} entrez → symbol entries")

    log(f"loading audit table {N_TABLE_PATH}")
    audit = pd.read_csv(N_TABLE_PATH, sep="\t")
    audit = audit[audit["status"] == "OK"]
    log(f"  {len(audit)} GWAS with OK status")

    gs_lines = ["TRAIT\tGENESET"]
    summary_rows = []
    for _, row in audit.iterrows():
        gwas_id = row["gwas_id"]
        genes_out = GENES_OUT_DIR / f"{gwas_id}.genes.out"
        if not genes_out.exists():
            log(f"  WARNING: missing {genes_out}")
            summary_rows.append({"gwas_id": gwas_id, "family": row["family"],
                                 "n_genes_in_magma": 0, "n_top": 0,
                                 "status": "missing_genes_out"})
            continue
        log(f"reading {genes_out.name}")
        # MAGMA --gene-analysis output: header starts with "GENE", whitespace-delimited
        df = pd.read_csv(genes_out, sep=r"\s+", comment="#")
        # Expected cols include GENE, CHR, START, STOP, NSNPS, NPARAM, N, ZSTAT, P
        if "ZSTAT" not in df.columns or "GENE" not in df.columns:
            log(f"  WARNING: missing expected columns; got {list(df.columns)}")
            continue
        # Map entrez GENE id → symbol
        df["symbol"] = df["GENE"].astype(str).map(id_to_sym)
        df = df.dropna(subset=["symbol", "ZSTAT"])
        df = df[df["ZSTAT"] > 0]  # only positive z (canonical scDRS use)
        df = df.sort_values("ZSTAT", ascending=False).drop_duplicates("symbol")
        n_total = len(df)
        top = df.head(TOP_N)
        gs_str = ",".join([f"{s}:{z:.6f}" for s, z in zip(top["symbol"], top["ZSTAT"])])
        gs_lines.append(f"{gwas_id}\t{gs_str}")
        summary_rows.append({"gwas_id": gwas_id, "family": row["family"],
                             "n_genes_in_magma": int(n_total),
                             "n_top": int(len(top)),
                             "z_max": float(top["ZSTAT"].max()),
                             "z_min_top": float(top["ZSTAT"].min()),
                             "status": "OK"})
        log(f"  {gwas_id}: {n_total} genes (z>0), top {len(top)} kept "
            f"(z range {top['ZSTAT'].min():.2f}–{top['ZSTAT'].max():.2f})")

    log(f"writing {GS_OUT}")
    with open(GS_OUT, "w") as f:
        f.write("\n".join(gs_lines) + "\n")
    log(f"  wrote {len(gs_lines)-1} traits")

    pd.DataFrame(summary_rows).to_csv(SUMMARY_OUT, index=False)
    log(f"wrote {SUMMARY_OUT}")
    log("done.")


if __name__ == "__main__":
    main()
