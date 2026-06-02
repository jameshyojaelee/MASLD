#!/usr/bin/env python3
"""
Phase 0.5 — Step 3 / v2 atlas smoke test
========================================

Quick sanity check of the v2 atlas vs v1:
  - Total cell count (alert if v2 drops >30 % of v1 cells)
  - Per-dataset breakdown
  - Hepatocyte subset count

Outputs:
    Analysis/SingleCell/results_gpu_v2_phase05/atlas_v2_smoke.txt
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc
import anndata as ad

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SC_DIR = BASE / "Analysis" / "SingleCell"
V1 = SC_DIR / "integration/output/human/scalesc_human_annotated_celltypist.h5ad"
V2_FULL = SC_DIR / "results_gpu_v2_phase05" / "atlas" / "scalesc_human_annotated_celltypist_v2.h5ad"
V2_HEP = SC_DIR / "results_gpu_v2_phase05" / "atlas" / "hepatocyte_atlas_v2.h5ad"
AUDIT = SC_DIR / "results_gpu_v2_phase05" / "atlas" / "qc_cell_audit_v2.tsv"
OUT = SC_DIR / "results_gpu_v2_phase05" / "atlas_v2_smoke.txt"


def main() -> int:
    lines: list[str] = []
    lines.append("Phase 0.5 v2 atlas smoke test")
    lines.append("=" * 60)

    # v1
    log.info("Reading v1 backed ...")
    a1 = sc.read_h5ad(V1, backed="r")
    n1 = a1.n_obs
    ds1 = a1.obs["dataset"].astype(str).value_counts()
    hep1 = int((a1.obs["cell_type"].astype(str) == "Hepatocytes").sum())
    a1.file.close()
    lines.append(f"v1 atlas: {n1:,} cells; hepatocytes: {hep1:,}")
    lines.append("v1 per-dataset:\n" + ds1.to_string())
    lines.append("")

    # v2
    if not V2_FULL.exists():
        lines.append(f"!! v2 atlas missing at {V2_FULL}")
        log.error("v2 atlas missing — writing smoke anyway")
        OUT.parent.mkdir(parents=True, exist_ok=True)
        with open(OUT, "w") as f:
            f.write("\n".join(lines))
        return 1

    log.info("Reading v2 backed ...")
    a2 = sc.read_h5ad(V2_FULL, backed="r")
    n2 = a2.n_obs
    ds2 = a2.obs["dataset"].astype(str).value_counts()
    if "cell_type" in a2.obs.columns:
        hep2 = int((a2.obs["cell_type"].astype(str) == "Hepatocytes").sum())
    else:
        hep2 = -1
    a2.file.close()

    lines.append(f"v2 atlas: {n2:,} cells; hepatocytes: {hep2:,}")
    lines.append("v2 per-dataset:\n" + ds2.to_string())
    lines.append("")

    pct_kept = (n2 / n1) * 100.0 if n1 > 0 else 0.0
    lines.append(f"Total retention (v2/v1): {pct_kept:.1f}%")
    if pct_kept < 70.0:
        lines.append("**ALERT**: v2 dropped >30% of v1 cells — investigate per-dataset audit.")

    pct_kept_hep = (hep2 / hep1) * 100.0 if hep1 > 0 else 0.0
    lines.append(f"Hepatocyte retention (v2/v1): {pct_kept_hep:.1f}%")
    if hep1 > 0 and pct_kept_hep < 70.0:
        lines.append("**ALERT**: hepatocyte fraction dropped >30%")

    # Cross-check with audit TSV
    if AUDIT.exists():
        lines.append("")
        lines.append("Per-dataset audit (from qc_cell_audit_v2.tsv):")
        adf = pd.read_csv(AUDIT, sep="\t")
        lines.append(adf.to_string(index=False))

    # Hepatocyte subset file check
    lines.append("")
    if V2_HEP.exists():
        log.info("Reading hepatocyte subset ...")
        ah = sc.read_h5ad(V2_HEP, backed="r")
        lines.append(f"hepatocyte_atlas_v2.h5ad: {ah.shape}")
        ah.file.close()
    else:
        lines.append(f"!! Hepatocyte subset h5ad missing at {V2_HEP}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    log.info("Wrote %s", OUT)
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
