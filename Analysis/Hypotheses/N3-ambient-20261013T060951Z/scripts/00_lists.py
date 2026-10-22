#!/usr/bin/env python
"""N3 ambient check, step 0: library list, gene axis and report genes. Idempotent; never overwrites.

Usage: python -I 00_lists.py <REPO> <N3_DIR> <RESULT_DIR>
"""
import subprocess
import sys
from pathlib import Path

import h5py
import pandas as pd

REPO, N3, RES = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
CONTROLS = ["GNMT", "MAT1A", "CYP2C19", "ALB", "APOA1"]
IMMUNE = ["PTPRC", "CD3E", "CD2", "CD247", "MS4A1", "CD74"]
STROMAL = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"]

if not (RES / "genes.txt").exists():
    with h5py.File(REPO / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad", "r") as f:
        genes = [x.decode() for x in f["raw/var/_index"][:]]
    assert len(set(genes)) == len(genes)
    (RES / "genes.txt").write_text("\n".join(genes) + "\n")
genes = (RES / "genes.txt").read_text().split("\n")[:-1]
if not (RES / "report_genes.txt").exists():
    tl = pd.read_csv(N3 / "target_list.tsv", sep="\t")
    rep = list(dict.fromkeys(list(tl.symbol) + CONTROLS + IMMUNE + STROMAL))
    missing = [g for g in rep if g not in set(genes)]
    assert not missing, missing
    (RES / "report_genes.txt").write_text("\n".join(rep) + "\n")
if not (RES / "library_list.tsv").exists():
    dm = pd.read_csv(N3 / "donor_meta.tsv", sep="\t")
    lm = pd.read_csv(N3 / "library_meta.tsv", sep="\t")
    L = lm[lm.donor.isin(dm.loc[dm.included, "donor"])].sort_values(["n_cells_all", "library"])
    assert len(L) == 208 and L.donor.nunique() == 70
    L[["library", "dataset"]].to_csv(RES / "library_list.tsv", sep=" ", index=False, header=False)
if not (RES / "pip_freeze.txt").exists():
    (RES / "pip_freeze.txt").write_text(subprocess.run([sys.executable, "-m", "pip", "freeze"],
                                                       capture_output=True, text=True).stdout)
print("lists ready", flush=True)
