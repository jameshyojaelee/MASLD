#!/usr/bin/env python3
"""H24 step 3: GSE244832 donor-by-lineage RNA pseudobulk for the 18 donors with
a paired snATAC sample (data/GSE244832/metadata/donor_pairing.csv).

RNA and ATAC come from different aliquots of the same donor: donor-matched, not
nucleus-matched. RNA counts are raw UMIs from each run's Cell Ranger
(refdata-gex-GRCh38-2024-A) filtered matrix; nuclei are labelled with the
atlas `cell_type` (scalesc_human_annotated_celltypist.h5ad obs; index
`<SRR>_<barcode>`). ATAC counts are read later from the ATAC Context v3
candidate (lineage-specific 500-bp consensus peaks).

Usage: 03_prepare_gse244832.py <OUT_DIR>
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
PAIR = ROOT / "data/GSE244832/metadata/donor_pairing.csv"
CR = ROOT / "data/GSE244832/cellranger"
ATLAS = ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
LINEAGE = {"Hepatocytes": "hepatocyte", "Fibroblasts": "stellate", "Macrophages": "macrophage",
           "Mono+mono derived cells": "macrophage", "Cholangiocytes": "cholangiocyte",
           "T cells": "t_nk", "Circulating NK/NKT": "t_nk", "Resident NK": "t_nk",
           "Endothelial cells": "endothelial_cell", "Plasma cells": "plasma_b"}


def decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def main(out: Path) -> None:
    out = out / "gse244832"
    out.mkdir(parents=True, exist_ok=True)
    pairs = list(csv.DictReader(open(PAIR)))
    srr_donor = {s: p["donor_id"] for p in pairs for s in p["rna_srrs"].split(";")}

    with h5py.File(ATLAS, "r") as f:
        obs = f["obs"]
        idx = decode(obs["_index"][:])
        g = obs["dataset"]
        ds = decode(g["categories"][:])[g["codes"][:]]
        g = obs["cell_type"]
        ct = decode(g["categories"][:])[g["codes"][:]]
    m = ds == "GSE244832"
    label = dict(zip(idx[m], ct[m]))

    genes = None
    sums: dict[tuple[str, str], np.ndarray] = {}
    nuclei: Counter = Counter()
    srr_report = []
    for srr, donor in sorted(srr_donor.items()):
        h5 = CR / srr / "outs/filtered_feature_bc_matrix.h5"
        if not h5.exists():
            srr_report.append({"srr": srr, "donor": donor, "status": "missing_h5"})
            continue
        with h5py.File(h5, "r") as f:
            mat = f["matrix"]
            bc = decode(mat["barcodes"][:])
            ids = decode(mat["features/id"][:])
            ftype = decode(mat["features/feature_type"][:])
            shape = tuple(mat["shape"][:])
            x = sparse.csc_matrix((mat["data"][:], mat["indices"][:], mat["indptr"][:]), shape=shape)
        keep_g = ftype == "Gene Expression"
        if genes is None:
            genes = ids[keep_g]
        elif not np.array_equal(genes, ids[keep_g]):
            raise SystemExit(f"feature order differs in {srr}")
        lab = np.array([LINEAGE.get(label.get(f"{srr}_{b}", ""), "") for b in bc])
        n_lab = int((lab != "").sum())
        srr_report.append({"srr": srr, "donor": donor, "status": "ok", "cells_filtered": len(bc),
                           "cells_labelled": n_lab})
        for lin in sorted(set(lab) - {""}):
            cols = np.where(lab == lin)[0]
            v = np.asarray(x[keep_g][:, cols].sum(axis=1)).ravel().astype(np.float64)
            key = (donor, lin)
            sums[key] = sums.get(key, 0) + v
            nuclei[key] += len(cols)

    keys = sorted(sums)
    mat = np.column_stack([sums[k] for k in keys]).astype(np.float32)
    if mat.max() >= 2 ** 24:
        raise SystemExit("a count exceeds exact float32 range")
    np.ascontiguousarray(mat.T, dtype="<f4").tofile(out / "rna_gene_by_unit.f32")
    (out / "dims.json").write_text(json.dumps({"rna": [int(mat.shape[0]), len(keys)]}))
    with open(out / "genes.tsv", "w") as fh:
        fh.write("ensembl_id\n" + "\n".join(genes) + "\n")
    with open(out / "rna_units.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["unit_index", "donor_id", "lineage_id", "nuclei", "rna_umis"])
        for i, k in enumerate(keys):
            w.writerow([i, k[0], k[1], nuclei[k], int(mat[:, i].sum())])
    with open(out / "srr_report.tsv", "w") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["srr", "donor", "status", "cells_filtered", "cells_labelled"])
        w.writeheader()
        w.writerows(srr_report)
    rep = {"donors": len({k[0] for k in keys}), "units": len(keys), "genes": int(mat.shape[0]),
           "srr_ok": sum(r["status"] == "ok" for r in srr_report), "srr_total": len(srr_report),
           "note": "runs of one donor are summed; if runs are re-sequencing of one library, UMIs can be counted once per run"}
    (out / "prepare_report.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
