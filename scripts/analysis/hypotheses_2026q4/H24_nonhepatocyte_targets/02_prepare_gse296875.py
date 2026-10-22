#!/usr/bin/env python3
"""H24 step 2: GSE296875 donor-by-lineage pseudobulk for RNA and ATAC.

Same-nucleus multiome, 39 donors, 7 lineages, 68,398 nuclei. Inputs:
- RNA: frozen raw-UMI donor-by-lineage pseudobulk (job 21107249), built from the
  eight raw 10x HDF5 files (not the decontX matrix).
- ATAC: the authors' `peaks` assay counts (306,706 peaks x 68,398 nuclei) dumped
  as dgCMatrix slots in the crossmodal rectangle; summed here over nuclei of each
  (donor, lineage) unit using the frozen fragment-membership roster, joined on
  (well_id, raw_barcode) through the well-namespaced cell_id.
No association between modalities is computed here.

Usage: 02_prepare_gse296875.py <OUT_DIR>
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

EXE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
RECT = EXE / "gse296875-crossmodal-recurrence-20260830T193000Z/rectangle"
ROSTER = EXE / "gse296875-fragment-membership-21063829/cell_membership.tsv.gz"
PB = EXE / "gse296875-donor-lineage-pseudobulk-21107249/pseudobulk"


def main(out: Path) -> None:
    out = out / "gse296875"
    out.mkdir(parents=True, exist_ok=True)
    nuclei = [l.strip() for l in open(RECT / "nuclei.txt")]
    peaks = [l.strip() for l in open(RECT / "atac.features.txt")]
    n_nuc, n_peak = len(nuclei), len(peaks)
    if n_nuc != 68398 or n_peak != 306706:
        raise SystemExit("rectangle shape differs from the recorded 306,706 x 68,398")

    with gzip.open(ROSTER, "rt") as fh:
        roster = {r["cell_id"]: (r["donor_id"], r["lineage_id"]) for r in csv.DictReader(fh, delimiter="\t")}
    units = list(csv.DictReader(open(PB / "units.tsv"), delimiter="\t"))
    unit_key = {(u["donor_id"], u["lineage_id"]): int(u["unit_index"]) for u in units}
    col_unit = np.array([unit_key[roster[c]] for c in nuclei], dtype=np.int64)

    data = np.fromfile(RECT / "atac.data.i32", dtype=np.int32)
    ind = np.fromfile(RECT / "atac.indices.i32", dtype=np.int32)
    ptr = np.fromfile(RECT / "atac.indptr.i32", dtype=np.int32).astype(np.int64)
    if ptr.size != n_nuc + 1 or ptr[-1] != data.size:
        raise SystemExit("ATAC slot dump inconsistent")
    atac = sparse.csc_matrix((data.astype(np.float32), ind, ptr), shape=(n_peak, n_nuc))
    del data, ind
    agg = sparse.csr_matrix((np.ones(n_nuc, dtype=np.float32), (np.arange(n_nuc), col_unit)),
                            shape=(n_nuc, len(units)))
    pb = (atac @ agg).toarray().astype(np.float64)  # peaks x units
    # Exact-integer check: every unit total equals the per-nucleus column sums.
    nuc_tot = np.asarray(atac.sum(axis=0)).ravel()
    unit_tot = np.bincount(col_unit, weights=nuc_tot, minlength=len(units))
    if not np.allclose(unit_tot, pb.sum(axis=0)):
        raise SystemExit("ATAC pseudobulk totals do not reproduce")
    nuc_per_unit = np.bincount(col_unit, minlength=len(units))

    z = np.load(PB / "donor_lineage_counts.npz")
    rna = z["counts"].T.astype(np.float64)  # genes x units
    if not np.array_equal(z["nuclei"], nuc_per_unit):
        raise SystemExit("RNA and ATAC unit nuclei counts disagree")

    if max(pb.max(), rna.max()) >= 2 ** 24:
        raise SystemExit("a count exceeds exact float32 range")
    # Column-major (feature x unit) little-endian float32 for R readBin.
    np.ascontiguousarray(pb.T, dtype="<f4").tofile(out / "atac_peak_by_unit.f32")
    np.ascontiguousarray(rna.T, dtype="<f4").tofile(out / "rna_gene_by_unit.f32")
    (out / "dims.json").write_text(json.dumps({"atac": [n_peak, len(units)], "rna": [int(rna.shape[0]), len(units)]}))
    with open(out / "peaks.tsv", "w") as fh:
        fh.write("peak_id\n" + "\n".join(peaks) + "\n")
    with open(out / "units.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["unit_index", "donor_id", "well_id", "lineage_id", "nuclei", "atac_in_peak_fragments", "rna_umis"])
        for u in units:
            i = int(u["unit_index"])
            w.writerow([i, u["donor_id"], u["well_id"], u["lineage_id"], int(nuc_per_unit[i]),
                        int(pb[:, i].sum()), int(rna[:, i].sum())])
    rep = {"units": len(units), "donors": len({u["donor_id"] for u in units}),
           "peaks": n_peak, "genes": int(rna.shape[0]), "nuclei": n_nuc,
           "atac_unit_totals_reproduced": True, "rna_atac_nuclei_agree": True}
    (out / "prepare_report.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
