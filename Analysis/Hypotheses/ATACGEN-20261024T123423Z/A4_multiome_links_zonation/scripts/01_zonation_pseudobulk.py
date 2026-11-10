#!/usr/bin/env python3
"""A4 step 1: hepatocyte zonation score and donor x zone-bin pseudobulk (GSE296875).

Raw RNA UMIs per hepatocyte nucleus from the eight raw 10x HDF5 files; ATAC peak
counts per nucleus from the crossmodal rectangle. Score and bins as in PRESPEC.md:
  x_nm = UMI_nm / total GE UMI_n, scaled by the mean of x_m over hepatocyte nuclei;
  PC_n, PP_n = sums over the scoring markers; s_n = PC_n / (PC_n + PP_n);
  within-donor tertiles of s (ties broken by a seeded random order), nuclei with
  PC_n + PP_n = 0 unbinned. Bin 0 = PP, 1 = mid, 2 = PC.
No association with ATAC or phenotype is computed here.

Usage: 01_zonation_pseudobulk.py <A4_DIR>
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 20261024
TIE_OFFSET = 7
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
EXE = ROOT / "Analysis/MASLD_Model_Benchmark/executions"
RECT = EXE / "gse296875-crossmodal-recurrence-20260830T193000Z/rectangle"
ROSTER = EXE / "gse296875-fragment-membership-21063829/cell_membership.tsv.gz"
GENES = EXE / "gse296875-donor-lineage-pseudobulk-21107249/pseudobulk/genes.tsv"
H24 = ROOT / "Analysis/Hypotheses/H24-20261008T211750Z/gse296875"
REF = ROOT / "Analysis/SingleCell/continual_integration/reference"
LOCK = REF / "gse296875_source_lock_v44.json"
PC = ["CYP2E1", "GLUL", "CYP1A2", "OAT"]
PP = ["HAL", "SDS", "CPS1", "PCK1", "ASS1"]


def read_well(path: Path, barcodes: list[str], gene_ids: np.ndarray) -> sparse.csc_matrix:
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        shape = tuple(g["shape"][()])
        ftype = g["features/feature_type"][()]
        fid = np.array([x.decode() for x in g["features/id"][()]])
        ge = np.where(ftype == b"Gene Expression")[0]
        if not np.array_equal(fid[ge], gene_ids):
            raise SystemExit(f"{path.name}: gene axis differs from genes.tsv")
        bc = pd.Index([x.decode() for x in g["barcodes"][()]])
        cols = bc.get_indexer(barcodes)
        if (cols < 0).any():
            raise SystemExit(f"{path.name}: {(cols < 0).sum()} roster barcodes missing")
        m = sparse.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), shape=shape)
    m = m[:, cols]
    return m[ge, :].tocsc()


def main(a4: Path) -> None:
    out = a4 / "work" / "zone"
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(GENES, sep="\t")
    gene_ids = genes["ensembl_id"].to_numpy()
    sym = genes["symbol"].to_numpy()
    for m in PC + PP:
        if (sym == m).sum() != 1:
            raise SystemExit(f"marker {m} is not a unique symbol")
    roster = pd.read_csv(ROSTER, sep="\t", dtype=str)
    hep = roster[roster["lineage_id"] == "hepatocyte"].reset_index(drop=True)
    print("hepatocyte nuclei:", len(hep), "donors:", hep["donor_id"].nunique(), flush=True)

    lock = json.loads(LOCK.read_text())["raw_rna_sources"]
    blocks, order = [], []
    for well, src in sorted(lock.items()):
        sub = hep[hep["well_id"] == well]
        if not len(sub):
            continue
        blocks.append(read_well((REF / src["path"]).resolve(), sub["raw_barcode"].tolist(), gene_ids))
        order.append(sub.index.to_numpy())
        print(well, "nuclei", len(sub), flush=True)
    rna = sparse.hstack(blocks).tocsc()
    perm = np.argsort(np.concatenate(order))
    rna = rna[:, perm].tocsc()  # genes x hepatocyte nuclei, roster order
    if rna.shape[1] != len(hep):
        raise SystemExit("RNA nuclei count mismatch")

    # Check against the H24 hepatocyte donor pseudobulk (same raw UMIs).
    units = pd.read_csv(H24 / "units.tsv", sep="\t", dtype={"donor_id": str})
    dims = json.loads((H24 / "dims.json").read_text())
    h24_rna = np.fromfile(H24 / "rna_gene_by_unit.f32", dtype="<f4").reshape(dims["rna"][1], dims["rna"][0]).T
    h24_atac = np.memmap(H24 / "atac_peak_by_unit.f32", dtype="<f4", mode="r",
                         shape=(dims["atac"][1], dims["atac"][0]))
    hu = units[units["lineage_id"] == "hepatocyte"].set_index("donor_id")
    donors = sorted(hep["donor_id"].unique(), key=int)
    d_idx = pd.Index(donors).get_indexer(hep["donor_id"])
    agg_d = sparse.csr_matrix((np.ones(len(hep)), (np.arange(len(hep)), d_idx)), shape=(len(hep), len(donors)))
    rna_d = (rna @ agg_d).toarray()
    ref = h24_rna[:, hu.loc[donors, "unit_index"].to_numpy()]
    if not np.array_equal(rna_d, ref):
        raise SystemExit("hepatocyte RNA donor sums do not reproduce the H24 pseudobulk")
    print("RNA donor sums reproduce H24 exactly", flush=True)

    total = np.asarray(rna.sum(axis=0)).ravel()
    def frac(ms):
        rows = [int(np.where(sym == m)[0][0]) for m in ms]
        x = rna[rows, :].toarray() / np.maximum(total, 1)
        return x / x.mean(axis=1, keepdims=True)
    pc = frac(PC).sum(axis=0)
    pp = frac(PP).sum(axis=0)
    den = pc + pp
    s = np.where(den > 0, pc / np.where(den > 0, den, 1), np.nan)

    rng = np.random.default_rng(SEED + TIE_OFFSET)
    tie = rng.random(len(hep))
    zbin = np.full(len(hep), -1, dtype=np.int64)
    for d in range(len(donors)):
        idx = np.where((d_idx == d) & np.isfinite(s))[0]
        o = idx[np.lexsort((tie[idx], s[idx]))]
        zbin[o] = (3 * np.arange(len(o))) // len(o)

    # ATAC per nucleus
    nuclei = [l.strip() for l in open(RECT / "nuclei.txt")]
    col = pd.Index(nuclei).get_indexer(hep["cell_id"])
    if (col < 0).any():
        raise SystemExit("hepatocyte nuclei missing from the ATAC rectangle")
    data = np.fromfile(RECT / "atac.data.i32", dtype=np.int32)
    ind = np.fromfile(RECT / "atac.indices.i32", dtype=np.int32)
    ptr = np.fromfile(RECT / "atac.indptr.i32", dtype=np.int32).astype(np.int64)
    n_peak = dims["atac"][0]
    atac = sparse.csc_matrix((data.astype(np.float32), ind, ptr), shape=(n_peak, len(nuclei)))
    del data, ind
    atac = atac[:, col].tocsc()
    atac_tot = np.asarray(atac.sum(axis=0)).ravel()
    atac_d = (atac @ agg_d).toarray()
    ref_a = np.asarray(h24_atac[hu.loc[donors, "unit_index"].to_numpy(), :]).T
    if not np.allclose(atac_d, ref_a):
        raise SystemExit("hepatocyte ATAC donor sums do not reproduce the H24 pseudobulk")
    print("ATAC donor sums reproduce H24", flush=True)

    keep = zbin >= 0
    u_key = d_idx[keep] * 3 + zbin[keep]
    n_unit = len(donors) * 3
    agg_u = sparse.csr_matrix((np.ones(keep.sum()), (np.arange(keep.sum()), u_key)), shape=(keep.sum(), n_unit))
    atac_u = (atac[:, np.where(keep)[0]] @ agg_u).toarray()
    rna_u = (rna[:, np.where(keep)[0]] @ agg_u).toarray()
    if max(atac_u.max(), rna_u.max()) >= 2 ** 24:
        raise SystemExit("a count exceeds exact float32 range")
    np.ascontiguousarray(atac_u.T, dtype="<f4").tofile(out / "atac_peak_by_unit.f32")
    np.ascontiguousarray(rna_u.T, dtype="<f4").tofile(out / "rna_gene_by_unit.f32")
    (out / "dims.json").write_text(json.dumps({"atac": [n_peak, n_unit], "rna": [int(rna.shape[0]), n_unit]}))
    well_of = hep.groupby("donor_id")["well_id"].first()
    nuc_u = np.bincount(u_key, minlength=n_unit)
    with open(out / "units.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["unit_index", "donor_id", "well_id", "zone_bin", "zone", "nuclei",
                    "atac_in_peak_fragments", "rna_umis"])
        for d, donor in enumerate(donors):
            for b, name in enumerate(["PP", "mid", "PC"]):
                i = d * 3 + b
                w.writerow([i, donor, well_of[donor], b, name, int(nuc_u[i]),
                            int(atac_u[:, i].sum()), int(rna_u[:, i].sum())])
    sc = pd.DataFrame({"cell_id": hep["cell_id"], "donor_id": hep["donor_id"], "well_id": hep["well_id"],
                       "rna_total_umi": total.astype(np.int64), "atac_in_peak": atac_tot.astype(np.int64),
                       "pc_sum": pc, "pp_sum": pp, "score": s, "zone_bin": zbin})
    with gzip.open(out / "nuclei_scores.tsv.gz", "wt") as fh:
        sc.to_csv(fh, sep="\t", index=False)
    rep = {"hepatocyte_nuclei": int(len(hep)), "binned": int(keep.sum()), "unbinned": int((~keep).sum()),
           "donors": len(donors), "units": n_unit, "min_nuclei_per_unit": int(nuc_u.min()),
           "rna_donor_sums_reproduce_H24": True, "atac_donor_sums_reproduce_H24": True}
    (out / "prepare_report.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
