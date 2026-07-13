#!/usr/bin/env python
"""
GLP-1RA decontX export: subset the integrated human atlas to Hepatocytes +
Endothelial cells, split by dataset, and write raw integer counts (genes x cells,
CSC) as compact binary components that the decontx R env (only `Matrix` available,
no HDF5 reader) reconstructs into a dgCMatrix via readBin().

Why Hep+Endo TOGETHER, per dataset:
  decontX estimates ONE ambient profile per run. Keeping hepatocytes (ALB-high,
  ~85% of the pool) with endothelial cells makes the ambient pool hepatocyte-
  flavored, so any ALB / hepatocyte bleed into endothelial cells is modeled as
  contamination. Splitting by dataset handles batch (different ambient pools /
  chemistries) by looping in R.

Transport: per dataset, column-chunked CSC binary (chunk size <= 250k cells so
each chunk's nnz stays < 2^31, R readBin limit), plus a per-cell meta table and a
checksum for the axis genes so the R side can VERIFY reconstruction.

Raw counts live in `.raw.X` (verified integer; `.X` is lognorm).
"""
import os
import json
import numpy as np
import scipy.sparse as sp
import anndata as ad

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
H5 = os.path.join(ROOT, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")
WORK = os.path.join(ROOT, "RNA-seq/results/glp1ra/scrna_incretin_axis/decontx_work")
os.makedirs(WORK, exist_ok=True)

CELLTYPES = ["Hepatocytes", "Endothelial cells"]
# z for decontX = cell_type mapped to consecutive integers (conservative, well
# populated, directly the biological unit). leiden also exported for sensitivity.
Z_MAP = {"Hepatocytes": 1, "Endothelial cells": 2}
CHUNK = 250_000  # cells per chunk -> keeps per-chunk nnz < 2^31

# genes whose reconstruction we verify + report on downstream
AXIS = ["GLP1R", "GIPR", "GCGR", "DPP4", "GCG", "GLP2R"]
CTRL = ["ALB", "PECAM1", "STAB2"]
CHECK_GENES = AXIS + CTRL

# datasets with < MIN_ENDO endothelial or < MIN_HEP hepatocytes are degenerate
# for the Hep-vs-Endo ambient test and are skipped (reported in manifest).
MIN_ENDO = 50
MIN_HEP = 50


def main():
    print(f"[export] reading {H5} (backed)")
    a = ad.read_h5ad(H5, backed="r")
    assert a.raw is not None, "raw counts required"
    var_names = np.asarray(a.raw.var_names if a.raw.var_names is not None else a.var_names)
    ngenes = len(var_names)
    # write gene order once (full transcriptome; decontX uses all genes for ambient)
    with open(os.path.join(WORK, "genes.txt"), "w") as fh:
        fh.write("\n".join(map(str, var_names)) + "\n")
    gene_pos = {g: i for i, g in enumerate(var_names)}
    check_idx = {g: gene_pos[g] for g in CHECK_GENES if g in gene_pos}
    missing = [g for g in CHECK_GENES if g not in gene_pos]
    print(f"[export] ngenes={ngenes}  check-genes present={list(check_idx)}  missing={missing}")

    obs = a.obs
    ct = obs["cell_type"].astype(str).values
    ds = obs["dataset"].astype(str).values
    leiden = obs["leiden"].astype(str).values
    sample = obs["sample"].astype(str).values
    cond = obs["condition"].astype(str).values

    in_scope = np.isin(ct, CELLTYPES)
    manifest = {"datasets": {}, "skipped": {}, "ngenes": int(ngenes),
                "check_genes": list(check_idx), "chunk": CHUNK}

    only = os.environ.get("ONLY_DATASETS", "").strip()
    only_set = set(x for x in only.split(",") if x) if only else None
    if only_set:
        print(f"[export] ONLY_DATASETS filter active: {sorted(only_set)}")

    for dataset in sorted(set(ds[in_scope])):
        if only_set is not None and dataset not in only_set:
            continue
        sel = in_scope & (ds == dataset)
        idx = np.where(sel)[0]
        ct_sel = ct[idx]
        n_hep = int((ct_sel == "Hepatocytes").sum())
        n_endo = int((ct_sel == "Endothelial cells").sum())
        if n_hep < MIN_HEP or n_endo < MIN_ENDO:
            print(f"[export] SKIP {dataset}: Hep={n_hep} Endo={n_endo} (degenerate)")
            manifest["skipped"][dataset] = {"n_hep": n_hep, "n_endo": n_endo}
            continue
        dsdir = os.path.join(WORK, dataset)
        os.makedirs(dsdir, exist_ok=True)
        ncells = len(idx)
        print(f"[export] {dataset}: {ncells} cells (Hep={n_hep} Endo={n_endo}) "
              f"-> {int(np.ceil(ncells/CHUNK))} chunk(s)")

        # per-cell meta (column order == reconstruction order)
        import pandas as pd
        meta = pd.DataFrame({
            "barcode": obs.index.values[idx],
            "cell_type": ct_sel,
            "leiden": leiden[idx],
            "sample": sample[idx],
            "condition": cond[idx],
            "z": [Z_MAP[c] for c in ct_sel],
        })
        meta.to_csv(os.path.join(dsdir, "meta.csv.gz"), index=False)

        chunks = []
        check_rowsum = {g: 0.0 for g in check_idx}
        total_sum = 0.0
        for ci, start in enumerate(range(0, ncells, CHUNK)):
            cidx = idx[start:start + CHUNK]
            Xsub = a.raw.X[cidx, :]                 # csr, cells x genes
            if not sp.isspmatrix_csr(Xsub):
                Xsub = sp.csr_matrix(Xsub)
            Xgc = Xsub.T.tocsc()                    # genes x cells
            Xgc.sort_indices()
            Xgc.eliminate_zeros()
            nnz = int(Xgc.nnz)
            assert nnz < 2**31, f"chunk nnz {nnz} exceeds int32 (raise CHUNK split)"
            assert Xgc.shape[0] == ngenes
            Xgc.data.astype("<f8").tofile(os.path.join(dsdir, f"chunk{ci}_data.bin"))
            Xgc.indices.astype("<i4").tofile(os.path.join(dsdir, f"chunk{ci}_indices.bin"))
            Xgc.indptr.astype("<i4").tofile(os.path.join(dsdir, f"chunk{ci}_indptr.bin"))
            chunks.append({"chunk": ci, "ncells": int(Xgc.shape[1]), "nnz": nnz})
            total_sum += float(Xgc.data.sum())
            for g, gi in check_idx.items():
                check_rowsum[g] += float(Xsub[:, gi].sum())
            del Xsub, Xgc

        info = {"ncells": int(ncells), "ngenes": int(ngenes),
                "n_hep": n_hep, "n_endo": n_endo, "chunks": chunks,
                "total_count_sum": total_sum,
                "check_rowsum": {g: float(v) for g, v in check_rowsum.items()}}
        with open(os.path.join(dsdir, "dims.json"), "w") as fh:
            json.dump(info, fh, indent=2)
        manifest["datasets"][dataset] = {"ncells": ncells, "n_hep": n_hep,
                                         "n_endo": n_endo, "n_chunks": len(chunks)}
        print(f"[export]   done {dataset}: total_count_sum={total_sum:.0f}")

    with open(os.path.join(WORK, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[export] manifest written. datasets={list(manifest['datasets'])} "
          f"skipped={list(manifest['skipped'])}")


if __name__ == "__main__":
    main()
