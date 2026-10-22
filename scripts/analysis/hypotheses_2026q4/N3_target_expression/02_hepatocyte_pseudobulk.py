#!/usr/bin/env python
"""N3 step 2: donor-level hepatocyte pseudobulk from the integrated atlas.

Implements the data section of PRESPEC_N3.md. One streaming pass over raw UMI
counts. Writes donor_meta.tsv, library_meta.tsv, hep_donor_counts.tsv.gz,
target_celltype_summary.tsv and qc.json into the result directory.

The library-to-donor maps follow H34 (h34_y_loss.py, library_to_donor).
Usage: python -I 02_hepatocyte_pseudobulk.py <REPO_ROOT> <RESULT_DIR>
"""
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

REPO = Path(sys.argv[1])
RES = Path(sys.argv[2])
ATLAS = REPO / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
ROSTER = REPO / ("Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/"
                 "results/current_source_diagnosis_donor_roster.tsv")
LIB_SEX = REPO / "Analysis/SingleCell/metadata/snrna_donor_sex.csv"
H34_SEX = REPO / "Analysis/Hypotheses/H34-20261008T015507Z/donor_sex_calls.tsv"
FSTAGE = REPO / ("Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk/"
                 "sample_celltype_metadata.tsv")
STAGE_CODE = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}
Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D"]
MIN_HEP_CELLS = 50
CHUNK = 100_000

log = lambda *a: print(*a, flush=True)


def decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


def library_to_donor():
    m = {}
    for ds in ["GSE244832", "GSE185477", "GSE202379", "GSE136103"]:
        dp = pd.read_csv(REPO / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
        for _, r in dp.iterrows():
            for s in r["rna_srrs"].split(";"):
                assert s not in m
                m[s] = f"{ds}_{r['donor_id']}"
    runs = pd.read_csv(REPO / "Liver_Atlas/metadata/SraRunTable.csv", dtype=str)
    smp = pd.read_csv(REPO / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv", sep="\t", dtype=str)
    j = runs[["Run", "shortfilename"]].merge(
        smp[["characteristics: shortFileName", "title"]].rename(
            columns={"characteristics: shortFileName": "shortfilename"}), on="shortfilename")
    j["donor_short"] = j["title"].map(lambda t: (re.search(r"H[0-9]+", t) or [None])[0])
    j = j.dropna(subset=["donor_short"])
    assert j.Run.nunique() == 48 and j.donor_short.nunique() == 19, "Liver Atlas donor map drift"
    for r, d in zip(j.Run, j.donor_short):
        m[r] = f"Liver_Atlas_{d}"
    return m


def iter_chunks(X, indptr, n_genes):
    n = len(indptr) - 1
    for s in range(0, n, CHUNK):
        e = min(s + CHUNK, n)
        a, b = indptr[s], indptr[e]
        data = X["data"][a:b].astype(np.float32)
        idx = X["indices"][a:b]
        ptr = indptr[s:e + 1] - a
        yield s, e, sparse.csr_matrix((data, idx, ptr), shape=(e - s, n_genes))


def main():
    tl = pd.read_csv(RES / "target_list.tsv", sep="\t")
    f = h5py.File(ATLAS, "r")
    X = f["raw/X"]
    indptr = X["indptr"][:].astype(np.int64)
    genes = decode(f["raw/var/_index"][:])
    gene_ids = decode(f["raw/var/gene_ids"][:])
    assert (gene_ids == decode(f["var/gene_ids"][:])).all()
    n_genes = len(genes)
    gidx = {g: i for i, g in enumerate(genes)}
    obs = f["obs"]
    lib_cats = decode(obs["sample/categories"][:]); lib_codes = obs["sample/codes"][:].astype(np.int64)
    ct_cats = decode(obs["cell_type/categories"][:]); ct_codes = obs["cell_type/codes"][:].astype(np.int64)
    ds_cats = decode(obs["dataset/categories"][:]); ds_codes = obs["dataset/codes"][:]
    n_lib, n_ct = len(lib_cats), len(ct_cats)
    hep = int(np.where(ct_cats == "Hepatocytes")[0][0])
    n_cells = len(lib_codes)
    assert n_cells == len(indptr) - 1
    missing = [g for g in Y_GENES + ["XIST"] + list(tl.symbol) if g not in gidx]
    assert not missing, f"genes absent from atlas: {missing}"
    sex_cols = [gidx[g] for g in Y_GENES + ["XIST"]]
    tcols = [gidx[g] for g in tl.symbol]
    assert (gene_ids[tcols] == tl.ensembl.values).all(), "target Ensembl IDs differ from atlas gene_ids"

    hep_lib = np.zeros((n_lib, n_genes), dtype=np.float64)
    sex_lib = np.zeros((n_lib, len(sex_cols)), dtype=np.float64)
    ct_tsum = np.zeros((n_ct, len(tcols)), dtype=np.float64)
    ct_tdet = np.zeros((n_ct, len(tcols)), dtype=np.float64)
    integral_checked = False
    for s, e, M in iter_chunks(X, indptr, n_genes):
        if not integral_checked:
            assert np.all(M.data == np.round(M.data)), "raw/X is not integer-valued"
            integral_checked = True
        n = e - s
        lc, cc = lib_codes[s:e], ct_codes[s:e]
        ones = np.ones(n, dtype=np.float32)
        G_all = sparse.csr_matrix((ones, (lc, np.arange(n))), shape=(n_lib, n))
        sex_lib += (G_all @ M[:, sex_cols]).toarray()
        G_ct = sparse.csr_matrix((ones, (cc, np.arange(n))), shape=(n_ct, n))
        Mt = M[:, tcols]
        ct_tsum += (G_ct @ Mt).toarray()
        Mb = Mt.copy(); Mb.data[:] = 1.0
        ct_tdet += (G_ct @ Mb).toarray()
        hm = cc == hep
        if hm.any():
            nh = int(hm.sum())
            G_h = sparse.csr_matrix((np.ones(nh, dtype=np.float32), (lc[hm], np.arange(nh))), shape=(n_lib, nh))
            hep_lib += (G_h @ M[hm]).toarray()
        log("pass", e, "of", n_cells)

    # ---------------------------------------------------------- libraries
    n_hep_lib = np.bincount(lib_codes[ct_codes == hep], minlength=n_lib)
    n_all_lib = np.bincount(lib_codes, minlength=n_lib)
    ds_of_lib = pd.Series(ds_cats[ds_codes]).groupby(lib_codes).first().reindex(range(n_lib)).values
    l2d = library_to_donor()
    donor_of_lib = np.array([l2d.get(l, l) for l in lib_cats])
    lib_meta = pd.DataFrame({"library": lib_cats, "dataset": ds_of_lib, "donor": donor_of_lib,
                             "n_cells_all": n_all_lib, "n_hep_cells": n_hep_lib, "hep_umis": hep_lib.sum(1)})
    for k, g in enumerate(Y_GENES + ["XIST"]):
        lib_meta[g] = sex_lib[:, k]
    lib_meta.to_csv(RES / "library_meta.tsv", sep="\t", index=False)

    # ------------------------------------------------------------- donors
    ros = pd.read_csv(ROSTER, sep="\t", dtype={"donor": str})
    ros["stage_ord_check"] = ros["disease_stage_coarse"].map(STAGE_CODE)
    live = ros[~ros["exclude"].astype(str).str.lower().eq("true") &
               ~ros["exclude_stage_analysis"].astype(str).str.lower().eq("true")]
    # The roster leaves stage_ordinal blank for Cirrhosis; the 0-3 coding comes from disease_stage_coarse.
    lv = live[live["stage_ordinal"].notna()]
    assert (lv["stage_ordinal"].astype(float) == lv["stage_ord_check"].astype(float)).all()
    assert live["stage_ord_check"].notna().all()
    donors = sorted(set(donor_of_lib))
    d_idx = {d: i for i, d in enumerate(donors)}
    A = sparse.csr_matrix((np.ones(n_lib), (np.arange(n_lib), [d_idx[d] for d in donor_of_lib])),
                          shape=(n_lib, len(donors)))
    don_hep = (A.T @ sparse.csr_matrix(hep_lib)).toarray()          # donors x genes
    don_sex = A.T @ sex_lib                                           # donors x (4 Y + XIST)
    dm = pd.DataFrame({"donor": donors})
    dm["n_libraries"] = np.asarray(A.sum(0)).ravel().astype(int)
    dm["dataset"] = [lib_meta.loc[lib_meta.donor == d, "dataset"].iloc[0] for d in donors]
    dm["n_hep_cells"] = (A.T @ n_hep_lib).astype(int)
    dm["hep_umis"] = don_hep.sum(1)
    dm["Y"] = np.asarray(don_sex[:, :len(Y_GENES)].sum(1)).ravel()
    dm["XIST"] = np.asarray(don_sex[:, len(Y_GENES)]).ravel()
    dm["log2_Y_over_XIST"] = np.log2((dm.Y + 0.5) / (dm.XIST + 0.5))
    dm["sex_call"] = np.where(dm.Y + dm.XIST < 50, "uncallable",
                              np.where(dm.log2_Y_over_XIST > 2, "Male",
                                       np.where(dm.log2_Y_over_XIST < -2, "Female", "ambiguous")))
    dm = dm.merge(ros[["donor", "disease_stage_coarse", "exclude", "exclude_stage_analysis"]], on="donor", how="left")
    dm["in_roster"] = dm["disease_stage_coarse"].notna()
    dm["stage_ord"] = dm["disease_stage_coarse"].map(STAGE_CODE)
    h34 = pd.read_csv(H34_SEX, sep="\t")[["donor", "sex"]].rename(columns={"sex": "h34_sex"})
    dm = dm.merge(h34, on="donor", how="left")
    ls = pd.read_csv(LIB_SEX, dtype=str)[["sample", "inferred_sex"]]
    ls["donor"] = ls["sample"].map(lambda s: l2d.get(s, s))
    rs = ls.groupby("donor")["inferred_sex"].agg(lambda x: x.iloc[0] if x.nunique() == 1 else "mixed").rename("repo_library_sex")
    dm = dm.merge(rs, on="donor", how="left")
    fs = pd.read_csv(FSTAGE, sep="\t", dtype=str)
    fs = fs[(fs.cell_type == "Hepatocytes") & fs.F_stage_documented.notna() & (fs.F_stage_documented != "")].copy()
    fs["donor"] = fs["sample"].map(lambda s: l2d.get(s, s))
    fd = fs.groupby("donor")["F_stage_documented"].agg(lambda x: float(x.iloc[0]) if x.nunique() == 1 else np.nan).rename("F_documented")
    dm = dm.merge(fd, on="donor", how="left")
    ex = dm["exclude"].astype(str).str.lower().eq("true") | dm["exclude_stage_analysis"].astype(str).str.lower().eq("true")
    reason = np.where(~dm.in_roster, "not_in_roster",
              np.where(ex, "roster_excluded",
               np.where(dm.n_hep_cells < MIN_HEP_CELLS, f"fewer_than_{MIN_HEP_CELLS}_hepatocytes",
                np.where(~dm.sex_call.isin(["Male", "Female"]), "no_sex_call", "included"))))
    dm["inclusion"] = reason
    dm["included"] = dm.inclusion == "included"
    dm.to_csv(RES / "donor_meta.tsv", sep="\t", index=False)

    keep = dm.in_roster & ~ex
    cols = [d for d, k in zip(donors, keep) if k]
    out = pd.DataFrame(don_hep[keep.values].T.astype(np.int64), columns=cols)
    out.insert(0, "ensembl", gene_ids); out.insert(0, "symbol", genes)
    out.to_csv(RES / "hep_donor_counts.tsv.gz", sep="\t", index=False)

    # ------------------------------------------------ target x cell type
    rows = []
    for i, c in enumerate(ct_cats):
        for j, g in enumerate(tl.symbol):
            rows.append((c, g, ct_tsum[i, j], ct_tdet[i, j], int((ct_codes == i).sum())))
    pd.DataFrame(rows, columns=["cell_type", "symbol", "umis", "cells_detected", "cells_total"]).to_csv(
        RES / "target_celltype_summary.tsv", sep="\t", index=False)

    inc = dm[dm.included]
    qc = {"atlas": str(ATLAS.resolve()), "n_cells": int(n_cells), "n_libraries": int(n_lib), "n_genes": int(n_genes),
          "n_hepatocytes": int((ct_codes == hep).sum()), "min_hep_cells": MIN_HEP_CELLS,
          "n_donors_all": len(donors), "n_donors_included": int(len(inc)),
          "included_by_stage_dataset": {f"{k[0]}|{k[1]}": int(v) for k, v in
                                        inc.groupby(["disease_stage_coarse", "dataset"]).size().items()},
          "inclusion_counts": dm.inclusion.value_counts().to_dict(),
          "sex_agree_h34": int((inc.sex_call == inc.h34_sex).sum()), "sex_disagree_h34": int((inc.sex_call != inc.h34_sex).sum()),
          "sex_agree_repo_library": int((inc.sex_call == inc.repo_library_sex).sum()),
          "sex_disagree_repo_library": int((inc.sex_call != inc.repo_library_sex).sum()),
          "libraries_not_mapped_to_roster_donor": sorted(set(lib_meta.loc[~lib_meta.donor.isin(ros.donor), "library"]) )[:50]}
    json.dump(qc, open(RES / "qc.json", "w"), indent=2)
    log(json.dumps({k: qc[k] for k in ["n_donors_all", "n_donors_included", "inclusion_counts", "sex_disagree_h34",
                                        "sex_disagree_repo_library"]}))
    log(inc.groupby(["disease_stage_coarse", "dataset"]).size())


if __name__ == "__main__":
    main()
