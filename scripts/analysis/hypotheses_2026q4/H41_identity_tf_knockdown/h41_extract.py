"""H41 step A: extract HEPG2 shRNA landmark signatures from LINCS GSE92742 level 5.

Usage: python3 -I h41_extract.py <downloads_untrusted_dir> <decompressed_dir> <out_dir>

Selects GSE92742 sig_info rows with cell_id == HEPG2 and pert_type in
{trt_sh, trt_sh.cgs}; reads the 978 landmark rows of those columns from the
level 5 GCTX. No scoring happens here.
Stop rule S3 (PRESPEC_H41.md): fewer than 978 landmark rows, or any selected
sig_id missing from the GCTX columns.
"""
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

GCTX = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
SIG_INFO = "GSE92742_Broad_LINCS_sig_info.txt.gz"
GENE_INFO = "GSE92742_Broad_LINCS_gene_info.txt.gz"
CELL = "HEPG2"
TYPES = ["trt_sh", "trt_sh.cgs"]
BLOCK = 8000


def stop(out_dir, msg):
    with open(os.path.join(out_dir, "STOPPED.txt"), "w") as fh:
        fh.write(msg + "\n")
    sys.exit("STOP: " + msg)


def decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def main(dl_dir, gctx_dir, out_dir):
    os.makedirs(out_dir, exist_ok=False)
    gi = pd.read_csv(os.path.join(dl_dir, GENE_INFO), sep="\t", dtype=str)
    lm = gi.loc[gi.pr_is_lm == "1", ["pr_gene_id", "pr_gene_symbol"]].reset_index(drop=True)
    if len(lm) != 978:
        stop(out_dir, f"gene_info has {len(lm)} landmarks, expected 978")
    si = pd.read_csv(os.path.join(dl_dir, SIG_INFO), sep="\t", dtype=str)
    sel = si[(si.cell_id == CELL) & si.pert_type.isin(TYPES)].copy()

    with h5py.File(os.path.join(gctx_dir, GCTX), "r") as h5:
        row_ids = decode(h5["/0/META/ROW/id"][:])
        col_ids = decode(h5["/0/META/COL/id"][:])
        dset = h5["/0/DATA/0/matrix"]
        shape = dset.shape
        if shape != (len(col_ids), len(row_ids)):
            stop(out_dir, f"matrix shape {shape} does not match (cols, rows)")
        row_pos = {r: i for i, r in enumerate(row_ids)}
        missing_lm = [g for g in lm.pr_gene_id if g not in row_pos]
        if missing_lm:
            stop(out_dir, f"{len(missing_lm)} landmark rows missing")
        lm_idx = np.array([row_pos[g] for g in lm.pr_gene_id])
        col_pos = {c: i for i, c in enumerate(col_ids)}
        n_missing = int((~sel.sig_id.isin(col_pos)).sum())
        if n_missing:
            stop(out_dir, f"{n_missing} of {len(sel)} selected sig_ids missing from GCTX")
        sel["col_idx"] = sel.sig_id.map(col_pos).astype(np.int64)
        sel = sel.sort_values("col_idx").reset_index(drop=True)
        want = sel.col_idx.to_numpy()
        out = np.empty((len(want), len(lm_idx)), dtype=np.float32)
        filled = 0
        for start in range(0, shape[0], BLOCK):
            stop_i = min(start + BLOCK, shape[0])
            lo = np.searchsorted(want, start)
            hi = np.searchsorted(want, stop_i)
            if hi == lo:
                continue
            block = dset[start:stop_i, :]
            out[lo:hi, :] = block[want[lo:hi] - start][:, lm_idx]
            filled += hi - lo
        assert filled == len(want)

    report = {"matrix_shape": list(shape), "nan_values": int(np.isnan(out).sum())}
    for t in TYPES:
        m = (sel.pert_type == t).to_numpy()
        np.save(os.path.join(out_dir, f"z_{t}.npy"), out[m])
        sel.loc[m].drop(columns="col_idx").to_csv(os.path.join(out_dir, f"meta_{t}.tsv"), sep="\t", index=False)
        report[t] = {"signatures": int(m.sum()),
                     "hairpins": int(sel.loc[m, "pert_id"].nunique()),
                     "genes": int(sel.loc[m, "pert_iname"].nunique())}
    lm.to_csv(os.path.join(out_dir, "landmark_order.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "extract_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
