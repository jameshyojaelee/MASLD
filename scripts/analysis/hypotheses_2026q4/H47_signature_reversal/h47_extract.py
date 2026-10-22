"""H47 step A: extract landmark x liver-line compound signatures from LINCS level 5 GCTX.

Usage: python3 -I h47_extract.py <downloads_untrusted_dir> <decompressed_dir> <out_dir>

Reads only metadata and the level 5 matrices. No scoring happens here.
Stop rule S2 (PRESPEC_H47.md): fewer than 978 landmark rows, or more than 1%
of selected sig_ids missing from the GCTX columns, stops the run.
"""
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

LIVER_LINES = None  # decided from cell_info (primary_site == liver)
BLOCK = 8000

SERIES = {
    "GSE92742": {
        "gctx": "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
        "sig_info": "GSE92742_Broad_LINCS_sig_info.txt.gz",
        "sig_metrics": "GSE92742_Broad_LINCS_sig_metrics.txt.gz",
        "cell_info": "GSE92742_Broad_LINCS_cell_info.txt.gz",
    },
    "GSE70138": {
        "gctx": "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
        "sig_info": "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
        "sig_metrics": "GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz",
        "cell_info": "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz",
    },
}


def stop(out_dir, msg):
    with open(os.path.join(out_dir, "STOPPED.txt"), "w") as fh:
        fh.write(msg + "\n")
    sys.exit("STOP: " + msg)


def decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def main(dl_dir, gctx_dir, out_dir):
    os.makedirs(out_dir, exist_ok=False)
    gi = pd.read_csv(os.path.join(dl_dir, "GSE92742_Broad_LINCS_gene_info.txt.gz"), sep="\t", dtype=str)
    lm_ids = gi.loc[gi.pr_is_lm == "1", "pr_gene_id"].tolist()
    if len(lm_ids) != 978:
        stop(out_dir, f"gene_info has {len(lm_ids)} landmarks, expected 978")
    report = {}
    metas = []
    for series, f in SERIES.items():
        ci = pd.read_csv(os.path.join(dl_dir, f["cell_info"]), sep="\t", dtype=str)
        lines = sorted(ci.loc[ci.primary_site == "liver", "cell_id"])
        si = pd.read_csv(os.path.join(dl_dir, f["sig_info"]), sep="\t", dtype=str)
        sel = si[(si.pert_type == "trt_cp") & si.cell_id.isin(lines)].copy()
        sm = pd.read_csv(os.path.join(dl_dir, f["sig_metrics"]), sep="\t", dtype=str)
        sm = sm[["sig_id", "distil_cc_q75", "distil_nsample", "tas"]]
        sel = sel.merge(sm, on="sig_id", how="left")
        if "pert_idose" not in sel.columns:
            sel["pert_idose"] = ""
        sel["series"] = series

        with h5py.File(os.path.join(gctx_dir, f["gctx"]), "r") as h5:
            row_ids = decode(h5["/0/META/ROW/id"][:])
            col_ids = decode(h5["/0/META/COL/id"][:])
            dset = h5["/0/DATA/0/matrix"]
            shape = dset.shape
            # GCTX stores the matrix as (n_columns, n_rows) = (signatures, genes)
            if shape != (len(col_ids), len(row_ids)):
                stop(out_dir, f"{series}: matrix shape {shape} does not match (cols, rows) = ({len(col_ids)}, {len(row_ids)})")
            row_pos = {r: i for i, r in enumerate(row_ids)}
            missing_lm = [g for g in lm_ids if g not in row_pos]
            if missing_lm:
                stop(out_dir, f"{series}: {len(missing_lm)} landmark rows missing")
            lm_idx = np.array([row_pos[g] for g in lm_ids])
            col_pos = {c: i for i, c in enumerate(col_ids)}
            present = sel.sig_id.isin(col_pos)
            n_missing = int((~present).sum())
            if n_missing > 0.01 * len(sel):
                stop(out_dir, f"{series}: {n_missing} of {len(sel)} selected sig_ids missing from GCTX")
            sel = sel[present].copy()
            sel["col_idx"] = sel.sig_id.map(col_pos).astype(np.int64)
            sel = sel.sort_values("col_idx").reset_index(drop=True)
            want = sel.col_idx.to_numpy()
            out = np.empty((len(want), len(lm_idx)), dtype=np.float32)
            # sequential block reads; pick selected signatures from each block
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
        np.save(os.path.join(out_dir, f"z_landmark_{series}.npy"), out)
        sel.drop(columns="col_idx").to_csv(os.path.join(out_dir, f"sig_meta_{series}.tsv"), sep="\t", index=False)
        report[series] = {
            "liver_lines": lines,
            "selected_signatures": int(len(sel) + n_missing),
            "missing_from_gctx": n_missing,
            "extracted": int(len(sel)),
            "per_line": sel.cell_id.value_counts().to_dict(),
            "matrix_shape": list(shape),
            "nan_values": int(np.isnan(out).sum()),
        }
        metas.append(sel)
    pd.Series(lm_ids, name="pr_gene_id").to_csv(os.path.join(out_dir, "landmark_order.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "extract_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
