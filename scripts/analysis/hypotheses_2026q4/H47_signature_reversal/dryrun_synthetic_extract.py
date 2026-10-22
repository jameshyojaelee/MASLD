"""Build a synthetic extract (real LINCS metadata, random z-scores) to dry-run h47_score.py.

Usage: python3 -I dryrun_synthetic_extract.py <downloads_untrusted_dir> <out_dir>
No L1000 expression value is read. Outputs are code checks only.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h47_extract import SERIES  # noqa: E402


def main(dl_dir, out_dir):
    os.makedirs(out_dir, exist_ok=False)
    rng = np.random.default_rng(7)
    gi = pd.read_csv(os.path.join(dl_dir, "GSE92742_Broad_LINCS_gene_info.txt.gz"), sep="\t", dtype=str)
    lm = gi.loc[gi.pr_is_lm == "1", "pr_gene_id"]
    lm.to_csv(os.path.join(out_dir, "landmark_order.tsv"), sep="\t", index=False)
    for series, f in SERIES.items():
        ci = pd.read_csv(os.path.join(dl_dir, f["cell_info"]), sep="\t", dtype=str)
        lines = ci.loc[ci.primary_site == "liver", "cell_id"]
        si = pd.read_csv(os.path.join(dl_dir, f["sig_info"]), sep="\t", dtype=str)
        sel = si[(si.pert_type == "trt_cp") & si.cell_id.isin(lines)].copy()
        sm = pd.read_csv(os.path.join(dl_dir, f["sig_metrics"]), sep="\t", dtype=str)
        sel = sel.merge(sm[["sig_id", "distil_cc_q75", "distil_nsample", "tas"]], on="sig_id", how="left")
        if "pert_idose" not in sel.columns:
            sel["pert_idose"] = ""
        sel["series"] = series
        sel.to_csv(os.path.join(out_dir, f"sig_meta_{series}.tsv"), sep="\t", index=False)
        np.save(os.path.join(out_dir, f"z_landmark_{series}.npy"),
                rng.normal(size=(len(sel), len(lm))).astype(np.float32))
    print("synthetic extract written")


if __name__ == "__main__":
    main(*sys.argv[1:])
