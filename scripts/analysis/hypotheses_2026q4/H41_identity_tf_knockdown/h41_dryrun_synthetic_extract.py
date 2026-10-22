"""H41 dry run: build a synthetic extract with the real HEPG2 shRNA metadata and random z-scores.

Usage: python3 -I h41_dryrun_synthetic_extract.py <downloads_untrusted_dir> <out_dir>

No L1000 expression value is read. z ~ N(0, 1), numpy default_rng(1) (dry run only).
"""
import json
import os
import sys

import numpy as np
import pandas as pd


def main(dl_dir, out_dir):
    os.makedirs(out_dir, exist_ok=False)
    gi = pd.read_csv(os.path.join(dl_dir, "GSE92742_Broad_LINCS_gene_info.txt.gz"), sep="\t", dtype=str)
    lm = gi.loc[gi.pr_is_lm == "1", ["pr_gene_id", "pr_gene_symbol"]].reset_index(drop=True)
    si = pd.read_csv(os.path.join(dl_dir, "GSE92742_Broad_LINCS_sig_info.txt.gz"), sep="\t", dtype=str)
    rng = np.random.default_rng(1)
    rep = {"synthetic": True}
    for t in ("trt_sh", "trt_sh.cgs"):
        m = si[(si.cell_id == "HEPG2") & (si.pert_type == t)].reset_index(drop=True)
        np.save(os.path.join(out_dir, f"z_{t}.npy"), rng.normal(size=(len(m), 978)).astype(np.float32))
        m.to_csv(os.path.join(out_dir, f"meta_{t}.tsv"), sep="\t", index=False)
        rep[t] = len(m)
    lm.to_csv(os.path.join(out_dir, "landmark_order.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "extract_report.json"), "w") as fh:
        json.dump(rep, fh, indent=2)
    print(rep)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
