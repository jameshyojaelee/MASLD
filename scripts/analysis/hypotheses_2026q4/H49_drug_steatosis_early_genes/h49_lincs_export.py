"""H49: export LINCS L1000 level 5 PHH signatures (landmark genes) from the H47 extract.

Usage: python3 -I h49_lincs_export.py <h47_extract_dir> <out_tsv_gz>

Reads (read-only) sig_meta_GSE92742.tsv, z_landmark_GSE92742.npy and
landmark_order.tsv from the H47 extract. Rows of the z matrix follow the rows
of sig_meta; columns follow landmark_order (Entrez IDs). Writes PHH rows only:
signature metadata columns followed by one column per landmark Entrez ID.
"""
import os
import sys

import numpy as np
import pandas as pd


def main(extract_dir: str, out_path: str) -> None:
    if os.path.exists(out_path):
        raise SystemExit(f"refusing to overwrite {out_path}")
    meta = pd.read_csv(os.path.join(extract_dir, "sig_meta_GSE92742.tsv"), sep="\t", dtype=str)
    z = np.load(os.path.join(extract_dir, "z_landmark_GSE92742.npy"), mmap_mode="r")
    genes = pd.read_csv(os.path.join(extract_dir, "landmark_order.tsv"), sep="\t", dtype=str)["pr_gene_id"]
    if z.shape != (len(meta), len(genes)):
        raise SystemExit(f"shape mismatch: z {z.shape}, meta {len(meta)}, genes {len(genes)}")
    keep = np.flatnonzero((meta["cell_id"] == "PHH").to_numpy())
    cols = ["sig_id", "pert_id", "pert_iname", "pert_type", "cell_id", "pert_dose", "pert_dose_unit",
            "pert_time", "pert_time_unit", "distil_cc_q75", "distil_nsample"]
    out = meta.iloc[keep][cols].reset_index(drop=True)
    zz = pd.DataFrame(np.asarray(z[keep, :], dtype=np.float64), columns=["g" + g for g in genes])
    out = pd.concat([out, zz], axis=1)
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"PHH signatures: {len(out)}; landmarks: {len(genes)}; nan: {int(np.isnan(zz.to_numpy()).sum())}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
