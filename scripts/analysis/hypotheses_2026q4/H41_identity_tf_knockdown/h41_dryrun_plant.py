"""H41 dry run only: plant the disease direction into synthetic HNF4A and HNF1A hairpin z-scores.

Usage: python3 -I h41_dryrun_plant.py <synthetic_extract_dir> <deg_results.csv>

Refuses to run unless extract_report.json says synthetic. Adds 1.0 x standardized
refit logFC (mapped landmarks) to every HNF4A and HNF1A hairpin signature, so the
call path and the residual test beyond HNF4A are exercised.
"""
import json
import os
import sys

import numpy as np
import pandas as pd


def main(ext_dir, deg_path):
    rep = json.load(open(os.path.join(ext_dir, "extract_report.json")))
    if not rep.get("synthetic"):
        sys.exit("refusing: extract is not synthetic")
    lm = pd.read_csv(os.path.join(ext_dir, "landmark_order.tsv"), sep="\t", dtype=str)
    mp = pd.read_csv(os.path.join(ext_dir, "landmark_map.tsv"), sep="\t", dtype=str).set_index("pr_gene_id")
    ens = mp.reindex(lm.pr_gene_id).ensembl
    deg = pd.read_csv(deg_path)
    deg.index = deg.gene.str.replace(r"\.[0-9]+(_PAR_Y)?$", "", regex=True)
    lfc = ens.map(deg.logFC).fillna(0.0).to_numpy()
    lfc = lfc / lfc.std()
    meta = pd.read_csv(os.path.join(ext_dir, "meta_trt_sh.tsv"), sep="\t", dtype=str)
    z = np.load(os.path.join(ext_dir, "z_trt_sh.npy"))
    rows = meta.pert_iname.isin(["HNF4A", "HNF1A"]).to_numpy()
    z[rows] += lfc.astype(np.float32)
    np.save(os.path.join(ext_dir, "z_trt_sh.npy"), z)
    rep["planted_rows"] = int(rows.sum())
    json.dump(rep, open(os.path.join(ext_dir, "extract_report.json"), "w"), indent=2)
    print("planted", int(rows.sum()), "rows")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
