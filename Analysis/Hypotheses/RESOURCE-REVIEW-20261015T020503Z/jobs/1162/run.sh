#!/bin/bash
# Job 1162 (review-R16, K08): did any dataset in the single-cell atlas receive
# all-zero Scrublet scores (the exception fallback or the <30-cell skip)?
# Read-only. Reads only obs columns with h5py; never loads X.
set -euo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$PWD
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python - "$BASE" "$OUT" <<'PY'
import sys, os, glob, h5py, numpy as np
base, out = sys.argv[1], sys.argv[2]

def read_col(obs, name):
    o = obs[name]
    if isinstance(o, h5py.Group):  # categorical (anndata >=0.8)
        cats = o["categories"][()]
        codes = o["codes"][()]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats], dtype=object)
        res = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
        return res
    v = o[()]
    if v.dtype.kind in "SO":
        v = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v], dtype=object)
    return v

rows = []
files = [os.path.join(base, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")]
files += sorted(glob.glob(os.path.join(base, "Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset/*_v2.h5ad")))
files += [os.path.join(base, "Analysis/SingleCell/results_gpu_v2_phase05/atlas/scalesc_human_annotated_celltypist_v2.h5ad")]
log = open(os.path.join(out, "k08_obs_columns.txt"), "w")
for f in files:
    real = os.path.realpath(f)
    with h5py.File(f, "r") as h:
        obs = h["obs"]
        keys = list(obs.keys())
        log.write(f"{f}\t{real}\t{len(keys)} obs keys\t" + ",".join(keys) + "\n")
        score_cols = [k for k in keys if "scrublet" in k.lower() or "doublet" in k.lower()]
        ds_col = next((k for k in ["dataset", "dataset_id", "source_dataset", "batch", "study"] if k in keys), None)
        n = None
        ds = read_col(obs, ds_col) if ds_col else None
        if ds is not None:
            n = len(ds)
        for sc in score_cols:
            v = read_col(obs, sc)
            if v.dtype.kind not in "fiub":
                try:
                    v = v.astype(float)
                except Exception:
                    rows.append((os.path.basename(f), sc, "ALL", len(v), "non_numeric", "", "", ""))
                    continue
            v = v.astype(float)
            groups = np.unique(ds) if ds is not None else np.array(["ALL"])
            for g in groups:
                sel = (ds == g) if ds is not None else np.ones(len(v), bool)
                x = v[sel]
                rows.append((os.path.basename(f), sc, str(g), int(sel.sum()),
                             f"{np.mean(x == 0):.4f}", f"{np.nanmax(x):.4f}" if len(x) else "NA",
                             f"{np.nanmedian(x):.4f}" if len(x) else "NA",
                             int(np.sum(x > 0.25))))
        if not score_cols:
            rows.append((os.path.basename(f), "NO_SCRUBLET_OR_DOUBLET_COLUMN", ds_col or "no_dataset_col", n if n is not None else -1, "", "", "", ""))
log.close()
with open(os.path.join(out, "k08_scrublet_by_dataset.tsv"), "w") as fh:
    fh.write("file\tcolumn\tdataset\tn_cells\tfraction_exact_zero\tmax\tmedian\tn_gt_0.25\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(open(os.path.join(out, "k08_scrublet_by_dataset.tsv")).read())
PY
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
echo DONE
