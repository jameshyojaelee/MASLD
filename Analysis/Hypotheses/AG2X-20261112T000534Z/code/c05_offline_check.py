#!/usr/bin/env python3
"""Offline checks of the Stage C code path (no hosted call, no key loaded)."""
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import c01_stageC_lib as S  # noqa: E402
import c10_score_stageC as C10  # noqa: E402

out = {}
# 1. request inputs exactly as built for sending
for s in S.SETS:
    for st in S.STATES:
        a, nb = C10.load_cells(s, st)
        out[f"{s}_{st}"] = {"n_cells": int(a.n_obs), "n_genes": int(a.n_vars), "obs_columns": list(a.obs.columns),
                            "obs_index_head": a.obs_names[:2].tolist(), "request_bytes": int(nb)}
# 2. readout track rule on a synthetic tidy table
nm = "set1_otherhep"
rows = []
for p, mode, val in [("SHARE-seq-cellranger-nuclear", "imputed_supported", 1.0), ("scATAC-seq 10X-bowtie2-nuclear", "imputed_supported", 2.0),
                     ("sciATAC-seq-bowtie2-nuclear", "imputed_exploratory", 50.0), ("EasySci-ATAC-easy_sci-nuclear", "imputed_supported", 9.0),
                     ("scATAC-seq 10X-cellranger-nuclear", "imputed_supported", 9.0)]:
    rows.append({"variant_id": "chr1:10:A>G", "output_type": "SC_ATAC", "track_name": f"{p}_{nm}", "imputation_mode": mode, "raw_score": val})
t = pd.DataFrame(rows)
assert S.readout(t, nm).loc["chr1:10:A:G"] == 1.5
assert S.readout(t, nm, all_tracks=True).loc["chr1:10:A:G"] == np.mean([1, 2, 9, 9])
# 3. Stage B references
v = pd.read_csv(S.C / "inputs/variants_stageC.tsv", sep="\t")
rh, rs = S.reference(S.REF_HEP), S.reference(S.REF_STELLATE)
out["reference_hep_variants_covered"] = int(rh.reindex(v.uid).notna().sum())
out["reference_stellate_variants_covered"] = int(rs.reindex(v.uid).notna().sum())
# 4. LOO nearest centroid on separable and on label-shuffled data
r = np.random.default_rng(S.SEED)
E = np.vstack([r.normal(0, 1, (50, 8)) + 5 * np.eye(8)[0], r.normal(0, 1, (40, 8)) + 5 * np.eye(8)[1]])
y = np.array(["a"] * 50 + ["b"] * 40)
out["loo_separable"] = S.loo_nearest_centroid_bal_acc(E, y)
out["loo_shuffled"] = S.loo_nearest_centroid_bal_acc(E, r.permutation(y))
assert out["loo_separable"] > 0.95 and out["loo_shuffled"] < 0.75
out["planned_calls"] = 3 + int(v.pilot.sum()) + len(v) * 3 * 2
out["ledger_spent"] = S.spent()
json.dump(out, open(S.TAB / "offline_check.json", "w"), indent=1)
print(json.dumps(out, indent=1))
