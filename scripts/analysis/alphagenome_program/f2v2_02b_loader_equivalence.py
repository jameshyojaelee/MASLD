#!/usr/bin/env python3
"""f2-haplotype-v2 step 2b: prove this package's scorer is the same loader as v1's. 4 API calls.

v2 re-implements score_pair rather than importing it from f2_02_score.py, because it drops the determinism
probe and records the two readout widths. A null must regenerate its statistic through the SAME loader as the
observed value, so that re-implementation has to be checked rather than asserted.

The check: rescore one of v1's own central-20kb chr6 null draws with THIS package's scorer and compare all
sixteen summed values with the numbers v1 deposited. v1 measured the REF arm as bit-identical on a repeated
call at this window, so any difference here is the loader, not the model.
"""

from __future__ import annotations

import json
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import f2_recipe as R                              # noqa: E402
from f2v2_02_score_genespan_null import Scorer     # noqa: E402

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
V1 = PROJECT / "GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z"
PROBE_TAG = "gnmtnull0"


def main() -> None:
    ref = None
    with open(V1 / "raw" / "scored_rows.jsonl") as handle:
        for line in handle:
            r = json.loads(line)
            if r["tag"] == PROBE_TAG:
                ref = r
                break
    if ref is None:
        raise SystemExit(f"{PROBE_TAG} not in v1's checkpoint")

    sc = Scorer()
    got = sc.score_pair(ref["chrom"], int(ref["pos1"]), ref["ref1"], ref["alt1"],
                        int(ref["pos2"]), ref["ref2"], ref["alt2"], ref["ensembl"], f"probe:{PROBE_TAG}")
    keys = [f"{ch}_{arm}_sum" for ch in R.CHANNELS for arm in ("ref", "v1", "v2", "joint")]
    cmp = {k: {"v1": float(ref[k]), "v2": float(got[k]),
               "identical": float(ref[k]) == float(got[k]),
               "relative_difference": (abs(float(got[k]) - float(ref[k])) / abs(float(ref[k]))
                                      if float(ref[k]) else None)} for k in keys}
    res = {"probe_tag": PROBE_TAG, "readout_branch_v1": ref["rna_readout"], "readout_branch_v2": got["rna_readout"],
           "readout_branch_identical": ref["rna_readout"] == got["rna_readout"],
           "n_sums_compared": len(keys),
           "n_identical": int(sum(1 for v in cmp.values() if v["identical"])),
           "max_relative_difference": max((v["relative_difference"] or 0.0) for v in cmp.values()),
           "residuals": {ch: {"v1": float(ref[f"{ch}_residual"]), "v2": float(got[f"{ch}_residual"])}
                         for ch in R.CHANNELS},
           "per_sum": cmp,
           "verdict": ("v2's scorer reproduces v1's deposited sums exactly, so the repaired null is regenerated "
                       "through the same loader"
                       if all(v["identical"] for v in cmp.values()) else
                       "v2's scorer does NOT reproduce v1's sums; the repaired null is not comparable to v1's "
                       "and every verdict below must be read as a different loader")}
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    json.dump(res, (OUT / "tables" / "f2v2_loader_equivalence.json").open("w"), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("probe_tag", "readout_branch_identical", "n_identical",
                                          "n_sums_compared", "max_relative_difference", "verdict")}, indent=1))


if __name__ == "__main__":
    main()
