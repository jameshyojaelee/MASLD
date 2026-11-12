#!/usr/bin/env python3
"""Stage D checks. No hosted call.
(a) Local change from saved tracks against the Stage B scorer value: per observed track,
    log2(sum ALT + c) - log2(sum REF + c) over 501 bp, c in {0, 1}, mean over tracks.
(b) SORT1 expectation 3, second half: mean |importance| inside the gained CEBPA hit in REF ISM and in
    ALT-background ISM, and the inside/outside ratio in REF ISM restricted to the same central 32 bp.
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402
import d20_panels_stageD as P  # noqa: E402


def main():
    loc = pd.read_csv(P.TAB / "stageD_local_change.tsv", sep="\t")
    rows = []
    for _, x in loc.iterrows():
        z, meta = P.load_tracks(x.uid, x.condition)
        m = meta["outputs"][x.output]
        res, vb = m["resolution"], m["variant_bin"]
        h = (P.LOCAL_W // 2) // res
        sr = z[f"{x.output}_ref"][vb - h:vb + h + 1].sum(axis=0)
        sa = z[f"{x.output}_alt"][vb - h:vb + h + 1].sum(axis=0)
        r = {"uid": x.uid, "condition": x.condition, "output": x.output, "stageB_raw": x.stageB_raw,
             "mean_sum_ref_501bp": float(sr.mean()), "mean_sum_alt_501bp": float(sa.mean())}
        for c in (0.0, 1.0):
            r[f"local_pc{int(c)}"] = float(np.mean(np.log2(sa + c) - np.log2(sr + c)))
            r[f"absdiff_vs_stageB_pc{int(c)}"] = abs(r[f"local_pc{int(c)}"] - x.stageB_raw)
        r["candidate_not_adopted"] = True
        rows.append(r)
    a = pd.DataFrame(rows)
    a.to_csv(P.TAB / "stageD_local_vs_scorer_check.tsv", sep="\t", index=False)
    L.log(f"local check: max absdiff pc0 {a.absdiff_vs_stageB_pc0.max():.4f}; pc1 {a.absdiff_vs_stageB_pc1.max():.4f}")

    u = "chr1:109274968:G:T"
    hits = pd.read_csv(P.D / "motifs/motif_hits_64bp.tsv.gz", sep="\t")
    h = hits[(hits.uid == u) & (hits.allele == "alt") & (hits.motif_id == "MA0102.5") & hits.covers_variant]
    b = h.sort_values("score", ascending=False).iloc[0]
    out = []
    for step in ("ism_ref", "ism_alt"):
        g = P.load_ism(u, step, "hepatocyte", "SC_ATAC")
        M, off = P.ism_matrix(g, u, step)
        imp = np.abs(np.nanmean(M, axis=0))
        lo, hi = off - 15, off + 16  # central 32 bp in both windows ([-15, +16] relative to the variant)
        idx = np.arange(lo, hi + 1)
        inside = (idx - off >= int(b.start_rel)) & (idx - off <= int(b.end_rel))
        w = imp[lo:hi + 1]
        out.append({"uid": u, "ism": step, "cebpa_hit_rel": f"{int(b.start_rel):+d}..{int(b.end_rel):+d} ({b.strand})",
                    "mean_abs_importance_inside": float(w[inside].mean()),
                    "mean_abs_importance_outside_central32": float(w[~inside].mean()),
                    "ratio_inside_vs_outside_central32": float(w[inside].mean() / w[~inside].mean()),
                    "candidate_not_adopted": True})
    o = pd.DataFrame(out)
    o.to_csv(P.TAB / "stageD_sort1_cebpa_span_ref_vs_alt.tsv", sep="\t", index=False)
    L.log(o.to_string())


if __name__ == "__main__":
    main()
