#!/usr/bin/env python3
"""Stage D 2026-11-18 checks: verify < 1e-8 agreement between Stage D matched-window slice
(w = 1201 if out == "SC_ATAC" else 501, np.log2(sa + 1.0) - np.log2(sr + 1.0)) and Stage B
ChromatinScorer raw parquet scores across ATAC, DNASE, and SC_ATAC.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402
import d20_panels_stageD_20261118 as P  # noqa: E402


def load_stageB_raw_map(uids: set[str]) -> dict[tuple[str, str, str], float]:
    cols = ["variant_id", "output_type", "imputation_mode", "raw_score"]
    rows = []
    for d in sorted((L.PKG / "stageB/raw").glob("t?_*")):
        cond = d.name.split("_", 1)[1]
        if cond not in ("bulk_liver", "hepatocyte", "monocyte_ctrl"):
            continue
        for f in sorted(d.glob("chunk_*.parquet")):
            x = pd.read_parquet(f, columns=cols)
            x = x[(x.imputation_mode == "observed") & x.output_type.isin(["ATAC", "DNASE", "SC_ATAC"])]
            if x.empty:
                continue
            x["uid"] = x.variant_id.str.replace(">", ":", regex=False)
            x = x[x.uid.isin(uids)]
            if not x.empty:
                x["condition"] = cond
                rows.append(x[["uid", "condition", "output_type", "raw_score"]])
    df = pd.concat(rows, ignore_index=True)
    g = df.groupby(["uid", "condition", "output_type"]).raw_score.mean()
    return {k: float(v) for k, v in g.items()}


def main():
    loc = pd.read_csv(P.TAB / "stageD_local_change.tsv", sep="\t")
    stageB_exact = load_stageB_raw_map(set(loc.uid.unique()))

    rows = []
    for _, x in loc.iterrows():
        z, meta = P.load_tracks(x.uid, x.condition)
        m = meta["outputs"][x.output]
        res, vb = m["resolution"], m["variant_bin"]
        w = P.window_bp(x.output)
        h = (w // 2) // res
        ref_slice = z[f"{x.output}_ref"][vb - h : vb + h + 1]
        alt_slice = z[f"{x.output}_alt"][vb - h : vb + h + 1]

        sr64 = ref_slice.astype(np.float64).sum(axis=0)
        sa64 = alt_slice.astype(np.float64).sum(axis=0)
        local_f64 = float(np.mean(np.log2(sa64 + 1.0) - np.log2(sr64 + 1.0)))

        sr32 = ref_slice.astype(np.float32).sum(axis=0)
        sa32 = alt_slice.astype(np.float32).sum(axis=0)
        local_f32 = float(np.mean(np.log2(sa32.astype(np.float64) + 1.0) - np.log2(sr32.astype(np.float64) + 1.0)))

        sb_exact = stageB_exact.get((x.uid, x.condition, x.output), float(x.stageB_raw))
        rows.append(
            {
                "uid": x.uid,
                "condition": x.condition,
                "output": x.output,
                "window_bp": w,
                "stageB_raw_tsv": x.stageB_raw,
                "stageB_raw_parquet": sb_exact,
                "mean_sum_ref": float(sr64.mean()),
                "mean_sum_alt": float(sa64.mean()),
                "local_pc1_f64": local_f64,
                "local_pc1_f32": local_f32,
                "absdiff_f64_vs_parquet": abs(local_f64 - sb_exact),
                "absdiff_f32_vs_parquet": abs(local_f32 - sb_exact),
                "absdiff_f64_vs_tsv": abs(local_f64 - x.stageB_raw),
                "candidate_not_adopted": True,
            }
        )

    a = pd.DataFrame(rows)
    a.to_csv(P.TAB / "stageD_local_vs_scorer_check.tsv", sep="\t", index=False)

    per_output_max = a.groupby("output")[["absdiff_f64_vs_parquet", "absdiff_f32_vs_parquet", "absdiff_f64_vs_tsv"]].max().to_dict("index")
    summary = {
        "max_absdiff_f64_vs_parquet": float(a.absdiff_f64_vs_parquet.max()),
        "max_absdiff_f32_vs_parquet": float(a.absdiff_f32_vs_parquet.max()),
        "max_absdiff_f64_vs_tsv": float(a.absdiff_f64_vs_tsv.max()),
        "per_output_max": per_output_max,
    }
    (P.TAB / "stageD_local_vs_scorer_summary.json").write_text(json.dumps(summary, indent=1))
    L.log(f"Stage D 20261118 check summary: {json.dumps(summary)}")

    # SORT1 CEBPA span check (reading from stageD/motifs & stageD/raw/ism_*)
    u = "chr1:109274968:G:T"
    hits = pd.read_csv(P.D / "motifs/motif_hits_64bp.tsv.gz", sep="\t")
    h = hits[(hits.uid == u) & (hits.allele == "alt") & (hits.motif_id == "MA0102.5") & hits.covers_variant]
    b = h.sort_values("score", ascending=False).iloc[0]
    out = []
    for step in ("ism_ref", "ism_alt"):
        g = P.load_ism(u, step, "hepatocyte", "SC_ATAC")
        M, off = P.ism_matrix(g, u, step)
        imp = np.abs(np.nanmean(M, axis=0))
        lo, hi = off - 15, off + 16
        idx = np.arange(lo, hi + 1)
        inside = (idx - off >= int(b.start_rel)) & (idx - off <= int(b.end_rel))
        w = imp[lo : hi + 1]
        out.append(
            {
                "uid": u,
                "ism": step,
                "cebpa_hit_rel": f"{int(b.start_rel):+d}..{int(b.end_rel):+d} ({b.strand})",
                "mean_abs_importance_inside": float(w[inside].mean()),
                "mean_abs_importance_outside_central32": float(w[~inside].mean()),
                "ratio_inside_vs_outside_central32": float(w[inside].mean() / w[~inside].mean()),
                "candidate_not_adopted": True,
            }
        )
    o = pd.DataFrame(out)
    o.to_csv(P.TAB / "stageD_sort1_cebpa_span_ref_vs_alt.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
