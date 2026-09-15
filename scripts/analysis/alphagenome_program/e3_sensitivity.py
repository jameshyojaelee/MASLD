#!/usr/bin/env python3
"""E3 post-hoc sensitivity checks. Nothing here changes a prespecified verdict.

1. Degeneracy: are any predicted effects exactly zero or missing (a Spearman on a constant is
   undefined and would silently drop out of the bootstrap).
2. One block, chr6_lr0003, holds 26 of the 200 random and 14 of the 100 top-effect elements. Point
   estimates are recomputed without it, because a block bootstrap on 93 blocks is sensitive to a
   block that large.
3. Range expansion: the top-effect stratum is selected on |measured allele effect| in HepG2 control,
   so its Spearman is not an unbiased estimate of performance on a randomly drawn element. The
   label's spread in each stratum is reported so the size of that selection is visible.

Output (tables/): e3_sensitivity.tsv
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

BIG_BLOCK = "chr6_lr0003"
LABEL = "measured_HepG2_control"


def sp(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return float("nan")
    return float(stats.spearmanr(x[ok], y[ok]).statistic)


def main() -> None:
    out_root = pathlib.Path(sys.argv[1])
    tables = out_root / "tables"
    d = pd.read_csv(tables / "e3_predicted_effects.tsv", sep="\t")
    d = d[d["readout_window"] == "primary"]

    rows = []
    for arm in ("native", "npad"):
        for ch in ("atac", "dnase", "rna"):
            v = d[f"{arm}_{ch}"].to_numpy(float)
            rows.append({"check": "degeneracy", "stratum": "all", "arm": arm, "channel": ch,
                         "n": int(v.size), "value": float(np.nanmedian(np.abs(v))),
                         "detail": f"exact_zero={int((v == 0).sum())} nan={int(np.isnan(v).sum())}"})

    for stratum, col in (("random", "in_random_stratum"), ("topeffect", "in_topeffect_stratum")):
        s = d[d[col]]
        s2 = s[s["long_range_block"] != BIG_BLOCK]
        for ch in ("atac", "dnase", "rna"):
            for arm in ("native", "npad"):
                a = sp(s[f"{arm}_{ch}"].to_numpy(float), s[LABEL].to_numpy(float))
                b = sp(s2[f"{arm}_{ch}"].to_numpy(float), s2[LABEL].to_numpy(float))
                rows.append({"check": f"drop_{BIG_BLOCK}", "stratum": stratum, "arm": arm,
                             "channel": ch, "n": int(s2.shape[0]), "value": b,
                             "detail": f"with_all_blocks_n={s.shape[0]}_rho={a:+.4f}"})
        v = s[LABEL]
        rows.append({"check": "label_spread", "stratum": stratum, "arm": "", "channel": "",
                     "n": int(s.shape[0]), "value": float(v.std()),
                     "detail": (f"IQR={float(v.quantile(.75) - v.quantile(.25)):.3f} "
                                f"median_abs={float(v.abs().median()):.3f}")})

    pd.DataFrame(rows).to_csv(tables / "e3_sensitivity.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
