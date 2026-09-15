#!/usr/bin/env python3
"""C2 post-hoc diagnostic (NOT prespecified): how precise the caQTL label itself is.

The Currin lead file releases beta, the nominal p and the estimated residual degrees of freedom but
no standard error.  For a FastQTL nominal test, p comes from a t statistic on `estimated_DF` degrees
of freedom, so SE = |beta| / t with t = isf(p/2, df).  The measurement-error attenuation ceiling on
the Pearson correlation any predictor can reach against a noisy label is sqrt(1 - E[SE^2]/Var(beta)).

Caveat carried into the report: the lead set is selected at q < 0.05, so |beta| is winner's-cursed
and Var(beta) is inflated, which makes this ceiling optimistic rather than conservative.  It is a
rough bound on a Pearson scale, reported beside a Spearman statistic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    frame = pd.read_csv(arguments.labels, sep="\t")
    source = pd.read_csv(
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/data/"
        "seqfunc_external/currin2025_caqtl_v1/"
        "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
        sep="\t",
        usecols=["lead_variant_ID", "peak_ID", "estimated_DF"],
    ).rename(columns={"lead_variant_ID": "lead_variant_id", "peak_ID": "peak_id"})
    frame = frame.merge(source, on=["lead_variant_id", "peak_id"], how="left", validate="one_to_one")
    if frame["estimated_DF"].isna().any():
        raise SystemExit("estimated_DF join is incomplete")
    beta = frame["beta_alt"].to_numpy(np.float64)
    p = frame["p_nominal"].to_numpy(np.float64)
    df = frame["estimated_DF"].to_numpy(np.float64)
    t = stats.t.isf(p / 2.0, df)
    se = np.abs(beta) / t
    usable = np.isfinite(se) & (se > 0)
    mean_error_variance = float(np.mean(se[usable] ** 2))
    observed_variance = float(np.var(beta[usable], ddof=1))
    reliability = 1.0 - mean_error_variance / observed_variance
    report = {
        "schema_version": "agp-c2-label-precision-v1",
        "prespecified": False,
        "variants": int(usable.sum()),
        "median_abs_beta": float(np.median(np.abs(beta[usable]))),
        "median_se": float(np.median(se[usable])),
        "median_abs_t": float(np.median(t[usable])),
        "estimated_df_median": float(np.median(df[usable])),
        "mean_error_variance": mean_error_variance,
        "observed_beta_variance": observed_variance,
        "implied_label_reliability": reliability,
        "attenuation_ceiling_pearson": float(np.sqrt(max(reliability, 0.0))),
        "caveat": (
            "leads are selected at q<0.05 so Var(beta) is winner's-cursed and inflated; this ceiling "
            "is optimistic. Pearson scale, reported beside a Spearman statistic."
        ),
    }
    arguments.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
