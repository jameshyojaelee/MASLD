#!/usr/bin/env python3
"""E2 step 4: the measured-label ceiling on transfer.

Added after the prespecification was sealed, as a diagnostic, not as a new prediction.  It is
model-free: the agreement between the two MEASURED labels on the same variants - the GSE281364
reporter allele effect and the Currin caQTL beta - bounds how much any model trained on one can
predict about the other.  A model cannot transfer a relation the measurements do not share.

Both labels are already oriented to the same allele (the MPRA `genomic_alt`), so a positive Spearman
means the allele that raises reporter activity also raises accessibility.  Uncertainty is the same
10,000-draw 1-Mb block bootstrap at seed 20260914 used by the transfer contrast, resampled within
the same five chromosome-group folds.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from c2_03_fit_heads import fast_spearman, interval  # noqa: E402
from e2_02_fit_transfer import STRATA, block_bootstrap, bootstrap_p, macro_over_folds  # noqa: E402

REPORTER_COLUMNS = {
    "hepg2_mean": "d_hepg2_mean",
    "hepg2_control": "d_HepG2_control",
    "hepg2_paoa": "d_HepG2_PAOA",
    "all_context_mean": "d_all_context_mean",
}
ENDOGENOUS_COLUMNS = {"nearest_peak": "beta_nearest", "minp_peak": "beta_minp"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resamples", type=int, default=10000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260914)
    arguments = parser.parse_args()
    (arguments.output / "tables").mkdir(parents=True, exist_ok=True)

    units = pd.read_csv(arguments.inputs / "e2_units.tsv.gz", sep="\t").sort_values("row_index")
    folds = units.chrom_fold.to_numpy(np.int64)
    blocks = units.block_1mb.to_numpy(dtype=str)
    absolute = units.abs_distance_nearest.to_numpy(np.float64)

    rows_out = []
    for endogenous, endogenous_column in ENDOGENOUS_COLUMNS.items():
        beta = units[endogenous_column].to_numpy(np.float64)
        for stratum, threshold in STRATA:
            keep = np.flatnonzero(np.isfinite(beta) & (absolute <= threshold))
            if keep.size == 0:
                continue
            arms = {
                name: units[column].to_numpy(np.float64)[keep]
                for name, column in REPORTER_COLUMNS.items()
            }
            samples = block_bootstrap(
                beta[keep], arms, folds[keep], blocks[keep],
                resamples=arguments.resamples, seed=arguments.bootstrap_seed,
            )
            for name, values in arms.items():
                point = macro_over_folds(beta[keep], values, folds[keep])
                low, high, standard_error = interval(samples[name])
                rows_out.append({
                    "endogenous_label": endogenous,
                    "reporter_label": name,
                    "stratum_max_abs_distance_bp": stratum,
                    "variants": int(keep.size),
                    "blocks_1mb": int(len(set(blocks[keep].tolist()))),
                    "macro_spearman_measured_vs_measured": point,
                    "ci_low": low,
                    "ci_high": high,
                    "bootstrap_se": standard_error,
                    "bootstrap_p_two_sided": bootstrap_p(samples[name]),
                    "pooled_spearman": fast_spearman(beta[keep], values),
                    "pearson": float(np.corrcoef(beta[keep], values)[0, 1]),
                })
            print(f"{endogenous}/{stratum} done", flush=True)
    frame = pd.DataFrame(rows_out)
    frame.to_csv(arguments.output / "tables" / "measured_label_ceiling.tsv", sep="\t", index=False)
    print(frame.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
