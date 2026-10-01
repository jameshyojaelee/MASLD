#!/usr/bin/env python3
"""Reconstruct high-impact historical estimates directly from deposited rows."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]


def macro(y, p, folds):
    values = [spearmanr(y[folds == f], p[folds == f]).statistic for f in sorted(set(folds))]
    if not np.isfinite(values).all():
        raise ValueError("undefined fold correlation")
    return float(np.tanh(np.mean(np.arctanh(values))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("run numerical scientific reconstruction on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    prog = ROOT / "GWAS/finemapping/results/alphagenome_program"
    path = prog / "ad1-gate-20260915T114615Z/raw/gate_1048576.tsv"
    raw = pd.read_csv(path, sep="\t")
    assert len(raw) == 3627 and not raw.key.duplicated().any()
    old = pd.read_csv(prog/"c2-endogenous-head-20260914T194000Z/tables/oof_predictions.tsv.gz", sep="\t")
    old["key"] = old.variant_id.str.removeprefix("chr")
    joined = raw.merge(old[["key", "hyenadna__delta_ridge__ensemble"]], on="key", validate="one_to_one")
    assert len(joined) == 3627
    y, f = joined.beta_alt.to_numpy(), joined.heldout_fold.to_numpy()
    native = macro(y, joined.local_atac_liver.to_numpy(), f)
    hyena = macro(y, joined["hyenadna__delta_ridge__ensemble"].to_numpy(), f)
    assert abs(native-0.7441117690812528) < 1e-12
    assert abs(hyena-0.23441130981965008) < 1e-12
    splice_path = ROOT/"GWAS/finemapping/results/alphagenome_atlas/p3a-splice-direction-20260915T141415Z/tables/splice_direction_pairs.tsv"
    sp = pd.read_csv(splice_path, sep="\t")
    good = np.isfinite(sp.excision_ratio_log2) & np.isfinite(sp.slope) & sp.excision_ratio_log2.ne(0)
    sp = sp.loc[good]
    p, m = sp.excision_ratio_log2.to_numpy(), sp.slope.to_numpy()
    concordance = float(np.mean((p > 0) == (m > 0)))
    marginal = float(np.mean(p > 0)*np.mean(m > 0) + np.mean(p <= 0)*np.mean(m <= 0))
    magnitude = float(spearmanr(np.abs(p), np.abs(m)).statistic)
    assert len(sp) == 637
    assert abs(concordance-0.5777080062794349) < 1e-12
    assert abs(marginal-0.5004325120943197) < 1e-12
    assert abs(magnitude-0.03569716222835225) < 1e-12
    report = dict(status="pass", independent_reconstruction="scipy_spearman_and_direct_sign_arithmetic_no_producer_metric_import",
                  local_native=dict(n=3627, macro_spearman=native, matched_hyenadna_macro=hyena,
                                    population="selected_subset_not_full32322", source=str(path)),
                  splice=dict(n=637, concordance=concordance, marginal=marginal,
                              magnitude_spearman=magnitude, source=str(splice_path),
                              original_thresholds="P1>0.65_and_p<0.01;D1<=0.60_or_p>=0.05_inconsistent;observed_below_both"),
                  seed="not_applicable_deterministic_point_estimates", protected_outcomes_read=False)
    (args.out/"rederived_predictions.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
