#!/usr/bin/env python3
"""Independently reconstruct adaptation errors with resampled group totals."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node execution required")
    args.out.mkdir(parents=True, exist_ok=False)
    table = pd.read_csv(args.review / "ten_paired_contrasts.tsv", sep="\t")
    assert len(table) == 20
    arms = sorted(set(table.arm) | set(table.comparator))
    errors, correlations, hashes = {}, {}, {}
    observed = None
    ids = None
    for arm in arms:
        path = args.sweep / arm / "validation_predictions.tsv"
        records = pd.read_csv(path, sep="\t").sort_values("variant_id")
        assert len(records) == 6834 and records.variant_id.is_unique and records.validation_fold.eq(1).all()
        if ids is None:
            ids = records.variant_id.to_numpy()
            observed = records.observed_beta.to_numpy()
        np.testing.assert_array_equal(records.variant_id, ids)
        np.testing.assert_array_equal(records.observed_beta, observed)
        errors[arm] = (records.predicted_beta.to_numpy() - observed)**2
        correlations[arm] = spearmanr(observed, records.predicted_beta).statistic
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    error_matrix = np.column_stack([errors[arm] for arm in arms])
    point_mse = error_matrix.mean(axis=0)
    identities = [value.split(":") for value in ids]
    chromosomes = np.array([x[0] for x in identities])
    historical_bins = np.array([f"{x[0]}:{int(x[1])//1000000}" for x in identities])
    outputs = []
    for name, labels, seed in (("historical_1mb_bins", historical_bins, 20260915),
                               ("whole_chromosomes", chromosomes, 20260916)):
        family = table.loc[table.uncertainty_unit.eq(name)].reset_index(drop=True)
        assert len(family) == 10 and family.arm.is_unique
        units = np.unique(labels)
        assert len(units) == (524 if name == "historical_1mb_bins" else 5)
        # Unlike the producer, sum squared errors within each unit once.
        # Resampling these totals reproduces the paired variant-weighted errors.
        totals = np.array([error_matrix[labels == unit].sum(axis=0) for unit in units])
        counts = np.array([(labels == unit).sum() for unit in units])
        selected = np.random.default_rng(seed).integers(len(units), size=(2000,len(units)))
        sampled_mse = totals[selected].sum(axis=1) / counts[selected].sum(axis=1)[:, None]
        p_values = []
        for row in family.itertuples(index=False):
            a, c = arms.index(row.arm), arms.index(row.comparator)
            pool = "symmetric" if row.arm.endswith("symmetric") else "variant"
            assert row.comparator == f"frozen_2048_{pool}_shared_mlp64"
            np.testing.assert_allclose([row.native_RMSE, row.comparator_native_RMSE], np.sqrt(point_mse[[a,c]]), rtol=0, atol=1e-13)
            np.testing.assert_allclose([row.signed_spearman, row.comparator_signed_spearman],
                [correlations[row.arm], correlations[row.comparator]], rtol=0, atol=1e-13)
            gains = sampled_mse[:,c] - sampled_mse[:,a]
            rmse_gains = np.sqrt(sampled_mse[:,c]) - np.sqrt(sampled_mse[:,a])
            np.testing.assert_allclose([row.MSE_improvement_beta_squared_low95,row.MSE_improvement_beta_squared_high95],
                np.quantile(gains,[.025,.975]), rtol=0, atol=1e-12)
            np.testing.assert_allclose([row.RMSE_improvement_beta_low95,row.RMSE_improvement_beta_high95],
                np.quantile(rmse_gains,[.025,.975]), rtol=0, atol=1e-12)
            point = point_mse[c] - point_mse[a]
            np.testing.assert_allclose(point,row.MSE_improvement_beta_squared,rtol=0,atol=1e-13)
            p_values.append((1 + np.count_nonzero(np.abs(gains-point) >= abs(point)))/2001)
        p_values = np.array(p_values)
        order = np.argsort(p_values)
        q_values = np.ones(10)
        for position, index in enumerate(order):
            q_values[index] = min(1, min(p_values[order[j]]*10/(j+1) for j in range(position,10)))
        np.testing.assert_allclose(p_values,family.nominal_centered_bootstrap_MSE_p,rtol=0,atol=1e-14)
        np.testing.assert_allclose(q_values,family.exploratory_BH_q,rtol=0,atol=1e-14)
        supported = family.arm.str.startswith("adapter_") & family.MSE_improvement_beta_squared_low95.gt(0) & family.exploratory_BH_q.lt(.05)
        assert supported.sum() == 6
        outputs.append(dict(unit=name, resampling_units=len(units), supported_adapter_contrasts=family.loc[supported,"arm"].tolist()))
    report = dict(status="pass", rows=6834, arms=12, contrasts=20,
        checks="independent_saved_prediction_RMSE_rho_group_total_MSE_RMSE_intervals_centered_p_BH",
        rho_intervals="not_reconstructed; point_rho_checked", result=outputs, sources=hashes,
        interpretation="conditional_single_split_single_seed_development; five_chromosomes; not_native_model_superiority")
    (args.out / "checks.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep",type=Path,required=True)
    parser.add_argument("--review",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
