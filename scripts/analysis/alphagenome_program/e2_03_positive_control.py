#!/usr/bin/env python3
"""E2 step 3: a positive control that bounds what the transfer design could have detected.

A null transfer result is only informative if the same pipeline would have found transfer had it
been there.  This script keeps every real component - the same 4,359-element HyenaDNA embeddings,
the same chromosome-group folds, the same 1-Mb blocks, the same head recipe, the same five seeds and
the same bootstrap - and replaces only the two labels with simulated ones that share a known,
linearly learnable sequence component:

    signal   s = standardised (mean of the REF/ALT delta embedding) @ w,  w fixed at seed 20260914
    reporter y_rep  = s + noise_level * e_rep
    endogenous y_endo = s + noise_level * e_endo

Both labels are generated from the same `s`, so a perfect transfer opportunity exists by
construction.  The recovered gain of arm B over arm A at each noise level is the design's
sensitivity.  Noise levels are chosen so that arm A's macro Spearman brackets whatever the real run
produced; the comparison is made in RESULTS.md, not here.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from c2_03_fit_heads import allele_identity_features, fast_spearman, interval  # noqa: E402
from e2_02_fit_transfer import (  # noqa: E402
    SEEDS,
    ProjectionCache,
    TransferError,
    block_bootstrap,
    bootstrap_p,
    build_reporter_oof,
    macro_over_folds,
    run_campaign,
)

SIGNAL_SEED = 20260914
PRIMARY_STRATUM_BP = 5000.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resamples", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260914)
    parser.add_argument("--noise-levels", type=float, nargs="+", default=(1.0, 2.0, 4.0))
    arguments = parser.parse_args()
    out = arguments.output
    (out / "tables").mkdir(parents=True, exist_ok=True)

    units = pd.read_csv(arguments.inputs / "e2_units.tsv.gz", sep="\t")
    units = units.sort_values("row_index").reset_index(drop=True)
    n = len(units)
    folds = units.chrom_fold.to_numpy(np.int64)
    blocks = units.block_1mb.to_numpy(dtype=str)
    absolute = units.abs_distance_nearest.to_numpy(np.float64)
    pool = units.has_endogenous.to_numpy(bool)

    with np.load(arguments.raw / "hyenadna_allele_embeddings.npz", allow_pickle=False) as data:
        if not np.array_equal(data["row_index"], units.row_index.to_numpy(np.int64)):
            raise TransferError("embedding row alignment differs")
        embeddings = data["embeddings"]

    raw_delta = (
        (embeddings[:, 1].astype(np.float64) + embeddings[:, 3].astype(np.float64)) / 2.0
        - (embeddings[:, 0].astype(np.float64) + embeddings[:, 2].astype(np.float64)) / 2.0
    )
    generator = np.random.default_rng(SIGNAL_SEED)
    weights = generator.normal(size=raw_delta.shape[1])
    signal = raw_delta @ weights
    signal = (signal - signal.mean()) / signal.std()

    cache = ProjectionCache(embeddings)
    identity = allele_identity_features(
        units.mpra_ref.to_numpy(dtype=str), units.mpra_alt.to_numpy(dtype=str)
    )

    rows_primary = np.flatnonzero(pool & (absolute <= PRIMARY_STRATUM_BP))
    summary: list[dict[str, object]] = []
    for level in arguments.noise_levels:
        noise = np.random.default_rng([SIGNAL_SEED, int(level * 1000)])
        reporter_y = signal + level * noise.normal(size=n)
        endogenous_y = signal + level * noise.normal(size=n)
        audit: list[dict[str, object]] = []
        reporter_oof: dict[int, dict[int, np.ndarray]] = {}
        self_eval: dict[int, np.ndarray] = {}
        for seed in SEEDS:
            per_eval, own = build_reporter_oof(
                cache=cache, reporter_y=reporter_y, reporter_pool=np.ones(n, bool),
                folds=folds, blocks=blocks, seed=seed, audit=audit,
            )
            reporter_oof[seed] = per_eval
            self_eval[seed] = own
        reporter_ensemble = np.mean(np.vstack([self_eval[s] for s in SEEDS]), axis=0)
        per_seed = run_campaign(
            cache=cache, identity=identity, endo_y=endogenous_y, endo_pool=pool,
            reporter_oof=reporter_oof, folds=folds, blocks=blocks, audit=audit,
            tag=f"positive_control_noise_{level}",
        )
        ensemble = {
            arm: np.mean(np.vstack([values[s] for s in SEEDS]), axis=0)
            for arm, values in per_seed.items()
        }
        alone = np.full(n, np.nan)
        for held in range(5):
            rows = np.flatnonzero(pool & (folds == held))
            alone[rows] = np.mean(
                np.vstack([reporter_oof[s][held][rows] for s in SEEDS]), axis=0
            )
        ensemble["reporter_oof_alone"] = alone
        arms = {name: values[rows_primary] for name, values in ensemble.items()}
        observed = endogenous_y[rows_primary]
        samples = block_bootstrap(
            observed, arms, folds[rows_primary], blocks[rows_primary],
            resamples=arguments.resamples, seed=arguments.bootstrap_seed,
        )
        reference = samples["endogenous_only"]
        base_point = macro_over_folds(observed, arms["endogenous_only"], folds[rows_primary])
        for name, values in arms.items():
            point = macro_over_folds(observed, values, folds[rows_primary])
            gain = samples[name] - reference
            low, high, _ = interval(gain)
            summary.append({
                "noise_level": level,
                "arm": name,
                "evaluation_variants": int(rows_primary.size),
                "evaluation_blocks_1mb": int(len(set(blocks[rows_primary].tolist()))),
                "macro_spearman": point,
                "gain_vs_endogenous_only": point - base_point,
                "gain_ci_low": low,
                "gain_ci_high": high,
                "gain_bootstrap_p_two_sided": bootstrap_p(gain) if name != "endogenous_only" else float("nan"),
                "reporter_head_macro_on_its_own_label": macro_over_folds(
                    reporter_y, reporter_ensemble, folds
                ),
                "reporter_head_pooled_on_its_own_label": fast_spearman(reporter_y, reporter_ensemble),
            })
        print(json.dumps(summary[-len(arms):], indent=2, sort_keys=True), flush=True)

    frame = pd.DataFrame(summary)
    frame.to_csv(out / "tables" / "positive_control.tsv", sep="\t", index=False)
    print(frame.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
