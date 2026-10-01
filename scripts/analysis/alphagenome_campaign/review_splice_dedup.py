#!/usr/bin/env python3
"""Reconstruct the historical loose-splice comparison after exact row deduplication.

Original effect definitions, single-linkage blocks, seed, bootstrap and null
are retained. No new sequence inference, event definition or threshold tuning.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT/"GWAS/finemapping/results/alphagenome_atlas/p3a-splice-direction-20260915T141415Z/tables"
SEED, DRAWS = 123, 2000


def blocks(frame):
    result = np.empty(len(frame), dtype=object)
    for chrom, group in frame.groupby("chrom", sort=False):
        component, previous = 1, None
        for position, index in sorted((int(frame.iloc[j].position), int(j)) for j in group.index):
            if previous is not None and position-previous > 1_000_000:
                component += 1
            result[index] = f"chr{str(chrom).replace('chr', '')}:component{component:04d}"
            previous = position
    return result


def interval(values, unit):
    unique = pd.unique(unit)
    sums = np.array([values[unit == key].sum() for key in unique])
    counts = np.array([np.sum(unit == key) for key in unique])
    rng = np.random.default_rng(SEED)
    draws = []
    for _ in range(DRAWS):
        # Equivalent to the historical expanded-row block bootstrap.
        selected = rng.choice(len(unique), size=len(unique), replace=True)
        draws.append(sums[selected].sum()/counts[selected].sum())
    return np.quantile(draws, [.025, .975]), len(unique)


def direction(predicted, measured, unit):
    valid = np.isfinite(predicted) & np.isfinite(measured) & (predicted != 0)
    p, m, b = predicted[valid], measured[valid], unit[valid]
    agree = ((p > 0) == (m > 0)).astype(float)
    ci, nblocks = interval(agree, b)
    unique = pd.unique(b)
    balance = np.array([2*agree[b == key].sum()-np.sum(b == key) for key in unique])
    above, decided = int(np.sum(balance > 0)), int(np.sum(balance != 0))
    pp, pm = np.mean(p > 0), np.mean(m > 0)
    marginal = pp*pm+(1-pp)*(1-pm)
    return dict(n=len(p), n_blocks=nblocks, concordance=float(agree.mean()),
        marginal_expectation=float(marginal), above_marginal=float(agree.mean()-marginal),
        share_pred_positive=float(pp), share_meas_positive=float(pm),
        blocks_above_half=above, blocks_decided=decided,
        block_binomial_p=float(binomtest(above, decided, .5).pvalue) if decided else None,
        block_lo=float(ci[0]), block_hi=float(ci[1]))


def population(frame):
    scored = frame.loc[frame.state == "scored"].reset_index(drop=True)
    unit = blocks(scored)
    measured = scored.slope.to_numpy(float)
    directions = {name: direction(scored[name].to_numpy(float), measured, unit)
                  for name in ("excision_ratio_log2", "target_usage_log2", "site_usage_delta")}
    detection = {}
    for name in ("rank_abs_excision_ratio", "rank_abs_usage_log2", "rank_abs_atlas_quantile", "rank_proximity"):
        values = scored[name].to_numpy(float)
        valid = np.isfinite(values)
        ci, nb = interval(values[valid], unit[valid])
        detection[name] = dict(n=int(valid.sum()), n_blocks=nb, mean=float(values[valid].mean()),
                               block_lo=float(ci[0]), block_hi=float(ci[1]))
    paired = {}
    for comparator in ("atlas", "proximity"):
        other = "rank_abs_atlas_quantile" if comparator == "atlas" else "rank_proximity"
        a, b = scored.rank_abs_excision_ratio.to_numpy(float), scored[other].to_numpy(float)
        valid = np.isfinite(a) & np.isfinite(b)
        ci, nb = interval((a-b)[valid], unit[valid])
        paired[comparator] = dict(n=int(valid.sum()), n_blocks=nb, difference=float((a-b)[valid].mean()),
                                  block_lo=float(ci[0]), block_hi=float(ci[1]))
    magnitude = float(spearmanr(np.abs(measured), np.abs(scored.excision_ratio_log2.to_numpy(float))).statistic)
    ratio = directions["excision_ratio_log2"]
    return dict(total_rows=len(frame), states=frame.state.value_counts().to_dict(), scored_rows=len(scored),
        unique_scored_variant_event_pairs=len(scored[["variant_uid", "phenotype_id"]].drop_duplicates()),
        direction=directions, detection=detection, detection_model_minus_atlas=paired["atlas"],
        detection_model_minus_proximity=paired["proximity"],
        magnitude=dict(spearman_abs_slope_vs_abs_ratio=magnitude),
        original_P1_pass=bool(ratio["concordance"] > .65 and ratio["above_marginal"] > .05 and ratio["block_binomial_p"] < .01),
        below_original_D1_direction_cutoff=bool(ratio["concordance"] <= .60))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run statistical reconstruction on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    source_path = SOURCE/"splice_direction_pairs.tsv"
    with source_path.open() as handle:
        raw = list(csv.DictReader(handle, delimiter="\t"))
    keep, copies, seen = [], [], {}
    missing_identity_rows = []
    for index, row in enumerate(raw):
        if not row["phenotype_id"]:
            keep.append(index)
            missing_identity_rows.append(index)
            continue
        identity = (row["variant_uid"], row["phenotype_id"])
        if identity not in seen:
            seen[identity] = index
            keep.append(index)
        else:
            original = seen[identity]
            if row != raw[original]:
                raise ValueError(f"Conflicting repeated variant/event rows {original}, {index}; exact deduplication not justified")
            copies.append(dict(original_source_row0=original, duplicate_source_row0=index,
                variant_uid=identity[0], phenotype_id=identity[1], state=row["state"],
                all_deposited_fields_identical=True))
    pd.DataFrame(copies).to_csv(args.out/"exact_duplicate_rows.tsv", sep="\t", index=False)
    frame = pd.read_csv(source_path, sep="\t")
    historical = population(frame)
    archived = json.loads((SOURCE/"splice_direction.json").read_text())
    assert historical["total_rows"] == archived["n_pairs_carried"] == 657
    assert historical["scored_rows"] == archived["n_scored"] == 637
    for name, values in historical["direction"].items():
        for field, value in values.items():
            np.testing.assert_allclose(value, archived["direction"][name][field], rtol=0, atol=1e-12)
    for name, values in historical["detection"].items():
        for field, value in values.items():
            np.testing.assert_allclose(value, archived["detection"][name][field], rtol=0, atol=1e-12)
    for comparator in ("atlas", "proximity"):
        key = "detection_model_minus_"+comparator
        for field, value in historical[key].items():
            np.testing.assert_allclose(value, archived[key][field], rtol=0, atol=1e-12)
    np.testing.assert_allclose(historical["magnitude"]["spearman_abs_slope_vs_abs_ratio"],
                               archived["magnitude"]["spearman_abs_slope_vs_abs_ratio"], rtol=0, atol=1e-12)
    unique = frame.iloc[keep].copy().reset_index(drop=True)
    deduplicated = population(unique)
    # A second exact duplicate removal must change nothing.
    identified = unique.phenotype_id.notna()
    assert not unique.loc[identified].duplicated(["variant_uid", "phenotype_id"]).any()
    result = dict(status="pass", exact_extra_copies_removed=len(copies),
        missing_event_identity_rows_preserved=len(missing_identity_rows),
        historical_637_rows_reconstructed=historical, exact_unique_rows=deduplicated,
        statistical_procedure=dict(bootstrap_draws=DRAWS, seed=SEED,
            blocks="original_1Mb_single_linkage_components_by_variant_position",
            null="two_sided_exact_binomial_on_block_majorities_against_one_half; tied_blocks_omitted",
            marginal="printed_row_matched_sign_skew_expectation_not_the_binomial_null",
            multiplicity="original_nominal_p_values_retained; no_new_confirmatory_family_or_threshold_tuning"),
        interpretation_limits=[
            "Historical loose donor_or_acceptor_neighbor definition is not the exact source cluster.",
            "Model REF_intersection_ALT proposed_junctions omit gained_lost union and ignore strand.",
            "GTEx normalized intron_excision phenotype slope per ALT dosage is not deltaPSI; predicted log2 ratio is a distinct unit.",
            "GTEx_v8 overlaps foundation_model_training; no_external_generalization_claim.",
            "Removing exact duplicate records corrects weighting, not independent participant replication.",
            "P1 requires concordance>.65,p<.01; written D1 uses concordance<=.60,p>=.05; inconsistency remains recorded.",
            "Below original ranking criteria does not mean zero directional information."],
        historical_outputs_unchanged=True, new_sequence_inference=False, new_event_tuning=False,
        protected_outcomes_read=False, source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest())
    (args.out/"comparison.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
