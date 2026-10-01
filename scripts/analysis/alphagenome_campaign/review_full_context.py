#!/usr/bin/env python3
"""Independent source identity and numerical review of full context comparisons."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from review_model_results import review_reporter, review_figures, review_sweep

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915"


def weights(elements):
    _, codes, counts = np.unique(elements, return_inverse=True, return_counts=True)
    return 1/counts[codes]


def check_transfer(directory, features, output):
    output.mkdir()
    receipt = json.loads((directory/"complete.json").read_text())
    frame = pd.read_csv(directory/"endogenous_oof.tsv.gz", sep="\t")
    assert not frame.duplicated(["element_id", "peak"]).any()
    source = pd.read_csv(BASE/"data/endpoints-revised-21772544/reporter_caqtl_all_nominal_pairs.tsv.gz", sep="\t")
    source = source.loc[source.eligible.astype(str).str.lower().eq("true")]
    assert not source.duplicated(["element_id", "peak"]).any()
    lookup = source.set_index(["element_id", "peak"])
    keys = pd.MultiIndex.from_frame(frame[["element_id", "peak"]])
    original = lookup.loc[keys].reset_index()
    direct = original.ref.eq(original.genomic_ref) & original.alt.eq(original.genomic_alt)
    swapped = original.alt.eq(original.genomic_ref) & original.ref.eq(original.genomic_alt)
    assert (direct | swapped).all()
    signs = np.where(direct, 1, -1)
    y = original.beta.to_numpy()*signs
    np.testing.assert_allclose(frame.effect, y, rtol=1e-10, atol=1e-12)
    np.testing.assert_array_equal(frame.orientation_multiplier, signs)
    np.testing.assert_array_equal(frame.fold, original.outer_fold)
    np.testing.assert_array_equal(frame.borzoi_long_range_group_id, original.borzoi_long_range_group_id)
    assert frame.groupby("borzoi_long_range_group_id").fold.nunique().max() == 1
    assert original.effect_unit.eq("source_normalized_accessibility_slope_per_reporter_ALT_dosage").all()
    manifest = pd.read_csv(features/"reporter_manifest.tsv", sep="\t")
    assert not manifest.duplicated(["contig", "variant_pos1", "genomic_ref", "genomic_alt"]).any()
    target_labels = pd.read_csv(features/"target_labels.tsv.gz", sep="\t")
    with np.load(features/"reporters_16384.npz", allow_pickle=False) as archive:
        np.testing.assert_array_equal(archive["element_ids"], manifest.element_id)
        np.testing.assert_array_equal(archive["allele_order"], ["REF", "ALT", "REF_RC", "ALT_RC"])
        sequence_available = pd.Series(archive["extracted"], index=manifest.element_id)
    with np.load(features/"target_pools_16384.npz", allow_pickle=False) as archive:
        np.testing.assert_array_equal(archive["pair_index"], target_labels.pair_index)
        target_available = archive["extracted"] & sequence_available.loc[target_labels.element_id].to_numpy()
    expected = target_labels.loc[target_available, ["element_id", "peak"]]
    assert pd.MultiIndex.from_frame(expected).equals(keys)
    # Reconstruct physical correspondence directly from coordinate boundaries.
    position = original.variant_pos1.to_numpy()-1
    start, end = original.peak_start0.to_numpy(), original.peak_end0.to_numpy()
    inside = (position >= start) & (position < end)
    overlap = np.maximum(0, np.minimum(original.reference_match_end0, end)-np.maximum(original.reference_match_start0, start))
    np.testing.assert_array_equal(original.physical_variant_peak_overlap, inside)
    np.testing.assert_array_equal(original.physical_construct_peak_overlap_bp, overlap)
    assert ((start >= position-8192) & (end <= position+8192) & (start < end)).all()
    distance = np.abs(original.distance_from_peakCenter.to_numpy())
    distance_names = ("0-500", "501-1000", "1001-5000", "5001-10000", "10001-100000", ">100000")
    distance_masks = ((distance <= 500), (distance > 500)&(distance <= 1000), (distance > 1000)&(distance <= 5000),
        (distance > 5000)&(distance <= 10000), (distance > 10000)&(distance <= 100000), distance > 100000)
    for name, mask in zip(distance_names, distance_masks):
        np.testing.assert_array_equal(original.distance_stratum.eq(name), mask)
    activity = original.LX2_control_mean_log2_activity.to_numpy()
    activity_masks = ((activity <= 0), (activity > 0), ~np.isfinite(activity))
    activity_names = ("RNA_over_DNA_le1", "RNA_over_DNA_gt1", "missing")
    for name, mask in zip(activity_names, activity_masks):
        np.testing.assert_array_equal(original.activity_stratum.eq(name), mask)
    strata = {"all_fully_observed": np.ones(len(frame), bool), "variant_inside_peak": inside, "construct_overlaps_peak": overlap > 0}
    strata.update({"distance_"+name: mask for name, mask in zip(distance_names, distance_masks)})
    strata.update({"activity_"+name: mask for name, mask in zip(activity_names, activity_masks)})
    recorded = pd.read_csv(directory/"performance_by_fixed_stratum.tsv", sep="\t").set_index("stratum")
    assert set(recorded.index) == set(strata) and len(recorded) == 12
    elements = original.element_id.to_numpy()
    global_weights = weights(elements)
    np.testing.assert_allclose(frame.variant_weight, global_weights, atol=1e-14, rtol=0)
    # Independently recover the separate reporter-unit means and training RMS.
    reporters = pd.read_csv(directory/"reporter_oof.tsv.gz", sep="\t").set_index("element_id")
    raw_reporters = pd.read_csv(BASE/"data/endpoints-revised-21772544/reporter_treatment_replicates.tsv.gz", sep="\t")
    raw_reporters = raw_reporters.loc[raw_reporters.cell_line.eq("LX2") & raw_reporters.eligible.astype(str).str.lower().eq("true")]
    replicates = raw_reporters.pivot(index="element_id", columns="experimental_replicate", values="effect_control").loc[reporters.index]
    assert replicates.shape[1] == 4
    np.testing.assert_allclose(reporters.measured_reporter_log2_effect, replicates.mean(axis=1), rtol=1e-12, atol=1e-12)
    expected_reporter = sequence_available.loc[reporters.index].to_numpy() & np.isfinite(replicates.to_numpy()).all(axis=1)
    np.testing.assert_array_equal(reporters.eligible, expected_reporter)
    choices = pd.read_csv(directory/"inner_selection.tsv", sep="\t")
    fold = original.outer_fold.to_numpy()
    for choice in choices.itertuples():
        assert choice.inner_fold == (choice.held_fold+1)%5
        train = fold != choice.held_fold
        rt = reporters.fold.ne(choice.held_fold) & reporters.eligible
        assert choice.training_pairs == train.sum() and choice.training_reporters == rt.sum()
        np.testing.assert_allclose(choice.training_caQTL_RMS, np.sqrt(np.average(y[train]**2, weights=global_weights[train])), atol=1e-12, rtol=0)
        np.testing.assert_allclose(choice.training_reporter_RMS, np.sqrt(np.mean(reporters.loc[rt, "measured_reporter_log2_effect"]**2)), atol=1e-12, rtol=0)
    saved_margins = pd.read_csv(directory/"training_only_operational_margins.tsv", sep="\t", index_col=0)
    block = original.borzoi_long_range_group_id.to_numpy()
    a, b = frame.endogenous_specialist.to_numpy(), frame.shared_private.to_numpy()
    assert np.isfinite(np.r_[a,b]).all()
    ea, eb = (a-y)**2, (b-y)**2
    rng = np.random.default_rng(20260915)
    p_values, records = [], []
    for name, mask in strata.items():
        row = recorded.loc[name]
        w = np.zeros(len(frame)); w[mask] = weights(elements[mask])
        per_fold_blocks = [np.unique(block[mask & (fold == f)]) for f in range(5)]
        support = min(len(x) for x in per_fold_blocks)
        assert row.pairs == mask.sum() and row.variants == len(set(elements[mask]))
        assert row.blocks == len(set(block[mask])) and row.minimum_blocks_per_fold == support
        margin_fold = np.full(5, np.nan)
        for held in range(5):
            train = mask & (fold != held)
            if train.any():
                margin_fold[held] = .05*np.average(y[train]**2, weights=w[train])
        np.testing.assert_allclose(saved_margins[name], margin_fold, atol=1e-12, rtol=0, equal_nan=True)
        if not mask.any() or support < 2:
            assert row.status == "insufficient_fold_support_or_outside_observable_targets"
            assert row.p_for_family_accounting == 1
            p_values.append(1.0)
            records.append(dict(stratum=name, pairs=int(mask.sum()), status="not_testable_preserved_in_12_family"))
            continue
        point = float(np.average(ea[mask]-eb[mask], weights=w[mask]))
        margin = float(np.average(margin_fold[fold[mask]], weights=w[mask]))
        unique = np.unique(block[mask]); code = {key: i for i,key in enumerate(unique)}
        total_weight = np.array([w[mask & (block == key)].sum() for key in unique])
        error_sum = np.array([np.sum(w[mask & (block == key)]*(ea-eb)[mask & (block == key)]) for key in unique])
        fold_codes = [np.array([code[key] for key in keys]) for keys in per_fold_blocks]
        draws = np.empty(receipt["bootstrap"])
        for draw in range(len(draws)):
            chosen = np.concatenate([rng.choice(indices, len(indices), replace=True) for indices in fold_codes])
            draws[draw] = error_sum[chosen].sum()/total_weight[chosen].sum()
        ci = np.quantile(draws, [.025,.975])
        p = float((1+sum(np.abs(draws-point) >= abs(point)))/(len(draws)+1))
        np.testing.assert_allclose([row.MSE_improvement,row.low95,row.high95,row.training_only_operational_margin], [point,*ci,margin], atol=1e-11, rtol=0)
        np.testing.assert_allclose([row.specialist_RMSE_beta,row.shared_private_RMSE_beta],
            [np.sqrt(np.average(ea[mask],weights=w[mask])),np.sqrt(np.average(eb[mask],weights=w[mask]))], atol=1e-12, rtol=0)
        assert abs(row.p_for_family_accounting-p) <= 1/(len(draws)+1)+1e-12
        p_values.append(float(row.p_for_family_accounting))
        records.append(dict(stratum=name, pairs=int(mask.sum()), variants=len(set(elements[mask])), blocks=len(unique),
            MSE_improvement=point, low95=ci[0], high95=ci[1], conditional_margin=margin,
            identical_predictions=bool(np.array_equal(a[mask],b[mask])),
            both_all_zero=bool(np.all(a[mask] == 0) and np.all(b[mask] == 0)),
            status="conditional_identical_predictions_not_information_limit" if np.array_equal(a[mask],b[mask]) else "descriptive_fixed_prediction_comparison"))
    # Independently evaluate the BH definition over the complete 12-member set.
    order = np.argsort(p_values)
    sorted_p = np.asarray(p_values)[order]
    q = np.empty(12)
    for index, original_index in enumerate(order):
        q[original_index] = min(1.0, min(12*sorted_p[j]/(j+1) for j in range(index,12)))
    np.testing.assert_allclose(recorded.loc[list(strata), "BH_q_complete_fixed_stratum_family"], q, atol=1e-12, rtol=0)
    if "endogenous_selected_with_null" in frame:
        old = pd.read_csv(BASE/"model/transfer_21773114/endogenous_oof.tsv.gz", sep="\t")
        np.testing.assert_array_equal(frame[["element_id","peak"]], old[["element_id","peak"]])
        np.testing.assert_array_equal(frame.endogenous_selected_with_null, old.endogenous_specialist)
        np.testing.assert_array_equal(frame.shared_private_selected_with_null, old.shared_private)
        assert choices.specialist_alpha.notna().all() and choices.joint_alpha.notna().all()
    pd.DataFrame(records).to_csv(output/"strata_reconstructed.tsv",sep="\t",index=False)
    return dict(status="pass", source_pairs=len(source), whole_target_eligible_pairs=len(target_labels), evaluated_pairs=len(frame),
        variants=len(set(elements)), locus_blocks=len(set(block)), family_members=12, estimable=sum(r["status"] != "not_testable_preserved_in_12_family" for r in records),
        physical_variant_peak_rows=int(inside.sum()), physical_construct_overlap_rows=int((overlap>0).sum()),
        signed_source_beta_and_target_identity=True, separate_training_unit_scales=True,
        source_effects_and_reporter_measurement_uncertainty="conditional_fixed_labels; not_propagated",
        identical_returned_models=bool(np.array_equal(a,b)), all_zero_returned_models=bool(np.all(a == 0) and np.all(b == 0)))


def check_sweep_summary(directory):
    source = pd.read_csv(directory/"matched_predictions.tsv.gz", sep="\t")
    metrics = pd.read_csv(directory/"performance.tsv", sep="\t").set_index("arm")
    assert len(source) == 6834 and not source.variant_id.duplicated().any()
    assert len(metrics) == 14
    y = source.observed_beta.to_numpy()
    for arm, row in metrics.iterrows():
        predicted = source[arm].to_numpy()
        assert np.isfinite(predicted).all()
        np.testing.assert_allclose(row.RMSE_beta_units, np.sqrt(np.mean((predicted-y)**2)), atol=1e-12,rtol=0)
        np.testing.assert_allclose(row.MAE_beta_units, np.mean(np.abs(predicted-y)), atol=1e-12,rtol=0)
        rho = float(spearmanr(y,predicted).statistic) if np.ptp(predicted) else np.nan
        np.testing.assert_allclose(row.signed_spearman,rho,atol=1e-12,rtol=0,equal_nan=True)
    return dict(arms=14, rows=6834, native_errors_and_ranks_reconstructed=True)


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run scientific review on a compute node")
    args.out.mkdir(parents=True,exist_ok=False)
    result = dict(status="pass", protected_outcomes_read=False, producer_metric_functions_imported=False)
    result["full_B1"] = review_reporter(args.lx2,args.out)
    result["figures"] = review_figures(args.figures,args.lx2,args.donor)
    result["adaptation"] = review_sweep(args.sweep,args.out,True)
    result["sweep_summary"] = check_sweep_summary(args.sweep_summary)
    result["transfer"] = {directory.name:check_transfer(directory,args.features,args.out/directory.name) for directory in args.transfer}
    (args.out/"checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lx2",type=Path,required=True)
    parser.add_argument("--figures",type=Path,required=True)
    parser.add_argument("--donor",type=Path,required=True)
    parser.add_argument("--sweep",type=Path,required=True)
    parser.add_argument("--sweep-summary",type=Path,required=True)
    parser.add_argument("--transfer",type=Path,action="append",required=True)
    parser.add_argument("--features",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
