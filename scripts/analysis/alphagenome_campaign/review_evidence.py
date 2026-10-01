#!/usr/bin/env python3
"""Independently read deposited estimates and expose interpretation limits."""
import argparse
import json
import os
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("run executable review checks on a compute node")
    args.output.mkdir(parents=True, exist_ok=False)
    ag = ROOT / "GWAS/finemapping/results/alphagenome_program"
    atlas = ROOT / "GWAS/finemapping/results/alphagenome_atlas"
    zero = pd.read_csv(ag / "ad1-gate-20260915T114615Z/tables/zeroshot_arms.tsv", sep="\t")
    local = zero.set_index("arm").loc["alphagenome_local_zeroshot_atac"]
    matched = zero.set_index("arm").loc["hyenadna_delta_ridge"]
    assert int(local.variants) == int(matched.variants) == 3627
    assert abs(float(local.macro_spearman) - 0.7441117690812528) < 1e-12
    splice_path = atlas / "p3a-splice-direction-20260915T141415Z/tables/splice_direction.json"
    splice = json.loads(splice_path.read_text())
    direction = splice["direction"]["excision_ratio_log2"]
    assert direction["block_lo"] > direction["marginal_expectation"]
    assert direction["concordance"] < 0.60
    reliability_path = ag / "a1-mpra-reliability-v2-20260914T204500Z/tables/reliability_by_context.tsv"
    rel = pd.read_csv(reliability_path, sep="\t")
    relrow = rel[(rel.context == "HepG2_PAOA") & (rel.quantity == "allele_effect") &
                 (rel.stratum == "all") & (rel.partition == "mean_of_three")].iloc[0]
    assert relrow.pearson_lo < 0.70 < relrow.pearson_hi
    reporter_path = ag / "ad-arm2-mpra-20260915T114407Z/tables/paired_contrasts.tsv"
    reporter = pd.read_csv(reporter_path, sep="\t")
    reporter = reporter[(reporter.scope == "own") & (reporter.bootstrap_seed == 20260914) &
                        (reporter.against == "hyenadna")].set_index("cell_id")
    primary, secondary = reporter.loc["len2048_poolregion"], reporter.loc["len2048_centrebin"]
    assert primary.paired_ci_low < 0 < primary.paired_ci_high
    assert secondary.paired_ci_low > 0
    # Preserve the entire input header and selected source row; differences in
    # precision or future changes must fail instead of silently reinterpreting.
    text = (ag / "a1-mpra-reliability-v2-20260914T204500Z/RESULTS.md").read_text()
    assert "0.6938 [0.6356, 0.7383]" in text
    assert 0.6356 < 0.70 < 0.7383
    source_script = ROOT / "Analysis/MASLD_Model_Benchmark/executions/corgi-film-h3k27ac-finetune-20260907T191151Z/score_film.py"
    scorer = source_script.read_text()
    assert "al = num / den" in scorer and "Pc, Mc = centre(Pm[:, cols]), centre(Mm[:, cols])" in scorer
    report = {
        "state": "interpretation_corrections_no_named_release_modified",
        "currin": {"macro_spearman":float(local.macro_spearman), "n":3627,"blocks":1676,
                   "matched_hyenadna":float(matched.macro_spearman),"full_universe_n":32322,
                   "full_population_alphagenome_performance":"not_established_by_this_subset"},
        "reporter": {"primary_macro":float(primary.alphagenome_score),"comparator_macro":float(primary.reference_score),
                     "primary_difference":float(primary.paired_difference),"primary_ci95":[float(primary.paired_ci_low),float(primary.paired_ci_high)],
                     "secondary_variant_bin_macro":float(secondary.alphagenome_score),"secondary_difference":float(secondary.paired_difference),
                     "secondary_ci95":[float(secondary.paired_ci_low),float(secondary.paired_ci_high)],
                     "interpretation":"development_signal_primary_superiority_not_established","source":str(reporter_path)},
        "splice": {"source":str(splice_path),"direction":direction,
                   "magnitude_spearman":splice["magnitude"]["spearman_abs_slope_vs_abs_ratio"],
                   "interpretation":"weak_direction_below_original_0.60_ranking_criterion; retain_failed_detection_and_threshold_inconsistency",
                   "old_decision_text":splice["verdicts"]["decision"]},
        "reliability": {"source":str(reliability_path),"deposited_columns":rel.columns.tolist(),
                       "HepG2_PAOA_point":float(relrow.pearson),"ci95":[float(relrow.pearson_lo),float(relrow.pearson_hi)],
                       "interval_contains_0.70":True,"point_rule_failed":True},
        "donor": {"source":str(source_script),
                  "problem":"reported_skill_fits_scalar_to_evaluation_outcomes_after_evaluation_centering",
                  "comparator_problem":"RNA_lane_uses_all_output_regions_and_is_not_region_held",
                  "action":"new_comparisons_fit_native_logCPM_from_training_donors_and_training_regions_only"},
        "source_use": {"sequence":"counted_interval_permitted_with_1bp_residual; called_peak_correspondence_blocked",
                       "weights":"derivative_weights_permitted_by_2026_09_07_owner_decision",
                       "data":"raw_and_processed_source_data_not_rehosted"},
        "protected_outcomes_read":False,
    }
    (args.output/"evidence_corrections.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
