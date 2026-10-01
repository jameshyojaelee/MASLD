#!/usr/bin/env python3
"""Check source specialist variant-target labels before final comparisons."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--verify-production", action="store_true")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run scientific identity check on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    tier = ROOT/"GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z/tables/endpoint1_matched_tierA4.tsv.gz"
    c2 = ROOT/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
    original = pd.read_csv(tier, sep="\t", usecols=["key", "peak_id", "beta_alt", "heldout_fold", "model_peak_overlap"])
    selected = original.loc[original.model_peak_overlap.eq(1)].copy()
    labels = pd.read_csv(c2, sep="\t", usecols=["lead_variant_id", "peak_id", "beta_alt", "heldout_fold"])
    labels["key"] = labels.lead_variant_id.str.removeprefix("chr")
    assert not labels.key.duplicated().any() and not selected.key.duplicated().any()
    checked = selected.merge(labels, on="key", how="left", suffixes=("_specialist", "_C2"), validate="one_to_one", indicator=True)
    checked["exact_target"] = checked.peak_id_specialist.eq(checked.peak_id_C2)
    checked["same_beta"] = np.isclose(checked.beta_alt_specialist, checked.beta_alt_C2, rtol=1e-7, atol=1e-8)
    checked["same_fold"] = checked.heldout_fold_specialist.eq(checked.heldout_fold_C2)
    present = checked._merge.eq("both")
    common_target = present & checked.exact_target
    eligible = common_target & checked.same_beta & checked.same_fold
    mismatch = common_target & ~eligible
    disagreement = checked.loc[mismatch]
    disagreement.to_csv(args.out/"source_disagreements.tsv", sep="\t", index=False)
    checked.loc[~present].to_csv(args.out/"C2_eligibility_exclusions.tsv", sep="\t", index=False)
    checked.loc[present & ~checked.exact_target].to_csv(args.out/"different_measured_target_exclusions.tsv", sep="\t", index=False)
    result = dict(status="fail" if mismatch.any() else "pass", all_TierA4_rows=len(original),
        source_model_peak_overlap1_rows=len(selected), exact_variant_target_beta_fold_rows=int(eligible.sum()),
        C2_population_eligibility_excluded_rows=int((~present).sum()),
        eligibility_exclusion_reason="variant_absent_from_original_C2_eligible_population; not_a_label_mismatch",
        different_measured_target_excluded_rows=int((present & ~checked.exact_target).sum()),
        target_exclusion_reason="same_variant_associated_with_a_different_source_peak; distinct_valid_molecular_targets_not_source_mislabeling",
        same_target_beta_mismatch_rows=int((common_target & ~checked.same_beta).sum()),
        same_target_fold_mismatch_rows=int((common_target & ~checked.same_fold).sum()),
        source_paths=[str(tier), str(c2)], effect_signs_changed=False, protected_outcomes_read=False)
    if args.verify_production:
        from model_final import verified_specialists
        producer_out = args.out/"production_target_guard"
        producer_out.mkdir()
        admitted = verified_specialists(producer_out)
        expected = checked.loc[eligible].set_index("key").sort_index()
        actual = admitted.set_index("key").sort_index()
        np.testing.assert_array_equal(actual.index, expected.index)
        np.testing.assert_array_equal(actual.peak_id, expected.peak_id_specialist)
        np.testing.assert_allclose(actual.beta_alt, expected.beta_alt_specialist, atol=1e-12, rtol=0)
        np.testing.assert_array_equal(actual.heldout_fold, expected.heldout_fold_specialist)
        exclusions = pd.read_csv(producer_out/"specialist_other_target_or_population_exclusions.tsv", sep="\t")
        assert set(exclusions.key) == set(checked.loc[~common_target, "key"])
        result["production_population_and_exclusions_independently_verified"] = True
    (args.out/"checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    if mismatch.any():
        raise SystemExit("Specialist target identity differs; variant-only comparison is invalid")


if __name__ == "__main__":
    main()
