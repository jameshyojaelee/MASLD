from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))
SPEC = importlib.util.spec_from_file_location(
    "plan45_target_adjudication", SCRIPT_ROOT / "27_adjudicate_and_freeze_stage_a_targets.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def candidate(
    index: int,
    direction: str,
    stratum: str,
    passing: bool = True,
) -> dict[str, object]:
    return {
        "target_uid": f"target-{index}",
        "coarse_locus_uid": f"locus-{index}",
        "ensembl_id": f"ENSG{index}",
        "n_concordant_gwas": 2 if index == 1 else 1,
        "n_concordant_gwas_evidence_families": 2 if index == 1 else 1,
        "best_orientation_consensus": 0.95 - index / 100,
        "best_susie_pp4": 0.90 - index / 100,
        "selected_exact_shared_posterior": 0.80 - index / 100,
        "oriented_risk_effect": direction,
        "phenotype_strata_set": {stratum},
        "candidate_gate_pass": passing,
    }


class TargetAdjudicationTests(unittest.TestCase):
    def test_balanced_screen_requires_both_directions_and_strata(self) -> None:
        rows = [
            candidate(1, "risk_increases_expression", "direct_masld_mash_diagnosis"),
            candidate(2, "risk_decreases_expression", "mri_pdff_or_histologic_steatosis"),
            candidate(3, "risk_increases_expression", "mri_pdff_or_histologic_steatosis"),
            candidate(4, "risk_decreases_expression", "direct_masld_mash_diagnosis"),
            candidate(5, "risk_increases_expression", "direct_masld_mash_diagnosis"),
            candidate(6, "risk_decreases_expression", "mri_pdff_or_histologic_steatosis"),
            candidate(7, "risk_increases_expression", "direct_masld_mash_diagnosis"),
        ]
        mode, selected = MODULE.choose_targets(rows)
        self.assertEqual(mode, "balanced_four_to_six_locus_screen")
        self.assertEqual(len(selected), 6)
        self.assertEqual(len({row["coarse_locus_uid"] for row in selected}), 6)
        self.assertEqual(len({row["oriented_risk_effect"] for row in selected}), 2)

    def test_incomplete_architecture_falls_back_to_one_locus(self) -> None:
        rows = [
            candidate(i, "risk_increases_expression", "direct_masld_mash_diagnosis")
            for i in range(1, 6)
        ]
        mode, selected = MODULE.choose_targets(rows)
        self.assertEqual(mode, "single_locus_mechanism_only")
        self.assertEqual([row["target_uid"] for row in selected], ["target-1"])

    def test_no_passing_candidate_freezes_nothing(self) -> None:
        mode, selected = MODULE.choose_targets(
            [candidate(1, "risk_increases_expression", "direct_masld_mash_diagnosis", False)]
        )
        self.assertEqual(mode, "no_target_freeze")
        self.assertEqual(selected, [])


if __name__ == "__main__":
    unittest.main()
