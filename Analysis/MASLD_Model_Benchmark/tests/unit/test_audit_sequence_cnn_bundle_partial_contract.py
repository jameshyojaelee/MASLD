from __future__ import annotations

import unittest

from scripts.audit_sequence_cnn_bundle_partial_contract import (
    SequenceCNNContractError,
    validate_expected_units,
)


class SequenceCNNBundlePartialContractTests(unittest.TestCase):
    def _units(self) -> list[dict[str, object]]:
        return [
            {
                "logical_work_id": f"seq__sequence_cnn_control__lineage_{index}__d{index % 5}g{index % 5}__s{20260824 + index}",
                "lineage_id": f"lineage_{index}",
                "outer_fold": index % 5,
                "split_id": f"donor{index % 5}_genomic{index % 5}",
                "seed": 20260824 + index,
            }
            for index in range(8)
        ]

    def test_exact_eight_unit_roster_preserves_receipt_seeds(self) -> None:
        units = self._units()
        result = validate_expected_units(units)
        self.assertEqual([row["seed"] for row in result], [20260824 + index for index in range(8)])

    def test_incomplete_coverage_is_rejected(self) -> None:
        with self.assertRaises(SequenceCNNContractError):
            validate_expected_units(self._units()[:7])

    def test_non_diagonal_donor_genomic_split_is_rejected(self) -> None:
        units = self._units()
        units[0]["split_id"] = "donor0_genomic1"
        with self.assertRaises(SequenceCNNContractError):
            validate_expected_units(units)


if __name__ == "__main__":
    unittest.main()
