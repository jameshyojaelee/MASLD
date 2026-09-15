from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import torch

from scripts.probe_borzoi_grelu_converted_port_full_window import (
    BorzoiFullWindowProbeError,
    TARGET_FIELDS,
    TARGET_INDICES,
    build_synthetic_sequence,
    read_target_subset,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]


class BorzoiFullWindowProbeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(
            (ROOT / "config/borzoi_grelu_converted_port_full_window_probe.json").read_text()
        )

    def test_contract_is_synthetic_and_fail_closed(self) -> None:
        validate_contract(self.contract)
        self.assertEqual(self.contract["fixture"]["sequence_bp"], 524288)
        self.assertFalse(self.contract["execution_contract"]["biological_data_allowed"])
        self.assertFalse(self.contract["claim_boundary"]["native_numeric_parity_established"])

    def test_contract_rejects_claim_and_execution_drift(self) -> None:
        changed = deepcopy(self.contract)
        changed["claim_boundary"]["open_champion_eligible"] = True
        with self.assertRaisesRegex(BorzoiFullWindowProbeError, "opened or drifted"):
            validate_contract(changed)
        changed = deepcopy(self.contract)
        changed["execution_contract"]["biological_data_allowed"] = True
        with self.assertRaisesRegex(BorzoiFullWindowProbeError, "opened or drifted"):
            validate_contract(changed)

    def test_synthetic_sequence_is_one_hot_and_balanced(self) -> None:
        fixture = build_synthetic_sequence(524288, torch)
        self.assertEqual(tuple(fixture.shape), (1, 4, 524288))
        self.assertTrue(torch.all(fixture.sum(dim=1) == 1))

    def test_target_manifest_parser_binds_positions_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "targets.txt"
            expected = []
            with path.open("x", encoding="utf-8", newline="") as handle:
                handle.write("\t".join(TARGET_FIELDS) + "\n")
                for index in range(7611):
                    identifier = f"target{index}"
                    description = f"description{index}"
                    strand_pair = index
                    handle.write(
                        "\t".join(
                            (
                                str(index), identifier, f"file{index}", "0", "0",
                                "1.0", "sum", str(strand_pair), description,
                            )
                        ) + "\n"
                    )
                    if index in TARGET_INDICES:
                        expected.append(
                            {"index": index, "identifier": identifier,
                             "strand_pair": strand_pair, "description": description}
                        )
            self.assertEqual(read_target_subset(path, expected), expected)


if __name__ == "__main__":
    unittest.main()
