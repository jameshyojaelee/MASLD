from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

from scripts.acquire_scooby_epicardioids_checkpoint import (
    LocalRangeReader,
    ScoobyCheckpointAcquisitionError,
    compare_header,
    validate_task_boundary,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/scooby_epicardioids_checkpoint_acquisition.json"


class ScoobyCheckpointAcquisitionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_acquisition_preserves_observed_atac_boundary(self) -> None:
        validate_task_boundary(self.config)
        mutated = copy.deepcopy(self.config)
        mutated["task_scope"]["rna_conditioned_atac_eligible"] = True
        with self.assertRaises(ScoobyCheckpointAcquisitionError):
            validate_task_boundary(mutated)

    def test_acquisition_cannot_deserialize_or_forward(self) -> None:
        mutated = copy.deepcopy(self.config)
        mutated["gate_policy"]["checkpoint_deserialization_allowed"] = True
        with self.assertRaises(ScoobyCheckpointAcquisitionError):
            validate_task_boundary(mutated)

    def test_header_comparison_rejects_tensor_change(self) -> None:
        frozen = {
            "file_size_bytes": 12,
            "header_size_bytes": 3,
            "tensor_payload_size_bytes": 1,
            "tensor_count": 1,
            "parameter_count": 1,
            "dtypes": ["U8"],
            "metadata": {},
            "tensors": [{"name": "x"}],
        }
        compare_header(frozen, frozen)
        observed = copy.deepcopy(frozen)
        observed["tensors"][0]["name"] = "y"
        with self.assertRaises(ScoobyCheckpointAcquisitionError):
            compare_header(observed, frozen)

    def test_local_range_reader_returns_exact_requested_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.safetensors"
            path.write_bytes(b"0123456789")
            reader = LocalRangeReader(path)
            self.assertEqual(reader.read(2, 4), b"2345")
            self.assertEqual(reader.bytes_downloaded, 4)


if __name__ == "__main__":
    unittest.main()
