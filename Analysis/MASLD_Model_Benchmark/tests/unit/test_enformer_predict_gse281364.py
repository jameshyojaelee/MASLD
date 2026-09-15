from __future__ import annotations

import tempfile
from pathlib import Path
import json
import unittest

import numpy as np

from scripts.enformer_predict_gse281364 import (
    TRACKS,
    EnformerPredictionError,
    _load_chunk,
    _prepare_cache,
    _write_chunk,
    reverse_complement,
    score_sums,
)


class EnformerPredictionTests(unittest.TestCase):
    def test_sad_sar_and_rc_ensemble(self) -> None:
        values = np.ones((4, TRACKS), dtype=np.float32)
        values[1] = 3.0
        values[2] = 5.0
        values[3] = 9.0
        result = score_sums(values)
        self.assertEqual(result.shape, (4, TRACKS))
        np.testing.assert_allclose(result[0], 2.0)
        np.testing.assert_allclose(result[2], 3.0)
        np.testing.assert_allclose(result[3], np.log2(7.0) - np.log2(4.0), rtol=1e-6)
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")

    def test_append_only_chunk_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            array = np.zeros((2, 4, TRACKS), dtype=np.float32)
            path = _write_chunk(root, 0, 2, array, deterministic_repeat_max_abs=0.0)
            np.testing.assert_array_equal(_load_chunk(path, 0, 2), array)
            with self.assertRaises(EnformerPredictionError):
                _write_chunk(root, 0, 2, array, deterministic_repeat_max_abs=0.0)

    def test_cached_chunk_checksum_and_bounds_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            array = np.zeros((2, 4, TRACKS), dtype=np.float32)
            path = _write_chunk(root, 0, 2, array, deterministic_repeat_max_abs=0.0)
            with (path / "features.npy").open("r+b") as handle:
                handle.seek(-1, 2)
                handle.write(b"x")
            with self.assertRaises(EnformerPredictionError):
                _load_chunk(path, 0, 2)
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            array = np.zeros((2, 4, TRACKS), dtype=np.float32)
            path = _write_chunk(root, 0, 2, array, deterministic_repeat_max_abs=0.0)
            receipt = json.loads((path / "receipt.json").read_text())
            receipt["end"] = 3
            (path / "receipt.json").write_text(json.dumps(receipt, sort_keys=True) + "\n")
            with self.assertRaises(EnformerPredictionError):
                _load_chunk(path, 0, 2)

    def test_cache_contract_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            contract = {"cache_id": "fixed", "manifest_sha256": "a", "outcomes_read": False}
            _prepare_cache(root, contract)
            changed = dict(contract, manifest_sha256="b")
            with self.assertRaises(EnformerPredictionError):
                _prepare_cache(root, changed)

    def test_negative_sum_fails(self) -> None:
        values = np.ones((4, TRACKS), dtype=np.float32)
        values[0, 0] = -1
        with self.assertRaises(EnformerPredictionError):
            score_sums(values)


if __name__ == "__main__":
    unittest.main()
