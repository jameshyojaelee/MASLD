from __future__ import annotations

import unittest

import numpy as np
import torch
from torch import nn

from scripts.mpralegnet_predict_gse281364 import aggregate_allele_scores, score_sequences


class _SumModel(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.scale = nn.Parameter(torch.tensor(1.0))

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return values.sum(dim=(1, 2)) * self.scale


class MPRALegNetGSE281364PredictionTests(unittest.TestCase):
    def test_cpu_scoring_uses_model_device_and_is_deterministic(self) -> None:
        names = ["a", "b", "c", "d"]
        records = {name: "A" * 230 for name in names}
        scores, repeat_max_abs = score_sequences(
            _SumModel().eval(), records, names, batch_size=2
        )
        np.testing.assert_array_equal(scores, np.full(4, 230.0, dtype=np.float32))
        self.assertEqual(repeat_max_abs, 0.0)

    def test_orientation_average_preserves_alt_minus_ref(self) -> None:
        reference, alternative, delta = aggregate_allele_scores(
            np.asarray([1.0, 3.0, 5.0, 9.0], dtype=np.float32)
        )
        self.assertEqual(reference, 3.0)
        self.assertEqual(alternative, 6.0)
        self.assertEqual(delta, 3.0)


if __name__ == "__main__":
    unittest.main()
