from __future__ import annotations

import unittest

import torch

from scripts.probe_corgi_regular_film_head_l40s import (
    CorgiFilmProbeError,
    poisson_multinomial_loss,
    tensor_digest,
)


class CorgiFilmHeadL40SProbeTests(unittest.TestCase):
    def test_poisson_multinomial_loss_is_finite_and_differentiable(self) -> None:
        prediction = torch.full((1, 1, 16), 0.75, requires_grad=True)
        target = torch.linspace(0.25, 1.25, 16).reshape(1, 1, 16)
        loss = poisson_multinomial_loss(prediction, target)
        self.assertTrue(torch.isfinite(loss))
        loss.backward()
        self.assertIsNotNone(prediction.grad)
        self.assertTrue(torch.isfinite(prediction.grad).all())

    def test_loss_rejects_shape_mismatch(self) -> None:
        with self.assertRaises(CorgiFilmProbeError):
            poisson_multinomial_loss(
                torch.ones((1, 1, 16)), torch.ones((1, 1, 15))
            )

    def test_tensor_digest_binds_shape_and_values(self) -> None:
        left = torch.arange(12, dtype=torch.float32).reshape(3, 4)
        right = left.reshape(2, 6)
        changed = left.clone()
        changed[0, 0] = 1.0
        self.assertNotEqual(tensor_digest(left), tensor_digest(right))
        self.assertNotEqual(tensor_digest(left), tensor_digest(changed))


if __name__ == "__main__":
    unittest.main()
