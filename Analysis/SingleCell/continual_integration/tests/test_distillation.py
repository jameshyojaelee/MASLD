from __future__ import annotations

import unittest

import torch
from torch.distributions import Normal
from scvi import REGISTRY_KEYS
from scvi.module._constants import MODULE_KEYS

from masld_cl.distillation import (
    LatentDistillationError, LatentDistillationRegularizer,
    posterior_means_at_anchor,
)


class TinyEncoder(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.encoder = torch.nn.Linear(2, 2, bias=False)
        with torch.no_grad():
            self.encoder.weight.copy_(torch.tensor([[1.0, 0.5], [-0.25, 2.0]]))

    def _get_inference_input(self, tensors):
        return {"x": tensors[REGISTRY_KEYS.X_KEY]}

    def inference(self, x):
        mean = self.encoder(x)
        return {MODULE_KEYS.QZ_KEY: Normal(mean, torch.ones_like(mean))}


def batch(x, indices):
    return {
        REGISTRY_KEYS.X_KEY: torch.as_tensor(x, dtype=torch.float32),
        REGISTRY_KEYS.INDICES_KEY: torch.as_tensor(indices, dtype=torch.int64),
    }


class TestLatentDistillation(unittest.TestCase):
    def setUp(self):
        self.module = TinyEncoder()
        self.x = torch.tensor([[1.0, 2.0], [3.0, 1.0], [2.0, 4.0]])
        with torch.no_grad():
            self.targets = self.module.encoder(self.x)
        self.regularizer = LatentDistillationRegularizer(
            self.targets, torch.tensor([True, False, True])
        )

    def test_zero_at_anchor_and_exact_permuted_index_matching(self):
        value = self.regularizer.penalty(
            self.module, batch(self.x[[2, 0, 1]], [2, 0, 1])
        )
        self.assertEqual(float(value), 0.0)
        maximum = self.regularizer.assert_zero_at_anchor(
            self.module,
            [batch(self.x[:2], [0, 1]), batch(self.x[2:], [2])],
        )
        self.assertEqual(maximum, 0.0)

    def test_pretraining_anchor_reader_uses_qz_mean_and_restores_mode(self):
        self.module.train()
        observed = posterior_means_at_anchor(
            self.module,
            [batch(self.x[:2], [0, 1]), batch(self.x[2:], [2])],
        )
        torch.testing.assert_close(observed, self.targets)
        self.assertTrue(self.module.training)

    def test_positive_encoder_perturbation_and_explicit_gradient_fixture(self):
        with torch.no_grad():
            self.module.encoder.weight[0, 0].add_(0.2)
        observed = self.regularizer.penalty(
            self.module, batch(self.x[[2, 0]], [2, 0])
        )
        current = self.module.encoder(self.x[[2, 0]])
        expected = ((current - self.targets[[2, 0]]) ** 2).sum(dim=1).mean()
        torch.testing.assert_close(observed, expected, rtol=1e-5, atol=1e-7)
        observed.backward()
        observed_gradient = self.module.encoder.weight.grad.detach().clone()

        comparison = TinyEncoder()
        with torch.no_grad():
            comparison.encoder.weight[0, 0].add_(0.2)
        manual = ((comparison.encoder(self.x[[2, 0]]) - self.targets[[2, 0]]) ** 2).sum(1).mean()
        manual.backward()
        torch.testing.assert_close(
            observed_gradient, comparison.encoder.weight.grad, rtol=1e-5, atol=1e-7
        )
        self.assertGreater(float(observed), 0.0)

    def test_query_only_batch_has_zero_penalty_and_zero_gradient(self):
        value = self.regularizer.penalty(self.module, batch(self.x[1:2], [1]))
        value.backward()
        self.assertEqual(float(value), 0.0)
        torch.testing.assert_close(
            self.module.encoder.weight.grad,
            torch.zeros_like(self.module.encoder.weight), rtol=0, atol=0,
        )

    def test_semisupervised_tuple_and_training_mode_restored(self):
        self.module.train()
        value = self.regularizer.penalty(
            self.module, (batch(self.x[[0, 2]], [0, 2]), batch(self.x[:1], [0]))
        )
        self.assertEqual(float(value), 0.0)
        self.assertTrue(self.module.training)

    def test_missing_duplicate_and_out_of_range_indices_fail_closed(self):
        with self.assertRaisesRegex(LatentDistillationError, "lacks exact"):
            self.regularizer.penalty(
                self.module, {REGISTRY_KEYS.X_KEY: self.x[:1]}
            )
        with self.assertRaisesRegex(LatentDistillationError, "duplicate"):
            self.regularizer.penalty(
                self.module, batch(self.x[:2], [0, 0])
            )
        with self.assertRaisesRegex(LatentDistillationError, "outside"):
            self.regularizer.penalty(
                self.module, batch(self.x[:1], [3])
            )

    def test_invalid_target_contracts_fail_closed(self):
        with self.assertRaises(LatentDistillationError):
            LatentDistillationRegularizer(
                torch.ones(2, 2), torch.zeros(2, dtype=torch.bool)
            )
        with self.assertRaises(LatentDistillationError):
            LatentDistillationRegularizer(
                torch.tensor([[float("nan")]]), torch.tensor([True])
            )

    @unittest.skipUnless(torch.cuda.is_available(), "CUDA-specific device invariant")
    def test_cuda_anchor_visit_audit_keeps_cpu_bookkeeping(self):
        module = TinyEncoder().cuda()
        targets = module.encoder(self.x.cuda()).detach().cpu()
        regularizer = LatentDistillationRegularizer(
            targets, torch.tensor([True, False, True])
        ).cuda()
        maximum = regularizer.assert_zero_at_anchor(
            module,
            [batch(self.x[:2], [0, 1]), batch(self.x[2:], [2])],
        )
        self.assertEqual(maximum, 0.0)


if __name__ == "__main__":
    unittest.main()
