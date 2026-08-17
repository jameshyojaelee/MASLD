from __future__ import annotations

import copy
import unittest

import torch
from torch import nn

from masld_cl.ewc import (
    EWCError,
    EWCRegularizer,
    build_anchor_and_masks,
    compute_empirical_fisher,
)


class ExpandedModel(nn.Module):
    def __init__(self, width=2):
        super().__init__()
        self.weight = nn.Parameter(torch.tensor([[1.0, -2.0] + [0.5] * (width - 2)]))
        self.bias = nn.Parameter(torch.tensor([0.25]))

    def forward(self, x):
        return x @ self.weight.t() + self.bias


class TestEWC(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(17)

    def test_anchor_requires_reference_loaded_first(self):
        reference = ExpandedModel(2)
        state = {name: value.detach().clone() for name, value in reference.named_parameters()}
        query = ExpandedModel(3)
        with torch.no_grad():
            query.weight[:, :2].zero_()
        with self.assertRaisesRegex(EWCError, "not the loaded reference"):
            build_anchor_and_masks(query, state)

    def test_new_category_column_is_masked(self):
        reference = ExpandedModel(2)
        state = {name: value.detach().clone() for name, value in reference.named_parameters()}
        query = ExpandedModel(3)
        anchors, masks = build_anchor_and_masks(query, state)
        self.assertTrue(torch.equal(masks["weight"], torch.tensor([[1.0, 1.0, 0.0]])))
        ones = {name: torch.ones_like(value) for name, value in anchors.items()}
        regularizer = EWCRegularizer(anchors, masks, ones, ones)
        regularizer.assert_zero_at_anchor(query)
        with torch.no_grad():
            query.weight[0, 2] += 10
        self.assertEqual(float(regularizer.penalty(query)), 0.0)
        with torch.no_grad():
            query.weight[0, 0] += 2
        self.assertGreater(float(regularizer.penalty(query)), 0.0)

    def test_empirical_fisher_is_finite_nonnegative(self):
        model = ExpandedModel(2)
        batches = [torch.randn(4, 2), torch.randn(4, 2), torch.randn(2, 2)]

        def loss_fn(current, batch):
            return current(batch).square().mean()

        fisher, summary = compute_empirical_fisher(
            model, batches, loss_fn, expected_batch_size=4
        )
        self.assertEqual(summary.n_batches, 2)
        self.assertEqual(summary.skipped_partial_batches, 1)
        for value in fisher.values():
            self.assertTrue(torch.isfinite(value).all())
            self.assertTrue(torch.all(value >= 0))

    def test_empirical_fisher_can_consume_exact_partial_sample(self):
        model = ExpandedModel(2)
        batches = [torch.ones(4, 2), torch.full((2, 2), 2.0)]

        fisher, summary = compute_empirical_fisher(
            model, batches, lambda current, batch: current(batch).square().mean(),
            expected_batch_size=4, include_partial=True,
        )
        self.assertEqual(summary.n_batches, 2)
        self.assertEqual(summary.n_observations, 6)
        self.assertEqual(summary.skipped_partial_batches, 0)
        self.assertEqual(
            summary.weighting,
            "observations_over_squared_minibatch_mean_gradients",
        )
        self.assertTrue(all(torch.isfinite(value).all() for value in fisher.values()))
        torch.testing.assert_close(
            fisher["weight"], torch.full((1, 2), 107.0 / 6.0)
        )
        torch.testing.assert_close(
            fisher["bias"], torch.tensor([33.5 / 6.0])
        )

    def test_lambda_zero_matches_unregularized_step(self):
        left = ExpandedModel(2)
        right = copy.deepcopy(left)
        state = {name: value.detach().clone() for name, value in left.named_parameters()}
        anchors, masks = build_anchor_and_masks(left, state)
        ones = {name: torch.ones_like(value) for name, value in anchors.items()}
        regularizer = EWCRegularizer(anchors, masks, ones, ones)
        x = torch.tensor([[1.0, 2.0], [-1.0, 0.5]])
        for model, use_regularizer in ((left, False), (right, True)):
            optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
            optimizer.zero_grad()
            loss = model(x).square().mean()
            if use_regularizer:
                loss = loss + 0.0 * regularizer.penalty(model)
            loss.backward()
            optimizer.step()
        for a, b in zip(left.parameters(), right.parameters()):
            self.assertTrue(torch.equal(a, b))

    def test_explicit_loss_and_gradient_fixture(self):
        model = nn.Linear(2, 1, bias=False)
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[1.0, 2.0]]))
        anchor = {"weight": model.weight.detach().clone()}
        mask = {"weight": torch.tensor([[1.0, 0.0]])}
        reference = {"weight": torch.tensor([[2.0, 3.0]])}
        control = {"weight": torch.tensor([[5.0, 7.0]])}
        regularizer = EWCRegularizer(anchor, mask, reference, control)
        with torch.no_grad():
            model.weight[0, 0] += 0.4
        penalty = regularizer.penalty(model)
        objective = 0.5 * penalty
        objective.backward()
        torch.testing.assert_close(
            penalty, torch.tensor(1.6), rtol=1e-5, atol=1e-7
        )
        torch.testing.assert_close(
            model.weight.grad, torch.tensor([[4.0, 0.0]]), rtol=1e-5, atol=1e-7
        )

    def test_block_mean_product_normalization_is_scale_invariant(self):
        model = ExpandedModel(2)
        state = {name: value.detach().clone() for name, value in model.named_parameters()}
        anchors, masks = build_anchor_and_masks(model, state)
        reference = {name: torch.full_like(value, 2.0) for name, value in anchors.items()}
        control = {name: torch.full_like(value, 5.0) for name, value in anchors.items()}
        regularizer = EWCRegularizer(
            anchors, masks, reference, control, normalization="block_mean_product"
        )
        scaled = EWCRegularizer(
            anchors, masks,
            {name: value * 1e6 for name, value in reference.items()},
            {name: value * 1e-3 for name, value in control.items()},
            normalization="block_mean_product",
        )
        with torch.no_grad():
            model.weight.add_(0.5)
            model.bias.add_(0.25)
        torch.testing.assert_close(regularizer.penalty(model), scaled.penalty(model))


if __name__ == "__main__":
    unittest.main()
