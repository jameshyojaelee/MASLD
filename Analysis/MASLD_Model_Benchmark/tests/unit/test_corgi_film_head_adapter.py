from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import sys
import unittest

import torch
from torch import nn

ADAPTER_PATH = (
    Path(__file__).resolve().parents[2]
    / "src/masld_bench/adapters/corgi_film_head.py"
)
SPEC = importlib.util.spec_from_file_location("corgi_film_head_adapter", ADAPTER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot construct Corgi FiLM-head adapter source spec")
ADAPTER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ADAPTER
SPEC.loader.exec_module(ADAPTER)
CorgiFilmAtacModel = ADAPTER.CorgiFilmAtacModel
CorgiFilmHeadError = ADAPTER.CorgiFilmHeadError
adaptation_delta = ADAPTER.adaptation_delta
build_optimizer = ADAPTER.build_optimizer
load_adaptation_delta = ADAPTER.load_adaptation_delta
trainable_parameter_manifest = ADAPTER.trainable_parameter_manifest


class FakeReleasedCorgi(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.config = {
            "input_trans_regulators": 3,
            "output_channels": 22,
            "dim": 2,
            "output_central_bins": 4,
            "film_dimensions_conv": [1, 1],
            "film_dimensions_transformer": [1, 1],
        }
        self.film_mlp_conv = nn.Sequential(nn.BatchNorm1d(3), nn.Linear(3, 4))
        self.film_mlp_transformer = nn.Sequential(
            nn.BatchNorm1d(3), nn.Linear(3, 4)
        )
        self.film_layers = nn.ModuleList([nn.Identity() for _ in range(4)])
        self.max_pool = nn.Identity()
        self.conv_0 = nn.Identity()
        self.conv_1 = nn.Identity()
        self.conv_2 = nn.Identity()
        self.conv_3 = nn.Identity()
        self.conv_4 = nn.Identity()
        self.conv_5 = nn.Identity()
        self.transformer_layers = nn.ModuleList([nn.Identity() for _ in range(9)])
        self.crop = nn.Identity()
        self.final_conv = nn.Identity()
        self.output_head = nn.Conv1d(1_920, 22, kernel_size=1)
        with torch.no_grad():
            self.output_head.weight.copy_(
                torch.arange(22 * 1_920, dtype=torch.float32).reshape(22, 1_920, 1)
            )
            self.output_head.bias.copy_(torch.arange(22, dtype=torch.float32))


class CorgiFilmHeadAdapterTests(unittest.TestCase):
    def test_atac_row_is_copied_and_only_contract_prefixes_train(self) -> None:
        released = FakeReleasedCorgi()
        expected_weight = released.output_head.weight[1].detach().clone()
        expected_bias = released.output_head.bias[1].detach().clone()
        adapted = CorgiFilmAtacModel(released, require_exact=False)
        self.assertTrue(torch.equal(adapted.masld_atac_head.weight[0], expected_weight))
        self.assertTrue(torch.equal(adapted.masld_atac_head.bias[0], expected_bias))
        manifest = trainable_parameter_manifest(adapted)
        self.assertGreater(manifest["trainable_parameter_numel"], 0)
        self.assertTrue(manifest["batchnorm_running_statistics_frozen_by_eval_mode"])
        trainable = [name for name, value in adapted.named_parameters() if value.requires_grad]
        self.assertTrue(all(name.startswith(tuple(manifest["trainable_prefixes"])) for name in trainable))
        self.assertFalse(any(name.startswith("conv_0.") for name in trainable))

    def test_optimizer_has_exactly_film_and_head_groups(self) -> None:
        adapted = CorgiFilmAtacModel(FakeReleasedCorgi(), require_exact=False)
        optimizer = build_optimizer(adapted)
        self.assertEqual(
            [group["group_name"] for group in optimizer.param_groups],
            ["film_mlp", "atac_head"],
        )
        optimized = {id(parameter) for group in optimizer.param_groups for parameter in group["params"]}
        expected = {id(parameter) for parameter in adapted.parameters() if parameter.requires_grad}
        self.assertEqual(optimized, expected)

    def test_delta_round_trip_restores_only_trainable_state(self) -> None:
        adapted = CorgiFilmAtacModel(FakeReleasedCorgi(), require_exact=False)
        frozen_before = {
            name: value.detach().clone()
            for name, value in adapted.named_parameters()
            if not value.requires_grad
        }
        delta = adaptation_delta(adapted)
        original = copy.deepcopy(delta["trainable_state_dict"])
        with torch.no_grad():
            for parameter in adapted.parameters():
                if parameter.requires_grad:
                    parameter.add_(1.0)
        load_adaptation_delta(adapted, delta)
        state = dict(adapted.named_parameters())
        for name, expected in original.items():
            self.assertTrue(torch.equal(state[name], expected))
        for name, expected in frozen_before.items():
            self.assertTrue(torch.equal(state[name], expected))

    def test_delta_rejects_changed_source_identity(self) -> None:
        adapted = CorgiFilmAtacModel(FakeReleasedCorgi(), require_exact=False)
        delta = adaptation_delta(adapted)
        delta["source_checkpoint_sha256"] = "0" * 64
        with self.assertRaises(CorgiFilmHeadError):
            load_adaptation_delta(adapted, delta)

    def test_rejects_sequential_or_softplus_released_head(self) -> None:
        released = FakeReleasedCorgi()
        released.output_head = nn.Sequential(nn.Conv1d(1_920, 22, 1), nn.Softplus())
        with self.assertRaises(CorgiFilmHeadError):
            CorgiFilmAtacModel(released, require_exact=False)


if __name__ == "__main__":
    unittest.main()
