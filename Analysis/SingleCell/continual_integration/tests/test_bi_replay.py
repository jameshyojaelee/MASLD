from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import torch
from torch import nn

from masld_cl.bi_replay import compute_bi_scores, load_bi_replay_unlock
from masld_cl.config import canonical_json_bytes
from masld_cl.contracts import ContractError, sha256_path


class FakeEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.scale = nn.Parameter(torch.tensor(1.0))
        self.batch_shapes: list[tuple[int, int]] = []
        self.zero_counts: list[torch.Tensor] = []

    def _get_inference_input(self, tensors):
        from scvi import REGISTRY_KEYS

        return {"x": tensors[REGISTRY_KEYS.X_KEY]}

    def inference(self, x):
        from scvi.module._constants import MODULE_KEYS

        self.batch_shapes.append(tuple(x.shape))
        self.zero_counts.append((x == 0).sum(dim=1).detach().cpu())
        location = torch.stack(
            (x[:, 0], x[:, 1] + 0.25 * x[:, 2]), dim=1
        ) * self.scale
        return {
            MODULE_KEYS.Z_KEY: location + torch.randn_like(location),
            MODULE_KEYS.QZ_KEY: torch.distributions.Normal(
                location, torch.ones_like(location)
            )
        }


class TestBIReplay(unittest.TestCase):
    def test_bi_scores_use_exactly_200_half_masked_augmentations(self):
        from scvi import REGISTRY_KEYS

        module = FakeEncoder()
        batches = [
            {REGISTRY_KEYS.X_KEY: torch.ones(2, 10)},
            {REGISTRY_KEYS.X_KEY: torch.ones(1, 10)},
        ]
        scores = compute_bi_scores(module, batches, seed=17)
        self.assertEqual(scores.shape, (3,))
        self.assertTrue(np.isfinite(scores).all())
        self.assertTrue(np.all(scores >= 0))
        self.assertEqual(module.batch_shapes, [(400, 10), (200, 10)])
        self.assertTrue(all(torch.all(value == 5) for value in module.zero_counts))
        repeated = compute_bi_scores(FakeEncoder(), batches, seed=17)
        np.testing.assert_array_equal(scores, repeated)

    def test_bi_unlock_is_hash_bound_to_promotion_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifacts = []
            for name in ("promotion.json", "metrics.tsv", "metrics.lock.json"):
                path = root / name
                path.write_text(name)
                artifacts.append(path)
            config = {"_config_sha256": "config"}
            selection = {
                "lock_sha256": "selection",
                "selected": {"ewc_lambda": 1.0, "replay_fraction": 0.2},
            }
            unlock = {
                "schema_version": "masld-cl-bi-replay-unlock-v1",
                "config_sha256": "config",
                "selection_lock_sha256": "selection",
                "selected_setting": selection["selected"],
                "promotion_decision_realpath": str(artifacts[0]),
                "promotion_decision_sha256": sha256_path(artifacts[0]),
                "metrics_realpath": str(artifacts[1]),
                "metrics_sha256": sha256_path(artifacts[1]),
                "metrics_lock_realpath": str(artifacts[2]),
                "metrics_lock_sha256": sha256_path(artifacts[2]),
                "metric_bundle_lock_sha256": "bundle",
                "sensitivity_only": True,
            }
            unlock["lock_sha256"] = hashlib.sha256(
                canonical_json_bytes(unlock)
            ).hexdigest()
            path = root / "bi_unlock.json"
            path.write_text(json.dumps(unlock))
            loaded = load_bi_replay_unlock(path, config, selection)
            self.assertEqual(loaded["lock_sha256"], unlock["lock_sha256"])
            artifacts[0].write_text("changed")
            with self.assertRaises(ContractError):
                load_bi_replay_unlock(path, config, selection)


if __name__ == "__main__":
    unittest.main()
