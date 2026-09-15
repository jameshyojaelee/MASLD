from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.reconcile_borzoi_grelu_to_local_port_mapping import (
    BorzoiMappingError,
    map_key,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]


class BorzoiSemanticMappingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(
            (ROOT / "config/borzoi_grelu_to_local_port_mapping_probe.json").read_text()
        )

    def test_contract_is_fail_closed(self) -> None:
        validate_contract(self.contract)
        self.assertFalse(
            self.contract["historical_source_candidate"][
                "checkpoint_producing_revision_proven"
            ]
        )
        self.assertFalse(self.contract["mapping_contract"]["model_forward_allowed"])

    def test_contract_rejects_state_dict_order_mapping(self) -> None:
        changed = deepcopy(self.contract)
        changed["mapping_contract"]["state_dict_order_used_for_mapping"] = True
        with self.assertRaisesRegex(BorzoiMappingError, "contract opened"):
            validate_contract(changed)

    def test_representative_semantic_mappings(self) -> None:
        cases = {
            "embedding.conv_tower.blocks.0.conv.weight": "conv_dna.conv_layer.weight",
            "embedding.conv_tower.blocks.6.norm.layer.running_mean": "unet1.1.norm.running_mean",
            "embedding.transformer_tower.blocks.3.mha.to_pos_k.weight": "transformer.3.0.fn.1.to_rel_k.weight",
            "embedding.transformer_tower.blocks.7.ffn.dense2.linear.bias": "transformer.7.1.fn.4.bias",
            "embedding.unet_tower.blocks.0.sconv.depthwise.weight": "separable1.conv_layer.0.weight",
            "embedding.unet_tower.blocks.1.channel_transform.conv.layer.bias": "horizontal_conv0.conv_layer.bias",
            "embedding.pointwise_conv.norm.layer.weight": "final_joined_convs.0.norm.weight",
            "head.channel_transform.conv.layer.weight": "human_head.weight",
        }
        for source, expected in cases.items():
            self.assertEqual(map_key(source)[0], expected)

    def test_unknown_key_fails_closed(self) -> None:
        with self.assertRaisesRegex(BorzoiMappingError, "unmapped"):
            map_key("embedding.unknown.weight")


if __name__ == "__main__":
    unittest.main()
