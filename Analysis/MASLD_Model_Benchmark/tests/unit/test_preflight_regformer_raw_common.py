#!/usr/bin/env python3
"""Unit checks for the RegFormer raw-common strict-restore gate."""

from __future__ import annotations

from pathlib import Path
import unittest

from scripts import preflight_regformer_raw_common as preflight


ROOT = Path(__file__).parents[2]
ACQUISITION = ROOT / "executions/regformer-acquisition-21066043"
INVENTORY = (
    ROOT
    / "executions/regformer-safe-inventory-21066069/checkpoint_inventory.json"
)
INSPECTION = (
    ROOT
    / "executions/regformer-weights-only-inspection-21066098/weights_only_inspection.json"
)
WRAPPER = ROOT / "slurm/preflight_regformer_raw_common_cpu.sbatch"


class RegFormerRawCommonPreflightTests(unittest.TestCase):
    def test_exact_small_asset_evidence_and_state_schema(self) -> None:
        self.assertEqual(
            preflight.sha256_file(ACQUISITION / "vocab.json"),
            preflight.VOCAB_SHA256,
        )
        self.assertEqual(
            preflight.sha256_file(ACQUISITION / "output_graph_no_cycles.dgl"),
            preflight.GRAPH_SHA256,
        )
        self.assertEqual(
            preflight.sha256_file(ACQUISITION / "pretrain.toml"),
            preflight.CONFIG_SHA256,
        )
        self.assertEqual(preflight.sha256_file(INVENTORY), preflight.INVENTORY_SHA256)
        schema, audit = preflight.inspection_schema(INSPECTION)
        self.assertEqual(len(schema), 139)
        self.assertEqual(schema["encoder.embedding.weight"]["shape"], [60_698, 512])
        self.assertEqual(schema["lm_head.weight"]["shape"], [60_698, 512])
        self.assertEqual(schema["encoder.layer_emb.pe"]["shape"], [100, 1, 512])
        self.assertEqual(audit["state_tensor_count"], 139)

    def test_released_vocab_append_and_config_are_exact(self) -> None:
        vocab = preflight.validate_vocab(ACQUISITION / "vocab.json")
        self.assertEqual(vocab["released_size"], 60_697)
        self.assertEqual(vocab["runtime_size_after_released_mask_append"], 60_698)
        self.assertEqual(vocab["mask_id"], 60_697)
        self.assertEqual(
            preflight.validate_config(ACQUISITION / "pretrain.toml"),
            dict(sorted(preflight.REQUIRED_CONFIG_LITERALS.items())),
        )

    def test_released_label_refinement_is_detected_and_prohibited(self) -> None:
        source = """
def run(model, batch, out, values, mask):
    raw = model._get_cell_emb_from_layer(out, values, src_key_padding_mask=mask)
    labels = batch["celltype_labels"]
    return refine_embedding(raw, batch_ids, labels, 2)
"""
        audit = preflight.audit_released_embedding_source(source)
        self.assertLess(audit["raw_pool_line"], audit["forbidden_refine_embedding_line"])
        self.assertEqual(
            audit["benchmark_action"],
            "do_not_import_or_execute_released_embedding_driver",
        )
        with self.assertRaises(preflight.RegFormerPreflightError):
            preflight.audit_released_embedding_source(
                source.replace("refine_embedding", "identity")
            )

    def test_assignment_parser_ignores_comments_but_not_values(self) -> None:
        observed = preflight.parse_literal_assignments(
            "layer_size = 512 # fixed\ncell_emb_style='avg-pool'\n"
        )
        self.assertEqual(
            observed, {"layer_size": "512", "cell_emb_style": "'avg-pool'"}
        )

    def test_cli_accepts_no_labels_or_outcomes(self) -> None:
        options = {
            option
            for action in preflight.parser()._actions
            for option in action.option_strings
        }
        self.assertFalse(
            any(
                marker in option
                for option in options
                for marker in ("label", "cell-type", "outcome", "histology")
            )
        )

    def test_wrapper_is_cpu_nslab_generic_and_non_submitting(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(
            line for line in text.splitlines() if line.startswith("#SBATCH")
        )
        self.assertIn("--job-name=model-check-031", header)
        self.assertIn("--partition=cpu", header)
        self.assertIn("--account=nslab", header)
        self.assertIn("--qos=nslab", header)
        self.assertNotIn("--gres", header)
        self.assertNotIn("--array", header)
        self.assertNotIn("sbatch ", text)
        self.assertIn(preflight.SOURCE_COMMIT, text)
        self.assertIn(preflight.RUNTIME_ARTIFACTS_SHA256, text)
        self.assertIn("strict_restore_receipt.json", text)


if __name__ == "__main__":
    unittest.main()
