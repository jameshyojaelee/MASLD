from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import tomllib
import unittest

from scripts.prepare_gse281364_sequence_context_refit_dag import (
    COMPONENT_FIELDS,
    MODEL_IDS,
    RefitDagError,
    materialize_contexts,
    prepare,
    validate_authorities,
    validate_config,
    validate_context_mappers,
    validate_inventory,
    validate_registry_states,
    validate_stages,
)
from scripts.project_gse281364_long_range_embeddings import (
    LongRangeProjectionError,
    build_long_range_index,
    read_tsv,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/gse281364_sequence_context_refit_dag.toml"


class GSE281364SequenceContextRefitDagTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = tomllib.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        validate_config(cls.config)
        cls.authorities = validate_authorities(ROOT, cls.config)

    def test_exact_model_inventory_registry_and_common_schema(self) -> None:
        inventory = validate_inventory(ROOT, self.config)
        states = validate_registry_states(self.authorities)
        self.assertEqual({row["model_id"] for row in inventory}, MODEL_IDS)
        inventory_to_registry = (
            MODEL_IDS
            - {
                "enformer_crested_restricted_port",
                "zero_training_mean_allele_identity_ridge",
            }
        ) | {"enformer"}
        self.assertEqual(set(states), inventory_to_registry)
        self.assertEqual(
            tuple(self.config["common_prediction_schema"]["fields"]),
            COMPONENT_FIELDS,
        )
        open_models = {
            row["model_id"]
            for row in inventory
            if row["open_after_refit_and_external_gates"] and not row["diagnostic_only"]
        }
        self.assertEqual(
            open_models,
            {
                "caduceus",
                "dnabert2",
                "hyenadna",
                "zero_training_mean_allele_identity_ridge",
            },
        )

    def test_long_range_row_binding_replaces_every_embedded_old_fold(self) -> None:
        rows = read_tsv(self.authorities["row_tsv"])
        index = build_long_range_index(rows)
        self.assertEqual(len(index), 1033)
        self.assertEqual(len({row["long_range_block_id"] for row in index.values()}), 239)
        self.assertEqual(
            {row["outer_fold"] for row in index.values()},
            {"fold-0", "fold-1", "fold-2", "fold-3", "fold-4"},
        )
        changed = [dict(row) for row in rows]
        changed[0]["outer_fold"] = "fold-9"
        with self.assertRaises(LongRangeProjectionError):
            build_long_range_index(changed)

    def test_dag_is_topological_and_ready_nodes_have_real_entrypoints(self) -> None:
        stages = validate_stages(ROOT, self.config)
        ready = [stage for stage in stages if stage["command"] != "BLOCKED"]
        self.assertEqual(
            {stage["stage_id"] for stage in ready},
            {
                "authority_reconciliation",
                "long_range_projection_open_sequence",
                "long_range_projection_restricted_sequence",
                "materialize_corgi_context_roster",
            },
        )
        self.assertTrue(all(not stage["outcomes_required"] for stage in ready))
        self.assertFalse(any("residual" in stage["stage_id"] for stage in stages))

    def test_context_path_is_executable_but_not_condition_matched_or_open(self) -> None:
        mappers = validate_context_mappers(ROOT, self.config)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "contexts"
            receipt = materialize_contexts(self.config, mappers, output)
            self.assertEqual(receipt["contexts"], 5)
            self.assertEqual(receipt["context_width"], 2891)
            self.assertFalse(receipt["condition_matched"])
            self.assertFalse(receipt["outcomes_read"])
            self.assertFalse(receipt["model_forward_executed"])
            self.assertTrue((output / "contexts.npz").is_file())
            self.assertTrue((output / "context_roster.tsv").is_file())

    def test_firewall_mutations_are_rejected(self) -> None:
        changed = deepcopy(self.config)
        changed["outcome_access_authorized"] = True
        with self.assertRaises(RefitDagError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["row_binding"]["original_6kb_predictions_reusable"] = True
        with self.assertRaises(RefitDagError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["eligibility"]["cross_cohort_hepatocyte_context_is_condition_matched"] = True
        with self.assertRaises(RefitDagError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["eligibility"]["current_open_context_paths"] = ["corgi_regular"]
        with self.assertRaises(RefitDagError):
            validate_config(changed)

    def test_integrated_reconciliation_does_not_open_outcomes_or_fit_models(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "dag"
            receipt = prepare(ROOT, CONFIG_PATH, output)
            self.assertEqual(receipt["current_trigger_compatible_pair_count"], 0)
            self.assertFalse(receipt["outcome_member_opened"])
            self.assertFalse(receipt["outcomes_read"])
            self.assertFalse(receipt["prediction_values_read"])
            self.assertFalse(receipt["predictions_generated"])
            self.assertFalse(receipt["model_fit"])
            self.assertFalse(receipt["stack_fit"])
            self.assertFalse(receipt["residual_correlation_calculated"])
            self.assertFalse(receipt["conditional_model_built_or_fit"])
            dag = json.loads((output / "refit_dag.json").read_text(encoding="utf-8"))
            self.assertEqual(dag["current_trigger_compatible_pair_count"], 0)
            self.assertFalse(dag["context_condition_matched"])
            self.assertTrue((output / "artifact_inventory.tsv").is_file())
            self.assertTrue((output / "common_prediction_schema.tsv").is_file())
            with self.assertRaises(RefitDagError):
                prepare(ROOT, CONFIG_PATH, output)


if __name__ == "__main__":
    unittest.main()
