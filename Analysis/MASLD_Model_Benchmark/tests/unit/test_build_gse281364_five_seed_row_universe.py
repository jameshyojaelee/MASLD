from __future__ import annotations

from collections import Counter
from copy import deepcopy
from pathlib import Path
import tempfile
import tomllib
import unittest

from scripts.build_gse281364_five_seed_row_universe import (
    CONTEXTS,
    ROW_FIELDS,
    SEEDS,
    STATUS,
    RowUniverseError,
    build_rows,
    hash_id,
    preflight,
    validate_authorities,
    validate_config,
    validate_prediction_identities,
    validate_registry_states,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/gse281364_five_seed_row_universe.toml"


class GSE281364FiveSeedRowUniverseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = tomllib.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        validate_config(cls.config)
        cls.authorities = validate_authorities(ROOT, cls.config)
        cls.rows, cls.metadata = build_rows(ROOT, cls.authorities, cls.config)

    def test_exact_common_row_namespace_and_five_seed_roster(self) -> None:
        self.assertEqual(tuple(self.config["row_contract"]["row_universe_fields"]), ROW_FIELDS)
        self.assertEqual(len(self.rows), 10_330)
        self.assertEqual({row["seed"] for row in self.rows}, set(SEEDS))
        self.assertEqual({row["assay_context_id"] for row in self.rows}, set(CONTEXTS))
        per_seed = Counter(row["seed"] for row in self.rows)
        self.assertEqual(per_seed, Counter({seed: 2066 for seed in SEEDS}))
        row_sets = {
            seed: {row["row_hash"] for row in self.rows if row["seed"] == seed}
            for seed in SEEDS
        }
        self.assertTrue(all(values == row_sets[SEEDS[0]] for values in row_sets.values()))
        first = self.rows[0]
        self.assertEqual(
            first["row_hash"],
            hash_id(
                "gse281364",
                "mpra_allelic_direction",
                str(first["element_id"]),
                str(first["source_locus_group_id"]),
                str(first["assay_context_id"]),
            ),
        )

    def test_long_range_fold_map_replaces_the_unsafe_original_map(self) -> None:
        self.assertEqual(self.metadata["long_range_blocks"], 239)
        self.assertEqual(
            self.metadata["fold_source_locus_groups"],
            {"fold-0": 208, "fold-1": 208, "fold-2": 207, "fold-3": 208, "fold-4": 202},
        )
        self.assertEqual(
            self.metadata["source_groups_reassigned_from_original_6kb_fold_map"],
            830,
        )
        group_folds = {}
        block_folds = {}
        for row in self.rows:
            group_folds.setdefault(row["source_locus_group_id"], row["outer_fold"])
            self.assertEqual(group_folds[row["source_locus_group_id"]], row["outer_fold"])
            block_folds.setdefault(row["long_range_block_id"], row["outer_fold"])
            self.assertEqual(block_folds[row["long_range_block_id"]], row["outer_fold"])

    def test_prediction_and_registry_sources_are_not_prematurely_eligible(self) -> None:
        identities = validate_prediction_identities(ROOT, self.authorities)
        registry = validate_registry_states(ROOT, self.authorities)
        self.assertEqual(identities["common_head_seed_ids"], [20260824])
        self.assertEqual(identities["mpralegnet_construct_rows"], 4359)
        self.assertFalse(identities["borzoi_prediction_available"])
        self.assertFalse(identities["context_prediction_authority_bound"])
        self.assertIn("candidate|target_label_unexposed", registry["corgi_regular"])
        self.assertTrue(
            all(
                not source["current_open_trigger_eligible"]
                for source in self.config["source_status"]
            )
        )

    def test_fold_denominator_or_conditional_source_injection_is_rejected(self) -> None:
        changed = deepcopy(self.config)
        changed["fold_contract"]["original_6kb_fold_predictions_directly_stackable"] = True
        with self.assertRaises(RowUniverseError):
            validate_config(changed)

        changed = deepcopy(self.config)
        changed["outcome_contract"]["all_4359_constructs_are_rows"] = True
        with self.assertRaises(RowUniverseError):
            validate_config(changed)

        changed = deepcopy(self.config)
        blocked = next(
            source
            for source in changed["source_status"]
            if source["source_id"] == "blocked_context_candidates"
        )
        blocked["eligible_after_gates"] = True
        with self.assertRaises(RowUniverseError):
            validate_config(changed)

    def test_integrated_build_is_outcome_blind_and_immutable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "row_contract"
            receipt = preflight(root=ROOT, config_path=CONFIG_PATH, output=output)
            self.assertEqual(receipt["status"], STATUS)
            self.assertEqual(receipt["base_rows_per_seed"], 2066)
            self.assertEqual(receipt["seeded_rows"], 10_330)
            self.assertEqual(receipt["current_trigger_compatible_pair_rows"], 0)
            self.assertFalse(receipt["outcomes_read"])
            self.assertFalse(receipt["prediction_values_read"])
            self.assertFalse(receipt["reporter_counts_read"])
            self.assertFalse(receipt["sealed_assets_read"])
            self.assertFalse(receipt["predictions_generated"])
            self.assertFalse(receipt["stack_fit"])
            self.assertFalse(receipt["residual_correlation_calculated"])
            self.assertFalse(receipt["conditional_model_built_or_fit"])
            self.assertTrue((output / "row_universe.tsv").is_file())
            self.assertTrue((output / "source_status.tsv").is_file())
            with self.assertRaises(RowUniverseError):
                preflight(root=ROOT, config_path=CONFIG_PATH, output=output)


if __name__ == "__main__":
    unittest.main()
