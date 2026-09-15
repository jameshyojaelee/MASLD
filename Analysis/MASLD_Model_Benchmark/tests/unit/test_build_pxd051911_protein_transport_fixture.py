"""Contract tests for the PXD051911 protein_transport fixture builder."""

from __future__ import annotations

import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ACTIVATION = ROOT / "executions" / "pxd051911-activation-readiness-21109461"
SOURCE = ACTIVATION / "source"
JOIN = ACTIVATION / "audit" / "readiness" / "participant_join.tsv"

_spec = importlib.util.spec_from_file_location(
    "build_pxd051911_protein_transport_fixture",
    ROOT / "scripts" / "build_pxd051911_protein_transport_fixture.py",
)
builder = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(builder)


class FoldIndexMirrorsRepository(unittest.TestCase):
    def test_fold_index_matches_the_repository_implementation(self) -> None:
        from masld_bench.evaluators.rna_atac_development import fold_index as canonical

        for unit in ("pxd051911:ID1", "pxd051911:ID93", "pxd051911:ID55"):
            self.assertEqual(
                builder.fold_index(unit, seed=20260821, outer_folds=5),
                canonical(unit, seed=20260821, outer_folds=5),
            )

    def test_state_codes_map_into_the_missing_state_vocabulary(self) -> None:
        from masld_bench.topology import MISSING_STATES

        for value in builder.STATE_TO_MISSING_STATE.values():
            self.assertIn(value, MISSING_STATES)
        self.assertEqual(set(builder.STATE_CODES), set(builder.STATE_TO_MISSING_STATE))
        self.assertEqual(len(set(builder.STATE_CODES.values())), 4)


class FixtureBuild(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = Path(tempfile.mkdtemp(prefix="protein-transport-fixture-"))
        cls.output = cls._tmp / "fixture"
        cls.receipt = builder.build(SOURCE, JOIN, cls.output)

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls._tmp, ignore_errors=True)

    def test_participant_count_and_donor_safety(self) -> None:
        self.assertEqual(self.receipt["participants"], 58)
        self.assertEqual(self.receipt["unit_of_replication"], "participant")
        self.assertFalse(self.receipt["plasma_arm_used"])
        self.assertFalse(self.receipt["rna_protein_pair_inferred"])

    def test_feature_axis_is_the_union_and_panels_are_source_native(self) -> None:
        axis = self.receipt["feature_axis"]
        self.assertEqual(axis["liver_panel"], 7096)
        self.assertEqual(axis["scwat_panel"], 6988)
        self.assertEqual(axis["owat_panel"], 6988)
        self.assertEqual(axis["union_proteins"], 8361)

    def test_cross_tissue_scopes_match_the_activation_artifact(self) -> None:
        self.assertEqual(self.receipt["cross_tissue_scopes"]["S-SCWAT"], 56)
        self.assertEqual(self.receipt["cross_tissue_scopes"]["S-OWAT"], 25)
        self.assertEqual(
            self.receipt["cross_tissue_scopes"]["pairing_level"], "same_donor_different_tissue"
        )
        self.assertTrue(self.receipt["cross_tissue_scopes"]["never_same_sample"])

    def test_fold_sizes_are_frozen_and_partition_the_cohort(self) -> None:
        sizes = self.receipt["split"]["fold_sizes"]
        self.assertEqual(sizes, [10, 17, 10, 14, 7])
        self.assertEqual(sum(sizes), 58)
        self.assertFalse(self.receipt["split"]["stratified"])

    def test_every_participant_sits_in_exactly_one_outer_fold(self) -> None:
        import csv

        with (self.output / "folds" / "participant_outer_folds.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 58)
        self.assertEqual(len({r["person_key"] for r in rows}), 58)
        by_key: dict[str, set[str]] = {}
        for row in rows:
            by_key.setdefault(row["person_key"], set()).add(row["outer_fold"])
        self.assertTrue(all(len(folds) == 1 for folds in by_key.values()))

    def test_no_identity_group_crosses_an_outer_fold(self) -> None:
        import csv

        from masld_bench.leakage import assert_group_integrity

        with (self.output / "folds" / "participant_outer_folds.tsv").open() as handle:
            rows = [dict(r) for r in csv.DictReader(handle, delimiter="\t")]
        assert_group_integrity(
            rows, group_keys=["cohort_family_id", "person_key"], fold_key="outer_fold"
        )

    def test_missingness_is_typed_and_never_zero(self) -> None:
        self.assertFalse(self.receipt["missingness"]["missing_encoded_as_zero"])
        self.assertTrue(self.receipt["missingness"]["explicit_typed_state_matrix"])
        for tissue in builder.TISSUES:
            values = np.load(self.output / "molecular" / f"{tissue}_values.npy")
            states = np.load(self.output / "molecular" / f"{tissue}_state.npy")
            observed = np.load(self.output / "molecular" / f"{tissue}_observed_mask.npy")
            self.assertEqual(values.dtype, np.float64)
            self.assertEqual(states.dtype, np.uint8)
            self.assertEqual(observed.dtype, np.bool_)
            self.assertEqual(values.shape, (58, 8361))
            self.assertEqual(states.shape, values.shape)
            np.testing.assert_array_equal(observed, states == builder.STATE_OBSERVED)
            self.assertFalse(bool(np.any(values[observed] == 0.0)))
            self.assertFalse(bool(np.any(np.isnan(values[observed]))))
            self.assertTrue(bool(np.all(np.isnan(values[~observed]))))

    def test_structural_and_detection_missingness_are_kept_apart(self) -> None:
        states = np.load(self.output / "molecular" / "owat_state.npy")
        present = {int(code) for code in np.unique(states)}
        self.assertIn(builder.STATE_TISSUE_NOT_COLLECTED, present)
        self.assertIn(builder.STATE_BELOW_QC, present)
        self.assertIn(builder.STATE_PROTEIN_NOT_IN_PANEL, present)
        collected = np.load(self.output / "molecular" / "owat_participant_collected_mask.npy")
        self.assertEqual(int(collected.sum()), 25)
        # A participant with no omental tissue carries one uniform structural state.
        for row in np.flatnonzero(~collected):
            self.assertEqual(
                set(np.unique(states[row])), {builder.STATE_TISSUE_NOT_COLLECTED}
            )

    def test_fat_pct_is_masked_not_zero_filled(self) -> None:
        self.assertEqual(self.receipt["missingness"]["fat_pct_observed"], 44)
        self.assertEqual(self.receipt["missingness"]["fat_pct_structurally_missing"], 14)

    def test_nas_equals_its_component_sum(self) -> None:
        self.assertTrue(self.receipt["endpoint_integrity"]["nas_equals_component_sum"])
        self.assertEqual(self.receipt["endpoint_integrity"]["checked_participants"], 58)

    def test_values_are_source_native(self) -> None:
        values = self.receipt["values"]
        self.assertTrue(values["source_native_raw_intensity"])
        self.assertFalse(values["normalisation_applied"])
        self.assertFalse(values["log_transform_applied"])
        self.assertFalse(values["imputation_applied"])

    def test_builder_refuses_to_overwrite_an_existing_fixture(self) -> None:
        with self.assertRaises(FileExistsError):
            builder.build(SOURCE, JOIN, self.output)


if __name__ == "__main__":
    unittest.main()
