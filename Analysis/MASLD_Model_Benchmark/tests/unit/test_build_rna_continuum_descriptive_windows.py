from __future__ import annotations

import math
import unittest

from scripts.build_rna_continuum_descriptive_windows import (
    ARM_PRIMARY,
    ARM_QUARANTINED,
    AXIS_LEAKAGE,
    HISTOLOGY_ANCHOR_STATUS,
    PRIMARY_COHORTS,
    QUARANTINED_COHORTS,
    WINDOW_CENTERS,
    WINDOW_WIDTH,
    FeatureScores,
    Substrate,
    build_windows,
    design_effect_se_multiplier,
    expected_row_count,
    feature_shape,
    in_window,
    load_hotspot_scores,
    load_system_scores,
    reproduce_frozen_windows,
    summarize_window,
    window_bounds,
)


class WindowGeometryTests(unittest.TestCase):
    def test_nine_fixed_windows_of_width_one_fifth(self) -> None:
        self.assertEqual(len(WINDOW_CENTERS), 9)
        self.assertAlmostEqual(WINDOW_WIDTH, 0.2)

    def test_intervals_are_lower_closed_and_upper_open(self) -> None:
        lower, upper, final = window_bounds(1)
        self.assertAlmostEqual(lower, 0.0)
        self.assertAlmostEqual(upper, 0.2)
        self.assertFalse(final)
        self.assertTrue(in_window(0.0, 1))
        self.assertFalse(in_window(0.2, 1))

    def test_final_window_closes_its_upper_bound(self) -> None:
        _, upper, final = window_bounds(9)
        self.assertTrue(final)
        self.assertAlmostEqual(upper, 1.0)
        self.assertTrue(in_window(1.0, 9))

    def test_windows_overlap_so_a_sample_can_occupy_two(self) -> None:
        occupied = [i for i in range(1, 10) if in_window(0.25, i)]
        self.assertEqual(occupied, [2, 3])


class SummaryTests(unittest.TestCase):
    def test_single_sample_window_has_a_mean_but_no_standard_error(self) -> None:
        stat = summarize_window([1.5])
        self.assertEqual(stat.n_samples, 1)
        self.assertAlmostEqual(stat.mean, 1.5)
        self.assertIsNone(stat.se)

    def test_empty_window_is_reported_as_zero_not_as_a_zero_score(self) -> None:
        stat = summarize_window([])
        self.assertEqual(stat.n_samples, 0)
        self.assertIsNone(stat.mean)
        self.assertIsNone(stat.se)

    def test_standard_error_uses_the_sample_standard_deviation(self) -> None:
        stat = summarize_window([1.0, 2.0, 3.0])
        self.assertAlmostEqual(stat.mean, 2.0)
        self.assertAlmostEqual(stat.se, 1.0 / math.sqrt(3))

    def test_design_effect_multiplier_matches_the_gse193066_ratio(self) -> None:
        self.assertAlmostEqual(design_effect_se_multiplier(), math.sqrt(164 / 106))
        self.assertGreater(design_effect_se_multiplier(), 1.0)


class ArmSeparationTests(unittest.TestCase):
    def test_the_two_arms_are_disjoint_and_cover_all_five_cohorts(self) -> None:
        self.assertEqual(set(PRIMARY_COHORTS) & set(QUARANTINED_COHORTS), set())
        self.assertEqual(len(set(PRIMARY_COHORTS) | set(QUARANTINED_COHORTS)), 5)

    def test_quarantined_arm_declares_the_leakage_on_every_row(self) -> None:
        self.assertEqual(AXIS_LEAKAGE[ARM_QUARANTINED],
                         "source_overlap_with_signature_discovery")
        self.assertNotEqual(AXIS_LEAKAGE[ARM_PRIMARY], AXIS_LEAKAGE[ARM_QUARANTINED])

    def test_gse126848_anchor_is_never_recorded_not_merely_missing(self) -> None:
        # "missing" invites someone to go looking for data that was never collected.
        self.assertEqual(HISTOLOGY_ANCHOR_STATUS["GSE126848"], "never_recorded_in_source")
        self.assertIn("nas_never_recorded", HISTOLOGY_ANCHOR_STATUS["GSE213621"])


def _substrate() -> Substrate:
    axis = {
        ("GSE162694", "S1"): 0.05,
        ("GSE162694", "S2"): 0.15,
        ("GSE162694", "S3"): 0.95,
    }
    return Substrate(
        axis=axis,
        stage_complete={("GSE162694", "S1")},
        scale_native={"GSE162694": "kleiner_0_4_plus_normal_histology"},
        stage_comparable={"GSE162694": True},
        sex_source={"GSE162694": "annotated"},
        cohort_n={"GSE162694": 3},
    )


def _layer(testable_scores: bool = True) -> FeatureScores:
    scores = {
        "prog_a": {
            ("GSE162694", "S1"): 1.0,
            ("GSE162694", "S2"): 3.0,
            ("GSE162694", "S3"): 5.0,
        }
    }
    if not testable_scores:
        scores["prog_b"] = {}
    return FeatureScores(layer="hotspot_program", aggregation="within_cohort_z",
                         display={"prog_a": "Program A", "prog_b": "Program B"},
                         scores=scores)


class BuildWindowTests(unittest.TestCase):
    def test_every_row_names_the_sample_unit_and_bars_inference(self) -> None:
        rows = build_windows(_substrate(), [_layer()], ("GSE162694",), ARM_PRIMARY, {})
        self.assertEqual(len(rows), 9)
        for row in rows:
            self.assertEqual(row["unit_of_observation"], "sequencing_sample")
            self.assertIsNone(row["n_donors"])
            self.assertEqual(row["n_donors_status"], "untestable_no_donor_key_on_disk")
            self.assertFalse(row["inference_use"])
            self.assertTrue(row["visualization_only"])
            self.assertTrue(row["longitudinal_interpretation_barred"])

    def test_stage_complete_and_missing_are_reported_per_window(self) -> None:
        rows = build_windows(_substrate(), [_layer()], ("GSE162694",), ARM_PRIMARY, {})
        first = next(row for row in rows if row["window_id"] == 1)
        # S1 (0.05) and S2 (0.15) fall in window 1; only S1 has a recorded stage.
        self.assertEqual(first["n_samples_in_window"], 2)
        self.assertEqual(first["n_samples_stage_complete"], 1)
        self.assertEqual(first["n_samples_stage_missing"], 1)

    def test_stage_counts_always_sum_to_the_window_occupancy(self) -> None:
        rows = build_windows(_substrate(), [_layer()], ("GSE162694",), ARM_PRIMARY, {})
        for row in rows:
            self.assertEqual(
                row["n_samples_stage_complete"] + row["n_samples_stage_missing"],
                row["n_samples_in_window"],
            )

    def test_design_effect_column_is_a_sensitivity_beside_the_naive_value(self) -> None:
        rows = build_windows(_substrate(), [_layer()], ("GSE162694",), ARM_PRIMARY, {})
        row = next(r for r in rows if r["se_score_naive_sample_level"] is not None)
        self.assertGreater(row["se_score_design_effect_sensitivity"],
                           row["se_score_naive_sample_level"])
        self.assertIn("out_of_substrate", row["design_effect_source"])

    def test_untestable_program_zero_occupancy_names_the_gate_not_a_defect(self) -> None:
        gate = {("prog_b", "GSE162694"): {"testable": False,
                                          "retained_l1_fraction": 0.72004694346637,
                                          "minimum_retained_l1_fraction": 0.8}}
        rows = build_windows(_substrate(), [_layer(testable_scores=False)],
                             ("GSE162694",), ARM_PRIMARY, gate)
        empty = [row for row in rows if row["feature_uid"] == "prog_b"]
        self.assertEqual(len(empty), 9)
        for row in empty:
            self.assertEqual(row["n_samples_in_window"], 0)
            self.assertIn("not_testable_after_signature_exclusion",
                          row["zero_occupancy_reason"])
            self.assertFalse(row["feature_testable"])

    def test_quarantined_arm_stamps_leakage_on_every_row(self) -> None:
        substrate = _substrate()
        axis = dict(substrate.axis)
        axis[("GSE126848", "Q1")] = 0.5
        quarantined = Substrate(
            axis=axis, stage_complete=set(),
            scale_native={"GSE126848": "not_recorded"},
            stage_comparable={"GSE126848": False},
            sex_source={"GSE126848": "annotated"},
            cohort_n={"GSE126848": 1},
        )
        layer = FeatureScores(layer="hotspot_program", aggregation="within_cohort_z",
                              display={"prog_a": "Program A"},
                              scores={"prog_a": {("GSE126848", "Q1"): 2.0}})
        rows = build_windows(quarantined, [layer], ("GSE126848",), ARM_QUARANTINED, {})
        for row in rows:
            self.assertEqual(row["axis_leakage"],
                             "source_overlap_with_signature_discovery")
            self.assertEqual(row["histology_anchor_status"], "never_recorded_in_source")


class UnscoredFeatureTests(unittest.TestCase):
    """The 21114163 defect: an all-NA feature must appear at zero, not disappear.

    hotspot_fibroblasts_8584ba834d31e608 (Activated stellate, PDGFRA) carries 844 source
    rows whose outcome_z is NA in every one, because the signature-exclusion testability
    gate leaves it unscored. Registering features by parsed score dropped it from the
    layer, which is worse than an unexplained zero: a reader comparing 116 emitted
    against 117 in the source has no row to inspect.
    """

    def _rows(self) -> list[dict[str, str]]:
        rows = []
        for sample in ("S1", "S2"):
            rows.append({"program_uid": "prog_scored", "sample_id": sample,
                         "dataset": "GSE162694", "outcome_z": "1.0",
                         "module_name": "Scored"})
            rows.append({"program_uid": "prog_all_na", "sample_id": sample,
                         "dataset": "GSE162694", "outcome_z": "NA",
                         "module_name": "Activated stellate (PDGFRA)"})
        return rows

    def test_all_na_program_is_registered_not_dropped(self) -> None:
        layers = load_hotspot_scores(self._rows())
        self.assertEqual(set(layers[0].scores), {"prog_scored", "prog_all_na"})
        self.assertEqual(layers[0].scores["prog_all_na"], {})

    def test_all_na_program_keeps_its_display_name(self) -> None:
        layers = load_hotspot_scores(self._rows())
        self.assertEqual(layers[0].display["prog_all_na"],
                         "Activated stellate (PDGFRA)")

    def test_all_na_program_emits_a_full_set_of_zero_occupancy_windows(self) -> None:
        layers = load_hotspot_scores(self._rows())
        rows = build_windows(_substrate(), layers, ("GSE162694",), ARM_PRIMARY, {})
        emitted = [row for row in rows if row["feature_uid"] == "prog_all_na"]
        self.assertEqual(len(emitted), 9)
        self.assertTrue(all(row["n_samples_in_window"] == 0 for row in emitted))

    def test_all_na_system_is_registered_under_every_aggregation(self) -> None:
        rows = [
            {"community_id": "system_01", "sample_id": "S1", "dataset": "GSE162694",
             "aggregation_method": "collection_balanced", "system_score_z": "1.0",
             "system_display": "system_01"},
            {"community_id": "system_02", "sample_id": "S1", "dataset": "GSE162694",
             "aggregation_method": "collection_balanced", "system_score_z": "NA",
             "system_display": "system_02"},
            {"community_id": "system_01", "sample_id": "S1", "dataset": "GSE162694",
             "aggregation_method": "feature_balanced_median", "system_score_z": "2.0",
             "system_display": "system_01"},
            {"community_id": "system_02", "sample_id": "S1", "dataset": "GSE162694",
             "aggregation_method": "feature_balanced_median", "system_score_z": "NA",
             "system_display": "system_02"},
        ]
        layers = load_system_scores(rows)
        self.assertEqual(len(layers), 2)
        for layer in layers:
            self.assertEqual(set(layer.scores), {"system_01", "system_02"})


class ShapeAssertionTests(unittest.TestCase):
    def test_expected_row_count_is_features_times_windows_times_cohorts(self) -> None:
        layer = FeatureScores(layer="hotspot_program", aggregation="within_cohort_z",
                              display={}, scores={f"p{i}": {} for i in range(117)})
        self.assertEqual(expected_row_count([layer], 2), 117 * 9 * 2)
        self.assertEqual(expected_row_count([layer], 3), 117 * 9 * 3)

    def test_expected_row_count_sums_across_layers_and_aggregations(self) -> None:
        hotspot = FeatureScores(layer="hotspot_program", aggregation="within_cohort_z",
                                display={}, scores={f"p{i}": {} for i in range(117)})
        systems = [
            FeatureScores(layer="molecular_system", aggregation=name, display={},
                          scores={f"s{i}": {} for i in range(43)})
            for name in ("collection_balanced", "feature_balanced_median")
        ]
        self.assertEqual(expected_row_count([hotspot, *systems], 2),
                         117 * 9 * 2 + 43 * 9 * 2 * 2)

    def test_emitted_rows_match_the_expected_shape(self) -> None:
        layers = load_hotspot_scores([
            {"program_uid": "a", "sample_id": "S1", "dataset": "GSE162694",
             "outcome_z": "1.0", "module_name": "A"},
            {"program_uid": "b", "sample_id": "S1", "dataset": "GSE162694",
             "outcome_z": "NA", "module_name": "B"},
        ])
        rows = build_windows(_substrate(), layers, ("GSE162694",), ARM_PRIMARY, {})
        self.assertEqual(len(rows), expected_row_count(layers, 1))

    def test_shape_receipt_pairs_expected_with_emitted(self) -> None:
        layers = load_hotspot_scores([
            {"program_uid": "a", "sample_id": "S1", "dataset": "GSE162694",
             "outcome_z": "1.0", "module_name": "A"},
            {"program_uid": "b", "sample_id": "S1", "dataset": "GSE162694",
             "outcome_z": "NA", "module_name": "B"},
        ])
        shape = feature_shape(layers)
        self.assertEqual(shape[0]["features_expected"], 2)
        self.assertEqual(shape[0]["features_emitted"], 2)
        self.assertEqual(shape[0]["features_with_no_finite_score"], 1)


class ReproductionControlTests(unittest.TestCase):
    def _computed(self) -> list[dict[str, object]]:
        return build_windows(_substrate(), [_layer()], ("GSE162694",), ARM_PRIMARY, {})

    def _frozen_from(self, computed) -> list[dict[str, str]]:
        return [
            {
                "program_uid": row["feature_uid"],
                "dataset": row["dataset"],
                "window_id": str(row["window_id"]),
                "n_participants": str(row["n_samples_in_window"]),
                "mean_score": ("NA" if row["mean_score_within_cohort_z"] is None
                               else repr(row["mean_score_within_cohort_z"])),
                "se_score": ("NA" if row["se_score_naive_sample_level"] is None
                             else repr(row["se_score_naive_sample_level"])),
            }
            for row in computed
        ]

    def test_matching_frozen_table_reproduces(self) -> None:
        computed = self._computed()
        control = reproduce_frozen_windows(computed, self._frozen_from(computed))
        self.assertEqual(control["n_checked"], 9)
        self.assertEqual(control["n_mismatched"], 0)

    def test_a_changed_frozen_mean_is_detected(self) -> None:
        computed = self._computed()
        frozen = self._frozen_from(computed)
        frozen[0]["mean_score"] = "99.0"
        control = reproduce_frozen_windows(computed, frozen)
        self.assertEqual(control["n_mismatched"], 1)

    def test_a_changed_frozen_occupancy_is_detected(self) -> None:
        computed = self._computed()
        frozen = self._frozen_from(computed)
        frozen[0]["n_participants"] = "999"
        control = reproduce_frozen_windows(computed, frozen)
        self.assertEqual(control["n_mismatched"], 1)


if __name__ == "__main__":
    unittest.main()
