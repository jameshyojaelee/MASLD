from __future__ import annotations

import unittest

from scripts.build_gse256398_terminal_disposition import (
    barcode_label_search,
    binomial_tail,
    donor_level_power,
    species_partition,
)


RELEVANT = {
    "healthy_control": 6,
    "masld_f0": 3,
    "mash_fibrosis": 4,
    "mash_cirrhosis": 4,
}


def sample(organism: str, *supplementary: str) -> dict[str, list[str]]:
    return {"organism": [organism], "supplementary": list(supplementary)}


class BarcodeLabelSearchTests(unittest.TestCase):
    def test_count_matrices_alone_are_not_cell_annotation(self) -> None:
        samples = {
            "GSM1": sample("Homo sapiens", "ftp://x/GSM1_CB_raw_feature_bc_matrix_filtered.h5"),
            "GSM2": sample("Homo sapiens", "ftp://x/GSM2_CB_raw_feature_bc_matrix_filtered.h5"),
        }
        found = barcode_label_search(samples)
        self.assertEqual(found["count_matrices_h5"], 2)
        self.assertEqual(found["candidate_annotation_files"], [])
        self.assertFalse(found["barcode_level_cell_annotation_released"])

    def test_an_annotation_file_would_be_detected(self) -> None:
        samples = {
            "GSM1": sample(
                "Homo sapiens",
                "ftp://x/GSM1_matrix.h5",
                "ftp://x/GSM1_cell_annotations.csv",
            )
        }
        found = barcode_label_search(samples)
        self.assertEqual(found["candidate_annotation_files"], ["GSM1_cell_annotations.csv"])
        self.assertTrue(found["barcode_level_cell_annotation_released"])

    def test_none_placeholder_is_not_counted_as_a_file(self) -> None:
        samples = {"GSM1": sample("Homo sapiens", "NONE")}
        self.assertEqual(barcode_label_search(samples)["sample_supplementary_files"], 0)


class SpeciesPartitionTests(unittest.TestCase):
    def test_mixed_species_series_is_flagged_and_mouse_listed(self) -> None:
        samples = {
            "GSM_h": sample("Homo sapiens"),
            "GSM_m": sample("Mus musculus"),
        }
        partition = species_partition(samples)
        self.assertTrue(partition["series_is_mixed_species"])
        self.assertEqual(partition["human_samples"], 1)
        self.assertEqual(partition["mouse_gsm"], ["GSM_m"])
        self.assertTrue(partition["human_only_filter_required"])


class DonorLevelPowerTests(unittest.TestCase):
    def test_binomial_tail_endpoints(self) -> None:
        self.assertAlmostEqual(binomial_tail(0, 5, 0.5), 1.0)
        self.assertAlmostEqual(binomial_tail(5, 5, 0.5), 0.5**5)

    def test_smallest_class_is_not_estimable_in_every_donor_safe_fold(self) -> None:
        power = donor_level_power(dict(RELEVANT))
        cv = power["donor_safe_grouped_cv"]
        self.assertEqual(cv["smallest_class_donors"], 3)
        self.assertFalse(cv["smallest_class_estimable_in_every_fold"])
        self.assertEqual(cv["outer_folds_guaranteed_without_smallest_class"], 2)

    def test_binary_contrast_is_underpowered_against_its_own_baseline(self) -> None:
        power = donor_level_power(dict(RELEVANT))
        binary = power["leave_one_donor_out_bounds"][
            "binary_healthy_vs_masld_spectrum"
        ]
        self.assertEqual(binary["correct_donors_required_at_alpha_0.05"], 15)
        self.assertLess(binary["power_at_true_accuracy_0.80"], 0.80)

    def test_permutation_floor_is_reported_for_the_binary_contrast(self) -> None:
        power = donor_level_power(dict(RELEVANT))
        binary = power["binary_healthy_vs_masld_spectrum"]
        self.assertEqual(binary["distinct_label_assignments"], 12376)
        self.assertLess(binary["minimum_attainable_two_sided_permutation_p"], 0.05)

    def test_group_sizes_are_reported_without_collapsing_ood_into_controls(self) -> None:
        power = donor_level_power(dict(RELEVANT))
        self.assertEqual(power["masld_relevant_donors"], 17)
        self.assertNotIn("alcohol_associated_cirrhosis", power["group_sizes"])
        self.assertNotIn("alcohol_associated_hepatitis", power["group_sizes"])


if __name__ == "__main__":
    unittest.main()
