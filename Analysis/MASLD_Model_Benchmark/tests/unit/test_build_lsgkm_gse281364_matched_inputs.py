from __future__ import annotations

from pathlib import Path
import gzip
import tempfile
import unittest

import numpy as np

from scripts.build_lsgkm_gse281364_matched_inputs import (
    MatchedInputError,
    Window,
    assign_decile,
    decile_boundaries,
    deduplicate_windows,
    exclude_cross_class_sequence_duplicates,
    interval_overlap_bases,
    load_repeat_index,
    match_row,
    match_windows,
    matching_gate_result,
    mean_with_missing_zero,
    merged_intervals,
    outer_fold_join_matches,
    tss_bin,
    validate_config_boundary,
    write_fasta_gzip,
)


def window(
    identifier: str,
    *,
    start: int,
    sequence: str,
    gc: float,
    mappability: float,
    repeat: float,
    decile: int = 4,
    distance_bin: int = 2,
) -> Window:
    return Window(
        identifier=identifier,
        contig="chr1",
        start=start,
        end=start + 300,
        sequence=sequence,
        gc=gc,
        mappability=mappability,
        log1p_atac=2.0,
        tss_distance=50_000,
        tss_bin=distance_bin,
        repeat_fraction=repeat,
        atac_decile=decile,
    )


class LSGKMMatchedInputTests(unittest.TestCase):
    def test_matching_gate_reports_each_prespecified_failure_exactly(self) -> None:
        result = matching_gate_result(
            matched_pairs=7_499,
            eligible_positive_windows=75_000,
            minimum_matched_pairs=10_000,
            minimum_positive_coverage=0.1,
        )
        self.assertFalse(result["pair_count_gate_passed"])
        self.assertFalse(result["positive_coverage_gate_passed"])
        self.assertFalse(result["passed"])
        self.assertAlmostEqual(result["positive_coverage"], 7_499 / 75_000)

    def test_row_universe_fold_label_joins_numeric_borzoi_fold(self) -> None:
        self.assertTrue(outer_fold_join_matches("fold-1", "1"))
        self.assertFalse(outer_fold_join_matches("1", "1"))
        self.assertFalse(outer_fold_join_matches("fold-5", "5"))
        self.assertFalse(outer_fold_join_matches("fold-1", "fold-1"))

    def test_interval_union_and_repeat_fraction_geometry(self) -> None:
        index = merged_intervals(
            [
                ("chr1", 100, 180),
                ("chr1", 160, 220),
                ("chr1", 260, 280),
                ("chrM", 0, 100),
            ]
        )
        self.assertEqual(interval_overlap_bases(index, "chr1", 150, 300), 90)
        self.assertEqual(interval_overlap_bases(index, "chr1", 220, 260), 0)
        self.assertEqual(interval_overlap_bases(index, "chrM", 0, 100), 0)

    def test_tss_bins_are_frozen_half_open_distance_classes(self) -> None:
        self.assertEqual([tss_bin(value) for value in (0, 999, 1000, 9999, 10_000, 99_999, 100_000, 999_999, 1_000_000)], [0, 0, 1, 1, 2, 2, 3, 3, 4])

    def test_mappability_mean_treats_uncovered_bases_as_zero(self) -> None:
        import pyBigWig

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "partial.bw"
            writer = pyBigWig.open(str(path), "w")
            writer.addHeader([("chr1", 300)])
            writer.addEntries(["chr1"], [0], ends=[150], values=[1.0])
            writer.close()
            reader = pyBigWig.open(str(path))
            try:
                total = reader.stats("chr1", 0, 300, type="sum", exact=True)[0]
            finally:
                reader.close()
            self.assertEqual(mean_with_missing_zero(total, 300), 0.5)
            self.assertEqual(mean_with_missing_zero(None, 300), 0.0)
            with self.assertRaises(MatchedInputError):
                mean_with_missing_zero(301.0, 300)

    def test_repeatmasker_one_based_inclusive_coordinates_are_converted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "hg38.fa.out.gz"
            with gzip.open(source, "wt", encoding="ascii") as handle:
                handle.write("SW perc perc perc query begin end left strand repeat class begin end left id\n")
                handle.write("100 0.0 0.0 0.0 chr1 11 20 (100) + Alu SINE/Alu 1 10 (0) 1\n")
                handle.write("100 0.0 0.0 0.0 chr1 21 25 (95) + Alu SINE/Alu 11 15 (0)\n")
            index = load_repeat_index(source)
            self.assertEqual(interval_overlap_bases(index, "chr1", 10, 25), 15)
            self.assertEqual(interval_overlap_bases(index, "chr1", 9, 10), 0)

    def test_common_training_deciles_are_deterministic(self) -> None:
        boundaries = decile_boundaries([float(value) for value in range(100)])
        self.assertEqual(len(boundaries), 9)
        self.assertEqual(assign_decile(-1.0, boundaries), 0)
        self.assertEqual(assign_decile(100.0, boundaries), 9)

    def test_exact_and_reverse_complement_duplicates_are_grouped(self) -> None:
        sequence = "A" * 150 + "C" * 150
        reverse = "G" * 150 + "T" * 150
        retained = deduplicate_windows(
            [
                window("forward", start=100, sequence=sequence, gc=0.5, mappability=0.9, repeat=0.1),
                window("reverse", start=500, sequence=reverse, gc=0.5, mappability=0.9, repeat=0.1),
            ]
        )
        self.assertEqual([value.identifier for value in retained], ["forward"])

    def test_identical_sequence_cannot_receive_both_class_labels(self) -> None:
        sequence = "A" * 150 + "C" * 150
        positive = window("positive", start=100, sequence=sequence, gc=0.5, mappability=0.9, repeat=0.1)
        duplicate = window("negative_duplicate", start=500, sequence="G" * 150 + "T" * 150, gc=0.5, mappability=0.9, repeat=0.1)
        retained = window("negative_retained", start=900, sequence="AC" * 150, gc=0.5, mappability=0.9, repeat=0.1)
        self.assertEqual(
            [value.identifier for value in exclude_cross_class_sequence_duplicates([positive], [duplicate, retained])],
            ["negative_retained"],
        )

    def test_matcher_is_row_order_invariant_without_replacement(self) -> None:
        positives = [
            window("p1", start=100, sequence="A" * 300, gc=0.40, mappability=0.90, repeat=0.10),
            window("p2", start=500, sequence="C" * 300, gc=0.42, mappability=0.91, repeat=0.11),
        ]
        negatives = [
            window("n1", start=1_000, sequence="AC" * 150, gc=0.405, mappability=0.90, repeat=0.10),
            window("n2", start=1_500, sequence="AG" * 150, gc=0.415, mappability=0.91, repeat=0.11),
            window("wrong_decile", start=2_000, sequence="AT" * 150, gc=0.41, mappability=0.90, repeat=0.10, decile=3),
            window("wrong_map", start=2_500, sequence="CG" * 150, gc=0.41, mappability=0.80, repeat=0.10),
        ]
        first = match_windows(positives, negatives, 1103)
        second = match_windows(list(reversed(positives)), list(reversed(negatives)), 1103)
        identity = lambda rows: [(p.identifier, n.identifier, round(d, 12)) for p, n, d in rows]
        self.assertEqual(identity(first), identity(second))
        self.assertEqual(len(first), 2)
        self.assertEqual(len({negative.identifier for _positive, negative, _distance in first}), 2)
        self.assertNotIn("wrong_decile", {negative.identifier for _positive, negative, _distance in first})
        self.assertNotIn("wrong_map", {negative.identifier for _positive, negative, _distance in first})

    def test_match_manifest_exposes_both_exact_strata(self) -> None:
        positive = window("p", start=100, sequence="A" * 300, gc=0.40, mappability=0.90, repeat=0.10)
        negative = window("n", start=500, sequence="C" * 300, gc=0.42, mappability=0.91, repeat=0.11)
        row = match_row(1103, "donor0_genomic0", 0, positive, negative, 1.0)
        self.assertEqual(row["positive_contig"], row["negative_contig"])
        self.assertEqual(row["positive_atac_decile"], row["negative_atac_decile"])
        self.assertEqual(row["positive_tss_distance_bin"], row["negative_tss_distance_bin"])

    def test_deterministic_gzip_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.fa.gz"
            second = Path(directory) / "second.fa.gz"
            records = [("one", "ACGT" * 75), ("two", "TGCA" * 75)]
            self.assertEqual(write_fasta_gzip(first, records), 2)
            self.assertEqual(write_fasta_gzip(second, records), 2)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_action_and_minimal_design_firewall(self) -> None:
        config = {
            "status": "prespecified_outcome_blind_input_materialization",
            "outcome_access_authorized": False,
            "sealed_asset_access_authorized": False,
            "production_training_authorized": False,
            "production_prediction_authorized": False,
            "design": {
                "active_state_roster": ["hepatocyte"],
                "active_assay_contexts": ["HepG2_control", "HepG2_PAOA"],
                "context_specific_prediction": False,
                "fixed_seeds": [1103, 2909, 4721, 6673, 8111],
                "superseded_legacy_seeds": [11, 29, 47, 71, 101],
                "diagonal_split_count": 5,
                "shared_fit_count": 25,
                "readouts_per_shared_fit": 2,
                "independent_model_families": 1,
            },
        }
        validate_config_boundary(config)
        config["production_training_authorized"] = True
        with self.assertRaises(MatchedInputError):
            validate_config_boundary(config)


if __name__ == "__main__":
    unittest.main()
