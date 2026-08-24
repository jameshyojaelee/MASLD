from __future__ import annotations

import unittest

from masld_bench.leakage import (
    LeakageError,
    SplitWindow,
    Variant,
    assert_group_integrity,
    assert_preprocessor_fit_scope,
    assert_sequence_window_isolation,
    merged_variant_blocks,
    salted_signature,
)
from masld_bench.topology import (
    MaskedModality,
    TopologyError,
    assert_biological_replication,
    assert_pairing,
    donor_pseudobulk,
)


class LeakageTopologyTests(unittest.TestCase):
    def test_grouped_split_and_training_only_preprocessing(self) -> None:
        with self.assertRaises(LeakageError):
            assert_group_integrity(
                [
                    {"cohort": "a", "donor": "d1", "outer_fold": 1},
                    {"cohort": "a", "donor": "d1", "outer_fold": 2},
                ],
                group_keys=["cohort", "donor"],
            )
        with self.assertRaises(LeakageError):
            assert_preprocessor_fit_scope(
                {"fit_partition_ids": ["train-1", "sealed-1"]},
                allowed_training_ids={"train-1"},
                sealed_ids={"sealed-1"},
            )

    def test_salted_fingerprint_and_merged_variant_blocks(self) -> None:
        signature = salted_signature("private-donor-vector", salt=b"0123456789abcdef")
        self.assertEqual(len(signature), 64)
        variants = [
            Variant("v1", "chr1", 100),
            Variant("v2", "chr1", 900_000),
            Variant("v3", "chr2", 100),
        ]
        blocks = merged_variant_blocks(
            variants,
            ld_edges_by_ancestry={"EUR": [("v2", "v3", 0.81)]},
        )
        self.assertEqual(len(set(blocks.values())), 1)

    def test_variant_block_inputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(LeakageError, "positive one-based"):
            merged_variant_blocks(
                [Variant("v1", "chr1", 0)], ld_edges_by_ancestry={}
            )
        with self.assertRaisesRegex(LeakageError, "between zero and one"):
            merged_variant_blocks(
                [Variant("v1", "chr1", 1)],
                ld_edges_by_ancestry={},
                ld_r2_threshold=1.1,
            )
        with self.assertRaisesRegex(LeakageError, "within \[0, 1\]"):
            merged_variant_blocks(
                [Variant("v1", "chr1", 1), Variant("v2", "chr2", 2)],
                ld_edges_by_ancestry={"EUR": [("v1", "v2", float("nan"))]},
            )
        with self.assertRaisesRegex(LeakageError, "unknown variant"):
            merged_variant_blocks(
                [Variant("v1", "chr1", 1)],
                ld_edges_by_ancestry={"EUR": [("v1", "v2", 0.9)]},
            )

    def test_missing_is_not_zero_and_false_pairing_is_rejected(self) -> None:
        with self.assertRaises(TopologyError):
            MaskedModality("structurally_missing", False, 0).validate()
        MaskedModality("structurally_missing", False, None).validate()
        with self.assertRaises(TopologyError):
            assert_pairing("GSE244832", "same_cell")
        for technical_unit in ("cell", "guide", "array", "section"):
            with self.subTest(technical_unit=technical_unit), self.assertRaises(
                TopologyError
            ):
                assert_biological_replication(
                    [{"unit_of_replication": technical_unit}]
                )
        assert_biological_replication(
            [{"unit_of_replication": "independent_experimental_batch"}]
        )
        counts = donor_pseudobulk(
            [
                {"donor": "d1", "gene": "g1", "count": 2},
                {"donor": "d1", "gene": "g1", "count": 3},
            ],
            donor_key="donor",
            feature_key="gene",
            count_key="count",
        )
        self.assertEqual(counts["d1"]["g1"], 5.0)
        reversed_counts = donor_pseudobulk(
            reversed(
                [
                    {"donor": "d1", "gene": "g1", "count": 2},
                    {"donor": "d1", "gene": "g1", "count": 3},
                ]
            ),
            donor_key="donor",
            feature_key="gene",
            count_key="count",
        )
        self.assertEqual(reversed_counts, counts)
        for invalid in (True, float("nan"), float("inf"), -1):
            with self.subTest(invalid=invalid), self.assertRaises(TopologyError):
                donor_pseudobulk(
                    [{"donor": "d1", "gene": "g1", "count": invalid}],
                    donor_key="donor",
                    feature_key="gene",
                    count_key="count",
                )

    def test_sequence_windows_respect_largest_receptive_field_buffer(self) -> None:
        isolated = [
            SplitWindow("train", "chr1", 0, 100, "train"),
            SplitWindow("test", "chr1", 1100, 1200, "withheld_sealed"),
        ]
        assert_sequence_window_isolation(isolated, boundary_buffer_bp=1000)
        with self.assertRaisesRegex(LeakageError, "receptive-field buffer"):
            assert_sequence_window_isolation(
                [
                    SplitWindow("train", "chr1", 0, 100, "train"),
                    SplitWindow("test", "chr1", 1099, 1200, "withheld_sealed"),
                ],
                boundary_buffer_bp=1000,
            )


if __name__ == "__main__":
    unittest.main()
