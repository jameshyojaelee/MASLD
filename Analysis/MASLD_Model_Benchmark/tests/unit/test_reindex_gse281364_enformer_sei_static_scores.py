#!/usr/bin/env python3
"""Unit tests for long-range-safe Enformer and Sei static-score binding."""

from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.reindex_gse281364_enformer_sei_static_scores import (
    ROW_FIELDS,
    StaticReindexError,
    collapse_row_universe,
    indexed_source_rows,
    validate_sei_scalarization,
    write_gzip_tsv,
)


def row(seed: int, context: str, *, fold: str = "fold-0") -> dict[str, str]:
    return {
        "seed": str(seed),
        "row_hash": f"row-{seed}-{context}",
        "unit_hash": "unit-a",
        "block_hash": "block-a",
        "stratum": "all",
        "outer_fold": fold,
        "study_id": "gse281364",
        "assay_context_id": context,
        "element_id": "chr1:10-20",
        "source_locus_group_id": "outer-a",
        "long_range_block_id": "long-a",
    }


class StaticReindexTests(unittest.TestCase):
    def test_row_universe_collapses_seed_context_replicates(self) -> None:
        rows = [
            row(seed, context)
            for seed in (11, 13)
            for context in ("HepG2_control", "HepG2_PAOA")
        ]
        collapsed = collapse_row_universe(
            rows,
            expected_seeds=(11, 13),
            expected_contexts=("HepG2_control", "HepG2_PAOA"),
        )
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(collapsed[0]["long_range_block_id"], "long-a")

    def test_row_universe_rejects_context_fold_drift(self) -> None:
        rows = [
            row(seed, context)
            for seed in (11, 13)
            for context in ("HepG2_control", "HepG2_PAOA")
        ]
        rows[-1]["outer_fold"] = "fold-1"
        with self.assertRaises(StaticReindexError):
            collapse_row_universe(
                rows,
                expected_seeds=(11, 13),
                expected_contexts=("HepG2_control", "HepG2_PAOA"),
            )

    def test_source_index_rejects_duplicate_element(self) -> None:
        fields = ("element_id", "outer_fold")
        rows = [
            {"element_id": "a", "outer_fold": "0"},
            {"element_id": "a", "outer_fold": "0"},
        ]
        with self.assertRaises(StaticReindexError):
            indexed_source_rows(rows, fields, label="toy")

    def test_sei_scalarization_rejects_misaligned_summary(self) -> None:
        matrix = np.asarray([[1.0, -2.0, 0.5], [3.0, 2.0, -1.0]], dtype=np.float32)
        rows = [
            {
                "max_abs_sequence_class_index0": "1",
                "signed_max_abs_score": "-2",
                "max_abs_score": "2",
                "mean_signed_score": str(float(matrix[0].mean())),
            },
            {
                "max_abs_sequence_class_index0": "0",
                "signed_max_abs_score": "3",
                "max_abs_score": "3",
                "mean_signed_score": "999",
            },
        ]
        with self.assertRaises(StaticReindexError):
            validate_sei_scalarization(matrix, rows)

    def test_deterministic_gzip_has_zero_timestamp(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.tsv.gz"
            second = Path(directory) / "second.tsv.gz"
            values = [{"a": "1", "b": "2"}]
            write_gzip_tsv(first, ("a", "b"), values)
            write_gzip_tsv(second, ("a", "b"), values)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with gzip.open(first, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(list(reader), values)
            self.assertEqual(first.read_bytes()[4:8], b"\x00\x00\x00\x00")


if __name__ == "__main__":
    unittest.main()
