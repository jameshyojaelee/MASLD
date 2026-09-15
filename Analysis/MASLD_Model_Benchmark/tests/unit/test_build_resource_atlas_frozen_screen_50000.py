#!/usr/bin/env python3
"""Unit checks for the deterministic 50,000-cell Atlas frozen screen."""

from __future__ import annotations

from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts import build_resource_atlas_frozen_screen_50000 as builder


ROOT = Path(__file__).parents[2]
WRAPPER = ROOT / "slurm/build_resource_atlas_frozen_screen_50000.sbatch"


class FrozenScreen50000Tests(unittest.TestCase):
    @staticmethod
    def _write_manifest(path: Path, *, small: int = 1, large: int = 7) -> None:
        fields = sorted(builder.REQUIRED_MANIFEST_FIELDS)
        rows: list[dict[str, str]] = []
        ordinal = 0
        for broad, source_types in builder.ONTOLOGY.items():
            for donor, count in ((f"{broad}-small", small), (f"{broad}-large", large)):
                for _ in range(count):
                    source_cell_id = f"cell-{ordinal:04d}"
                    library_id = f"library-{ordinal % 3}"
                    rows.append(
                        {
                            "analysis_eligible": "True",
                            "cell_id": f"{library_id}|{source_cell_id}",
                            "cell_type": source_types[0],
                            "dataset": "synthetic-development",
                            "donor_id": donor,
                            "library_id": library_id,
                            "source_cell_id": source_cell_id,
                        }
                    )
                    ordinal += 1
        rows.append(
            {
                "analysis_eligible": "False",
                "cell_id": "excluded-library|excluded-cell",
                "cell_type": "not-used",
                "dataset": "synthetic-development",
                "donor_id": "excluded-donor",
                "library_id": "excluded-library",
                "source_cell_id": "excluded-cell",
            }
        )
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_frozen_budget_and_disclosure_contract(self) -> None:
        self.assertEqual(builder.CELL_BUDGET, 50_000)
        self.assertEqual(builder.CELLS_PER_CLASS, 10_000)
        self.assertEqual(len(builder.ONTOLOGY), 5)
        self.assertIn("source_cell_type", builder.OBS_FIELDS)
        self.assertIn("broad_label", builder.OBS_FIELDS)
        forbidden = {"fibrosis", "nas", "histology", "sealed_outcome", "age", "sex"}
        self.assertTrue(forbidden.isdisjoint(builder.OBS_FIELDS))

    def test_capacity_aware_round_robin_is_exact(self) -> None:
        quotas = builder.allocate_round_robin(
            {"small": 1, "middle": 4, "large": 8},
            target=10,
            donor_order=("small", "middle", "large"),
        )
        self.assertEqual(quotas, {"small": 1, "middle": 4, "large": 5})
        with self.assertRaises(builder.FrozenScreenError):
            builder.allocate_round_robin(
                {"a": 1, "b": 1}, target=3, donor_order=("a", "b")
            )

    def test_two_pass_selection_is_balanced_unique_and_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "cell_manifest.tsv.gz"
            self._write_manifest(manifest)
            first, first_audit = builder.select_cells(manifest, cells_per_class=6)
            second, second_audit = builder.select_cells(manifest, cells_per_class=6)

        self.assertEqual(first, second)
        self.assertEqual(first_audit, second_audit)
        self.assertEqual(len(first), 30)
        self.assertEqual(
            Counter(str(row["broad_label"]) for row in first),
            Counter({broad: 6 for broad in builder.ONTOLOGY}),
        )
        self.assertEqual(
            [int(row["source_row_index"]) for row in first],
            sorted(int(row["source_row_index"]) for row in first),
        )
        row_ids = [str(row["row_id"]) for row in first]
        self.assertEqual(len(row_ids), len(set(row_ids)))
        for broad in builder.ONTOLOGY:
            self.assertEqual(first_audit["quotas"][broad][f"{broad}-small"], 1)
            self.assertEqual(first_audit["quotas"][broad][f"{broad}-large"], 5)
            class_ranks = sorted(
                int(row["selection_class_rank"])
                for row in first
                if row["broad_label"] == broad
            )
            self.assertEqual(class_ranks, list(range(6)))

    def test_manifest_identity_and_class_capacity_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "cell_manifest.tsv.gz"
            self._write_manifest(manifest, small=1, large=1)
            with self.assertRaises(builder.FrozenScreenError):
                builder.select_cells(manifest, cells_per_class=3)

            rows = list(builder._manifest_rows(manifest))
            fields = sorted(builder.REQUIRED_MANIFEST_FIELDS)
            rows[0][1]["cell_id"] = "not-the-compound-identity"
            with gzip.open(manifest, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(row for _, row in rows)
            with self.assertRaises(builder.FrozenScreenError):
                builder.scan_capacities(manifest)

    def test_contract_lock_checks_file_and_canonical_content_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            lock = Path(directory) / "contract_lock.json"
            payload = {
                "atlas_realpath": "/synthetic/atlas.h5ad",
                "atlas_size_bytes": 10,
                "cell_order_sha256": "c" * 64,
                "gene_order_sha256": "g" * 64,
                "manifest_sha256": {
                    "assay_manifest.tsv": "a" * 64,
                    "cell_manifest.tsv.gz": "m" * 64,
                    "donor_manifest.tsv": "d" * 64,
                    "library_manifest.tsv": "l" * 64,
                },
                "production_ready": True,
                "raw_counts": {"nnz": 1, "shape": [1, 1]},
                "raw_counts_sha256": "r" * 64,
                "schema_version": "masld-cl-contract-v1",
            }
            payload["lock_sha256"] = builder.canonical_sha256(payload)
            lock.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
            file_sha = builder.sha256_file(lock)
            observed = builder.load_contract(
                lock,
                expected_file_sha256=file_sha,
                expected_lock_sha256=str(payload["lock_sha256"]),
                expected_raw_counts_sha256="r" * 64,
                expected_cell_order_sha256="c" * 64,
                expected_gene_order_sha256="g" * 64,
            )
            self.assertTrue(observed["production_ready"])
            with self.assertRaises(builder.FrozenScreenError):
                builder.load_contract(
                    lock,
                    expected_file_sha256=file_sha,
                    expected_lock_sha256=str(payload["lock_sha256"]),
                    expected_raw_counts_sha256="x" * 64,
                    expected_cell_order_sha256="c" * 64,
                    expected_gene_order_sha256="g" * 64,
                )

    def test_interval_copy_preserves_values_and_hashes_full_source(self) -> None:
        source = np.arange(10, dtype=np.int32)
        output = np.empty(5, dtype=np.int32)
        digest = hashlib.sha256()
        builder._copy_selected_intervals(
            dataset=source,
            output=output,
            source_starts=np.asarray([1, 6]),
            source_ends=np.asarray([3, 9]),
            target_starts=np.asarray([0, 2]),
            digest=digest,
            component="data",
            n_features=20,
        )
        np.testing.assert_array_equal(output, np.asarray([1, 2, 6, 7, 8]))
        self.assertEqual(digest.hexdigest(), hashlib.sha256(source.tobytes()).hexdigest())

    def test_receipts_follow_independent_verifier_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "fixture"
            root.mkdir()
            (root / "payload.tsv").write_text("row_id\ncell-1\n", encoding="utf-8")
            builder.freeze_output_tree(root, manifest_sha256="f" * 64)
            try:
                artifacts = json.loads(
                    (root / "ARTIFACTS.json").read_text(encoding="utf-8")
                )
                complete = json.loads(
                    (root / "COMPLETE").read_text(encoding="utf-8")
                )
                self.assertEqual(artifacts["schema_version"], "masld-bench-artifacts-v1")
                self.assertEqual(artifacts["metadata"]["row_count"], 50_000)
                self.assertFalse(artifacts["metadata"]["sealed_outcomes_read"])
                self.assertEqual(
                    complete["manifest_sha256"],
                    builder.sha256_file(root / "ARTIFACTS.json"),
                )
                for record in artifacts["artifacts"]:
                    path = root / record["path"]
                    self.assertEqual(builder.sha256_file(path), record["sha256"])
                    self.assertEqual(path.stat().st_size, record["size_bytes"])
            finally:
                root.chmod(0o700)
                for path in root.iterdir():
                    path.chmod(0o600)

    def test_wrapper_is_generic_cpu_nslab_and_non_submitting(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(line for line in text.splitlines() if line.startswith("#SBATCH"))
        self.assertIn("--job-name=model-data-050", header)
        self.assertIn("--partition=cpu", header)
        self.assertIn("--account=nslab", header)
        self.assertIn("--qos=nslab", header)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("--array", header)
        self.assertNotIn("sbatch ", text)
        self.assertIn("test_build_resource_atlas_frozen_screen_50000", text)
        self.assertIn("ee67de1c51116285f52cbc2c31a69b468b9ebb585a36800311bdb9cc117ee4ed", text)
        self.assertIn("ef1ffcca70d396aaf9f6f9da0b59bb113ed521aa9f8c989cd12fb194c316805f", text)


if __name__ == "__main__":
    unittest.main()
