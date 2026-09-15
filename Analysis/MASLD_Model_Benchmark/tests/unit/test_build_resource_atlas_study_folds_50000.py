#!/usr/bin/env python3
"""Unit checks for the 50,000-cell study-held-out split contract."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file
from scripts import build_resource_atlas_study_folds_50000 as builder


ROOT = Path(__file__).parents[2]
WRAPPER = ROOT / "slurm/build_resource_atlas_study_folds_50000.sbatch"


class StudyFold50000Tests(unittest.TestCase):
    def test_outer_assignment_depends_only_on_study_size_and_id(self) -> None:
        counts = {"s7": 4, "s6": 5, "s5": 6, "s4": 7, "s3": 8, "s2": 9, "s1": 10}
        expected = {"s1": 0, "s2": 1, "s3": 2, "s4": 3, "s5": 4, "s6": 4, "s7": 4}
        self.assertEqual(builder.assign_outer_folds(counts), expected)
        self.assertEqual(builder.assign_outer_folds(dict(reversed(list(counts.items())))), expected)

    def test_inner_assignment_is_balanced_and_outer_specific(self) -> None:
        donors = [f"donor-{index:02d}" for index in range(23)]
        first = builder.assign_inner_donors(donors, outer_fold=0)
        second = builder.assign_inner_donors(reversed(donors), outer_fold=0)
        other = builder.assign_inner_donors(donors, outer_fold=1)
        self.assertEqual(first, second)
        counts = [list(first.values()).count(fold) for fold in range(5)]
        self.assertLessEqual(max(counts) - min(counts), 1)
        self.assertNotEqual(first, other)

    @staticmethod
    def _source(root: Path) -> tuple[dict[str, int], dict[str, int]]:
        counts = {"s1": 10, "s2": 9, "s3": 8, "s4": 7, "s5": 6, "s6": 5, "s7": 4}
        donor_counts: dict[str, int] = {}
        labels = list(builder.EXPECTED_LABELS)
        with (root / "selection.tsv").open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("row_id", "donor_id", "dataset", "broad_label"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            ordinal = 0
            for dataset, count in counts.items():
                donor_counts[dataset] = 2
                for index in range(count):
                    writer.writerow(
                        {
                            "row_id": f"row-{ordinal}",
                            "donor_id": f"donor-{dataset}-{index % 2}",
                            "dataset": dataset,
                            "broad_label": labels[index % len(labels)],
                        }
                    )
                    ordinal += 1
        freeze_tree(
            root,
            {
                "artifact_class": "resource_atlas_frozen_screen_fixture",
                "subset_id": builder.DATASET_VIEW_ID,
                "row_count": sum(counts.values()),
                "sealed_outcomes_read": False,
                "histology_read": False,
            },
        )
        return counts, donor_counts

    def test_build_freezes_donor_safe_study_folds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source"
            output = Path(directory) / "output"
            source.mkdir()
            counts, donor_counts = self._source(source)
            digest = sha256_file(source / "ARTIFACTS.json")
            builder.build(
                source=source,
                output=output,
                expected_source_artifacts_sha256=digest,
                expected_rows=sum(counts.values()),
                expected_donors=14,
                expected_dataset_counts=counts,
                expected_dataset_donors=donor_counts,
            )
            verified = verify_frozen_tree(output)
            manifest = json.loads((output / "split_manifest.json").read_text())
            with (output / "row_outer_folds.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            with (output / "inner_donor_folds.tsv").open(newline="") as handle:
                inner = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(verified["metadata"]["rows"], sum(counts.values()))
        self.assertFalse(verified["metadata"]["target_labels_used_for_assignment"])
        self.assertTrue(manifest["post_assignment_validation"]["target_labels_read"])
        dataset_folds: dict[str, set[str]] = {}
        for row in rows:
            dataset_folds.setdefault(row["dataset"], set()).add(row["outer_fold"])
        self.assertTrue(all(len(values) == 1 for values in dataset_folds.values()))
        donor_outer = {row["donor_id"]: row["outer_fold"] for row in rows}
        self.assertTrue(
            all(donor_outer[row["donor_id"]] != row["held_outer_fold"] for row in inner)
        )

    def test_production_contract_and_wrapper_are_fixed(self) -> None:
        self.assertEqual(sum(builder.EXPECTED_DATASET_COUNTS.values()), 50_000)
        self.assertEqual(sum(builder.EXPECTED_DATASET_DONORS.values()), 102)
        self.assertEqual(
            builder.assign_outer_folds(builder.EXPECTED_DATASET_COUNTS),
            {
                "GSE202379": 0,
                "GSE244832": 1,
                "Liver_Atlas": 2,
                "GSE136103": 3,
                "GSE174748": 4,
                "GSE185477": 4,
                "GSE189600": 4,
            },
        )
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(line for line in text.splitlines() if line.startswith("#SBATCH"))
        self.assertIn("--job-name=model-data-052", header)
        self.assertIn("--partition=cpu", header)
        self.assertIn("--qos=nslab", header)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("--array", header)
        self.assertFalse(
            any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
        )
        self.assertIn("672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509", text)


if __name__ == "__main__":
    unittest.main()
