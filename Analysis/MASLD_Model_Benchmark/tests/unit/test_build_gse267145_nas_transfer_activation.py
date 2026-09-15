from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from scripts.build_gse267145_nas_transfer_activation import (
    NasActivationError,
    project_target_roster,
    read_feature_axis,
)


def write_axis(path: Path, ids) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(("feature_index", "stable_gene_id"))
        for i, v in enumerate(ids):
            w.writerow((i, v))


def write_roster(path: Path, extra=None) -> None:
    fields = ["participant_index", "participant_id", "rna_source_sample_accession",
              "pairing", "rna_observation_state"]
    row = {"participant_index": "0", "participant_id": "ABJ468",
           "rna_source_sample_accession": "GSM8314797",
           "pairing": "same_sample_different_aliquot",
           "rna_observation_state": "observed"}
    if extra:
        fields += list(extra)
        row.update(extra)
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(row)


class FeatureAxisTests(unittest.TestCase):
    def test_reads_stable_gene_ids(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "axis.tsv"
            write_axis(p, ["ENSG2", "ENSG1"])
            self.assertEqual(read_feature_axis(p), ["ENSG2", "ENSG1"])

    def test_versioned_id_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "axis.tsv"
            write_axis(p, ["ENSG1.4"])
            with self.assertRaises(NasActivationError):
                read_feature_axis(p)

    def test_duplicate_id_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "axis.tsv"
            write_axis(p, ["ENSG1", "ENSG1"])
            with self.assertRaises(NasActivationError):
                read_feature_axis(p)


class TargetRosterFirewallTests(unittest.TestCase):
    """The NASH-CRN component sum is the evaluator's; it must not be projectable."""

    def test_clean_roster_projects(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "axis.tsv"
            write_roster(p)
            rows = project_target_roster(p)
            self.assertEqual(rows[0]["participant_id"], "ABJ468")
            self.assertNotIn("nash_crn_component_sum", rows[0])

    def test_each_forbidden_outcome_column_raises(self) -> None:
        for column in ("nash_crn_component_sum", "steatosis", "ballooning",
                       "lobular_inflammation", "fibrosis", "stage3", "nas_score",
                       "recorded_sex"):
            with self.subTest(column=column):
                with tempfile.TemporaryDirectory() as d:
                    p = Path(d) / "axis.tsv"
                    write_roster(p, extra={column: "3"})
                    with self.assertRaises(NasActivationError):
                        project_target_roster(p)

    def test_undeclared_extra_column_is_dropped_not_carried(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "axis.tsv"
            write_roster(p, extra={"some_future_column": "x"})
            rows = project_target_roster(p)
            self.assertNotIn("some_future_column", rows[0])


if __name__ == "__main__":
    unittest.main()
