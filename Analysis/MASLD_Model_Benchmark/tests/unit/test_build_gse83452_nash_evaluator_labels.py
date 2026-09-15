from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts.build_gse83452_nash_evaluator_labels import (
    EvaluatorLabelError,
    build_labels,
    verify_predictions_are_frozen,
)

COLUMNS = [
    "cohort_family_id",
    "series",
    "participant_id",
    "sample_accession",
    "timepoint",
    "repeat_topology",
    "nash_status",
    "age",
    "sex",
    "intervention",
]


def write_records(
    path: Path,
    *,
    nash: int = 104,
    no_nash: int = 44,
    undefined: int = 4,
    followup: int = 79,
) -> None:
    rows = []
    index = 0
    for status, count in (("nash", nash), ("no_nash", no_nash), ("undefined", undefined)):
        for _ in range(count):
            rows.append(
                {
                    "cohort_family_id": "antwerp_inserm_shared",
                    "series": "GSE83452",
                    "participant_id": f"gse83452::GSM{2203254 + index}",
                    "sample_accession": f"GSM{2203254 + index}",
                    "timepoint": "baseline",
                    "repeat_topology": "baseline_only_or_pair_anchor",
                    "nash_status": status,
                    "age": "55",
                    "sex": "female",
                    "intervention": "BS",
                }
            )
            index += 1
    for _ in range(followup):
        rows.append(
            {
                "cohort_family_id": "antwerp_inserm_shared",
                "series": "GSE83452",
                "participant_id": f"gse83452::GSM{2203254 + index}",
                "sample_accession": f"GSM{2203254 + index}",
                "timepoint": "follow-up",
                "repeat_topology": "followup_only_unpaired",
                "nash_status": "nash",
                "age": "55",
                "sex": "male",
                "intervention": "Diet",
            }
        )
        index += 1
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


class LabelTests(unittest.TestCase):
    def test_baseline_selection_and_undefined_handling(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "records.tsv"
            write_records(path)
            labels, receipt = build_labels(path)
            self.assertEqual(len(labels), 148)
            self.assertEqual(receipt["baseline_records"], 152)
            self.assertEqual(receipt["endpoint_evaluable"], 148)
            self.assertEqual(receipt["class_counts"], {"no_nash": 44, "nash": 104})
            self.assertEqual(receipt["undefined_retained_as_missing"], 4)
            self.assertFalse(receipt["undefined_mapped_to_no_nash"])

    def test_undefined_never_appears_in_the_label_table(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "records.tsv"
            write_records(path)
            labels, _ = build_labels(path)
            self.assertNotIn("undefined", {row["nash_status"] for row in labels})

    def test_age_sex_and_intervention_are_never_projected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "records.tsv"
            write_records(path)
            labels, _ = build_labels(path)
            for row in labels:
                self.assertEqual(sorted(row), ["nash_status", "row_id"])

    def test_a_shifted_undefined_census_is_rejected(self) -> None:
        """Mapping one undefined to no-NASH would produce exactly this census."""

        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "records.tsv"
            write_records(path, no_nash=45, undefined=3)
            with self.assertRaises(EvaluatorLabelError):
                build_labels(path)

    def test_a_followup_record_never_reaches_the_endpoint(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "records.tsv"
            write_records(path)
            labels, receipt = build_labels(path)
            self.assertEqual(receipt["timepoint"], "baseline")
            self.assertEqual(len(labels), 148)

    def test_an_off_roster_status_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "records.tsv"
            write_records(path)
            text = path.read_text(encoding="utf-8").replace(
                "\tno_nash\t", "\tborderline\t", 1
            )
            path.write_text(text, encoding="utf-8")
            with self.assertRaises(EvaluatorLabelError):
                build_labels(path)


class OrderingLockTests(unittest.TestCase):
    def _bundle(self, base: Path, *, label_blind: bool = True) -> Path:
        root = base / "campaign" / "predictions" / "model"
        root.mkdir(parents=True)
        (root / "prediction_bundle.json").write_text(
            json.dumps(
                {
                    "model_id": "per_array_rank_elastic_net",
                    "metadata": {
                        "external_labels_read": False if label_blind else True,
                        "prediction_frozen_before_evaluator_label_join": True,
                    },
                }
            ),
            encoding="utf-8",
        )
        (root / "ARTIFACTS.json").write_text("{}", encoding="utf-8")
        (root / "COMPLETE").write_text("{}", encoding="utf-8")
        return root

    def test_labels_are_refused_when_a_committed_sha_differs(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = self._bundle(Path(value))
            with self.assertRaises(EvaluatorLabelError):
                verify_predictions_are_frozen([root], ["0" * 64])

    def test_labels_are_refused_when_a_bundle_is_not_label_blind(self) -> None:
        import hashlib

        with tempfile.TemporaryDirectory() as value:
            root = self._bundle(Path(value), label_blind=False)
            digest = hashlib.sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()
            with self.assertRaises(EvaluatorLabelError):
                verify_predictions_are_frozen([root], [digest])

    def test_an_unfrozen_bundle_is_refused(self) -> None:
        import hashlib

        with tempfile.TemporaryDirectory() as value:
            root = self._bundle(Path(value))
            digest = hashlib.sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()
            (root / "COMPLETE").unlink()
            with self.assertRaises(EvaluatorLabelError):
                verify_predictions_are_frozen([root], [digest])


    def test_two_arms_of_one_model_are_recorded_separately(self) -> None:
        import hashlib

        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            roots = []
            for arm in ("primary", "sensitivity"):
                root = base / arm / "predictions" / "model"
                root.mkdir(parents=True)
                (root / "prediction_bundle.json").write_text(
                    json.dumps(
                        {
                            "model_id": "per_array_rank_elastic_net",
                            "metadata": {
                                "external_labels_read": False,
                                "prediction_frozen_before_evaluator_label_join": True,
                            },
                        }
                    ),
                    encoding="utf-8",
                )
                (root / "ARTIFACTS.json").write_text(f'{{"arm": "{arm}"}}', encoding="utf-8")
                (root / "COMPLETE").write_text("{}", encoding="utf-8")
                roots.append(root)
            digests = [
                hashlib.sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()
                for root in roots
            ]
            observed = verify_predictions_are_frozen(roots, digests)
            self.assertEqual(len(observed), 2)


if __name__ == "__main__":
    unittest.main()
