from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from scripts.build_gse267145_histology_fit_views import (
    QUERY_FIELDS,
    TRAINING_FIELDS,
    build_fit_views,
)


ROOT = Path(__file__).resolve().parents[2]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class GSE267145FitViewTests(unittest.TestCase):
    def test_real_fold_views_exclude_held_labels_and_sex(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "views"
            result = build_fit_views(
                participant_axis=ROOT / "executions/model-data-064-21079902/fixture/molecular/participant_axis.tsv",
                outer_folds=ROOT / "executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv",
                endpoints=ROOT / "executions/model-data-061-21079623/activation/participant_endpoints.tsv",
                inner_roster=ROOT / "executions/model-data-066-21080119/validation/inner_fold_roster.tsv",
                output=output,
            )
            self.assertEqual(result["participants"], 99)
            for outer_fold in range(5):
                training = rows(output / f"outer_{outer_fold}/training_endpoints.tsv")
                query = rows(output / f"outer_{outer_fold}/query_participants.tsv")
                self.assertEqual(tuple(training[0]), TRAINING_FIELDS)
                self.assertEqual(tuple(query[0]), QUERY_FIELDS)
                self.assertFalse({row["participant_id"] for row in training} & {row["participant_id"] for row in query})
                self.assertEqual(len(training) + len(query), 99)
                self.assertNotIn("recorded_sex", training[0])
                self.assertNotIn("stage5", training[0])
                self.assertNotIn("stage3", query[0])


if __name__ == "__main__":
    unittest.main()
