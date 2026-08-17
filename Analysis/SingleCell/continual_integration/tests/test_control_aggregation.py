from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from masld_cl.control_aggregation import aggregate_control_metrics
from masld_cl.contracts import sha256_path


class TestControlAggregation(unittest.TestCase):
    def test_six_models_reduce_by_prespecified_rules(self):
        models = ["all_lineage", "H", "M", "F", "C", "T"]
        metrics = [
            "reference_macro_f1_change_ci_low", "reference_neighborhood_jaccard_loss",
            "shift_control", "shift_control_standard_error",
        ]
        with tempfile.TemporaryDirectory() as directory:
            inputs = []
            for index, model in enumerate(models):
                path = Path(directory) / f"{index}.tsv"
                with path.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=[
                        "setting_id", "ewc_lambda", "replay_fraction", "model_kind",
                        "metric", "scope", "value",
                    ], delimiter="\t")
                    writer.writeheader()
                    values = [-0.001 * index, 0.01 * index, 1.0 + index, 0.02 + index * 0.01]
                    for metric, value in zip(metrics, values):
                        writer.writerow({
                            "setting_id": "s", "ewc_lambda": 1, "replay_fraction": .2,
                            "model_kind": model, "metric": metric, "scope": "x", "value": value,
                        })
                reference = Path(directory) / f"reference_{index}.json"
                candidate = Path(directory) / f"candidate_{index}.json"
                reference_embedding = Path(directory) / f"reference_embedding_{index}.json"
                candidate_embedding = Path(directory) / f"candidate_embedding_{index}.json"
                reference_record = Path(directory) / f"reference_record_{index}.json"
                candidate_record = Path(directory) / f"candidate_record_{index}.json"
                reference.write_text("{}\n")
                candidate.write_text("{}\n")
                reference_embedding.write_text("{}\n")
                candidate_embedding.write_text("{}\n")
                reference_record.write_text("{}\n")
                candidate_record.write_text("{}\n")
                details = {
                    "schema_version": "masld-cl-control-evaluation-v1",
                    "config_sha256": "cfg",
                    "metrics_realpath": str(path.resolve()),
                    "metrics_sha256": sha256_path(path),
                    "candidate_seed": 17,
                    "reference_run_manifest": str(reference.resolve()),
                    "reference_run_manifest_sha256": sha256_path(reference),
                    "candidate_run_manifest": str(candidate.resolve()),
                    "candidate_run_manifest_sha256": sha256_path(candidate),
                    "reference_embedding": str(reference_embedding.resolve()),
                    "candidate_embedding": str(candidate_embedding.resolve()),
                    "reference_execution_record": str(reference_record.resolve()),
                    "reference_execution_record_sha256": sha256_path(reference_record),
                    "candidate_execution_record": str(candidate_record.resolve()),
                    "candidate_execution_record_sha256": sha256_path(candidate_record),
                }
                path.with_suffix(path.suffix + ".details.json").write_text(
                    json.dumps(details)
                )
                inputs.append(path)
            output = Path(directory) / "aggregate.tsv"
            config_path = Path(directory) / "config.json"
            config_path.write_text("{}\n")
            with patch("masld_cl.control_aggregation.verify_execution_record", return_value={}), patch(
                "masld_cl.control_aggregation.require_execution_ownership"
            ):
                rows = aggregate_control_metrics(
                    {"lineages": models[1:], "_config_sha256": "cfg", "_config_path": str(config_path), "screen": {"seed": 17}},
                    inputs, output
                )
            values = {row["metric"]: row["value"] for row in rows}
            self.assertEqual(values["reference_macro_f1_change_ci_low"], -0.005)
            self.assertEqual(values["reference_neighborhood_jaccard_loss"], 0.05)
            self.assertEqual(values["shift_control"], 3.5)
            self.assertEqual(values["shift_control_standard_error"], 0.07)
