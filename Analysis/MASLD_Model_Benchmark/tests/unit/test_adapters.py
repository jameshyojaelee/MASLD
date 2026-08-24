from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

from masld_bench.adapters.base import (
    AdapterAction,
    AdapterError,
    AdapterReceipt,
    SubprocessAdapter,
    withheld_input_roles_for_action,
)
from masld_bench.artifacts import (
    canonical_hash,
    sha256_file,
    verify_frozen_tree,
    write_json_exclusive,
)


def _fixture_run_spec() -> dict[str, object]:
    inputs = [
        {
            "path": "/fixture/environment.lock",
            "sha256": "e" * 64,
            "size_bytes": 1,
            "media_type": "text/plain",
            "role": "environment:fixture_cpu",
        }
    ]
    return {
        "runtime_id": "fixture_cpu",
        "stage": "smoke",
        "dataset_ids": ["development_fixture"],
        "action": ["prepare", "predict"],
        "inputs": inputs,
        "metadata": {
            "input_owner_by_role": {"environment:fixture_cpu": None},
            "full_input_inventory_sha256": canonical_hash(inputs),
        },
    }


def _fixture_dataset_access() -> dict[str, object]:
    run_spec = _fixture_run_spec()
    return {
        "fit_dataset_ids": [],
        "development_prediction_dataset_ids": ["development_fixture"],
        "prediction_first_stress_dataset_ids": [],
        "sealed_prediction_dataset_ids": [],
        "action_dataset_ids": ["development_fixture"],
        "dataset_routing": {
            "development_fixture": {
                "task_partition": "development_prediction",
                "fit_allowed": False,
                "prediction_first_policy": None,
            }
        },
        "full_input_inventory_sha256": run_spec["metadata"][  # type: ignore[index]
            "full_input_inventory_sha256"
        ],
        "action_input_inventory_sha256": canonical_hash(run_spec["inputs"]),
    }


class AdapterReceiptTests(unittest.TestCase):
    def test_action_specific_outcome_role_is_withheld_and_environment_is_not(self) -> None:
        roles = frozenset(
            {
                "environment:fixture_cpu",
                "dataset_view_data:paired_fixture",
                "dataset_view_topology:paired_fixture",
            }
        )
        request = {
            "withheld_input_roles_by_action": {
                "fit": ["dataset_view_data:paired_fixture"],
                "predict": ["dataset_view_data:paired_fixture"],
            }
        }
        self.assertEqual(
            withheld_input_roles_for_action(
                request=request,
                action=AdapterAction.PREDICT,
                available_roles=roles,
            ),
            frozenset({"dataset_view_data:paired_fixture"}),
        )
        frozen_request = {
            "withheld_input_roles_by_action": {
                "fit": ("dataset_view_data:paired_fixture",),
            }
        }
        self.assertEqual(
            withheld_input_roles_for_action(
                request=frozen_request,
                action=AdapterAction.FIT,
                available_roles=roles,
            ),
            frozenset({"dataset_view_data:paired_fixture"}),
        )
        with self.assertRaisesRegex(AdapterError, "environment lock"):
            withheld_input_roles_for_action(
                request={
                    "withheld_input_roles_by_action": {
                        "predict": ["environment:fixture_cpu"]
                    }
                },
                action=AdapterAction.PREDICT,
                available_roles=roles,
            )

    def test_fixture_runs_through_the_real_subprocess_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            request = root / "request.json"
            output = root / "output"
            run_id = "f" * 64
            write_json_exclusive(
                request,
                {
                    "schema_version": "masld-bench-adapter-request-v1",
                    "action": "predict",
                    "run_id": run_id,
                    "row_ids": ["donor-1", "donor-2"],
                    "run_spec": _fixture_run_spec(),
                    **_fixture_dataset_access(),
                },
            )
            adapter = SubprocessAdapter(
                [sys.executable, "-m", "masld_bench.adapters.stub"],
                timeout_seconds=60,
            )
            receipt = adapter.predict(request, output)
            self.assertEqual(receipt.run_id, run_id)
            self.assertEqual(receipt.status, "complete")
            self.assertEqual(len(receipt.artifacts), 1)
            manifest = verify_frozen_tree(output)
            self.assertEqual(
                {item["path"] for item in manifest["artifacts"]},
                {
                    "adapter_receipt.json",
                    "adapter_stderr.log",
                    "adapter_stdout.log",
                    "predictions.json",
                },
            )

    def test_receipt_requires_full_run_hash_and_exact_artifact_size(self) -> None:
        with self.assertRaisesRegex(AdapterError, "full lowercase SHA-256"):
            AdapterReceipt.from_mapping(
                {
                    "schema_version": "masld-bench-adapter-receipt-v1",
                    "action": "probe",
                    "run_id": "short",
                    "status": "complete",
                    "artifacts": [],
                }
            )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / "embedding.jsonl"
            artifact.write_text("{}\n", encoding="utf-8")
            receipt = AdapterReceipt.from_mapping(
                {
                    "schema_version": "masld-bench-adapter-receipt-v1",
                    "action": "predict",
                    "run_id": "a" * 64,
                    "status": "complete",
                    "artifacts": [
                        {
                            "path": artifact.name,
                            "sha256": sha256_file(artifact),
                            "size_bytes": artifact.stat().st_size + 1,
                        }
                    ],
                }
            )
            with self.assertRaisesRegex(AdapterError, "size mismatch"):
                receipt.validate_artifacts(root)

    def test_receipt_rejects_unknown_and_evaluation_fields(self) -> None:
        base = {
            "schema_version": "masld-bench-adapter-receipt-v1",
            "action": "probe",
            "run_id": "a" * 64,
            "status": "complete",
            "artifacts": [],
        }
        with self.assertRaisesRegex(AdapterError, "unknown fields"):
            AdapterReceipt.from_mapping({**base, "surprise": True})
        with self.assertRaisesRegex(AdapterError, "evaluation fields"):
            AdapterReceipt.from_mapping({**base, "metadata": {"macro_f1": 0.9}})
        with self.assertRaisesRegex(AdapterError, "invalid adapter receipt action"):
            AdapterReceipt.from_mapping({**base, "action": "score"})

    def test_symlink_artifact_cannot_escape_output_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            root = parent / "output"
            root.mkdir()
            outside = parent / "outside.bin"
            outside.write_bytes(b"not an adapter artifact")
            os.symlink(outside, root / "link.bin")
            receipt = AdapterReceipt.from_mapping(
                {
                    "schema_version": "masld-bench-adapter-receipt-v1",
                    "action": "export",
                    "run_id": "b" * 64,
                    "status": "complete",
                    "artifacts": [
                        {
                            "path": "link.bin",
                            "sha256": sha256_file(outside),
                            "size_bytes": outside.stat().st_size,
                        }
                    ],
                }
            )
            with self.assertRaisesRegex(AdapterError, "symlink"):
                receipt.validate_artifacts(root)

    def test_subprocess_rejects_symlinked_output_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            request = root / "request.json"
            write_json_exclusive(
                request,
                {
                    "schema_version": "masld-bench-adapter-request-v1",
                    "action": "predict",
                    "run_id": "c" * 64,
                    "row_ids": ["donor-1"],
                    "run_spec": _fixture_run_spec(),
                    **_fixture_dataset_access(),
                },
            )
            actual = root / "actual"
            linked = root / "linked"
            os.symlink(actual, linked)
            adapter = SubprocessAdapter(
                [sys.executable, "-m", "masld_bench.adapters.stub"],
                timeout_seconds=60,
            )
            with self.assertRaisesRegex(AdapterError, "invalid"):
                adapter.predict(request, linked)

    def test_fit_request_rejects_development_or_prediction_first_data(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            request = root / "request.json"
            access = _fixture_dataset_access()
            access["action_dataset_ids"] = ["development_fixture"]
            write_json_exclusive(
                request,
                {
                    "schema_version": "masld-bench-adapter-request-v1",
                    "action": "fit",
                    "run_id": "d" * 64,
                    "run_spec": _fixture_run_spec(),
                    **access,
                },
            )
            adapter = SubprocessAdapter(
                [sys.executable, "-m", "masld_bench.adapters.stub"],
                timeout_seconds=60,
            )
            with self.assertRaisesRegex(AdapterError, "frozen action scope"):
                adapter.fit(request, root / "output")

    def test_fit_request_cannot_expose_development_artifact_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            request = root / "request.json"
            inputs = [
                {
                    "path": "/fixture/environment.lock",
                    "sha256": "e" * 64,
                    "size_bytes": 1,
                    "role": "environment:fixture_cpu",
                },
                {
                    "path": "/private/train.tsv",
                    "sha256": "a" * 64,
                    "size_bytes": 1,
                    "role": "dataset:train",
                },
                {
                    "path": "/private/development.tsv",
                    "sha256": "b" * 64,
                    "size_bytes": 1,
                    "role": "dataset:development",
                },
            ]
            owners = {
                "environment:fixture_cpu": None,
                "dataset:train": "train",
                "dataset:development": "development",
            }
            run_spec = {
                "runtime_id": "fixture_cpu",
                "stage": "smoke",
                "dataset_ids": ["train", "development"],
                "inputs": inputs,
                "metadata": {
                    "input_owner_by_role": owners,
                    "full_input_inventory_sha256": canonical_hash(inputs),
                },
            }
            write_json_exclusive(
                request,
                {
                    "schema_version": "masld-bench-adapter-request-v1",
                    "action": "fit",
                    "run_id": "e" * 64,
                    "run_spec": run_spec,
                    "fit_dataset_ids": ["train"],
                    "development_prediction_dataset_ids": ["development"],
                    "prediction_first_stress_dataset_ids": [],
                    "sealed_prediction_dataset_ids": [],
                    "action_dataset_ids": ["train"],
                    "dataset_routing": {
                        "train": {
                            "task_partition": "fit",
                            "fit_allowed": True,
                            "prediction_first_policy": None,
                        },
                        "development": {
                            "task_partition": "development_prediction",
                            "fit_allowed": False,
                            "prediction_first_policy": None,
                        },
                    },
                    "full_input_inventory_sha256": canonical_hash(inputs),
                    "action_input_inventory_sha256": canonical_hash(inputs),
                },
            )
            adapter = SubprocessAdapter(
                [sys.executable, "-m", "masld_bench.adapters.stub"],
                timeout_seconds=60,
            )
            with self.assertRaisesRegex(AdapterError, "incorrect action inventory"):
                adapter.fit(request, root / "output")

    def test_prior_action_output_is_recursively_verified(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            attempt = Path(temporary) / "attempt-001"
            requests = attempt / "requests"
            requests.mkdir(parents=True)
            first_request = requests / "001-prepare.json"
            run_id = "f" * 64
            write_json_exclusive(
                first_request,
                {
                    "schema_version": "masld-bench-adapter-request-v1",
                    "action": "prepare",
                    "run_id": run_id,
                    "run_spec": _fixture_run_spec(),
                    "prior_action_outputs": [],
                    **_fixture_dataset_access(),
                },
            )
            adapter = SubprocessAdapter(
                [sys.executable, "-m", "masld_bench.adapters.stub"],
                timeout_seconds=60,
            )
            first_output = attempt / "adapter_actions" / "001-prepare"
            adapter.prepare(first_request, first_output)
            prior = {
                "action": "prepare",
                "status": "complete",
                "output_path": "adapter_actions/001-prepare",
                "adapter_receipt_sha256": sha256_file(
                    first_output / "adapter_receipt.json"
                ),
                "output_manifest_sha256": sha256_file(
                    first_output / "ARTIFACTS.json"
                ),
            }
            second_request = requests / "002-predict.json"
            write_json_exclusive(
                second_request,
                {
                    "schema_version": "masld-bench-adapter-request-v1",
                    "action": "predict",
                    "run_id": run_id,
                    "row_ids": ["donor-1"],
                    "run_spec": _fixture_run_spec(),
                    "prior_action_outputs": [prior],
                    **_fixture_dataset_access(),
                },
            )
            second_output = attempt / "adapter_actions" / "002-predict"
            self.assertEqual(
                adapter.predict(second_request, second_output).status, "complete"
            )

            forged = dict(prior)
            forged["output_manifest_sha256"] = "0" * 64
            second_request.write_text(
                json.dumps(
                    {
                        "schema_version": "masld-bench-adapter-request-v1",
                        "action": "predict",
                        "run_id": run_id,
                        "row_ids": ["donor-1"],
                        "run_spec": _fixture_run_spec(),
                        "prior_action_outputs": [forged],
                        **_fixture_dataset_access(),
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(AdapterError, "binding changed"):
                adapter.predict(
                    second_request, attempt / "adapter_actions" / "forged-output"
                )


if __name__ == "__main__":
    unittest.main()
