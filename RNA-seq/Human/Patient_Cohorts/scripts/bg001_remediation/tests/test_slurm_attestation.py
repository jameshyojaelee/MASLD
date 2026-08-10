#!/usr/bin/env python3
from __future__ import annotations

import copy
import csv
import importlib.util
import io
import shlex
import subprocess
import tempfile
import unittest
from datetime import datetime, timedelta
from pathlib import Path
from unittest import mock


SCRIPT_DIR = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


attestation = load_module("bg001_slurm_attestation", SCRIPT_DIR / "slurm_attestation.py")


def scheduler_timestamp(value: datetime) -> str:
    return value.strftime("%Y-%m-%dT%H:%M:%S")


def ledger_timestamp(value: datetime) -> str:
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


def encode_sacct(rows: list[dict[str, str]]) -> str:
    output = io.StringIO()
    writer = csv.writer(output, delimiter="|", lineterminator="\n")
    for row in rows:
        writer.writerow([row[field] for field in attestation.SACCT_FIELDS])
    return output.getvalue()


class ExactSlurmFixture:
    """Build the exact six-stage, 20-logical-job donor execution contract."""

    def __init__(self, base: Path):
        self.root = (
            base / "bg001-fragment-v211-gencode49-20260806T055119Z"
        ).resolve()
        (self.root / "contract").mkdir(parents=True)
        (self.root / "logs").mkdir()
        (self.root / ".bg001_candidate_root").write_text("candidate\n")
        for arm in ("R0", "F_locked", "F_legacy", "F_five"):
            (self.root / "arms" / arm).mkdir(parents=True)

        self.ledger_rows = self._make_ledger_rows()
        self.write_ledger()
        self.accounting_rows = self._make_accounting_rows()

    def _make_ledger_rows(self) -> list[dict[str, str]]:
        submitted = datetime(2026, 8, 6, 0, 0, 0)
        job_ids = {
            str(spec["stage"]): str(19_560_000 + index)
            for index, spec in enumerate(attestation.STAGE_SPECS, start=1)
        }
        rows = []
        for index, spec in enumerate(attestation.STAGE_SPECS):
            dependency_stage = str(spec["dependency_stage"])
            dependency = (
                f"afterok:{job_ids[dependency_stage]}" if dependency_stage else ""
            )
            rows.append(
                {
                    "stage": str(spec["stage"]),
                    "job_id": job_ids[str(spec["stage"])],
                    "dependency": dependency,
                    "job_name": str(spec["job_name"]),
                    "partition": str(spec["partition"]),
                    "qos": str(spec["qos"]),
                    "array": str(spec["array"]),
                    "cpus_per_task": str(spec["cpus"]),
                    "memory": str(spec["memory"]),
                    "time_limit": attestation.LEDGER_TIME_LIMIT,
                    "submitted_utc": ledger_timestamp(submitted + timedelta(minutes=index)),
                }
            )
        return rows

    def write_ledger(self) -> None:
        with (self.root / "contract/slurm_submission.tsv").open(
            "w", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=attestation.LEDGER_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(self.ledger_rows)

    def _make_accounting_rows(self) -> list[dict[str, str]]:
        ledger_by_stage = {row["stage"]: row for row in self.ledger_rows}
        workflow_start = datetime(2026, 8, 6, 1, 0, 0)
        raw_job_id = 90_000_000
        rows = []
        for stage_index, spec in enumerate(attestation.STAGE_SPECS):
            ledger = ledger_by_stage[str(spec["stage"])]
            ledger_submit = datetime.strptime(
                ledger["submitted_utc"], "%Y-%m-%dT%H:%M:%SZ"
            )
            stage_start = workflow_start + timedelta(minutes=60 * stage_index)
            stage_end = stage_start + timedelta(minutes=20)
            task_count = int(spec["tasks"])
            for task_index in range(task_count):
                task_start = stage_start + timedelta(
                    minutes=21 if task_index >= 4 else 0
                )
                task_end = task_start + timedelta(minutes=20)
                logical_id = (
                    f"{ledger['job_id']}_{task_index}"
                    if task_count > 1
                    else ledger["job_id"]
                )
                executed = bool(spec["executed"])
                rows.append(
                    {
                        "JobID": logical_id,
                        "JobIDRaw": str(raw_job_id),
                        "Cluster": attestation.EXPECTED_CLUSTER,
                        "User": attestation.EXPECTED_USER,
                        "UID": str(self.root.stat().st_uid),
                        "Account": attestation.EXPECTED_ACCOUNT,
                        "AssocID": attestation.EXPECTED_ASSOC_ID,
                        "DBIndex": str(raw_job_id + 100_000_000),
                        "JobName": str(spec["job_name"]),
                        "State": "COMPLETED" if executed else "CANCELLED by 12345",
                        "ExitCode": "0:0",
                        "Partition": str(spec["partition"]),
                        "QOS": str(spec["qos"]),
                        "ReqCPUS": str(spec["cpus"]),
                        "AllocCPUS": str(spec["cpus"]) if executed else "0",
                        "ReqMem": str(spec["memory"]),
                        "Timelimit": attestation.TIME_LIMIT,
                        "Submit": scheduler_timestamp(
                            ledger_submit - timedelta(seconds=9)
                        ),
                        "Start": scheduler_timestamp(task_start) if executed else "None",
                        "End": scheduler_timestamp(task_end if executed else stage_end),
                        "WorkDir": str(self.root),
                        "NodeList": "io001" if executed else "None assigned",
                        "StdOut": str(
                            self.root / "logs" / f"{spec['log_stem']}.out"
                        ),
                        "StdErr": str(
                            self.root / "logs" / f"{spec['log_stem']}.err"
                        ),
                        "SubmitLine": shlex.join(
                            attestation.expected_submit_tokens(self.root, ledger, spec)
                        ),
                    }
                )
                raw_job_id += 1
        return rows

    def ledger(self) -> list[dict[str, str]]:
        return attestation.read_submission_ledger(self.root)

    def record(self, logical_id: str) -> dict[str, str]:
        return next(row for row in self.accounting_rows if row["JobID"] == logical_id)


class SlurmAttestationTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.fixture = ExactSlurmFixture(Path(self.temporary_directory.name))

    def tearDown(self):
        self.temporary_directory.cleanup()

    def validate(self, rows: list[dict[str, str]] | None = None):
        accounting = self.fixture.accounting_rows if rows is None else rows
        return attestation.validate_accounting(
            self.fixture.root,
            self.fixture.ledger(),
            accounting,
        )

    def test_exact_six_stage_twenty_record_fixture_passes(self):
        raw = encode_sacct(self.fixture.accounting_rows)
        parsed = attestation.parse_sacct(raw)
        normalized = attestation.validate_accounting(
            self.fixture.root,
            self.fixture.ledger(),
            parsed,
        )

        self.assertEqual(len(normalized), 20)
        self.assertEqual(sum(row["State"] == "COMPLETED" for row in normalized), 19)
        self.assertEqual(
            sum(row["State"].startswith("CANCELLED") for row in normalized),
            1,
        )
        self.assertEqual(
            {row["stage"] for row in normalized},
            {str(spec["stage"]) for spec in attestation.STAGE_SPECS},
        )

    def test_missing_array_task_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        rows.remove(self.fixture.record("19560001_7"))

        with self.assertRaisesRegex(SystemExit, "expected 20 logical jobs, found 19"):
            self.validate(rows)

    def test_extra_array_task_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        extra = copy.deepcopy(self.fixture.record("19560001_7"))
        extra["JobID"] = "19560001_8"
        extra["JobIDRaw"] = "99999999"
        rows.append(extra)

        with self.assertRaisesRegex(SystemExit, "expected 20 logical jobs, found 21"):
            self.validate(rows)

    def test_array_task_substitution_is_rejected_at_cardinality_gate(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        record = next(row for row in rows if row["JobID"] == "19560001_7")
        record["JobID"] = "19560001_8"

        with self.assertRaisesRegex(SystemExit, "array/cardinality differs"):
            self.validate(rows)

    def test_nonzero_exit_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        next(row for row in rows if row["JobID"] == "19560004_3")[
            "ExitCode"
        ] = "1:0"

        with self.assertRaisesRegex(SystemExit, "did not complete cleanly"):
            self.validate(rows)

    def test_scheduler_submit_time_mismatch_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        rows[0]["Submit"] = "2026-08-06T05:00:00"

        with self.assertRaisesRegex(SystemExit, "ledger/scheduler submit times differ"):
            self.validate(rows)

    def test_duplicate_raw_job_id_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        rows[1]["JobIDRaw"] = rows[0]["JobIDRaw"]

        with self.assertRaisesRegex(SystemExit, "duplicate raw job IDs"):
            self.validate(rows)

    def test_scheduler_identity_mismatch_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        rows[0]["User"] = "another-user"

        with self.assertRaisesRegex(SystemExit, "resources/paths differ"):
            self.validate(rows)

    def test_executed_job_submit_start_end_chronology_is_required(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        rows[0]["Start"] = "2026-08-05T23:00:00"

        with self.assertRaisesRegex(SystemExit, "executed-job chronology differs"):
            self.validate(rows)

    def test_cancelled_analysis_requires_an_unassigned_node_and_finite_end(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        analysis = next(row for row in rows if row["JobID"] == "19560006")
        analysis["NodeList"] = "cpu001"

        with self.assertRaisesRegex(SystemExit, "not cancelled before execution"):
            self.validate(rows)

        rows = copy.deepcopy(self.fixture.accounting_rows)
        analysis = next(row for row in rows if row["JobID"] == "19560006")
        analysis["End"] = "Unknown"
        with self.assertRaisesRegex(SystemExit, "Missing scheduler timestamp"):
            self.validate(rows)

    def test_observed_array_concurrency_above_four_is_rejected(self):
        normalized = self.validate()
        for row in normalized:
            if row["stage"] == "bam_manifest_hash":
                row["Start"] = "2026-08-06T01:00:00"
                row["End"] = "2026-08-06T01:20:00"
        with self.assertRaisesRegex(SystemExit, "exceeded the four-task cap"):
            attestation.observed_array_concurrency(normalized)

    def test_wrong_requested_memory_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        next(row for row in rows if row["JobID"] == "19560004_0")[
            "ReqMem"
        ] = "31G"

        with self.assertRaisesRegex(SystemExit, "resources/paths differ"):
            self.validate(rows)

    def test_wrong_frozen_wrapper_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        record = next(row for row in rows if row["JobID"] == "19560004_0")
        tokens = shlex.split(record["SubmitLine"])
        tokens[-1] = str(self.fixture.root / "unfrozen/recount_array.sbatch")
        record["SubmitLine"] = shlex.join(tokens)

        with self.assertRaisesRegex(SystemExit, "SubmitLine differs from frozen wrapper"):
            self.validate(rows)

    def test_wrong_submitline_dependency_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        record = next(row for row in rows if row["JobID"] == "19560004_0")
        tokens = shlex.split(record["SubmitLine"])
        dependency_index = next(
            index for index, token in enumerate(tokens) if token.startswith("--dependency=")
        )
        tokens[dependency_index] = "--dependency=afterok:99999999"
        record["SubmitLine"] = shlex.join(tokens)

        with self.assertRaisesRegex(SystemExit, "SubmitLine differs from frozen wrapper"):
            self.validate(rows)

    def test_wrong_ledger_dependency_is_rejected_before_accounting(self):
        merge_row = next(
            row
            for row in self.fixture.ledger_rows
            if row["stage"] == "merge_and_validate"
        )
        merge_row["dependency"] = "afterok:99999999"
        self.fixture.write_ledger()

        with self.assertRaisesRegex(SystemExit, "resources/dependency differ"):
            self.fixture.ledger()

    def test_cancelled_analysis_with_start_timestamp_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        analysis = next(row for row in rows if row["JobID"] == "19560006")
        analysis["Start"] = "2026-08-06T03:30:00"

        with self.assertRaisesRegex(SystemExit, "not cancelled before execution"):
            self.validate(rows)

    def test_cancelled_analysis_with_arm_artifact_is_rejected(self):
        (self.fixture.root / "arms/F_locked/deg_results.csv").write_text(
            "gene,logFC\nENSG000001,0.5\n"
        )

        with self.assertRaisesRegex(SystemExit, "nevertheless wrote an arm artifact"):
            self.validate()

    def test_temporal_dependency_violation_is_rejected(self):
        rows = copy.deepcopy(self.fixture.accounting_rows)
        finalize = next(row for row in rows if row["JobID"] == "19560002")
        finalize["Start"] = "2026-08-06T01:10:00"

        with self.assertRaisesRegex(SystemExit, "temporal dependency order differs"):
            self.validate(rows)

    def test_capture_binds_exact_ledger_and_raw_accounting_hashes(self):
        raw = encode_sacct(self.fixture.accounting_rows)
        fake_sacct = Path("/usr/bin/true")

        def fake_run(command, **kwargs):
            self.assertEqual(kwargs["env"]["TZ"], "UTC")
            self.assertEqual(kwargs["env"]["LC_ALL"], "C")
            if command[-1] == "--version":
                return subprocess.CompletedProcess(
                    command,
                    0,
                    stdout="slurm 23.11.10\n",
                    stderr="",
                )
            return subprocess.CompletedProcess(command, 0, stdout=raw, stderr="")

        with mock.patch.object(attestation.shutil, "which", return_value=str(fake_sacct)), mock.patch.object(
            attestation.subprocess, "run", side_effect=fake_run
        ):
            result = attestation.capture(self.fixture.root)

        validation = result["validation"]
        self.assertEqual(result["raw_sacct_psv"], raw)
        self.assertEqual(
            validation["ledger_sha256"],
            attestation.sha256(self.fixture.root / "contract/slurm_submission.tsv"),
        )
        self.assertEqual(
            validation["raw_sacct_sha256"],
            attestation.sha256_bytes(raw.encode()),
        )
        self.assertEqual(validation["logical_jobs"], 20)
        self.assertEqual(validation["completed_jobs"], 19)
        self.assertEqual(validation["cancelled_unstarted_jobs"], 1)
        self.assertEqual(
            validation["observed_array_max_concurrency"],
            {"bam_manifest_hash": 4, "fragment_recount": 4},
        )

    def test_capture_rejects_ledger_mutation_during_scheduler_query(self):
        raw = encode_sacct(self.fixture.accounting_rows)
        fake_sacct = Path("/usr/bin/true")
        calls = 0

        def fake_run(command, **_kwargs):
            nonlocal calls
            calls += 1
            if calls == 1:
                with (self.fixture.root / "contract/slurm_submission.tsv").open("a") as handle:
                    handle.write("mutated\n")
                return subprocess.CompletedProcess(command, 0, stdout=raw, stderr="")
            return subprocess.CompletedProcess(
                command, 0, stdout="slurm 23.11.10\n", stderr=""
            )

        with mock.patch.object(
            attestation.shutil, "which", return_value=str(fake_sacct)
        ), mock.patch.object(attestation.subprocess, "run", side_effect=fake_run):
            with self.assertRaisesRegex(SystemExit, "changed during attestation"):
                attestation.capture(self.fixture.root)


if __name__ == "__main__":
    unittest.main()
