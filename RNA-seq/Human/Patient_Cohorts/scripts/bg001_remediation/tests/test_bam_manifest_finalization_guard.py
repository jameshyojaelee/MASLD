#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))


def load_guard():
    spec = importlib.util.spec_from_file_location(
        "bg001_bam_manifest_finalization_guard",
        SCRIPT_DIR / "bam_manifest_finalization_guard.py",
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


guard = load_guard()


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    path.chmod(0o600)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TransactionFixture:
    def __init__(self, base: Path):
        self.root = base / "bg001-fragment-v211-gencode49-20260806T120000Z"
        self.root.mkdir(mode=0o700)
        self.manifests = self.root / "manifests"
        self.manifests.mkdir(mode=0o700)
        self.shards = self.manifests / "bam_hash_shards"
        self.shards.mkdir(mode=0o700)
        self.bams = self.root / "fixture_bams"
        self.bams.mkdir(mode=0o700)
        for directory in (self.root, self.manifests, self.shards, self.bams):
            directory.chmod(0o700)
        self.rows: list[dict[str, str]] = []
        serial = 0
        for dataset, count in guard.EXPECTED.items():
            for index in range(count):
                sample = f"{dataset}_S{index:03d}"
                bam = self.bams / f"{sample}.bam"
                bam.write_bytes(f"{sample}\n".encode())
                bam.chmod(0o600)
                status = bam.stat()
                self.rows.append(
                    {
                        "dataset": dataset,
                        "sample_id": sample,
                        "bam_path": str(bam),
                        "layout": "paired",
                        "strandedness": "2",
                        "expected_status": "included",
                        "exclusion_reason": "",
                        "bam_size": str(status.st_size),
                        "bam_mtime": str(int(status.st_mtime)),
                        "bam_sha256": "",
                    }
                )
                serial += 1
        for (dataset, sample), reason in guard.DOCUMENTED_UNAVAILABLE.items():
            self.rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample,
                    "bam_path": "",
                    "layout": "paired",
                    "strandedness": "2",
                    "expected_status": "documented_unavailable",
                    "exclusion_reason": reason,
                    "bam_size": "",
                    "bam_mtime": "",
                    "bam_sha256": "",
                }
            )
        write_tsv(self.manifests / "bam_manifest.draft.tsv", guard.MANIFEST_FIELDS, self.rows)
        for task in range(8):
            rows = []
            for index, source in enumerate(guard.expected_task_rows(self.rows, task)):
                rows.append(
                    {
                        "task": str(task),
                        "dataset": source["dataset"],
                        "sample_id": source["sample_id"],
                        "bam_path": source["bam_path"],
                        "bam_size": source["bam_size"],
                        "bam_mtime": source["bam_mtime"],
                        "bam_sha256": hashlib.sha256(
                            f"{task}:{index}:{source['sample_id']}".encode()
                        ).hexdigest(),
                    }
                )
            write_tsv(self.shards / f"task_{task}.tsv", guard.SHARD_FIELDS, rows)

    def publish_unchecked(self) -> None:
        hashes: dict[tuple[str, str], str] = {}
        shard_hashes: dict[str, str] = {}
        for task in range(8):
            path = self.shards / f"task_{task}.tsv"
            shard_hashes[path.name] = guard.sha256(path)
            for row in read_tsv(path):
                hashes[(row["dataset"], row["sample_id"])] = row["bam_sha256"]
        final_rows = [dict(row) for row in self.rows]
        for row in final_rows:
            if row["expected_status"] == "included":
                row["bam_sha256"] = hashes[(row["dataset"], row["sample_id"])]
        manifest = self.manifests / "bam_manifest.tsv"
        write_tsv(manifest, guard.MANIFEST_FIELDS, final_rows)
        record = {
            "schema_version": "1.0",
            "run_id": self.root.name,
            "draft_sha256": guard.sha256(self.manifests / "bam_manifest.draft.tsv"),
            "manifest_sha256": guard.sha256(manifest),
            "included_bams": 820,
            "documented_unavailable": 2,
            "shard_sha256": shard_hashes,
            "finalized_utc": "2026-08-06T12:30:00+00:00",
        }
        freeze = self.manifests / "bam_manifest.freeze.json"
        freeze.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        freeze.chmod(0o600)
        marker = self.root / "BAM_MANIFEST_FROZEN"
        marker.write_text(record["manifest_sha256"] + "\n")
        marker.chmod(0o600)


class FinalizationGuardTests(unittest.TestCase):
    def test_valid_pre_and_postflight_bind_exact_transaction(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            pre = guard.preflight(fixture.root, invoke_verifier=False)
            fixture.publish_unchecked()
            post = guard.postflight(fixture.root, invoke_verifier=False)

            self.assertEqual(pre["status"], "PASS")
            self.assertEqual(post["status"], "PASS")
            self.assertEqual(post["included_bams"], 820)
            self.assertEqual(post["documented_unavailable"], 2)
            self.assertEqual(post["task_cardinalities"], [78, 216, 93, 66, 92, 92, 92, 91])
            self.assertEqual(post["gse213621_chunk_cardinalities"], [92, 92, 92, 91])
            verified = guard.verify_artifacts(fixture.root, invoke_verifier=False)
            self.assertEqual(verified, post)
            for relative in (guard.PRE_ARTIFACT, guard.POST_ARTIFACT):
                status = os.lstat(fixture.root / relative)
                self.assertTrue(stat.S_ISREG(status.st_mode))
                self.assertEqual(stat.S_IMODE(status.st_mode), 0o600)
                self.assertEqual(status.st_nlink, 1)

    def test_verify_rejects_postflight_content_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            guard.preflight(fixture.root, invoke_verifier=False)
            fixture.publish_unchecked()
            guard.postflight(fixture.root, invoke_verifier=False)
            post_path = fixture.root / guard.POST_ARTIFACT
            post = json.loads(post_path.read_text())
            post["task_cardinalities"][0] = 77
            post_path.write_text(json.dumps(post, indent=2, sort_keys=True) + "\n")
            post_path.chmod(0o600)
            with self.assertRaisesRegex(guard.GuardError, "independent reconstruction"):
                guard.verify_artifacts(fixture.root, invoke_verifier=False)

    def test_preflight_rejects_group_writable_candidate_parent(self):
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            fixture = TransactionFixture(parent)
            parent.chmod(0o770)
            try:
                with self.assertRaisesRegex(guard.GuardError, "candidate parent"):
                    guard.preflight(fixture.root, invoke_verifier=False)
            finally:
                parent.chmod(0o700)

    def test_preflight_rejects_predictable_temp_symlink_without_touching_target(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            outside = Path(temporary) / "outside.tsv"
            outside.write_text("must remain unchanged\n")
            temporary_path = fixture.manifests / "bam_manifest.tsv.tmp"
            temporary_path.symlink_to(outside)

            with self.assertRaisesRegex(guard.GuardError, "predictable.*temporaries"):
                guard.preflight(fixture.root, invoke_verifier=False)
            self.assertEqual(outside.read_text(), "must remain unchanged\n")
            self.assertFalse((fixture.manifests / "bam_manifest.tsv").exists())

    def test_preflight_rejects_non_owner_only_manifests_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            fixture.manifests.chmod(0o770)
            with self.assertRaisesRegex(guard.GuardError, "owner-only"):
                guard.preflight(fixture.root, invoke_verifier=False)

    def test_preflight_rejects_existing_final_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            (fixture.root / "BAM_MANIFEST_FROZEN").write_text("premature\n")
            with self.assertRaisesRegex(guard.GuardError, "final outputs"):
                guard.preflight(fixture.root, invoke_verifier=False)

    def test_postflight_rejects_cross_chunk_attribution_even_when_outputs_agree(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            chunk0_path = fixture.shards / "task_4.tsv"
            chunk1_path = fixture.shards / "task_5.tsv"
            chunk0 = read_tsv(chunk0_path)
            chunk1 = read_tsv(chunk1_path)
            chunk0[0], chunk1[0] = chunk1[0], chunk0[0]
            chunk0[0]["task"] = "4"
            chunk1[0]["task"] = "5"
            write_tsv(chunk0_path, guard.SHARD_FIELDS, chunk0)
            write_tsv(chunk1_path, guard.SHARD_FIELDS, chunk1)
            guard.preflight(fixture.root, invoke_verifier=False)
            fixture.publish_unchecked()

            with self.assertRaisesRegex(guard.GuardError, "provenance differs"):
                guard.postflight(fixture.root, invoke_verifier=False)

    def test_postflight_rejects_manifest_byte_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            guard.preflight(fixture.root, invoke_verifier=False)
            fixture.publish_unchecked()
            manifest = fixture.manifests / "bam_manifest.tsv"
            manifest.write_bytes(manifest.read_bytes() + b"\n")
            manifest.chmod(0o600)

            with self.assertRaisesRegex(guard.GuardError, "byte-identical"):
                guard.postflight(fixture.root, invoke_verifier=False)

    def test_postflight_rejects_marker_hash_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            guard.preflight(fixture.root, invoke_verifier=False)
            fixture.publish_unchecked()
            marker = fixture.root / "BAM_MANIFEST_FROZEN"
            marker.write_text("0" * 64 + "\n")
            marker.chmod(0o600)

            with self.assertRaisesRegex(guard.GuardError, "exact final-manifest"):
                guard.postflight(fixture.root, invoke_verifier=False)

    def test_postflight_rejects_freeze_json_hash_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            guard.preflight(fixture.root, invoke_verifier=False)
            fixture.publish_unchecked()
            freeze = fixture.manifests / "bam_manifest.freeze.json"
            record = json.loads(freeze.read_text())
            record["shard_sha256"]["task_7.tsv"] = "0" * 64
            freeze.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            freeze.chmod(0o600)

            with self.assertRaisesRegex(guard.GuardError, "hashes or counts"):
                guard.postflight(fixture.root, invoke_verifier=False)

    def test_postflight_rejects_hardlinked_final_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            guard.preflight(fixture.root, invoke_verifier=False)
            fixture.publish_unchecked()
            manifest = fixture.manifests / "bam_manifest.tsv"
            os.link(manifest, fixture.manifests / "manifest_second_link.tsv")

            with self.assertRaisesRegex(guard.GuardError, "exactly one hard link"):
                guard.postflight(fixture.root, invoke_verifier=False)

    def test_postflight_rejects_input_mutation_after_preflight(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            guard.preflight(fixture.root, invoke_verifier=False)
            shard = fixture.shards / "task_0.tsv"
            shard.write_bytes(shard.read_bytes() + b"\n")
            shard.chmod(0o600)
            fixture.publish_unchecked()

            with self.assertRaisesRegex(guard.GuardError, "changed after preflight"):
                guard.postflight(fixture.root, invoke_verifier=False)

    def test_frozen_verifier_is_called_with_required_frozen_flag(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = TransactionFixture(Path(temporary))
            verifier = fixture.root / (
                "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/"
                "bg001_remediation/verify_run_contract.py"
            )
            verifier.parent.mkdir(parents=True)
            verifier.write_text(
                "import sys\n"
                "assert '--run-root' in sys.argv\n"
                "assert '--require-bam-frozen' in sys.argv\n"
                "print('fixture verifier PASS')\n"
            )
            verifier.chmod(0o600)
            record = guard.run_contract_verifier(fixture.root, require_frozen=True)
            self.assertEqual(record["stdout"], "fixture verifier PASS")
            self.assertIn("--require-bam-frozen", record["command"])


if __name__ == "__main__":
    unittest.main()
