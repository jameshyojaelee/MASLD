#!/usr/bin/env python3
"""Synthetic fail-closed tests for GSE105127 revision-v2 gzip admission."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file
from scripts import preflight_gse105127_gzip_inputs_v2 as preflight
from scripts import quantify_gse105127_rna_bundle_v2 as quantifier_v2


class GSE105127GzipPreflightV2Tests(unittest.TestCase):
    def admit(self, path: Path) -> dict[str, object]:
        payload = path.read_bytes()
        return preflight.strict_single_member_gzip(
            path,
            expected_sha256=hashlib.sha256(payload).hexdigest(),
            expected_bytes=len(payload),
            expected_md5=hashlib.md5(payload, usedforsecurity=False).hexdigest(),
        )

    def test_exactly_one_complete_member_passes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "valid.tsv.gz"
            source = b"header\nvalue\n"
            path.write_bytes(gzip.compress(source, mtime=0))
            observed = self.admit(path)
        self.assertEqual(observed["gzip_member_count"], 1)
        self.assertEqual(observed["trailing_bytes"], 0)
        self.assertEqual(observed["uncompressed_sha256"], hashlib.sha256(source).hexdigest())
        self.assertTrue(observed["gzip_test_passed"])

    def test_truncated_member_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "truncated.tsv.gz"
            path.write_bytes(gzip.compress(b"payload\n", mtime=0)[:-5])
            with self.assertRaisesRegex(
                preflight.GSE105127GzipPreflightError, "truncated|trailer"
            ):
                self.admit(path)

    def test_trailing_garbage_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "garbage.tsv.gz"
            path.write_bytes(gzip.compress(b"payload\n", mtime=0) + b"not-a-member")
            with self.assertRaisesRegex(
                preflight.GSE105127GzipPreflightError, "trailing_garbage"
            ):
                self.admit(path)

    def test_second_member_fails_wrong_member_count(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "two-members.tsv.gz"
            path.write_bytes(
                gzip.compress(b"first\n", mtime=0)
                + gzip.compress(b"second\n", mtime=0)
            )
            with self.assertRaisesRegex(
                preflight.GSE105127GzipPreflightError, "additional_gzip_member"
            ):
                self.admit(path)

    def test_wrong_compressed_identity_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "wrong-sha.tsv.gz"
            path.write_bytes(gzip.compress(b"payload\n", mtime=0))
            with self.assertRaisesRegex(
                preflight.GSE105127GzipPreflightError, "SHA-256 differs"
            ):
                preflight.strict_single_member_gzip(path, expected_sha256="0" * 64)

    def test_invalid_expected_digest_is_rejected_before_admission(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "valid.tsv.gz"
            path.write_bytes(gzip.compress(b"payload\n", mtime=0))
            with self.assertRaisesRegex(
                preflight.GSE105127GzipPreflightError, "expected compressed SHA-256"
            ):
                preflight.strict_single_member_gzip(path, expected_sha256="not-a-digest")

    def test_every_consumer_wrapper_runs_same_job_preflight_first(self) -> None:
        project_root = Path(__file__).resolve().parents[2]
        contracts = {
            "slurm/gse105127_v2_rsem_reference_cpu.sbatch": (
                "run_preflight rsem_reference",
                "scripts/build_gse105127_rsem_reference.py",
            ),
            "slurm/gse105127_v2_rrbs_collapse_io.sbatch": (
                "run_preflight rrbs_collapse",
                "scripts/collapse_gse105127_rrbs_bundle.py",
            ),
            "slurm/gse105127_v2_rna_quant_cpu.sbatch": (
                "run_preflight rna_quantification",
                "scripts/quantify_gse105127_rna_bundle_v2.py",
            ),
            "slurm/gse105127_v2_crosswalk_io.sbatch": (
                "run_preflight cpg_crosswalk",
                "scripts/build_gse105127_cpg_crosswalk.py",
            ),
        }
        for relative, (preflight_call, consumer_call) in contracts.items():
            source = (project_root / relative).read_text(encoding="utf-8")
            self.assertLess(source.index(preflight_call), source.index(consumer_call))
            self.assertIn("freeze_stage_link", source)

    def test_preflight_and_linker_surface_excludes_outcome_assets(self) -> None:
        project_root = Path(__file__).resolve().parents[2]
        for relative in (
            "scripts/preflight_gse105127_gzip_inputs_v2.py",
            "scripts/freeze_gse105127_gzip_stage_link_v2.py",
            "scripts/prepare_gse105127_gzip_preflight_revision_v2.py",
            "slurm/gse105127_gzip_v2_common.sh",
        ):
            source = (project_root / relative).read_text(encoding="utf-8").lower()
            self.assertNotIn("evaluator_only", source)
            self.assertNotIn("outcomes.tsv", source)
            self.assertNotIn("predictions.tsv", source)

    def test_quantifier_v2_retains_raw_path_through_receipt_hash(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw_member = root / "raw"
            reference = root / "reference"
            stage = root / "stage"
            raw_member.mkdir()
            reference.mkdir()
            stage.mkdir()
            (raw_member / "fixture.txt").write_text(
                "synthetic raw identity\n", encoding="utf-8"
            )
            (reference / "fixture.txt").write_text(
                "synthetic reference identity\n", encoding="utf-8"
            )
            freeze_tree(raw_member, {"fixture": "raw"})
            freeze_tree(reference, {"fixture": "reference"})
            (stage / "rsem.genes.results").write_text(
                "gene_id\ttranscript_id(s)\tlength\teffective_length\t"
                "expected_count\tTPM\tFPKM\n"
                "ENSG_FIXTURE\tENST_FIXTURE\t1000\t900\t1.5\t1000000\t1\n",
                encoding="utf-8",
            )
            row = {
                "row_id": "synthetic-row",
                "participant_group_id": "synthetic-participant",
                "zone": "CV",
            }
            receipt = quantifier_v2.finalize_rsem_stage(
                stage=stage,
                raw_member=raw_member,
                reference_root=reference,
                row=row,
                seed=105127,
            )
            verify_frozen_tree(stage)
            with gzip.open(
                stage / "rsem.genes.results.gz", "rt", encoding="utf-8"
            ) as handle:
                compressed_text = handle.read()
            expected_raw_sha = sha256_file(raw_member / "ARTIFACTS.json")
        self.assertEqual(
            receipt["schema_version"], "masld-bench-gse105127-rna-rsem-v2"
        )
        self.assertEqual(receipt["raw_artifacts_sha256"], expected_raw_sha)
        self.assertTrue(receipt["gene_results_gzip_test_passed"])
        self.assertIn("ENSG_FIXTURE", compressed_text)


if __name__ == "__main__":
    unittest.main()
