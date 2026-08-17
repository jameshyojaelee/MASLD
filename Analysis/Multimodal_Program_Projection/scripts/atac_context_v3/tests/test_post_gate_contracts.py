#!/usr/bin/env python3
"""Post-gate fixtures for contracts that do not alter sealed producers."""

from __future__ import annotations

import csv
import gzip
import hashlib
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import sys


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

import atac_context_v3_lib as contracts  # noqa: E402
from genetics_context import (  # noqa: E402
    aggregate_lineage_mass,
    normalize_against_reference,
)


class PostGateContractTests(unittest.TestCase):
    def test_candidate_guard_rejects_canonical_traversal_and_symlink_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            fake_project = Path(temp) / "project"
            expected_parent = (
                fake_project / "Analysis/Multimodal_Program_Projection/candidates"
            )
            expected_parent.mkdir(parents=True)
            outside = Path(temp) / "outside"
            outside.mkdir()
            escaped = expected_parent / contracts.RELEASE_ID
            escaped.symlink_to(outside, target_is_directory=True)
            with mock.patch.object(contracts, "project_root", return_value=fake_project):
                with self.assertRaises(contracts.ContractError):
                    contracts.candidate_root(
                        fake_project / "RNA-seq/results" / contracts.RELEASE_ID
                    )
                with self.assertRaises(contracts.ContractError):
                    contracts.candidate_root(expected_parent / ".." / contracts.RELEASE_ID)
                with self.assertRaises(contracts.ContractError):
                    contracts.candidate_root(escaped)

    def test_new_path_guard_rejects_existing_and_dangling_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            existing = root / "existing.tsv"
            existing.write_text("sealed\n", encoding="utf-8")
            dangling = root / "dangling.tsv"
            dangling.symlink_to(root / "absent.tsv")
            with self.assertRaises(contracts.ContractError):
                contracts.require_new_path(existing)
            with self.assertRaises(contracts.ContractError):
                contracts.require_new_path(dangling)

    def test_schema_round_trip_and_hash_are_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "round_trip.tsv"
            rows = [{"b": "two", "a": "one"}, {"b": "four", "a": "three"}]
            contracts.write_tsv(path, ("a", "b"), rows)
            columns, observed = contracts.read_tsv(path, required=("a", "b"))
            self.assertEqual(columns, ["a", "b"])
            self.assertEqual(observed, [
                {"a": "one", "b": "two"},
                {"a": "three", "b": "four"},
            ])
            digest = contracts.sha256_file(path)
            self.assertEqual(digest, contracts.sha256_file(path))

    def test_complete_dynamic_and_genetic_state_vocabularies(self) -> None:
        dynamic = {
            contracts.classify_da(0.01, 0.01, 1.0, 1.0),
            contracts.classify_da(0.01, 0.01, 1.0, -1.0),
            contracts.classify_da(0.01, 0.50, 1.0, 1.0),
            contracts.classify_da(0.50, 0.50, 1.0, 1.0),
            "untestable",
        }
        self.assertEqual(dynamic, contracts.EVIDENCE_STATES)
        genetic = {
            contracts.classify_genetic(0.94, 0.80, 0.80, 0.80),
            contracts.classify_genetic(0.99, 0.60, 0.70, 0.70),
            contracts.classify_genetic(0.99, 0.20, 0.60, 0.20),
            contracts.classify_genetic(0.99, 0.10, 0.20, 0.10),
            contracts.classify_genetic(0.99, 0.00, 0.00, 0.00),
        }
        self.assertEqual(genetic, contracts.GENETIC_STATES)

    def test_unmapped_allele_and_posterior_loss_are_fail_closed(self) -> None:
        reference = {"chr1": "AACCGGTTAA"}

        def fetch(chrom: str, start: int, end: int) -> str:
            return reference[chrom][start:end]

        self.assertIsNone(
            normalize_against_reference("chr1", 3, "A", "T", fetch)
        )
        variants = [
            {"posterior": 0.94, "mapped": True, "chrom": "chr1", "position_1based": 3},
            {"posterior": 0.06, "mapped": False, "chrom": "", "position_1based": 0},
        ]
        result = aggregate_lineage_mass(
            variants,
            lambda chrom, position: True,
            lambda chrom, position: True,
        )
        self.assertEqual(result["state"], "untestable")
        self.assertAlmostEqual(float(result["lost_mass"]), 0.06)

    def test_liftover_wiring_requires_unique_mapping_and_preserves_mass(self) -> None:
        liftover = (SCRIPT_DIR / "08_liftover_variants.R").read_text(encoding="utf-8")
        aggregate = (SCRIPT_DIR / "09_aggregate_genetic_context.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("n_liftover_mappings == 1L", liftover)
        self.assertIn('status = "unmapped_or_multimapped"', aggregate)
        self.assertIn("SNP.PP.H4 does not sum to one", aggregate)
        self.assertIn('"lost_posterior_mass": 1.0 - mapped', aggregate)
        self.assertIn("mapped_posterior_mass_below_0.95", aggregate)
        self.assertNotIn("renormal", liftover.lower())
        self.assertNotIn("renormal", aggregate.lower())

    def test_no_primary_masl_design_label_in_model_qc(self) -> None:
        candidate = contracts.default_candidate_root()
        with (candidate / "da/da_model_qc.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        primary = [
            row for row in rows
            if row["cohort"] == "GSE244832"
            and row["model"] == "primary_condition_only"
        ]
        self.assertTrue(primary)
        self.assertFalse(any(row["condition"] == "MASL" for row in primary))

    def test_promoted_replay_preparation_is_atomic_and_manifest_pinned(self) -> None:
        def digest(path: Path) -> str:
            return hashlib.sha256(path.read_bytes()).hexdigest()

        with tempfile.TemporaryDirectory() as temp:
            fixture = Path(temp)
            candidate = fixture / "candidate"
            candidate.mkdir()
            (candidate / "consensus_peak_manifest.tsv").write_text(
                "release_id\tlineage\nfixture\thepatocyte\n", encoding="utf-8"
            )
            registry = fixture / "registry.tsv"
            tier = fixture / "tier.tsv"
            aggregate = fixture / "aggregate.csv"
            upstream_inputs = fixture / "upstream_inputs.tsv"
            eqtl = fixture / "eqtl/chr1/ENSG_FIXTURE_susie.rds"
            eqtl.parent.mkdir(parents=True)
            eqtl.write_bytes(b"fixture-rds")
            studies = [f"study_{index:02d}" for index in range(35)]
            with registry.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=("study_name", "ancestry"), delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows({"study_name": study, "ancestry": "EUR"} for study in studies)
            with tier.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("study_name", "trait", "tier", "tier_label", "placement"),
                    delimiter="\t", lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows({
                    "study_name": study,
                    "trait": "PDFF" if index == 0 else "ALT",
                    "tier": "1" if index == 0 else "2",
                    "tier_label": "direct_MASLD" if index == 0 else "liver_enzyme",
                    "placement": "main",
                } for index, study in enumerate(studies))
            with aggregate.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("gwas_name", "gene", "ensembl", "chr", "PP.H4.susie", "method"),
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow({
                    "gwas_name": studies[0], "gene": "FIXTURE", "ensembl": "ENSG_FIXTURE",
                    "chr": 1, "PP.H4.susie": 0.75, "method": "susie",
                })
            with upstream_inputs.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=("path", "sha256"), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow({"path": str(tier), "sha256": digest(tier)})
            runner = contracts.project_root() / "GWAS/finemapping/src/06_susie_coloc.R"
            runner_before = digest(runner)
            promoted = fixture / "promoted.tsv"
            fields = (
                "release_id", "status", "aggregate_path", "aggregate_sha256",
                "runner_path", "runner_sha256", "gwas_registry_path",
                "gwas_registry_sha256", "gwas_tier_path", "gwas_tier_sha256",
                "eqtl_susie_dir", "upstream_input_manifest_path",
                "upstream_input_manifest_sha256",
            )
            with promoted.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow({
                    "release_id": "fixture-coloc", "status": "PROMOTED",
                    "aggregate_path": aggregate, "aggregate_sha256": digest(aggregate),
                    "runner_path": runner, "runner_sha256": digest(runner),
                    "gwas_registry_path": registry, "gwas_registry_sha256": digest(registry),
                    "gwas_tier_path": tier, "gwas_tier_sha256": digest(tier),
                    "eqtl_susie_dir": eqtl.parents[1],
                    "upstream_input_manifest_path": upstream_inputs,
                    "upstream_input_manifest_sha256": digest(upstream_inputs),
                })
            subprocess.run(
                [
                    sys.executable, str(SCRIPT_DIR / "07_prepare_genetic_replay.py"),
                    "--candidate-root", str(candidate), "--fixture-mode",
                    "--promoted-manifest", str(promoted),
                ],
                check=True, text=True, capture_output=True,
            )
            prepared = candidate / "genetics/prepared"
            self.assertTrue((prepared / "PROMOTED_UPSTREAM_GATE.tsv").is_file())
            self.assertTrue((prepared / "replay_input_manifest.tsv").is_file())
            self.assertEqual(len(list((prepared / "replay_batches").glob("batch_*.tsv"))), 5)
            with (prepared / "replay_input_manifest.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                replay_inputs = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(replay_inputs)
            self.assertTrue(all(Path(row["path"]).is_file() for row in replay_inputs))
            self.assertTrue(all(digest(Path(row["path"])) == row["sha256"] for row in replay_inputs))
            self.assertEqual(digest(runner), runner_before)
            self.assertFalse(list((candidate / "genetics").glob(".prepared.pending.*")))

    def test_replay_batch_validator_checks_complete_posterior_family(self) -> None:
        def digest(path: Path) -> str:
            return hashlib.sha256(path.read_bytes()).hexdigest()

        def write(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
            path.parent.mkdir(parents=True, exist_ok=True)
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)

        with tempfile.TemporaryDirectory() as temp:
            candidate = Path(temp) / "candidate"
            prepared = candidate / "genetics/prepared"
            frozen = candidate / "frozen.txt"
            frozen.parent.mkdir(parents=True)
            frozen.write_text("frozen\n", encoding="utf-8")
            input_manifest = prepared / "replay_input_manifest.tsv"
            write(
                input_manifest,
                ("release_id", "role", "path", "bytes", "sha256"),
                [{
                    "release_id": contracts.RELEASE_ID, "role": "fixture",
                    "path": frozen, "bytes": frozen.stat().st_size,
                    "sha256": digest(frozen),
                }],
            )
            write(
                prepared / "PROMOTED_UPSTREAM_GATE.tsv",
                ("status", "replay_input_manifest_sha256"),
                [{"status": "PROMOTED", "replay_input_manifest_sha256": digest(input_manifest)}],
            )
            plan_fields = (
                "release_id", "batch_id", "gwas_name", "gene", "ensembl", "chr",
                "trait", "trait_class", "ancestry", "promoted_pp_h4_susie",
            )
            write(
                prepared / "replay_batches/batch_1.tsv", plan_fields,
                [{
                    "release_id": contracts.RELEASE_ID, "batch_id": 1,
                    "gwas_name": "study", "gene": "GENE", "ensembl": "ENSG1",
                    "chr": 1, "trait": "PDFF", "trait_class": "direct_MASLD",
                    "ancestry": "EUR", "promoted_pp_h4_susie": 0.75,
                }],
            )
            pending = candidate / "genetics/.batch_1.pending.fixture"
            export = pending / "exports/study/chr1/ENSG1"
            write(
                export / "signal_pairs.tsv",
                (
                    "gwas_name", "ensembl", "signal_pair_index", "posterior_column",
                    "PP.H4.abf", "gwas_signal", "eqtl_signal",
                ),
                [{
                    "gwas_name": "study", "ensembl": "ENSG1", "signal_pair_index": 1,
                    "posterior_column": "SNP.PP.H4.abf", "PP.H4.abf": 0.75,
                    "gwas_signal": "L1", "eqtl_signal": "L1",
                }],
            )
            write(
                export / "variant_posteriors.tsv.gz",
                (
                    "gwas_name", "ensembl", "signal_pair_index", "posterior_column",
                    "snp", "SNP.PP.H4", "hg19_position", "allele1", "allele2",
                ),
                [{
                    "gwas_name": "study", "ensembl": "ENSG1", "signal_pair_index": 1,
                    "posterior_column": "SNP.PP.H4.abf", "snp": "1:100",
                    "SNP.PP.H4": 1.0, "hg19_position": 100, "allele1": "A", "allele2": "C",
                }],
            )
            subprocess.run(
                [
                    sys.executable, str(SCRIPT_DIR / "22_validate_replay_batch.py"),
                    "--batch", "1", "--pending", str(pending),
                    "--candidate-root", str(candidate), "--fixture-mode",
                ],
                check=True, text=True, capture_output=True,
            )
            manifest = pending / "batch_manifest.tsv"
            self.assertTrue(manifest.is_file())
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 2)


if __name__ == "__main__":
    unittest.main()
