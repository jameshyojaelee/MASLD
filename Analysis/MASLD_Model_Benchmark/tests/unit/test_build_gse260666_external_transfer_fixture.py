from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import tomllib
import unittest

import numpy as np

from scripts.build_gse260666_external_transfer_fixture import (
    ExternalTransferFixtureError,
    build_fixture,
    sha256_file,
)


class GSE260666ExternalFixtureTests(unittest.TestCase):
    def _write_tsv(self, path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def _root(self, path: Path) -> str:
        (path / "ARTIFACTS.json").write_text("{}\n", encoding="utf-8")
        return sha256_file(path / "ARTIFACTS.json")

    def _fixture(self, root: Path, *, overlap: bool = False) -> dict[str, object]:
        benchmark = root / "Analysis" / "MASLD_Model_Benchmark"
        activation = benchmark / "activation"
        matrix = benchmark / "matrix"
        training = benchmark / "training"
        training_join = benchmark / "training_join"
        for path in (activation, matrix, training, training_join):
            path.mkdir(parents=True)
        groups = ["healthy_control"] * 6 + ["non-alcoholic_fatty_liver_disease_(nafld)"] * 6 + ["non-alcoholic_steatohepatitis_(nash)"] * 4
        participants = []
        for index, group in enumerate(groups):
            participants.append(
                {
                    "participant_id": f"E{index:02d}",
                    "sample_accession": f"GSME{index:02d}",
                    "biosample_accession": f"SAMNE{index:02d}",
                    "sra_experiment": f"SRXE{index:02d}",
                    "group": group,
                    "fibrosis": "not_reported",
                    "sex": "not_reported",
                    "biological_unit": "participant",
                }
            )
        self._write_tsv(
            activation / "participant_join.tsv",
            tuple(participants[0]),
            participants,
        )
        genes = [f"GENE{index:03d}" for index in range(30)]
        crosswalk = [
            {
                "source_feature_id": gene,
                "mapping_state": "unique_symbol",
                "gencode_v49_stable_id": f"ENSG{index:011d}",
            }
            for index, gene in enumerate(genes)
        ]
        self._write_tsv(activation / "gene_crosswalk.tsv", tuple(crosswalk[0]), crosswalk)
        rng = np.random.default_rng(91)
        external = rng.poisson(20, size=(16, len(genes))).astype(int)
        raw = matrix / "raw" / "GSE260nnn" / "GSE260666_raw_counts.txt.gz"
        raw.parent.mkdir(parents=True)
        with gzip.open(raw, "wt", encoding="utf-8") as handle:
            handle.write("gene\t" + "\t".join(row["participant_id"] for row in participants) + "\n")
            for gene_index, gene in enumerate(genes):
                handle.write(gene + "\t" + "\t".join(str(value) for value in external[:, gene_index]) + "\n")
        training_participants = [
            {
                "participant_id": f"T{index:03d}",
                "rna_source_sample_accession": f"GSMT{index:03d}",
                "h3k27ac_source_sample_accession": f"GSMH{index:03d}",
            }
            for index in range(99)
        ]
        self._write_tsv(
            training / "participant_axis.tsv", tuple(training_participants[0]), training_participants
        )
        training_features = [
            {"stable_gene_id": f"ENSG{index:011d}"} for index in range(len(genes))
        ]
        self._write_tsv(
            training / "rna_feature_axis.tsv", ("stable_gene_id",), training_features
        )
        np.save(training / "rna_values.npy", rng.gamma(2.0, 3.0, size=(99, len(genes))))
        soft_alias = "GSE260666" if overlap else "GSE267145"
        (training_join / "raw").mkdir()
        with gzip.open(training_join / "raw" / "GSE267145_family.soft.gz", "wt", encoding="utf-8") as handle:
            handle.write(f"^SERIES = {soft_alias}\n!Series_relation = BioProject: PRJNA999999\n")
            for index in range(99):
                handle.write(f"^SAMPLE = GSMT{index:03d}\n!Sample_relation = SRA: SRXT{index:03d}\n")
        task = benchmark / "task.toml"
        task.write_text('task_id = "gse260666_external_histology_transfer"\n', encoding="utf-8")
        salt = benchmark / "private-salt.bin"
        salt.write_bytes(bytes(range(32)))
        salt.chmod(0o600)
        return {
            "activation": activation,
            "matrix": matrix,
            "training_molecular": training,
            "training_join": training_join,
            "task_spec": task,
            "project_root": root,
            "benchmark_root": benchmark,
            "fingerprint_salt_file": salt,
            "activation_artifacts_sha256": self._root(activation),
            "matrix_artifacts_sha256": self._root(matrix),
            "training_molecular_artifacts_sha256": self._root(training),
            "training_join_artifacts_sha256": self._root(training_join),
            "task_spec_sha256": sha256_file(task),
        }

    def test_outcome_separation_mapping_and_contamination(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            arguments = self._fixture(root)
            output = root / "output"
            result = build_fixture(**arguments, output=output)
            self.assertEqual(result["status"], "passed")
            model_text = (output / "model_input" / "participant_axis.tsv").read_text()
            model_text += (output / "model_input" / "receipt.json").read_text()
            self.assertNotIn("healthy_control", model_text)
            self.assertNotIn("NAFL", model_text)
            labels = (output / "evaluator_only" / "labels.tsv").read_text()
            self.assertIn("healthy_control", labels)
            self.assertIn("\tNOR\n", labels)
            counts = np.load(output / "model_input" / "rna_counts.npy", allow_pickle=False)
            self.assertEqual(counts.shape, (16, 30))
            self.assertEqual(counts.dtype, np.uint64)
            audit = json.loads((output / "contamination_audit" / "audit.json").read_text())
            self.assertEqual(audit["status"], "passed")
            self.assertEqual(audit["prior_project_use"]["matched_file_count"], 0)
            self.assertTrue(audit["salted_expression_fingerprints_only"])
            self.assertFalse(audit["champion_claim_eligible"])

    def test_source_alias_overlap_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            arguments = self._fixture(root, overlap=True)
            with self.assertRaises(ExternalTransferFixtureError):
                build_fixture(**arguments, output=root / "output")

    def test_task_spec_freezes_one_to_one_mapping_and_claim_boundary(self) -> None:
        task = Path(__file__).parents[2] / "config" / "evaluation" / "gse260666_external_histology_transfer_task.toml"
        with task.open("rb") as handle:
            config = tomllib.load(handle)
        self.assertEqual(config["source_to_evaluation_mapping"]["healthy_control"], "NOR")
        self.assertEqual(config["source_to_evaluation_mapping"]["non-alcoholic_fatty_liver_disease_(nafld)"], "NAFL")
        self.assertEqual(config["source_to_evaluation_mapping"]["non-alcoholic_steatohepatitis_(nash)"], "NASH")
        self.assertEqual(config["datasets_sealed"], [])
        self.assertIn("not project-sealed", config["claim_gate"])


if __name__ == "__main__":
    unittest.main()
