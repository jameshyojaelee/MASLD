from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from scripts.activate_alphagenome_gse281364_task_native import (
    AlphaGenomeActivationError,
    AUTOSOMES,
    compare_primary_assemblies,
    credential_presence,
    parse_assembly_report,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/alphagenome_gse281364_task_native_activation.json"


def assembly_report(*, changed_chr22: bool = False) -> str:
    rows = ["# synthetic NCBI assembly report"]
    for index, contig in enumerate(AUTOSOMES, start=1):
        accession_index = index + (100 if changed_chr22 and contig == "chr22" else 0)
        rows.append(
            "\t".join(
                (
                    str(index),
                    "assembled-molecule",
                    str(index),
                    "Chromosome",
                    f"CM{accession_index:06d}.2",
                    "=",
                    f"NC_{accession_index:06d}.11",
                    "Primary Assembly",
                    str(1000000 + index),
                    contig,
                )
            )
        )
    return "\n".join(rows) + "\n"


class AlphaGenomeGse281364TaskNativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_activation_keeps_outcomes_and_model_execution_closed(self) -> None:
        disposition = self.config["activation_disposition"]
        self.assertFalse(disposition["api_metadata_call_allowed"])
        self.assertFalse(disposition["api_prediction_call_allowed"])
        self.assertFalse(disposition["checkpoint_download_allowed"])
        self.assertFalse(disposition["outcome_access_allowed"])
        self.assertFalse(disposition["sealed_access_allowed"])
        self.assertFalse(disposition["open_champion_eligible"])

    def test_native_plan_prespecifies_cell_matched_and_fitted_lanes(self) -> None:
        task = self.config["task_contract"]
        self.assertEqual(task["contexts"], ["HepG2_control", "HepG2_PAOA"])
        self.assertTrue(task["sequence_prediction_reused_identically_across_contexts"])
        self.assertEqual(task["input_length_bp"], 524288)
        self.assertEqual(task["variant_index0"], 262144)
        self.assertEqual(task["native_scorers"][0]["output_type"], "DNASE")
        self.assertEqual(
            [row["comparator_id"] for row in task["prespecified_comparators"]],
            [
                "alphagenome_hepg2_dnase_zero_shot",
                "alphagenome_all_dnase_lasso",
                "alphagenome_all_multimodal_lasso",
            ],
        )

    def test_assembly_report_requires_all_autosomes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "report.txt"
            path.write_text(assembly_report(), encoding="utf-8")
            self.assertEqual(set(parse_assembly_report(path)), set(AUTOSOMES))
            path.write_text("\n".join(assembly_report().splitlines()[:-1]) + "\n", encoding="utf-8")
            with self.assertRaises(AlphaGenomeActivationError):
                parse_assembly_report(path)

    def test_primary_assembly_comparison_rejects_accession_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            p13 = root / "p13.txt"
            p14 = root / "p14.txt"
            p13.write_text(assembly_report(), encoding="utf-8")
            p14.write_text(assembly_report(changed_chr22=True), encoding="utf-8")
            with self.assertRaises(AlphaGenomeActivationError):
                compare_primary_assemblies(p13, p14, root / "crosswalk.tsv")

    def test_credential_audit_returns_booleans_without_secret_material(self) -> None:
        policy = self.config["credential_policy"]
        with mock.patch.dict(
            os.environ,
            {"ALPHA_GENOME_API_KEY": "do-not-record", "HF_TOKEN": "also-do-not-record"},
            clear=True,
        ):
            observed = credential_presence(self.config)
        self.assertEqual(
            observed,
            {"api_credential_present": True, "huggingface_credential_present": True},
        )
        self.assertNotIn("do-not-record", json.dumps(observed))
        self.assertTrue(policy["values_hashes_lengths_and_prefixes_forbidden"])


if __name__ == "__main__":
    unittest.main()
