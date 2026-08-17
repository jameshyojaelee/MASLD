import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from masld_cl.contracts import ContractError, sha256_path
from masld_cl.external_reference_evaluation import (
    _external26_gate_decision,
    _held_external_scores,
    _load_owned_reference,
    _require_identical_common_universe,
    lock_external26_reference_decision,
)


class TestExternalReferenceEvaluation(unittest.TestCase):
    def test_all_locked_gates_select_expanded_reference(self):
        gates, decision = _external26_gate_decision(
            {"ci_low": -0.019}, 0.05, -0.02,
            {"Hepatocytes": -0.05, "T cells": 0.01},
            {"improvement": 0.01, "ci_low": 0.001},
        )
        self.assertTrue(all(item["pass"] for item in gates))
        self.assertEqual(decision, "select_external_clean26")

    def test_any_failed_gate_retains_strict_reference(self):
        gates, decision = _external26_gate_decision(
            {"ci_low": -0.02}, 0.01, 0.0,
            {"Hepatocytes": 0.0, "T cells": 0.0},
            {"improvement": 0.01, "ci_low": 0.001},
        )
        self.assertFalse(gates[0]["pass"])
        self.assertEqual(decision, "retain_common_strict7")

    def test_control_ci_must_exclude_no_improvement(self):
        gates, decision = _external26_gate_decision(
            {"ci_low": 0.0}, 0.01, 0.0,
            {"Hepatocytes": 0.0, "T cells": 0.0},
            {"improvement": 0.01, "ci_low": 0.0},
        )
        self.assertFalse(gates[-1]["pass"])
        self.assertEqual(decision, "retain_common_strict7")

    def test_held_external_scores_require_exact_19_donors(self):
        cells = pd.DataFrame({
            "cell_id": [f"HLiCA|cell{i}" for i in range(19)],
            "donor_id": [f"HLiCA|donor{i}" for i in range(19)],
            "audit_cell_type": ["Hepatocytes"] * 19,
            "predicted_cell_type": ["Hepatocytes"] * 19,
        })
        result = _held_external_scores(cells, ["Hepatocytes"])
        self.assertEqual(result["n_donors"], 19)
        self.assertEqual(result["mean_donor_macro_f1"], 1.0)
        with self.assertRaisesRegex(ContractError, "requires 19 donors"):
            _held_external_scores(cells.iloc[:-1], ["Hepatocytes"])

    def test_owned_reference_rejects_wrong_roster(self):
        embedding = {"schema_version": "masld-cl-embedding-v1"}
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = directory / "embedding_manifest.json"
            path.write_text("{}")
            (directory / "reference_manifest.json").write_text(json.dumps({
                "schema_version": "masld-cl-external26-common-reference-v21",
                "config_sha256": "config",
                "roster_name": "wrong",
                "embedding": embedding,
                "query_labels_hidden": True,
                "case_stage_program_hero_gene_umap_cas13_used": False,
            }))
            bundle = (embedding, np.zeros((1, 2)), pd.DataFrame())
            with patch(
                "masld_cl.external_reference_evaluation.load_embedding",
                return_value=bundle,
            ):
                with self.assertRaisesRegex(ContractError, "provenance differs"):
                    _load_owned_reference(
                        {"_config_sha256": "config"}, path, "common_strict7"
                    )

    def test_common_universe_lock_must_be_identical(self):
        common = {"common_prepared_lock": {"lock_sha256": "same"}}
        expanded = {"common_prepared_lock": {"lock_sha256": "same"}}
        self.assertEqual(_require_identical_common_universe(common, expanded), "same")
        expanded["common_prepared_lock"]["lock_sha256"] = "different"
        with self.assertRaisesRegex(ContractError, "common-universe lock"):
            _require_identical_common_universe(common, expanded)

    def test_reference_decision_lock_selects_strict7_after_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            common_embedding = directory / "common_embedding.json"
            common_manifest = directory / "common_manifest.json"
            expanded_embedding = directory / "expanded_embedding.json"
            expanded_manifest = directory / "expanded_manifest.json"
            for path in (
                common_embedding, common_manifest,
                expanded_embedding, expanded_manifest,
            ):
                path.write_text("{}")
            source = lambda path: {  # noqa: E731
                "path": str(path.resolve()), "sha256": sha256_path(path)
            }
            decision = {
                "schema_version": "masld-cl-external26-reference-decision-v21",
                "config_sha256": "config",
                "control_only": True,
                "case_stage_program_hero_gene_umap_cas13_read": False,
                "decision": "retain_common_strict7",
                "pilot_pass": False,
                "identical_common_universe": True,
                "common_universe_lock_sha256": "common-lock",
                "gates": [{"gate": "retention", "pass": False}],
                "sources": {
                    "common_strict7_embedding": source(common_embedding),
                    "common_strict7_manifest": source(common_manifest),
                    "external_clean26_embedding": source(expanded_embedding),
                    "external_clean26_manifest": source(expanded_manifest),
                },
            }
            decision_path = directory / "decision.json"
            decision_path.write_text(json.dumps(decision))
            lock = lock_external26_reference_decision(
                {"_config_sha256": "config"}, decision_path,
                directory / "lock.json",
            )
            self.assertEqual(lock["selected_roster"], "common_strict7")
            self.assertFalse(lock["expanded_reference_allowed_downstream"])
            self.assertEqual(lock["failed_gates"], ["retention"])


if __name__ == "__main__":
    unittest.main()
