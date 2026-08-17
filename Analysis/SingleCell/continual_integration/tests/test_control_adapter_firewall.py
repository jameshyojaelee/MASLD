from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from masld_cl.config import canonical_json_bytes
from masld_cl.firewall import FirewallError, load_control_adapter_selection_lock


class TestControlAdapterFirewall(unittest.TestCase):
    def _lock(self, root: Path):
        config = {"_config_sha256": "config", "lineages": ["L1"]}
        policy = root / "policy.json"
        policy.write_text("{}")
        results = []
        for weight in (0.25, 0.5, 0.75):
            for kind in ("all_lineage", "L1"):
                manifest = root / f"manifest-{weight}-{kind}.json"
                manifest.write_text(json.dumps({
                    "adapter_source_identity": {"identity_sha256": "source"}
                }))
                result = root / f"result-{weight}-{kind}.json"
                result.write_text(json.dumps({
                    "schema_version": "masld-cl-control-adapter-pilot-v7",
                    "config_sha256": "config", "control_only": True,
                    "outcomes_unlocked": False, "global_reference_weight": weight,
                    "model_kind": kind, "sources": {
                        "adapter_manifest": str(manifest),
                        "adapter_manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
                    },
                }))
                results.append({"path": str(result), "sha256": hashlib.sha256(result.read_bytes()).hexdigest()})
        lock = {
            "schema_version": "masld-cl-control-adapter-selection-v7",
            "config_sha256": "config", "selection_frozen": True,
            "outcomes_unlocked": True, "selected_global_reference_weight": 0.75,
            "policy_realpath": str(policy),
            "policy_sha256": hashlib.sha256(policy.read_bytes()).hexdigest(),
            "control_results": results,
            "source_identity": {"identity_sha256": "source"},
        }
        lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
        path = root / "selection.json"
        path.write_text(json.dumps(lock))
        return config, path

    def test_accepts_exact_frozen_grid_and_fails_on_mutation(self):
        with tempfile.TemporaryDirectory() as value:
            config, path = self._lock(Path(value))
            loaded = load_control_adapter_selection_lock(path, config)
            self.assertEqual(loaded["selected_global_reference_weight"], 0.75)
            source = Path(loaded["control_results"][0]["path"])
            source.write_text("{}")
            with self.assertRaises(FirewallError):
                load_control_adapter_selection_lock(path, config)
