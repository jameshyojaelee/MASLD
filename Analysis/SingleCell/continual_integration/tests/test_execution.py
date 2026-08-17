from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from masld_cl.config import canonical_json_bytes
from masld_cl.contracts import ContractError, sha256_tree
from masld_cl.execution import (
    _load_and_verify_execution_lock, _load_and_verify_historical_execution_lock,
    verify_source_identity_payload, write_execution_lock,
)


class TestExecutionLock(unittest.TestCase):
    def test_historical_source_identity_payload_detects_mutation(self):
        entries = [{"path": "x.py", "sha256": "abc", "size_bytes": 1}]
        identity = {
            "entries": entries,
            "identity_sha256": hashlib.sha256(canonical_json_bytes(entries)).hexdigest(),
        }
        self.assertEqual(verify_source_identity_payload(identity), identity)
        identity["entries"][0]["size_bytes"] = 2
        with self.assertRaises(ContractError):
            verify_source_identity_payload(identity)

    def test_source_and_input_tampering_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            project = Path(directory) / "project"
            pipeline = project / "Analysis" / "SingleCell" / "continual_integration"
            for relative in ("masld_cl", "scripts", "tests", "slurm", "reference"):
                (pipeline / relative).mkdir(parents=True)
            (pipeline / "masld_cl/core.py").write_text("VALUE = 1\n")
            (pipeline / "scripts/train_reference.py").write_text("pass\n")
            (pipeline / "tests/test_x.py").write_text("pass\n")
            (pipeline / "slurm/x.sbatch").write_text("#!/bin/bash\n")
            (pipeline / "reference/x.tsv").write_text("x\n")
            config = pipeline / "config_v1.json"
            config.write_text(json.dumps({"schema_version": "masld-cl-v1"}))
            (pipeline / "README.md").write_text("candidate\n")
            prepared = project / "prepared.h5ad"
            prepared.write_bytes(b"counts")
            output = project / "result"
            spec_path = project / "run.json"
            spec = {
                "schema_version": "masld-cl-execution-spec-v2",
                "resource_class": "large",
                "script": "train_reference.py",
                "arguments": [
                    "--config", str(config), "--prepared", str(prepared),
                    "--output", str(output),
                ],
            }
            spec["spec_sha256"] = hashlib.sha256(canonical_json_bytes(spec)).hexdigest()
            spec_path.write_text(json.dumps(spec))
            lock_path = spec_path.with_suffix(".json.source-lock.json")
            write_execution_lock(spec_path, pipeline, lock_path)
            _load_and_verify_execution_lock(spec, spec_path, pipeline)
            (pipeline / "masld_cl/core.py").write_text("VALUE = 2\n")
            _load_and_verify_historical_execution_lock(spec, spec_path, pipeline)
            with self.assertRaisesRegex(ContractError, "source or input changed"):
                _load_and_verify_execution_lock(spec, spec_path, pipeline)
            prepared.write_bytes(b"changed")
            with self.assertRaisesRegex(ContractError, "historical execution input changed"):
                _load_and_verify_historical_execution_lock(spec, spec_path, pipeline)

    def test_directory_bundle_hash_detects_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "weights.pt").write_bytes(b"one")
            before, _ = sha256_tree(root)
            (root / "weights.pt").write_bytes(b"two")
            after, _ = sha256_tree(root)
            self.assertNotEqual(before, after)


if __name__ == "__main__":
    unittest.main()
