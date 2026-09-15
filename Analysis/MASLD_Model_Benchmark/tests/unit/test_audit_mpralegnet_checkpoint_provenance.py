from __future__ import annotations

from collections import defaultdict
import os
import pickle
import pickletools
import socket
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

from safetensors.torch import load_file, save_file
import torch

from scripts.audit_mpralegnet_checkpoint_provenance import (
    ProvenanceAuditError,
    _extract_tensor_state,
    _state_dict_from_author_checkpoint,
    build_evaluator_contract,
    compare_checkpoint_states,
    inspect_pickle_program,
    load_config,
    probe_official_unsafe_global_scanner,
    validate_inner_member_names,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/gse281364_mpralegnet_provenance_audit.toml"


class MPRALegNetProvenanceAuditTests(unittest.TestCase):
    def test_config_freezes_outcome_and_seed_firewalls(self) -> None:
        config = load_config(CONFIG)
        validate_config(config)
        self.assertFalse(config["outcome_access_authorized"])
        self.assertFalse(config["metric_calculation_authorized"])
        self.assertEqual(
            config["rekeyed_prediction"]["independent_fitted_predictions"], 1
        )
        self.assertFalse(
            config["rekeyed_prediction"]["schema_seed_repeats_are_independent"]
        )

    def test_tensor_equivalence_is_exact_not_approximate(self) -> None:
        author = {
            "model.weight": torch.tensor([[1.0, 2.0]], dtype=torch.float32),
            "model.count": torch.tensor(3, dtype=torch.int64),
        }
        equal, rows, summary = compare_checkpoint_states(
            author, {key: value.clone() for key, value in author.items()}
        )
        self.assertTrue(equal)
        self.assertTrue(summary["all_tensors_bitwise_equal"])
        self.assertTrue(all(row["bitwise_equal"] == "true" for row in rows))
        changed = {key: value.clone() for key, value in author.items()}
        changed["model.weight"][0, 0] += torch.finfo(torch.float32).eps
        equal, _, summary = compare_checkpoint_states(author, changed)
        self.assertFalse(equal)
        self.assertFalse(summary["all_tensors_bitwise_equal"])

    def test_blocked_scanner_never_calls_torch_load(self) -> None:
        with (
            patch.object(
                torch.serialization,
                "get_unsafe_globals_in_checkpoint",
                side_effect=pickle.UnpicklingError("Unsupported operand 48"),
            ),
            patch.object(torch, "load") as load,
        ):
            receipt = probe_official_unsafe_global_scanner(b"checkpoint")
        self.assertEqual(
            receipt["status"],
            "blocked_official_weights_only_scanner_unsupported_POP_opcode",
        )
        self.assertFalse(receipt["torch_load_invoked"])
        self.assertFalse(receipt["weights_only_load_attempted"])
        self.assertFalse(receipt["weights_only_false_used"])
        self.assertFalse(receipt["custom_deserializer_used"])
        self.assertFalse(receipt["unrestricted_deserializer_used"])
        load.assert_not_called()

    def test_undeclared_pickle_global_is_rejected(self) -> None:
        payload = pickle.dumps(defaultdict(list), protocol=2)
        opcodes = sorted(
            {operation.name for operation, _, _ in pickletools.genops(payload)}
        )
        with self.assertRaisesRegex(ProvenanceAuditError, "global set differs"):
            inspect_pickle_program(
                payload, declared_globals=(), allowed_opcodes=opcodes
            )

    def test_undeclared_pickle_opcode_is_rejected(self) -> None:
        payload = pickle.dumps({"value": 1}, protocol=2)
        with self.assertRaisesRegex(ProvenanceAuditError, "opcode set differs"):
            inspect_pickle_program(
                payload, declared_globals=(), allowed_opcodes=("PROTO", "STOP")
            )

    def test_malicious_code_reducer_is_rejected_without_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / "executed"

            class Malicious:
                def __reduce__(self):
                    return os.system, (f"touch {marker}",)

            payload = pickle.dumps(Malicious(), protocol=2)
            opcodes = sorted(
                {operation.name for operation, _, _ in pickletools.genops(payload)}
            )
            with self.assertRaisesRegex(ProvenanceAuditError, "global set differs"):
                inspect_pickle_program(
                    payload, declared_globals=(), allowed_opcodes=opcodes
                )
            self.assertFalse(marker.exists())

    def test_malicious_network_reducer_is_rejected_without_network(self) -> None:
        class Malicious:
            def __reduce__(self):
                return socket.create_connection, (("127.0.0.1", 9),)

        payload = pickle.dumps(Malicious(), protocol=2)
        opcodes = sorted(
            {operation.name for operation, _, _ in pickletools.genops(payload)}
        )
        with self.assertRaisesRegex(ProvenanceAuditError, "global set differs"):
            inspect_pickle_program(
                payload, declared_globals=(), allowed_opcodes=opcodes
            )

    def test_non_tensor_state_dict_leaf_is_rejected(self) -> None:
        with self.assertRaisesRegex(ProvenanceAuditError, "non-tensor leaves"):
            _extract_tensor_state(
                {"state_dict": {"model.weight": torch.ones(1), "epoch": 24}}
            )

    def test_unexpected_inner_archive_member_is_rejected(self) -> None:
        names = {
            "archive/data.pkl",
            "archive/version",
            "archive/data/0",
            "archive/evil.py",
        }
        with self.assertRaisesRegex(ProvenanceAuditError, "member names differ"):
            validate_inner_member_names(names, 1)

    def test_tensor_schema_drift_is_rejected(self) -> None:
        reference = {"model.weight": torch.ones((2, 2), dtype=torch.float32)}
        shape_drift = {"model.weight": torch.ones((4,), dtype=torch.float32)}
        key_drift = {"model.bias": torch.ones((2, 2), dtype=torch.float32)}
        self.assertFalse(compare_checkpoint_states(reference, shape_drift)[0])
        self.assertFalse(compare_checkpoint_states(reference, key_drift)[0])

    def test_safetensors_round_trip_matches_comparator(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.safetensors"
            state = {"model.weight": torch.arange(4, dtype=torch.float32)}
            save_file(state, path)
            equal, _, _ = compare_checkpoint_states(state, load_file(path))
            self.assertTrue(equal)

    def test_evaluator_contract_forbids_execution_and_deduplicates_seeds(self) -> None:
        contract = build_evaluator_contract(load_config(CONFIG))
        self.assertFalse(contract["execution"]["authorized_in_this_audit"])
        self.assertFalse(contract["execution"]["outcomes_read"])
        self.assertEqual(contract["prediction_input"]["expected_rows"], 1033)
        self.assertEqual(
            contract["evaluation"]["schema_seed_rows"],
            "deduplicated_not_independent_replications",
        )


if __name__ == "__main__":
    unittest.main()
