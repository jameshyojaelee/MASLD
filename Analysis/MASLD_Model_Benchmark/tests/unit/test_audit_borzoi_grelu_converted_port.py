from __future__ import annotations

from collections import OrderedDict
from copy import deepcopy
import json
from pathlib import Path
import pickle
import unittest

from scripts.audit_borzoi_grelu_converted_port import (
    ConvertedPortAuditError,
    scan_pickle,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]


class BorzoiGreluConvertedPortAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(
            (ROOT / "config/borzoi_grelu_converted_port_acquisition.json").read_text()
        )

    def test_contract_is_fail_closed_and_records_format_incident(self) -> None:
        validate_contract(self.contract)
        self.assertFalse(self.contract["model_identity"]["native_borzoi_substitute"])
        self.assertEqual(
            self.contract["source_format_incident"]["observed_file_magic_hex"],
            "504b0304",
        )
        self.assertFalse(
            self.contract["admission_disposition"]["open_champion_eligible"]
        )

    def test_contract_rejects_weights_only_load_authorization(self) -> None:
        changed = deepcopy(self.contract)
        changed["admission_disposition"][
            "torch_weights_only_load_allowed_by_this_contract"
        ] = True
        with self.assertRaisesRegex(ConvertedPortAuditError, "contract opened"):
            validate_contract(changed)

    def test_static_pickle_scan_does_not_unpickle(self) -> None:
        payload = pickle.dumps(OrderedDict(), protocol=2)
        receipt = scan_pickle(
            payload,
            allowed_globals={"collections OrderedDict"},
            forbidden_opcodes={"STACK_GLOBAL", "INST", "OBJ"},
        )
        self.assertFalse(receipt["pickle_payload_unpickled"])
        self.assertEqual(receipt["pickle_globals"], ["collections OrderedDict"])

    def test_static_pickle_scan_rejects_unexpected_global(self) -> None:
        payload = b"\x80\x02cos\nsystem\n."
        with self.assertRaisesRegex(ConvertedPortAuditError, "unexpected"):
            scan_pickle(
                payload,
                allowed_globals={"collections OrderedDict"},
                forbidden_opcodes={"STACK_GLOBAL", "INST", "OBJ"},
            )


if __name__ == "__main__":
    unittest.main()
