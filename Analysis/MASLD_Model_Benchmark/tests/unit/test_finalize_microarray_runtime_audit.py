from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.finalize_microarray_runtime_audit import (
    BASELINE_RUNTIME_IDS,
    MicroarrayRuntimeAuditError,
    RUNTIME_FIXED_KEYS,
    read_runtime_receipt,
    read_single_array_capability,
    require_runtime_inventory,
)

PASSING_CAPABILITY_PROBE = {
    "status": "pass_single_array_runtime_capability",
    "arrays_read_per_platform": 1,
    "capabilities": {
        "calvin_header_parse_both_platforms": True,
        "GPL570_single_array_summarization": True,
        "GPL16686_single_array_summarization": True,
        "probe_to_entrez_to_ensembl_both_platforms": True,
    },
    "all_sample_RMA_run": False,
    "across_array_quantile_normalization_run": False,
    "GEO_series_matrix_read": False,
    "expression_values_exported": False,
    "labels_read": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
}


class MicroarrayRuntimeReceiptTests(unittest.TestCase):
    def write_receipt(self, path: Path, values: list[str]) -> None:
        rows = ["key\tvalue"]
        rows.extend(
            f"{key}\t{value}" for key, value in zip(RUNTIME_FIXED_KEYS, values, strict=True)
        )
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    def test_complete_fixed_schema_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            receipt = Path(temporary) / "runtime.tsv"
            self.write_receipt(receipt, ["explicit_state"] * len(RUNTIME_FIXED_KEYS))
            observed = read_runtime_receipt(receipt)
        self.assertEqual(tuple(observed), RUNTIME_FIXED_KEYS)
        self.assertEqual(len(observed), 55)

    def test_blank_missing_state_is_rejected(self) -> None:
        values = ["explicit_state"] * len(RUNTIME_FIXED_KEYS)
        values[RUNTIME_FIXED_KEYS.index("reader.GPL570.error")] = ""
        with tempfile.TemporaryDirectory() as temporary:
            receipt = Path(temporary) / "runtime.tsv"
            self.write_receipt(receipt, values)
            with self.assertRaisesRegex(MicroarrayRuntimeAuditError, "explicit"):
                read_runtime_receipt(receipt)

    def test_dropped_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            receipt = Path(temporary) / "runtime.tsv"
            rows = ["key\tvalue"]
            rows.extend(f"{key}\texplicit_state" for key in RUNTIME_FIXED_KEYS[:-1])
            receipt.write_text("\n".join(rows) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(MicroarrayRuntimeAuditError, "schema"):
                read_runtime_receipt(receipt)


class MicroarrayRuntimeInventoryTests(unittest.TestCase):
    def test_activated_runtime_joins_the_baseline_inventory(self) -> None:
        require_runtime_inventory(set(BASELINE_RUNTIME_IDS) | {"microarray_r_bioc_portable"})

    def test_dropped_baseline_runtime_is_rejected(self) -> None:
        reduced = set(BASELINE_RUNTIME_IDS) - {"module_R_4.4.3"}
        with self.assertRaisesRegex(MicroarrayRuntimeAuditError, "incomplete"):
            require_runtime_inventory(reduced)


class MicroarraySingleArrayCapabilityTests(unittest.TestCase):
    def write_probe(self, directory: str, probe: dict) -> Path:
        path = Path(directory) / "capability_probe.json"
        path.write_text(json.dumps(probe), encoding="utf-8")
        return path

    def test_absent_probe_reports_not_provided(self) -> None:
        self.assertEqual(read_single_array_capability(None), "not_provided")

    def test_passing_probe_is_summarized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = self.write_probe(temporary, PASSING_CAPABILITY_PROBE)
            observed = read_single_array_capability(path)
        self.assertIsInstance(observed, dict)
        self.assertEqual(observed["arrays_read_per_platform"], 1)
        self.assertEqual(len(observed["sha256"]), 64)

    def test_all_sample_rma_probe_is_rejected(self) -> None:
        probe = dict(PASSING_CAPABILITY_PROBE, all_sample_RMA_run=True)
        with tempfile.TemporaryDirectory() as temporary:
            path = self.write_probe(temporary, probe)
            with self.assertRaisesRegex(MicroarrayRuntimeAuditError, "capability probe"):
                read_single_array_capability(path)

    def test_unsummarized_platform_is_rejected(self) -> None:
        capabilities = dict(
            PASSING_CAPABILITY_PROBE["capabilities"],
            GPL16686_single_array_summarization=False,
        )
        probe = dict(PASSING_CAPABILITY_PROBE, capabilities=capabilities)
        with tempfile.TemporaryDirectory() as temporary:
            path = self.write_probe(temporary, probe)
            with self.assertRaisesRegex(MicroarrayRuntimeAuditError, "capability probe"):
                read_single_array_capability(path)


if __name__ == "__main__":
    unittest.main()
