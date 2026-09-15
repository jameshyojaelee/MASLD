from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.build_numeric_text_serialization_finding import (
    SerializationFindingError,
    build_finding,
)


def roundtrip_receipt(
    *, default_ok: bool = False, fifteen_ok: bool = False, seventeen_ok: bool = True,
    perturbed: bool = True, status: str = "pass_repaired_writer_round_trips_and_default_does_not",
) -> dict[str, object]:
    return {
        "status": status,
        "test_values_perturbed_to_restore_mantissa_entropy": perturbed,
        "values_tested": 3_936_600,
        "values_already_clean_at_15_significant_digits": 336_481,
        "default_write_table": {
            "bitwise_identical": default_ok,
            "differing_elements": 3_600_119,
            "maximum_absolute_difference": 4.974e-14,
        },
        "fifteen_significant_digits": {
            "bitwise_identical": fifteen_ok,
            "differing_elements": 3_600_119,
            "maximum_absolute_difference": 4.974e-14,
        },
        "seventeen_significant_digits": {
            "bitwise_identical": seventeen_ok,
            "differing_elements": 0,
            "maximum_absolute_difference": 0.0,
        },
    }


def summarization_verdict(status: str = "pass_label_blind_single_array_summarization") -> dict[str, object]:
    return {"status": status}


class BuildFindingTests(unittest.TestCase):
    def _write(self, base: Path, receipt: dict[str, object], verdict: dict[str, object]) -> dict[str, Path]:
        r = base / "roundtrip.json"
        v = base / "verdict.json"
        r.write_text(json.dumps(receipt), encoding="utf-8")
        v.write_text(json.dumps(verdict), encoding="utf-8")
        return {"roundtrip_receipt": r, "summarization_verdict": v}

    def test_clean_evidence_produces_the_finding(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(base, roundtrip_receipt(), summarization_verdict())
            finding = build_finding(output=base / "finding", **paths)
            self.assertEqual(
                finding["status"],
                "resolved_by_writing_17_significant_digits_and_asserting_the_round_trip",
            )
            self.assertAlmostEqual(
                finding["fraction_differing_under_default_writer"], 3_600_119 / 3_936_600
            )
            self.assertIs(finding["seventeen_digit_writer_bitwise_identical"], True)
            self.assertEqual(len(finding["bound_evidence"]), 2)

    def test_vacuous_negative_control_is_rejected(self) -> None:
        """If the default writer 'passes', the test proved nothing."""
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(
                base, roundtrip_receipt(default_ok=True), summarization_verdict()
            )
            with self.assertRaisesRegex(SerializationFindingError, "vacuous"):
                build_finding(output=base / "finding", **paths)

    def test_unperturbed_test_data_is_rejected_as_circular(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(base, roundtrip_receipt(perturbed=False), summarization_verdict())
            with self.assertRaisesRegex(SerializationFindingError, "circular"):
                build_finding(output=base / "finding", **paths)

    def test_failing_repaired_writer_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(
                base, roundtrip_receipt(seventeen_ok=False), summarization_verdict()
            )
            with self.assertRaisesRegex(SerializationFindingError, "does not round-trip"):
                build_finding(output=base / "finding", **paths)

    def test_inconclusive_receipt_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(
                base,
                roundtrip_receipt(status="fail_round_trip_evidence_is_inconclusive"),
                summarization_verdict(),
            )
            with self.assertRaisesRegex(SerializationFindingError, "inconclusive"):
                build_finding(output=base / "finding", **paths)

    def test_failing_summarization_verdict_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            paths = self._write(base, roundtrip_receipt(), summarization_verdict("fail"))
            with self.assertRaisesRegex(SerializationFindingError, "not a pass"):
                build_finding(output=base / "finding", **paths)


if __name__ == "__main__":
    unittest.main()
