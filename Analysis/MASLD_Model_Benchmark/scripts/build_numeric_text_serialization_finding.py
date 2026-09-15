#!/usr/bin/env python3
"""Freeze the reusable finding that R's default text serialization loses bits.

``utils::write.table`` converts doubles at 15 significant digits.  An IEEE-754
double needs 17.  On a log2 expression scale roughly 91 percent of values lose
their low bits, silently, while the pipeline reports success.

This matters more here than it usually would.  The GSE49541 and GSE83452 lanes
rest on bitwise reproducibility: the separation fixtures certify that an array's
summary is bit-identical across contexts, and a lossy writer would discard those
exact bits on the way to disk, so a downstream re-read would not reproduce the
digests the separation certified.

The record also carries the methodological trap that the first verification
attempt fell into, because it is the more transferable lesson: the only
preserved values were *post*-truncation, so re-writing them at 15 digits
reproduced the same decimals and every writer appeared lossless.  Test data
drawn from the output of the defect cannot detect the defect.  A negative
control that must fail is what makes a pass mean anything.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


class SerializationFindingError(RuntimeError):
    """Raised when the evidence does not support the finding."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _scalar(value: Any) -> Any:
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def build_finding(
    *, roundtrip_receipt: Path, summarization_verdict: Path, output: Path
) -> dict[str, Any]:
    evidence = json.loads(roundtrip_receipt.read_text(encoding="utf-8"))
    if _scalar(evidence.get("status")) != "pass_repaired_writer_round_trips_and_default_does_not":
        raise SerializationFindingError("round-trip evidence is inconclusive")
    if _scalar(evidence.get("test_values_perturbed_to_restore_mantissa_entropy")) is not True:
        raise SerializationFindingError(
            "test values were not perturbed; the check would be circular"
        )

    default = {key: _scalar(value) for key, value in evidence["default_write_table"].items()}
    fifteen = {
        key: _scalar(value) for key, value in evidence["fifteen_significant_digits"].items()
    }
    seventeen = {
        key: _scalar(value) for key, value in evidence["seventeen_significant_digits"].items()
    }
    # The negative controls must fail, or a pass by the repaired writer is empty.
    if default["bitwise_identical"] is not False or fifteen["bitwise_identical"] is not False:
        raise SerializationFindingError("negative control did not fail; the test is vacuous")
    if seventeen["bitwise_identical"] is not True:
        raise SerializationFindingError("the repaired writer does not round-trip")

    tested = int(_scalar(evidence["values_tested"]))
    differing = int(default["differing_elements"])

    verdict = json.loads(summarization_verdict.read_text(encoding="utf-8"))
    if verdict.get("status") != "pass_label_blind_single_array_summarization":
        raise SerializationFindingError("the summarization verdict is not a pass")

    finding = {
        "schema_version": "masld-bench-numeric-text-serialization-finding-v1",
        "status": "resolved_by_writing_17_significant_digits_and_asserting_the_round_trip",
        "finding_id": "r_write_table_lossy_for_doubles_20260825",
        "finding": (
            "utils::write.table serialises doubles at 15 significant digits; an "
            "IEEE-754 double needs 17. Written matrices were a lossy "
            "representation of what the summarizer computed."
        ),
        "detected_by": "a bitwise round-trip assertion after every numeric matrix write",
        "would_not_have_been_detected_by": (
            "any tolerance-based comparison, any summary statistic, or the job "
            "exit status; the pipeline reported success while writing rounded values"
        ),
        "measured_on": "real GSE49541 fRMA values at genuine log2 magnitudes",
        "values_tested": tested,
        "values_differing_under_default_writer": differing,
        "fraction_differing_under_default_writer": differing / tested,
        "maximum_absolute_error_log2": float(default["maximum_absolute_difference"]),
        "default_writer_bitwise_identical": False,
        "fifteen_digit_writer_bitwise_identical": False,
        "seventeen_digit_writer_bitwise_identical": True,
        "resolution": {
            "action": 'format explicitly with sprintf("%.17g", x) before write.table',
            "guard_retained": True,
            "guard": (
                "read the file back and require identical() on dim, dimnames and "
                "the full numeric vector; never relax to all.equal or a tolerance"
            ),
            "applied_to": [
                "gse49541_gpl570_platform_feature_matrix.tsv",
                "gse49541_gpl570_gene_matrix.tsv",
            ],
        },
        "methodological_trap": {
            "what_happened": (
                "the first verification reused the preserved matrix directly and "
                "reported every writer as lossless"
            ),
            "why": (
                "those values had already been written at 15 significant digits, so "
                "each was the nearest double to a 15-digit decimal; re-writing them "
                "at 15 digits reproduced the same decimal and the same double"
            ),
            "lesson": (
                "test data drawn from the output of a defect cannot detect that "
                "defect; a negative control that must fail is what makes a pass "
                "meaningful"
            ),
            "values_already_clean_at_15_significant_digits": int(
                _scalar(evidence["values_already_clean_at_15_significant_digits"])
            ),
        },
        "generalization": {
            "applies_to": "every numeric matrix this campaign writes as text",
            "named_exposures": [
                "GSE83452 GPL16686 summarization of 231 arrays",
                "any downstream gene matrix, score table or prediction file",
            ],
            "why_it_matters_here": (
                "the firewall fixtures certify bitwise reproducibility of an "
                "array's summary; a lossy writer discards exactly those bits, so a "
                "re-read would not reproduce the certified digests"
            ),
        },
        "bound_evidence": {
            "roundtrip_receipt": {
                "path": str(roundtrip_receipt),
                "sha256": sha256_file(roundtrip_receipt),
            },
            "summarization_verdict": {
                "path": str(summarization_verdict),
                "sha256": sha256_file(summarization_verdict),
            },
        },
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }

    output.mkdir(parents=True, exist_ok=False)
    (output / "numeric_text_serialization_finding.json").write_text(
        json.dumps(finding, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return finding


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--roundtrip-receipt", type=Path, required=True)
    parser.add_argument("--summarization-verdict", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    finding = build_finding(
        roundtrip_receipt=arguments.roundtrip_receipt,
        summarization_verdict=arguments.summarization_verdict,
        output=arguments.output,
    )
    print(json.dumps({k: finding[k] for k in ("status", "finding_id", "finding")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
