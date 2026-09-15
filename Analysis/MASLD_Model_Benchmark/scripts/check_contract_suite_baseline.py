#!/usr/bin/env python3
"""Check a campaign on the full contract suite without excluding any test.

The registry contract suite currently carries pre-existing failures owned by
other in-flight model-registry campaigns.  Dropping those tests from a check
would leave the check with a hole that outlives the defect, so this runs the
whole suite and instead requires two things:

* every failure is one of the known pre-existing ones, so this campaign has
  introduced no new failure, and
* the tests this campaign is actually responsible for ran and passed.

A pre-existing failure that gets fixed elsewhere simply disappears from the
observed set, which still passes.  A new failure anywhere fails the check.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys


# Failures observed on the unmodified tree before this campaign changed
# anything. All are model-registry bundles and census bindings; none is a
# dataset, split, pairing, or schema requirements.
# Two entries were removed on 2026-08-25 as the debt was paid:
# test_every_campaign_binds_the_frozen_model_census (census-binding scope was a
# test defect, not a campaign defect) and
# test_uce_upstream_preflight_and_exposure_are_frozen (the assertion was stale;
# the weight bodies are on disk and re-derive byte-exact). The register shrinks
# test_scgpt_upstream_preflight_and_exposure_are_frozen was removed on
# 2026-08-25 once the bundle was corrected to match the bodies on disk; see
# config/artifacts/incidents/scgpt_weight_rights_disposition_20260825.json.
# test_model_census_maps_exactly_to_audit_bundles was removed on 2026-08-25:
# the two orphaned directories are cross-cutting disposition records, not census
# members, and registering them was provably self-defeating (104 bundles / 2
# orphans becomes 89 / 17). The audit root is now scoped on species via
# non_census_record_directories, guarded by
# tests/unit/test_non_census_record_directories.py.
# test_every_model_audit_bundle_is_complete_and_hash_bound was removed on
# 2026-08-25: its original completeness cause was resolved by the census
# scoping, and the scimilarity -v2 schemas it then failed on are sanctioned
# evolution, not drift. All three v2 files are strict supersets of the shape
# shared by every v1 peer (checked against all 62 / 103 / 61), and schema
# versioning past v1 is established practice here. The accepted sets were
# widened; test_versioned_audit_schemas_only_ever_add_fields stops a bump ever
# being used to SHED a field.
# test_scimilarity_preflight_and_exposure_are_frozen was removed on 2026-08-25:
# two assertions added that same day required admission_blocking False and
# missing pieces [], while config/models/cell_foundation.toml was rewritten the same
# day to keep the model blocked. The registry won because
# admission_blocking == bool(missing pieces) holds for 136 of 136 models with zero
# violations. Guarded by
# tests/unit/test_scimilarity_admission_and_schema_correction.py.
# The register shrinks as debt is paid; it started at 8.
PREEXISTING_FAILURES = frozenset(
    {
    }
)

_OUTCOME = re.compile(r"^(FAIL|ERROR): (\w+)")
_RAN = re.compile(r"^Ran (\d+) tests?")
_VERBOSE = re.compile(r"^\S+ \(([\w.]+)\)")


class ContractBaselineError(RuntimeError):
    """Raised when the suite regressed or the required tests did not run."""


def parse_failures(text: str) -> set[str]:
    return {m.group(2) for line in text.splitlines() if (m := _OUTCOME.match(line))}


def parse_ran(text: str) -> int:
    for line in text.splitlines():
        if match := _RAN.match(line):
            return int(match.group(1))
    raise ContractBaselineError("no unittest summary line found")


def parse_passed(text: str) -> set[str]:
    """Names of tests that reported ok, from unittest -v output.

    A verbose line reads ``method (dotted.path.Class.method) ... ok``.  The
    dotted path inside the parentheses is what carries the module name, so it
    is what gets recorded; the bare leading method name would not.
    """

    passed: set[str] = set()
    for line in text.splitlines():
        if " ... ok" not in line:
            continue
        if match := _VERBOSE.match(line):
            passed.add(match.group(1))
        else:
            passed.add(line.split(" ", 1)[0])
    return passed


def check(text: str, required_modules: list[str], minimum_tests: int) -> dict:
    ran = parse_ran(text)
    if ran < minimum_tests:
        raise ContractBaselineError(
            f"suite ran only {ran} tests, expected at least {minimum_tests}; "
            "a collection error may have hidden tests"
        )
    observed = parse_failures(text)
    regressions = sorted(observed - PREEXISTING_FAILURES)
    if regressions:
        raise ContractBaselineError(
            "new contract failures introduced by this campaign: "
            + ", ".join(regressions)
        )
    passed = parse_passed(text)
    for module in required_modules:
        hits = [name for name in passed if module in name]
        if not hits:
            raise ContractBaselineError(
                f"no passing test found for required module {module}"
            )
    return {
        "tests_run": ran,
        "preexisting_failures_still_present": sorted(
            observed & PREEXISTING_FAILURES
        ),
        "preexisting_failures_since_fixed": sorted(
            PREEXISTING_FAILURES - observed
        ),
        "new_failures": [],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unittest-output", type=Path, required=True)
    parser.add_argument("--require-module", action="append", default=[])
    parser.add_argument("--minimum-tests", type=int, default=70)
    args = parser.parse_args()
    text = args.unittest_output.read_text(encoding="utf-8", errors="replace")
    summary = check(text, args.require_module, args.minimum_tests)
    print("contract suite baseline check passed")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ContractBaselineError as error:
        print(f"contract suite baseline check FAILED: {error}", file=sys.stderr)
        raise SystemExit(1)
