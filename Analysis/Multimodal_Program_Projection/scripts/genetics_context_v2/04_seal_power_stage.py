#!/usr/bin/env python3
"""Seal the base-R power preflight products."""

from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    assert_candidate_root,
    write_stage_seal,
)


def main() -> None:
    candidate = assert_candidate_root(PROJECT_ROOT, DEFAULT_CANDIDATE_ROOT)
    outputs = [
        candidate / name
        for name in (
            "power_stratified_interface.tsv",
            "power_match_balance.tsv",
            "sensitivity.tsv",
            "gate_status.tsv",
        )
    ]
    missing = [str(path) for path in outputs if not path.is_file()]
    if missing:
        raise ContractError("cannot seal missing power products: " + ", ".join(missing))
    upstream = [
        candidate / "work" / "stage_seals" / "03_build_eqtl_observability.json"
    ]
    write_stage_seal(candidate, "04_power_interface", outputs, upstream)
    print("PASS: sealed 04_power_interface")


if __name__ == "__main__":
    main()

