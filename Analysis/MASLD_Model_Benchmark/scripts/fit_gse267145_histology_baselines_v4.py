#!/usr/bin/env python3
"""Launch the validated v3 fitter from the production child environment."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if ROOT.as_posix() not in sys.path:
    sys.path.insert(0, ROOT.as_posix())

import scripts.fit_gse267145_histology_baselines_v3 as v3


SOFTWARE_REVISION_ID = "model-training-069-production-v4"
SOFTWARE_CONTRACT_PATH = (
    ROOT / "config/evaluation/gse267145_histology_production_v4.json"
)
SOFTWARE_CONTRACT_SHA256 = (
    "df1f77f82c3657e6786be9a67100b3ce138384f8b28eb6a22b48bc3e0cf97de2"
)


def _verify_software_contract() -> None:
    if (
        not SOFTWARE_CONTRACT_PATH.is_file()
        or SOFTWARE_CONTRACT_PATH.is_symlink()
        or v3.v2.sha256_file(SOFTWARE_CONTRACT_PATH) != SOFTWARE_CONTRACT_SHA256
    ):
        raise v3.v2.HistologyBaselineError("v4 software contract SHA-256 differs")
    contract = json.loads(SOFTWARE_CONTRACT_PATH.read_text(encoding="utf-8"))
    correction = contract.get("software_correction", {})
    if (
        contract.get("revision_id") != SOFTWARE_REVISION_ID
        or correction.get("v3_fitter_sha256")
        != "0cf28366d86040ce18fec73e3a5cde862ab33dad71cf0719bdd1d15a14508d3e"
        or correction.get("modeling_code_changed") is not False
        or contract.get("scoring_authorized") is not False
    ):
        raise v3.v2.HistologyBaselineError("v4 software contract differs")


def main() -> int:
    _verify_software_contract()
    if "--software-self-test" in sys.argv:
        if sys.argv != [sys.argv[0], "--software-self-test"]:
            raise v3.v2.HistologyBaselineError("v4 self-test arguments differ")
        v3._verify_contract()
        print(
            json.dumps(
                {
                    "status": "passed_child_import_contract",
                    "software_revision_id": SOFTWARE_REVISION_ID,
                    "repository_root_on_sys_path": ROOT.as_posix() in sys.path,
                    "v3_fallback_contract_sha256": v3.CONTRACT_SHA256,
                    "outcomes_read": False,
                    "metrics_calculated": False,
                    "scorer_called": False,
                },
                sort_keys=True,
            )
        )
        return 0
    return v3.main()


if __name__ == "__main__":
    raise SystemExit(main())
