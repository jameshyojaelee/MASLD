from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest

from scripts.run_lsgkm_gse281364_fit import LSGKMFitError, validate_fit_request


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/lsgkm_gse281364_dinucleotide_full_campaign.json"


def authorization(config: dict) -> dict:
    return {
        "status": "authorized_remaining_24_outcome_blind_fits",
        "authorized_run_fits": [
            {"split_id": split_id, "model_seed": seed}
            for split_id, split in config["fit_grid"]["splits"].items()
            for seed in split["run_seeds"]
        ],
        "action_firewall": copy.deepcopy(config["action_firewall"]),
    }


class LSGKMFullFitTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))
        self.authorization = authorization(self.config)

    def test_remaining_fit_is_admitted(self) -> None:
        validate_fit_request(
            self.config, self.authorization, "donor4_genomic4", 8111
        )

    def test_reused_probe_fit_is_rejected(self) -> None:
        with self.assertRaises(LSGKMFitError):
            validate_fit_request(
                self.config, self.authorization, "donor0_genomic0", 1103
            )

    def test_unknown_seed_is_rejected(self) -> None:
        with self.assertRaises(LSGKMFitError):
            validate_fit_request(
                self.config, self.authorization, "donor1_genomic1", 9999
            )

    def test_expanded_authority_is_rejected(self) -> None:
        changed = copy.deepcopy(self.authorization)
        changed["authorized_run_fits"].append(
            {"split_id": "donor0_genomic0", "model_seed": 1103}
        )
        with self.assertRaises(LSGKMFitError):
            validate_fit_request(
                self.config, changed, "donor1_genomic1", 1103
            )

    def test_metric_or_outcome_permission_is_rejected(self) -> None:
        for field in (
            "outcome_access_authorized",
            "benchmark_metric_calculation_authorized",
            "control_prediction_value_access_authorized",
        ):
            with self.subTest(field=field):
                changed = copy.deepcopy(self.authorization)
                changed["action_firewall"][field] = True
                with self.assertRaises(LSGKMFitError):
                    validate_fit_request(
                        self.config, changed, "donor1_genomic1", 1103
                    )


if __name__ == "__main__":
    unittest.main()
