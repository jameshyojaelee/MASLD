from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts.prepare_corgi_regular_film_head_fixture import (
    CONTEXT_ARMS,
    CorgiFilmFixtureError,
    build_fixture,
    stable_tiles,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/corgi_regular_film_plus_head_smoke.toml"
CONTRACT_PREFLIGHT = ROOT / "executions/model-check-260-21096725"
TILE_ROSTER = ROOT / "executions/corgi-tile-smoke-subset-21068454"
MAPPERS = {
    0: ROOT / "executions/corgi-gse296875-mapper-fold0-21066334",
    1: ROOT / "executions/corgi-gse296875-mapper-fold1-21066335",
    2: ROOT / "executions/corgi-gse296875-mapper-fold2-21066337",
    3: ROOT / "executions/corgi-gse296875-mapper-fold3-21066336",
    4: ROOT / "executions/corgi-gse296875-mapper-fold4-21066338",
}


class CorgiFilmHeadFixtureTests(unittest.TestCase):
    def test_stable_tile_sampling_is_order_invariant(self) -> None:
        tiles = [
            {"tile_id": f"tile-{index}", "genomic_fold": "2"}
            for index in range(20)
        ]
        forward = stable_tiles(
            tiles,
            limit=4,
            phase="inner_train",
            evaluation_fold=0,
            epoch=1,
            donor_id="donor-a",
        )
        reverse = stable_tiles(
            list(reversed(tiles)),
            limit=4,
            phase="inner_train",
            evaluation_fold=0,
            epoch=1,
            donor_id="donor-a",
        )
        self.assertEqual(forward, reverse)
        self.assertEqual(len({row["tile_id"] for row in forward}), 4)

    def test_real_fixture_is_outcome_blind_and_crossed_fold_safe(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "fixture"
            receipt = build_fixture(
                contract_path=CONTRACT,
                contract_preflight=CONTRACT_PREFLIGHT,
                mapper_roots=MAPPERS,
                tile_roster=TILE_ROSTER,
                output=output,
            )
            self.assertEqual(
                receipt["status"], "pass_outcome_blind_crossed_fold_fixture"
            )
            self.assertEqual(receipt["contexts"], 39 * len(CONTEXT_ARMS))
            self.assertEqual(receipt["donors"], 39)
            self.assertIs(receipt["development_atac_outcomes_read"], False)
            self.assertIs(receipt["model_checkpoint_loaded"], False)
            self.assertIs(receipt["training_execution_authorized"], False)
            with (output / "schedule.tsv").open(encoding="utf-8", newline="") as handle:
                schedule = list(csv.DictReader(handle, delimiter="\t"))
            for evaluation_fold in range(5):
                rows = [
                    row
                    for row in schedule
                    if int(row["evaluation_fold"]) == evaluation_fold
                ]
                fit = [row for row in rows if row["phase"] != "held_predict"]
                held = [row for row in rows if row["phase"] == "held_predict"]
                self.assertFalse(
                    {row["donor_id"] for row in fit}
                    & {row["donor_id"] for row in held}
                )
                self.assertTrue(
                    all(int(row["genomic_fold"]) != evaluation_fold for row in fit)
                )
                self.assertEqual({row["context_arm"] for row in held}, set(CONTEXT_ARMS))

    def test_rejects_mapper_authority_substitution(self) -> None:
        changed = dict(MAPPERS)
        changed[0] = MAPPERS[1]
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(CorgiFilmFixtureError):
                build_fixture(
                    contract_path=CONTRACT,
                    contract_preflight=CONTRACT_PREFLIGHT,
                    mapper_roots=changed,
                    tile_roster=TILE_ROSTER,
                    output=Path(temporary) / "fixture",
                )

    def test_contract_preflight_remains_execution_blocked(self) -> None:
        receipt = json.loads(
            (CONTRACT_PREFLIGHT / "contract/contract_receipt.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIs(receipt["training_execution_authorized"], False)
        self.assertIs(receipt["prediction_execution_authorized"], False)
        self.assertIs(receipt["evaluation_execution_authorized"], False)


if __name__ == "__main__":
    unittest.main()
