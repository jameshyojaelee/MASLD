from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.build_saf_partition_substrate import (
    SAF_ACTIVITY_COMPONENTS,
    build_arm_a_saf_triple,
)
from scripts.build_showcase_aspect_lineage_substrate import SubstrateError

HEADER = [
    "participant_id", "outer_fold", "steatosis", "ballooning",
    "lobular_inflammation", "nas_sum_deposited", "lobular_necrosis", "fibrosis",
]


def write_arm_a(root: Path, rows: list[list[str]]) -> None:
    path = root / "arm_a" / "aspect_endpoints.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(["\t".join(HEADER)] + ["\t".join(r) for r in rows]) + "\n",
        encoding="utf-8",
    )


class SafActivityDerivationTests(unittest.TestCase):
    def test_activity_is_ballooning_plus_lobular_inflammation(self) -> None:
        self.assertEqual(
            SAF_ACTIVITY_COMPONENTS, ("ballooning", "lobular_inflammation")
        )
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            write_arm_a(root, [["P1", "0", "2", "1", "2", "5", "0", "1"]])
            out = root / "out"
            receipt = build_arm_a_saf_triple(root, out)
            body = (out / "arm_a_saf_triple.tsv").read_text().splitlines()
            columns = body[0].split("\t")
            values = dict(zip(columns, body[1].split("\t")))
            self.assertEqual(values["saf_activity_A"], "3")
            self.assertEqual(values["saf_steatosis_S"], "2")
            self.assertEqual(values["saf_fibrosis_F"], "1")
            self.assertEqual(receipt["saf_triple_state"], "derived")

    def test_lobular_necrosis_never_enters_the_activity_grade(self) -> None:
        """Necrosis is deposited beside the components and is not one of them."""

        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            # necrosis 2 would push A to 5 if it were wrongly included
            write_arm_a(root, [["P1", "0", "1", "1", "1", "3", "2", "0"]])
            out = root / "out"
            build_arm_a_saf_triple(root, out)
            body = (out / "arm_a_saf_triple.tsv").read_text().splitlines()
            values = dict(zip(body[0].split("\t"), body[1].split("\t")))
            self.assertEqual(values["saf_activity_A"], "2")
            self.assertEqual(values["lobular_necrosis"], "2")

    def test_a_reconstruction_that_misses_the_nas_sum_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            # S + A = 4 but the deposited sum says 5
            write_arm_a(root, [["P1", "0", "2", "1", "1", "5", "0", "1"]])
            with self.assertRaises(SubstrateError):
                build_arm_a_saf_triple(root, root / "out")

    def test_activity_realises_its_full_zero_to_four_range(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = Path(value)
            rows = []
            for index, (b, i) in enumerate(
                [(0, 0), (0, 1), (1, 1), (2, 1), (2, 2)]
            ):
                rows.append(
                    [f"P{index}", "0", "1", str(b), str(i), str(1 + b + i), "0", "0"]
                )
            out = root / "out"
            write_arm_a(root, rows)
            receipt = build_arm_a_saf_triple(root, out)
            self.assertEqual(
                receipt["census"]["activity_A"], {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
            )
            self.assertEqual(receipt["levels_realised"]["activity_A"], 5)


if __name__ == "__main__":
    unittest.main()
