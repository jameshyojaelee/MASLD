from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from masld_bench.artifacts import verify_frozen_tree


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "build_gse296875_rna_atac_smoke.py"
)
R_EXPORTER = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "export_gse296875_atac_smoke.R"
)
SPEC = importlib.util.spec_from_file_location("gse296875_smoke_builder", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load GSE296875 smoke builder")
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


class GSE296875SelectionTests(unittest.TestCase):
    def test_peak_coordinates_accept_exact_scientific_notation_only(self) -> None:
        self.assertEqual(
            builder._parse_exact_integer("2.6e+07", label="peak start"),
            26_000_000,
        )
        for value in ("2.5", "NaN", "Inf", "not-a-coordinate"):
            with self.subTest(value=value):
                with self.assertRaises(builder.GSE296875SmokeError):
                    builder._parse_exact_integer(value, label="peak start")

    def test_python310_builder_freeze_is_control_plane_verifiable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "artifact"
            root.mkdir()
            (root / "payload.txt").write_text("same nucleus\n", encoding="utf-8")
            manifest_sha256 = builder._freeze_smoke_tree(
                root,
                {"artifact_class": "multimodal_smoke_subset"},
            )
            verified = verify_frozen_tree(root)
            self.assertEqual(
                builder.sha256_file(root / "ARTIFACTS.json"), manifest_sha256
            )
            self.assertEqual(len(verified["artifacts"]), 1)

    def test_r_exporter_does_not_require_signac_class_introspection(self) -> None:
        source = R_EXPORTER.read_text(encoding="utf-8")
        self.assertIn("peak_storage <- attributes(peak_assay)", source)
        self.assertNotIn("slotNames(peak_assay)", source)
        self.assertNotIn("library(Signac)", source)

    @staticmethod
    def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
        with path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_selection_is_deterministic_donor_balanced_and_atac_blind(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mapping_path = root / "mapping.tsv"
            labels_path = root / "labels.tsv"
            mapping_rows: list[dict[str, str]] = []
            label_rows: list[dict[str, str]] = []
            source_labels = [
                "Cholangiocytes",
                "Hepatocytes",
                "Kupffer",
                "Mesenchymal",
                "NK-T",
                "B cells",
                "LSEC",
            ]
            counter = 0
            for source_label in source_labels:
                for donor_index in range(3):
                    for replicate in range(3):
                        counter += 1
                        well = f"well{1 + donor_index % 2}"
                        barcode = f"{'ACGT' * 4}{counter:08b}".replace("0", "A").replace("1", "C")
                        cell_id = f"{well}_{barcode}-1"
                        mapping_rows.append(
                            {
                                "cell_id": cell_id,
                                "donor_id": f"donor{donor_index}",
                                "well_id": well,
                                "raw_barcode": cell_id,
                            }
                        )
                        label_rows.append(
                            {"cell_id": cell_id, "author_label": source_label}
                        )
            self.write_tsv(
                mapping_path,
                ("cell_id", "donor_id", "well_id", "raw_barcode"),
                mapping_rows,
            )
            self.write_tsv(
                labels_path, ("cell_id", "author_label"), label_rows
            )
            first = root / "first.tsv"
            second = root / "second.tsv"
            with patch.multiple(
                builder,
                EXPECTED_CELLS=len(mapping_rows),
                EXPECTED_DONORS=3,
                EXPECTED_WELLS=("well1", "well2"),
                CELLS_PER_CLASS=6,
                MAPPING_METADATA_SHA256=builder.sha256_file(mapping_path),
                AUTHOR_LABELS_SHA256=builder.sha256_file(labels_path),
            ):
                builder.select_cells(
                    mapping_metadata=mapping_path,
                    author_labels=labels_path,
                    output=first,
                )
                builder.select_cells(
                    mapping_metadata=mapping_path,
                    author_labels=labels_path,
                    output=second,
                )
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with first.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 30)
            self.assertEqual(len({row["cell_id"] for row in rows}), 30)
            self.assertEqual({row["donor_id"] for row in rows}, {
                "donor0",
                "donor1",
                "donor2",
            })
            self.assertEqual(
                {row["broad_label"] for row in rows},
                set(builder.LABEL_MAP.values()),
            )
            self.assertFalse(
                {row["source_label"] for row in rows} & {"B cells", "LSEC"}
            )
            receipt = json.loads(first.with_suffix(".json").read_text())
            self.assertFalse(receipt["selection_outcomes_used"])
            self.assertFalse(receipt["atac_values_or_qc_used_for_selection"])


if __name__ == "__main__":
    unittest.main()
