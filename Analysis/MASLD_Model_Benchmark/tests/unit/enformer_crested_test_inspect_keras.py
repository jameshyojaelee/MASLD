from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
import zipfile


MODULE = Path(__file__).parents[2] / "scripts" / "enformer_crested_inspect_keras.py"
SPEC = importlib.util.spec_from_file_location("enformer_crested_inspect_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
inspector = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inspector)


class EnformerCREstedInspectTests(unittest.TestCase):
    def test_safe_name_rejects_traversal(self) -> None:
        with self.assertRaisesRegex(inspector.KerasInspectionError, "traversal"):
            inspector._safe_name("../model.keras")

    def test_serialized_class_collection_and_lambda_gate(self) -> None:
        rows: set[tuple[str, str, str]] = set()
        inspector._collect_serialized_classes(
            {
                "module": "keras.layers",
                "class_name": "Dense",
                "registered_name": None,
                "config": [{"class_name": "Activation", "module": "keras.layers"}],
            },
            rows,
        )
        self.assertIn(("keras.layers", "Dense", "None"), rows)
        self.assertIn(("keras.layers", "Activation", ""), rows)

    def test_keras_inventory_rejects_lambda(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.keras"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(
                    "config.json",
                    json.dumps({"class_name": "Lambda", "module": "keras.layers"}),
                )
                archive.writestr("metadata.json", json.dumps({"keras_version": "3"}))
                archive.writestr("model.weights.h5", b"fixture")
            with self.assertRaisesRegex(inspector.KerasInspectionError, "Lambda"):
                inspector._inspect_keras(path)

    def test_outer_link_contract_is_not_accepted_as_regular_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.tar.gz"
            with tarfile.open(path, "w:gz") as archive:
                member = tarfile.TarInfo(inspector.KERAS_NAME)
                member.type = tarfile.SYMTYPE
                member.linkname = "elsewhere"
                archive.addfile(member)
            with tarfile.open(path, "r:gz") as archive:
                self.assertFalse(archive.getmembers()[0].isfile())


if __name__ == "__main__":
    unittest.main()
