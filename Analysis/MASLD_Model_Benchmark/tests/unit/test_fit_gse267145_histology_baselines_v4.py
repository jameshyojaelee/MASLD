from __future__ import annotations

import json
from pathlib import Path
import unittest

import scripts.fit_gse267145_histology_baselines_v4 as v4


ROOT = Path(__file__).resolve().parents[2]


class GSE267145HistologyFitterV4Tests(unittest.TestCase):
    def test_software_contract_preserves_v3_algorithm_and_unscored_state(self) -> None:
        contract = json.loads(v4.SOFTWARE_CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["revision_id"], v4.SOFTWARE_REVISION_ID)
        self.assertFalse(contract["software_correction"]["modeling_code_changed"])
        self.assertFalse(contract["software_correction"]["fallback_code_changed"])
        self.assertEqual(
            contract["preserved_contract"]["v4_units_to_compute"],
            ["1/1721", "1/1723", "1/1733"],
        )
        self.assertFalse(contract["scoring_authorized"])

    def test_repository_root_is_inserted_before_v3_import(self) -> None:
        source = (
            ROOT / "scripts/fit_gse267145_histology_baselines_v4.py"
        ).read_text(encoding="utf-8")
        insert_at = source.index("sys.path.insert(0, ROOT.as_posix())")
        import_at = source.index(
            "import scripts.fit_gse267145_histology_baselines_v3 as v3"
        )
        self.assertLess(insert_at, import_at)
        self.assertIn(ROOT.as_posix(), v4.sys.path)

    def test_v4_delegates_to_the_validated_v3_entrypoint(self) -> None:
        source = (
            ROOT / "scripts/fit_gse267145_histology_baselines_v4.py"
        ).read_text(encoding="utf-8")
        self.assertIn("return v3.main()", source)
        self.assertEqual(
            v4.v3.CONTRACT_SHA256,
            "8be11a742def0f80465bca92fd02baf626d2d207ee30a87cdb25bcf59e47e29f",
        )
        self.assertEqual(
            v4.v3.FALLBACK_ID,
            "training_only_fibrosis_group3_class_prior_v1",
        )


if __name__ == "__main__":
    unittest.main()
