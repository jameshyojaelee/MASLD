from __future__ import annotations

import copy
import io
import json
from pathlib import Path
import unittest
import zipfile

from scripts.acquire_scooby_native_context_release import (
    ScoobyNativeContextReleaseError,
    extract_member,
    validate_boundary,
)
from scripts.scooby_remote_zip_inventory import BytesRangeReader, inventory_reader


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/scooby_epicardioids_native_context_release.json"


class ScoobyNativeContextReleaseTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_boundary_rejects_coordinate_compatible_liver_claim(self) -> None:
        validate_boundary(self.config)
        mutated = copy.deepcopy(self.config)
        mutated["release_boundary"]["coordinate_compatible_liver_query_supported"] = True
        with self.assertRaises(ScoobyNativeContextReleaseError):
            validate_boundary(mutated)

    def test_selective_deflate_extraction(self) -> None:
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr("released/embedding.pq", b"PAR1fixturePAR1")
        payload = stream.getvalue()
        inventory = inventory_reader(BytesRangeReader(payload))
        member = dict(inventory["members"][0])
        member["destination"] = "embedding.pq"
        extracted, receipt = extract_member(BytesRangeReader(payload), member)
        self.assertEqual(extracted, b"PAR1fixturePAR1")
        self.assertEqual(receipt["uncompressed_size"], len(extracted))

    def test_rna_conditioned_lane_remains_closed(self) -> None:
        mutated = copy.deepcopy(self.config)
        mutated["release_boundary"]["rna_conditioned_atac_eligible"] = True
        with self.assertRaises(ScoobyNativeContextReleaseError):
            validate_boundary(mutated)


if __name__ == "__main__":
    unittest.main()
