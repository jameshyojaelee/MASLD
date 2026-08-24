from __future__ import annotations

import unittest

from scripts.audit_gse267145_reference_crosswalk import build_gene_crosswalk


def record(stable: str, version: str) -> dict[str, str]:
    return {
        "stable_id": stable,
        "versioned_id": f"{stable}.{version}",
        "contig": "1",
        "start_1based": "1",
        "end_1based": "2",
        "strand": "+",
        "gene_name": stable,
        "gene_type": "protein_coding",
    }


class GSE267145ReferenceCrosswalkTests(unittest.TestCase):
    def test_stable_id_mapping_does_not_require_equal_versions(self) -> None:
        stable = "ENSG00000000001"
        rows, states = build_gene_crosswalk(
            [stable], {stable: record(stable, "1")}, {stable: record(stable, "9")}
        )
        self.assertEqual(rows[0]["mapping_state"], "stable_id_exact_v98_and_v49")
        self.assertEqual(states["stable_id_exact_v98_and_v49"], 1)

    def test_retired_gene_is_masked(self) -> None:
        stable = "ENSG00000000002"
        rows, _states = build_gene_crosswalk(
            [stable], {stable: record(stable, "1")}, {}
        )
        self.assertEqual(rows[0]["allowed_project_input"], "false")


if __name__ == "__main__":
    unittest.main()
