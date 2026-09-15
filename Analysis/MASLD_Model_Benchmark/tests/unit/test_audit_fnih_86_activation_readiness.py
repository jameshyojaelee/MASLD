from __future__ import annotations

import copy
from pathlib import Path
import tomllib
import unittest

from scripts.audit_fnih_86_activation_readiness import (
    FNIHReadinessError,
    discover_portal_assets,
    evaluate_gate,
    inventory_candidate_data_paths,
    summarize_public_evidence,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/evaluation/fnih_86_activation_readiness.toml"


class FNIHActivationReadinessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = tomllib.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_topology_counts_and_pairing_are_explicit(self) -> None:
        _, _ = validate_contract(self.contract)
        by_id = {row["assay_id"]: row for row in self.contract["assay_topology"]}
        self.assertEqual(by_id["tenx_multiome_rna_atac"]["donor_n"], 86)
        self.assertEqual(
            by_id["tenx_multiome_rna_atac"]["within_assay_pairing"],
            "same_nucleus",
        )
        self.assertEqual(by_id["droplet_hic"]["donor_n"], 85)
        self.assertEqual(by_id["visium_hd"]["donor_n"], 7)
        self.assertTrue(
            all(not row["public_donor_barcode_join"] for row in by_id.values())
        )

    def test_activation_is_fail_closed(self) -> None:
        _, criteria = validate_contract(self.contract)
        gate = evaluate_gate(criteria)
        self.assertFalse(gate["activation_gate_pass"])
        self.assertEqual(gate["current_activation_state"], "blocked")
        self.assertGreater(gate["nonpassing_criterion_count"], 0)
        self.assertFalse(gate["eligible_as_external_seal_now"])

    def test_all_required_criteria_must_pass_before_manual_revision(self) -> None:
        _, criteria = validate_contract(self.contract)
        changed = copy.deepcopy(criteria)
        for row in changed:
            row["current_state"] = "pass"
        gate = evaluate_gate(changed)
        self.assertTrue(gate["activation_gate_pass"])
        self.assertEqual(
            gate["current_activation_state"], "eligible_for_manual_manifest_revision"
        )
        self.assertFalse(gate["automatic_activation"])
        self.assertFalse(gate["eligible_as_external_seal_now"])

    def test_firewall_flags_cannot_be_relaxed(self) -> None:
        for key in (
            "automatic_activation",
            "sealed_outcomes_read",
            "biological_matrices_downloaded",
            "controlled_data_accessed",
        ):
            changed = copy.deepcopy(self.contract)
            changed[key] = True
            with self.assertRaises(FNIHReadinessError):
                validate_contract(changed)

    def test_portal_discovery_fetches_only_same_origin_static_code(self) -> None:
        html = (
            '<script src="/MASLD/assets/index-abc.js"></script>'
            '<link href="/MASLD/assets/index-def.css" rel="stylesheet">'
        )
        self.assertEqual(
            discover_portal_assets("https://epigenome.wustl.edu/MASLD/", html),
            [
                "https://epigenome.wustl.edu/MASLD/assets/index-abc.js",
                "https://epigenome.wustl.edu/MASLD/assets/index-def.css",
            ],
        )
        with self.assertRaises(FNIHReadinessError):
            discover_portal_assets(
                "https://epigenome.wustl.edu/MASLD/",
                '<script src="https://example.org/other.js"></script>',
            )

    def test_candidate_biological_paths_are_inventory_only(self) -> None:
        rows = inventory_candidate_data_paths(
            "portal_asset_01",
            b'const x="data/example.h5ad"; const y="tracks/sample.bw";',
        )
        self.assertEqual({row["suffix"] for row in rows}, {".h5ad", ".bw"})
        self.assertTrue(all(row["requested"] is False for row in rows))

    def test_code_host_unavailability_is_explicitly_fail_closed(self) -> None:
        sources, _ = validate_contract(self.contract)
        github = next(row for row in sources if row["source_id"] == "github_repository")
        self.assertTrue(github["allow_unavailable"])
        self.assertIn("activation stays blocked", github["unavailable_disposition"])

    def test_public_summary_distinguishes_aggregate_tracks_from_donor_joins(self) -> None:
        payloads = {
            "zenodo_record": b'{"id":15298484,"metadata":{"access_right":"restricted","license":{"id":"cc-by-4.0"}},"files":[]}',
            "gds_by_pmid": b'{"esearchresult":{"count":"0","idlist":[]}}',
            "gds_by_exact_title": b'{"esearchresult":{"count":"0","idlist":[]}}',
        }
        paths = [
            {
                "literal": "https://example.org/FNIH.Liver/ATAC/disease/Hepatocytes_MASH.bw",
                "suffix": ".bw",
            },
            {
                "literal": "https://example.org/FNIH.Liver/H3K27ac/celltypes/HSC.bw",
                "suffix": ".bw",
            },
        ]
        summary = summarize_public_evidence([], payloads, paths)
        self.assertEqual(summary["portal_candidate_suffix_counts"], {".bw": 2})
        self.assertEqual(summary["portal_aggregate_celltype_or_disease_path_count"], 2)
        self.assertEqual(summary["portal_identity_keyword_path_count"], 0)


if __name__ == "__main__":
    unittest.main()
