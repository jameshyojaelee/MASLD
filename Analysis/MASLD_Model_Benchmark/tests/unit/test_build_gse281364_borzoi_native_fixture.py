from __future__ import annotations

import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from contextlib import redirect_stdout

from scripts import build_gse281364_borzoi_native_fixture as fixture
from scripts.alphagenome_sei_build_fixture import IndexedFasta


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/gse281364_borzoi_native_screen.json"


def source_group(identifier: str, start: int, end: int, contig: str = "chr1") -> dict[str, str]:
    return {
        "outer_locus_sequence_group_id": identifier,
        "outer_fold": "0",
        "contig": contig,
        "elements": "1",
        "minimum_variant_pos0": str(start),
        "maximum_variant_pos0": str(end),
    }


def write_fasta(path: Path, sequence: str) -> Path:
    header = b">chr1\n"
    with path.open("wb") as handle:
        handle.write(header)
        for start in range(0, len(sequence), 60):
            handle.write((sequence[start : start + 60] + "\n").encode())
    fai = Path(f"{path}.fai")
    fai.write_text(f"chr1\t{len(sequence)}\t{len(header)}\t60\t61\n", encoding="ascii")
    return fai


class BorzoiNativeFixtureTests(unittest.TestCase):
    def test_cli_maps_contract_flag_to_build_contract_path(self) -> None:
        with patch.object(fixture, "build", return_value={"status": "ok"}) as build:
            with redirect_stdout(io.StringIO()):
                status = fixture.main(
                    [
                        "--project-root",
                        "/project",
                        "--contract",
                        "/contract.json",
                        "--output",
                        "/output",
                    ]
                )
        self.assertEqual(status, 0)
        build.assert_called_once_with(
            project_root=Path("/project"),
            contract_path=Path("/contract.json"),
            output=Path("/output"),
        )

    def test_registered_contract_is_four_member_native_and_fail_closed(self) -> None:
        contract = fixture.load_contract(CONTRACT)
        ensemble = contract["ensemble_contract"]
        self.assertEqual(
            ensemble["member_order"],
            ["replicate_0", "replicate_1", "replicate_2", "replicate_3"],
        )
        self.assertEqual(ensemble["orientations"], ["forward", "reverse_complement"])
        self.assertEqual(ensemble["shift_bp"], [0])
        self.assertIn("8_restored", ensemble["aggregation"])
        self.assertIn("separately_named_ports", ensemble["generic_port_rule"])
        self.assertFalse(contract["execution_gate"]["checkpoint_download_allowed"])
        self.assertFalse(contract["execution_gate"]["model_forward_allowed"])

    def test_long_range_groups_merge_overlapping_receptive_fields(self) -> None:
        rows = [
            source_group("outer_00000000000000000001", 1_000_000, 1_000_100),
            source_group("outer_00000000000000000002", 1_400_000, 1_400_100),
            source_group("outer_00000000000000000003", 2_100_000, 2_100_100),
            source_group("outer_00000000000000000004", 3_000_000, 3_000_100, "chr2"),
            source_group("outer_00000000000000000005", 4_000_000, 4_000_100, "chr3"),
            source_group("outer_00000000000000000006", 5_000_000, 5_000_100, "chr4"),
        ]
        components, mapping = fixture.build_long_range_groups(
            rows, buffer_bp=524_288, folds=5, seed="test"
        )
        self.assertEqual(len(components), 5)
        self.assertEqual(mapping[rows[0]["outer_locus_sequence_group_id"]][0], mapping[rows[1]["outer_locus_sequence_group_id"]][0])
        self.assertNotEqual(mapping[rows[1]["outer_locus_sequence_group_id"]][0], mapping[rows[2]["outer_locus_sequence_group_id"]][0])
        self.assertEqual({fold for _, fold in mapping.values()}, set(range(5)))
        fixture.verify_long_range_firewall(components, 524_288)

    def test_ref_coordinate_alt_and_reverse_complement_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reference.fa"
            fai = write_fasta(path, "ACGT" * 400_000)
            reference = IndexedFasta(path, fai)
            position0 = 700_000
            ref = reference.fetch("chr1", position0, position0 + 1)
            alt = "A" if ref != "A" else "C"
            row = {
                "contig": "chr1",
                "variant_pos0": str(position0),
                "variant_pos1": str(position0 + 1),
                "genomic_ref": ref,
                "genomic_alt": alt,
            }
            reference_sequence, alternative, start, end = fixture.allele_window(
                reference,
                row,
                window=524_288,
                center=262_144,
                allowed_contigs={"chr1"},
            )
            self.assertEqual(end - start, 524_288)
            self.assertEqual(reference_sequence[262_144], ref)
            self.assertEqual(alternative[262_144], alt)
            self.assertEqual(
                fixture.reverse_complement(fixture.reverse_complement(reference_sequence)),
                reference_sequence,
            )
            wrong = dict(row, variant_pos1=str(position0 + 2))
            with self.assertRaises(fixture.BorzoiFixtureError):
                fixture.allele_window(
                    reference,
                    wrong,
                    window=524_288,
                    center=262_144,
                    allowed_contigs={"chr1"},
                )

    def test_native_track_map_requires_strand_closed_primary_set(self) -> None:
        rows = [
            ["0", "plus", "x", "1", "1", "1", "sum", "1", "CAGE:liver +"],
            ["1", "minus", "x", "1", "1", "1", "sum", "0", "CAGE:liver -"],
            ["2", "atac", "x", "1", "1", "1", "sum", "2", "ATAC:Hepatocyte"],
            ["3", "rna", "x", "1", "1", "1", "sum", "3", "RNA:liver"],
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "targets.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write("\tidentifier\tfile\tclip\tclip_soft\tscale\tsum_stat\tstrand_pair\tdescription\n")
                for row in rows:
                    handle.write("\t".join(row) + "\n")
            contract = json.loads(CONTRACT.read_text())
            output = contract["native_output_contract"]
            output["human_tracks"] = 4
            output["assay_track_counts"] = {"CAGE": 2, "ATAC": 1, "RNA": 1}
            output["primary_track_indices"] = [0, 1, 2, 3]
            output["primary_track_roles"] = {
                "0": "expression_plus",
                "1": "expression_minus",
                "2": "accessibility",
                "3": "expression",
            }
            output["primary_track_identifiers"] = {
                "0": "plus",
                "1": "minus",
                "2": "atac",
                "3": "rna",
            }
            rendered = fixture.load_native_tracks(path, contract)
            self.assertEqual(len(rendered), 4)
            self.assertEqual(rendered[1]["strand_pair"], "0")
            output["primary_track_indices"] = [0, 2, 3]
            with self.assertRaises(fixture.BorzoiFixtureError):
                fixture.load_native_tracks(path, contract)


if __name__ == "__main__":
    unittest.main()
