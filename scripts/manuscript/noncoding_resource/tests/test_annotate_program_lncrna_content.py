#!/usr/bin/env python3
"""Focused tests for frozen-program lncRNA annotation."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
from collections import Counter
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCRIPT_DIR = HERE.parent
PROJECT = SCRIPT_DIR.parents[2]
SCRIPT = SCRIPT_DIR / "annotate_program_lncrna_content.py"
SPEC = importlib.util.spec_from_file_location("program_lncrna", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


MEMBERSHIP_FIELDS = list(MODULE.MEMBERSHIP_REQUIRED)
IDENTITY_FIELDS = list(MODULE.IDENTITY_REQUIRED)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def member(
    uid: str, module: str, gene: str, mapped: str, status: str, weight: str
) -> dict[str, str]:
    return {
        "cell_type": "hepatocytes",
        "module": module,
        "source_gene": gene,
        "canonical_gene": gene,
        "source_weight": weight,
        "mapped_symbol": mapped,
        "mapped_symbol_status": status,
        "original_l1_weight": weight,
        "canonical_weight_text": weight,
        "membership_sha256": ("a" if uid == "p1" else "b") * 64,
        "program_uid": uid,
    }


def identity(
    gene_id: str,
    symbol: str,
    biotype: str,
    n_symbol_ids: int = 1,
    symbol_status: str = "unique",
) -> dict[str, str]:
    return {
        "annotation_release": "GENCODE v49",
        "gene_id_versioned": gene_id,
        "gene_id_base": gene_id.split(".")[0],
        "gene_version": gene_id.rsplit(".", 1)[1],
        "gene_name": symbol,
        "gene_type": biotype,
        "chromosome": "chr1",
        "start_1based": "1",
        "end_1based": "100",
        "strand": "+",
        "source": "HAVANA",
        "level": "2",
        "is_canonical_chromosome": "true",
        "n_versions_for_base_id": "1",
        "base_id_mapping_status": "unique",
        "n_gene_ids_for_symbol": str(n_symbol_ids),
        "symbol_mapping_status": symbol_status,
    }


class ProgramLncRNAAnnotationTests(unittest.TestCase):
    def fixture(self, root: Path) -> tuple[Path, Path]:
        membership = root / "membership.tsv"
        identities = root / "identity.tsv"
        write_tsv(
            membership,
            MEMBERSHIP_FIELDS,
            [
                member(
                    "p1",
                    "1",
                    "LNC1",
                    "LNC1",
                    "gencode_v49_unique_symbol_confirmed",
                    "0.03",
                ),
                member(
                    "p1",
                    "1",
                    "PC1",
                    "PC1",
                    "gencode_v49_unique_symbol_confirmed",
                    "0.77",
                ),
                member("p1", "1", "DUP", "NA", "unmapped_or_ambiguous", "0.20"),
                member(
                    "p2",
                    "2",
                    "LNC1",
                    "LNC1",
                    "gencode_v49_unique_symbol_confirmed",
                    "0.06",
                ),
                member(
                    "p2",
                    "2",
                    "ENSG00000000005",
                    "LNC2",
                    "gencode_v49_unambiguous_ensembl_to_symbol",
                    "0.04",
                ),
                member("p2", "2", "MISSING", "NA", "unmapped_or_ambiguous", "0.90"),
            ],
        )
        write_tsv(
            identities,
            IDENTITY_FIELDS,
            [
                identity("ENSG00000000001.1", "LNC1", "lncRNA"),
                identity("ENSG00000000002.1", "PC1", "protein_coding"),
                identity("ENSG00000000003.1", "DUP", "lncRNA", 2, "duplicated"),
                identity("ENSG00000000004.1", "DUP", "protein_coding", 2, "duplicated"),
                identity("ENSG00000000005.1", "LNC2", "lncRNA"),
            ],
        )
        return membership, identities

    def args(self, membership: Path, identities: Path, output: Path):
        return type(
            "Args",
            (),
            {
                "membership": membership,
                "identity": identities,
                "output": output,
                "expected_programs": 2,
                "expected_rows": 6,
                "expected_membership_sha256": None,
                "lncrna_weight_threshold": 0.05,
            },
        )()

    def test_preserves_membership_and_summarizes_lncrna_weight(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, identities = self.fixture(root)
            output = root / "output"
            MODULE.run(self.args(membership, identities, output))
            source = read_tsv(membership)
            annotated = read_tsv(output / "program_membership_lncrna_annotation.tsv")
            self.assertEqual(len(source), len(annotated))
            for before, after in zip(source, annotated):
                self.assertEqual(
                    before, {field: after[field] for field in MEMBERSHIP_FIELDS}
                )
            summary = {
                row["program_uid"]: row
                for row in read_tsv(output / "program_lncrna_content.tsv")
            }
            self.assertEqual(summary["p1"]["n_uniquely_mapped_lncrna"], "1")
            self.assertAlmostEqual(float(summary["p1"]["lncrna_l1_fraction"]), 0.03)
            self.assertEqual(
                summary["p1"]["requires_leave_all_lncrna_out_sensitivity"], "false"
            )
            self.assertEqual(summary["p1"]["n_ambiguous_members"], "1")
            self.assertEqual(summary["p2"]["n_uniquely_mapped_lncrna"], "2")
            self.assertAlmostEqual(float(summary["p2"]["lncrna_l1_fraction"]), 0.10)
            self.assertEqual(
                summary["p2"]["requires_leave_all_lncrna_out_sensitivity"], "true"
            )
            self.assertEqual(summary["p2"]["n_missing_members"], "1")

    def test_identity_drift_rejects_frozen_unique_symbol(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, identities = self.fixture(root)
            rows = read_tsv(identities)
            rows[0]["n_gene_ids_for_symbol"] = "2"
            rows[0]["symbol_mapping_status"] = "duplicated"
            rows.append(
                identity("ENSG00000000006.1", "LNC1", "lncRNA", 2, "duplicated")
            )
            write_tsv(identities, IDENTITY_FIELDS, rows)
            with self.assertRaisesRegex(
                MODULE.AnnotationError, "unique-symbol mapping drift"
            ):
                MODULE.run(self.args(membership, identities, root / "output"))

    def test_refuses_existing_output_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, identities = self.fixture(root)
            output = root / "output"
            output.mkdir()
            with self.assertRaisesRegex(
                MODULE.AnnotationError, "refusing to overwrite"
            ):
                MODULE.run(self.args(membership, identities, output))

    def test_rejects_membership_hash_drift(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            membership, identities = self.fixture(root)
            args = self.args(membership, identities, root / "output")
            args.expected_membership_sha256 = "0" * 64
            with self.assertRaisesRegex(
                MODULE.AnnotationError, "membership SHA256 drift"
            ):
                MODULE.run(args)

    def test_real_frozen_registry_smoke(self) -> None:
        membership = (
            PROJECT / "Analysis/Multimodal_Program_Projection/candidates/"
            "program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
        )
        identity_path = PROJECT / "data/gencode_v49_gene_metadata.tsv.gz"
        if not membership.is_file() or not identity_path.is_file():
            self.skipTest("repository integration inputs unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            normalized_identity = root / "gene_identity.tsv"
            with gzip.open(identity_path, "rt", encoding="utf-8", newline="") as handle:
                legacy = list(csv.DictReader(handle, delimiter="\t"))
            cardinality = Counter(row["gene_name"] for row in legacy)
            normalized_rows = [
                identity(
                    row["gene_id"],
                    row["gene_name"],
                    row["gene_biotype"],
                    cardinality[row["gene_name"]],
                    "unique" if cardinality[row["gene_name"]] == 1 else "duplicated",
                )
                | {
                    "chromosome": row["chromosome"],
                    "is_canonical_chromosome": "true"
                    if row["chromosome"]
                    in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
                    else "false",
                }
                for row in legacy
            ]
            write_tsv(normalized_identity, IDENTITY_FIELDS, normalized_rows)
            args = type(
                "Args",
                (),
                {
                    "membership": membership,
                    "identity": normalized_identity,
                    "output": root / "output",
                    "expected_programs": 117,
                    "expected_rows": 7093,
                    "expected_membership_sha256": MODULE.FROZEN_MEMBERSHIP_SHA256,
                    "lncrna_weight_threshold": 0.05,
                },
            )()
            MODULE.run(args)
            summary = read_tsv(args.output / "program_lncrna_content.tsv")
            annotated = read_tsv(
                args.output / "program_membership_lncrna_annotation.tsv"
            )
            self.assertEqual(len(summary), 117)
            self.assertEqual(
                sum(row["is_uniquely_mapped_lncrna"] == "true" for row in annotated),
                409,
            )
            self.assertEqual(
                sum(row["identity_mapping_status"] == "ambiguous" for row in annotated),
                458,
            )
            self.assertEqual(
                sum(row["identity_mapping_status"] == "missing" for row in annotated), 2
            )
            self.assertEqual(
                sum(
                    row["requires_leave_all_lncrna_out_sensitivity"] == "true"
                    for row in summary
                ),
                49,
            )
            self.assertTrue(
                all(
                    abs(float(row["total_original_l1_weight"]) - 1.0) <= 1e-12
                    for row in summary
                )
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
