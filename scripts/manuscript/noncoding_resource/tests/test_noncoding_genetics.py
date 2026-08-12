#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest


TEST_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = TEST_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from noncoding_genetics import (  # noqa: E402
    INPUT_ROLES,
    NoncodingGeneticsError,
    build_release,
)


def write_table(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fields, delimiter=delimiter, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def identity_row(
    gene_id: str,
    base_id: str,
    name: str,
    gene_type: str,
    start: int,
    end: int,
    strand: str,
) -> dict[str, object]:
    return {
        "annotation_release": "GENCODE v49",
        "gene_id_versioned": gene_id,
        "gene_id_base": base_id,
        "gene_name": name,
        "gene_type": gene_type,
        "chromosome": "chr1",
        "start_1based": start,
        "end_1based": end,
        "strand": strand,
        "is_canonical_chromosome": "true",
        "base_id_mapping_status": "unique",
    }


class Fixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.paths = {role: root / f"{role}.tsv" for role in INPUT_ROLES}
        self.paths["credible_sets"] = root / "credible_sets.csv"
        self.paths["consequence_annotation"] = root / "consequence.csv"
        self.paths["coloc"] = root / "coloc.csv"
        self.paths["variant_liftover"] = root / "variant_liftover.csv"
        self.paths["abc_context"] = root / "abc.csv"
        self.paths["atac_context"] = root / "atac.csv"
        self.promotion_manifest = root / "promotion_manifest.tsv"
        self._write_inputs()
        self.write_manifest()

    def _write_inputs(self) -> None:
        write_table(
            self.paths["trait_registry"],
            ["study_name", "trait", "tier", "tier_label", "placement"],
            [
                {
                    "study_name": "DIRECT",
                    "trait": "PDFF",
                    "tier": 1,
                    "tier_label": "direct_MASLD",
                    "placement": "main",
                },
                {
                    "study_name": "ENZYME",
                    "trait": "ALT",
                    "tier": 2,
                    "tier_label": "liver_enzyme",
                    "placement": "main",
                },
                {
                    "study_name": "SUPP",
                    "trait": "Cirrhosis",
                    "tier": 3,
                    "tier_label": "sequela_mixed_etiology",
                    "placement": "supp",
                },
            ],
        )
        write_table(
            self.paths["gene_identity"],
            [
                "annotation_release",
                "gene_id_versioned",
                "gene_id_base",
                "gene_name",
                "gene_type",
                "chromosome",
                "start_1based",
                "end_1based",
                "strand",
                "is_canonical_chromosome",
                "base_id_mapping_status",
            ],
            [
                identity_row(
                    "ENSG00000000001.4",
                    "ENSG00000000001",
                    "PC1",
                    "protein_coding",
                    1000,
                    2000,
                    "+",
                ),
                identity_row(
                    "ENSG00000000002.2",
                    "ENSG00000000002",
                    "LNC1",
                    "lncRNA",
                    10000,
                    11000,
                    "+",
                ),
                identity_row(
                    "ENSG00000000003.1",
                    "ENSG00000000003",
                    "OTH1",
                    "pseudogene",
                    20000,
                    21000,
                    "-",
                ),
            ],
        )
        cs_fields = [
            "chromosome",
            "position",
            "allele1",
            "allele2",
            "trait",
            "locus",
            "study",
            "ancestry",
            "susie_pip",
            "susie_cs",
            "susie_converged",
            "susie_reliable",
            "variant_id",
        ]
        cs_rows = []
        direct_specs = (
            (100, 0.20),
            (200, 0.20),
            (300, 0.20),
            (400, 0.25),
            (500, 0.10),
        )
        for position, pip in direct_specs:
            cs_rows.append(
                {
                    "chromosome": 1,
                    "position": position,
                    "allele1": "A",
                    "allele2": "G",
                    "trait": "PDFF",
                    "locus": "1.100",
                    "study": "DIRECT",
                    "ancestry": "EUR",
                    "susie_pip": pip,
                    "susie_cs": 1,
                    "susie_converged": "TRUE",
                    "susie_reliable": "TRUE",
                    "variant_id": f"1:{position}:A:G",
                }
            )
        for position, pip in ((600, 0.4), (700, 0.4)):
            cs_rows.append(
                {
                    "chromosome": 1,
                    "position": position,
                    "allele1": "C",
                    "allele2": "T",
                    "trait": "ALT",
                    "locus": "1.600",
                    "study": "ENZYME",
                    "ancestry": "EUR",
                    "susie_pip": pip,
                    "susie_cs": 2,
                    "susie_converged": "TRUE",
                    "susie_reliable": "TRUE",
                    "variant_id": f"1:{position}:C:T",
                }
            )
        write_table(self.paths["credible_sets"], cs_fields, cs_rows)
        write_table(
            self.paths["consequence_annotation"],
            ["variant_id", "Consequence", "class"],
            [
                {
                    "variant_id": "1:100:A:G",
                    "Consequence": "missense_variant",
                    "class": "coding_protein_altering",
                },
                {
                    "variant_id": "1:200:A:G",
                    "Consequence": "splice_donor_variant",
                    "class": "splice_region",
                },
                {
                    "variant_id": "1:300:A:G",
                    "Consequence": "3_prime_UTR_variant",
                    "class": "noncoding",
                },
                {
                    "variant_id": "1:400:A:G",
                    "Consequence": "intron_variant",
                    "class": "noncoding",
                },
                {
                    "variant_id": "1:600:C:T",
                    "Consequence": "intergenic_variant",
                    "class": "noncoding",
                },
                {
                    "variant_id": "1:700:C:T",
                    "Consequence": "intergenic_variant",
                    "class": "noncoding",
                },
            ],
        )
        lift_rows = []
        for position in (100, 200, 300, 400, 600, 700):
            allele1, allele2 = ("A", "G") if position < 600 else ("C", "T")
            lift_rows.append(
                {
                    "chromosome": 1,
                    "position": position,
                    "allele1": allele1,
                    "allele2": allele2,
                    "chr_hg38": "chr1",
                    "pos_hg38": 2500 if position == 200 else 50000 + position,
                    "overlaps_any_peak": "TRUE" if position == 100 else "FALSE",
                    "cell_types_overlapping": "Hepatocyte" if position == 100 else "",
                }
            )
        write_table(
            self.paths["variant_liftover"],
            [
                "chromosome",
                "position",
                "allele1",
                "allele2",
                "chr_hg38",
                "pos_hg38",
                "overlaps_any_peak",
                "cell_types_overlapping",
            ],
            lift_rows,
        )
        write_table(
            self.paths["abc_context"],
            ["variant_id", "abc_is_self_promoter"],
            [
                {"variant_id": "1:300:A:G", "abc_is_self_promoter": "FALSE"},
                {"variant_id": "1:400:A:G", "abc_is_self_promoter": "TRUE"},
            ],
        )
        write_table(
            self.paths["atac_context"],
            ["variant_id", "cell_type"],
            [{"variant_id": "1:100:A:G", "cell_type": "Hepatocyte"}],
        )
        write_table(
            self.paths["coloc"],
            [
                "gwas_name",
                "gene",
                "ensembl",
                "PP.H4.abf",
                "PP.H4.susie",
                "ancestry",
                "top_snp",
            ],
            [
                {
                    "gwas_name": "DIRECT",
                    "gene": "LNC1",
                    "ensembl": "ENSG00000000002",
                    "PP.H4.abf": 0.4,
                    "PP.H4.susie": 0.8,
                    "ancestry": "EUR",
                    "top_snp": "1:100",
                },
                {
                    "gwas_name": "ENZYME",
                    "gene": "PC1",
                    "ensembl": "ENSG00000000001.4",
                    "PP.H4.abf": 0.9,
                    "PP.H4.susie": 0.7,
                    "ancestry": "EUR",
                    "top_snp": "1:600",
                },
                {
                    "gwas_name": "DIRECT",
                    "gene": "OTH1",
                    "ensembl": "ENSG00000000003",
                    "PP.H4.abf": 0.1,
                    "PP.H4.susie": "",
                    "ancestry": "EUR",
                    "top_snp": "1:300",
                },
            ],
        )
        write_table(
            self.paths["evidence_classes"],
            ["gene_id_versioned", "primary_evidence_class", "joint_testable"],
            [
                {
                    "gene_id_versioned": "ENSG00000000001.4",
                    "primary_evidence_class": "genetic_only",
                    "joint_testable": "true",
                },
                {
                    "gene_id_versioned": "ENSG00000000002.2",
                    "primary_evidence_class": "disease_state_only",
                    "joint_testable": "true",
                },
                {
                    "gene_id_versioned": "ENSG00000000003.1",
                    "primary_evidence_class": "neither",
                    "joint_testable": "false",
                },
            ],
        )

    def write_manifest(self, terminal_tasks: int = 1100) -> None:
        rows = []
        for role in INPUT_ROLES:
            path = self.paths[role]
            rows.append(
                {
                    "release_id": "fixture-release",
                    "input_role": role,
                    "source_path": str(path.resolve()),
                    "source_sha256": sha256(path),
                    "promotion_status": "promoted",
                    "n_expected_tasks": 1100 if role == "coloc" else "",
                    "n_terminal_tasks": terminal_tasks if role == "coloc" else "",
                }
            )
        write_table(
            self.promotion_manifest,
            [
                "release_id",
                "input_role",
                "source_path",
                "source_sha256",
                "promotion_status",
                "n_expected_tasks",
                "n_terminal_tasks",
            ],
            rows,
        )

    def build(self, output_name: str = "candidate") -> Path:
        output = self.root / output_name
        build_release(
            output_dir=output,
            promotion_manifest=self.promotion_manifest,
            credible_sets=self.paths["credible_sets"],
            consequence_annotation=self.paths["consequence_annotation"],
            trait_registry=self.paths["trait_registry"],
            gene_identity=self.paths["gene_identity"],
            coloc=self.paths["coloc"],
            evidence_classes=self.paths["evidence_classes"],
            variant_liftover=self.paths["variant_liftover"],
            abc_context=self.paths["abc_context"],
            atac_context=self.paths["atac_context"],
        )
        return output


class NoncodingGeneticsTests(unittest.TestCase):
    def test_build_separates_architecture_context_and_target_biotype(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            output = fixture.build()

            architecture = read_tsv(output / "credible_set_pip_architecture.tsv")
            self.assertEqual(len(architecture), 2)
            direct = next(row for row in architecture if row["study"] == "DIRECT")
            enzyme = next(row for row in architecture if row["study"] == "ENZYME")
            self.assertEqual(direct["pip_sum_gate_passed"], "true")
            self.assertAlmostEqual(
                float(direct["protein_altering_pip_mass"]), 0.20 / 0.95
            )
            self.assertAlmostEqual(
                float(direct["canonical_splice_pip_mass"]), 0.20 / 0.95
            )
            self.assertAlmostEqual(
                float(direct["synonymous_or_utr_pip_mass"]), 0.20 / 0.95
            )
            self.assertAlmostEqual(
                float(direct["other_noncoding_pip_mass"]), 0.25 / 0.95
            )
            self.assertAlmostEqual(float(direct["unresolved_pip_mass"]), 0.10 / 0.95)
            self.assertAlmostEqual(float(direct["mass_sum_check"]), 1.0)
            self.assertEqual(enzyme["pip_sum_gate_passed"], "false")
            self.assertEqual(enzyme["protein_altering_pip_mass"], "")

            context = next(
                row
                for row in read_tsv(output / "credible_set_regulatory_context.tsv")
                if row["study"] == "DIRECT"
            )
            self.assertAlmostEqual(
                float(context["promoter_proximal_pip_mass"]), 0.20 / 0.95
            )
            self.assertAlmostEqual(float(context["abc_enhancer_pip_mass"]), 0.20 / 0.95)
            self.assertAlmostEqual(
                float(context["lineage_accessible_pip_mass"]), 0.20 / 0.95
            )
            self.assertAlmostEqual(
                float(context["unresolved_context_pip_mass"]), 0.10 / 0.95
            )
            self.assertEqual(context["lineage_cell_types"], "Hepatocyte")

            targets = read_tsv(output / "coloc_target_biotype.tsv")
            self.assertEqual(len(targets), 3)
            lnc = next(row for row in targets if row["gene_biotype"] == "lncRNA")
            self.assertEqual(lnc["method"], "susie")
            self.assertEqual(lnc["gene_id_versioned"], "ENSG00000000002.2")
            self.assertEqual(
                {row["trait_scope"] for row in targets},
                {"tier1_direct_masld_pdff", "tier2_liver_enzyme"},
            )

            links = read_tsv(output / "noncoding_dna_lncrna_link_status.tsv")
            self.assertTrue(all(row["credible_set_id"] == "" for row in links))
            self.assertTrue(
                all(
                    row["credible_set_to_gene_link_status"]
                    == "not_established_no_explicit_signal_pair"
                    for row in links
                )
            )
            self.assertNotIn(
                "top_snp", (output / "noncoding_dna_lncrna_link_status.tsv").read_text()
            )

            verdicts = read_tsv(output / "genetics_noncoding_verdict.tsv")
            pip_gate = next(
                row
                for row in verdicts
                if row["metric"] == "reliable_credible_set_pip_sum_coverage"
            )
            self.assertEqual(pip_gate["value"], "0.5")
            self.assertEqual(pip_gate["passed"], "false")

    def test_promotion_gate_requires_all_1100_terminal_tasks(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.write_manifest(terminal_tasks=1099)
            with self.assertRaisesRegex(
                NoncodingGeneticsError, "all corrected COLOC tasks must be terminal"
            ):
                fixture.build()

    def test_manifest_hash_and_nonoverwriting_output_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            output = fixture.build()
            with self.assertRaisesRegex(NoncodingGeneticsError, "already exists"):
                fixture.build()
            self.assertTrue((output / "checksum_manifest.tsv").is_file())
            fixture.paths["coloc"].write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(NoncodingGeneticsError, "SHA256 differs"):
                fixture.build("second_candidate")


if __name__ == "__main__":
    unittest.main()
