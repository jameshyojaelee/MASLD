from __future__ import annotations

import csv
import gzip
from pathlib import Path
import tempfile
import unittest
import unittest.mock

from scripts.build_gpl16686_gencode_v49_crosswalk import (
    CONTRACT_STATES,
    MAPPING_STATE_TO_CONTRACT_STATE,
    CrosswalkError,
    build_crosswalk,
    classify_feature,
    digest_pairs,
    parse_gencode_genes,
    read_platform_axis,
    stable_gene_id,
)


GENCODE = {
    "ENSG1": {"versioned_id": "ENSG1.4", "gene_name": "AAA", "gene_type": "protein_coding", "contig": "chr1"},
    "ENSG2": {"versioned_id": "ENSG2.9", "gene_name": "BBB", "gene_type": "lncRNA", "contig": "chr2"},
}


def classify(feature: str, **overrides: object) -> dict[str, str]:
    arguments: dict[str, object] = {
        "annotation_absent": set(),
        "probe_entrez": {},
        "entrez_absent": set(),
        "entrez_ensembl": {},
        "gencode": GENCODE,
    }
    arguments.update(overrides)
    return classify_feature(feature, **arguments)  # type: ignore[arg-type]


class ClassifyFeatureTests(unittest.TestCase):
    def test_one_to_one_feature_is_eligible(self) -> None:
        row = classify(
            "p1", probe_entrez={"p1": {"100"}}, entrez_ensembl={"100": {"ENSG1"}}
        )
        self.assertEqual(row["mapping_state"], "one_to_one_gencode_v49")
        self.assertEqual(row["contract_state"], "one_to_one")
        self.assertEqual(row["eligible_for_gene_matrix"], "true")
        self.assertEqual(row["ensembl_gene_id"], "ENSG1")
        self.assertEqual(row["gencode_v49_gene_id_versioned"], "ENSG1.4")
        self.assertEqual(row["gencode_v49_gene_name"], "AAA")

    def test_absent_from_annotation_package_is_unmapped(self) -> None:
        row = classify("p1", annotation_absent={"p1"})
        self.assertEqual(
            row["mapping_state"], "platform_feature_absent_from_annotation_package"
        )
        self.assertEqual(row["contract_state"], "unmapped_or_control")
        self.assertEqual(row["eligible_for_gene_matrix"], "false")

    def test_control_feature_without_entrez_is_unmapped(self) -> None:
        row = classify("AFFX-control")
        self.assertEqual(row["mapping_state"], "no_entrez_id_in_annotation_package")
        self.assertEqual(row["contract_state"], "unmapped_or_control")

    def test_multiple_entrez_ids_are_join_unresolved(self) -> None:
        row = classify(
            "p1",
            probe_entrez={"p1": {"100", "200"}},
            entrez_ensembl={"100": {"ENSG1"}, "200": {"ENSG2"}},
        )
        self.assertEqual(row["mapping_state"], "multiple_entrez_ids")
        self.assertEqual(row["contract_state"], "join_unresolved")
        self.assertEqual(row["entrez_id_count"], "2")
        self.assertEqual(row["entrez_id"], "")

    def test_multiple_ensembl_genes_are_join_unresolved(self) -> None:
        row = classify(
            "p1", probe_entrez={"p1": {"100"}}, entrez_ensembl={"100": {"ENSG1", "ENSG2"}}
        )
        self.assertEqual(row["mapping_state"], "multiple_ensembl_genes")
        self.assertEqual(row["contract_state"], "join_unresolved")
        self.assertEqual(row["ensembl_gene_id_count"], "2")
        self.assertEqual(row["ensembl_gene_id"], "")

    def test_entrez_without_ensembl_is_join_unresolved(self) -> None:
        row = classify("p1", probe_entrez={"p1": {"100"}})
        self.assertEqual(row["mapping_state"], "entrez_without_ensembl_gene")
        self.assertEqual(row["contract_state"], "join_unresolved")

    def test_entrez_absent_from_org_db_is_join_unresolved(self) -> None:
        row = classify("p1", probe_entrez={"p1": {"100"}}, entrez_absent={"100"})
        self.assertEqual(row["mapping_state"], "entrez_absent_from_org_hs_eg_db")
        self.assertEqual(row["contract_state"], "join_unresolved")

    def test_retired_gene_is_below_qc_not_join_unresolved(self) -> None:
        row = classify(
            "p1", probe_entrez={"p1": {"100"}}, entrez_ensembl={"100": {"ENSG_RETIRED"}}
        )
        self.assertEqual(row["mapping_state"], "ensembl_gene_absent_from_gencode_v49")
        self.assertEqual(row["contract_state"], "below_qc")
        self.assertEqual(row["ensembl_gene_id"], "ENSG_RETIRED")
        self.assertEqual(row["gencode_v49_gene_id_versioned"], "")

    def test_every_mapping_state_lands_in_a_declared_contract_state(self) -> None:
        self.assertEqual(
            set(MAPPING_STATE_TO_CONTRACT_STATE.values()) - set(CONTRACT_STATES), set()
        )


class DigestTests(unittest.TestCase):
    def test_digest_is_order_sensitive_and_content_sensitive(self) -> None:
        first = digest_pairs([("p1", "ENSG1"), ("p2", "ENSG2")])
        self.assertNotEqual(first, digest_pairs([("p2", "ENSG2"), ("p1", "ENSG1")]))
        self.assertNotEqual(first, digest_pairs([("p1", "ENSG1"), ("p2", "ENSG3")]))
        self.assertEqual(first, digest_pairs([("p1", "ENSG1"), ("p2", "ENSG2")]))

    def test_stable_gene_id_removes_only_version(self) -> None:
        self.assertEqual(stable_gene_id("ENSG1.12"), "ENSG1")
        self.assertEqual(stable_gene_id("ENSG1"), "ENSG1")


class FixtureBuildTests(unittest.TestCase):
    """Exercise the whole builder against a miniature frozen axis."""

    def _write(self, path: Path, header: list[str], rows: list[list[str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)

    def _fixture(self, base: Path, features: list[str]) -> dict[str, Path]:
        axis = base / "axis.tsv"
        self._write(
            axis,
            ["platform_feature_id"],
            [[feature] for feature in features],
        )
        probe_entrez = base / "probe_entrez.tsv"
        self._write(
            probe_entrez,
            ["platform_feature_id", "entrez_id"],
            [["p1", "100"], ["p2", "100"], ["p3", "200"], ["p3", "300"], ["p4", "400"]],
        )
        entrez_ensembl = base / "entrez_ensembl.tsv"
        self._write(
            entrez_ensembl,
            ["entrez_id", "ensembl_gene_id"],
            [["100", "ENSG1"], ["200", "ENSG2"], ["300", "ENSG2"], ["400", "ENSG_GONE"]],
        )
        annotation_absent = base / "absent.tsv"
        self._write(annotation_absent, ["platform_feature_id", "state"], [["p6", "x"]])
        entrez_absent = base / "entrez_absent.tsv"
        self._write(entrez_absent, ["entrez_id", "state"], [])
        versions = base / "versions.tsv"
        self._write(
            versions,
            ["package", "version"],
            [["hgu133plus2.db", "3.13.0"], ["org.Hs.eg.db", "3.22.0"]],
        )
        gtf = base / "genes.gtf.gz"
        with gzip.open(gtf, "wt", encoding="utf-8") as handle:
            handle.write("##description: test\n")
            handle.write(
                'chr1\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "ENSG1.4"; '
                'gene_type "protein_coding"; gene_name "AAA";\n'
            )
            handle.write(
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t+\t.\tgene_id "ENSG1.4"; '
                'gene_type "protein_coding"; gene_name "AAA";\n'
            )
            handle.write(
                'chr2\tHAVANA\tgene\t30\t40\t.\t-\t.\tgene_id "ENSG2.9"; '
                'gene_type "lncRNA"; gene_name "BBB";\n'
            )
        return {
            "platform_axis": axis,
            "probe_entrez_path": probe_entrez,
            "entrez_ensembl_path": entrez_ensembl,
            "annotation_absent_path": annotation_absent,
            "entrez_absent_path": entrez_absent,
            "gencode_gtf": gtf,
            "annotation_versions_path": versions,
        }

    def test_builder_states_counts_and_collapse_map(self) -> None:
        features = ["p1", "p2", "p3", "p4", "p5", "p6"]
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            arguments = self._fixture(base, features)
            output = base / "crosswalk"
            with unittest.mock.patch(
                "scripts.build_gpl16686_gencode_v49_crosswalk.EXPECTED_PLATFORM_FEATURES",
                len(features),
            ):
                receipt = build_crosswalk(output=output, **arguments)

            self.assertEqual(receipt["platform_features"], 6)
            self.assertEqual(
                receipt["contract_state_counts"],
                {
                    "one_to_one": 2,
                    "join_unresolved": 1,
                    "below_qc": 1,
                    "unmapped_or_control": 2,
                },
            )
            self.assertEqual(receipt["eligible_platform_features"], 2)
            self.assertEqual(receipt["one_to_one_genes"], 1)
            self.assertEqual(receipt["genes_with_multiple_platform_features"], 1)
            self.assertEqual(receipt["maximum_platform_features_per_gene"], 2)
            self.assertIs(receipt["depends_on_labels"], False)
            self.assertIs(receipt["depends_on_CEL_intensity"], False)
            self.assertIs(receipt["depends_on_held_cohort_arrays"], False)
            self.assertIs(receipt["coordinate_or_GB_ACC_mapping_inference_used"], False)

            with (output / "gpl16686_gencode_v49_crosswalk.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["platform_feature_id"] for row in rows], features)
            by_feature = {row["platform_feature_id"]: row for row in rows}
            self.assertEqual(by_feature["p3"]["contract_state"], "join_unresolved")
            self.assertEqual(by_feature["p4"]["contract_state"], "below_qc")
            self.assertEqual(by_feature["p5"]["mapping_state"], "no_entrez_id_in_annotation_package")
            self.assertEqual(
                by_feature["p6"]["mapping_state"],
                "platform_feature_absent_from_annotation_package",
            )

            with (output / "gpl16686_gene_collapse_map.tsv").open(encoding="utf-8") as handle:
                collapse = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(collapse), 1)
            self.assertEqual(collapse[0]["ensembl_gene_id"], "ENSG1")
            self.assertEqual(collapse[0]["platform_feature_ids"], "p1;p2")
            self.assertEqual(collapse[0]["entrez_id_count"], "1")

    def test_builder_refuses_to_overwrite_an_existing_output(self) -> None:
        features = ["p1", "p2", "p3", "p4", "p5", "p6"]
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            arguments = self._fixture(base, features)
            output = base / "crosswalk"
            output.mkdir()
            with unittest.mock.patch(
                "scripts.build_gpl16686_gencode_v49_crosswalk.EXPECTED_PLATFORM_FEATURES",
                len(features),
            ):
                with self.assertRaises(FileExistsError):
                    build_crosswalk(output=output, **arguments)

    def test_axis_of_unexpected_size_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            arguments = self._fixture(base, ["p1", "p2"])
            with self.assertRaisesRegex(CrosswalkError, "platform axis differs"):
                read_platform_axis(arguments["platform_axis"])

    def test_duplicate_gencode_stable_ids_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "dup.gtf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(
                    'chr1\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "ENSG1.4"; '
                    'gene_type "protein_coding"; gene_name "AAA";\n'
                )
                handle.write(
                    'chrX\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "ENSG1.5"; '
                    'gene_type "protein_coding"; gene_name "AAA";\n'
                )
            with self.assertRaisesRegex(CrosswalkError, "duplicated"):
                parse_gencode_genes(path)


if __name__ == "__main__":
    unittest.main()
