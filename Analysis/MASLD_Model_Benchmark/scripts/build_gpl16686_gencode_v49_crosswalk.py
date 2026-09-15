#!/usr/bin/env python3
"""Freeze the GPL16686 platform feature to GENCODE v49 gene crosswalk.

The crosswalk is the training-frozen feature identity object named by
``config/evaluation/microarray_transfer_preprocessing.toml``.  It is built from
the frozen GEO platform axis and the exact annotation package versions only.
No CEL intensity, no GEO series matrix, no label, and no outcome participates,
so the object is identical whichever arrays are later held out.

Mapping authority, in the requirements' order:

    hgu133plus2.db exact version -> org.Hs.eg.db exact version -> GENCODE v49

Ambiguity is never broken by a heuristic.  A platform feature reaches the gene
matrix only when it carries exactly one Entrez ID and that Entrez ID carries
exactly one Ensembl gene that is present in GENCODE v49.  Everything else keeps
an explicit audit state: ``join_unresolved`` for an unresolved join and
``below_qc`` for a gene that GENCODE v49 retired or never carried.  Coordinates
and GB_ACC are never used to infer a mapping.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path

EXPECTED_PLATFORM_FEATURES = 53_617
PLATFORM_ID = "GPL16686"
ANNOTATION_PACKAGE = "hugene20sttranscriptcluster.db"

CONTRACT_STATES = ("one_to_one", "join_unresolved", "below_qc", "unmapped_or_control")

MAPPING_STATE_TO_CONTRACT_STATE = {
    "one_to_one_gencode_v49": "one_to_one",
    "platform_feature_absent_from_annotation_package": "unmapped_or_control",
    "no_entrez_id_in_annotation_package": "unmapped_or_control",
    "multiple_entrez_ids": "join_unresolved",
    "entrez_absent_from_org_hs_eg_db": "join_unresolved",
    "entrez_without_ensembl_gene": "join_unresolved",
    "multiple_ensembl_genes": "join_unresolved",
    "ensembl_gene_absent_from_gencode_v49": "below_qc",
}

CROSSWALK_COLUMNS = (
    "platform_feature_id",
    "mapping_state",
    "contract_state",
    "eligible_for_gene_matrix",
    "entrez_id_count",
    "entrez_id",
    "ensembl_gene_id_count",
    "ensembl_gene_id",
    "gencode_v49_gene_id_versioned",
    "gencode_v49_gene_name",
    "gencode_v49_gene_type",
    "gencode_v49_contig",
)


class CrosswalkError(RuntimeError):
    """Raised when a frozen input or a derived state differs from the requirement."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def digest_pairs(pairs: list[tuple[str, str]]) -> str:
    """Digest an ordered feature-to-gene relation independently of file layout."""
    value = sha256(b"masld-bench-gpl16686-gencode-v49-crosswalk-v1\0")
    for left, right in pairs:
        value.update(left.encode("utf-8"))
        value.update(b"\t")
        value.update(right.encode("utf-8"))
        value.update(b"\n")
    return value.hexdigest()


def stable_gene_id(value: str) -> str:
    return value.split(".", 1)[0]


def parse_attributes(value: str) -> dict[str, str]:
    output: dict[str, str] = {}
    for token in value.rstrip(";").split(";"):
        token = token.strip()
        if not token:
            continue
        key, separator, observed = token.partition(" ")
        if not separator:
            raise CrosswalkError("GTF attribute differs")
        output[key] = observed.strip().strip('"')
    return output


def parse_gencode_genes(path: Path) -> dict[str, dict[str, str]]:
    """Return GENCODE gene records keyed by unversioned Ensembl gene ID."""
    genes: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise CrosswalkError(f"GTF width differs at line {line_number}")
            if fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            versioned = attributes.get("gene_id")
            if not versioned:
                raise CrosswalkError("GENCODE gene lacks gene_id")
            stable = stable_gene_id(versioned)
            record = {
                "versioned_id": versioned,
                "gene_name": attributes.get("gene_name", ""),
                "gene_type": attributes.get("gene_type", ""),
                "contig": fields[0],
            }
            if stable in genes and genes[stable] != record:
                raise CrosswalkError(f"GENCODE stable gene ID is duplicated: {stable}")
            genes[stable] = record
    if not genes:
        raise CrosswalkError("GENCODE GTF contains no genes")
    return genes


def read_platform_axis(path: Path) -> list[str]:
    """Read the transcript-cluster axis emitted by the annotation dump.

    The GEO GPL16686 table carries 53981 rows, but the core summarizer returns
    53617 transcript clusters and the matrix is keyed on what the summarizer
    produces.  The annotation package's PROBEID keys are that axis, so the axis
    is taken from there rather than from the GEO table.
    """
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "platform_feature_id" not in reader.fieldnames:
            raise CrosswalkError("platform axis lacks a platform_feature_id column")
        features = [row["platform_feature_id"] for row in reader]
    if len(features) != EXPECTED_PLATFORM_FEATURES or len(set(features)) != len(features):
        raise CrosswalkError("GPL16686 platform axis differs")
    return features


def read_relation(path: Path, left: str, right: str) -> dict[str, set[str]]:
    relation: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or left not in reader.fieldnames or right not in reader.fieldnames:
            raise CrosswalkError(f"relation {path.name} lacks {left} or {right}")
        for row in reader:
            key = row[left].strip()
            value = row[right].strip()
            if not key or not value:
                raise CrosswalkError(f"relation {path.name} carries an empty identifier")
            relation[key].add(value)
    return dict(relation)


def read_single_column(path: Path, column: str) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise CrosswalkError(f"{path.name} lacks {column}")
        return {row[column].strip() for row in reader if row[column].strip()}


def classify_feature(
    feature: str,
    *,
    annotation_absent: set[str],
    probe_entrez: dict[str, set[str]],
    entrez_absent: set[str],
    entrez_ensembl: dict[str, set[str]],
    gencode: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Apply the frozen state machine to one platform feature."""
    record = {
        "platform_feature_id": feature,
        "mapping_state": "",
        "contract_state": "",
        "eligible_for_gene_matrix": "false",
        "entrez_id_count": "0",
        "entrez_id": "",
        "ensembl_gene_id_count": "0",
        "ensembl_gene_id": "",
        "gencode_v49_gene_id_versioned": "",
        "gencode_v49_gene_name": "",
        "gencode_v49_gene_type": "",
        "gencode_v49_contig": "",
    }

    def finish(state: str) -> dict[str, str]:
        record["mapping_state"] = state
        record["contract_state"] = MAPPING_STATE_TO_CONTRACT_STATE[state]
        record["eligible_for_gene_matrix"] = "true" if state == "one_to_one_gencode_v49" else "false"
        return record

    if feature in annotation_absent:
        return finish("platform_feature_absent_from_annotation_package")
    entrez_ids = sorted(probe_entrez.get(feature, set()))
    record["entrez_id_count"] = str(len(entrez_ids))
    if not entrez_ids:
        return finish("no_entrez_id_in_annotation_package")
    if len(entrez_ids) > 1:
        return finish("multiple_entrez_ids")
    entrez = entrez_ids[0]
    record["entrez_id"] = entrez
    if entrez in entrez_absent:
        return finish("entrez_absent_from_org_hs_eg_db")
    ensembl_ids = sorted(entrez_ensembl.get(entrez, set()))
    record["ensembl_gene_id_count"] = str(len(ensembl_ids))
    if not ensembl_ids:
        return finish("entrez_without_ensembl_gene")
    if len(ensembl_ids) > 1:
        return finish("multiple_ensembl_genes")
    ensembl = ensembl_ids[0]
    record["ensembl_gene_id"] = ensembl
    gene = gencode.get(ensembl)
    if gene is None:
        return finish("ensembl_gene_absent_from_gencode_v49")
    record["gencode_v49_gene_id_versioned"] = gene["versioned_id"]
    record["gencode_v49_gene_name"] = gene["gene_name"]
    record["gencode_v49_gene_type"] = gene["gene_type"]
    record["gencode_v49_contig"] = gene["contig"]
    return finish("one_to_one_gencode_v49")


def build_crosswalk(
    *,
    platform_axis: Path,
    probe_entrez_path: Path,
    entrez_ensembl_path: Path,
    annotation_absent_path: Path,
    entrez_absent_path: Path,
    gencode_gtf: Path,
    annotation_versions_path: Path,
    output: Path,
) -> dict[str, object]:
    features = read_platform_axis(platform_axis)
    probe_entrez = read_relation(probe_entrez_path, "platform_feature_id", "entrez_id")
    entrez_ensembl = read_relation(entrez_ensembl_path, "entrez_id", "ensembl_gene_id")
    annotation_absent = read_single_column(annotation_absent_path, "platform_feature_id")
    entrez_absent = read_single_column(entrez_absent_path, "entrez_id")
    gencode = parse_gencode_genes(gencode_gtf)

    unknown = set(probe_entrez) - set(features)
    if unknown:
        raise CrosswalkError("annotation dump carries a feature outside the frozen axis")
    if annotation_absent & set(probe_entrez):
        raise CrosswalkError("a feature is both absent from and present in the annotation package")

    with annotation_versions_path.open(encoding="utf-8", newline="") as handle:
        versions = {
            row["package"]: row["version"]
            for row in csv.DictReader(handle, delimiter="\t")
        }

    rows = [
        classify_feature(
            feature,
            annotation_absent=annotation_absent,
            probe_entrez=probe_entrez,
            entrez_absent=entrez_absent,
            entrez_ensembl=entrez_ensembl,
            gencode=gencode,
        )
        for feature in features
    ]

    mapping_counts: dict[str, int] = defaultdict(int)
    contract_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        mapping_counts[row["mapping_state"]] += 1
        contract_counts[row["contract_state"]] += 1
    if sum(contract_counts.values()) != EXPECTED_PLATFORM_FEATURES:
        raise CrosswalkError("state assignment lost a platform feature")
    if set(contract_counts) - set(CONTRACT_STATES):
        raise CrosswalkError("state assignment produced an undeclared contract state")

    eligible = [row for row in rows if row["eligible_for_gene_matrix"] == "true"]
    gene_features: dict[str, list[str]] = defaultdict(list)
    gene_entrez: dict[str, set[str]] = defaultdict(set)
    for row in eligible:
        gene_features[row["ensembl_gene_id"]].append(row["platform_feature_id"])
        gene_entrez[row["ensembl_gene_id"]].add(row["entrez_id"])

    pairs = [(row["platform_feature_id"], row["ensembl_gene_id"]) for row in eligible]
    crosswalk_sha256 = digest_pairs(pairs)

    output.mkdir(parents=True, exist_ok=False)
    crosswalk_path = output / "gpl16686_gencode_v49_crosswalk.tsv"
    with crosswalk_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(CROSSWALK_COLUMNS), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    collapse_path = output / "gpl16686_gene_collapse_map.tsv"
    with collapse_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "ensembl_gene_id",
                "gencode_v49_gene_id_versioned",
                "gencode_v49_gene_name",
                "gencode_v49_gene_type",
                "entrez_id_count",
                "platform_feature_count",
                "platform_feature_ids",
            ]
        )
        for gene in sorted(gene_features):
            record = gencode[gene]
            members = sorted(gene_features[gene])
            writer.writerow(
                [
                    gene,
                    record["versioned_id"],
                    record["gene_name"],
                    record["gene_type"],
                    len(gene_entrez[gene]),
                    len(members),
                    ";".join(members),
                ]
            )

    multi_entrez_genes = sum(1 for gene in gene_entrez if len(gene_entrez[gene]) > 1)
    feature_counts = [len(members) for members in gene_features.values()]
    receipt = {
        "schema_version": "masld-bench-gpl16686-gencode-v49-crosswalk-v1",
        "status": "pass_gpl16686_gencode_v49_crosswalk_frozen",
        "platform_id": PLATFORM_ID,
        "cohort_family_id": "antwerp_inserm_shared",
        "series": "GSE83452",
        "annotation_authority_order": [
            "hugene20sttranscriptcluster.db_exact_version",
            "org.Hs.eg.db_exact_version",
            "GENCODE_v49",
        ],
        "annotation_package": ANNOTATION_PACKAGE,
        "package_versions": versions,
        "genome_build": "GRCh38.p14",
        "annotation": "GENCODE_v49",
        "gencode_gtf": str(gencode_gtf),
        "gencode_gtf_sha256": sha256_file(gencode_gtf),
        "gencode_v49_gene_records": len(gencode),
        "platform_axis_sha256": sha256_file(platform_axis),
        "platform_features": len(features),
        "mapping_state_counts": dict(sorted(mapping_counts.items())),
        "contract_state_counts": {
            state: contract_counts.get(state, 0) for state in CONTRACT_STATES
        },
        "ambiguous_mapping_state": "join_unresolved",
        "retired_or_absent_GENCODE_state": "below_qc",
        "eligible_platform_features": len(eligible),
        "one_to_one_genes": len(gene_features),
        "genes_with_multiple_platform_features": sum(1 for count in feature_counts if count > 1),
        "maximum_platform_features_per_gene": max(feature_counts) if feature_counts else 0,
        "genes_reached_by_multiple_entrez_ids": multi_entrez_genes,
        "gene_collapse_rule": "training_frozen_median_of_log_scale_features_after_one_to_one_mapping",
        "thermo_release_36_asset": "unavailable_permission_login_gated",
        "geo_platform_table_rows_not_used_as_axis": 53_981,
        "gene_collapse_is_within_array": True,
        "crosswalk_sha256": crosswalk_sha256,
        "coordinate_or_GB_ACC_mapping_inference_used": False,
        "outcome_selected_probe_mapping": False,
        "depends_on_CEL_intensity": False,
        "depends_on_held_cohort_arrays": False,
        "depends_on_labels": False,
        "CEL_expression_values_read": False,
        "GEO_series_matrix_read": False,
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }
    (output / "crosswalk_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform-axis", type=Path, required=True)
    parser.add_argument("--probe-entrez", type=Path, required=True)
    parser.add_argument("--entrez-ensembl", type=Path, required=True)
    parser.add_argument("--annotation-absent", type=Path, required=True)
    parser.add_argument("--entrez-absent", type=Path, required=True)
    parser.add_argument("--gencode-gtf", type=Path, required=True)
    parser.add_argument("--annotation-versions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = build_crosswalk(
        platform_axis=arguments.platform_axis,
        probe_entrez_path=arguments.probe_entrez,
        entrez_ensembl_path=arguments.entrez_ensembl,
        annotation_absent_path=arguments.annotation_absent,
        entrez_absent_path=arguments.entrez_absent,
        gencode_gtf=arguments.gencode_gtf,
        annotation_versions_path=arguments.annotation_versions,
        output=arguments.output,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
