#!/usr/bin/env python3
"""Map the frozen GSE256398 10x gene axis to GENCODE v49 without outcomes."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path

import h5py
import numpy as np


EXPECTED_FEATURES = 36_601


class GSE256398CrosswalkError(RuntimeError):
    """Raised when a source or reference differs from the frozen requirements."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def digest_values(values: list[str]) -> str:
    value = sha256()
    for item in values:
        value.update(item.encode("utf-8"))
        value.update(b"\n")
    return value.hexdigest()


def decode_many(values: np.ndarray) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


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
            raise GSE256398CrosswalkError("GTF attribute differs")
        output[key] = observed.strip().strip('"')
    return output


def parse_gencode_genes(path: Path) -> dict[str, dict[str, str]]:
    genes: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise GSE256398CrosswalkError(f"GTF width differs at line {line_number}")
            if fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            versioned = attributes.get("gene_id")
            if not versioned:
                raise GSE256398CrosswalkError("GENCODE gene lacks gene_id")
            stable = stable_gene_id(versioned)
            record = {
                "stable_id": stable,
                "versioned_id": versioned,
                "gene_name": attributes.get("gene_name", ""),
                "gene_type": attributes.get("gene_type", attributes.get("gene_biotype", "")),
                "contig": fields[0],
                "start_1based": fields[3],
                "end_1based": fields[4],
                "strand": fields[6],
            }
            if stable in genes and genes[stable] != record:
                raise GSE256398CrosswalkError(f"GENCODE stable gene ID is duplicated: {stable}")
            genes[stable] = record
    if not genes:
        raise GSE256398CrosswalkError("GENCODE GTF contains no genes")
    return genes


def read_source_features(path: Path) -> tuple[list[str], list[str], list[str], list[str]]:
    with h5py.File(path, "r") as handle:
        try:
            features = handle["matrix/features"]
            ids = decode_many(features["id"][...])
            names = decode_many(features["name"][...])
            feature_types = decode_many(features["feature_type"][...])
            genomes = decode_many(features["genome"][...])
            shape = tuple(int(value) for value in handle["matrix/shape"][...])
        except KeyError as error:
            raise GSE256398CrosswalkError("10x feature schema differs") from error
    if (
        shape[0] != EXPECTED_FEATURES
        or len(ids) != EXPECTED_FEATURES
        or len(names) != EXPECTED_FEATURES
        or len(feature_types) != EXPECTED_FEATURES
        or len(genomes) != EXPECTED_FEATURES
        or len(set(ids)) != len(ids)
        or set(feature_types) != {"Gene Expression"}
        or set(genomes) != {"GRCh38"}
    ):
        raise GSE256398CrosswalkError("GSE256398 frozen feature axis differs")
    return ids, names, feature_types, genomes


def build_crosswalk(
    source_ids: list[str], source_names: list[str], gencode: dict[str, dict[str, str]]
) -> tuple[list[dict[str, object]], Counter[str]]:
    if len(source_ids) != len(source_names) or len(set(source_ids)) != len(source_ids):
        raise GSE256398CrosswalkError("source feature IDs are not a unique aligned axis")
    stable_ids = [stable_gene_id(value) for value in source_ids]
    if len(set(stable_ids)) != len(stable_ids):
        raise GSE256398CrosswalkError("source stable gene IDs are duplicated")
    rows: list[dict[str, object]] = []
    states: Counter[str] = Counter()
    for index, (source_id, stable, source_name) in enumerate(
        zip(source_ids, stable_ids, source_names, strict=True)
    ):
        current = gencode.get(stable)
        if current is None:
            state = "absent_or_retired_from_gencode_v49"
            allowed = "false"
        else:
            state = "stable_id_exact_gencode_v49"
            allowed = "true"
        states[state] += 1
        rows.append(
            {
                "source_feature_index": index,
                "source_feature_id": source_id,
                "source_stable_gene_id": stable,
                "source_feature_name": source_name,
                "mapping_state": state,
                "allowed_project_input": allowed,
                "source_and_v49_name_equal": (
                    "" if current is None else str(source_name == current["gene_name"]).lower()
                ),
                "gencode49_gene_id": "" if current is None else current["versioned_id"],
                "gencode49_gene_name": "" if current is None else current["gene_name"],
                "gencode49_gene_type": "" if current is None else current["gene_type"],
                "gencode49_contig": "" if current is None else current["contig"],
                "gencode49_start_1based": "" if current is None else current["start_1based"],
                "gencode49_end_1based": "" if current is None else current["end_1based"],
                "gencode49_strand": "" if current is None else current["strand"],
            }
        )
    return rows, states


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--gencode49-gtf", type=Path, required=True)
    parser.add_argument("--gencode49-gtf-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sha256_file(args.source / "ARTIFACTS.json") != args.source_artifacts_sha256:
        raise GSE256398CrosswalkError("source ARTIFACTS SHA-256 differs")
    if sha256_file(args.gencode49_gtf) != args.gencode49_gtf_sha256:
        raise GSE256398CrosswalkError("GENCODE v49 GTF SHA-256 differs")
    source_receipt = json.loads((args.source / "audit/receipt.json").read_text(encoding="utf-8"))
    h5_paths = sorted((args.source / "h5").glob("*.h5"))
    if len(h5_paths) != 26:
        raise GSE256398CrosswalkError("source donor H5 census differs")
    ids, names, _, _ = read_source_features(h5_paths[0])
    if (
        digest_values(ids) != source_receipt["feature_ids_sha256"]
        or digest_values(names) != source_receipt["feature_names_sha256"]
    ):
        raise GSE256398CrosswalkError("source feature roster hash differs")
    gencode = parse_gencode_genes(args.gencode49_gtf)
    rows, states = build_crosswalk(ids, names, gencode)
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "gene_crosswalk.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    exact = states["stable_id_exact_gencode_v49"]
    name_equal = sum(row["source_and_v49_name_equal"] == "true" for row in rows)
    summary = {
        "schema_version": "masld-bench-gse256398-gencode49-crosswalk-v1",
        "status": "pass_with_unmapped_features_masked",
        "dataset_id": "gse256398",
        "source_features": len(rows),
        "gencode49_genes": len(gencode),
        "mapping_states": dict(sorted(states.items())),
        "exact_stable_id_fraction": exact / len(rows),
        "exact_stable_id_and_name_count": name_equal,
        "source_feature_ids_sha256": source_receipt["feature_ids_sha256"],
        "source_feature_names_sha256": source_receipt["feature_names_sha256"],
        "source_artifacts_sha256": args.source_artifacts_sha256,
        "gencode49_gtf_sha256": args.gencode49_gtf_sha256,
        "project_input_rule": "retain_only_stable_id_exact_gencode_v49_and_explicitly_mask_all_other_features",
        "reference_mapping_is_outcome_independent": True,
        "normalization_or_feature_selection_fitted": False,
        "model_training_activated": False,
        "sealed_or_final_outcomes_read": False,
        "exact_cell_ranger_reference": "unresolved",
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
