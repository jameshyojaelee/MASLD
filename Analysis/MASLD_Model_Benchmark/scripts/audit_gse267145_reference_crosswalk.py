#!/usr/bin/env python3
"""Crosswalk GSE267145 RNA genes and H3K27ac regions to project references."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import urllib.request


class ReferenceCrosswalkError(RuntimeError):
    """Raised when a source or project reference differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path, maximum_bytes: int = 100_000_000) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "masld-bench/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, destination.open(
        "xb"
    ) as handle:
        digest = sha256()
        size = 0
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > maximum_bytes:
                raise ReferenceCrosswalkError("reference download exceeds size ceiling")
            handle.write(block)
            digest.update(block)
        headers = {
            key.lower(): value
            for key, value in response.headers.items()
            if key.lower() in {"content-length", "last-modified", "etag"}
        }
    return {"url": url, "bytes": size, "sha256": digest.hexdigest(), "headers": headers}


def parse_attributes(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in value.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, observed = item.partition(" ")
        if not separator:
            raise ReferenceCrosswalkError("GTF attribute differs")
        result[key] = observed.strip().strip('"')
    return result


def parse_gtf_genes(path: Path) -> dict[str, dict[str, str]]:
    genes: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ReferenceCrosswalkError(f"GTF width differs at {line_number}")
            if fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            raw_gene = attributes.get("gene_id")
            if raw_gene is None:
                raise ReferenceCrosswalkError("GTF gene lacks gene_id")
            stable = raw_gene.split(".", 1)[0]
            record = {
                "stable_id": stable,
                "versioned_id": raw_gene,
                "contig": fields[0],
                "start_1based": fields[3],
                "end_1based": fields[4],
                "strand": fields[6],
                "gene_name": attributes.get("gene_name", ""),
                "gene_type": attributes.get(
                    "gene_type", attributes.get("gene_biotype", "")
                ),
            }
            if stable in genes and genes[stable] != record:
                raise ReferenceCrosswalkError(f"GTF stable ID is duplicated: {stable}")
            genes[stable] = record
    if not genes:
        raise ReferenceCrosswalkError("GTF has no genes")
    return genes


def read_matrix_genes(path: Path) -> list[str]:
    genes: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "ensembl_gene_id":
            raise ReferenceCrosswalkError("RNA matrix header differs")
        for row in reader:
            if len(row) != len(header):
                raise ReferenceCrosswalkError("RNA matrix width differs")
            genes.append(row[0])
    if len(genes) != 43_285 or len(set(genes)) != len(genes):
        raise ReferenceCrosswalkError("RNA matrix feature axis differs")
    return genes


def build_gene_crosswalk(
    matrix_genes: list[str],
    ensembl98: dict[str, dict[str, str]],
    gencode49: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], Counter[str]]:
    rows: list[dict[str, str]] = []
    states: Counter[str] = Counter()
    for matrix_gene in matrix_genes:
        stable = matrix_gene.split(".", 1)[0]
        old = ensembl98.get(stable)
        current = gencode49.get(stable)
        if old is not None and current is not None:
            state = "stable_id_exact_v98_and_v49"
            allowed = "true"
        elif old is not None:
            state = "ensembl98_only_retired_or_absent_v49"
            allowed = "false"
        elif current is not None:
            state = "gencode49_only_not_in_frozen_v98"
            allowed = "false"
        else:
            state = "absent_from_both_frozen_annotations"
            allowed = "false"
        states[state] += 1
        rows.append(
            {
                "matrix_gene_id": matrix_gene,
                "stable_gene_id": stable,
                "mapping_state": state,
                "allowed_project_input": allowed,
                "ensembl98_gene_name": "" if old is None else old["gene_name"],
                "ensembl98_gene_type": "" if old is None else old["gene_type"],
                "gencode49_gene_id": "" if current is None else current["versioned_id"],
                "gencode49_gene_name": "" if current is None else current["gene_name"],
                "gencode49_gene_type": "" if current is None else current["gene_type"],
            }
        )
    return rows, states


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            contig, size = line.rstrip("\n").split("\t")
            sizes[contig] = int(size)
    if not sizes:
        raise ReferenceCrosswalkError("chromosome sizes are empty")
    return sizes


def audit_h3_regions(path: Path, chrom_sizes: dict[str, int]) -> dict:
    pattern = re.compile(r'"?(chr(?:[0-9]+|X|Y|M)):(\d+)-(\d+)"?')
    regions: set[tuple[str, int, int]] = set()
    contigs: Counter[str] = Counter()
    minimum_width = math.inf
    maximum_width = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline()
        if not header:
            raise ReferenceCrosswalkError("H3K27ac matrix is empty")
        for line_number, line in enumerate(handle, start=2):
            token = line.split(maxsplit=1)[0]
            match = pattern.fullmatch(token)
            if match is None:
                raise ReferenceCrosswalkError(f"H3K27ac region differs at {line_number}")
            contig, start_text, end_text = match.groups()
            start, end = int(start_text), int(end_text)
            if contig not in chrom_sizes or start < 0 or end <= start or end > chrom_sizes[contig]:
                raise ReferenceCrosswalkError("H3K27ac region is outside project primary contigs")
            region = (contig, start, end)
            if region in regions:
                raise ReferenceCrosswalkError("H3K27ac region is duplicated")
            regions.add(region)
            contigs[contig] += 1
            minimum_width = min(minimum_width, end - start)
            maximum_width = max(maximum_width, end - start)
    if len(regions) != 96_460:
        raise ReferenceCrosswalkError("H3K27ac region count differs")
    return {
        "regions": len(regions),
        "contig_counts": dict(sorted(contigs.items())),
        "minimum_interval_width": minimum_width,
        "maximum_interval_width": maximum_width,
        "all_regions_on_GRCh38p14_primary_contigs": True,
        "zero_based_half_open_bounds_compatible": True,
        "one_based_inclusive_bounds_compatible": True,
        "coordinate_base_resolved": False,
        "sequence_extraction_allowed": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-audit", type=Path, required=True)
    parser.add_argument("--matrix-artifacts-sha256", required=True)
    parser.add_argument("--sequence-reference", type=Path, required=True)
    parser.add_argument("--sequence-artifacts-sha256", required=True)
    parser.add_argument("--gencode49-gtf", type=Path, required=True)
    parser.add_argument("--gencode49-gtf-sha256", required=True)
    parser.add_argument("--ensembl98-url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for root, expected in (
        (args.matrix_audit, args.matrix_artifacts_sha256),
        (args.sequence_reference, args.sequence_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise ReferenceCrosswalkError("input ARTIFACTS SHA-256 differs")
    if sha256_file(args.gencode49_gtf) != args.gencode49_gtf_sha256:
        raise ReferenceCrosswalkError("GENCODE v49 GTF SHA-256 differs")
    args.output.mkdir(parents=True, exist_ok=False)
    ensembl_path = args.output / "Homo_sapiens.GRCh38.98.gtf.gz"
    receipt = download(args.ensembl98_url, ensembl_path)
    ensembl98 = parse_gtf_genes(ensembl_path)
    gencode49 = parse_gtf_genes(args.gencode49_gtf)
    genes = read_matrix_genes(args.matrix_audit / "raw" / "GSE269412_RNA.txt.gz")
    crosswalk, states = build_gene_crosswalk(genes, ensembl98, gencode49)
    fields = tuple(crosswalk[0])
    with (args.output / "rna_gene_crosswalk.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(crosswalk)
    h3 = audit_h3_regions(
        args.matrix_audit / "raw" / "GSE267119_H3K27ac.txt.gz",
        read_chrom_sizes(args.sequence_reference / "validation" / "source.chrom.sizes"),
    )
    exact = states["stable_id_exact_v98_and_v49"]
    summary = {
        "schema_version": "masld-bench-gse267145-reference-crosswalk-v1",
        "status": "pass_with_coordinate_base_blocker",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "rna_matrix_genes": len(genes),
        "ensembl98_genes": len(ensembl98),
        "gencode49_genes": len(gencode49),
        "gene_mapping_states": dict(sorted(states.items())),
        "exact_gene_fraction": exact / len(genes),
        "rna_project_input_rule": "retain_only_stable_id_exact_v98_and_v49_and_mask_all_other_features",
        "h3k27ac": h3,
        "ensembl98_receipt": receipt,
        "matrix_artifacts_sha256": args.matrix_artifacts_sha256,
        "sequence_artifacts_sha256": args.sequence_artifacts_sha256,
        "gencode49_gtf_sha256": args.gencode49_gtf_sha256,
    }
    (args.output / "crosswalk_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
