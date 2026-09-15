#!/usr/bin/env python3
"""Build and admit a GSE268273 GENCODE-v36-to-v49 RSEM reference."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

from scripts.plan_gse268273_raw_campaign import (
    GENCODE_V36_GTF_BYTES,
    GENCODE_V36_GTF_MD5,
    GENCODE_V49_GTF_SHA256,
    GENOME_SHA256,
    PROHIBITED_REFERENCE_PATHS,
    sha256_file,
)


EXPECTED_PRIMARY_SCAFFOLD_SOURCE_FEATURES = frozenset(
    {
        "ENSG00000271254.6",
        "ENSG00000274847.1",
        "ENSG00000276256.1",
        "ENSG00000277196.4",
        "ENSG00000277400.1",
        "ENSG00000278384.1",
        "ENSG00000278817.1",
    }
)


class GSE268273ReferenceError(RuntimeError):
    """Raised when a source or built RSEM reference fails inclusion."""


def tool_environment() -> dict[str, str]:
    path = os.environ.get("GSE268273_TOOL_PATH", "")
    if not path:
        raise GSE268273ReferenceError("isolated RSEM/STAR tool PATH is absent")
    environment = dict(os.environ)
    environment["PATH"] = path
    library_path = os.environ.get("GSE268273_TOOL_LD_LIBRARY_PATH", "")
    if library_path:
        environment["LD_LIBRARY_PATH"] = library_path
    else:
        environment.pop("LD_LIBRARY_PATH", None)
    environment.pop("PYTHONPATH", None)
    return environment


def tool_executable(name: str, environment: Mapping[str, str]) -> str:
    executable = shutil.which(name, path=environment["PATH"])
    if executable is None:
        raise GSE268273ReferenceError(f"isolated tool executable is absent: {name}")
    return executable


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_test(path: Path) -> None:
    result = subprocess.run(
        ["gzip", "-t", str(path)],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise GSE268273ReferenceError(
            f"gzip integrity failed: {path}: {result.stderr.strip()}"
        )


def gtf_gene_contract(path: Path) -> tuple[set[str], set[str], Counter[str]]:
    exact: set[str] = set()
    stable: set[str] = set()
    contigs: Counter[str] = Counter()
    pattern = re.compile(r'(?:^|;\s*)gene_id "([^"]+)";')
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise GSE268273ReferenceError("GENCODE GTF row width differs")
            contigs[fields[0]] += 1
            if fields[2] != "gene":
                continue
            match = pattern.search(fields[8])
            if match is None:
                raise GSE268273ReferenceError("GENCODE gene row lacks gene_id")
            gene_id = match.group(1)
            if gene_id in exact:
                raise GSE268273ReferenceError("GENCODE exact gene_id is duplicated")
            exact.add(gene_id)
            stable.add(gene_id.split(".", 1)[0])
    if not exact:
        raise GSE268273ReferenceError("GENCODE GTF contains no genes")
    return exact, stable, contigs


def fasta_contigs(path: Path) -> set[str]:
    observed: set[str] = set()
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            contig = line[1:].split(None, 1)[0]
            if contig in observed:
                raise GSE268273ReferenceError("genome FASTA contig is duplicated")
            observed.add(contig)
    if not observed:
        raise GSE268273ReferenceError("genome FASTA contains no contigs")
    return observed


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE268273ReferenceError(f"TSV lacks a header: {path}")
        return tuple(reader.fieldnames), list(reader)


def admit_reference_inputs(
    *,
    genome: Path,
    source_gtf: Path,
    target_gtf: Path,
    activation_root: Path,
    model_input_root: Path,
) -> dict[str, Any]:
    """Validate the exact source and target gene-axis inputs before indexing."""

    if (
        str(genome) in PROHIBITED_REFERENCE_PATHS
        or str(target_gtf) in PROHIBITED_REFERENCE_PATHS
    ):
        raise GSE268273ReferenceError("a prohibited truncated reference was supplied")
    verify_frozen_tree(activation_root)
    verify_frozen_tree(model_input_root)
    for path in (genome, source_gtf, target_gtf):
        gzip_test(path)
    if sha256_file(genome) != GENOME_SHA256:
        raise GSE268273ReferenceError("GRCh38.p14 compressed SHA-256 differs")
    if (
        source_gtf.stat().st_size != GENCODE_V36_GTF_BYTES
        or md5_file(source_gtf) != GENCODE_V36_GTF_MD5
    ):
        raise GSE268273ReferenceError(
            "GENCODE v36 primary-assembly byte/MD5 contract differs"
        )
    if sha256_file(target_gtf) != GENCODE_V49_GTF_SHA256:
        raise GSE268273ReferenceError("GENCODE v49 compressed SHA-256 differs")
    source_exact, _, source_contigs = gtf_gene_contract(source_gtf)
    _, target_stable, target_contigs = gtf_gene_contract(target_gtf)
    genome_contigs = fasta_contigs(genome)
    if not set(source_contigs) <= genome_contigs:
        missing = sorted(set(source_contigs) - genome_contigs)
        raise GSE268273ReferenceError(
            "GENCODE v36 primary-assembly contig is absent from GRCh38.p14: "
            + ",".join(missing[:10])
        )
    if not set(target_contigs) <= genome_contigs:
        missing = sorted(set(target_contigs) - genome_contigs)
        raise GSE268273ReferenceError(
            "GENCODE v49 annotation contig is absent from GRCh38.p14: "
            + ",".join(missing[:10])
        )
    crosswalk_fields, crosswalk = read_tsv(activation_root / "gene_crosswalk.tsv")
    candidate_fields, candidates = read_tsv(
        model_input_root / "candidate_gene_axis.tsv"
    )
    if (
        not {
            "source_feature_id",
            "mapping_state",
            "gencode_v49_stable_id",
        }
        <= set(crosswalk_fields)
        or not {
            "stable_gene_id",
            "source_feature_count",
            "source_feature_ids",
            "raw_count_aggregation",
        }
        <= set(candidate_fields)
        or len(crosswalk) != 14_149
        or len(candidates) != 14_078
    ):
        raise GSE268273ReferenceError("source-to-target gene contract differs")
    missing_source = sorted(
        row["source_feature_id"]
        for row in crosswalk
        if row["source_feature_id"] not in source_exact
    )
    missing_target = sorted(
        row["stable_gene_id"]
        for row in candidates
        if row["stable_gene_id"] not in target_stable
    )
    duplicate_candidates = [
        row for row in candidates if int(row["source_feature_count"]) > 1
    ]
    wrong_aggregation = sorted(
        row["stable_gene_id"]
        for row in duplicate_candidates
        if row["raw_count_aggregation"]
        != "deterministic_sum_before_normalization"
    )
    mapping_states = Counter(row["mapping_state"] for row in crosswalk)
    missing_scaffold_features = sorted(
        EXPECTED_PRIMARY_SCAFFOLD_SOURCE_FEATURES - source_exact
    )
    mismatches: dict[str, Any] = {}
    if missing_source:
        mismatches["source_ids_absent_from_v36_primary"] = {
            "count": len(missing_source),
            "examples": missing_source[:10],
        }
    if missing_target:
        mismatches["target_ids_absent_from_v49"] = {
            "count": len(missing_target),
            "examples": missing_target[:10],
        }
    if len(duplicate_candidates) != 11:
        mismatches["duplicate_sum_groups"] = len(duplicate_candidates)
    if wrong_aggregation:
        mismatches["wrong_duplicate_aggregation"] = wrong_aggregation[:10]
    expected_states = Counter({"stable_id": 14_089, "unmapped_stable_id": 60})
    if mapping_states != expected_states:
        mismatches["mapping_states"] = dict(sorted(mapping_states.items()))
    if missing_scaffold_features:
        mismatches["required_primary_scaffold_features_absent"] = (
            missing_scaffold_features
        )
    if len(source_exact) != 60_719 or len(source_contigs) != 47:
        mismatches["v36_primary_annotation_census"] = {
            "gene_ids": len(source_exact),
            "contigs": len(source_contigs),
        }
    if mismatches:
        raise GSE268273ReferenceError(
            "GENCODE v36-primary-to-v49 admission census differs: "
            + json.dumps(mismatches, sort_keys=True)
        )
    return {
        "genome_contigs": len(genome_contigs),
        "source_gtf_contigs": len(source_contigs),
        "source_annotation_gene_ids": len(source_exact),
        "source_crosswalk_features_admitted": len(crosswalk),
        "source_primary_scaffold_features_admitted": len(
            EXPECTED_PRIMARY_SCAFFOLD_SOURCE_FEATURES
        ),
        "target_gtf_contigs": len(target_contigs),
        "target_stable_genes_admitted": len(candidates),
        "target_duplicate_sum_groups": len(duplicate_candidates),
        "source_features_masked_absent_v49": mapping_states[
            "unmapped_stable_id"
        ],
        "status": "passed",
    }


def build_reference(
    *,
    genome: Path,
    source_gtf: Path,
    target_gtf: Path,
    activation_root: Path,
    model_input_root: Path,
    work: Path,
    output: Path,
    threads: int,
) -> dict[str, Any]:
    if output.exists() or work.exists():
        raise GSE268273ReferenceError("refusing to overwrite reference output or work")
    admission = admit_reference_inputs(
        genome=genome,
        source_gtf=source_gtf,
        target_gtf=target_gtf,
        activation_root=activation_root,
        model_input_root=model_input_root,
    )
    output.mkdir(parents=True)
    work.mkdir(parents=True)
    source_root = output / "source_evidence"
    reference_root = output / "reference"
    source_root.mkdir()
    reference_root.mkdir()
    shutil.copyfile(source_gtf, source_root / source_gtf.name)
    plain_genome = work / "GRCh38.p14.genome.fa"
    plain_source_gtf = work / "gencode.v36.primary_assembly.annotation.gtf"
    with gzip.open(genome, "rb") as source, plain_genome.open("xb") as target:
        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    with gzip.open(source_gtf, "rb") as source, plain_source_gtf.open("xb") as target:
        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    tools = tool_environment()
    star_executable = tool_executable("STAR", tools)
    rsem_prepare = tool_executable("rsem-prepare-reference", tools)
    rsem_calculate = tool_executable("rsem-calculate-expression", tools)
    prefix = reference_root / "gse268273_gencode_v36"
    command = [
        rsem_prepare,
        "--gtf",
        str(plain_source_gtf),
        "--star",
        "--star-path",
        str(Path(star_executable).parent),
        "--star-sjdboverhang",
        "100",
        "--num-threads",
        str(threads),
        str(plain_genome),
        str(prefix),
    ]
    with (output / "rsem_prepare.stdout.txt").open("x", encoding="utf-8") as stdout, (
        output / "rsem_prepare.stderr.txt"
    ).open("x", encoding="utf-8") as stderr:
        result = subprocess.run(
            command,
            check=False,
            stdout=stdout,
            stderr=stderr,
            text=True,
            env=tools,
        )
    if result.returncode != 0:
        raise GSE268273ReferenceError("rsem-prepare-reference failed")
    required_core = [
        prefix.with_suffix(suffix)
        for suffix in (".grp", ".ti", ".seq", ".chrlist", ".transcripts.fa")
    ]
    required_star = [reference_root / value for value in ("Genome", "SA", "SAindex")]
    if any(not path.is_file() or path.stat().st_size == 0 for path in required_core + required_star):
        raise GSE268273ReferenceError("built RSEM/STAR reference is incomplete")
    plain_genome.unlink()
    plain_source_gtf.unlink()
    work.rmdir()
    version_outputs = {}
    for name, command_value in (
        ("rsem", [rsem_calculate, "--version"]),
        ("star", [star_executable, "--version"]),
    ):
        version_result = subprocess.run(
            command_value,
            check=False,
            capture_output=True,
            text=True,
            env=tools,
        )
        if version_result.returncode != 0:
            raise GSE268273ReferenceError(f"{name} version probe failed")
        version = "\n".join(
            value.strip()
            for value in (version_result.stdout, version_result.stderr)
            if value.strip()
        )
        if not version:
            raise GSE268273ReferenceError(f"{name} version probe returned no text")
        version_outputs[name] = version
    if (
        "RSEM v1.3.1" not in version_outputs["rsem"]
        or "2.7.10b" not in version_outputs["star"]
    ):
        raise GSE268273ReferenceError("observed RSEM or STAR binary version differs")
    receipt = {
        "schema_version": "masld-bench-gse268273-rsem-reference-build-v2",
        "status": "passed",
        "genome": "GRCh38.p14",
        "genome_compressed_sha256": GENOME_SHA256,
        "source_annotation": "GENCODE_v36_primary_assembly",
        "reference_contract_revision": (
            "gse268273_reference_primary_assembly_revision_v2"
        ),
        "superseded_failed_reference_job": 21082744,
        "source_annotation_md5": GENCODE_V36_GTF_MD5,
        "target_annotation": "GENCODE_v49",
        "target_annotation_compressed_sha256": GENCODE_V49_GTF_SHA256,
        **admission,
        "source_features_admitted": admission[
            "source_crosswalk_features_admitted"
        ],
        "rsem_version": version_outputs["rsem"],
        "rsem_runtime_module_label": "RSEM/1.3.3-foss-2022a",
        "rsem_module_label_differs_from_binary_reported_version": True,
        "star_version": version_outputs["star"],
        "sjdb_overhang": 100,
        "alignment_parameters": "RSEM_integrated_STAR_defaults",
        "source_exact_ENCODE_alignment_invocation_available": False,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "reference_prefix": "reference/gse268273_gencode_v36",
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse268273_rsem_star_reference",
            "source_features": 14_149,
            "target_stable_genes": 14_078,
            "status": "passed",
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--source-gtf", type=Path, required=True)
    parser.add_argument("--target-gtf", type=Path, required=True)
    parser.add_argument("--activation-root", type=Path, required=True)
    parser.add_argument("--model-input-root", type=Path, required=True)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--admission-only", action="store_true")
    args = parser.parse_args()
    if args.admission_only:
        receipt = admit_reference_inputs(
            genome=args.genome,
            source_gtf=args.source_gtf,
            target_gtf=args.target_gtf,
            activation_root=args.activation_root,
            model_input_root=args.model_input_root,
        )
        print(json.dumps(receipt, sort_keys=True))
        return 0
    if args.work is None or args.output is None or args.threads is None:
        parser.error("--work, --output, and --threads are required for a build")
    receipt = build_reference(
        genome=args.genome,
        source_gtf=args.source_gtf,
        target_gtf=args.target_gtf,
        activation_root=args.activation_root,
        model_input_root=args.model_input_root,
        work=args.work,
        output=args.output,
        threads=args.threads,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
