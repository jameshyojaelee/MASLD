#!/usr/bin/env python3
"""Build a label-free, no-fit GSE105127 RNA/RRBS activation fixture."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


EXPECTED_GEO_ARTIFACTS_SHA256 = (
    "bac433336705b74dcb2a6ef6b6594af10fef362769d564adfea9788fe546e88e"
)
EXPECTED_JOIN_ARTIFACTS_SHA256 = (
    "574eeee1119054a2d7b55066804accf29fb45f23b591dd26f9d61d907b27a75f"
)
EXPECTED_RRBS_ARTIFACTS_SHA256 = (
    "5e3f698a5172d612b38bfe0a78b57f7199c10031ae9740190ead1835ef90e216"
)
EXPECTED_SOFT_SHA256 = (
    "9b9c0b53712fcf18d3b6f48cd3902d4cf2c13541634c5d1bc86ea6cb17922b0e"
)
EXPECTED_SAMPLE_TABLE_SHA256 = (
    "28d5442df8c5e7eb52b4b517398f4b08f9dd6259edc12da14f13f2d2d7fe2a6f"
)
EXPECTED_RRBS_FILES = 57
EXPECTED_RNA_FILES = 57
EXPECTED_CYTOSINE_ROWS = 800_713_381
EXPECTED_RRBS_COMPRESSED_BYTES = 8_192_644_943
EXPECTED_RNA_FASTQ_BYTES = 79_576_393_854
EXPECTED_RNA_READS = 1_140_779_295
EXPECTED_RNA_BASES = 86_699_226_420
TITLE = re.compile(r"^(?P<participant>[0-9]+)_(?P<zone>CV|IZ|PP)_(?P<assay>RNA|RRBS)$")
BIOSAMPLE = re.compile(r"/biosample/(?P<value>SAMN[0-9]+)/?$")
EXPERIMENT = re.compile(r"[?&]term=(?P<value>SRX[0-9]+)$")
MD5 = re.compile(r"^[0-9a-f]{32}$")
FORBIDDEN_MODEL_FIELDS = {
    "phenotype",
    "phenotype_name",
    "label",
    "outcome",
    "disease",
    "fibrosis",
    "nas",
    "sex",
    "age",
    "bmi",
    "outer_fold",
}


class GSE105127ActivationError(RuntimeError):
    """Raised when a GSE105127 no-fit activation invariant differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def row_id(participant: str, zone: str) -> str:
    return "g105_" + sha256(
        b"masld-bench-gse105127-row-v1\0"
        + participant.encode("ascii")
        + b"\0"
        + zone.encode("ascii")
    ).hexdigest()[:20]


def participant_group_id(participant: str) -> str:
    return "g105p_" + sha256(
        b"masld-bench-gse105127-participant-v1\0" + participant.encode("ascii")
    ).hexdigest()[:20]


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE105127ActivationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise GSE105127ActivationError(f"refusing to write empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = tuple(rows[0])
    if FORBIDDEN_MODEL_FIELDS & set(fields):
        raise GSE105127ActivationError("label-bearing field entered model-facing fixture")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_soft_samples(payload: bytes) -> list[dict[str, str]]:
    try:
        text = gzip.decompress(payload).decode("utf-8", errors="strict")
    except (OSError, UnicodeDecodeError) as error:
        raise GSE105127ActivationError("GEO SOFT is not strict UTF-8 gzip") from error
    samples: list[dict[str, str]] = []
    current: dict[str, Any] | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            if current is not None:
                samples.append(_finish_soft_sample(current))
            current = {
                "sample_accession": line.split(" = ", 1)[1],
                "relations": [],
                "supplementary": [],
            }
        elif current is not None and line.startswith("!Sample_") and " = " in line:
            key, value = line.split(" = ", 1)
            if key == "!Sample_relation":
                current["relations"].append(value)
            elif key.startswith("!Sample_supplementary_file_"):
                current["supplementary"].append(value)
            else:
                current[key] = value
    if current is not None:
        samples.append(_finish_soft_sample(current))
    if len(samples) != 114 or len({row["sample_accession"] for row in samples}) != 114:
        raise GSE105127ActivationError("GEO SOFT sample census differs")
    return samples


def _finish_soft_sample(raw: Mapping[str, Any]) -> dict[str, str]:
    title = str(raw.get("!Sample_title", ""))
    match = TITLE.fullmatch(title)
    if match is None:
        raise GSE105127ActivationError("GEO sample title violates participant-zone-assay contract")
    relations = [str(value) for value in raw.get("relations", [])]
    biosamples = [match["value"] for value in relations if (match := BIOSAMPLE.search(value))]
    experiments = [match["value"] for value in relations if (match := EXPERIMENT.search(value))]
    if len(biosamples) != 1 or len(experiments) != 1:
        raise GSE105127ActivationError("GEO BioSample/SRA-experiment relation differs")
    assay = match_title["assay"] if (match_title := TITLE.fullmatch(title)) else ""
    expected = {
        "RNA": ("RNA-Seq", "transcriptomic", "cDNA", "total RNA"),
        "RRBS": ("Bisulfite-Seq", "genomic", "Reduced Representation", "genomic DNA"),
    }[assay]
    observed = tuple(
        str(raw.get(key, ""))
        for key in (
            "!Sample_library_strategy",
            "!Sample_library_source",
            "!Sample_library_selection",
            "!Sample_molecule_ch1",
        )
    )
    if observed != expected or raw.get("!Sample_instrument_model") != "Illumina HiSeq 2500":
        raise GSE105127ActivationError("GEO assay or instrument semantics differ")
    supplementary = [str(value) for value in raw.get("supplementary", [])]
    expected_suffix = "_RNA.coverage.bw" if assay == "RNA" else ".CG.bed.gz"
    primary = [value for value in supplementary if value.endswith(expected_suffix)]
    if len(primary) != 1:
        raise GSE105127ActivationError("GEO primary processed-file identity differs")
    return {
        "participant_id": match_title["participant"],
        "zone": match_title["zone"],
        "assay": assay,
        "sample_accession": str(raw["sample_accession"]),
        "biosample_accession": biosamples[0],
        "experiment_accession": experiments[0],
        "library_strategy": expected[0],
        "primary_processed_file": primary[0],
    }


def parse_ena_runinfo(payload: bytes) -> list[dict[str, str]]:
    try:
        text = payload.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise GSE105127ActivationError("ENA run manifest is not strict UTF-8") from error
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    required = {
        "run_accession",
        "study_accession",
        "sample_accession",
        "experiment_accession",
        "scientific_name",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
    }
    if not rows or not required <= set(rows[0]):
        raise GSE105127ActivationError("ENA run manifest schema differs")
    if (
        len(rows) != 114
        or len({row["run_accession"] for row in rows}) != 114
        or len({row["experiment_accession"] for row in rows}) != 114
        or len({row["sample_accession"] for row in rows}) != 114
    ):
        raise GSE105127ActivationError("ENA run/experiment/BioSample census differs")
    strategies = Counter(row["library_strategy"] for row in rows)
    if strategies != Counter(
        {"RNA-Seq": EXPECTED_RNA_FILES, "Bisulfite-Seq": EXPECTED_RRBS_FILES}
    ):
        raise GSE105127ActivationError("ENA assay census differs")
    for row in rows:
        expected_source, expected_selection = {
            "RNA-Seq": ("TRANSCRIPTOMIC", "cDNA"),
            "Bisulfite-Seq": ("GENOMIC", "Reduced Representation"),
        }[row["library_strategy"]]
        try:
            reads = int(row["read_count"])
            bases = int(row["base_count"])
            size = int(row["fastq_bytes"])
        except ValueError as error:
            raise GSE105127ActivationError("ENA numerical field differs") from error
        path = PurePosixPath(row["fastq_ftp"])
        if (
            row["study_accession"] != "PRJNA414905"
            or row["scientific_name"] != "Homo sapiens"
            or row["library_source"] != expected_source
            or row["library_selection"] != expected_selection
            or row["library_layout"] != "SINGLE"
            or row["instrument_platform"] != "ILLUMINA"
            or row["instrument_model"] != "Illumina HiSeq 2500"
            or reads <= 0
            or bases <= 0
            or size <= 0
            or not row["fastq_ftp"].startswith("ftp.sra.ebi.ac.uk/")
            or path.name != f"{row['run_accession']}.fastq.gz"
            or MD5.fullmatch(row["fastq_md5"].lower()) is None
        ):
            raise GSE105127ActivationError("ENA assay-native run semantics differ")
    rna = [row for row in rows if row["library_strategy"] == "RNA-Seq"]
    if (
        sum(int(row["read_count"]) for row in rna) != EXPECTED_RNA_READS
        or sum(int(row["base_count"]) for row in rna) != EXPECTED_RNA_BASES
        or sum(int(row["fastq_bytes"]) for row in rna) != EXPECTED_RNA_FASTQ_BYTES
        or any(int(row["base_count"]) != 76 * int(row["read_count"]) for row in rna)
    ):
        raise GSE105127ActivationError("RNA FASTQ byte/read/base census differs")
    return rows


def canonical_cpg_interval(start0: int, end0: int, strand: str) -> tuple[int, int]:
    """Return the two-base reference CpG interval for one BisSNP cytosine row."""
    if start0 < 0 or end0 != start0 + 1 or strand not in {"+", "-"}:
        raise GSE105127ActivationError("BisSNP coordinate or strand differs")
    canonical_start = start0 if strand == "+" else start0 - 1
    if canonical_start < 0:
        raise GSE105127ActivationError("negative-strand BisSNP row underflows contig")
    return canonical_start, canonical_start + 2


def methylated_read_count(percent: float, coverage: int) -> int:
    """Recover the integer BisSNP methylated count from a 2-decimal percentage."""
    if not 0 <= percent <= 100 or coverage < 1:
        raise GSE105127ActivationError("BisSNP percent or coverage differs")
    implied = percent * coverage / 100.0
    observed = int(implied + 0.5)
    if abs(implied - observed) > coverage * 0.00005 + 0.00001:
        raise GSE105127ActivationError("BisSNP percentage cannot recover an integer count")
    return observed


def collapse_cpg_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Collapse one or two strand-specific BisSNP rows without double counting."""
    if len(rows) not in {1, 2}:
        raise GSE105127ActivationError("canonical CpG has more than two source rows")
    intervals = {
        canonical_cpg_interval(int(row["start0"]), int(row["end0"]), str(row["strand"]))
        for row in rows
    }
    strands = [str(row["strand"]) for row in rows]
    if len(intervals) != 1 or len(strands) != len(set(strands)):
        raise GSE105127ActivationError("BisSNP strand rows do not identify one CpG")
    coverage = sum(int(row["coverage"]) for row in rows)
    methylated = sum(
        methylated_read_count(float(row["percent"]), int(row["coverage"]))
        for row in rows
    )
    if methylated > coverage:
        raise GSE105127ActivationError("collapsed methylated count exceeds coverage")
    start0, end0 = next(iter(intervals))
    return {
        "start0": start0,
        "end0": end0,
        "methylated_reads": methylated,
        "coverage": coverage,
        "beta": methylated / coverage,
        "strand_state": "both_strands" if len(rows) == 2 else f"{strands[0]}_only",
        "auxiliary_columns_used": False,
    }


def liftover_state(
    *,
    source_is_cpg: bool,
    source_contig_recognized: bool,
    mapped_intervals: Sequence[tuple[str, int, int]],
    target_primary: bool = False,
    target_is_cpg: bool = False,
    roundtrip_exact: bool = False,
) -> str:
    """Classify one CpG crosswalk disposition in fail-closed order."""
    if not source_contig_recognized:
        return "source_contig_unrecognized"
    if not source_is_cpg:
        return "source_not_cpg"
    if not mapped_intervals:
        return "unmapped_chain"
    if len(mapped_intervals) != 1:
        return "mapped_nonunique"
    _, start0, end0 = mapped_intervals[0]
    if start0 < 0 or end0 - start0 != 2:
        return "mapped_length_changed"
    if not target_primary:
        return "mapped_nonprimary"
    if not target_is_cpg:
        return "target_not_cpg"
    if not roundtrip_exact:
        return "roundtrip_failed"
    return "mapped_unique_cpg"


def _verify_sources(geo: Path, join: Path, rrbs: Path) -> None:
    for root, expected in (
        (geo, EXPECTED_GEO_ARTIFACTS_SHA256),
        (join, EXPECTED_JOIN_ARTIFACTS_SHA256),
        (rrbs, EXPECTED_RRBS_ARTIFACTS_SHA256),
    ):
        verify_frozen_tree(root)
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise GSE105127ActivationError(f"frozen source identity differs: {root}")
    if (
        sha256_file(geo / "raw/GSE105127_family.soft.gz") != EXPECTED_SOFT_SHA256
        or sha256_file(geo / "samples/GSE105127.tsv") != EXPECTED_SAMPLE_TABLE_SHA256
    ):
        raise GSE105127ActivationError("GEO source member identity differs")
    join_audit = json.loads((join / "join_audit.json").read_text(encoding="utf-8"))
    if (
        join_audit.get("status") != "pass"
        or join_audit.get("participants") != 19
        or join_audit.get("assay_records") != 114
        or join_audit.get("pairing") != "adjacent_section"
        or join_audit.get("same_cell_or_same_section_pairing") is not False
        or join_audit.get("model_training_activated") is not False
    ):
        raise GSE105127ActivationError("frozen adjacent-section join authority differs")


def _validate_rrbs_audit(rrbs: Path) -> list[dict[str, str]]:
    summary = json.loads((rrbs / "qc_summary.json").read_text(encoding="utf-8"))
    artifacts = json.loads((rrbs / "ARTIFACTS.json").read_text(encoding="utf-8"))
    artifact_sha = {
        row["path"]: row["sha256"] for row in artifacts.get("artifacts", [])
    }
    _, rows = read_tsv(rrbs / "rrbs_file_qc.tsv")
    if (
        summary.get("status") != "pass"
        or summary.get("rrbs_files") != EXPECTED_RRBS_FILES
        or summary.get("cytosine_rows") != EXPECTED_CYTOSINE_ROWS
        or summary.get("compressed_bytes") != EXPECTED_RRBS_COMPRESSED_BYTES
        or summary.get("hard_qc_failures") != 0
        or len(rows) != EXPECTED_RRBS_FILES
        or len({row["sample_accession"] for row in rows}) != EXPECTED_RRBS_FILES
        or any(row["hard_qc_state"] != "observed" or row["error_total"] != "0" for row in rows)
        or sum(int(row["rows"]) for row in rows) != EXPECTED_CYTOSINE_ROWS
        or sum(int(row["bytes"]) for row in rows) != EXPECTED_RRBS_COMPRESSED_BYTES
    ):
        raise GSE105127ActivationError("frozen exhaustive RRBS audit differs")
    for row in rows:
        path = rrbs / row["local_file"]
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or artifact_sha.get(row["local_file"]) != row["sha256"]
        ):
            raise GSE105127ActivationError("frozen RRBS source file identity differs")
    return rows


def build_activation(
    *,
    geo_root: Path,
    join_root: Path,
    rrbs_root: Path,
    ena_runinfo: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise GSE105127ActivationError(f"refusing to overwrite activation fixture: {output}")
    _verify_sources(geo_root, join_root, rrbs_root)
    rrbs_qc = _validate_rrbs_audit(rrbs_root)
    soft = parse_soft_samples((geo_root / "raw/GSE105127_family.soft.gz").read_bytes())
    ena_payload = ena_runinfo.read_bytes()
    ena = parse_ena_runinfo(ena_payload)
    by_experiment = {row["experiment_accession"]: row for row in ena}
    if len(by_experiment) != 114:
        raise GSE105127ActivationError("ENA experiment axis differs")
    rrbs_by_gsm = {row["sample_accession"]: row for row in rrbs_qc}
    sample_by_key = {(row["participant_id"], row["zone"], row["assay"]): row for row in soft}
    expected_keys = {
        (participant, zone, assay)
        for participant in {row["participant_id"] for row in soft}
        for zone in ("CV", "IZ", "PP")
        for assay in ("RNA", "RRBS")
    }
    if len({row["participant_id"] for row in soft}) != 19 or set(sample_by_key) != expected_keys:
        raise GSE105127ActivationError("participant-zone-assay cross-product differs")

    row_axis: list[dict[str, Any]] = []
    rna_runs: list[dict[str, Any]] = []
    rrbs_sources: list[dict[str, Any]] = []
    for participant, zone in sorted({(row["participant_id"], row["zone"]) for row in soft}):
        rna = sample_by_key[(participant, zone, "RNA")]
        rrbs = sample_by_key[(participant, zone, "RRBS")]
        identifier = row_id(participant, zone)
        rna_run = by_experiment.get(rna["experiment_accession"])
        rrbs_run = by_experiment.get(rrbs["experiment_accession"])
        rrbs_receipt = rrbs_by_gsm.get(rrbs["sample_accession"])
        if (
            rna_run is None
            or rrbs_run is None
            or rrbs_receipt is None
            or rna_run["sample_accession"] != rna["biosample_accession"]
            or rrbs_run["sample_accession"] != rrbs["biosample_accession"]
            or rna_run["library_strategy"] != "RNA-Seq"
            or rrbs_run["library_strategy"] != "Bisulfite-Seq"
            or rrbs_receipt["participant_id"] != participant
            or rrbs_receipt["zone"] != zone
            or rrbs_receipt["pairing"] != "adjacent_section"
        ):
            raise GSE105127ActivationError("GEO/ENA/RRBS sample join differs")
        row_axis.append(
            {
                "row_id": identifier,
                "participant_group_id": participant_group_id(participant),
                "zone": zone,
                "biological_unit": "participant",
                "repeated_measurement_unit": "liver_zone",
                "pairing_topology": "adjacent_section",
                "same_section": False,
                "same_cell": False,
                "same_molecule": False,
                "rna_state": "derivable_not_processed",
                "rrbs_state": "observed_native_hg19_remap_pending",
            }
        )
        rna_runs.append(
            {
                "row_id": identifier,
                "sample_accession": rna["sample_accession"],
                "biosample_accession": rna["biosample_accession"],
                "experiment_accession": rna["experiment_accession"],
                "run_accession": rna_run["run_accession"],
                "library_layout": "SINGLE",
                "read_length": 76,
                "read_count": int(rna_run["read_count"]),
                "base_count": int(rna_run["base_count"]),
                "fastq_url": "https://" + rna_run["fastq_ftp"],
                "fastq_bytes": int(rna_run["fastq_bytes"]),
                "fastq_md5": rna_run["fastq_md5"].lower(),
                "processed_coverage_bigwig_is_expression_matrix": False,
                "model_input_state": "derivable_not_processed",
            }
        )
        rrbs_sources.append(
            {
                "row_id": identifier,
                "sample_accession": rrbs["sample_accession"],
                "biosample_accession": rrbs["biosample_accession"],
                "experiment_accession": rrbs["experiment_accession"],
                "run_accession": rrbs_run["run_accession"],
                "processed_bed_source_url": rrbs_receipt["source_url"],
                "processed_bed_sha256": rrbs_receipt["sha256"],
                "processed_bed_bytes": int(rrbs_receipt["bytes"]),
                "cytosine_rows": int(rrbs_receipt["rows"]),
                "hard_qc_state": "observed",
                "coordinate_state": "native_hg19_remap_pending",
            }
        )

    if (
        len(row_axis) != EXPECTED_RNA_FILES
        or len({row["row_id"] for row in row_axis}) != EXPECTED_RNA_FILES
    ):
        raise GSE105127ActivationError("label-free participant-zone row axis differs")
    if len({row["participant_group_id"] for row in row_axis}) != 19:
        raise GSE105127ActivationError("label-free participant census differs")

    output.mkdir(parents=True)
    write_tsv(output / "model_inputs/row_axis.tsv", row_axis)
    write_tsv(output / "authorities/rna_run_manifest.tsv", rna_runs)
    write_tsv(output / "authorities/rrbs_source_manifest.tsv", rrbs_sources)

    source_semantics = {
        "schema_version": "masld-bench-gse105127-source-semantics-v1",
        "rrbs_files_exhaustively_audited": 57,
        "rrbs_cytosine_rows_exhaustively_audited": EXPECTED_CYTOSINE_ROWS,
        "rrbs_hard_format_failures": 0,
        "bed_style": "BisSNP_BED9_plus_2",
        "coordinate_semantics": "zero_based_half_open_single_reference_cytosine",
        "negative_strand_semantics": "reported_cytosine_is_shifted_plus_one_from_reference_C_of_CpG",
        "canonical_cpg_interval": {
            "plus": "[start0,start0+2)",
            "minus": "[start0-1,start0+1)",
        },
        "methylation_percent_column": 4,
        "coverage_column": 5,
        "strand_column": 6,
        "methylated_read_reconstruction": "round(percent_times_coverage_divided_by_100)_after_two_decimal_tolerance_check",
        "strand_collapse": "sum_reconstructed_methylated_reads_and_column5_coverage_once_per_observed_strand",
        "one_strand_only_sites": "retain_with_observed_coverage_and_explicit_strand_state",
        "columns_10_and_11": "BisSNP_bedDetail_auxiliary_fields_excluded_from_model_measurement_and_never_double_counted",
        "beta": "methylated_reads_divided_by_coverage",
        "missing_site": "not_observed_or_below_qc_never_zero_methylation",
        "likelihood": "coverage_aware_binomial_or_training_only_dispersion_beta_binomial",
        "rna_count_likelihood_for_rrbs_forbidden": True,
        "evidence": [
            "https://www.bioconductor.org/packages/devel/bioc/vignettes/RnBeads/inst/doc/RnBeads.pdf",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3491382/",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC6175862/",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105127",
        ],
    }
    remap = {
        "schema_version": "masld-bench-gse105127-cpg-remap-contract-v1",
        "status": "designed_not_executed",
        "source_reference": {
            "assembly": "1000_Genomes_GRCh37",
            "url": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz",
            "declared_bytes": 892_331_003,
            "sha256": "derivable_not_processed",
            "required_checks": ["gzip_test", "sha256", "contig_inventory", "fai_rebuild_and_compare"],
        },
        "forward_chain": {
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
            "upstream_md5": "35887f73fe5e2231656504d1f6430900",
            "sha256": "derivable_not_processed",
        },
        "reverse_chain": {
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
            "upstream_md5": "ff3031d93792f4cbb86af44055efd903",
            "sha256": "derivable_not_processed",
        },
        "target_reference": {
            "assembly": "GRCh38.p14",
            "path": "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz",
            "sha256": "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215",
            "accepted_contigs": "chr1_to_chr22_chrX_chrY_chrM_only",
        },
        "unit": "unique_two_base_canonical_CpG_interval_after_within_sample_strand_collapse",
        "mapping_tool": "CrossMap_0.7.0_bed_after_runtime_probe",
        "forward_command_shape": "CrossMap bed --chromid l --unmap-file FORWARD_UNMAP FORWARD_CHAIN SOURCE_CPG_BED FORWARD_BED",
        "reverse_command_shape": "CrossMap bed --chromid l --unmap-file REVERSE_UNMAP REVERSE_CHAIN TARGET_CPG_BED REVERSE_BED",
        "record_identity": "opaque_source_CpG_ID_retained_in_column4_and_required_exactly_once_after_each_direction",
        "source_contig_aliases": "1_to_chr1_through_22_to_chr22_X_to_chrX_Y_to_chrY_MT_to_chrM_only",
        "accepted_state": "mapped_unique_cpg",
        "failed_or_unmapped_states": [
            "source_contig_unrecognized",
            "source_not_cpg",
            "strand_pair_inconsistent",
            "unmapped_chain",
            "mapped_nonunique",
            "mapped_length_changed",
            "mapped_nonprimary",
            "target_not_cpg",
            "roundtrip_failed",
        ],
        "acceptance_checks": [
            "source_two_base_reference_sequence_is_CG",
            "exactly_one_forward_mapping",
            "target_interval_length_is_two",
            "target_contig_is_primary",
            "target_two_base_reference_sequence_is_CG",
            "exactly_one_reverse_mapping",
            "reverse_mapping_equals_source_contig_start_end",
        ],
        "split_output_policy": "a_split_or_repeated_source_CpG_ID_is_mapped_nonunique_or_mapped_length_changed_never_accepted",
        "failure_policy": "retain_every_source_CpG_in_crosswalk_with_state_and_never_silently_drop",
        "bed_coordinate_convention": "zero_based_half_open",
    }
    rna_contract = {
        "schema_version": "masld-bench-gse105127-rna-preprocessing-contract-v1",
        "status": "raw_reads_available_quantification_not_run",
        "libraries": 57,
        "participants": 19,
        "zones_per_participant": 3,
        "layout": "single_end",
        "read_length": 76,
        "reads": EXPECTED_RNA_READS,
        "bases": EXPECTED_RNA_BASES,
        "fastq_bytes": EXPECTED_RNA_FASTQ_BYTES,
        "processed_bigwig_role": "coverage_visualization_only_not_gene_expression_input",
        "required_processing": [
            "verify_ENA_MD5_and_gzip_integrity",
            "adapter_and_quality_trim_with_exact_frozen_runtime",
            "quantify_each_RNA_GSM_independently_against_GRCh38p14_GENCODE_v49",
            "retain_continuous_RSEM_expected_counts_effective_lengths_and_TPM_without_rounding",
            "fit_filtering_normalization_batch_adjustment_and_count_to_context_mapping_inside_outer_training_participants_only",
        ],
        "recommended_runtime": "RSEM_1.3.3_with_STAR_2.7.10b_after_runtime_probe",
        "reference_fasta_sha256": "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215",
        "gencode_v49_gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
        "source_pipeline_for_provenance_only": "TrimGalore_0.4.2_q20_STAR_2.5.2a_per_sample_2pass_featureCounts_1.5.0p3_GENCODE_v19_primary_alignments",
        "same_section_pairing_allowed": False,
    }
    topology = {
        "schema_version": "masld-bench-gse105127-topology-v1",
        "participants": 19,
        "participant_zone_rows": 57,
        "rna_libraries": 57,
        "rrbs_libraries": 57,
        "pairing": "adjacent_section",
        "same_section": False,
        "same_cell": False,
        "same_molecule": False,
        "inference_unit": "participant",
        "zones_sections_libraries_runs_are_independent_replicates": False,
        "split_rule": "all_three_zones_and_both_assays_from_one_participant_remain_in_one_outer_fold",
    }
    missingness = {
        "schema_version": "masld-bench-gse105127-missingness-v1",
        "states": [
            "observed",
            "structurally_missing",
            "not_applicable",
            "below_qc",
            "unavailable_permission",
            "join_unresolved",
            "withheld_sealed",
            "derivable_not_processed",
            "unmapped_reference",
        ],
        "missing_as_zero_forbidden": True,
        "rrbs_uncovered_cpg_state": "structurally_missing_or_below_qc_never_zero",
    }
    for name, payload in (
        ("authorities/source_semantics.json", source_semantics),
        ("authorities/remap_contract.json", remap),
        ("authorities/rna_preprocessing_contract.json", rna_contract),
        ("authorities/topology.json", topology),
        ("authorities/missingness.json", missingness),
    ):
        write_json(output / name, payload)
    audit = {
        "schema_version": "masld-bench-gse105127-no-fit-activation-v1",
        "status": "design_ready_no_fit_remap_and_rna_quantification_pending",
        "participants": 19,
        "participant_zone_rows": 57,
        "rrbs_files": 57,
        "rrbs_cytosine_rows_exhaustively_audited": EXPECTED_CYTOSINE_ROWS,
        "rna_runs": 57,
        "adjacent_section_only": True,
        "model_input_labels_present": False,
        "phenotypes_or_histology_read_for_feature_construction": False,
        "normalization_run": False,
        "model_fit_run": False,
        "model_scoring_run": False,
        "rrbs_remap_run": False,
        "rna_quantification_run": False,
        "ena_runinfo_sha256": sha256(ena_payload).hexdigest(),
        "source_artifacts": {
            "geo": EXPECTED_GEO_ARTIFACTS_SHA256,
            "join": EXPECTED_JOIN_ARTIFACTS_SHA256,
            "rrbs": EXPECTED_RRBS_ARTIFACTS_SHA256,
        },
        "remaining_activation_gates": [
            "freeze_and_verify_exact_source_FASTA_and_both_chain_files",
            "run_exhaustive_strand_pair_auxiliary_consistency_audit",
            "materialize_failure_aware_roundtrip_CpG_crosswalk",
            "download_and_quantify_57_RNA_FASTQs_without_labels",
            "freeze_task_evaluator_and_fold_authority_before_any_fit",
        ],
    }
    write_json(output / "activation_audit.json", audit)
    freeze_tree(
        output,
        {
            "artifact_class": "gse105127_label_free_assay_native_activation_design",
            "status": audit["status"],
            "fit_or_score_performed": False,
        },
    )
    verify_frozen_tree(output)
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geo-root", type=Path, required=True)
    parser.add_argument("--join-root", type=Path, required=True)
    parser.add_argument("--rrbs-root", type=Path, required=True)
    parser.add_argument("--ena-runinfo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit = build_activation(
        geo_root=args.geo_root,
        join_root=args.join_root,
        rrbs_root=args.rrbs_root,
        ena_runinfo=args.ena_runinfo,
        output=args.output,
    )
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
