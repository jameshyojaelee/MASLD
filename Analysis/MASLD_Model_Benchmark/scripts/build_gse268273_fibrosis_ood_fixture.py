#!/usr/bin/env python3
"""Build a no-fit GSE268273 fibrosis/etiology OOD transfer fixture."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import hmac
import io
import json
from pathlib import Path
import re
import subprocess
import tomllib
from typing import Any, Mapping, Sequence
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
import zipfile

import numpy as np

from masld_bench.contracts import TaskSpec


ARTICLE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11474426/fullTextXML"
ENA_PARAMETERS = {
    "accession": "PRJNA1116068",
    "result": "read_run",
    "fields": (
        "study_accession,sample_accession,experiment_accession,run_accession,"
        "scientific_name,library_layout,fastq_ftp,fastq_bytes,fastq_md5"
    ),
    "format": "tsv",
    "download": "false",
}
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urlencode(
    ENA_PARAMETERS
)
EXPECTED_SOURCE_ARTIFACTS_SHA256 = (
    "883f870e2a76f78c39d257c24f8e795d2e74f1bd32aba785251b52aa7cb8229e"
)
EXPECTED_MATRIX_ARTIFACTS_SHA256 = (
    "3f934ab9c6ee730c6bd0887b10d963f8210f274c36dc0878cb814c85026dc038"
)
EXPECTED_ACTIVATION_ARTIFACTS_SHA256 = (
    "d9a9e8b1a99739908587f0b854bb3adc472a3543b769eeeb0813544e78d57c55"
)
EXPECTED_TRAINING_MOLECULAR_ARTIFACTS_SHA256 = (
    "2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456"
)
EXPECTED_TRAINING_JOIN_ARTIFACTS_SHA256 = (
    "e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580"
)
EXPECTED_GSE260666_ACTIVATION_ARTIFACTS_SHA256 = (
    "734d45aee49221ae9b72c2b79de6b4959ae1d3d01ac85420c89135c05469e8f7"
)
EXPECTED_MICROARRAY_PARTICIPANTS_ARTIFACTS_SHA256 = (
    "1ba06484d9d064198c6c85d480c6ae2c0dadea712d0b53f9f15d81c6aa35fddb"
)
EXPECTED_MICROARRAY_OVERLAP_ARTIFACTS_SHA256 = (
    "e577e1f581a484095748a888cfaebcfb83548d452bab523eab8ef19dbb0c89eb"
)
EXPECTED_TASK_SHA256 = (
    "0408cf4335326494261265c7730e1aa246d86adf5299fa40bb875b5735bbce5c"
)
EXPECTED_PROMOTION_SHA256 = (
    "1288cd1c190448a0749cd64182849fca1f54dbaedf0248925140956c7f8da072"
)
EXPECTED_PARTICIPANTS = 109
EXPECTED_ENA_RUNS = 824
EXPECTED_GROUPS = {"classic-nafld": 40, "imid-nafld": 69}
EXPECTED_FIBROSIS = {"F0": 36, "F1": 41, "F2": 13, "F3": 14, "F4": 5}
EXPECTED_SEX = {"female": 52, "male": 57}
FIVE_COHORTS = (
    "GSE126848",
    "GSE130970",
    "GSE135251",
    "GSE162694",
    "GSE213621",
)


class GSE268273FixtureError(RuntimeError):
    """Raised when GSE268273 cannot support the frozen no-fit fixture."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
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


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE268273FixtureError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def download(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "masld-bench-gse268273-audit/1"})
    with urlopen(request, timeout=120) as response:
        payload = response.read()
    if not payload:
        raise GSE268273FixtureError(f"official source returned no bytes: {url}")
    return payload


def characteristics(row: Mapping[str, str]) -> dict[str, str]:
    output: dict[str, str] = {}
    for value in json.loads(row["characteristics_json"]):
        if ": " not in value:
            raise GSE268273FixtureError("GEO characteristic lacks key-value separator")
        key, item = value.split(": ", 1)
        normalized = key.strip().lower()
        if normalized in output:
            raise GSE268273FixtureError(f"duplicated GEO characteristic: {normalized}")
        output[normalized] = item.strip()
    return output


def _row_id(sample_accession: str) -> str:
    return "g268_" + hashlib.sha256(
        f"gse268273_ood_v1\0{sample_accession}".encode("utf-8")
    ).hexdigest()[:20]


def _normalize_boolean(value: str) -> str:
    normalized = value.strip().lower()
    if normalized not in {"yes", "no"}:
        raise GSE268273FixtureError(f"source yes/no value differs: {value!r}")
    return normalized


def _numeric(value: str | None) -> tuple[str, str]:
    if value is None or not value.strip():
        return "", "structurally_missing"
    normalized = value.strip().replace(",", ".")
    try:
        number = float(normalized)
    except ValueError as error:
        raise GSE268273FixtureError(f"source numeric value differs: {value!r}") from error
    if not np.isfinite(number):
        raise GSE268273FixtureError("source numeric value is nonfinite")
    return format(number, ".17g"), "observed"


def parse_article_xml(payload: bytes) -> dict[str, Any]:
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as error:
        raise GSE268273FixtureError("Europe PMC article XML is invalid") from error
    text = " ".join(" ".join(root.itertext()).split())
    required = (
        "PMC11474426",
        "10.1016/j.jhepr.2024.101167",
        "cross-sectional, case-control study",
        "TSC of 109 histologically characterized cases",
        "Total RNA was extracted from freshly frozen biopsies",
        "human reference genome (GRCh38)",
        "gencode release 36",
        "RSEM (version 1.3.0)",
        "limma-voom",
        "adjusting for BMI, sex, and fibrosis severity",
    )
    missing = [value for value in required if value not in text]
    if missing:
        raise GSE268273FixtureError(
            "source article evidence differs: " + "; ".join(missing)
        )
    license_evidence = " ".join(
        " ".join(
            [*element.attrib.values(), *element.itertext()]
        )
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "license"
    )
    normalized_license = " ".join(license_evidence.lower().split())
    if not (
        "creativecommons.org/licenses/by-nc-nd/4.0" in normalized_license
        or ("cc by-nc-nd" in normalized_license and "4.0" in normalized_license)
    ):
        raise GSE268273FixtureError(
            "source article CC BY-NC-ND 4.0 license evidence differs"
        )
    return {
        "pmcid": "PMC11474426",
        "pmid": "39411649",
        "doi": "10.1016/j.jhepr.2024.101167",
        "design": "cross-sectional_case_control",
        "transcriptomic_participant_wording": "109_histologically_characterized_cases",
        "biopsy_selection": "suspected_significant_fibrosis_by_liver_stiffness",
        "group_matching": "source_reports_age_and_fibrosis_matching",
        "source_reference": "GRCh38_GENCODE_release_36",
        "source_quantification": "STAR_2.5.3a_RSEM_1.3.0",
        "source_de_model": "limma_voom_adjusted_BMI_sex_fibrosis",
        "article_license": "CC_BY_NC_ND_4_0",
        "article_license_does_not_assign_dataset_terms": True,
    }


def parse_ena_tsv(
    payload: bytes,
    expected_biosamples: set[str],
    expected_experiments: set[str],
    *,
    expected_runs: int = EXPECTED_ENA_RUNS,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    try:
        rows = list(
            csv.DictReader(io.StringIO(payload.decode("utf-8")), delimiter="\t")
        )
    except UnicodeDecodeError as error:
        raise GSE268273FixtureError("ENA run manifest is not UTF-8") from error
    required = {
        "run_accession",
        "study_accession",
        "sample_accession",
        "experiment_accession",
        "scientific_name",
        "library_layout",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
    }
    if not rows or not required <= set(rows[0]):
        raise GSE268273FixtureError("ENA run manifest schema differs")
    if len(rows) != expected_runs:
        raise GSE268273FixtureError("ENA run count differs")
    if (
        len({row["run_accession"] for row in rows}) != len(rows)
        or {row["sample_accession"] for row in rows} != expected_biosamples
        or {row["experiment_accession"] for row in rows} != expected_experiments
        or any(
            row["study_accession"] != "PRJNA1116068"
            or row["scientific_name"] != "Homo sapiens"
            or row["library_layout"] != "SINGLE"
            for row in rows
        )
    ):
        raise GSE268273FixtureError("ENA topology or library identity differs")
    experiments_by_sample: dict[str, set[str]] = defaultdict(set)
    runs_per_experiment = Counter()
    total_bytes = 0
    for row in rows:
        experiments_by_sample[row["sample_accession"]].add(
            row["experiment_accession"]
        )
        runs_per_experiment[row["experiment_accession"]] += 1
        try:
            sizes = [int(value) for value in row["fastq_bytes"].split(";")]
        except ValueError as error:
            raise GSE268273FixtureError("ENA FASTQ byte field differs") from error
        paths = row["fastq_ftp"].split(";")
        md5s = row["fastq_md5"].split(";")
        if (
            not row["fastq_ftp"]
            or len(paths) != len(sizes)
            or len(md5s) != len(sizes)
            or any(value <= 0 for value in sizes)
            or any(re.fullmatch(r"[0-9a-fA-F]{32}", value) is None for value in md5s)
        ):
            raise GSE268273FixtureError("ENA FASTQ path, byte size, or MD5 differs")
        total_bytes += sum(sizes)
    if any(len(values) != 1 for values in experiments_by_sample.values()):
        raise GSE268273FixtureError("one BioSample maps to multiple SRA experiments")
    run_count_distribution = Counter(runs_per_experiment.values())
    if expected_runs == EXPECTED_ENA_RUNS and run_count_distribution != Counter(
        {4: 36, 8: 49, 12: 24}
    ):
        raise GSE268273FixtureError("technical-run distribution differs")
    receipt = {
        "runs": len(rows),
        "biosamples": len(experiments_by_sample),
        "experiments": len(runs_per_experiment),
        "layout": "SINGLE",
        "runs_per_experiment": {
            str(key): value for key, value in sorted(run_count_distribution.items())
        },
        "fastq_bytes": total_bytes,
        "fastq_gib": total_bytes / (1024**3),
        "biological_unit": "participant",
        "run_role": "technical_partition_not_biological_replicate",
    }
    return rows, receipt


_XLSX_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_XLSX_DOCUMENT_REL_NS = (
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
)
_XLSX_PACKAGE_REL_NS = (
    "http://schemas.openxmlformats.org/package/2006/relationships"
)


def _xlsx_column_index(cell_reference: str) -> int:
    match = re.fullmatch(r"([A-Z]+)[1-9][0-9]*", cell_reference)
    if match is None:
        raise GSE268273FixtureError("XLSX cell reference differs")
    index = 0
    for letter in match.group(1):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def _iter_xlsx_rows(workbook_path: Path):
    """Yield values from the single deposited worksheet using only stdlib XML."""
    with zipfile.ZipFile(workbook_path) as archive:
        try:
            workbook = ET.parse(archive.open("xl/workbook.xml")).getroot()
            relationships = ET.parse(
                archive.open("xl/_rels/workbook.xml.rels")
            ).getroot()
        except (KeyError, ET.ParseError) as error:
            raise GSE268273FixtureError("normalized XLSX workbook structure differs") from error
        sheets = workbook.findall(f".//{{{_XLSX_MAIN_NS}}}sheet")
        if len(sheets) != 1:
            raise GSE268273FixtureError(
                "normalized expression workbook sheet count differs"
            )
        relationship_id = sheets[0].get(f"{{{_XLSX_DOCUMENT_REL_NS}}}id")
        targets = {
            item.get("Id"): item.get("Target")
            for item in relationships.findall(
                f".//{{{_XLSX_PACKAGE_REL_NS}}}Relationship"
            )
        }
        sheet_target = targets.get(relationship_id)
        if (
            sheets[0].get("name") != "RNAseq_expression_values_limma_"
            or sheet_target != "worksheets/sheet1.xml"
        ):
            raise GSE268273FixtureError("normalized expression worksheet identity differs")
        try:
            shared_root = ET.parse(archive.open("xl/sharedStrings.xml")).getroot()
        except (KeyError, ET.ParseError) as error:
            raise GSE268273FixtureError("normalized XLSX shared strings differ") from error
        shared_strings = [
            "".join(
                value.text or ""
                for value in item.iter(f"{{{_XLSX_MAIN_NS}}}t")
            )
            for item in shared_root.findall(f"{{{_XLSX_MAIN_NS}}}si")
        ]
        try:
            sheet_handle = archive.open("xl/worksheets/sheet1.xml")
        except KeyError as error:
            raise GSE268273FixtureError("normalized XLSX worksheet is absent") from error
        with sheet_handle:
            try:
                for _, row_element in ET.iterparse(sheet_handle, events=("end",)):
                    if row_element.tag != f"{{{_XLSX_MAIN_NS}}}row":
                        continue
                    row: list[Any] = []
                    for cell in row_element.findall(f"{{{_XLSX_MAIN_NS}}}c"):
                        reference = cell.get("r", "")
                        column = _xlsx_column_index(reference)
                        if column < len(row):
                            raise GSE268273FixtureError(
                                "normalized XLSX cell order differs"
                            )
                        row.extend([None] * (column + 1 - len(row)))
                        cell_type = cell.get("t")
                        value = cell.find(f"{{{_XLSX_MAIN_NS}}}v")
                        if value is None or value.text is None:
                            row[column] = None
                        elif cell_type == "s":
                            try:
                                row[column] = shared_strings[int(value.text)]
                            except (ValueError, IndexError) as error:
                                raise GSE268273FixtureError(
                                    "normalized XLSX shared-string index differs"
                                ) from error
                        elif cell_type in {"str", "e"}:
                            row[column] = value.text
                        elif cell_type == "b":
                            row[column] = value.text == "1"
                        elif cell_type is None or cell_type == "n":
                            try:
                                row[column] = float(value.text)
                            except ValueError as error:
                                raise GSE268273FixtureError(
                                    "normalized XLSX numeric cell differs"
                                ) from error
                        else:
                            raise GSE268273FixtureError(
                                "normalized XLSX cell type differs"
                            )
                    yield row
                    row_element.clear()
            except ET.ParseError as error:
                raise GSE268273FixtureError(
                    "normalized XLSX worksheet XML differs"
                ) from error


def load_deposited_expression(
    workbook_path: Path,
    crosswalk_path: Path,
    participant_order: Sequence[str],
) -> tuple[np.ndarray, list[dict[str, str]], dict[str, Any]]:
    crosswalk_fields, crosswalk = read_tsv(crosswalk_path)
    required_crosswalk = {
        "source_feature_id",
        "mapping_state",
        "gencode_v49_stable_id",
        "gencode_v49_gene_name",
        "gencode_v49_gene_type",
    }
    if not required_crosswalk <= set(crosswalk_fields) or len(crosswalk) != 14149:
        raise GSE268273FixtureError("GSE268273 gene crosswalk differs")
    iterator = _iter_xlsx_rows(workbook_path)
    header = list(next(iterator))
    sample_axis = [str(value) for value in header[3:] if value is not None]
    if len(sample_axis) != len(participant_order) or set(sample_axis) != set(
        participant_order
    ):
        raise GSE268273FixtureError("normalized expression participant axis differs")
    lookup = {value: index for index, value in enumerate(sample_axis)}
    reorder = [lookup[value] for value in participant_order]
    feature_rows: list[dict[str, str]] = []
    values: list[list[float]] = []
    for index, row in enumerate(iterator):
        if index >= len(crosswalk) or row[0] is None:
            raise GSE268273FixtureError("normalized expression feature axis differs")
        cross = crosswalk[index]
        if str(row[0]) != cross["source_feature_id"] or len(row) < 3 + len(sample_axis):
            raise GSE268273FixtureError("workbook and crosswalk feature order differs")
        try:
            numeric = [float(row[3 + source_index]) for source_index in reorder]
        except (TypeError, ValueError) as error:
            raise GSE268273FixtureError("normalized expression contains nonnumeric data") from error
        values.append(numeric)
        feature_rows.append(cross)
    matrix = np.asarray(values, dtype=np.float64).T
    if matrix.shape != (len(participant_order), 14149) or not np.all(
        np.isfinite(matrix)
    ):
        raise GSE268273FixtureError("normalized expression matrix differs")
    mapped = [
        row["gencode_v49_stable_id"]
        for row in feature_rows
        if row["mapping_state"] == "stable_id"
    ]
    target_counts = Counter(mapped)
    unique_rows = [
        index
        for index, row in enumerate(feature_rows)
        if row["mapping_state"] == "stable_id"
        and target_counts[row["gencode_v49_stable_id"]] == 1
    ]
    duplicate_targets = sorted(
        stable for stable, count in target_counts.items() if count > 1
    )
    if (
        len(mapped) != 14089
        or len(target_counts) != 14078
        or len(duplicate_targets) != 11
        or len(unique_rows) != 14067
    ):
        raise GSE268273FixtureError("GENCODE v49 mapping multiplicity differs")
    receipt = {
        "participants": matrix.shape[0],
        "source_features": matrix.shape[1],
        "v49_mapped_rows": len(mapped),
        "v49_unique_stable_ids": len(target_counts),
        "v49_one_to_one_rows": len(unique_rows),
        "v49_duplicate_target_stable_ids": len(duplicate_targets),
        "v49_unmapped_rows": sum(
            row["mapping_state"] == "unmapped_stable_id" for row in feature_rows
        ),
        "minimum": float(matrix.min()),
        "maximum": float(matrix.max()),
        "measurement": "source_global_limma_voom_transformed_expression",
        "scored_model_input_allowed": False,
    }
    return matrix, feature_rows, {**receipt, "one_to_one_row_indices": unique_rows}


def _rank_rows(matrix: np.ndarray) -> np.ndarray:
    output = np.empty(matrix.shape, dtype=np.float64)
    for row_index, row in enumerate(matrix):
        order = np.argsort(row, kind="mergesort")
        sorted_values = row[order]
        ranks = np.empty(len(row), dtype=np.float64)
        start = 0
        while start < len(row):
            stop = start + 1
            while stop < len(row) and sorted_values[stop] == sorted_values[start]:
                stop += 1
            ranks[order[start:stop]] = (start + 1 + stop) / 2.0
            start = stop
        output[row_index] = ranks
    return output


def _row_unit(matrix: np.ndarray) -> np.ndarray:
    centered = np.asarray(matrix, dtype=np.float64) - np.mean(
        matrix, axis=1, keepdims=True
    )
    scale = np.sqrt(np.sum(centered * centered, axis=1, keepdims=True))
    if np.any(scale <= np.finfo(np.float64).eps):
        raise GSE268273FixtureError("expression fingerprint row is constant")
    return centered / scale


def _log_cpm(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or np.any(values < 0) or not np.all(np.isfinite(values)):
        raise GSE268273FixtureError("count fingerprint input is invalid")
    totals = values.sum(axis=1)
    if np.any(totals <= 0):
        raise GSE268273FixtureError("count fingerprint input has an empty library")
    return np.log1p(values * (1_000_000.0 / totals[:, None]))


def _signature(key: bytes, row_id: str, values: np.ndarray) -> str:
    quantized = np.round(values, 6).astype("<f8", copy=False)
    return hmac.new(
        key,
        row_id.encode("utf-8") + b"\0" + quantized.tobytes(),
        hashlib.sha256,
    ).hexdigest()


def load_featurecounts(
    path: Path, target_genes: set[str]
) -> tuple[np.ndarray, list[str], list[str]]:
    with path.open(encoding="utf-8", errors="strict") as handle:
        program = handle.readline()
        header = handle.readline().rstrip("\n").split("\t")
        if not program.startswith("# Program:featureCounts") or header[:6] != [
            "Geneid",
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
        ]:
            raise GSE268273FixtureError(f"featureCounts header differs: {path}")
        sample_ids = [
            Path(value).name.split(".Aligned", 1)[0] for value in header[6:]
        ]
        if len(sample_ids) != len(set(sample_ids)):
            raise GSE268273FixtureError(f"featureCounts sample axis is duplicated: {path}")
        by_gene: dict[str, np.ndarray] = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                raise GSE268273FixtureError(f"featureCounts row width differs: {path}")
            stable = fields[0].split(".", 1)[0]
            if stable not in target_genes:
                continue
            if stable in by_gene:
                raise GSE268273FixtureError(f"featureCounts stable gene duplicated: {path}")
            try:
                counts = np.asarray([int(value) for value in fields[6:]], dtype=np.float64)
            except ValueError as error:
                raise GSE268273FixtureError("featureCounts value is not integer") from error
            if np.any(counts < 0):
                raise GSE268273FixtureError("featureCounts value is negative")
            by_gene[stable] = counts
    genes = sorted(by_gene)
    if len(genes) < 10000:
        raise GSE268273FixtureError(f"too few common featureCounts genes: {path}")
    matrix = np.column_stack([by_gene[gene] for gene in genes])
    return matrix, genes, sample_ids


def load_gse267145(
    molecular: Path, target_genes: set[str]
) -> tuple[np.ndarray, list[str], list[str]]:
    fields, feature_rows = read_tsv(molecular / "rna_feature_axis.tsv")
    participant_fields, participant_rows = read_tsv(molecular / "participant_axis.tsv")
    if "stable_gene_id" not in fields or "participant_id" not in participant_fields:
        raise GSE268273FixtureError("GSE267145 molecular axes differ")
    genes_all = [row["stable_gene_id"] for row in feature_rows]
    if len(genes_all) != len(set(genes_all)):
        raise GSE268273FixtureError("GSE267145 stable-gene axis is duplicated")
    lookup = {gene: index for index, gene in enumerate(genes_all)}
    genes = sorted(target_genes & set(genes_all))
    if len(genes) < 10000:
        raise GSE268273FixtureError("too few common GSE267145 genes")
    values = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    if values.shape != (len(participant_rows), len(genes_all)):
        raise GSE268273FixtureError("GSE267145 molecular matrix differs")
    matrix = np.asarray(values[:, [lookup[gene] for gene in genes]], dtype=np.float64)
    if np.any(matrix < 0) or not np.all(np.isfinite(matrix)):
        raise GSE268273FixtureError("GSE267145 fingerprint matrix differs")
    return matrix, genes, [row["participant_id"] for row in participant_rows]


def compare_expression(
    *,
    target_values: np.ndarray,
    target_genes: Sequence[str],
    target_row_ids: Sequence[str],
    source_values: np.ndarray,
    source_genes: Sequence[str],
    source_row_ids: Sequence[str],
    source_measurement: str,
    key: bytes,
) -> dict[str, Any]:
    target_lookup = {gene: index for index, gene in enumerate(target_genes)}
    source_lookup = {gene: index for index, gene in enumerate(source_genes)}
    common = sorted(set(target_lookup) & set(source_lookup))
    if len(common) < 10000:
        raise GSE268273FixtureError("expression fingerprint common axis is too small")
    target = target_values[:, [target_lookup[gene] for gene in common]]
    source_raw = source_values[:, [source_lookup[gene] for gene in common]]
    source = (
        _log_cpm(source_raw)
        if source_measurement in {"raw_counts", "continuous_count_estimates"}
        else source_raw
    )
    pearson = _row_unit(target) @ _row_unit(source).T
    spearman = _row_unit(_rank_rows(target)) @ _row_unit(_rank_rows(source)).T
    joint = np.minimum(pearson, spearman)
    flat = int(np.argmax(joint))
    target_index, source_index = np.unravel_index(flat, joint.shape)
    maximum_pearson = float(pearson[target_index, source_index])
    maximum_spearman = float(spearman[target_index, source_index])
    near_duplicate = maximum_pearson >= 0.999 and maximum_spearman >= 0.999
    target_signatures = sorted(
        _signature(key, target_row_ids[index], target[index])
        for index in range(target.shape[0])
    )
    source_signatures = sorted(
        _signature(key, source_row_ids[index], source[index])
        for index in range(source.shape[0])
    )
    return {
        "common_gencode_v49_stable_genes": len(common),
        "target_participants": target.shape[0],
        "source_profiles": source.shape[0],
        "source_measurement": source_measurement,
        "maximum_cross_source_pearson": maximum_pearson,
        "maximum_cross_source_spearman": maximum_spearman,
        "maximum_target_row_hmac_sha256": hmac.new(
            key, target_row_ids[target_index].encode("utf-8"), hashlib.sha256
        ).hexdigest(),
        "maximum_source_row_hmac_sha256": hmac.new(
            key, source_row_ids[source_index].encode("utf-8"), hashlib.sha256
        ).hexdigest(),
        "near_duplicate_rule": "Pearson >= 0.999 and tied-rank Spearman >= 0.999",
        "near_duplicate_detected": near_duplicate,
        "target_signature_set_sha256": hashlib.sha256(
            "\n".join(target_signatures).encode("utf-8")
        ).hexdigest(),
        "source_signature_set_sha256": hashlib.sha256(
            "\n".join(source_signatures).encode("utf-8")
        ).hexdigest(),
        "raw_expression_fingerprints_retained": False,
    }


def within_target_duplicate_audit(
    values: np.ndarray, row_ids: Sequence[str], key: bytes
) -> dict[str, Any]:
    pearson = _row_unit(values) @ _row_unit(values).T
    spearman = _row_unit(_rank_rows(values)) @ _row_unit(_rank_rows(values)).T
    joint = np.minimum(pearson, spearman)
    joint[np.triu_indices_from(joint, 0)] = -np.inf
    flat = int(np.argmax(joint))
    left, right = np.unravel_index(flat, joint.shape)
    maximum_pearson = float(pearson[left, right])
    maximum_spearman = float(spearman[left, right])
    return {
        "profiles": values.shape[0],
        "maximum_distinct_profile_pearson": maximum_pearson,
        "maximum_distinct_profile_spearman": maximum_spearman,
        "maximum_left_row_hmac_sha256": hmac.new(
            key, row_ids[left].encode("utf-8"), hashlib.sha256
        ).hexdigest(),
        "maximum_right_row_hmac_sha256": hmac.new(
            key, row_ids[right].encode("utf-8"), hashlib.sha256
        ).hexdigest(),
        "near_duplicate_detected": maximum_pearson >= 0.999
        and maximum_spearman >= 0.999,
    }


def explicit_comparator_aliases(
    training_join: Path,
    gse260666_activation: Path,
    microarray_participants: Path,
    microarray_overlap: Path,
) -> dict[str, set[str]]:
    output: dict[str, set[str]] = {
        "gse267145_znf469_human_liver": {
            "GSE267145",
            "GSE267119",
            "GSE269412",
        },
        "gse260666_bulk_rna": {"GSE260666", "PRJNA1082656"},
        "gse31803_gse49541_fibrosis_array": {"GSE31803", "GSE49541"},
        "antwerp_inserm_shared": {"GSE106737", "GSE83452"},
    }
    fields, rows = read_tsv(training_join / "participant_join.tsv")
    for field in ("participant_id", "rna_gsm", "h3k27ac_gsm"):
        if field not in fields:
            raise GSE268273FixtureError("GSE267145 alias axis differs")
        output["gse267145_znf469_human_liver"].update(
            row[field] for row in rows if row[field]
        )
    alias_pattern = re.compile(r"\b(?:GSE|GSM|SAMN|SRX|PRJNA)[0-9]+\b")
    with gzip.open(
        training_join / "raw" / "GSE267145_family.soft.gz",
        "rt",
        encoding="utf-8",
        errors="strict",
    ) as handle:
        for line in handle:
            output["gse267145_znf469_human_liver"].update(alias_pattern.findall(line))
    fields, rows = read_tsv(gse260666_activation / "participant_join.tsv")
    for field in (
        "participant_id",
        "sample_accession",
        "biosample_accession",
        "sra_experiment",
    ):
        if field not in fields:
            raise GSE268273FixtureError("GSE260666 alias axis differs")
        output["gse260666_bulk_rna"].update(row[field] for row in rows if row[field])
    fields, rows = read_tsv(microarray_participants / "gse49541_participants.tsv")
    for field in ("participant_id", "sample_accession", "source_sample_key"):
        if field not in fields:
            raise GSE268273FixtureError("GSE49541 alias axis differs")
        output["gse31803_gse49541_fibrosis_array"].update(
            row[field] for row in rows if row[field]
        )
    fields, rows = read_tsv(microarray_participants / "gse83452_records.tsv")
    for field in (
        "participant_id",
        "sample_accession",
        "source_sample_number",
        "paired_baseline_accession",
    ):
        if field not in fields:
            raise GSE268273FixtureError("GSE83452 alias axis differs")
        output["antwerp_inserm_shared"].update(
            row[field]
            for row in rows
            if row[field] and row[field] not in {"not_applicable", "structurally_missing"}
        )
    for filename in (
        "cross_accession_participant_aliases.tsv",
        "cross_accession_sample_aliases.tsv",
    ):
        fields, rows = read_tsv(microarray_overlap / filename)
        for row in rows:
            for field in fields:
                value = row[field]
                if value and value not in {"not_applicable", "structurally_missing"}:
                    output["antwerp_inserm_shared"].add(value)
    return output


def prior_project_use(
    project_root: Path,
    benchmark_root: Path,
    aliases: Sequence[str],
    key: bytes,
) -> dict[str, Any]:
    excluded = benchmark_root.relative_to(project_root).as_posix()
    command = [
        "rg",
        "-l",
        "-F",
        "--hidden",
        "--glob",
        f"!{excluded}/**",
        "--glob",
        "!.git/**",
    ]
    for alias in sorted(set(aliases)):
        command.extend(("-e", alias))
    command.append(".")
    completed = subprocess.run(
        command,
        cwd=project_root,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode not in {0, 1}:
        raise GSE268273FixtureError("prior project-use alias scan failed")
    paths = sorted(set(completed.stdout.splitlines()))
    return {
        "search_scope": "project_root_excluding_benchmark_package_and_git_metadata",
        "aliases_searched": len(set(aliases)),
        "matched_file_count": len(paths),
        "matched_file_hmac_sha256": [
            hmac.new(key, path.encode("utf-8"), hashlib.sha256).hexdigest()
            for path in paths
        ],
        "raw_paths_retained": False,
    }


def build_fixture(
    *,
    source: Path,
    matrix: Path,
    activation: Path,
    training_molecular: Path,
    training_join: Path,
    gse260666_activation: Path,
    microarray_participants: Path,
    microarray_overlap: Path,
    task_spec: Path,
    promotion_gate: Path,
    expansion_config: Path,
    project_root: Path,
    benchmark_root: Path,
    fingerprint_salt_file: Path,
    output: Path,
    article_payload: bytes,
    ena_payload: bytes,
    expected_participants: int = EXPECTED_PARTICIPANTS,
    expected_ena_runs: int = EXPECTED_ENA_RUNS,
) -> dict[str, Any]:
    if output.exists():
        raise GSE268273FixtureError(f"refusing to overwrite fixture: {output}")
    expected_hashes = (
        (source / "ARTIFACTS.json", EXPECTED_SOURCE_ARTIFACTS_SHA256),
        (matrix / "ARTIFACTS.json", EXPECTED_MATRIX_ARTIFACTS_SHA256),
        (activation / "ARTIFACTS.json", EXPECTED_ACTIVATION_ARTIFACTS_SHA256),
        (
            training_molecular / "ARTIFACTS.json",
            EXPECTED_TRAINING_MOLECULAR_ARTIFACTS_SHA256,
        ),
        (
            training_join / "ARTIFACTS.json",
            EXPECTED_TRAINING_JOIN_ARTIFACTS_SHA256,
        ),
        (
            gse260666_activation / "ARTIFACTS.json",
            EXPECTED_GSE260666_ACTIVATION_ARTIFACTS_SHA256,
        ),
        (
            microarray_participants / "ARTIFACTS.json",
            EXPECTED_MICROARRAY_PARTICIPANTS_ARTIFACTS_SHA256,
        ),
        (
            microarray_overlap / "ARTIFACTS.json",
            EXPECTED_MICROARRAY_OVERLAP_ARTIFACTS_SHA256,
        ),
        (task_spec, EXPECTED_TASK_SHA256),
        (promotion_gate, EXPECTED_PROMOTION_SHA256),
    )
    for path, expected in expected_hashes:
        if sha256_file(path) != expected:
            raise GSE268273FixtureError(f"frozen input SHA-256 differs: {path}")
    if fingerprint_salt_file.stat().st_mode & 0o077:
        raise GSE268273FixtureError("fingerprint salt is group/world accessible")
    salt = fingerprint_salt_file.read_bytes()
    if len(salt) < 32:
        raise GSE268273FixtureError("fingerprint salt is shorter than 32 bytes")
    task = TaskSpec.from_toml(task_spec.read_text(encoding="utf-8"))
    if (
        task.task_id != "gse268273_fibrosis_etiology_ood_transport"
        or task.status.value != "blocked"
        or task.primary_metric != "participant_spearman_fibrosis_f0_f4"
        or task.promotion_gate_config_sha256 != EXPECTED_PROMOTION_SHA256
    ):
        raise GSE268273FixtureError("TaskSpec identity or blocked gate differs")
    source_fields, source_rows = read_tsv(source / "samples" / "GSE268273.tsv")
    join_fields, join_rows = read_tsv(activation / "participant_join.tsv")
    required_source = {
        "accession",
        "description_json",
        "characteristics_json",
        "relations_json",
    }
    required_join = {
        "participant_id",
        "sample_accession",
        "biosample_accession",
        "sra_experiment",
        "group",
        "fibrosis",
        "sex",
        "biological_unit",
    }
    if (
        not required_source <= set(source_fields)
        or not required_join <= set(join_fields)
        or len(source_rows) != expected_participants
        or len(join_rows) != expected_participants
    ):
        raise GSE268273FixtureError("GEO source or participant join schema differs")
    for field in (
        "participant_id",
        "sample_accession",
        "biosample_accession",
        "sra_experiment",
    ):
        values = [row[field] for row in join_rows]
        if len(values) != len(set(values)):
            raise GSE268273FixtureError(f"participant identifier is duplicated: {field}")
    if any(row["biological_unit"] != "participant" for row in join_rows):
        raise GSE268273FixtureError("GSE268273 biological unit differs")
    if expected_participants == EXPECTED_PARTICIPANTS and (
        Counter(row["group"] for row in join_rows) != Counter(EXPECTED_GROUPS)
        or Counter(row["fibrosis"] for row in join_rows) != Counter(EXPECTED_FIBROSIS)
        or Counter(row["sex"] for row in join_rows) != Counter(EXPECTED_SEX)
    ):
        raise GSE268273FixtureError("GSE268273 phenotype census differs")
    source_by_accession = {row["accession"]: row for row in source_rows}
    if set(source_by_accession) != {row["sample_accession"] for row in join_rows}:
        raise GSE268273FixtureError("GEO and participant-join accessions differ")
    article = parse_article_xml(article_payload)
    ena_rows, ena_receipt = parse_ena_tsv(
        ena_payload,
        {row["biosample_accession"] for row in join_rows},
        {row["sra_experiment"] for row in join_rows},
        expected_runs=expected_ena_runs,
    )
    participant_order = [row["participant_id"] for row in join_rows]
    deposited, feature_rows, matrix_receipt = load_deposited_expression(
        matrix
        / "raw"
        / "GSE268nnn"
        / "GSE268273_EGN_RNAseq_Normalized_data.xlsx",
        activation / "gene_crosswalk.tsv",
        participant_order,
    )
    if expected_participants != EXPECTED_PARTICIPANTS:
        raise GSE268273FixtureError("production fixture requires the full source cohort")
    row_ids = [_row_id(row["sample_accession"]) for row in join_rows]
    row_id_by_biosample = {
        row["biosample_accession"]: row_id
        for row, row_id in zip(join_rows, row_ids, strict=True)
    }
    run_count_by_biosample = Counter(row["sample_accession"] for row in ena_rows)
    byte_count_by_biosample: Counter[str] = Counter()
    for row in ena_rows:
        byte_count_by_biosample[row["sample_accession"]] += sum(
            int(value) for value in row["fastq_bytes"].split(";")
        )
    output.mkdir(parents=True)
    source_root = output / "source_evidence"
    preprocessing_root = output / "preprocessing_input"
    model_root = output / "model_input"
    evaluator_root = output / "evaluator_only"
    audit_root = output / "contamination_audit"
    for path in (source_root, preprocessing_root, model_root, evaluator_root, audit_root):
        path.mkdir()
    (source_root / "article.xml").write_bytes(article_payload)
    (source_root / "ena_run_manifest.tsv").write_bytes(ena_payload)
    source_receipt = {
        "schema_version": "masld-bench-gse268273-source-evidence-v1",
        "status": "passed",
        "article": article,
        "article_url": ARTICLE_URL,
        "article_sha256": hashlib.sha256(article_payload).hexdigest(),
        "ena_url": ENA_URL,
        "ena_sha256": hashlib.sha256(ena_payload).hexdigest(),
        "ena": ena_receipt,
        "geo_records": expected_participants,
        "unique_participants": expected_participants,
        "repeated_biopsy_or_specimen_detected": False,
        "multiple_runs_are_technical_partitions": True,
        "source_record_to_participant_is_one_to_one": True,
    }
    (source_root / "receipt.json").write_text(
        json.dumps(source_receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_tsv(
        preprocessing_root / "raw_run_manifest.tsv",
        (
            "row_id",
            "run_accession",
            "experiment_accession",
            "biosample_accession",
            "library_layout",
            "fastq_ftp",
            "fastq_bytes",
            "fastq_md5",
            "replicate_role",
        ),
        [
            {
                "row_id": row_id_by_biosample[row["sample_accession"]],
                "run_accession": row["run_accession"],
                "experiment_accession": row["experiment_accession"],
                "biosample_accession": row["sample_accession"],
                "library_layout": row["library_layout"],
                "fastq_ftp": row["fastq_ftp"],
                "fastq_bytes": row["fastq_bytes"],
                "fastq_md5": row["fastq_md5"].lower(),
                "replicate_role": "technical_run_partition",
            }
            for row in sorted(
                ena_rows,
                key=lambda value: (
                    row_id_by_biosample[value["sample_accession"]],
                    value["run_accession"],
                ),
            )
        ],
    )
    preprocessing_receipt = {
        "schema_version": "masld-bench-gse268273-preprocessing-input-v1",
        "status": "raw_manifest_frozen_materialization_pending",
        "participants": expected_participants,
        "technical_runs": len(ena_rows),
        "raw_fastq_bytes": ena_receipt["fastq_bytes"],
        "raw_fastq_materialized": False,
        "technical_runs_collapsed_to_participant": False,
        "quantification_complete": False,
        "training_frozen_preprocessing_available": False,
        "labels_included": False,
        "source_global_voom_allowed": False,
    }
    (preprocessing_root / "receipt.json").write_text(
        json.dumps(preprocessing_receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_tsv(
        model_root / "participant_axis.tsv",
        (
            "row_index",
            "row_id",
            "cohort_family_id",
            "rna_observation_state",
            "raw_technical_run_count",
            "raw_fastq_bytes",
        ),
        [
            {
                "row_index": index,
                "row_id": row_id,
                "cohort_family_id": "gse268273_imid_masld",
                "rna_observation_state": "derivable_not_processed",
                "raw_technical_run_count": run_count_by_biosample[
                    row["biosample_accession"]
                ],
                "raw_fastq_bytes": byte_count_by_biosample[
                    row["biosample_accession"]
                ],
            }
            for index, (row, row_id) in enumerate(
                zip(join_rows, row_ids, strict=True)
            )
        ],
    )
    mapped_counts = Counter(
        row["gencode_v49_stable_id"]
        for row in feature_rows
        if row["mapping_state"] == "stable_id"
    )
    rows_by_stable: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in feature_rows:
        if row["mapping_state"] == "stable_id":
            rows_by_stable[row["gencode_v49_stable_id"]].append(row)
    admitted_features = []
    for index, stable in enumerate(sorted(rows_by_stable)):
        members = sorted(rows_by_stable[stable], key=lambda row: row["source_feature_id"])
        names = {row["gencode_v49_gene_name"] for row in members}
        types = {row["gencode_v49_gene_type"] for row in members}
        if len(names) != 1 or len(types) != 1:
            raise GSE268273FixtureError("duplicate stable-gene annotations disagree")
        admitted_features.append(
            {
                "feature_index": index,
                "stable_gene_id": stable,
                "source_feature_count": len(members),
                "source_feature_ids": ";".join(
                    row["source_feature_id"] for row in members
                ),
                "gencode_v49_gene_name": next(iter(names)),
                "gencode_v49_gene_type": next(iter(types)),
                "raw_count_aggregation": (
                    "identity"
                    if len(members) == 1
                    else "deterministic_sum_before_normalization"
                ),
                "admission_state": "raw_count_requantification_required",
            }
        )
    write_tsv(
        model_root / "candidate_gene_axis.tsv",
        (
            "feature_index",
            "stable_gene_id",
            "source_feature_count",
            "source_feature_ids",
            "gencode_v49_gene_name",
            "gencode_v49_gene_type",
            "raw_count_aggregation",
            "admission_state",
        ),
        admitted_features,
    )
    np.save(
        model_root / "rna_observed_mask.npy",
        np.zeros(expected_participants, dtype=np.bool_),
        allow_pickle=False,
    )
    model_receipt = {
        "schema_version": "masld-bench-gse268273-model-input-v1",
        "status": "outcome_free_raw_reprocessing_blocked",
        "participants": expected_participants,
        "candidate_unique_v49_genes": len(admitted_features),
        "candidate_one_to_one_v49_genes": sum(
            int(row["source_feature_count"]) == 1 for row in admitted_features
        ),
        "candidate_aggregated_v49_genes": sum(
            int(row["source_feature_count"]) > 1 for row in admitted_features
        ),
        "duplicate_mapping_resolution": "deterministic_sum_raw_counts_before_normalization",
        "rna_values_present": False,
        "rna_observation_state": "derivable_not_processed",
        "missing_encoded_as_zero": False,
        "labels_included": False,
        "clinical_covariates_included": False,
        "deposited_expression_included": False,
        "deposited_expression_scored_use_allowed": False,
        "processed_differential_expression_opened": False,
        "external_fit_or_calibration_performed": False,
        "scored_prediction_allowed": False,
        "blocker": "raw_FASTQ_requantification_and_training_frozen_preprocessing_pending",
    }
    (model_root / "receipt.json").write_text(
        json.dumps(model_receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    evaluator_rows: list[dict[str, Any]] = []
    for row_id, joined in zip(row_ids, join_rows, strict=True):
        source_row = source_by_accession[joined["sample_accession"]]
        description = json.loads(source_row["description_json"])
        if description != [joined["participant_id"]]:
            raise GSE268273FixtureError("GEO participant description differs")
        chars = characteristics(source_row)
        group = chars["case/control"].lower()
        fibrosis_text = chars["fibrosis degree"].upper()
        sex = chars["sex"].lower()
        if (
            group != joined["group"]
            or fibrosis_text != joined["fibrosis"]
            or sex != joined["sex"]
        ):
            raise GSE268273FixtureError("source phenotype and frozen join differ")
        fibrosis = int(fibrosis_text[1:])
        numeric: dict[str, str] = {}
        states: dict[str, str] = {}
        for source_key, output_key in (
            ("fasting glycemia", "fasting_glycemia"),
            ("hba1c", "hba1c"),
            ("total cholesterol", "total_cholesterol"),
            ("hdl-cholesterol", "hdl_cholesterol"),
            ("ldl-cholesterol", "ldl_cholesterol"),
            ("tryglicerides", "triglycerides_source_spelling"),
        ):
            numeric[output_key], states[f"{output_key}_state"] = _numeric(
                chars.get(source_key)
            )
        evaluator_rows.append(
            {
                "row_id": row_id,
                "participant_id": joined["participant_id"],
                "sample_accession": joined["sample_accession"],
                "source_group": group,
                "source_disease": chars["disease"],
                "fibrosis_stage": fibrosis,
                "advanced_fibrosis_f3_f4": int(fibrosis >= 3),
                "recorded_sex": sex,
                "obesity": _normalize_boolean(chars["obesity"]),
                "t2d_ir": _normalize_boolean(chars["t2d/ir"]),
                "hypertension": _normalize_boolean(chars["hypertension"]),
                **numeric,
                **states,
                "age": "",
                "age_state": "structurally_missing",
                "bmi": "",
                "bmi_state": "structurally_missing",
                "nas": "",
                "nas_state": "structurally_missing",
                "steatosis_grade": "",
                "steatosis_grade_state": "structurally_missing",
                "ballooning_grade": "",
                "ballooning_grade_state": "structurally_missing",
                "lobular_inflammation_grade": "",
                "lobular_inflammation_grade_state": "structurally_missing",
            }
        )
    evaluator_fields = tuple(evaluator_rows[0])
    write_tsv(evaluator_root / "outcomes.tsv", evaluator_fields, evaluator_rows)
    evaluator_receipt = {
        "schema_version": "masld-bench-gse268273-evaluator-outcomes-v1",
        "status": "passed_evaluator_only_no_scoring",
        "participants": expected_participants,
        "fibrosis_counts": dict(
            sorted(Counter(str(row["fibrosis_stage"]) for row in evaluator_rows).items())
        ),
        "advanced_f3_f4": sum(
            int(row["advanced_fibrosis_f3_f4"]) for row in evaluator_rows
        ),
        "source_group_counts": dict(
            sorted(Counter(row["source_group"] for row in evaluator_rows).items())
        ),
        "source_disease_counts": dict(
            sorted(Counter(row["source_disease"] for row in evaluator_rows).items())
        ),
        "recorded_sex_counts": dict(
            sorted(Counter(row["recorded_sex"] for row in evaluator_rows).items())
        ),
        "clinical_missingness": {
            field: dict(sorted(Counter(row[field] for row in evaluator_rows).items()))
            for field in states
        },
        "age_state": "structurally_missing_participant_level_public_deposit",
        "bmi_state": "structurally_missing_participant_level_public_deposit",
        "nas_and_component_histology_state": "structurally_missing_participant_level_public_deposit",
        "source_group_is_healthy_control": False,
        "source_group_is_masld_negative_control": False,
        "mapped_to_gse267145_stage3": False,
        "model_environment_access": False,
        "outcomes_scored": False,
    }
    (evaluator_root / "receipt.json").write_text(
        json.dumps(evaluator_receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    aliases = {
        "GSE268273",
        "PRJNA1116068",
        "PMC11474426",
        "39411649",
        "10.1016/j.jhepr.2024.101167",
        "NCT03760172",
        "NCT02749864",
    }
    aliases.update(
        row[field]
        for row in join_rows
        for field in (
            "participant_id",
            "sample_accession",
            "biosample_accession",
            "sra_experiment",
        )
    )
    aliases.update(row["run_accession"] for row in ena_rows)
    with expansion_config.open("rb") as handle:
        expansion = tomllib.load(handle)
    registry_overlap: dict[str, list[str]] = {}
    for family in expansion.get("cohort_family", []):
        family_id = str(family["family_id"])
        if family_id == "gse268273_imid_masld":
            continue
        overlap = sorted(aliases & set(family.get("accessions", [])))
        if overlap:
            registry_overlap[family_id] = overlap
    comparator_alias_sets = explicit_comparator_aliases(
        training_join,
        gse260666_activation,
        microarray_participants,
        microarray_overlap,
    )
    explicit_alias_overlap = {
        family_id: sorted(aliases & comparator_aliases)
        for family_id, comparator_aliases in comparator_alias_sets.items()
        if aliases & comparator_aliases
    }
    comparator_alias_counts = {
        family_id: len(comparator_aliases)
        for family_id, comparator_aliases in sorted(comparator_alias_sets.items())
    }
    fingerprint_key = hmac.new(
        salt,
        (
            "gse268273-fibrosis-ood-fingerprint-v1\0"
            + EXPECTED_SOURCE_ARTIFACTS_SHA256
            + EXPECTED_MATRIX_ARTIFACTS_SHA256
            + EXPECTED_ACTIVATION_ARTIFACTS_SHA256
            + EXPECTED_TRAINING_MOLECULAR_ARTIFACTS_SHA256
            + EXPECTED_TRAINING_JOIN_ARTIFACTS_SHA256
            + EXPECTED_GSE260666_ACTIVATION_ARTIFACTS_SHA256
            + EXPECTED_MICROARRAY_PARTICIPANTS_ARTIFACTS_SHA256
            + EXPECTED_MICROARRAY_OVERLAP_ARTIFACTS_SHA256
        ).encode("utf-8"),
        hashlib.sha256,
    ).digest()
    one_to_one_indices = matrix_receipt.pop("one_to_one_row_indices")
    target_values = deposited[:, one_to_one_indices]
    target_genes = [
        feature_rows[index]["gencode_v49_stable_id"] for index in one_to_one_indices
    ]
    comparisons: dict[str, Any] = {}
    matrix_hashes: dict[str, str] = {}
    counts_root = (
        project_root / "RNA-seq" / "Human" / "Patient_Cohorts" / "results"
    )
    for cohort in FIVE_COHORTS:
        counts_path = counts_root / cohort / "counts" / "featurecounts" / "gene_counts.txt"
        source_values, source_genes, source_ids = load_featurecounts(
            counts_path, set(target_genes)
        )
        comparisons[cohort] = compare_expression(
            target_values=target_values,
            target_genes=target_genes,
            target_row_ids=row_ids,
            source_values=source_values,
            source_genes=source_genes,
            source_row_ids=source_ids,
            source_measurement="raw_counts",
            key=hmac.new(fingerprint_key, cohort.encode("utf-8"), hashlib.sha256).digest(),
        )
        matrix_hashes[cohort] = sha256_file(counts_path)
    training_values, training_genes, training_ids = load_gse267145(
        training_molecular, set(target_genes)
    )
    comparisons["gse267145_znf469_human_liver"] = compare_expression(
        target_values=target_values,
        target_genes=target_genes,
        target_row_ids=row_ids,
        source_values=training_values,
        source_genes=training_genes,
        source_row_ids=training_ids,
        source_measurement="continuous_count_estimates",
        key=hmac.new(
            fingerprint_key, b"gse267145_znf469_human_liver", hashlib.sha256
        ).digest(),
    )
    near_duplicate_sources = sorted(
        source_id
        for source_id, result in comparisons.items()
        if result["near_duplicate_detected"]
    )
    within_target = within_target_duplicate_audit(
        target_values,
        row_ids,
        hmac.new(fingerprint_key, b"within_gse268273", hashlib.sha256).digest(),
    )
    prior_use = prior_project_use(
        project_root.resolve(strict=True),
        benchmark_root.resolve(strict=True),
        sorted(aliases),
        fingerprint_key,
    )
    contamination_status = (
        "blocked_duplicate"
        if registry_overlap
        or explicit_alias_overlap
        or near_duplicate_sources
        or within_target["near_duplicate_detected"]
        else "passed"
    )
    contamination = {
        "schema_version": "masld-bench-gse268273-contamination-audit-v1",
        "status": contamination_status,
        "cohort_family_id": "gse268273_imid_masld",
        "alias_count": len(aliases),
        "registered_nonself_family_alias_overlap": registry_overlap,
        "explicit_comparator_alias_counts": comparator_alias_counts,
        "explicit_comparator_alias_overlap": explicit_alias_overlap,
        "explicit_alias_comparators": [
            "gse267145_znf469_human_liver",
            "gse260666_bulk_rna",
            "gse31803_gse49541_fibrosis_array",
            "antwerp_inserm_shared",
        ],
        "cohort_family_overlap_detected": bool(
            registry_overlap or explicit_alias_overlap
        ),
        "prior_project_use": prior_use,
        "prior_project_use_changes_role": bool(prior_use["matched_file_count"]),
        "evaluation_role": "project_exposed_external_development",
        "within_gse268273": within_target,
        "cross_source_expression_fingerprints": comparisons,
        "cross_source_matrix_sha256": matrix_hashes,
        "near_duplicate_source_ids": near_duplicate_sources,
        "fingerprint_method": "one-to-one GENCODE v49 genes; deposited voom values versus per-profile log1p-CPM; Pearson and tied-rank Spearman across genes",
        "fingerprint_salt_id_sha256": hashlib.sha256(salt).hexdigest(),
        "permissioned_random_salt_outside_release": True,
        "salted_expression_fingerprints_only": True,
        "raw_expression_fingerprints_retained": False,
        "noncomparable_registered_modalities": [
            "single_cell_or_single_nucleus_RNA_without_frozen_donor_pseudobulk",
            "RNA_microarray_without_platform_native_inductive_preprocessing",
            "protein_ATAC_histone_methylation_spatial_or_perturbation_assays",
        ],
        "champion_claim_eligible": False,
    }
    (audit_root / "audit.json").write_text(
        json.dumps(contamination, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if contamination_status != "passed":
        raise GSE268273FixtureError("GSE268273 contamination audit detected a duplicate")
    rights = {
        "schema_version": "masld-bench-gse268273-rights-v1",
        "access_tier": "public",
        "new_dua_or_controlled_access_required": False,
        "internal_nonclinical_research_use": "allowed_by_project_public_data_definition",
        "article_license": "CC_BY_NC_ND_4_0",
        "dataset_explicit_license_detected": False,
        "article_license_inherited_by_dataset": False,
        "raw_or_processed_data_redistribution": "prohibited_by_project_policy",
        "released_weight_redistribution": "requires_model_specific_terms_and_legal_review",
        "clinical_use": False,
    }
    (output / "rights_contract.json").write_text(
        json.dumps(rights, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {
        "schema_version": "masld-bench-gse268273-fibrosis-ood-fixture-v1",
        "status": "passed_fixture_scored_use_blocked_raw_reprocessing",
        "participants": expected_participants,
        "biological_unit": "participant",
        "geo_records_are_unique_participants": True,
        "technical_runs_are_replicates": False,
        "evaluation_role": "project_exposed_external_development",
        "project_sealed": False,
        "task_id": task.task_id,
        "model_input_path": "model_input",
        "evaluator_only_path": "evaluator_only",
        "preprocessing_input_path": "preprocessing_input",
        "contamination_audit_path": "contamination_audit",
        "model_input_contains_labels": False,
        "model_input_contains_expression_values": False,
        "deposited_voom_scored_input_allowed": False,
        "processed_differential_expression_opened": False,
        "external_fit_or_calibration_performed": False,
        "external_outcomes_scored": False,
        "scored_transfer_active": False,
        "next_gate": "materialize_and_requantify_824_raw_runs_then_apply_training_frozen_preprocessing",
    }
    (output / "receipt.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--activation", type=Path, required=True)
    parser.add_argument("--training-molecular", type=Path, required=True)
    parser.add_argument("--training-join", type=Path, required=True)
    parser.add_argument("--gse260666-activation", type=Path, required=True)
    parser.add_argument("--microarray-participants", type=Path, required=True)
    parser.add_argument("--microarray-overlap", type=Path, required=True)
    parser.add_argument("--task-spec", type=Path, required=True)
    parser.add_argument("--promotion-gate", type=Path, required=True)
    parser.add_argument("--expansion-config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--benchmark-root", type=Path, required=True)
    parser.add_argument("--fingerprint-salt-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    article_payload = download(ARTICLE_URL)
    ena_payload = download(ENA_URL)
    result = build_fixture(
        source=args.source,
        matrix=args.matrix,
        activation=args.activation,
        training_molecular=args.training_molecular,
        training_join=args.training_join,
        gse260666_activation=args.gse260666_activation,
        microarray_participants=args.microarray_participants,
        microarray_overlap=args.microarray_overlap,
        task_spec=args.task_spec,
        promotion_gate=args.promotion_gate,
        expansion_config=args.expansion_config,
        project_root=args.project_root,
        benchmark_root=args.benchmark_root,
        fingerprint_salt_file=args.fingerprint_salt_file,
        output=args.output,
        article_payload=article_payload,
        ena_payload=ena_payload,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
