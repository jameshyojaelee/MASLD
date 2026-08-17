#!/usr/bin/env python3
"""Independent fail-closed validation for the non-genetic ATAC v3 package."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import scipy.io


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
PROGRAM_RELEASE_ID = "program-context-v2-candidate-2026-08-07"
STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
PROJECT = Path(__file__).resolve().parents[4]
EXPECTED_ROOT = (
    PROJECT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()
PROGRAM_ROOT = (
    PROJECT / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "hotspot"
)


class ValidationFailure(RuntimeError):
    pass


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=EXPECTED_ROOT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=lambda index: (values[index], index))
    result = [1.0] * n
    running = 1.0
    for offset in range(n - 1, -1, -1):
        index = order[offset]
        running = min(running, values[index] * n / (offset + 1))
        result[index] = min(1.0, running)
    return result


def assert_bh(rows: list[dict[str, str]], family_fields: tuple[str, ...], label: str) -> int:
    families: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("pvalue", "") not in {"", "NA", "NaN"}:
            families[tuple(row[field] for field in family_fields)].append(row)
    total = 0
    for family, members in families.items():
        expected = bh([float(row["pvalue"]) for row in members])
        for row, value in zip(members, expected):
            if not math.isclose(float(row["qvalue"]), value, rel_tol=1e-8, abs_tol=1e-10):
                raise ValidationFailure(f"BH mismatch in {label} family {family}")
        total += len(members)
    return total


def validate_inputs(root: Path) -> dict[str, object]:
    identity = json.loads((root / "candidate_identity.json").read_text(encoding="utf-8"))
    if identity["release_id"] != RELEASE_ID or identity["seed"] != 42:
        raise ValidationFailure("candidate identity drift")
    manifest = read_tsv(root / "input_manifest.tsv")
    if not manifest or any(row["sha256"] == "deferred_large_file_hash" for row in manifest):
        raise ValidationFailure("input manifest contains deferred or absent hashes")
    for row in manifest:
        path = Path(row["relative_or_absolute_path"])
        if not path.is_absolute():
            path = PROJECT / path
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise ValidationFailure(f"frozen input hash mismatch: {path}")
    blindness = read_tsv(root / "condition_blindness/condition_blindness_audit.tsv")
    if len(blindness) != 10 or any(row["pass"] != "TRUE" for row in blindness):
        raise ValidationFailure("condition-label permutation audit failed or is incomplete")
    return {"n_frozen_inputs": len(manifest), "n_condition_blindness_checks": len(blindness)}


def validate_peaks_and_counts(root: Path) -> dict[str, object]:
    peaks = read_tsv(root / "consensus_peak_manifest.tsv")
    support = read_tsv(root / "peaks/native_peak_support.tsv")
    support_families: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in support:
        family = (row["cohort"], row["lineage"])
        support_families[family].append(row)
        retained = row["retained"] == "TRUE"
        n_support = int(row["n_supporting_donors"])
        if int(row["min_supporting_donors"]) != 2:
            raise ValidationFailure(f"native peak donor-support threshold drift: {family}")
        if retained != (n_support >= 2):
            raise ValidationFailure(f"native peak retention/support mismatch: {family}")
        donors = [value for value in row["supporting_donors"].split(",") if value]
        if len(donors) != n_support or len(set(donors)) != len(donors):
            raise ValidationFailure(f"native peak supporting-donor audit mismatch: {family}")
        if int(row["eligible_donors"]) < 2 or n_support > int(row["eligible_donors"]):
            raise ValidationFailure(f"invalid native peak donor denominator: {family}")
    expected_families = {
        (cohort, lineage)
        for cohort in ("GSE244832", "GSE281367")
        for lineage in ("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")
    }
    if set(support_families) != expected_families:
        raise ValidationFailure("native peak support audit is missing cohort-lineage families")
    for family, rows in support_families.items():
        cohort, lineage = family
        retained_coordinates = {
            (row["chrom"], int(row["start0"]), int(row["end"]))
            for row in rows if row["retained"] == "TRUE"
        }
        native = read_tsv(root / "peaks/native" / f"{cohort}.{lineage}.tsv")
        native_coordinates = {
            (row["chrom"], int(row["start0"]), int(row["end"])) for row in native
        }
        if not retained_coordinates or native_coordinates != retained_coordinates:
            raise ValidationFailure(f"native peak table/support audit mismatch: {family}")
    coordinates = set()
    by_lineage = Counter()
    expected_by_lineage: dict[str, list[str]] = defaultdict(list)
    blacklist_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed").open(
        "r", encoding="utf-8"
    ) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in STANDARD_CHROMS:
                blacklist_by_chrom[fields[0]].append((int(fields[1]), int(fields[2])))
    for row in peaks:
        coordinate = (row["lineage"], row["chrom"], int(row["start0"]), int(row["end"]))
        if coordinate in coordinates:
            raise ValidationFailure(f"duplicate consensus coordinate: {coordinate}")
        coordinates.add(coordinate)
        if row["chrom"] not in STANDARD_CHROMS or int(row["width"]) != 500:
            raise ValidationFailure(f"invalid consensus coordinate: {coordinate}")
        if row["blacklist_overlap"] != "FALSE":
            raise ValidationFailure(f"blacklisted consensus coordinate: {coordinate}")
        if any(int(row["start0"]) < end and int(row["end"]) > start for start, end in blacklist_by_chrom[row["chrom"]]):
            raise ValidationFailure(f"independent blacklist overlap: {coordinate}")
        if row["supported_both"] != str(
            row["supported_gse244832"] == "TRUE" and row["supported_gse281367"] == "TRUE"
        ).upper():
            raise ValidationFailure("consensus provenance flag mismatch")
        by_lineage[row["lineage"]] += 1
        expected_by_lineage[row["lineage"]].append(
            f'{row["chrom"]}:{row["start0"]}-{row["end"]}'
        )
    for lineage, values in expected_by_lineage.items():
        parsed = sorted(
            (
                value.split(":")[0],
                int(value.split(":")[1].split("-")[0]),
                int(value.rsplit("-", 1)[1]),
            )
            for value in values
        )
        for previous, current in zip(parsed, parsed[1:]):
            if previous[0] == current[0] and current[1] < previous[2]:
                raise ValidationFailure(f"overlapping consensus peaks within {lineage}")
    matrix_count = 0
    for cohort in ("GSE244832", "GSE281367"):
        for lineage in ("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"):
            count_root = root / "counts" / cohort
            donors = read_tsv(count_root / f"{lineage}.donors.tsv")
            peak_table = read_tsv(count_root / f"{lineage}.peaks.tsv")
            with gzip.open(count_root / f"{lineage}.mtx.gz", "rb") as handle:
                matrix = scipy.io.mmread(handle)
            if matrix.shape != (len(donors), len(peak_table)):
                raise ValidationFailure(f"count matrix dimension mismatch: {cohort}/{lineage}")
            observed_peaks = [row["peak_coordinate"] for row in peak_table]
            if observed_peaks != expected_by_lineage[lineage]:
                raise ValidationFailure(f"count peak order drift: {cohort}/{lineage}")
            if len({row["donor_id"] for row in donors}) != len(donors):
                raise ValidationFailure(f"duplicated donor row: {cohort}/{lineage}")
            row_sums = matrix.tocsr().sum(axis=1).A1
            for row, total in zip(donors, row_sums):
                if int(row["total_counted_fragments"]) != int(round(total)):
                    raise ValidationFailure(f"library total mismatch: {cohort}/{lineage}/{row['donor_id']}")
            matrix_count += 1
    for cohort in ("GSE244832", "GSE281367"):
        for lineage in ("hepatocyte", "stellate", "macrophage"):
            count_root = root / "counts_legacy_gse244832" / cohort
            donors = read_tsv(count_root / f"{lineage}.donors.tsv")
            peak_table = read_tsv(count_root / f"{lineage}.peaks.tsv")
            with gzip.open(count_root / f"{lineage}.mtx.gz", "rb") as handle:
                matrix = scipy.io.mmread(handle)
            if matrix.shape != (len(donors), len(peak_table)):
                raise ValidationFailure(f"legacy count dimension mismatch: {cohort}/{lineage}")
            matrix_count += 1
    recount = read_tsv(root / "recount/exact_fragment_recount.tsv")
    if len(recount) != 10 or any(
        row.get("counting_unit") != "deduplicated_fragment_record"
        or row.get("record_count_match") != "TRUE"
        for row in recount
    ):
        raise ValidationFailure("seed-fixed direct fragment recount is incomplete or failed")
    return {
        "n_consensus_peaks": len(peaks), "n_count_matrices": matrix_count,
        "n_exact_recounts": len(recount), "n_native_peak_support_rows": len(support),
        "peaks_by_lineage": dict(by_lineage),
    }


def validate_da(root: Path) -> dict[str, object]:
    da_root = root / "da"
    total = 0
    for cohort in ("GSE244832", "GSE281367"):
        for lineage in ("hepatocyte", "stellate", "macrophage"):
            rows = read_tsv(da_root / f"cohort_primary_{cohort}_{lineage}.tsv")
            assert_bh(rows, tuple(), f"DA {cohort}/{lineage}")
            total += len(rows)
    joined = read_tsv(da_root / "da_peak_results.tsv.gz")
    valid = {"supported", "discordant", "source_dependent", "indeterminate", "untestable"}
    if any(row["evidence_state"] not in valid for row in joined):
        raise ValidationFailure("invalid DA evidence state")
    summary = read_tsv(da_root / "da_lineage_summary.tsv")
    p_rows = [row for row in summary if row["directional_binomial_pvalue"] not in {"", "NA", "NaN"}]
    expected = bh([float(row["directional_binomial_pvalue"]) for row in p_rows])
    for row, value in zip(p_rows, expected):
        if not math.isclose(float(row["directional_binomial_qvalue"]), value, rel_tol=1e-8, abs_tol=1e-10):
            raise ValidationFailure("lineage directional BH mismatch")
    if len(summary) != 3:
        raise ValidationFailure("DA lineage family must contain three lineages")
    untestable = [row for row in summary if row["lineage_state"] == "untestable"]
    for row in untestable:
        if int(row["n_jointly_testable"]) != 0 or not row["untestable_reason"]:
            raise ValidationFailure("untestable DA lineage lacks an explicit zero-denominator reason")
    if any(row["lineage_state"] not in {"testable", "untestable"} for row in summary):
        raise ValidationFailure("invalid DA lineage testability state")
    model_qc = read_tsv(da_root / "da_model_qc.tsv")
    expected_model_rows = sum(
        int(row["n_normal_gse244832"]) + int(row["n_mash_gse244832"]) +
        int(row["n_normal_gse281367"]) + int(row["n_mash_gse281367"])
        for row in summary if row["lineage_state"] == "testable"
    )
    if len(model_qc) != expected_model_rows:
        raise ValidationFailure(
            f"expected {expected_model_rows} testable-lineage DA donor-model rows, observed {len(model_qc)}"
        )
    expected_families = {
        (cohort, row["lineage"])
        for row in summary if row["lineage_state"] == "testable"
        for cohort in ("GSE244832", "GSE281367")
    }
    observed_families = {(row["cohort"], row["lineage"]) for row in model_qc}
    if observed_families != expected_families:
        raise ValidationFailure("DA model families do not match testable lineage states")
    if any(row["condition"] == "MASL" for row in model_qc):
        raise ValidationFailure("GSE244832 MASL donor entered a primary DA model")
    if any(row["design_rank"] != row["n_model_columns"] for row in model_qc):
        raise ValidationFailure("DA design is not full rank")
    if any(row["sample_quality_weight"] in {"", "NA", "NaN"} for row in model_qc):
        raise ValidationFailure("DA sample-quality weight is missing")
    legacy = read_tsv(da_root / "legacy_gse244832_peak_space_peak_results.tsv.gz")
    for lineage in ("hepatocyte", "stellate", "macrophage"):
        family = [row for row in legacy if row["lineage"] == lineage]
        for suffix in ("gse244832", "gse281367"):
            expected = bh([float(row[f"pvalue_{suffix}"]) for row in family])
            for row, value in zip(family, expected):
                if not math.isclose(float(row[f"qvalue_{suffix}"]), value, rel_tol=1e-8, abs_tol=1e-10):
                    raise ValidationFailure(f"legacy peak-space BH mismatch: {lineage}/{suffix}")
    return {
        "n_joined_da_peaks": len(joined), "n_cohort_da_tests": total,
        "n_da_model_samples": len(model_qc), "n_legacy_da_peaks": len(legacy),
    }


def validate_programs(root: Path) -> dict[str, object]:
    rows = read_tsv(root / "programs/program_atac_results.tsv")
    if len(rows) != 234:
        raise ValidationFailure(f"program output must have 234 rows, observed {len(rows)}")
    counts = Counter(row["cohort"] for row in rows)
    if counts != {"GSE244832": 117, "GSE281367": 117}:
        raise ValidationFailure(f"program cohort completeness mismatch: {counts}")
    registry = read_tsv(PROGRAM_ROOT / "program_registry_v2.tsv")
    expected_ids = {row["program_uid"] for row in registry}
    expected_hashes = {row["program_uid"]: row["membership_sha256"] for row in registry}
    expected_names = {row["program_uid"]: row["module_name"] for row in registry}
    for cohort in counts:
        observed = {row["program_uid"] for row in rows if row["cohort"] == cohort}
        if observed != expected_ids:
            raise ValidationFailure(f"program registry mismatch in {cohort}")
    if any(row["membership_sha256"] != expected_hashes[row["program_uid"]] for row in rows):
        raise ValidationFailure("frozen program membership hash drift")
    if any(row["program_name"] != expected_names[row["program_uid"]] for row in rows):
        raise ValidationFailure("frozen program name drift")
    assert_bh(rows, ("cohort",), "program")
    coverage = read_tsv(root / "programs/program_measurement_coverage.tsv")
    if len(coverage) != 234:
        raise ValidationFailure("program coverage must contain 234 rows")
    model_qc = read_tsv(root / "programs/program_model_qc.tsv")
    if len(model_qc) != 10:
        raise ValidationFailure("program model QC must contain ten cohort-lineage rows")
    for row in model_qc:
        if row["contrast_testable"] == "TRUE" and row["design_rank"] != row["n_model_columns"]:
            raise ValidationFailure("program design is not full rank")
    return {"n_program_rows": len(rows), "n_program_models": len(model_qc)}


def validate_chromvar(root: Path) -> dict[str, object]:
    chromvar = root / "chromvar"
    session = (chromvar / "sessionInfo.txt").read_text(encoding="utf-8")
    required = ("chromVAR", "motifmatchr", "TFBSTools", "JASPAR2024", "BSgenome.Hsapiens.UCSC.hg38")
    if any(value not in session for value in required):
        raise ValidationFailure("official chromVAR session is incomplete")
    lock = root / "environments/chromvar.explicit.txt"
    if not lock.is_file():
        raise ValidationFailure("chromVAR explicit environment lock is missing")
    rows = read_tsv(chromvar / "chromvar_results.tsv.gz")
    motif_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        motif_sets[(row["cohort"], row["lineage"])].add(row["motif_id"])
    universes = list(motif_sets.values())
    if len(universes) != 10 or any(value != universes[0] for value in universes[1:]):
        raise ValidationFailure("JASPAR2024 motif universe differs across cohort-lineage families")
    tested = assert_bh(rows, ("cohort", "lineage"), "chromVAR")
    model_qc = read_tsv(chromvar / "chromvar_model_qc.tsv")
    if len(model_qc) != 10:
        raise ValidationFailure("chromVAR model QC must contain ten cohort-lineage rows")
    for row in model_qc:
        if row["contrast_testable"] == "TRUE" and row["design_rank"] != row["n_model_columns"]:
            raise ValidationFailure("chromVAR design is not full rank")
    peak_filter = read_tsv(chromvar / "chromvar_peak_filter_qc.tsv")
    if len(peak_filter) != 10:
        raise ValidationFailure("chromVAR peak-filter QC must contain ten cohort-lineage rows")
    for row in peak_filter:
        if row["filter_method"] != "condition_blind_nonzero_total_donor_count":
            raise ValidationFailure("chromVAR peak filter is not the frozen condition-blind rule")
        if int(row["n_chromvar_input_peaks"]) <= 0:
            raise ValidationFailure("chromVAR lineage retained no measurable peaks")
        if int(row["n_consensus_peaks"]) != (
            int(row["n_zero_library_peaks_excluded"]) + int(row["n_chromvar_input_peaks"])
        ):
            raise ValidationFailure("chromVAR peak-filter denominator mismatch")
    return {
        "n_chromvar_rows": len(rows), "n_chromvar_tests": tested,
        "n_chromvar_models": len(model_qc), "n_chromvar_peak_filter_families": len(peak_filter),
    }


def validate_figures(root: Path) -> dict[str, object]:
    manifest = read_tsv(root / "figures/figure_manifest.tsv")
    if len(manifest) != 4:
        raise ValidationFailure("expected four candidate-only non-genetic PDFs")
    for row in manifest:
        pdf = root / row["pdf"]
        source = root / row["source_table"]
        if sha256(pdf) != row["pdf_sha256"] or sha256(source) != row["source_sha256"]:
            raise ValidationFailure(f"figure/source checksum mismatch: {row['panel']}")
        info = subprocess.run(["pdfinfo", str(pdf)], text=True, capture_output=True, check=True).stdout
        if "Pages:           1" not in info or "Page size:" not in info:
            raise ValidationFailure(f"PDF page count or dimensions invalid: {row['panel']}")
        pdffonts = shutil.which("pdffonts")
        if pdffonts:
            fonts = subprocess.run(
                [pdffonts, str(pdf)], text=True, capture_output=True, check=True
            ).stdout
            has_type3 = "Type 3" in fonts or "Type3" in fonts
        else:
            raw_pdf = pdf.read_bytes()
            font_subtypes = re.findall(rb"/Subtype\s*/([A-Za-z0-9]+)", raw_pdf)
            if not font_subtypes:
                raise ValidationFailure(
                    f"cannot inspect PDF fonts without pdffonts: {row['panel']}"
                )
            has_type3 = b"Type3" in font_subtypes
        if has_type3:
            raise ValidationFailure(f"PDF contains a Type 3 font: {row['panel']}")
    return {"n_candidate_pdfs": len(manifest)}


def atomic_table(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise ValidationFailure(f"refusing to overwrite validator output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, encoding="utf-8", newline="") as handle:
        temp = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def main() -> None:
    args = arguments()
    root = args.candidate_root.resolve()
    if root != EXPECTED_ROOT:
        raise ValidationFailure(f"unsafe candidate root: {root}")
    validation = root / "validation"
    if validation.exists() or (root / "NON_GENETIC_READY").exists():
        raise ValidationFailure("refusing to overwrite validation or readiness seal")
    checks = {
        "inputs": validate_inputs(root),
        "peaks_and_counts": validate_peaks_and_counts(root),
        "da": validate_da(root),
        "programs": validate_programs(root),
        "chromvar": validate_chromvar(root),
        "figures": validate_figures(root),
    }
    rows = [
        {"release_id": RELEASE_ID, "check": name, "status": "PASS", "details_json": json.dumps(value, sort_keys=True)}
        for name, value in checks.items()
    ]
    atomic_table(validation / "validation_results.tsv", ("release_id", "check", "status", "details_json"), rows)
    hashed_outputs = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        relative = path.relative_to(root).as_posix()
        if relative in {"output_manifest.tsv", "NON_GENETIC_READY"} or relative.startswith("genetics/"):
            continue
        hashed_outputs.append({
            "release_id": RELEASE_ID,
            "artifact": relative,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    atomic_table(
        root / "output_manifest.tsv",
        ("release_id", "artifact", "bytes", "sha256"),
        hashed_outputs,
    )
    artifacts = (
        root / "consensus_peak_manifest.tsv",
        root / "da/da_peak_results.tsv.gz",
        root / "programs/program_atac_results.tsv",
        root / "chromvar/chromvar_results.tsv.gz",
        root / "figures/figure_manifest.tsv",
        validation / "validation_results.tsv",
        root / "output_manifest.tsv",
    )
    seal_rows = [{
        "release_id": RELEASE_ID,
        "gate": "NON_GENETIC_READY",
        "status": "READY",
        "artifact": path.relative_to(root).as_posix(),
        "sha256": sha256(path),
    } for path in artifacts]
    atomic_table(root / "NON_GENETIC_READY", ("release_id", "gate", "status", "artifact", "sha256"), seal_rows)
    print("NON_GENETIC_READY")


if __name__ == "__main__":
    main()
