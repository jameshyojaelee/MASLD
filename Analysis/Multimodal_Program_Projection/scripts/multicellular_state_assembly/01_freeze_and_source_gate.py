#!/usr/bin/env python3
"""Seal Plan 44 lineage hypotheses and audit public metadata without outcomes.

This script intentionally downloads GEO SOFT metadata only. It never downloads
or opens an expression matrix, and it must run before any Plan 44 atlas score is
calculated.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "multicellular-assembly-response-2026-08-09"
)
PLAN43 = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "chronic-state-risk-bridge-2026-08-09"
)
FROZEN = PLAN43 / "frozen_inputs"

PROGRAM_REGISTRY = FROZEN / "program_registry__program_registry_v2.tsv"
PROGRAM_MEMBERSHIP = FROZEN / "program_membership__program_membership_v2.tsv"
EVIDENCE_CLASSES = FROZEN / "evidence_classes__frozen_evidence_classes.tsv"
CORRECTED_LOCI = PLAN43 / "corrected_genetic_locus_registry.tsv"

SOFT_URLS = {
    "GSE168285": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/"
        "GSE168285/soft/GSE168285_family.soft.gz"
    ),
    "GSE175448": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/"
        "GSE175448/soft/GSE175448_family.soft.gz"
    ),
    "GSE238219": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/"
        "GSE238219/soft/GSE238219_family.soft.gz"
    ),
}

EXPECTED_LINEAGES = {
    "hepatocytes": 30,
    "fibroblasts": 29,
    "macrophages": 19,
    "cholangiocytes": 28,
    "tcells": 11,
}

SOURCE_TARGETS = (
    "C6orf106",
    "GPAM",
    "LYPLAL1",
    "NCKIPSD",
    "PNPLA3",
    "PPP1R3B",
    "RBM6",
    "TNKS",
    "TRIB1",
    "VKORC1",
    "WDR6",
)
SYMBOL_ALIASES = {"C6orf106": "ILRUN"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def download_once(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "MASLD-Plan44/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    # Validate the complete gzip payload before changing the candidate file.
    gzip.decompress(payload)
    if path.exists():
        if path.read_bytes() != payload:
            raise RuntimeError(f"public metadata changed after freeze: {path}")
        return
    descriptor, temporary_name = tempfile.mkstemp(dir=path.parent, prefix=".download-")
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def parse_soft(path: Path) -> list[dict[str, list[str]]]:
    records: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current = {"^SAMPLE": [line.split(" = ", 1)[1]]}
                records.append(current)
            elif current is not None and line.startswith("!") and " = " in line:
                key, value = line.split(" = ", 1)
                current.setdefault(key, []).append(value)
    return records


def one(record: dict[str, list[str]], key: str) -> str:
    values = record.get(key, [])
    if len(values) != 1:
        raise RuntimeError(f"expected one {key}, found {values!r}")
    return values[0]


def characteristic(record: dict[str, list[str]], prefix: str) -> str:
    values = [
        value.split(": ", 1)[1]
        for value in record.get("!Sample_characteristics_ch1", [])
        if value.startswith(prefix + ": ")
    ]
    if len(values) != 1:
        raise RuntimeError(f"expected one characteristic {prefix}, found {values!r}")
    return values[0]


def freeze_lineages() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    registry = read_tsv(PROGRAM_REGISTRY)
    if len(registry) != 117:
        raise RuntimeError(f"expected 117 programs, found {len(registry)}")
    counts = Counter(row["cell_type"] for row in registry)
    if dict(counts) != EXPECTED_LINEAGES:
        raise RuntimeError(f"lineage universe drift: {dict(counts)!r}")
    if any(
        row.get("external_outcomes_read", "").strip().lower()
        not in {"", "false", "0", "no"}
        for row in registry
    ):
        raise RuntimeError("program registry is not outcome-blind")

    denominators: dict[str, float] = defaultdict(float)
    for row in registry:
        denominators[row["cell_type"]] += abs(float(row["primary_beta"]))
    if any(value <= 0 for value in denominators.values()):
        raise RuntimeError(f"nonpositive lineage loading denominator: {denominators!r}")

    geometry: list[dict[str, object]] = []
    for row in sorted(registry, key=lambda item: (item["cell_type"], int(item["module"]))):
        beta = float(row["primary_beta"])
        geometry.append(
            {
                "program_uid": row["program_uid"],
                "cell_type": row["cell_type"],
                "module": row["module"],
                "module_name": row["module_name"],
                "stage_beta": f"{beta:.17g}",
                "within_lineage_loading": f"{beta / denominators[row['cell_type']]:.17g}",
                "membership_sha256": row["membership_sha256"],
                "program_release_id": row["release_id"],
                "external_outcomes_read_at_freeze": row.get("external_outcomes_read", "FALSE")
                or "FALSE",
            }
        )

    for lineage in EXPECTED_LINEAGES:
        total = sum(
            abs(float(row["within_lineage_loading"]))
            for row in geometry
            if row["cell_type"] == lineage
        )
        if not math.isclose(total, 1.0, rel_tol=0, abs_tol=1e-12):
            raise RuntimeError(f"lineage weights do not sum to one: {lineage}={total}")

    hypotheses = [
        {
            "hypothesis_id": "MPS_FAT_HEPATOCYTE_SPECIFICITY",
            "dataset": "GSE168285",
            "estimand": "fat x (hepatocyte - mean(fibroblast,macrophage))",
            "expected_direction": "positive",
            "multiplicity_family": "mps_three_primary_maxT",
            "promotion_role": "required",
        },
        {
            "hypothesis_id": "MPS_TGFB_REMODELING_SPECIFICITY",
            "dataset": "GSE168285",
            "estimand": "TGF_beta x (mean(fibroblast,macrophage) - hepatocyte)",
            "expected_direction": "positive",
            "multiplicity_family": "mps_three_primary_maxT",
            "promotion_role": "required",
        },
        {
            "hypothesis_id": "MPS_NPC_REMODELING_SPECIFICITY",
            "dataset": "GSE168285",
            "estimand": "high_NPC x (mean(fibroblast,macrophage) - hepatocyte)",
            "expected_direction": "positive",
            "multiplicity_family": "mps_three_primary_maxT",
            "promotion_role": "required",
        },
        {
            "hypothesis_id": "HUMAN_REMODELING_DISASSEMBLY",
            "dataset": "GSE175448",
            "estimand": "mean(delta(R-H)|improved)-mean(delta(R-H)|not_improved)",
            "expected_direction": "negative",
            "multiplicity_family": "single_primary_stratified_exact",
            "promotion_role": "required",
        },
        {
            "hypothesis_id": "HEPARG_OBSERVABILITY_BOUNDARY",
            "dataset": "GSE238219",
            "estimand": "absolute hepatocyte-axis effect - absolute remodeling-axis projection",
            "expected_direction": "positive",
            "multiplicity_family": "supplementary_complete_11_target_family",
            "promotion_role": "supplementary_only",
        },
    ]
    return geometry, hypotheses


def audit_gse168(records: list[dict[str, list[str]]]) -> dict[str, object]:
    if len(records) != 179:
        raise RuntimeError(f"GSE168285 sample drift: {len(records)}")
    combinations: Counter[tuple[str, str]] = Counter()
    experiments: set[str] = set()
    for record in records:
        title = one(record, "!Sample_title")
        experiment = title.split("_", 1)[0]
        treatment = characteristic(record, "treatment")
        experiments.add(experiment)
        combinations[(experiment, treatment)] += 1
    expected_experiments = {"NAFT18012", "NAFT18013", "NAFT18014", "NAFT18015"}
    if experiments != expected_experiments:
        raise RuntimeError(f"GSE168285 experiment drift: {sorted(experiments)!r}")
    if len(combinations) != 60:
        raise RuntimeError(f"GSE168285 condition-design drift: {len(combinations)}")
    replicate_counts = Counter(combinations.values())
    if replicate_counts != Counter({3: 59, 2: 1}):
        raise RuntimeError(f"GSE168285 replicate drift: {dict(replicate_counts)!r}")
    return {
        "dataset": "GSE168285",
        "source_gate": "mps_metadata_and_biological_unit",
        "status": "pass_model_specific_no_donor_key",
        "n_records": 179,
        "n_inferential_units": 60,
        "inferential_unit": "experiment_by_condition_mean",
        "n_human_donors": "unknown",
        "outcomes_accessed": "false",
        "details": "four experiment blocks; 59 triplicate and one duplicate condition combinations; donor-generalized inference prohibited",
    }


def audit_gse175(records: list[dict[str, list[str]]]) -> dict[str, object]:
    if len(records) != 38:
        raise RuntimeError(f"GSE175448 sample drift: {len(records)}")
    participants: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for record in records:
        description = one(record, "!Sample_description")
        participant, source_arm, *_ = description.split("_")
        arm = "cenicriviroc" if source_arm == "drugA" else "placebo"
        timing = characteristic(record, "timing of biopsy")
        improvement = characteristic(record, "fibrosis improvement")
        participants[participant].append((arm, timing, improvement))
    if len(participants) != 19:
        raise RuntimeError(f"GSE175448 participant drift: {len(participants)}")
    arm_counts: Counter[str] = Counter()
    improvement_counts: Counter[tuple[str, str]] = Counter()
    for participant, rows in participants.items():
        if len(rows) != 2:
            raise RuntimeError(f"GSE175448 incomplete pair: {participant}={rows!r}")
        arms = {row[0] for row in rows}
        timings = {row[1] for row in rows}
        if len(arms) != 1 or timings != {"pre-treatment", "post-treatment"}:
            raise RuntimeError(f"GSE175448 invalid pair: {participant}={rows!r}")
        arm = next(iter(arms))
        post = [row for row in rows if row[1] == "post-treatment"][0]
        arm_counts[arm] += 1
        improvement_counts[(arm, post[2])] += 1
    if arm_counts != Counter({"placebo": 10, "cenicriviroc": 9}):
        raise RuntimeError(f"GSE175448 arm drift: {dict(arm_counts)!r}")
    expected = Counter(
        {
            ("cenicriviroc", "improved"): 4,
            ("cenicriviroc", "not_improved"): 5,
            ("placebo", "improved"): 3,
            ("placebo", "not_improved"): 7,
        }
    )
    if improvement_counts != expected:
        raise RuntimeError(f"GSE175448 response drift: {dict(improvement_counts)!r}")
    return {
        "dataset": "GSE175448",
        "source_gate": "authoritative_participant_pairs",
        "status": "pass_authoritative_pairs",
        "n_records": 38,
        "n_inferential_units": 19,
        "inferential_unit": "participant_pair",
        "n_human_donors": 19,
        "outcomes_accessed": "false",
        "details": "9 cenicriviroc/10 placebo; 4/9 and 3/10 improved; 15120 stratified exact allocations",
    }


def audit_gse238(
    records: list[dict[str, list[str]]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    if len(records) != 5:
        raise RuntimeError(f"GSE238219 library drift: {len(records)}")
    evidence = {row["gene_symbol"]: row for row in read_tsv(EVIDENCE_CLASSES)}
    loci = {row["gene_symbol"]: row for row in read_tsv(CORRECTED_LOCI)}
    targets: list[dict[str, object]] = []
    for source_symbol in SOURCE_TARGETS:
        canonical_symbol = SYMBOL_ALIASES.get(source_symbol, source_symbol)
        evidence_row = evidence.get(canonical_symbol, {})
        locus_row = loci.get(canonical_symbol, {})
        targets.append(
            {
                "source_symbol": source_symbol,
                "canonical_symbol": canonical_symbol,
                "alias_resolution": (
                    f"{source_symbol}->{canonical_symbol}"
                    if source_symbol != canonical_symbol
                    else "identity"
                ),
                "evidence_class": evidence_row.get("static_class", "not_in_atlas"),
                "in_corrected_locus_registry": str(bool(locus_row)).lower(),
                "is_corrected_locus_representative": locus_row.get("is_representative", ""),
                "any_direct_tier1_support": locus_row.get("any_direct_tier1_support", ""),
                "tier12_driving_trait": locus_row.get("tier12_driving_trait", ""),
            }
        )
    genetic_only = [row for row in targets if row["evidence_class"] == "genetic_only"]
    if {row["canonical_symbol"] for row in genetic_only} != {
        "ILRUN",
        "GPAM",
        "PPP1R3B",
        "TNKS",
    }:
        raise RuntimeError(f"GSE238219 genetic-only overlap drift: {genetic_only!r}")
    direct = [row for row in genetic_only if row["any_direct_tier1_support"] == "true"]
    if [row["canonical_symbol"] for row in direct] != ["GPAM"]:
        raise RuntimeError(f"GSE238219 direct-support overlap drift: {direct!r}")
    gate = {
        "dataset": "GSE238219",
        "source_gate": "supplementary_perturbseq_observability",
        "status": "pass_supplementary_only",
        "n_records": 5,
        "n_inferential_units": 5,
        "inferential_unit": "source_replicate_library",
        "n_human_donors": "not_applicable_cell_line",
        "outcomes_accessed": "false",
        "details": "11 targets; four corrected genetic-only after C6orf106-to-ILRUN alias resolution; GPAM alone has direct Tier-1 support",
    }
    return gate, targets


def main() -> None:
    for required in (
        PROGRAM_REGISTRY,
        PROGRAM_MEMBERSHIP,
        EVIDENCE_CLASSES,
        CORRECTED_LOCI,
    ):
        if not required.is_file():
            raise FileNotFoundError(required)

    source_dir = CANDIDATE / "source_metadata"
    source_manifest: list[dict[str, object]] = []
    soft_records: dict[str, list[dict[str, list[str]]]] = {}
    for accession, url in SOFT_URLS.items():
        path = source_dir / f"{accession}_family.soft.gz"
        download_once(url, path)
        soft_records[accession] = parse_soft(path)
        source_manifest.append(
            {
                "source_id": accession,
                "source_url": url,
                "relative_path": path.relative_to(CANDIDATE).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "content_class": "metadata_only_no_expression_outcome",
                "retrieval_date": "2026-08-09",
            }
        )

    geometry, hypotheses = freeze_lineages()
    gse168_gate = audit_gse168(soft_records["GSE168285"])
    gse175_gate = audit_gse175(soft_records["GSE175448"])
    gse238_gate, target_registry = audit_gse238(soft_records["GSE238219"])

    write_tsv(
        CANDIDATE / "frozen_lineage_geometry.tsv",
        geometry,
        [
            "program_uid",
            "cell_type",
            "module",
            "module_name",
            "stage_beta",
            "within_lineage_loading",
            "membership_sha256",
            "program_release_id",
            "external_outcomes_read_at_freeze",
        ],
    )
    write_tsv(
        CANDIDATE / "frozen_hypotheses.tsv",
        hypotheses,
        [
            "hypothesis_id",
            "dataset",
            "estimand",
            "expected_direction",
            "multiplicity_family",
            "promotion_role",
        ],
    )
    write_tsv(
        CANDIDATE / "gse238219_target_registry.tsv",
        target_registry,
        [
            "source_symbol",
            "canonical_symbol",
            "alias_resolution",
            "evidence_class",
            "in_corrected_locus_registry",
            "is_corrected_locus_representative",
            "any_direct_tier1_support",
            "tier12_driving_trait",
        ],
    )
    gates = [gse168_gate, gse175_gate, gse238_gate]
    write_tsv(
        CANDIDATE / "source_gate_status.tsv",
        gates,
        [
            "dataset",
            "source_gate",
            "status",
            "n_records",
            "n_inferential_units",
            "inferential_unit",
            "n_human_donors",
            "outcomes_accessed",
            "details",
        ],
    )
    write_tsv(
        CANDIDATE / "source_manifest.tsv",
        source_manifest,
        [
            "source_id",
            "source_url",
            "relative_path",
            "bytes",
            "sha256",
            "content_class",
            "retrieval_date",
        ],
    )

    known_findings = [
        {
            "source_id": "PLAN41",
            "known_before_plan44": "true",
            "finding": "PCLS culture time increased the two externally frozen programs in 5/6 donors whereas acute combined lipid did not pass",
            "ownership": "prior_local_result_not_plan44_confirmation",
        },
        {
            "source_id": "PLAN43",
            "known_before_plan44": "true",
            "finding": "the frozen whole-registry global axis failed external human reversal",
            "ownership": "terminal_failed_predecessor",
        },
        {
            "source_id": "GSE168285_SOURCE",
            "known_before_plan44": "true",
            "finding": "the source reports that fat, TGF-beta, and nonparenchymal context shape its NASH/fibrosis phenotype",
            "ownership": "source_owned_general_finding",
        },
        {
            "source_id": "GSE175448_SOURCE",
            "known_before_plan44": "true",
            "finding": "the source reports suppression of fibrogenic pathways among fibrosis improvers",
            "ownership": "source_owned_general_finding",
        },
        {
            "source_id": "GSE238219_SOURCE",
            "known_before_plan44": "true",
            "finding": "the source reports target-specific lipid and transcriptional effects in HepaRG",
            "ownership": "source_owned_target_findings",
        },
    ]
    write_tsv(
        CANDIDATE / "source_known_findings.tsv",
        known_findings,
        ["source_id", "known_before_plan44", "finding", "ownership"],
    )

    upstream_rows = []
    for path in (
        PROGRAM_REGISTRY,
        PROGRAM_MEMBERSHIP,
        EVIDENCE_CLASSES,
        CORRECTED_LOCI,
        PLAN43 / "release_manifest.tsv",
        PLAN43 / "VALIDATED.json",
    ):
        if path.is_file():
            upstream_rows.append(
                {
                    "path": path.relative_to(ROOT).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    write_tsv(
        CANDIDATE / "upstream_immutability_baseline.tsv",
        upstream_rows,
        ["path", "bytes", "sha256"],
    )

    validation_rows = [
        {"check_id": "program_count", "status": "pass", "observed": 117},
        {
            "check_id": "lineage_counts",
            "status": "pass",
            "observed": json.dumps(EXPECTED_LINEAGES, sort_keys=True),
        },
        {"check_id": "gse168_samples", "status": "pass", "observed": 179},
        {"check_id": "gse168_condition_means", "status": "pass", "observed": 60},
        {"check_id": "gse175_pairs", "status": "pass", "observed": 19},
        {"check_id": "gse175_improved", "status": "pass", "observed": 7},
        {"check_id": "gse238_libraries", "status": "pass", "observed": 5},
        {"check_id": "gse238_targets", "status": "pass", "observed": 11},
        {
            "check_id": "new_expression_outcomes_accessed",
            "status": "pass",
            "observed": "false",
        },
    ]
    write_tsv(
        CANDIDATE / "preflight_validation.tsv",
        validation_rows,
        ["check_id", "status", "observed"],
    )

    active_paths = sorted(
        path
        for path in CANDIDATE.rglob("*")
        if path.is_file() and path.name not in {"release_manifest.tsv", "PREFLIGHT_READY.json"}
    )
    release_rows = [
        {
            "relative_path": path.relative_to(CANDIDATE).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in active_paths
    ]
    write_tsv(
        CANDIDATE / "release_manifest.tsv",
        release_rows,
        ["relative_path", "bytes", "sha256"],
    )
    ready = {
        "status": "PREFLIGHT_READY_OUTCOMES_UNOPENED",
        "release_id": "multicellular-assembly-response-2026-08-09",
        "n_programs": 117,
        "lineage_counts": EXPECTED_LINEAGES,
        "gse168285": {"samples": 179, "condition_experiment_means": 60},
        "gse175448": {"profiles": 38, "participants": 19, "improved": 7},
        "gse238219": {
            "libraries": 5,
            "targets": 11,
            "corrected_genetic_only_targets": 4,
        },
        "new_expression_outcomes_accessed": False,
        "release_manifest_sha256": sha256(CANDIDATE / "release_manifest.tsv"),
    }
    (CANDIDATE / "PREFLIGHT_READY.json").write_text(
        json.dumps(ready, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "PLAN44_PREFLIGHT_READY\t"
        f"programs=117\tgse168_samples=179\tgse175_pairs=19\t"
        f"gse238_genetic_only=4\tmanifest={ready['release_manifest_sha256']}"
    )


if __name__ == "__main__":
    main()
