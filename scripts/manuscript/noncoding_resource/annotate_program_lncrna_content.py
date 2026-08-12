#!/usr/bin/env python3
"""Annotate the frozen Hotspot registry with GENCODE lncRNA content.

This script is deliberately annotation-only. It copies every frozen membership
field unchanged, adds stable-ID/biotype mapping fields, and summarizes mapped
lncRNA membership by program. It never rescales weights or computes scores.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


MEMBERSHIP_REQUIRED = (
    "cell_type",
    "module",
    "source_gene",
    "canonical_gene",
    "source_weight",
    "mapped_symbol",
    "mapped_symbol_status",
    "original_l1_weight",
    "canonical_weight_text",
    "membership_sha256",
    "program_uid",
)
IDENTITY_REQUIRED = (
    "annotation_release",
    "gene_id_versioned",
    "gene_id_base",
    "gene_version",
    "gene_name",
    "gene_type",
    "chromosome",
    "start_1based",
    "end_1based",
    "strand",
    "source",
    "level",
    "is_canonical_chromosome",
    "n_versions_for_base_id",
    "base_id_mapping_status",
    "n_gene_ids_for_symbol",
    "symbol_mapping_status",
)
MAPPED_STATUSES = frozenset(
    {
        "gencode_v49_unique_symbol_confirmed",
        "gencode_v49_unambiguous_ensembl_to_symbol",
    }
)
UNMAPPED_STATUS = "unmapped_or_ambiguous"
ENSG_BASE = re.compile(r"^ENSG[0-9]+$")
FROZEN_MEMBERSHIP_SHA256 = (
    "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b"
)
ANNOTATION_RELEASE = "GENCODE v49"
OUTPUT_FILES = (
    "program_membership_lncrna_annotation.tsv",
    "program_lncrna_content.tsv",
    "program_lncrna_mapping_audit.tsv",
    "source_manifest.tsv",
    "execution_manifest.json",
)


class AnnotationError(RuntimeError):
    """Raised when a frozen input or identity contract is violated."""


@dataclass(frozen=True)
class Identity:
    gene_id: str
    ensembl_base: str
    gene_name: str
    gene_biotype: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AnnotationError(message)


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(
    path: Path, required: Iterable[str]
) -> tuple[list[str], list[dict[str, str]]]:
    require(path.is_file(), f"missing input: {path}")
    require(not path.is_symlink(), f"symlinked input is prohibited: {path}")
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = [field for field in required if field not in fields]
        require(not missing, f"{path} is missing columns: {missing}")
        rows = list(reader)
    require(rows, f"empty input: {path}")
    return fields, rows


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_positive_float(value: str, label: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise AnnotationError(f"non-numeric {label}: {value!r}") from error
    require(
        math.isfinite(number) and number > 0.0,
        f"non-positive or non-finite {label}: {value!r}",
    )
    return number


def normalize_base(value: str) -> str:
    return re.sub(r"\.[0-9]+$", "", value.strip())


def build_identity_indexes(
    rows: list[dict[str, str]],
) -> tuple[dict[str, Identity], dict[str, list[Identity]]]:
    by_base: dict[str, Identity] = {}
    by_symbol: dict[str, list[Identity]] = defaultdict(list)
    versioned_seen: set[str] = set()
    releases: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        release = row["annotation_release"].strip()
        gene_id = row["gene_id_versioned"].strip()
        base = normalize_base(row["gene_id_base"])
        symbol = row["gene_name"].strip()
        biotype = row["gene_type"].strip()
        require(release != "", f"blank annotation_release at identity row {row_number}")
        require(
            release == ANNOTATION_RELEASE,
            f"identity annotation release drift at row {row_number}: {release!r}",
        )
        require(gene_id != "", f"blank gene_id_versioned at identity row {row_number}")
        require(base != "", f"blank gene_id_base at identity row {row_number}")
        require(symbol != "", f"blank gene_name at identity row {row_number}")
        require(biotype != "", f"blank gene_type at identity row {row_number}")
        require(
            normalize_base(gene_id) == base,
            f"gene_id/base mismatch at row {row_number}",
        )
        require(
            gene_id not in versioned_seen, f"duplicate versioned gene_id: {gene_id}"
        )
        require(base not in by_base, f"duplicate ensembl_base: {base}")
        require(
            row["n_versions_for_base_id"] == "1"
            and row["base_id_mapping_status"] == "unique",
            f"base-ID compatibility mapping is not unique for {base}",
        )
        releases.add(release)
        versioned_seen.add(gene_id)
        identity = Identity(gene_id, base, symbol, biotype)
        by_base[base] = identity
        by_symbol[symbol].append(identity)
    require(
        len(releases) == 1,
        f"identity table mixes annotation releases: {sorted(releases)}",
    )
    for row_number, row in enumerate(rows, start=2):
        symbol = row["gene_name"].strip()
        try:
            declared_n = int(row["n_gene_ids_for_symbol"])
        except ValueError as error:
            raise AnnotationError(
                f"invalid n_gene_ids_for_symbol at identity row {row_number}"
            ) from error
        observed_n = len(by_symbol[symbol])
        expected_status = "unique" if observed_n == 1 else "duplicated"
        require(
            declared_n == observed_n,
            f"symbol cardinality metadata mismatch for {symbol}: declared={declared_n}, observed={observed_n}",
        )
        require(
            row["symbol_mapping_status"] == expected_status,
            f"symbol_mapping_status mismatch for {symbol}: {row['symbol_mapping_status']!r}",
        )
    return by_base, dict(by_symbol)


def validate_membership(
    rows: list[dict[str, str]], expected_programs: int, expected_rows: int
) -> dict[str, list[dict[str, str]]]:
    require(
        len(rows) == expected_rows,
        f"expected {expected_rows} membership rows, found {len(rows)}",
    )
    by_program: dict[str, list[dict[str, str]]] = defaultdict(list)
    member_keys: set[tuple[str, str]] = set()
    for row_number, row in enumerate(rows, start=2):
        uid = row["program_uid"].strip()
        gene = row["canonical_gene"].strip()
        require(uid != "", f"blank program_uid at membership row {row_number}")
        require(gene != "", f"blank canonical_gene at membership row {row_number}")
        key = (uid, gene)
        require(key not in member_keys, f"duplicate frozen program member: {key}")
        member_keys.add(key)
        parse_positive_float(row["source_weight"], f"source_weight row {row_number}")
        parse_positive_float(
            row["original_l1_weight"], f"original_l1_weight row {row_number}"
        )
        require(
            row["mapped_symbol_status"] in MAPPED_STATUSES | {UNMAPPED_STATUS},
            f"unknown mapped_symbol_status at row {row_number}: {row['mapped_symbol_status']!r}",
        )
        by_program[uid].append(row)
    require(
        len(by_program) == expected_programs,
        f"expected {expected_programs} programs, found {len(by_program)}",
    )
    for uid, program_rows in by_program.items():
        for field in ("cell_type", "module", "membership_sha256"):
            values = {row[field] for row in program_rows}
            require(len(values) == 1, f"{field} varies within {uid}: {sorted(values)}")
        total = sum(float(row["original_l1_weight"]) for row in program_rows)
        require(
            abs(total - 1.0) <= 1e-12,
            f"original L1 weights do not sum to one for {uid}: {total:.17g}",
        )
    return dict(by_program)


def map_member(
    row: dict[str, str],
    by_base: dict[str, Identity],
    by_symbol: dict[str, list[Identity]],
) -> tuple[str, str, Identity | None]:
    registry_status = row["mapped_symbol_status"]
    mapped = row["mapped_symbol"].strip()
    canonical = normalize_base(row["canonical_gene"])

    if registry_status == "gencode_v49_unique_symbol_confirmed":
        matches = by_symbol.get(mapped, [])
        require(
            len(matches) == 1,
            f"frozen unique-symbol mapping drift for {row['program_uid']}:{canonical}",
        )
        return "unique", "frozen_unique_symbol_confirmed", matches[0]

    if registry_status == "gencode_v49_unambiguous_ensembl_to_symbol":
        identity = by_base.get(canonical)
        require(
            identity is not None,
            f"frozen Ensembl mapping missing from identity table: {canonical}",
        )
        require(
            mapped in {identity.gene_name, identity.ensembl_base},
            f"frozen Ensembl/display mapping drift for {canonical}: {mapped!r}",
        )
        return "unique", "frozen_unambiguous_ensembl_confirmed", identity

    matches = (
        [by_base[canonical]]
        if ENSG_BASE.fullmatch(canonical) and canonical in by_base
        else by_symbol.get(canonical, [])
    )
    if not matches:
        return "missing", "no_gencode_identity_for_frozen_member", None
    if len(matches) > 1:
        return "ambiguous", "multiple_gencode_identities_for_frozen_member", None
    return "ambiguous", "frozen_registry_unresolved_despite_single_identity", None


def annotate(
    membership_rows: list[dict[str, str]],
    by_base: dict[str, Identity],
    by_symbol: dict[str, list[Identity]],
) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    resolved_keys: set[tuple[str, str]] = set()
    for source in membership_rows:
        status, reason, identity = map_member(source, by_base, by_symbol)
        row = dict(source)
        row.update(
            {
                "identity_mapping_status": status,
                "identity_mapping_reason": reason,
                "gene_id_versioned": identity.gene_id if identity else "",
                "gene_id_base": identity.ensembl_base if identity else "",
                "identity_gene_name": identity.gene_name if identity else "",
                "gene_type": identity.gene_biotype if identity else "",
                "is_uniquely_mapped_lncrna": "true"
                if identity and identity.gene_biotype == "lncRNA"
                else "false",
            }
        )
        if identity:
            key = (source["program_uid"], identity.ensembl_base)
            require(
                key not in resolved_keys,
                f"two frozen members resolve to one stable gene within a program: {key}",
            )
            resolved_keys.add(key)
        result.append(row)
    return result


def summarize(
    annotated: list[dict[str, str]], threshold: float
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_program: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in annotated:
        by_program[row["program_uid"]].append(row)

    program_rows: list[dict[str, object]] = []
    audit_counter: Counter[tuple[str, str]] = Counter()
    audit_weight: defaultdict[tuple[str, str], float] = defaultdict(float)
    for uid in sorted(by_program):
        rows = by_program[uid]
        total_weight = sum(float(row["original_l1_weight"]) for row in rows)
        unique = [row for row in rows if row["identity_mapping_status"] == "unique"]
        lncrna = [row for row in unique if row["gene_type"] == "lncRNA"]
        ambiguous = [
            row for row in rows if row["identity_mapping_status"] == "ambiguous"
        ]
        missing = [row for row in rows if row["identity_mapping_status"] == "missing"]
        lncrna_weight = sum(float(row["original_l1_weight"]) for row in lncrna)
        for row in rows:
            key = (row["identity_mapping_status"], row["identity_mapping_reason"])
            audit_counter[key] += 1
            audit_weight[key] += float(row["original_l1_weight"])
        program_rows.append(
            {
                "program_uid": uid,
                "cell_type": rows[0]["cell_type"],
                "module": rows[0]["module"],
                "membership_sha256": rows[0]["membership_sha256"],
                "n_frozen_members": len(rows),
                "total_original_l1_weight": f"{total_weight:.17g}",
                "n_uniquely_mapped_members": len(unique),
                "uniquely_mapped_l1_weight": f"{sum(float(row['original_l1_weight']) for row in unique):.17g}",
                "n_uniquely_mapped_lncrna": len(
                    {row["gene_id_base"] for row in lncrna}
                ),
                "lncrna_original_l1_weight": f"{lncrna_weight:.17g}",
                "lncrna_l1_fraction": f"{lncrna_weight / total_weight:.17g}",
                "n_ambiguous_members": len(ambiguous),
                "ambiguous_original_l1_weight": f"{sum(float(row['original_l1_weight']) for row in ambiguous):.17g}",
                "n_missing_members": len(missing),
                "missing_original_l1_weight": f"{sum(float(row['original_l1_weight']) for row in missing):.17g}",
                "sensitivity_trigger_threshold": f"{threshold:.17g}",
                "requires_leave_all_lncrna_out_sensitivity": "true"
                if lncrna_weight / total_weight >= threshold
                else "false",
            }
        )

    audit_rows = [
        {
            "identity_mapping_status": status,
            "identity_mapping_reason": reason,
            "n_members": audit_counter[(status, reason)],
            "original_l1_weight_sum_across_programs": f"{audit_weight[(status, reason)]:.17g}",
        }
        for status, reason in sorted(audit_counter)
    ]
    return program_rows, audit_rows


def source_manifest_rows(membership: Path, identity: Path) -> list[dict[str, object]]:
    producer = Path(__file__).resolve()
    return [
        {
            "role": role,
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for role, path in (
            ("frozen_program_membership", membership),
            ("stable_gene_identity", identity),
            ("producer_script", producer),
        )
    ]


def run(args: argparse.Namespace) -> None:
    membership = args.membership.resolve()
    identity = args.identity.resolve()
    output = args.output.resolve()
    require(
        math.isfinite(args.lncrna_weight_threshold)
        and 0.0 <= args.lncrna_weight_threshold <= 1.0,
        "lncRNA threshold must be finite and between zero and one",
    )
    expected_membership_sha256 = getattr(args, "expected_membership_sha256", None)
    if expected_membership_sha256:
        require(
            sha256(membership) == expected_membership_sha256,
            "frozen program-membership SHA256 drift",
        )
    membership_fields, membership_rows = read_tsv(membership, MEMBERSHIP_REQUIRED)
    _, identity_rows = read_tsv(identity, IDENTITY_REQUIRED)
    by_base, by_symbol = build_identity_indexes(identity_rows)
    validate_membership(membership_rows, args.expected_programs, args.expected_rows)
    annotated = annotate(membership_rows, by_base, by_symbol)
    program_rows, audit_rows = summarize(annotated, args.lncrna_weight_threshold)

    require(not output.exists(), f"refusing to overwrite output directory: {output}")
    output.mkdir(parents=True, exist_ok=False)
    annotation_fields = membership_fields + [
        "identity_mapping_status",
        "identity_mapping_reason",
        "gene_id_versioned",
        "gene_id_base",
        "identity_gene_name",
        "gene_type",
        "is_uniquely_mapped_lncrna",
    ]
    write_tsv(output / OUTPUT_FILES[0], annotation_fields, annotated)
    write_tsv(output / OUTPUT_FILES[1], list(program_rows[0]), program_rows)
    write_tsv(output / OUTPUT_FILES[2], list(audit_rows[0]), audit_rows)
    write_tsv(
        output / OUTPUT_FILES[3],
        ["role", "path", "size_bytes", "sha256"],
        source_manifest_rows(membership, identity),
    )
    execution = {
        "contract": "frozen_hotspot_lncrna_annotation_v1",
        "expected_programs": args.expected_programs,
        "expected_membership_rows": args.expected_rows,
        "lncrna_weight_threshold": args.lncrna_weight_threshold,
        "python": sys.version,
        "platform": platform.platform(),
        "program_scores_computed": False,
        "program_memberships_modified": False,
        "program_weights_modified": False,
    }
    with (output / OUTPUT_FILES[4]).open("x", encoding="utf-8") as handle:
        json.dump(execution, handle, indent=2, sort_keys=True)
        handle.write("\n")
    output_rows = [
        {
            "relative_path": name,
            "size_bytes": (output / name).stat().st_size,
            "sha256": sha256(output / name),
        }
        for name in OUTPUT_FILES
    ]
    write_tsv(
        output / "output_manifest.tsv",
        ["relative_path", "size_bytes", "sha256"],
        output_rows,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--identity", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-programs", type=int, default=117)
    parser.add_argument("--expected-rows", type=int, default=7093)
    parser.add_argument(
        "--expected-membership-sha256",
        default=FROZEN_MEMBERSHIP_SHA256,
        help="Expected byte-level SHA256 of the frozen membership TSV",
    )
    parser.add_argument("--lncrna-weight-threshold", type=float, default=0.05)
    return parser.parse_args()


if __name__ == "__main__":
    try:
        run(parse_args())
    except AnnotationError as error:
        raise SystemExit(f"ERROR: {error}") from error
